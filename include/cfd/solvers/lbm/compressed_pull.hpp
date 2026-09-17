#pragma once

#include "cfd/core/parallel.hpp"
#include "cfd/core/soa_field.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace cfd::lbm {

struct IeeeHalfCodec {
    static std::uint16_t encode(float value, int = 0) noexcept {
        if (std::isnan(value)) return 0x7e00U;
        const bool negative = std::signbit(value);
        const float a = std::abs(value);
        const std::uint16_t sign = negative ? 0x8000U : 0U;
        if (std::isinf(a)) return static_cast<std::uint16_t>(sign | 0x7c00U);
        if (a == 0.0F) return sign;
        if (a < std::ldexp(1.0F, -14)) {
            long mantissa = std::lround(std::ldexp(a, 24));
            if (mantissa <= 0) return sign;
            if (mantissa > 1023) mantissa = 1023;
            return static_cast<std::uint16_t>(sign | static_cast<std::uint16_t>(mantissa));
        }
        int exponent = static_cast<int>(std::floor(std::log2(a)));
        if (exponent > 15) return static_cast<std::uint16_t>(sign | 0x7c00U);
        float normalized = std::ldexp(a, -exponent);
        long mantissa = std::lround((normalized - 1.0F) * 1024.0F);
        if (mantissa == 1024) {
            mantissa = 0;
            ++exponent;
            if (exponent > 15) return static_cast<std::uint16_t>(sign | 0x7c00U);
        }
        const auto exp_bits = static_cast<std::uint16_t>(exponent + 15);
        return static_cast<std::uint16_t>(sign | (exp_bits << 10U) | static_cast<std::uint16_t>(mantissa));
    }

    static float decode(std::uint16_t bits, int = 0) noexcept {
        const float sign = (bits & 0x8000U) ? -1.0F : 1.0F;
        const unsigned exponent = (bits >> 10U) & 0x1fU;
        const unsigned mantissa = bits & 0x03ffU;
        if (exponent == 0x1fU) {
            return mantissa ? std::numeric_limits<float>::quiet_NaN()
                            : sign * std::numeric_limits<float>::infinity();
        }
        if (exponent == 0U) return sign * std::ldexp(static_cast<float>(mantissa), -24);
        return sign * std::ldexp(1.0F + static_cast<float>(mantissa) / 1024.0F,
                                 static_cast<int>(exponent) - 15);
    }
};

struct BFloat16Codec {
    static std::uint16_t encode(float value, int = 0) noexcept {
        std::uint32_t bits = std::bit_cast<std::uint32_t>(value);
        const std::uint32_t lsb = (bits >> 16U) & 1U;
        bits += 0x7fffU + lsb;
        return static_cast<std::uint16_t>(bits >> 16U);
    }
    static float decode(std::uint16_t bits, int = 0) noexcept {
        return std::bit_cast<float>(static_cast<std::uint32_t>(bits) << 16U);
    }
};

template<class Descriptor>
struct ShiftedHalfCodec {
    static std::uint16_t encode(float value, int direction) noexcept {
        return IeeeHalfCodec::encode(value - Descriptor::weight(direction));
    }
    static float decode(std::uint16_t bits, int direction) noexcept {
        return IeeeHalfCodec::decode(bits) + Descriptor::weight(direction);
    }
};

template<class Descriptor>
struct PackedDeviation16Codec {
    static constexpr float scale = 32768.0F;
    static std::uint16_t encode(float value, int direction) noexcept {
        const float deviation = value - Descriptor::weight(direction);
        const long quantized = std::lround(std::clamp(deviation * scale, -32768.0F, 32767.0F));
        return static_cast<std::uint16_t>(static_cast<std::int16_t>(quantized));
    }
    static float decode(std::uint16_t bits, int direction) noexcept {
        const auto quantized = static_cast<std::int16_t>(bits);
        return Descriptor::weight(direction) + static_cast<float>(quantized) / scale;
    }
};

template<class Descriptor, class Codec>
class CompressedPopulationPullSolver {
public:
    static constexpr int q = Descriptor::q;
    struct Config { std::size_t nx{64}, ny{64}, nz{64}; float tau{0.6F}; };
    struct Macroscopic { std::vector<float> rho, ux, uy, uz; };

    explicit CompressedPopulationPullSolver(Config config)
        : config_(normalize(config)), cells_(config_.nx * config_.ny * config_.nz),
          f_(cells_), next_(cells_) {
        if (!(config_.tau > 0.5F)) throw std::invalid_argument("compressed LBM tau must be > 0.5");
        initialize_uniform();
    }

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F, float uz = 0.0F) {
        initialize([=](std::size_t, std::size_t, std::size_t){ return std::array<float,4>{rho,ux,uy,uz}; });
    }
    void initialize_taylor_green(float amplitude = 0.03F) {
        constexpr float pi = 3.14159265358979323846F;
        initialize([=,this](std::size_t x,std::size_t y,std::size_t z){
            const float xf=(static_cast<float>(x)+0.5F)/static_cast<float>(config_.nx);
            const float yf=(static_cast<float>(y)+0.5F)/static_cast<float>(config_.ny);
            float zf=1.0F;
            if constexpr(Descriptor::dimensions==3) zf=std::cos(2*pi*(static_cast<float>(z)+0.5F)/static_cast<float>(config_.nz));
            return std::array<float,4>{1.0F, amplitude*std::sin(2*pi*xf)*std::cos(2*pi*yf)*zf,
                -amplitude*std::cos(2*pi*xf)*std::sin(2*pi*yf)*zf,0.0F};
        });
    }
    void step(std::size_t count=1){ while(count--) step_once(); }

    [[nodiscard]] Macroscopic compute_macroscopic() const {
        Macroscopic out{std::vector<float>(cells_),std::vector<float>(cells_),std::vector<float>(cells_),std::vector<float>(cells_)};
        cfd::core::parallel_for(cells_,[&](std::size_t n){
            float rho=0,mx=0,my=0,mz=0;
            for(int d=0;d<q;++d){const float v=load(f_,d,n);rho+=v;mx+=v*Descriptor::cx(d);my+=v*Descriptor::cy(d);mz+=v*Descriptor::cz(d);}
            out.rho[n]=rho;if(rho>0){out.ux[n]=mx/rho;out.uy[n]=my/rho;out.uz[n]=mz/rho;}
        });
        return out;
    }
    [[nodiscard]] double mass() const { double sum=0; for(std::size_t n=0;n<cells_;++n)for(int d=0;d<q;++d)sum+=load(f_,d,n); return sum; }
    [[nodiscard]] std::size_t population_bytes() const noexcept { return static_cast<std::size_t>(q)*cells_*sizeof(std::uint16_t); }

private:
    Config config_{}; std::size_t cells_{};
    cfd::core::StaticSoA<std::uint16_t,static_cast<std::size_t>(q)> f_,next_;
    static Config normalize(Config c){if(c.nx<2||c.ny<2)throw std::invalid_argument("compressed LBM nx/ny >=2");if constexpr(Descriptor::dimensions==2)c.nz=1;else if(c.nz<2)throw std::invalid_argument("compressed 3D LBM nz >=2");return c;}
    static float equilibrium(int d,float rho,float ux,float uy,float uz){const float cu=3.0F*(Descriptor::cx(d)*ux+Descriptor::cy(d)*uy+Descriptor::cz(d)*uz);const float uu=1.5F*(ux*ux+uy*uy+uz*uz);return Descriptor::weight(d)*rho*(1.0F+cu+0.5F*cu*cu-uu);}
    static float load(const cfd::core::StaticSoA<std::uint16_t,static_cast<std::size_t>(q)>& a,int d,std::size_t n){return Codec::decode(a(static_cast<std::size_t>(d),n),d);}
    static void store(cfd::core::StaticSoA<std::uint16_t,static_cast<std::size_t>(q)>& a,int d,std::size_t n,float v){a(static_cast<std::size_t>(d),n)=Codec::encode(v,d);}
    std::size_t source(std::size_t x,std::size_t y,std::size_t z,int d)const{auto shift=[](std::size_t c,int delta,std::size_t n){if(delta>0)return c+1==n?0U:c+1;if(delta<0)return c==0?n-1:c-1;return c;};const auto xs=shift(x,-Descriptor::cx(d),config_.nx),ys=shift(y,-Descriptor::cy(d),config_.ny),zs=shift(z,-Descriptor::cz(d),config_.nz);return(zs*config_.ny+ys)*config_.nx+xs;}
    template<class Init>void initialize(Init init){cfd::core::parallel_for(cells_,[&](std::size_t n){const auto x=n%config_.nx,y=(n/config_.nx)%config_.ny,z=n/(config_.nx*config_.ny);const auto s=init(x,y,z);for(int d=0;d<q;++d){const float v=equilibrium(d,s[0],s[1],s[2],s[3]);store(f_,d,n,v);store(next_,d,n,v);}});}
    void step_once(){const float omega=1.0F/config_.tau;cfd::core::parallel_for(cells_,[&](std::size_t n){const auto x=n%config_.nx,y=(n/config_.nx)%config_.ny,z=n/(config_.nx*config_.ny);std::array<float,q> fin{};float rho=0,mx=0,my=0,mz=0;for(int d=0;d<q;++d){const float v=load(f_,d,source(x,y,z,d));fin[static_cast<std::size_t>(d)]=v;rho+=v;mx+=v*Descriptor::cx(d);my+=v*Descriptor::cy(d);mz+=v*Descriptor::cz(d);}const float ux=rho>0?mx/rho:0,uy=rho>0?my/rho:0,uz=rho>0?mz/rho:0;for(int d=0;d<q;++d){const auto i=static_cast<std::size_t>(d);store(next_,d,n,fin[i]-omega*(fin[i]-equilibrium(d,rho,ux,uy,uz)));}});std::swap(f_,next_);}
};

} // namespace cfd::lbm
