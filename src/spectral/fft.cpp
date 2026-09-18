#include "cfd/spectral/fft.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::spectral {

bool is_power_of_two(std::size_t n) noexcept { return n!=0U && (n&(n-1U))==0U; }

void fft_inplace(std::span<Complex> a,bool inverse){
    const std::size_t n=a.size();
    if(!is_power_of_two(n)) throw std::invalid_argument("FFT length must be a non-zero power of two");
    for(std::size_t i=1,j=0;i<n;++i){
        std::size_t bit=n>>1U;
        for(;j&bit;bit>>=1U) j^=bit;
        j^=bit;
        if(i<j) std::swap(a[i],a[j]);
    }
    for(std::size_t len=2;len<=n;len<<=1U){
        const double angle=(inverse?2.0:-2.0)*std::numbers::pi/static_cast<double>(len);
        const Complex root{std::cos(angle),std::sin(angle)};
        for(std::size_t i=0;i<n;i+=len){
            Complex w{1.0,0.0};
            const std::size_t half=len>>1U;
            for(std::size_t j=0;j<half;++j){
                const Complex u=a[i+j];
                const Complex v=a[i+j+half]*w;
                a[i+j]=u+v;
                a[i+j+half]=u-v;
                w*=root;
            }
        }
    }
    if(inverse){
        const double scale=1.0/static_cast<double>(n);
        for(auto& v:a) v*=scale;
    }
}

void fft2_inplace(std::span<Complex> values,std::size_t nx,std::size_t ny,bool inverse){
    if(nx==0U||ny==0U||values.size()!=nx*ny||!is_power_of_two(nx)||!is_power_of_two(ny))
        throw std::invalid_argument("2-D FFT requires power-of-two dimensions matching the array size");
    std::vector<Complex> scratch(std::max(nx,ny));
    for(std::size_t j=0;j<ny;++j){
        auto row=values.subspan(j*nx,nx);
        fft_inplace(row,inverse);
    }
    for(std::size_t i=0;i<nx;++i){
        for(std::size_t j=0;j<ny;++j) scratch[j]=values[j*nx+i];
        fft_inplace(std::span<Complex>(scratch.data(),ny),inverse);
        for(std::size_t j=0;j<ny;++j) values[j*nx+i]=scratch[j];
    }
}

std::vector<double> wavenumbers(std::size_t n,double spacing){
    if(!is_power_of_two(n)||!(spacing>0.0)||!std::isfinite(spacing)) throw std::invalid_argument("invalid spectral grid");
    const double length=static_cast<double>(n)*spacing;
    const double factor=2.0*std::numbers::pi/length;
    std::vector<double> k(n);
    for(std::size_t i=0;i<n;++i){
        const auto mode=(i<=n/2U)?static_cast<long long>(i):static_cast<long long>(i)-static_cast<long long>(n);
        k[i]=factor*static_cast<double>(mode);
    }
    return k;
}

std::vector<double> derivative_periodic_1d(std::span<const double> field,double spacing){
    if(!is_power_of_two(field.size())) throw std::invalid_argument("periodic spectral derivative requires a power-of-two grid");
    std::vector<Complex> f(field.size());
    for(std::size_t i=0;i<field.size();++i) f[i]=field[i];
    fft_inplace(f,false);
    const auto k=wavenumbers(field.size(),spacing);
    const Complex imag{0.0,1.0};
    for(std::size_t i=0;i<f.size();++i) f[i]*=imag*((i==f.size()/2U)?0.0:k[i]);
    fft_inplace(f,true);
    std::vector<double> out(f.size());
    for(std::size_t i=0;i<f.size();++i) out[i]=f[i].real();
    return out;
}

std::vector<double> fractional_laplacian_periodic_1d(std::span<const double> field,double spacing,double order){
    if(!is_power_of_two(field.size())||!(order>=0.0)||!std::isfinite(order)) throw std::invalid_argument("invalid fractional Laplacian request");
    std::vector<Complex> f(field.size());
    for(std::size_t i=0;i<field.size();++i) f[i]=field[i];
    fft_inplace(f,false);
    const auto k=wavenumbers(field.size(),spacing);
    for(std::size_t i=0;i<f.size();++i){
        const double factor=(k[i]==0.0)?(order==0.0?1.0:0.0):std::pow(std::abs(k[i]),2.0*order);
        f[i]*=factor;
    }
    fft_inplace(f,true);
    std::vector<double> out(f.size());
    for(std::size_t i=0;i<f.size();++i) out[i]=f[i].real();
    return out;
}

} // namespace cfd::spectral
