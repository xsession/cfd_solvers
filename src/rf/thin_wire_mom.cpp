#include "cfd/rf/thin_wire_mom.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::rf {
namespace {
using Complex = std::complex<double>;
constexpr double eps0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;

struct Quadrature {
    std::vector<double> x;
    std::vector<double> w;
};

Quadrature gauss_legendre(std::size_t n) {
    if (n < 2U || n > 32U) throw std::invalid_argument("thin-wire MoM quadrature order must be 2..32");
    Quadrature q; q.x.resize(n); q.w.resize(n);
    const std::size_t half=(n+1U)/2U;
    for(std::size_t i=0;i<half;++i){
        double z=std::cos(std::numbers::pi*(static_cast<double>(i)+0.75)/(static_cast<double>(n)+0.5));
        double previous=0.0,derivative=0.0;
        for(std::size_t iteration=0;iteration<50U;++iteration){
            double p0=1.0,p1=z;
            if(n==1U){derivative=1.0;}
            else{
                for(std::size_t k=2;k<=n;++k){
                    const double kd=static_cast<double>(k);
                    const double p=((2.0*kd-1.0)*z*p1-(kd-1.0)*p0)/kd;
                    p0=p1;p1=p;
                }
                derivative=static_cast<double>(n)*(z*p1-p0)/(z*z-1.0);
                previous=z; z-=p1/derivative;
                if(std::abs(z-previous)<1e-15)break;
                continue;
            }
        }
        // Re-evaluate derivative at the converged root.
        double p0=1.0,p1=z;
        for(std::size_t k=2;k<=n;++k){
            const double kd=static_cast<double>(k);
            const double p=((2.0*kd-1.0)*z*p1-(kd-1.0)*p0)/kd;
            p0=p1;p1=p;
        }
        derivative=static_cast<double>(n)*(z*p1-p0)/(z*z-1.0);
        const double weight=2.0/((1.0-z*z)*derivative*derivative);
        q.x[i]=-z;q.x[n-1U-i]=z;q.w[i]=weight;q.w[n-1U-i]=weight;
    }
    return q;
}

template<class T>
std::vector<T> solve_dense(std::vector<T> a,std::vector<T> b,std::size_t n){
    for(std::size_t col=0;col<n;++col){
        std::size_t pivot=col;double best=std::abs(a[col*n+col]);
        for(std::size_t row=col+1U;row<n;++row){const double v=std::abs(a[row*n+col]);if(v>best){best=v;pivot=row;}}
        if(best<1e-24)throw std::runtime_error("singular thin-wire MoM matrix");
        if(pivot!=col){for(std::size_t j=col;j<n;++j)std::swap(a[pivot*n+j],a[col*n+j]);std::swap(b[pivot],b[col]);}
        const T d=a[col*n+col];for(std::size_t j=col;j<n;++j)a[col*n+j]/=d;b[col]/=d;
        for(std::size_t row=0;row<n;++row){
            if(row==col)continue;
            const T f=a[row*n+col];
            if(std::abs(f)==0.0)continue;
            for(std::size_t j=col;j<n;++j)a[row*n+j]-=f*a[col*n+j];
            b[row]-=f*b[col];
        }
    }
    return b;
}

Complex kernel_integral(double observation_z,double source_center,double dz,double radius,
                        double wave_number,double omega,const Quadrature& q){
    const Complex j{0.0,1.0};
    Complex sum{};
    for(std::size_t k=0;k<q.x.size();++k){
        const double zp=source_center+0.5*dz*q.x[k];
        const double axial=observation_z-zp;
        const double r=std::sqrt(radius*radius+axial*axial);
        const Complex phase=std::exp(-j*wave_number*r);
        const Complex bracket=(Complex{1.0,0.0}+j*wave_number*r)*(2.0*r*r-3.0*radius*radius)
            + wave_number*wave_number*radius*radius*r*r;
        sum+=q.w[k]*phase*bracket/(4.0*std::numbers::pi*std::pow(r,5));
    }
    return j/(omega*eps0)*(0.5*dz)*sum;
}
}

ThinWireMomResult solve_center_fed_thin_wire(const ThinWireMomConfig& config){
    if(!(config.length_m>0.0)||!(config.radius_m>0.0)||!(config.frequency_hz>0.0)
       ||!std::isfinite(config.length_m)||!std::isfinite(config.radius_m)||!std::isfinite(config.frequency_hz)
       ||config.radius_m>=0.1*config.length_m||config.segments<5U||(config.segments%2U)==0U
       ||!(config.feed_voltage_v!=0.0)||!std::isfinite(config.feed_voltage_v))
        throw std::invalid_argument("invalid thin-wire MoM configuration");
    const std::size_t n=config.segments;
    const double dz=config.length_m/static_cast<double>(n);
    if(config.radius_m>=0.5*dz)throw std::invalid_argument("thin-wire MoM requires radius < half segment length");
    const double c=1.0/std::sqrt(mu0*eps0);
    const double omega=2.0*std::numbers::pi*config.frequency_hz;
    const double k=omega/c;
    const auto quadrature=gauss_legendre(config.quadrature_order);
    ThinWireMomResult result;result.segment_center_z_m.resize(n);
    for(std::size_t i=0;i<n;++i)result.segment_center_z_m[i]=-0.5*config.length_m+(static_cast<double>(i)+0.5)*dz;
    std::vector<Complex> matrix(n*n),rhs(n,Complex{});
    for(std::size_t row=0;row<n;++row){
        const double zo=result.segment_center_z_m[row];
        for(std::size_t col=0;col<n;++col)
            matrix[row*n+col]=kernel_integral(zo,result.segment_center_z_m[col],dz,config.radius_m,k,omega,quadrature);
    }
    const std::size_t feed=n/2U;
    // PEC boundary condition: total tangential E vanishes. A delta-gap voltage
    // is represented by a uniform impressed field over the center segment.
    rhs[feed]=-config.feed_voltage_v/dz;
    result.current_a=solve_dense(std::move(matrix),std::move(rhs),n);
    const Complex feed_current=result.current_a[feed];
    if(std::abs(feed_current)<1e-30)throw std::runtime_error("thin-wire MoM zero feed current");
    result.feed_impedance_ohm=config.feed_voltage_v/feed_current;
    if(result.feed_impedance_ohm.real()<0.0){
        // The EFIE sign depends on the impressed-field convention. Normalize
        // current orientation so passive input resistance is positive.
        for(auto& current:result.current_a)current=-current;
        result.feed_impedance_ohm=-result.feed_impedance_ohm;
    }
    result.feed_current_a=std::abs(result.current_a[feed]);
    return result;
}

} // namespace cfd::rf
