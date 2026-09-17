#include "cfd/em/wave_port.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::em {
namespace {
constexpr double epsilon0=8.8541878128e-12;
constexpr double mu0=1.25663706212e-6;

double cos2_integral(double length,std::size_t index){return index==0U?length:0.5*length;}
double sin2_integral(double length,std::size_t index){return index==0U?0.0:0.5*length;}
}

std::vector<RectangularWaveguideMode> rectangular_waveguide_modes(double a,double b,double frequency,
                                                                   double eps_r,double mu_r,std::size_t count){
    if(!(a>0.0)||!(b>0.0)||!(frequency>0.0)||!(eps_r>0.0)||!(mu_r>0.0)||count==0U)throw std::invalid_argument("invalid rectangular wave-port controls");
    const double epsilon=epsilon0*eps_r,mu=mu0*mu_r,omega=2.0*std::numbers::pi*frequency,k=omega*std::sqrt(mu*epsilon),velocity=1.0/std::sqrt(mu*epsilon);
    const std::size_t limit=std::max<std::size_t>(4U,count+3U);std::vector<RectangularWaveguideMode> modes;
    for(std::size_t m=0;m<=limit;++m)for(std::size_t n=0;n<=limit;++n){
        if(m==0U&&n==0U)continue;
        const double kx=std::numbers::pi*static_cast<double>(m)/a,ky=std::numbers::pi*static_cast<double>(n)/b,kc=std::hypot(kx,ky),fc=velocity*kc/(2.0*std::numbers::pi);
        RectangularWaveguideMode te;te.family=RectangularWaveguideMode::Family::te;te.m=m;te.n=n;te.frequency_hz=frequency;te.cutoff_frequency_hz=fc;te.propagating=k>kc;
        if(te.propagating){const double beta=std::sqrt(k*k-kc*kc);te.propagation_constant_rad_per_m={beta,0.0};te.wave_impedance_ohm={omega*mu/beta,0.0};}else{const double alpha=std::sqrt(kc*kc-k*k);te.propagation_constant_rad_per_m={0.0,-alpha};te.wave_impedance_ohm=Complex{0.0,-omega*mu/alpha};}
        const double integral=ky*ky*cos2_integral(a,m)*sin2_integral(b,n)+kx*kx*sin2_integral(a,m)*cos2_integral(b,n);
        te.electric_scale=te.propagating?std::sqrt(2.0*te.wave_impedance_ohm.real()/integral):1.0/std::sqrt(integral);modes.push_back(te);
        if(m>0U&&n>0U){RectangularWaveguideMode tm;tm.family=RectangularWaveguideMode::Family::tm;tm.m=m;tm.n=n;tm.frequency_hz=frequency;tm.cutoff_frequency_hz=fc;tm.propagating=k>kc;if(tm.propagating){const double beta=std::sqrt(k*k-kc*kc);tm.propagation_constant_rad_per_m={beta,0.0};tm.wave_impedance_ohm={beta/(omega*epsilon),0.0};}else{const double alpha=std::sqrt(kc*kc-k*k);tm.propagation_constant_rad_per_m={0.0,-alpha};tm.wave_impedance_ohm=Complex{0.0,alpha/(omega*epsilon)};}const double ti=kx*kx*cos2_integral(a,m)*sin2_integral(b,n)+ky*ky*sin2_integral(a,m)*cos2_integral(b,n);tm.electric_scale=tm.propagating?std::sqrt(2.0*tm.wave_impedance_ohm.real()/ti):1.0/std::sqrt(ti);modes.push_back(tm);}
    }
    std::sort(modes.begin(),modes.end(),[](const auto& lhs,const auto& rhs){if(lhs.cutoff_frequency_hz!=rhs.cutoff_frequency_hz)return lhs.cutoff_frequency_hz<rhs.cutoff_frequency_hz;if(lhs.family!=rhs.family)return lhs.family==RectangularWaveguideMode::Family::te;return std::pair{lhs.m,lhs.n}<std::pair{rhs.m,rhs.n};});
    if(modes.size()>count)modes.resize(count);
    return modes;
}

TransverseModeFields rectangular_waveguide_transverse_fields(const RectangularWaveguideMode& mode,double a,double b,double x,double y){
    if(!(a>0.0)||!(b>0.0)||x<0.0||x>a||y<0.0||y>b)throw std::invalid_argument("invalid rectangular wave-port field coordinate");
    const double kx=std::numbers::pi*static_cast<double>(mode.m)/a,ky=std::numbers::pi*static_cast<double>(mode.n)/b;double ex=0.0,ey=0.0;
    if(mode.family==RectangularWaveguideMode::Family::te){ex=ky*std::cos(kx*x)*std::sin(ky*y);ey=-kx*std::sin(kx*x)*std::cos(ky*y);}else{ex=kx*std::cos(kx*x)*std::sin(ky*y);ey=ky*std::sin(kx*x)*std::cos(ky*y);}
    TransverseModeFields field;field.electric_v_per_m={mode.electric_scale*ex,mode.electric_scale*ey};if(std::abs(mode.wave_impedance_ohm)>0.0){field.magnetic_a_per_m[0]=-field.electric_v_per_m[1]/mode.wave_impedance_ohm;field.magnetic_a_per_m[1]=field.electric_v_per_m[0]/mode.wave_impedance_ohm;}return field;
}

double rectangular_mode_power_w(const RectangularWaveguideMode& mode,double a,double b,std::size_t nx,std::size_t ny){
    if(nx==0U||ny==0U)throw std::invalid_argument("wave-port power quadrature requires non-zero sampling");
    const double dx=a/static_cast<double>(nx),dy=b/static_cast<double>(ny);double power=0.0;for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i){const auto f=rectangular_waveguide_transverse_fields(mode,a,b,(static_cast<double>(i)+0.5)*dx,(static_cast<double>(j)+0.5)*dy);const Complex pz=f.electric_v_per_m[0]*std::conj(f.magnetic_a_per_m[1])-f.electric_v_per_m[1]*std::conj(f.magnetic_a_per_m[0]);power+=0.5*pz.real()*dx*dy;}return power;
}

} // namespace cfd::em
