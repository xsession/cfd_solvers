#include "cfd/solvers/optics/polarization.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::optics {

StokesVector stokes_from_jones(JonesVector j){
    const double ex2=std::norm(j.x),ey2=std::norm(j.y);
    const auto cross=j.x*std::conj(j.y);
    return {ex2+ey2,ex2-ey2,2.0*cross.real(),-2.0*cross.imag()};
}

FresnelCoefficients fresnel_dielectric(double n1,double n2,double cos_i){
    if(!(n1>0.0)||!(n2>0.0)) throw std::invalid_argument("Fresnel refractive index must be positive");
    cos_i=std::clamp(cos_i,0.0,1.0);
    const double eta=n1/n2;
    const double sin2_t=eta*eta*(1.0-cos_i*cos_i);
    if(sin2_t>1.0) return {0.0,1.0,1.0,0.0,0.0,1.0,1.0,true};
    const double cos_t=std::sqrt(std::max(0.0,1.0-sin2_t));
    const double rs=(n1*cos_i-n2*cos_t)/(n1*cos_i+n2*cos_t);
    const double rp=(n2*cos_i-n1*cos_t)/(n2*cos_i+n1*cos_t);
    const double ts=2.0*n1*cos_i/(n1*cos_i+n2*cos_t);
    const double tp=2.0*n1*cos_i/(n2*cos_i+n1*cos_t);
    return {cos_t,rs,rp,ts,tp,rs*rs,rp*rp,false};
}

double thin_film_reflectance_normal(double n0,double ns,double wavelength,const std::vector<ThinFilmLayer>& layers){
    if(!(n0>0.0)||!(ns>0.0)||!(wavelength>0.0)) throw std::invalid_argument("invalid thin-film optical constants");
    using C=std::complex<double>; const C I(0.0,1.0);
    C m00(1.0,0.0),m01(0.0,0.0),m10(0.0,0.0),m11(1.0,0.0);
    for(const auto&layer:layers){
        if(!(layer.refractive_index>0.0)||!(layer.thickness_nm>=0.0)) throw std::invalid_argument("invalid thin-film layer");
        const double delta=2.0*std::numbers::pi*layer.refractive_index*layer.thickness_nm/wavelength;
        const C a=std::cos(delta),b=I*std::sin(delta)/layer.refractive_index,c=I*layer.refractive_index*std::sin(delta),d=a;
        const C r00=m00*a+m01*c,r01=m00*b+m01*d,r10=m10*a+m11*c,r11=m10*b+m11*d;
        m00=r00;m01=r01;m10=r10;m11=r11;
    }
    const C numerator=n0*m00+n0*ns*m01-m10-ns*m11;
    const C denominator=n0*m00+n0*ns*m01+m10+ns*m11;
    return std::norm(numerator/denominator);
}

} // namespace cfd::optics
