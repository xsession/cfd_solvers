#include "cfd/fem/material_laws.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::fem {
LinearElastic1D::LinearElastic1D(double E):E_(E){if(!(E>0))throw std::invalid_argument("invalid Young modulus");}MaterialPoint1D LinearElastic1D::update(double e){return {E_*e,E_,0,0};}
BilinearPlasticity1D::BilinearPlasticity1D(double E,double y,double H):E_(E),yield_(y),H_(H){if(!(E>0)||!(y>0)||H<0)throw std::invalid_argument("invalid plastic material");}MaterialPoint1D BilinearPlasticity1D::update(double e){double trial=E_*(e-plastic_),f=std::abs(trial)-(yield_+H_*alpha_);double tangent=E_;if(f>0){double dg=f/(E_+H_),sgn=trial>=0?1.0:-1.0;plastic_+=dg*sgn;alpha_+=dg;trial-=E_*dg*sgn;tangent=E_*H_/(E_+H_);}return {trial,tangent,plastic_,alpha_};}
NeoHookeanResult compressible_neo_hookean(std::array<double,3> l,double mu,double lambda){if(!(mu>0)||lambda<0)throw std::invalid_argument("invalid Neo-Hookean moduli");double J=l[0]*l[1]*l[2];for(double v:l)if(!(v>0))throw std::invalid_argument("invalid principal stretch");double lnJ=std::log(J),i1=l[0]*l[0]+l[1]*l[1]+l[2]*l[2];NeoHookeanResult r;r.energy_density=0.5*mu*(i1-3.0)-mu*lnJ+0.5*lambda*lnJ*lnJ;for(int i=0;i<3;++i)r.nominal_stress[i]=mu*(l[i]-1.0/l[i])+lambda*lnJ/l[i];return r;}
double penalty_contact_pressure(double gap,double k){if(!(k>0)||!std::isfinite(gap))throw std::invalid_argument("invalid contact controls");return gap<0?-k*gap:0.0;}
}
