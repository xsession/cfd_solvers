#include "cfd/fvm/advanced_models.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
double spalart_allmaras_eddy_viscosity(double nu,double nt,double cv1){if(!(nu>0)||!(nt>=0)||!(cv1>0))throw std::invalid_argument("invalid SA state");const double chi=nt/nu,chi3=chi*chi*chi,cv3=cv1*cv1*cv1;return nt*chi3/(chi3+cv3);}
double k_epsilon_eddy_viscosity(double rho,double k,double eps,double cmu){if(!(rho>0)||!(k>=0)||!(eps>0)||!(cmu>0))throw std::invalid_argument("invalid k-epsilon state");return rho*cmu*k*k/eps;}
double k_omega_sst_eddy_viscosity(double rho,double k,double omega,double strain,double f2,double a1){if(!(rho>0)||!(k>=0)||!(omega>0)||!(a1>0)||!(f2>=0)||!(strain>=0))throw std::invalid_argument("invalid SST state");return rho*a1*k/std::max(a1*omega,strain*f2);}
double des_length_scale(double wall,double grid,double c){if(!(wall>=0)||!(grid>0)||!(c>0))throw std::invalid_argument("invalid DES scale");return std::min(wall,c*grid);}
ReynoldsStress2D boussinesq_reynolds_stress(double rho,double k,double mut,double dux,double dvy,double shear){if(!(rho>0)||!(k>=0)||!(mut>=0))throw std::invalid_argument("invalid RSM state");const double iso=2.0*rho*k/3.0;return {iso-2.0*mut*dux,iso-2.0*mut*dvy,-mut*shear};}
PhaseChangeState linear_mushy_phase_change(double t,double ts,double tl,double cps,double cpl,double latent){if(!(tl>ts)||!(cps>0)||!(cpl>0)||!(latent>=0))throw std::invalid_argument("invalid phase-change material");double f=std::clamp((t-ts)/(tl-ts),0.0,1.0);double cp=(1-f)*cps+f*cpl+(t>ts&&t<tl?latent/(tl-ts):0.0);double h=cps*std::min(t,ts);if(t>ts){double tm=std::min(t,tl)-ts;h+=0.5*(cps+cpl)*tm+latent*std::clamp(tm/(tl-ts),0.0,1.0);}if(t>tl)h+=cpl*(t-tl);return {f,cp,h};}
double drift_flux_slip_velocity(double rhoc,double rhod,double d,double mu,double g){if(!(rhoc>0)||!(rhod>0)||!(d>0)||!(mu>0))throw std::invalid_argument("invalid drift-flux state");return (rhoc-rhod)*g*d*d/(18.0*mu);}
double euler_euler_drag_source(double uc,double ud,double cd,double area){if(!(cd>=0)||!(area>=0))throw std::invalid_argument("invalid Euler-Euler drag");return cd*area*(ud-uc)*std::abs(ud-uc);}
double lubrication_thin_film_flux(double h,double dpdx,double mu,double us){if(!(h>=0)||!(mu>0))throw std::invalid_argument("invalid thin-film state");return -h*h*h*dpdx/(12.0*mu)+0.5*h*us;}
} // namespace cfd::fvm
