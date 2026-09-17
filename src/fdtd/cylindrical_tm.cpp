#include "cfd/solvers/fdtd/cylindrical_tm.hpp"
#include "cfd/core/parallel.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace cfd::fdtd {
namespace {constexpr double eps0=8.8541878128e-12,mu0=1.25663706212e-6,c0=299792458.0;}
CylindricalTM::CylindricalTM(CylindricalTMConfig c):c_(c){if(c.nr<3||c.nz<3||!(c.dr>0)||!(c.dz>0)||!(c.courant>0&&c.courant<=0.5)||!(c.epsilon_r>0)||!(c.mu_r>0))throw std::invalid_argument("invalid cylindrical FDTD config");eps_=eps0*c.epsilon_r;mu_=mu0*c.mu_r;dt_=c.courant/(c0/std::sqrt(c.epsilon_r*c.mu_r)*std::sqrt(1/(c.dr*c.dr)+1/(c.dz*c.dz)));er_.assign((c.nr-1)*c.nz,0.0);ez_.assign(c.nr*(c.nz-1),0.0);hphi_.assign((c.nr-1)*(c.nz-1),0.0);}
void CylindricalTM::initialize_gaussian_ez(double r0,double z0,double w,double a){if(!(w>0))throw std::invalid_argument("Gaussian width must be positive");for(std::size_t j=0;j<c_.nz-1;++j)for(std::size_t i=0;i<c_.nr;++i){double r=i*c_.dr,z=(j+0.5)*c_.dz,d2=(r-r0)*(r-r0)+(z-z0)*(z-z0);ez_[ez_idx(i,j)]=a*std::exp(-d2/(2*w*w));}std::fill(er_.begin(),er_.end(),0.0);std::fill(hphi_.begin(),hphi_.end(),0.0);}
void CylindricalTM::add_soft_ez_source(std::size_t i,std::size_t j,double v){if(i>=c_.nr||j>=c_.nz-1)throw std::out_of_range("cylindrical source index");ez_[ez_idx(i,j)]+=v;}
void CylindricalTM::step(std::size_t n){for(std::size_t q=0;q<n;++q)step_once();}
void CylindricalTM::step_once(){
    cfd::core::parallel_for((c_.nr-1)*(c_.nz-1),[&](std::size_t n){std::size_t i=n%(c_.nr-1),j=n/(c_.nr-1);double dez=(ez_[ez_idx(i+1,j)]-ez_[ez_idx(i,j)])/c_.dr;double der=(er_[er_idx(i,j+1)]-er_[er_idx(i,j)])/c_.dz;hphi_[n]+=dt_/mu_*(dez-der);});
    cfd::core::parallel_for((c_.nr-1)*(c_.nz-2),[&](std::size_t n){std::size_t i=n%(c_.nr-1),j=1+n/(c_.nr-1);er_[er_idx(i,j)]-=dt_/eps_*(hphi_[h_idx(i,j)]-hphi_[h_idx(i,j-1)])/c_.dz;});
    cfd::core::parallel_for(c_.nz-1,[&](std::size_t j){ez_[ez_idx(0,j)]+=dt_/eps_*4.0*hphi_[h_idx(0,j)]/c_.dr;for(std::size_t i=1;i<c_.nr-1;++i){double rp=(i+0.5)*c_.dr,rm=(i-0.5)*c_.dr,r=i*c_.dr;ez_[ez_idx(i,j)]+=dt_/eps_*(rp*hphi_[h_idx(i,j)]-rm*hphi_[h_idx(i-1,j)])/(r*c_.dr);}ez_[ez_idx(c_.nr-1,j)]=0.0;});
    for(std::size_t i=0;i<c_.nr-1;++i){er_[er_idx(i,0)]=0;er_[er_idx(i,c_.nz-1)]=0;}
}
double CylindricalTM::energy()const{double e=0;for(std::size_t j=0;j<c_.nz-1;++j)for(std::size_t i=0;i<c_.nr;++i){double r=std::max(i*c_.dr,0.5*c_.dr),v=ez_[ez_idx(i,j)];e+=0.5*eps_*v*v*2*std::numbers::pi*r*c_.dr*c_.dz;}for(std::size_t j=0;j<c_.nz;++j)for(std::size_t i=0;i<c_.nr-1;++i){double r=(i+0.5)*c_.dr,v=er_[er_idx(i,j)];e+=0.5*eps_*v*v*2*std::numbers::pi*r*c_.dr*c_.dz;}for(std::size_t j=0;j<c_.nz-1;++j)for(std::size_t i=0;i<c_.nr-1;++i){double r=(i+0.5)*c_.dr,v=hphi_[h_idx(i,j)];e+=0.5*mu_*v*v*2*std::numbers::pi*r*c_.dr*c_.dz;}return e;}
double CylindricalTM::max_field()const{double m=0;for(double v:ez_)m=std::max(m,std::abs(v));for(double v:er_)m=std::max(m,std::abs(v));return m;}
} // namespace cfd::fdtd
