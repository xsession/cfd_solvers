#include "cfd/solvers/fvm/compressible1d.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
namespace {
double minmod(double a,double b){return a*b<=0.0?0.0:(std::abs(a)<std::abs(b)?a:b);}
EulerState1D conservative(EulerPrimitive1D q,const IdealGasEquationOfState&eos){if(!(q.density>0.0)||!(q.pressure>0.0))throw std::invalid_argument("non-positive compressible primitive state");return {q.density,q.density*q.velocity,q.pressure/(eos.gamma-1.0)+0.5*q.density*q.velocity*q.velocity};}
EulerState1D flux(EulerState1D u,const IdealGasEquationOfState&eos){const double p=eos.pressure(u.density,u.momentum,u.energy),v=u.momentum/u.density;return {u.momentum,u.momentum*v+p,(u.energy+p)*v};}
EulerState1D add(EulerState1D a,EulerState1D b,double scale=1.0){return {a.density+scale*b.density,a.momentum+scale*b.momentum,a.energy+scale*b.energy};}
EulerState1D mul(EulerState1D a,double s){return {a.density*s,a.momentum*s,a.energy*s};}
}
double IdealGasEquationOfState::pressure(double rho,double mom,double e)const{if(!(gamma>1.0)||!(gas_constant>0.0)||!(rho>0.0))throw std::invalid_argument("invalid ideal-gas state");const double p=(gamma-1.0)*(e-0.5*mom*mom/rho);if(!(p>0.0)||!std::isfinite(p))throw std::runtime_error("non-positive ideal-gas pressure");return p;}
double IdealGasEquationOfState::sound_speed(double rho,double p)const{if(!(rho>0.0)||!(p>0.0))throw std::invalid_argument("invalid sound-speed state");return std::sqrt(gamma*p/rho);}
double IdealGasEquationOfState::temperature(double rho,double p)const{if(!(rho>0.0)||!(p>0.0))throw std::invalid_argument("invalid temperature state");return p/(rho*gas_constant);}
CompressibleEuler1D::CompressibleEuler1D(Compressible1DConfig c):config_(c),state_(c.cells),next_(c.cells){if(c.cells<8||!(c.length>0.0)||!(c.cfl>0.0&&c.cfl<1.0)||!(c.eos.gamma>1.0))throw std::invalid_argument("invalid compressible 1-D controls");initialize_uniform({1.0,0.0,1.0});}
void CompressibleEuler1D::initialize(EulerPrimitive1D left,EulerPrimitive1D right,double fraction){if(!(fraction>0.0&&fraction<1.0))throw std::invalid_argument("invalid Riemann interface");const auto ul=conservative(left,config_.eos),ur=conservative(right,config_.eos);for(std::size_t i=0;i<state_.size();++i)state_[i]=(static_cast<double>(i)+0.5)/static_cast<double>(state_.size())<fraction?ul:ur;}
void CompressibleEuler1D::initialize_uniform(EulerPrimitive1D q){const auto u=conservative(q,config_.eos);std::fill(state_.begin(),state_.end(),u);}
EulerPrimitive1D CompressibleEuler1D::primitive(std::size_t i)const{const auto&u=state_.at(i);return {u.density,u.momentum/u.density,config_.eos.pressure(u.density,u.momentum,u.energy)};}
double CompressibleEuler1D::stable_timestep()const{double speed=0.0;for(std::size_t i=0;i<state_.size();++i){auto q=primitive(i);speed=std::max(speed,std::abs(q.velocity)+config_.eos.sound_speed(q.density,q.pressure));}return config_.cfl*(config_.length/static_cast<double>(state_.size()))/speed;}
double CompressibleEuler1D::step(double dt){if(dt==0.0)dt=stable_timestep();if(!(dt>0.0))throw std::invalid_argument("invalid compressible timestep");const std::size_t n=state_.size();const double dx=config_.length/static_cast<double>(n);std::vector<EulerState1D>slope(n);for(std::size_t i=0;i<n;++i){if(!config_.periodic&&(i==0||i+1==n)){slope[i]={};continue;}const std::size_t im=i?i-1:n-1,ip=i+1<n?i+1:0;const auto&a=state_[im],&b=state_[i],&c=state_[ip];slope[i]={minmod(b.density-a.density,c.density-b.density),minmod(b.momentum-a.momentum,c.momentum-b.momentum),minmod(b.energy-a.energy,c.energy-b.energy)};}
    std::vector<EulerState1D>face(n+1);for(std::size_t f=0;f<=n;++f){std::size_t li=f==0?(config_.periodic?n-1:0):f-1,ri=f==n?(config_.periodic?0:n-1):f;EulerState1D ul=add(state_[li],slope[li],0.5),ur=add(state_[ri],slope[ri],-0.5);auto ql=[&]{double p=config_.eos.pressure(ul.density,ul.momentum,ul.energy);return EulerPrimitive1D{ul.density,ul.momentum/ul.density,p};}();auto qr=[&]{double p=config_.eos.pressure(ur.density,ur.momentum,ur.energy);return EulerPrimitive1D{ur.density,ur.momentum/ur.density,p};}();const double a=std::max(std::abs(ql.velocity)+config_.eos.sound_speed(ql.density,ql.pressure),std::abs(qr.velocity)+config_.eos.sound_speed(qr.density,qr.pressure));auto fl=flux(ul,config_.eos),fr=flux(ur,config_.eos);face[f]=add(mul(add(fl,fr),0.5),add(ur,ul,-1.0),-0.5*a);}
    for(std::size_t i=0;i<n;++i){next_[i]=add(state_[i],add(face[i+1],face[i],-1.0),-dt/dx);if(!(next_[i].density>0.0))throw std::runtime_error("compressible update produced non-positive density");(void)config_.eos.pressure(next_[i].density,next_[i].momentum,next_[i].energy);}state_.swap(next_);return dt;}
void CompressibleEuler1D::run(double duration){if(!(duration>=0.0))throw std::invalid_argument("invalid compressible duration");double t=0.0;while(t<duration){double dt=std::min(stable_timestep(),duration-t);t+=step(dt);}}
double CompressibleEuler1D::mass()const{double s=0;const double dx=config_.length/static_cast<double>(state_.size());for(auto u:state_)s+=u.density*dx;return s;}
double CompressibleEuler1D::total_energy()const{double s=0;const double dx=config_.length/static_cast<double>(state_.size());for(auto u:state_)s+=u.energy*dx;return s;}
} // namespace cfd::fvm
