#include "cfd/solvers/fvm/premixed_flame1d.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
PremixedFlame1D::PremixedFlame1D(PremixedFlame1DConfig c):cfg_(c),c_(c.cells),next_(c.cells){double dx=c.length/static_cast<double>(c.cells-1);if(c.cells<4||!(c.length>0)||!(c.diffusivity>0)||!(c.reaction_rate>0)||!(c.dt>0)||c.dt>0.45*dx*dx/c.diffusivity)throw std::invalid_argument("invalid flame controls");}
void PremixedFlame1D::initialize_front(double x0,double w){if(!(w>0))throw std::invalid_argument("invalid flame thickness");double dx=cfg_.length/static_cast<double>(cfg_.cells-1);for(std::size_t i=0;i<cfg_.cells;++i){double x=static_cast<double>(i)*dx;c_[i]=1.0/(1.0+std::exp((x-x0)/w));}}
void PremixedFlame1D::step(std::size_t n){double dx=cfg_.length/static_cast<double>(cfg_.cells-1),a=cfg_.diffusivity/(dx*dx);for(std::size_t q=0;q<n;++q){next_[0]=1.0;next_.back()=0.0;for(std::size_t i=1;i+1<cfg_.cells;++i){double diffusion=a*(c_[i-1]-2*c_[i]+c_[i+1]);double reaction=cfg_.reaction_rate*c_[i]*(1-c_[i]);next_[i]=std::clamp(c_[i]+cfg_.dt*(diffusion+reaction),0.0,1.0);}c_.swap(next_);}}
double PremixedFlame1D::front_location()const{double dx=cfg_.length/static_cast<double>(cfg_.cells-1);for(std::size_t i=1;i<c_.size();++i)if(c_[i]<=0.5&&c_[i-1]>0.5){double f=(c_[i-1]-0.5)/(c_[i-1]-c_[i]);return (static_cast<double>(i-1)+f)*dx;}return c_.back()>0.5?cfg_.length:0.0;}
double PremixedFlame1D::theoretical_speed()const{return 2.0*std::sqrt(cfg_.diffusivity*cfg_.reaction_rate);}
double PremixedFlame1D::temperature(std::size_t i)const{if(i>=c_.size())throw std::out_of_range("flame cell");return cfg_.unburned_temperature+(cfg_.burned_temperature-cfg_.unburned_temperature)*c_[i];}
} // namespace cfd::fvm
