#include "cfd/particle/wakefield.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::particle {
namespace { constexpr double c0=299792458.0; }

double resonator_longitudinal_wake_v_per_c(double distance,double rs,double frequency,double q){
    if(distance<0.0||!(rs>=0.0)||!(frequency>0.0)||!(q>0.5))throw std::invalid_argument("invalid resonator wake controls");
    const double omega=2.0*std::numbers::pi*frequency,alpha=omega/(2.0*q),omega_d=omega*std::sqrt(1.0-1.0/(4.0*q*q)),time=distance/c0;
    return (omega*rs/q)*std::exp(-alpha*time)*(std::cos(omega_d*time)-(alpha/omega_d)*std::sin(omega_d*time));
}

std::vector<double> bunch_wake_potential_v(std::span<const double> wake,std::span<const double> line_charge,double ds){
    if(wake.empty()||line_charge.empty()||!(ds>0.0))throw std::invalid_argument("invalid bunch-wake convolution inputs");
    std::vector<double> potential(line_charge.size(),0.0);for(std::size_t i=0;i<line_charge.size();++i){double value=0.0;const std::size_t count=std::min(i+1U,wake.size());for(std::size_t lag=0;lag<count;++lag)value-=wake[lag]*line_charge[i-lag]*ds;potential[i]=value;}return potential;
}

Complex wake_impedance_ohm(std::span<const double> wake,double ds,double frequency){
    if(wake.empty()||!(ds>0.0)||!(frequency>=0.0))throw std::invalid_argument("invalid wake-impedance controls");const double omega=2.0*std::numbers::pi*frequency,dt=ds/c0;Complex z{};for(std::size_t i=0;i<wake.size();++i){const double time=static_cast<double>(i)*dt,weight=(i==0U||i+1U==wake.size())?0.5:1.0;z+=weight*wake[i]*Complex{std::cos(-omega*time),std::sin(-omega*time)}*dt;}return z;
}

} // namespace cfd::particle
