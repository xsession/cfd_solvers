#include "cfd/solvers/fdtd/port.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fdtd {

WavePort1D::WavePort1D(double frequency_hz,double impedance_ohm,int propagation_direction)
    : frequency_hz_(frequency_hz),impedance_ohm_(impedance_ohm),direction_(propagation_direction),
      forward_(frequency_hz),backward_(frequency_hz) {
    if(!std::isfinite(frequency_hz)||frequency_hz<0.0||!std::isfinite(impedance_ohm)||!(impedance_ohm>0.0)||
       (propagation_direction!=1&&propagation_direction!=-1))
        throw std::invalid_argument("invalid 1-D wave-port configuration");
}

WaveAmplitudes1D WavePort1D::instantaneous(double electric_z,double magnetic_y) const noexcept{
    const double signed_zh=static_cast<double>(direction_)*impedance_ohm_*magnetic_y;
    return {0.5*(electric_z-signed_zh),0.5*(electric_z+signed_zh)};
}

void WavePort1D::sample(double time,double electric_z,double magnetic_y){
    if(!std::isfinite(electric_z)||!std::isfinite(magnetic_y)) throw std::invalid_argument("non-finite wave-port sample");
    const auto a=instantaneous(electric_z,magnetic_y);
    forward_.sample(time,a.forward);
    backward_.sample(time,a.backward);
}

void WavePort1D::clear() noexcept{forward_.clear();backward_.clear();}

SParameters1D s_parameters(const WavePort1D& reference,const WavePort1D& transmitted){
    const auto incident=reference.forward_spectrum();
    if(std::abs(incident)<=1.0e-30) throw std::runtime_error("cannot compute S-parameters with zero incident spectrum");
    return {reference.backward_spectrum()/incident,transmitted.forward_spectrum()/incident};
}

} // namespace cfd::fdtd
