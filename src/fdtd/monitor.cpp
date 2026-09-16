#include "cfd/solvers/fdtd/monitor.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::fdtd {

void TimeProbe::sample(double time,double value){
    if(!std::isfinite(time)||!std::isfinite(value)) throw std::invalid_argument("non-finite FDTD probe sample");
    if(!time_.empty()&&time<time_.back()) throw std::invalid_argument("FDTD probe time must be monotone");
    time_.push_back(time);value_.push_back(value);
}
void TimeProbe::clear() noexcept{time_.clear();value_.clear();}

DftMonitor::DftMonitor(double frequency_hz):frequency_(frequency_hz){
    if(!(frequency_hz>=0.0)&&std::isfinite(frequency_hz)) throw std::invalid_argument("invalid DFT frequency");
    if(!std::isfinite(frequency_hz)||frequency_hz<0.0) throw std::invalid_argument("invalid DFT frequency");
}
void DftMonitor::sample(double time,double value){
    if(!std::isfinite(time)||!std::isfinite(value)) throw std::invalid_argument("non-finite DFT sample");
    const double omega=2.0*std::numbers::pi*frequency_;
    const std::complex<double> phase=std::exp(std::complex<double>(0.0,-omega*time));
    if(have_previous_){
        const double dt=time-previous_time_;
        if(!(dt>=0.0)) throw std::invalid_argument("DFT sample time must be monotone");
        const std::complex<double> prev_phase=std::exp(std::complex<double>(0.0,-omega*previous_time_));
        integral_+=0.5*dt*(previous_value_*prev_phase+value*phase);
    }
    previous_time_=time;previous_value_=value;have_previous_=true;++samples_;
}
void DftMonitor::clear() noexcept{integral_={};previous_time_=0.0;previous_value_=0.0;have_previous_=false;samples_=0U;}
double DftMonitor::amplitude(double duration) const{
    if(!(duration>0.0)) throw std::invalid_argument("DFT duration must be positive");
    return 2.0*std::abs(integral_)/duration;
}

} // namespace cfd::fdtd
