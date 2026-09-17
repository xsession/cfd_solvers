#include "cfd/rf/antenna.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <limits>
#include <stdexcept>

namespace cfd::rf {
namespace {
constexpr double c0=299792458.0;
constexpr double eta0=376.730313668;
}

HertzianFarField hertzian_dipole_far_field(double length,double current,double frequency,double theta,double distance) {
    if (!(length>0.0) || !std::isfinite(length) || !std::isfinite(current) || !(frequency>0.0)
        || !std::isfinite(frequency) || !std::isfinite(theta) || !(distance>0.0) || !std::isfinite(distance))
        throw std::invalid_argument("invalid Hertzian-dipole controls");
    const double k=2.0*std::numbers::pi*frequency/c0;
    const std::complex<double> phase=std::exp(std::complex<double>{0.0,-k*distance});
    const auto e=std::complex<double>{0.0,1.0}*eta0*k*current*length*std::sin(theta)
        *phase/(4.0*std::numbers::pi*distance);
    return {e,e/eta0};
}

double shielding_effectiveness_db(double incident,double transmitted) {
    if (!(incident>0.0) || !(transmitted>0.0) || !std::isfinite(incident) || !std::isfinite(transmitted))
        throw std::invalid_argument("shielding field magnitudes must be finite and positive");
    return 20.0*std::log10(incident/transmitted);
}

SinusoidalDipoleSolver::SinusoidalDipoleSolver(double length,double frequency,double current)
    :length_(length),frequency_(frequency),current_(current) {
    if(!(length_>0.0)||!(frequency_>0.0)||!(current_>0.0)||!std::isfinite(length_)
       ||!std::isfinite(frequency_)||!std::isfinite(current_))
        throw std::invalid_argument("invalid dipole parameters");
}

double SinusoidalDipoleSolver::wavelength() const noexcept { return c0/frequency_; }
double SinusoidalDipoleSolver::electrical_length() const noexcept { return length_/wavelength(); }

double SinusoidalDipoleSolver::pattern_factor(double theta) const {
    const double k=2.0*std::numbers::pi/wavelength();
    const double s=std::sin(theta);
    if(std::abs(s)<1.0e-10) return 0.0;
    return (std::cos(0.5*k*length_*std::cos(theta))-std::cos(0.5*k*length_))/s;
}

std::complex<double> SinusoidalDipoleSolver::far_e_theta(double theta,double r) const {
    if(!(r>0.0)||!std::isfinite(r))throw std::invalid_argument("far-field distance must be positive");
    const double k=2.0*std::numbers::pi/wavelength();
    const std::complex<double> phase=std::exp(std::complex<double>(0.0,-k*r));
    return std::complex<double>(0.0,1.0)*eta0*current_/(2.0*std::numbers::pi*r)*pattern_factor(theta)*phase;
}

double SinusoidalDipoleSolver::radiated_power(std::size_t n) const {
    if(n<16U)throw std::invalid_argument("dipole integration requires at least 16 points");
    if(n%2U)++n;
    const double h=std::numbers::pi/static_cast<double>(n);
    double integral=0.0;
    for(std::size_t i=0;i<=n;++i){
        const double theta=h*static_cast<double>(i);
        const double f=pattern_factor(theta);
        const double value=f*f*std::sin(theta);
        const double weight=(i==0||i==n)?1.0:(i%2U?4.0:2.0);
        integral+=weight*value;
    }
    integral*=h/3.0;
    return eta0*current_*current_/(4.0*std::numbers::pi)*integral;
}

double SinusoidalDipoleSolver::radiation_resistance(std::size_t n) const {
    return 2.0*radiated_power(n)/(current_*current_);
}

double SinusoidalDipoleSolver::directivity(std::size_t n) const {
    if(n<16U)throw std::invalid_argument("dipole integration requires at least 16 points");
    const double p=radiated_power(n);
    double max_u=0.0;
    for(std::size_t i=0;i<=n;++i){
        const double theta=std::numbers::pi*static_cast<double>(i)/static_cast<double>(n);
        const double f=pattern_factor(theta);
        max_u=std::max(max_u,eta0*current_*current_/(8.0*std::numbers::pi*std::numbers::pi)*f*f);
    }
    return 4.0*std::numbers::pi*max_u/p;
}

std::vector<PatternSample> SinusoidalDipoleSolver::normalized_pattern(std::size_t samples) const {
    if(samples<2U)throw std::invalid_argument("dipole pattern requires at least two samples");
    std::vector<PatternSample> out(samples);double maximum=0.0;
    for(std::size_t i=0;i<samples;++i){
        const double theta=std::numbers::pi*static_cast<double>(i)/static_cast<double>(samples-1U);
        const double f=pattern_factor(theta);out[i]={theta,f*f};maximum=std::max(maximum,f*f);
    }
    if(maximum>0.0)for(auto& sample:out)sample.normalized_power/=maximum;
    return out;
}

std::complex<double> array_factor(std::span<const ArrayElement> elements,double frequency,double theta,double phi) {
    if(elements.empty()||!(frequency>0.0)||!std::isfinite(frequency)||!std::isfinite(theta)||!std::isfinite(phi))
        throw std::invalid_argument("invalid phased-array controls");
    const double k=2.0*std::numbers::pi*frequency/c0;
    const double sx=std::sin(theta)*std::cos(phi),sy=std::sin(theta)*std::sin(phi),sz=std::cos(theta);
    std::complex<double> result{};
    double norm=0.0;
    for(const auto& element:elements){
        if(!std::isfinite(element.x_m)||!std::isfinite(element.y_m)||!std::isfinite(element.z_m))
            throw std::invalid_argument("non-finite array-element coordinate");
        const double phase=k*(element.x_m*sx+element.y_m*sy+element.z_m*sz);
        result+=element.weight*std::exp(std::complex<double>(0.0,phase));
        norm+=std::abs(element.weight);
    }
    return norm>0.0?result/norm:std::complex<double>{};
}

PolarizationMetrics polarization_metrics(std::complex<double> et,std::complex<double> ep) {
    const double s0=std::norm(et)+std::norm(ep);
    if(!(s0>0.0)||!std::isfinite(s0)) throw std::invalid_argument("polarization field must be nonzero and finite");
    const double s1=std::norm(et)-std::norm(ep);
    const double s2=2.0*std::real(et*std::conj(ep));
    const double s3=-2.0*std::imag(et*std::conj(ep));
    const double linear=std::hypot(s1,s2);
    const double major2=0.5*(s0+linear),minor2=std::max(0.0,0.5*(s0-linear));
    const double axial=minor2>1.0e-30?std::sqrt(major2/minor2):std::numeric_limits<double>::infinity();
    const double tilt=0.5*std::atan2(s2,s1);
    const double ellipticity=0.5*std::asin(std::clamp(s3/s0,-1.0,1.0));
    return {axial,tilt,ellipticity,s3<0.0};
}

} // namespace cfd::rf
