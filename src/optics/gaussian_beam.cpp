#include "cfd/solvers/optics/gaussian_beam.hpp"
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfd::optics {
GaussianBeam::GaussianBeam(double wavelength,double waist,double distance,double index)
    :wavelength_(wavelength),index_(index) {
    if(!(wavelength>0.0)||!(waist>0.0)||!(index>0.0)||!std::isfinite(wavelength)
       ||!std::isfinite(waist)||!std::isfinite(distance)||!std::isfinite(index))
        throw std::invalid_argument("invalid Gaussian beam parameters");
    q_={distance,std::numbers::pi*index*waist*waist/wavelength};
    if(!(q_.imag()>0.0)||!std::isfinite(q_.imag()))throw std::invalid_argument("Gaussian beam scale out of range");
}
void GaussianBeam::transform(double a,double b,double c,double d,double output_index) {
    if(!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(c)||!std::isfinite(d)
       ||!(output_index>0.0)||!std::isfinite(output_index)||!(a*d-b*c>0.0))
        throw std::invalid_argument("invalid Gaussian beam ABCD transform");
    const double determinant=a*d-b*c,expected=index_/output_index;
    if(std::abs(determinant-expected)>1.0e-10*expected)
        throw std::invalid_argument("ABCD determinant must equal input/output index ratio");
    const auto next=(a*q_+b)/(c*q_+d);
    if(!(next.imag()>0.0)||!std::isfinite(next.real())||!std::isfinite(next.imag()))
        throw std::invalid_argument("singular Gaussian beam transform");
    q_=next;index_=output_index;
}
void GaussianBeam::propagate(double distance){transform(1.0,distance,0.0,1.0,index_);}
void GaussianBeam::thin_lens(double focal_length){
    if(focal_length==0.0||!std::isfinite(focal_length))throw std::invalid_argument("invalid lens focal length");
    transform(1.0,0.0,-1.0/focal_length,1.0,index_);
}
double GaussianBeam::radius() const {return std::sqrt(-wavelength_/(std::numbers::pi*index_*(1.0/q_).imag()));}
double GaussianBeam::curvature_radius() const {
    const double inverse=(1.0/q_).real();
    return inverse==0.0?std::numeric_limits<double>::infinity():1.0/inverse;
}
} // namespace cfd::optics
