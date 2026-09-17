#include "cfd/optics/birefringence.hpp"
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>
namespace cfd::optics {
double uniaxial_extraordinary_index(double no,double ne,double theta){if(!(no>0)||!(ne>0)||!std::isfinite(theta))throw std::invalid_argument("invalid uniaxial indices");double c=std::cos(theta),s=std::sin(theta);return no*ne/std::sqrt(ne*ne*c*c+no*no*s*s);}
double birefringent_retardance(double no,double neff,double d,double lambda){if(!(no>0)||!(neff>0)||d<0||!(lambda>0))throw std::invalid_argument("invalid retardance controls");return 2.0*std::numbers::pi*(neff-no)*d/lambda;}
JonesVector apply_linear_retarder(JonesVector in,double delta,double angle){const double c=std::cos(angle),s=std::sin(angle);auto fast=c*in.x+s*in.y,slow=-s*in.x+c*in.y;slow*=std::exp(std::complex<double>(0.0,delta));return {c*fast-s*slow,s*fast+c*slow};}
}
