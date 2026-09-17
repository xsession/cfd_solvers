#pragma once
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {
[[nodiscard]] inline double bdf2_derivative(double current,double old,double older,double dt){if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid BDF2 timestep");return (1.5*current-2.0*old+0.5*older)/dt;}
[[nodiscard]] inline double crank_nicolson_update(double old,double rhs_old,double rhs_new,double dt){if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid CN timestep");return old+0.5*dt*(rhs_old+rhs_new);}
} // namespace cfd::fvm
