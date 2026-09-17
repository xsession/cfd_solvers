#pragma once
#include <cmath>
#include <stdexcept>
namespace cfd::fvm {

enum class TemporalScheme { euler, backward_bdf2, crank_nicolson };

[[nodiscard]] inline double bdf2_derivative(double current,double old,double older,double dt){if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid BDF2 timestep");return (1.5*current-2.0*old+0.5*older)/dt;}
[[nodiscard]] inline double crank_nicolson_update(double old,double rhs_old,double rhs_new,double dt){if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid CN timestep");return old+0.5*dt*(rhs_old+rhs_new);}

// OpenFOAM-style off-centering: 0 -> Euler, 1 -> centered Crank-Nicolson.
// Values between zero and one retain some backward-Euler damping.
[[nodiscard]] inline double crank_nicolson_implicit_weight(double off_centering){
    if(!(off_centering>=0.0&&off_centering<=1.0)||!std::isfinite(off_centering))
        throw std::invalid_argument("Crank-Nicolson off-centering must be in [0,1]");
    return 1.0/(1.0+off_centering);
}
} // namespace cfd::fvm
