#pragma once
#include "cfd/solvers/optics/polarization.hpp"
namespace cfd::optics {
[[nodiscard]] double uniaxial_extraordinary_index(double ordinary_index,double extraordinary_index,double angle_to_optic_axis);
[[nodiscard]] double birefringent_retardance(double ordinary_index,double extraordinary_effective_index,double thickness,double vacuum_wavelength);
[[nodiscard]] JonesVector apply_linear_retarder(JonesVector input,double retardance,double fast_axis_angle=0.0);
}
