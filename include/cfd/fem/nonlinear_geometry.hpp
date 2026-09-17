#pragma once
namespace cfd::fem {
struct NonlinearTrussState{double stretch{},green_lagrange_strain{},second_piola_stress{},axial_force{},tangent{};};
[[nodiscard]] NonlinearTrussState green_lagrange_truss(double reference_length,double current_length,double area,double youngs_modulus);
} // namespace cfd::fem
