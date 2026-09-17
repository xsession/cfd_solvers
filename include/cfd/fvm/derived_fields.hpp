#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <span>
#include <vector>
namespace cfd::fvm {
struct DerivedFlowFields {std::vector<Vec3> vorticity;std::vector<double> q_criterion,lambda2,enstrophy;};
[[nodiscard]] DerivedFlowFields velocity_derived_fields(const PolyMesh& mesh,std::span<const Vec3> velocity);
[[nodiscard]] double wall_y_plus(double wall_distance,double tangential_shear_stress,double density,double dynamic_viscosity);
[[nodiscard]] Vec3 pressure_force_on_patch(const PolyMesh& mesh,std::span<const double> pressure,std::size_t patch);
[[nodiscard]] double sample_nearest_cell(const PolyMesh& mesh,std::span<const double> field,Vec3 point);
} // namespace cfd::fvm
