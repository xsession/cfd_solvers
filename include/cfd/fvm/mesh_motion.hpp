#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <span>
#include <vector>
namespace cfd::fvm {
[[nodiscard]] PolyMesh rigidly_moved_mesh(const PolyMesh& mesh,Vec3 translation,Vec3 rotation_axis,double angle_radians);
[[nodiscard]] std::vector<double> ale_mesh_flux(const PolyMesh& old_mesh,const PolyMesh& new_mesh,double dt);
[[nodiscard]] std::vector<double> conservative_nearest_remap(const PolyMesh& old_mesh,std::span<const double> old_field,const PolyMesh& new_mesh);
}
