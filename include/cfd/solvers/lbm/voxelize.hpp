#pragma once
#include "cfd/io/mesh_io.hpp"
#include "cfd/fem/reference_element.hpp"
#include <array>
#include <cstddef>
#include <vector>
namespace cfd::lbm {
struct VoxelGrid3D {std::size_t nx{},ny{},nz{};cfd::fem::Point3 minimum{},maximum{};};
[[nodiscard]] std::vector<unsigned char> voxelize_surface(const cfd::io::TriangleSurface& surface,const VoxelGrid3D& grid);
[[nodiscard]] cfd::io::TriangleSurface transform_surface(const cfd::io::TriangleSurface& surface,cfd::fem::Point3 translation,cfd::fem::Point3 rotation_axis,double angle_radians);
}
