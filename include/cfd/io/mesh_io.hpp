#pragma once
#include "cfd/fem/mesh2d.hpp"
#include "cfd/fem/mesh3d.hpp"
#include <array>
#include <filesystem>
#include <string>
#include <vector>
namespace cfd::io {
struct Triangle3 { cfd::fem::Point3 a,b,c; };
struct TriangleSurface { std::vector<Triangle3> triangles; };
[[nodiscard]] cfd::fem::Mesh2D read_gmsh_v2_tri2d(const std::filesystem::path& path);
[[nodiscard]] cfd::fem::Mesh3D read_gmsh_v2_tet3d(const std::filesystem::path& path);
[[nodiscard]] TriangleSurface read_obj_triangles(const std::filesystem::path& path);
[[nodiscard]] TriangleSurface read_ascii_stl(const std::filesystem::path& path);
[[nodiscard]] cfd::fem::Mesh2D read_vtu_tri2d(const std::filesystem::path& path);
void write_vtu(const std::filesystem::path& path,const cfd::fem::Mesh2D& mesh,const std::vector<double>& nodal_scalar={},std::string scalar_name="field");
void write_vtu(const std::filesystem::path& path,const cfd::fem::Mesh3D& mesh,const std::vector<double>& nodal_scalar={},std::string scalar_name="field");
void write_xdmf_inline(const std::filesystem::path& path,const cfd::fem::Mesh2D& mesh,const std::vector<double>& nodal_scalar={},std::string scalar_name="field");
[[nodiscard]] bool hdf5_output_available() noexcept;
void write_hdf5_xdmf(const std::filesystem::path& xdmf_path,const std::filesystem::path& hdf5_path,const cfd::fem::Mesh2D& mesh,const std::vector<double>& nodal_scalar={},std::string scalar_name="field");
} // namespace cfd::io
