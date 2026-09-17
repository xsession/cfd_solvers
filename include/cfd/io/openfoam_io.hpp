#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <filesystem>
#include <string>
#include <vector>
namespace cfd::io {
[[nodiscard]] cfd::fvm::PolyMesh read_openfoam_polymesh(const std::filesystem::path& poly_mesh_directory);
[[nodiscard]] std::vector<double> read_openfoam_vol_scalar_field(const std::filesystem::path& field_file,std::size_t expected_cells);
void write_openfoam_vol_scalar_field(const std::filesystem::path& field_file,const std::string& object_name,const std::vector<double>& values,const std::string& dimensions="[0 0 0 0 0 0 0]");
} // namespace cfd::io
