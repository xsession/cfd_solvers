#pragma once
#include <cstddef>
#include <filesystem>
#include <span>
#include <vector>
namespace cfd::lbm {
struct Vector2{double x{},y{};};
void write_vtk_structured_2d(const std::filesystem::path& path,std::size_t nx,std::size_t ny,std::span<const double> scalar,std::span<const Vector2> velocity={},const char* scalar_name="rho");
[[nodiscard]] std::vector<double> vorticity_2d(std::size_t nx,std::size_t ny,std::span<const Vector2> velocity,double dx=1.0,double dy=1.0);
[[nodiscard]] std::vector<double> horizontal_slice(std::size_t nx,std::size_t ny,std::span<const double> field,std::size_t y);
[[nodiscard]] std::vector<Vector2> streamline_2d(std::size_t nx,std::size_t ny,std::span<const Vector2> velocity,Vector2 start,double step,std::size_t max_points);
void write_pgm_scalar(const std::filesystem::path& path,std::size_t nx,std::size_t ny,std::span<const double> field);
} // namespace cfd::lbm
