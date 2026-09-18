#pragma once
#include <cstdint>
#include <cstddef>
#include <filesystem>
#include <span>
#include <vector>
namespace cfd::lbm {
struct Vector2 {
    double x{}, y{};
};
struct Vector3 {
    double x{}, y{}, z{};
};
void write_vtk_structured_2d(const std::filesystem::path& path, std::size_t nx, std::size_t ny,
                             std::span<const double> scalar, std::span<const Vector2> velocity = {},
                             const char* scalar_name = "rho");
void write_vtk_structured_3d(const std::filesystem::path& path, std::size_t nx, std::size_t ny, std::size_t nz,
                             std::span<const double> scalar, std::span<const Vector3> velocity = {},
                             std::span<const std::uint8_t> solid = {}, const char* scalar_name = "field");
void write_json_state_3d(const std::filesystem::path& path, std::size_t nx, std::size_t ny, std::size_t nz,
                         std::uint64_t step, std::span<const double> scalar, std::span<const std::uint8_t> solid = {},
                         const char* scalar_name = "field");
[[nodiscard]] std::vector<double> vorticity_2d(std::size_t nx, std::size_t ny, std::span<const Vector2> velocity,
                                               double dx = 1.0, double dy = 1.0);
[[nodiscard]] std::vector<double> q_criterion_3d(std::size_t nx, std::size_t ny, std::size_t nz,
                                                 std::span<const Vector3> velocity, double dx = 1.0, double dy = 1.0,
                                                 double dz = 1.0);
[[nodiscard]] std::vector<double> horizontal_slice(std::size_t nx, std::size_t ny, std::span<const double> field,
                                                   std::size_t y);
[[nodiscard]] std::vector<Vector2> streamline_2d(std::size_t nx, std::size_t ny, std::span<const Vector2> velocity,
                                                 Vector2 start, double step, std::size_t max_points);
void write_pgm_scalar(const std::filesystem::path& path, std::size_t nx, std::size_t ny, std::span<const double> field);
void write_html_scalar_viewer(const std::filesystem::path& path, std::size_t nx, std::size_t ny,
                              std::span<const double> field);
} // namespace cfd::lbm
