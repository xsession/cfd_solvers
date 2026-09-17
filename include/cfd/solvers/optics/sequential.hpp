#pragma once
#include "cfd/solvers/optics/ray.hpp"
#include <array>
#include <cstddef>
#include <vector>
namespace cfd::optics {
enum class SurfaceType { plane, sphere, conic, even_asphere };
struct SequentialSurface {
    SurfaceType type{SurfaceType::plane}; double vertex_z{}; double radius{};
    double aperture_radius{1.0e30}; double refractive_index_after{1.0}; double conic_constant{};
    std::array<double,4> even_coefficients{};
};
[[nodiscard]] double surface_sag(const SequentialSurface& surface,double radial_distance);
struct RayTraceResult { Ray ray{}; bool valid{}; std::size_t surfaces_traced{}; };
class SequentialOpticalSystem {
public:
    explicit SequentialOpticalSystem(double object_space_index = 1.0);
    void add_surface(SequentialSurface surface);
    [[nodiscard]] RayTraceResult trace(Ray ray) const;
    [[nodiscard]] std::vector<RayTraceResult> trace_many(const std::vector<Ray>& rays) const;
    [[nodiscard]] const std::vector<SequentialSurface>& surfaces() const noexcept { return surfaces_; }
    [[nodiscard]] double object_space_index() const noexcept { return object_space_index_; }
private:
    double object_space_index_{}; std::vector<SequentialSurface> surfaces_;
};
[[nodiscard]] double rms_spot_radius_at_plane(const std::vector<RayTraceResult>& rays,double image_z);
} // namespace cfd::optics
