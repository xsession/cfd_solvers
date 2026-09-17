#pragma once

#include "cfd/solvers/optics/ray.hpp"
#include <cstddef>
#include <vector>

namespace cfd::optics {

struct PolynomialTerm { unsigned x_power{}, y_power{}; double coefficient{}; };

class PolynomialFreeformSurface {
public:
    PolynomialFreeformSurface(double vertex_z, std::vector<PolynomialTerm> terms,
                              double aperture_radius = 1.0e30);
    [[nodiscard]] double sag(double x, double y) const;
    [[nodiscard]] Vec3 normal(double x, double y) const;
    [[nodiscard]] bool intersect(const Ray& ray, Vec3& point, Vec3& surface_normal,
                                 std::size_t max_iterations = 40,
                                 double tolerance = 1.0e-12) const;
private:
    double vertex_z_{};
    std::vector<PolynomialTerm> terms_;
    double aperture_radius_{};
};

} // namespace cfd::optics
