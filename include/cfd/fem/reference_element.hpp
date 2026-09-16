#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fem {

enum class ElementType {
    line2,
    tri3,
    quad4,
    tet4,
    hex8,
    prism6,
    pyramid5
};

struct Point3 {
    double x{};
    double y{};
    double z{};
};

struct ReferencePoint {
    double xi{};
    double eta{};
    double zeta{};
};

struct ShapeEvaluation {
    std::vector<double> value;
    std::vector<Point3> gradient_reference;
};

struct QuadraturePoint {
    ReferencePoint point{};
    double weight{};
};

struct IsoparametricEvaluation {
    Point3 point{};
    std::array<std::array<double, 3>, 3> jacobian{};
    double measure{};
    std::vector<Point3> gradient_physical;
};

[[nodiscard]] std::size_t element_dimension(ElementType type) noexcept;
[[nodiscard]] std::size_t element_node_count(ElementType type) noexcept;
[[nodiscard]] ShapeEvaluation evaluate_shape(ElementType type, ReferencePoint point);
[[nodiscard]] std::vector<QuadraturePoint> gaussian_quadrature(ElementType type,
                                                                unsigned order = 2U);

// 1-D elements may be embedded in 3-D. 2-D mappings currently use the x-y
// plane, which is the convention used by Mesh2D. 3-D elements use the full
// Jacobian. The returned measure is |dx/dxi| for line elements, |det J| for
// 2-D x-y mappings and |det J| for 3-D mappings.
[[nodiscard]] IsoparametricEvaluation evaluate_isoparametric(
    ElementType type,
    std::span<const Point3> physical_nodes,
    ReferencePoint point);

} // namespace cfd::fem
