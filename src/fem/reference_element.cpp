#include "cfd/fem/reference_element.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {
namespace {

Point3 add(Point3 a, Point3 b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Point3 scale(Point3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}

double determinant3(const std::array<std::array<double, 3>, 3>& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

std::array<std::array<double, 3>, 3> inverse3(
    const std::array<std::array<double, 3>, 3>& a,
    double det) {
    std::array<std::array<double, 3>, 3> r{};
    r[0][0] =  (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    r[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) / det;
    r[0][2] =  (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    r[1][0] = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]) / det;
    r[1][1] =  (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    r[1][2] = -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) / det;
    r[2][0] =  (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    r[2][1] = -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) / det;
    r[2][2] =  (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return r;
}

} // namespace

std::size_t element_dimension(ElementType type) noexcept {
    switch (type) {
        case ElementType::line2: return 1U;
        case ElementType::tri3:
        case ElementType::quad4: return 2U;
        case ElementType::tet4:
        case ElementType::hex8:
        case ElementType::prism6:
        case ElementType::pyramid5: return 3U;
    }
    return 0U;
}

std::size_t element_node_count(ElementType type) noexcept {
    switch (type) {
        case ElementType::line2: return 2U;
        case ElementType::tri3: return 3U;
        case ElementType::quad4: return 4U;
        case ElementType::tet4: return 4U;
        case ElementType::hex8: return 8U;
        case ElementType::prism6: return 6U;
        case ElementType::pyramid5: return 5U;
    }
    return 0U;
}

ShapeEvaluation evaluate_shape(ElementType type, ReferencePoint p) {
    ShapeEvaluation out;
    const double x = p.xi;
    const double y = p.eta;
    const double z = p.zeta;
    switch (type) {
        case ElementType::line2:
            out.value = {0.5 * (1.0 - x), 0.5 * (1.0 + x)};
            out.gradient_reference = {{-0.5, 0.0, 0.0}, {0.5, 0.0, 0.0}};
            break;
        case ElementType::tri3:
            out.value = {1.0 - x - y, x, y};
            out.gradient_reference = {{-1.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}};
            break;
        case ElementType::quad4: {
            static constexpr std::array<int, 4> sx{-1, 1, 1, -1};
            static constexpr std::array<int, 4> sy{-1, -1, 1, 1};
            out.value.resize(4U);
            out.gradient_reference.resize(4U);
            for (std::size_t i = 0; i < 4U; ++i) {
                const double ax = 1.0 + static_cast<double>(sx[i]) * x;
                const double ay = 1.0 + static_cast<double>(sy[i]) * y;
                out.value[i] = 0.25 * ax * ay;
                out.gradient_reference[i] = {
                    0.25 * static_cast<double>(sx[i]) * ay,
                    0.25 * static_cast<double>(sy[i]) * ax,
                    0.0};
            }
            break;
        }
        case ElementType::tet4:
            out.value = {1.0 - x - y - z, x, y, z};
            out.gradient_reference = {
                {-1.0, -1.0, -1.0}, {1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
            break;
        case ElementType::hex8: {
            static constexpr std::array<int, 8> sx{-1, 1, 1, -1, -1, 1, 1, -1};
            static constexpr std::array<int, 8> sy{-1, -1, 1, 1, -1, -1, 1, 1};
            static constexpr std::array<int, 8> sz{-1, -1, -1, -1, 1, 1, 1, 1};
            out.value.resize(8U);
            out.gradient_reference.resize(8U);
            for (std::size_t i = 0; i < 8U; ++i) {
                const double ax = 1.0 + static_cast<double>(sx[i]) * x;
                const double ay = 1.0 + static_cast<double>(sy[i]) * y;
                const double az = 1.0 + static_cast<double>(sz[i]) * z;
                out.value[i] = 0.125 * ax * ay * az;
                out.gradient_reference[i] = {
                    0.125 * static_cast<double>(sx[i]) * ay * az,
                    0.125 * static_cast<double>(sy[i]) * ax * az,
                    0.125 * static_cast<double>(sz[i]) * ax * ay};
            }
            break;
        }
        case ElementType::prism6: {
            const std::array<double, 3> t{1.0 - x - y, x, y};
            const std::array<Point3, 3> gt{{{-1.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}};
            out.value.resize(6U);
            out.gradient_reference.resize(6U);
            for (std::size_t i = 0; i < 3U; ++i) {
                const double lower = 0.5 * (1.0 - z);
                const double upper = 0.5 * (1.0 + z);
                out.value[i] = t[i] * lower;
                out.value[i + 3U] = t[i] * upper;
                out.gradient_reference[i] = {gt[i].x * lower, gt[i].y * lower, -0.5 * t[i]};
                out.gradient_reference[i + 3U] = {gt[i].x * upper, gt[i].y * upper, 0.5 * t[i]};
            }
            break;
        }
        case ElementType::pyramid5: {
            static constexpr std::array<int, 4> sx{-1, 1, 1, -1};
            static constexpr std::array<int, 4> sy{-1, -1, 1, 1};
            out.value.resize(5U);
            out.gradient_reference.resize(5U);
            const double one_minus_z = 1.0 - z;
            for (std::size_t i = 0; i < 4U; ++i) {
                const double ax = 1.0 + static_cast<double>(sx[i]) * x;
                const double ay = 1.0 + static_cast<double>(sy[i]) * y;
                out.value[i] = 0.25 * ax * ay * one_minus_z;
                out.gradient_reference[i] = {
                    0.25 * static_cast<double>(sx[i]) * ay * one_minus_z,
                    0.25 * static_cast<double>(sy[i]) * ax * one_minus_z,
                    -0.25 * ax * ay};
            }
            out.value[4] = z;
            out.gradient_reference[4] = {0.0, 0.0, 1.0};
            break;
        }
    }
    return out;
}

std::vector<QuadraturePoint> gaussian_quadrature(ElementType type, unsigned order) {
    if (order == 0U || order > 2U) throw std::invalid_argument("FEM quadrature currently supports order 1 or 2");
    if (order == 1U) {
        switch (type) {
            case ElementType::line2: return {{{0.0, 0.0, 0.0}, 2.0}};
            case ElementType::tri3: return {{{1.0 / 3.0, 1.0 / 3.0, 0.0}, 0.5}};
            case ElementType::quad4: return {{{0.0, 0.0, 0.0}, 4.0}};
            case ElementType::tet4: return {{{0.25, 0.25, 0.25}, 1.0 / 6.0}};
            case ElementType::hex8: return {{{0.0, 0.0, 0.0}, 8.0}};
            case ElementType::prism6: return {{{1.0 / 3.0, 1.0 / 3.0, 0.0}, 1.0}};
            case ElementType::pyramid5: return {{{0.0, 0.0, 0.25}, 4.0 / 3.0}};
        }
    }

    const double g = 1.0 / std::sqrt(3.0);
    switch (type) {
        case ElementType::line2:
            return {{{-g, 0.0, 0.0}, 1.0}, {{g, 0.0, 0.0}, 1.0}};
        case ElementType::tri3:
            return {{{1.0 / 6.0, 1.0 / 6.0, 0.0}, 1.0 / 6.0},
                    {{2.0 / 3.0, 1.0 / 6.0, 0.0}, 1.0 / 6.0},
                    {{1.0 / 6.0, 2.0 / 3.0, 0.0}, 1.0 / 6.0}};
        case ElementType::quad4: {
            std::vector<QuadraturePoint> q;
            for (double y : {-g, g}) for (double x : {-g, g}) q.push_back({{x, y, 0.0}, 1.0});
            return q;
        }
        case ElementType::tet4: {
            constexpr double a = 0.5854101966249685;
            constexpr double b = 0.1381966011250105;
            constexpr double w = 1.0 / 24.0;
            return {{{b, b, b}, w}, {{a, b, b}, w}, {{b, a, b}, w}, {{b, b, a}, w}};
        }
        case ElementType::hex8: {
            std::vector<QuadraturePoint> q;
            for (double z : {-g, g}) for (double y : {-g, g}) for (double x : {-g, g})
                q.push_back({{x, y, z}, 1.0});
            return q;
        }
        case ElementType::prism6: {
            const auto tri = gaussian_quadrature(ElementType::tri3, 2U);
            std::vector<QuadraturePoint> q;
            for (double z : {-g, g}) for (const auto& t : tri)
                q.push_back({{t.point.xi, t.point.eta, z}, t.weight});
            return q;
        }
        case ElementType::pyramid5: {
            // Product quadrature in a collapsed cube: x=(1-t)r, y=(1-t)s,
            // z=t. The transformation contributes (1-t)^2 to the weight.
            const std::array<double, 2> t{{0.5 * (1.0 - g), 0.5 * (1.0 + g)}};
            std::vector<QuadraturePoint> q;
            for (double tv : t) {
                const double scale_xy = 1.0 - tv;
                const double wt = 0.5 * scale_xy * scale_xy;
                for (double s : {-g, g}) for (double r : {-g, g})
                    q.push_back({{scale_xy * r, scale_xy * s, tv}, wt});
            }
            return q;
        }
    }
    throw std::invalid_argument("unsupported FEM element type");
}

IsoparametricEvaluation evaluate_isoparametric(ElementType type,
                                                std::span<const Point3> nodes,
                                                ReferencePoint point) {
    const std::size_t expected = element_node_count(type);
    if (nodes.size() != expected) throw std::invalid_argument("FEM isoparametric node count mismatch");
    const auto shape = evaluate_shape(type, point);
    const std::size_t dim = element_dimension(type);

    IsoparametricEvaluation out;
    out.gradient_physical.resize(expected);
    for (std::size_t i = 0; i < expected; ++i) out.point = add(out.point, scale(nodes[i], shape.value[i]));

    // J[a][b] = dx_a / dxi_b.
    for (std::size_t i = 0; i < expected; ++i) {
        const auto& g = shape.gradient_reference[i];
        out.jacobian[0][0] += nodes[i].x * g.x;
        out.jacobian[0][1] += nodes[i].x * g.y;
        out.jacobian[0][2] += nodes[i].x * g.z;
        out.jacobian[1][0] += nodes[i].y * g.x;
        out.jacobian[1][1] += nodes[i].y * g.y;
        out.jacobian[1][2] += nodes[i].y * g.z;
        out.jacobian[2][0] += nodes[i].z * g.x;
        out.jacobian[2][1] += nodes[i].z * g.y;
        out.jacobian[2][2] += nodes[i].z * g.z;
    }

    if (dim == 1U) {
        const Point3 tangent{out.jacobian[0][0], out.jacobian[1][0], out.jacobian[2][0]};
        const double length = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y + tangent.z * tangent.z);
        if (!(length > 1.0e-30)) throw std::runtime_error("singular FEM line mapping");
        out.measure = length;
        const double inv_l2 = 1.0 / (length * length);
        for (std::size_t i = 0; i < expected; ++i) {
            const double d = shape.gradient_reference[i].x * inv_l2;
            out.gradient_physical[i] = scale(tangent, d);
        }
        return out;
    }

    if (dim == 2U) {
        const double j00 = out.jacobian[0][0];
        const double j01 = out.jacobian[0][1];
        const double j10 = out.jacobian[1][0];
        const double j11 = out.jacobian[1][1];
        const double det = j00 * j11 - j01 * j10;
        if (std::abs(det) <= 1.0e-30) throw std::runtime_error("singular FEM 2-D mapping");
        out.measure = std::abs(det);
        for (std::size_t i = 0; i < expected; ++i) {
            const auto& g = shape.gradient_reference[i];
            // grad_x N = J^{-T} grad_xi N.
            out.gradient_physical[i] = {
                ( j11 * g.x - j10 * g.y) / det,
                (-j01 * g.x + j00 * g.y) / det,
                0.0};
        }
        return out;
    }

    const double det = determinant3(out.jacobian);
    if (std::abs(det) <= 1.0e-30) throw std::runtime_error("singular FEM 3-D mapping");
    out.measure = std::abs(det);
    const auto inv = inverse3(out.jacobian, det);
    for (std::size_t i = 0; i < expected; ++i) {
        const auto& g = shape.gradient_reference[i];
        out.gradient_physical[i] = {
            inv[0][0] * g.x + inv[1][0] * g.y + inv[2][0] * g.z,
            inv[0][1] * g.x + inv[1][1] * g.y + inv[2][1] * g.z,
            inv[0][2] * g.x + inv[1][2] * g.y + inv[2][2] * g.z};
    }
    return out;
}

} // namespace cfd::fem
