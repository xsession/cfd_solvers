#include "cfd/rf/surface_mom.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <stdexcept>

namespace cfd::rf {
namespace {

using Complex = std::complex<double>;
using Point = fem::Point3;

constexpr double eps0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;
constexpr std::size_t no_index = std::numeric_limits<std::size_t>::max();

Point add(Point a, Point b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
Point sub(Point a, Point b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Point scale(Point a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double dot(Point a, Point b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Point cross(Point a, Point b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm(Point a) {
    return std::sqrt(dot(a, a));
}
bool finite(Point a) {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}
Point normalize(Point a) {
    const double length = norm(a);
    if (!(length > 0.0) || !std::isfinite(length))
        throw std::invalid_argument("surface MoM vector must be finite and nonzero");
    return scale(a, 1.0 / length);
}

struct ComplexPoint {
    Complex x{};
    Complex y{};
    Complex z{};
};

ComplexPoint add(ComplexPoint a, ComplexPoint b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
ComplexPoint scale(ComplexPoint a, Complex s) {
    return {a.x * s, a.y * s, a.z * s};
}
ComplexPoint cross(Point a, ComplexPoint b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
ComplexPoint cross(ComplexPoint a, Point b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
double norm2(ComplexPoint a) {
    return std::norm(a.x) + std::norm(a.y) + std::norm(a.z);
}
Complex dot(Point a, ComplexPoint b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct TriangleGeometry {
    Point a{};
    Point b{};
    Point c{};
    Point normal{};
    Point centroid{};
    double area{};
};

TriangleGeometry geometry(const SurfaceTriangle& triangle) {
    if (!finite(triangle.a) || !finite(triangle.b) || !finite(triangle.c))
        throw std::invalid_argument("surface MoM triangle contains a non-finite vertex");
    const Point ab = sub(triangle.b, triangle.a);
    const Point ac = sub(triangle.c, triangle.a);
    const Point twice_normal = cross(ab, ac);
    const double twice_area = norm(twice_normal);
    if (!(twice_area > 1.0e-24) || !std::isfinite(twice_area))
        throw std::invalid_argument("surface MoM triangle is degenerate");
    return {triangle.a,
            triangle.b,
            triangle.c,
            scale(twice_normal, 1.0 / twice_area),
            scale(add(add(triangle.a, triangle.b), triangle.c), 1.0 / 3.0),
            0.5 * twice_area};
}

struct Sample {
    Point position{};
    double weight{};
};

std::array<Sample, 3> samples(const TriangleGeometry& triangle) {
    // Symmetric degree-two rule. The finite point rule is intentional: the
    // singular self term is regularized explicitly by the public configuration.
    return {{{add(add(scale(triangle.a, 2.0 / 3.0), scale(triangle.b, 1.0 / 6.0)), scale(triangle.c, 1.0 / 6.0)),
              triangle.area / 3.0},
             {add(add(scale(triangle.a, 1.0 / 6.0), scale(triangle.b, 2.0 / 3.0)), scale(triangle.c, 1.0 / 6.0)),
              triangle.area / 3.0},
             {add(add(scale(triangle.a, 1.0 / 6.0), scale(triangle.b, 1.0 / 6.0)), scale(triangle.c, 2.0 / 3.0)),
              triangle.area / 3.0}}};
}

struct EdgeUse {
    std::size_t triangle{};
    std::size_t free_vertex{};
    Point free_point{};
    int sign{};
};

struct EdgeBuild {
    std::size_t start{};
    std::size_t end{};
    std::vector<EdgeUse> uses;
};

struct BasisBuild {
    RwgBasisFunction public_basis{};
    EdgeUse plus{};
    EdgeUse minus{};
    bool has_minus{};
};

struct Support {
    std::size_t triangle{};
    std::size_t free_vertex{};
    int sign{};
};

std::size_t find_or_add_vertex(std::vector<Point>& vertices, Point point, double tolerance) {
    const double tolerance2 = tolerance * tolerance;
    for (std::size_t index = 0; index < vertices.size(); ++index)
        if (dot(sub(vertices[index], point), sub(vertices[index], point)) <= tolerance2)
            return index;
    vertices.push_back(point);
    return vertices.size() - 1U;
}

ComplexMatrix multiply_vector(const ComplexMatrix& matrix, std::span<const Complex> vector) {
    if (matrix.size() != vector.size())
        throw std::invalid_argument("surface MoM matrix/vector size mismatch");
    ComplexMatrix result(matrix.size(), Complex{});
    for (std::size_t row = 0; row < matrix.size(); ++row)
        for (std::size_t column = 0; column < matrix.size(); ++column)
            result(row, 0U) += matrix(row, column) * vector[column];
    return result;
}

std::vector<Complex> solve_vector(const ComplexMatrix& matrix, std::span<const Complex> rhs) {
    const auto inverse_matrix = inverse(matrix);
    std::vector<Complex> result(matrix.size(), Complex{});
    for (std::size_t row = 0; row < matrix.size(); ++row)
        for (std::size_t column = 0; column < matrix.size(); ++column)
            result[row] += inverse_matrix(row, column) * rhs[column];
    return result;
}

Point basis_value(const BasisBuild& basis, const TriangleGeometry& triangle, std::size_t triangle_index, Point position,
                  double edge_length) {
    const EdgeUse* use = nullptr;
    if (basis.plus.triangle == triangle_index)
        use = &basis.plus;
    else if (basis.has_minus && basis.minus.triangle == triangle_index)
        use = &basis.minus;
    if (use == nullptr)
        return {};
    return scale(sub(position, use->free_point), static_cast<double>(use->sign) * edge_length / (2.0 * triangle.area));
}

double basis_divergence(const BasisBuild& basis, const TriangleGeometry& triangle, std::size_t triangle_index,
                        double edge_length) {
    if (basis.plus.triangle == triangle_index)
        return static_cast<double>(basis.plus.sign) * edge_length / triangle.area;
    if (basis.has_minus && basis.minus.triangle == triangle_index)
        return static_cast<double>(basis.minus.sign) * edge_length / triangle.area;
    return 0.0;
}

} // namespace

SurfaceMomResult solve_pec_surface_mom(std::span<const SurfaceTriangle> input, const SurfaceMomConfig& config) {
    if (input.empty())
        throw std::invalid_argument("surface MoM requires at least one triangle");
    if (!(config.frequency_hz > 0.0) || !std::isfinite(config.frequency_hz) || !(config.relative_permittivity > 0.0) ||
        !std::isfinite(config.relative_permittivity) || !(config.relative_permeability > 0.0) ||
        !std::isfinite(config.relative_permeability) || !(config.singular_regularization_fraction > 0.0) ||
        !std::isfinite(config.singular_regularization_fraction) || !(config.mesh_tolerance_m > 0.0) ||
        !std::isfinite(config.mesh_tolerance_m))
        throw std::invalid_argument("invalid surface MoM configuration");

    const Point direction = normalize(config.plane_wave_direction);
    Point electric = config.plane_wave_electric;
    if (!finite(electric))
        throw std::invalid_argument("surface MoM plane-wave electric field is not finite");
    electric = sub(electric, scale(direction, dot(direction, electric)));
    const double electric_norm = norm(electric);
    if (!(electric_norm > 1.0e-14))
        throw std::invalid_argument("plane-wave electric field must have a transverse component");
    electric = scale(electric, 1.0 / electric_norm);
    if (config.excitation == SurfaceMomExcitation::delta_gap &&
        (config.feed_edge == no_index || std::abs(config.feed_voltage_v) == 0.0 ||
         !std::isfinite(config.feed_voltage_v.real()) || !std::isfinite(config.feed_voltage_v.imag())))
        throw std::invalid_argument("delta-gap surface MoM excitation requires a finite nonzero feed");

    std::vector<Point> vertices;
    std::vector<std::array<std::size_t, 3>> triangle_vertices;
    std::vector<TriangleGeometry> triangles;
    vertices.reserve(input.size() * 3U);
    triangle_vertices.reserve(input.size());
    triangles.reserve(input.size());
    for (const auto& triangle : input) {
        const auto a = find_or_add_vertex(vertices, triangle.a, config.mesh_tolerance_m);
        const auto b = find_or_add_vertex(vertices, triangle.b, config.mesh_tolerance_m);
        const auto c = find_or_add_vertex(vertices, triangle.c, config.mesh_tolerance_m);
        if (a == b || b == c || a == c)
            throw std::invalid_argument("surface MoM triangle has repeated vertices");
        triangle_vertices.push_back({a, b, c});
        triangles.push_back(geometry(triangle));
    }

    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge_indices;
    std::vector<EdgeBuild> edges;
    for (std::size_t triangle = 0; triangle < triangle_vertices.size(); ++triangle) {
        const auto& nodes = triangle_vertices[triangle];
        for (std::size_t local_edge = 0; local_edge < 3U; ++local_edge) {
            const auto local_start = nodes[local_edge];
            const auto local_end = nodes[(local_edge + 1U) % 3U];
            const auto key = std::minmax(local_start, local_end);
            const auto [found, inserted] = edge_indices.emplace(std::make_pair(key.first, key.second), edges.size());
            const auto edge = inserted ? (edges.emplace_back(EdgeBuild{key.first, key.second, {}}), edges.size() - 1U)
                                       : found->second;
            auto& uses = edges[edge].uses;
            if (uses.size() >= 2U)
                throw std::invalid_argument("surface MoM mesh has a non-manifold edge");
            uses.push_back({triangle, nodes[(local_edge + 2U) % 3U], vertices[nodes[(local_edge + 2U) % 3U]],
                            local_start == key.first ? 1 : -1});
        }
    }

    std::vector<BasisBuild> basis;
    basis.reserve(edges.size());
    SurfaceMomResult result;
    result.triangles.assign(input.begin(), input.end());
    for (std::size_t edge = 0; edge < edges.size(); ++edge) {
        const auto& build = edges[edge];
        const auto& plus = build.uses.front();
        const auto& minus = build.uses.size() == 2U ? build.uses.back() : EdgeUse{};
        const Point start = vertices[build.start];
        const Point end = vertices[build.end];
        const double length = norm(sub(end, start));
        RwgBasisFunction public_basis;
        public_basis.edge = edge;
        public_basis.edge_start = start;
        public_basis.edge_end = end;
        public_basis.plus_free_vertex = vertices[plus.free_vertex];
        public_basis.minus_free_vertex = build.uses.size() == 2U ? vertices[minus.free_vertex] : Point{};
        public_basis.plus_triangle = plus.triangle;
        public_basis.minus_triangle = build.uses.size() == 2U ? minus.triangle : no_index;
        public_basis.edge_length_m = length;
        public_basis.plus_sign = plus.sign;
        public_basis.minus_sign = build.uses.size() == 2U ? minus.sign : 0;
        basis.push_back({public_basis, plus, minus, build.uses.size() == 2U});
    }
    if (basis.empty())
        throw std::invalid_argument("surface MoM mesh produced no RWG basis functions");
    if (config.excitation == SurfaceMomExcitation::delta_gap && config.feed_edge >= basis.size())
        throw std::out_of_range("surface MoM feed edge out of range");

    const double omega = 2.0 * std::numbers::pi * config.frequency_hz;
    const double epsilon = eps0 * config.relative_permittivity;
    const double mu = mu0 * config.relative_permeability;
    const double wave_number = omega * std::sqrt(mu * epsilon);
    const double wave_impedance = std::sqrt(mu / epsilon);
    const Complex j{0.0, 1.0};
    const auto triangle_samples = [&triangles](std::size_t index) { return samples(triangles[index]); };
    const std::size_t n = basis.size();
    ComplexMatrix matrix(n, Complex{});
    std::vector<Complex> rhs(n, Complex{});
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) {
            Complex value{};
            const std::array<Support, 2> row_support{
                {{basis[row].plus.triangle, basis[row].plus.free_vertex, basis[row].plus.sign},
                 {basis[row].minus.triangle, basis[row].minus.free_vertex, basis[row].minus.sign}}};
            const std::array<Support, 2> column_support{
                {{basis[column].plus.triangle, basis[column].plus.free_vertex, basis[column].plus.sign},
                 {basis[column].minus.triangle, basis[column].minus.free_vertex, basis[column].minus.sign}}};
            const std::size_t row_count = basis[row].has_minus ? 2U : 1U;
            const std::size_t column_count = basis[column].has_minus ? 2U : 1U;
            for (std::size_t rs = 0; rs < row_count; ++rs) {
                const auto& row_triangle = triangles[row_support[rs].triangle];
                const auto row_quadrature = triangle_samples(row_support[rs].triangle);
                for (std::size_t cs = 0; cs < column_count; ++cs) {
                    const auto& column_triangle = triangles[column_support[cs].triangle];
                    const auto column_quadrature = triangle_samples(column_support[cs].triangle);
                    for (const auto& observation : row_quadrature) {
                        const Point test = basis_value(basis[row], row_triangle, row_support[rs].triangle,
                                                       observation.position, basis[row].public_basis.edge_length_m);
                        const double test_divergence = basis_divergence(
                            basis[row], row_triangle, row_support[rs].triangle, basis[row].public_basis.edge_length_m);
                        for (const auto& source : column_quadrature) {
                            const Point trial = basis_value(basis[column], column_triangle, column_support[cs].triangle,
                                                            source.position, basis[column].public_basis.edge_length_m);
                            const double trial_divergence =
                                basis_divergence(basis[column], column_triangle, column_support[cs].triangle,
                                                 basis[column].public_basis.edge_length_m);
                            double distance = norm(sub(observation.position, source.position));
                            const double regularization = config.singular_regularization_fraction *
                                                          std::sqrt(std::min(row_triangle.area, column_triangle.area));
                            if (distance < regularization)
                                distance = regularization;
                            const Complex green =
                                std::exp(-j * wave_number * distance) / (4.0 * std::numbers::pi * distance);
                            const Complex electric_term = -j * wave_number * dot(test, trial) +
                                                          j / wave_number * test_divergence * trial_divergence;
                            value += wave_impedance * green * electric_term * observation.weight * source.weight;
                        }
                    }
                }
            }
            matrix(row, column) = value;
        }
        if (config.excitation == SurfaceMomExcitation::plane_wave) {
            const std::array<Support, 2> supports{
                {{basis[row].plus.triangle, basis[row].plus.free_vertex, basis[row].plus.sign},
                 {basis[row].minus.triangle, basis[row].minus.free_vertex, basis[row].minus.sign}}};
            const std::size_t support_count = basis[row].has_minus ? 2U : 1U;
            for (std::size_t support = 0; support < support_count; ++support) {
                const auto& triangle = triangles[supports[support].triangle];
                for (const auto& quadrature : triangle_samples(supports[support].triangle)) {
                    const Point test = basis_value(basis[row], triangle, supports[support].triangle,
                                                   quadrature.position, basis[row].public_basis.edge_length_m);
                    const Complex phase = std::exp(-j * wave_number * dot(direction, quadrature.position));
                    rhs[row] -= dot(test, ComplexPoint{electric.x * phase, electric.y * phase, electric.z * phase}) *
                                quadrature.weight;
                }
            }
        }
    }
    if (config.excitation == SurfaceMomExcitation::delta_gap)
        rhs[config.feed_edge] = config.feed_voltage_v * basis[config.feed_edge].public_basis.edge_length_m;

    result.system_matrix = matrix;
    result.current_a = solve_vector(matrix, rhs);
    result.basis.reserve(basis.size());
    for (const auto& value : basis)
        result.basis.push_back(value.public_basis);
    result.feed_edge = config.excitation == SurfaceMomExcitation::delta_gap ? config.feed_edge : no_index;
    result.wave_number_per_m = wave_number;
    result.wave_impedance_ohm = wave_impedance;
    const auto reconstructed = multiply_vector(matrix, result.current_a);
    double residual2 = 0.0;
    double rhs2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        residual2 += std::norm(reconstructed(i, 0U) - rhs[i]);
        rhs2 += std::norm(rhs[i]);
    }
    result.relative_residual = std::sqrt(residual2) / std::max(std::sqrt(rhs2), 1.0e-30);
    if (config.excitation == SurfaceMomExcitation::delta_gap) {
        result.feed_current_a = result.current_a[config.feed_edge];
        if (std::abs(result.feed_current_a) > 1.0e-30) {
            result.feed_impedance_ohm = config.feed_voltage_v / result.feed_current_a;
            result.accepted_power_w = 0.5 * std::real(config.feed_voltage_v * std::conj(result.feed_current_a));
        }
    }
    return result;
}

SurfaceMomResult solve_pec_surface_mom(const io::TriangleSurface& surface, const SurfaceMomConfig& config) {
    if (surface.triangles.empty())
        throw std::invalid_argument("surface MoM triangle surface is empty");
    std::vector<SurfaceTriangle> triangles;
    triangles.reserve(surface.triangles.size());
    for (const auto& triangle : surface.triangles)
        triangles.push_back({triangle.a, triangle.b, triangle.c});
    return solve_pec_surface_mom(std::span<const SurfaceTriangle>(triangles), config);
}

SurfaceMomFarField surface_mom_far_field(const SurfaceMomResult& result, fem::Point3 observation_direction,
                                         double distance_m) {
    if (result.triangles.empty() || result.basis.size() != result.current_a.size() ||
        !(result.wave_number_per_m > 0.0) || !(result.wave_impedance_ohm > 0.0) || !(distance_m > 0.0) ||
        !std::isfinite(distance_m))
        throw std::invalid_argument("invalid surface MoM far-field result or observation point");
    const Point direction = normalize(observation_direction);
    std::vector<TriangleGeometry> triangles;
    triangles.reserve(result.triangles.size());
    for (const auto& triangle : result.triangles)
        triangles.push_back(geometry(triangle));
    ComplexPoint current_integral{};
    for (std::size_t index = 0; index < result.basis.size(); ++index) {
        const auto& basis = result.basis[index];
        const auto integrate = [&](std::size_t triangle_index, Point free_vertex, int sign) {
            if (triangle_index == no_index)
                return;
            const auto& triangle = triangles[triangle_index];
            for (const auto& quadrature : samples(triangle)) {
                const Point value = scale(sub(quadrature.position, free_vertex),
                                          static_cast<double>(sign) * basis.edge_length_m / (2.0 * triangle.area));
                const Complex phase =
                    std::exp(Complex{0.0, 1.0} * result.wave_number_per_m * dot(direction, quadrature.position));
                current_integral = add(current_integral, scale({result.current_a[index] * value.x * phase,
                                                                result.current_a[index] * value.y * phase,
                                                                result.current_a[index] * value.z * phase},
                                                               Complex{quadrature.weight, 0.0}));
            }
        };
        integrate(basis.plus_triangle, basis.plus_free_vertex, basis.plus_sign);
        if (basis.minus_triangle != no_index)
            integrate(basis.minus_triangle, basis.minus_free_vertex, basis.minus_sign);
    }
    const ComplexPoint transverse = cross(direction, cross(current_integral, direction));
    const Complex factor = Complex{0.0, 1.0} * result.wave_number_per_m * result.wave_impedance_ohm *
                           std::exp(Complex{0.0, -result.wave_number_per_m * distance_m}) /
                           (4.0 * std::numbers::pi * distance_m);
    const ComplexPoint field = scale(transverse, factor);
    return {direction, field.x, field.y, field.z, norm2(field) / (2.0 * result.wave_impedance_ohm)};
}

} // namespace cfd::rf
