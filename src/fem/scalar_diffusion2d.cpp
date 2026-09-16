#include "cfd/solvers/fem/scalar_diffusion2d.hpp"

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fem/reference_element.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::fem {
namespace {

Point3 as_point3(Node2 p) { return {p.x, p.y, 0.0}; }
Node2 as_node2(Point3 p) { return {p.x, p.y}; }

} // namespace

ScalarDiffusion2D::ScalarDiffusion2D(Mesh2D mesh, ScalarDiffusion2DConfig config)
    : mesh_(std::move(mesh)), config_(config), solution_(mesh_.node_count(), 0.0) {
    mesh_.validate();
    if (config_.max_iterations == 0U || !(config_.relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid ScalarDiffusion2D controls");
    }
}

void ScalarDiffusion2D::set_boundary(int patch, ScalarBoundaryCondition2D condition) {
    if (!condition.value) throw std::invalid_argument("scalar FEM boundary requires a value function");
    if (condition.type == ScalarBoundaryType::robin && !condition.alpha) {
        throw std::invalid_argument("Robin boundary requires alpha function");
    }
    boundary_[patch] = std::move(condition);
}

void ScalarDiffusion2D::solve(const std::function<double(Node2)>& diffusivity,
                              const std::function<double(Node2)>& reaction,
                              const std::function<double(Node2)>& source) {
    if (!diffusivity || !reaction || !source) throw std::invalid_argument("scalar FEM coefficient callback missing");
    const std::size_t n = mesh_.node_count();
    std::vector<double> prescribed(n, std::numeric_limits<double>::quiet_NaN());

    for (const auto& edge : mesh_.boundary_edges) {
        const auto it = boundary_.find(edge.patch);
        if (it == boundary_.end() || it->second.type != ScalarBoundaryType::dirichlet) continue;
        for (const auto node : edge.node) {
            const double v = it->second.value(mesh_.nodes[node]);
            if (std::isfinite(prescribed[node]) && std::abs(prescribed[node] - v) > 1.0e-10 * std::max(1.0, std::abs(v))) {
                throw std::runtime_error("conflicting Dirichlet values at FEM boundary corner");
            }
            prescribed[node] = v;
            solution_[node] = v;
        }
    }

    std::vector<std::size_t> free_index(n, static_cast<std::size_t>(-1));
    std::size_t free_count = 0U;
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(prescribed[i])) free_index[i] = free_count++;
    if (free_count == 0U) throw std::runtime_error("scalar FEM system has no free nodes");

    cfd::core::CsrBuilder builder(free_count, free_count);
    std::vector<double> rhs(free_count, 0.0);
    const auto volume_q = gaussian_quadrature(ElementType::tri3, 2U);

    for (const auto& tri : mesh_.triangles) {
        std::array<Point3, 3> nodes{{
            as_point3(mesh_.nodes[tri.node[0]]),
            as_point3(mesh_.nodes[tri.node[1]]),
            as_point3(mesh_.nodes[tri.node[2]])}};
        double local_a[3][3]{};
        double local_b[3]{};
        for (const auto& qp : volume_q) {
            const auto shape = evaluate_shape(ElementType::tri3, qp.point);
            const auto iso = evaluate_isoparametric(ElementType::tri3, nodes, qp.point);
            const Node2 x = as_node2(iso.point);
            const double k = diffusivity(x);
            const double c = reaction(x);
            const double f = source(x);
            if (!(k > 0.0) || !(c >= 0.0) || !std::isfinite(f)) throw std::runtime_error("invalid scalar FEM coefficient");
            const double w = qp.weight * iso.measure;
            for (std::size_t i = 0; i < 3U; ++i) {
                local_b[i] += f * shape.value[i] * w;
                for (std::size_t j = 0; j < 3U; ++j) {
                    const auto& gi = iso.gradient_physical[i];
                    const auto& gj = iso.gradient_physical[j];
                    local_a[i][j] += (k * (gi.x * gj.x + gi.y * gj.y)
                                      + c * shape.value[i] * shape.value[j]) * w;
                }
            }
        }
        for (std::size_t i = 0; i < 3U; ++i) {
            const std::size_t ni = tri.node[i];
            if (std::isfinite(prescribed[ni])) continue;
            const std::size_t row = free_index[ni];
            rhs[row] += local_b[i];
            for (std::size_t j = 0; j < 3U; ++j) {
                const std::size_t nj = tri.node[j];
                if (std::isfinite(prescribed[nj])) rhs[row] -= local_a[i][j] * prescribed[nj];
                else builder.add(row, free_index[nj], local_a[i][j]);
            }
        }
    }

    // Boundary integrals use 2-point Gauss on a physical straight edge.
    const double g = 1.0 / std::sqrt(3.0);
    for (const auto& edge : mesh_.boundary_edges) {
        const auto it = boundary_.find(edge.patch);
        if (it == boundary_.end() || it->second.type == ScalarBoundaryType::dirichlet) continue;
        const auto& bc = it->second;
        const Node2 a = mesh_.nodes[edge.node[0]];
        const Node2 b = mesh_.nodes[edge.node[1]];
        const double length = std::hypot(b.x - a.x, b.y - a.y);
        for (double xi : {-g, g}) {
            const double n0 = 0.5 * (1.0 - xi);
            const double n1 = 0.5 * (1.0 + xi);
            const Node2 x{n0 * a.x + n1 * b.x, n0 * a.y + n1 * b.y};
            const double beta = bc.value(x);
            const double alpha = bc.type == ScalarBoundaryType::robin ? bc.alpha(x) : 0.0;
            if (!std::isfinite(beta) || !(alpha >= 0.0)) throw std::runtime_error("invalid scalar FEM boundary coefficient");
            const double w = 0.5 * length;
            const std::array<double, 2> nn{{n0, n1}};
            for (std::size_t i = 0; i < 2U; ++i) {
                const std::size_t ni = edge.node[i];
                if (std::isfinite(prescribed[ni])) continue;
                const std::size_t row = free_index[ni];
                rhs[row] += beta * nn[i] * w;
                if (bc.type == ScalarBoundaryType::robin) {
                    for (std::size_t j = 0; j < 2U; ++j) {
                        const std::size_t nj = edge.node[j];
                        const double aij = alpha * nn[i] * nn[j] * w;
                        if (std::isfinite(prescribed[nj])) rhs[row] -= aij * prescribed[nj];
                        else builder.add(row, free_index[nj], aij);
                    }
                }
            }
        }
    }

    const auto matrix = builder.build();
    const auto diagonal = matrix.diagonal();
    cfd::core::JacobiPreconditioner jacobi(diagonal);
    cfd::core::KrylovWorkspace workspace;
    std::vector<double> x(free_count, 0.0);
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(prescribed[i])) x[free_index[i]] = solution_[i];
    linear_result_ = cfd::core::preconditioned_conjugate_gradient(
        rhs, x,
        [&](std::span<const double> v, std::span<double> out) { matrix.multiply(v, out); },
        [&](std::span<const double> r, std::span<double> z) { jacobi(r, z); },
        workspace, config_.max_iterations, config_.relative_tolerance);
    if (!linear_result_.converged) throw std::runtime_error("ScalarDiffusion2D PCG did not converge");
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(prescribed[i])) solution_[i] = x[free_index[i]];
}

} // namespace cfd::fem
