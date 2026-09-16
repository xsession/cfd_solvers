#include "cfd/solvers/fem/nonlinear_poisson2d.hpp"

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fem/reference_element.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>

namespace cfd::fem {
namespace {

Point3 p3(Node2 p) { return {p.x, p.y, 0.0}; }

} // namespace

NonlinearPoisson2D::NonlinearPoisson2D(Mesh2D mesh, NonlinearPoisson2DConfig config)
    : mesh_(std::move(mesh)), config_(config), solution_(mesh_.node_count(), 0.0) {
    mesh_.validate();
    if (!(config_.diffusivity > 0.0) || !(config_.cubic_coefficient >= 0.0)) {
        throw std::invalid_argument("invalid nonlinear Poisson coefficients");
    }
}

void NonlinearPoisson2D::solve(const std::function<double(Node2)>& source,
                               const std::function<double(Node2)>& dirichlet_value) {
    if (!source || !dirichlet_value) throw std::invalid_argument("nonlinear Poisson callback missing");
    const std::size_t n = mesh_.node_count();
    std::vector<std::size_t> free_index(n, static_cast<std::size_t>(-1));
    std::size_t free_count = 0U;
    for (std::size_t i = 0U; i < n; ++i) {
        if (mesh_.boundary_node[i]) solution_[i] = dirichlet_value(mesh_.nodes[i]);
        else free_index[i] = free_count++;
    }
    if (free_count == 0U) throw std::runtime_error("nonlinear Poisson mesh has no free nodes");

    std::vector<double> x(free_count, 0.0);
    for (std::size_t i = 0U; i < n; ++i) if (!mesh_.boundary_node[i]) x[free_index[i]] = solution_[i];
    const auto quadrature = gaussian_quadrature(ElementType::tri3, 2U);

    auto gather_full = [&](std::span<const double> free, std::vector<double>& full) {
        full = solution_;
        for (std::size_t i = 0U; i < n; ++i) if (!mesh_.boundary_node[i]) full[i] = free[free_index[i]];
    };

    auto residual = [&](std::span<const double> free, std::span<double> out) {
        if (out.size() != free_count) throw std::runtime_error("nonlinear residual size mismatch");
        std::fill(out.begin(), out.end(), 0.0);
        std::vector<double> full;
        gather_full(free, full);
        for (const auto& tri : mesh_.triangles) {
            std::array<Point3, 3> nodes{{p3(mesh_.nodes[tri.node[0]]), p3(mesh_.nodes[tri.node[1]]), p3(mesh_.nodes[tri.node[2]])}};
            double local[3]{};
            for (const auto& qp : quadrature) {
                const auto shape = evaluate_shape(ElementType::tri3, qp.point);
                const auto iso = evaluate_isoparametric(ElementType::tri3, nodes, qp.point);
                double uq = 0.0;
                Point3 grad{};
                for (std::size_t j = 0U; j < 3U; ++j) {
                    const double u = full[tri.node[j]];
                    uq += shape.value[j] * u;
                    grad.x += iso.gradient_physical[j].x * u;
                    grad.y += iso.gradient_physical[j].y * u;
                }
                const Node2 point{iso.point.x, iso.point.y};
                const double f = source(point);
                const double w = qp.weight * iso.measure;
                for (std::size_t i = 0U; i < 3U; ++i) {
                    const auto& gi = iso.gradient_physical[i];
                    local[i] += (config_.diffusivity * (gi.x * grad.x + gi.y * grad.y)
                                 + config_.cubic_coefficient * shape.value[i] * uq * uq * uq
                                 - shape.value[i] * f) * w;
                }
            }
            for (std::size_t i = 0U; i < 3U; ++i) {
                const auto node = tri.node[i];
                if (!mesh_.boundary_node[node]) out[free_index[node]] += local[i];
            }
        }
    };

    auto jacobian = [&](std::span<const double> free) {
        std::vector<double> full;
        gather_full(free, full);
        cfd::core::CsrBuilder builder(free_count, free_count);
        for (const auto& tri : mesh_.triangles) {
            std::array<Point3, 3> nodes{{p3(mesh_.nodes[tri.node[0]]), p3(mesh_.nodes[tri.node[1]]), p3(mesh_.nodes[tri.node[2]])}};
            double local[3][3]{};
            for (const auto& qp : quadrature) {
                const auto shape = evaluate_shape(ElementType::tri3, qp.point);
                const auto iso = evaluate_isoparametric(ElementType::tri3, nodes, qp.point);
                double uq = 0.0;
                for (std::size_t j = 0U; j < 3U; ++j) uq += shape.value[j] * full[tri.node[j]];
                const double w = qp.weight * iso.measure;
                for (std::size_t i = 0U; i < 3U; ++i) {
                    for (std::size_t j = 0U; j < 3U; ++j) {
                        const auto& gi = iso.gradient_physical[i];
                        const auto& gj = iso.gradient_physical[j];
                        local[i][j] += (config_.diffusivity * (gi.x * gj.x + gi.y * gj.y)
                            + 3.0 * config_.cubic_coefficient * shape.value[i] * shape.value[j] * uq * uq) * w;
                    }
                }
            }
            for (std::size_t i = 0U; i < 3U; ++i) {
                const auto ni = tri.node[i];
                if (mesh_.boundary_node[ni]) continue;
                for (std::size_t j = 0U; j < 3U; ++j) {
                    const auto nj = tri.node[j];
                    if (!mesh_.boundary_node[nj]) builder.add(free_index[ni], free_index[nj], local[i][j]);
                }
            }
        }
        return builder.build();
    };

    nonlinear_result_ = cfd::core::newton_solve(x, residual, jacobian, config_.newton);
    if (!nonlinear_result_.converged) throw std::runtime_error("NonlinearPoisson2D Newton solve did not converge");
    for (std::size_t i = 0U; i < n; ++i) if (!mesh_.boundary_node[i]) solution_[i] = x[free_index[i]];
}

} // namespace cfd::fem
