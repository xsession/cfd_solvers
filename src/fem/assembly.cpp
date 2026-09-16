#include "cfd/fem/assembly.hpp"

#include "cfd/fem/reference_element.hpp"

#include <array>
#include <algorithm>
#include <stdexcept>

namespace cfd::fem {
namespace {

Point3 p3(Node2 p) { return {p.x, p.y, 0.0}; }

std::array<std::array<double, 3>, 3> element_stiffness(
    const Mesh2D& mesh,
    const Tri3& tri,
    const std::function<double(Node2)>& diffusivity) {
    std::array<Point3, 3> nodes{{p3(mesh.nodes[tri.node[0]]), p3(mesh.nodes[tri.node[1]]), p3(mesh.nodes[tri.node[2]])}};
    std::array<std::array<double, 3>, 3> ke{};
    for (const auto& qp : gaussian_quadrature(ElementType::tri3, 2U)) {
        const auto iso = evaluate_isoparametric(ElementType::tri3, nodes, qp.point);
        const double k = diffusivity({iso.point.x, iso.point.y});
        if (!(k > 0.0)) throw std::runtime_error("invalid FEM diffusivity during assembly");
        const double w = qp.weight * iso.measure;
        for (std::size_t i = 0; i < 3U; ++i) {
            for (std::size_t j = 0; j < 3U; ++j) {
                const auto& gi = iso.gradient_physical[i];
                const auto& gj = iso.gradient_physical[j];
                ke[i][j] += k * (gi.x * gj.x + gi.y * gj.y) * w;
            }
        }
    }
    return ke;
}

} // namespace

cfd::core::CsrMatrix assemble_tri3_laplace_matrix(
    const Mesh2D& mesh,
    const std::function<double(Node2)>& diffusivity) {
    mesh.validate();
    if (!diffusivity) throw std::invalid_argument("FEM diffusivity callback missing");
    cfd::core::CsrBuilder builder(mesh.node_count(), mesh.node_count());
    for (const auto& tri : mesh.triangles) {
        const auto ke = element_stiffness(mesh, tri, diffusivity);
        for (std::size_t i = 0; i < 3U; ++i) {
            for (std::size_t j = 0; j < 3U; ++j) builder.add(tri.node[i], tri.node[j], ke[i][j]);
        }
    }
    return builder.build();
}

void apply_tri3_laplace_matrix_free(const Mesh2D& mesh,
                                    const std::function<double(Node2)>& diffusivity,
                                    std::span<const double> x,
                                    std::span<double> y) {
    mesh.validate();
    if (x.size() != mesh.node_count() || y.size() != mesh.node_count()) {
        throw std::invalid_argument("matrix-free FEM vector size mismatch");
    }
    std::fill(y.begin(), y.end(), 0.0);
    for (const auto& tri : mesh.triangles) {
        const auto ke = element_stiffness(mesh, tri, diffusivity);
        double local[3]{};
        for (std::size_t i = 0; i < 3U; ++i) {
            for (std::size_t j = 0; j < 3U; ++j) local[i] += ke[i][j] * x[tri.node[j]];
        }
        for (std::size_t i = 0; i < 3U; ++i) y[tri.node[i]] += local[i];
    }
}

} // namespace cfd::fem
