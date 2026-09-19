#include "cfd/solvers/fem/mixed_stokes.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {

MixedStokesTri3LocalSystem assemble_mixed_stokes_tri3(const Mesh2D& mesh, const Tri3& triangle, double viscosity,
                                                      double pressure_stabilization) {
    if (!(viscosity > 0.0) || !(pressure_stabilization >= 0.0))
        throw std::invalid_argument("invalid mixed Stokes controls");
    const Node2 a = mesh.nodes.at(triangle.node[0]), b = mesh.nodes.at(triangle.node[1]),
                c = mesh.nodes.at(triangle.node[2]);
    const double det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    const double area = 0.5 * std::abs(det);
    if (!(area > 0.0))
        throw std::invalid_argument("degenerate mixed Stokes triangle");
    const double inv_det = 1.0 / det;
    const std::array<std::array<double, 2>, 3> grad{{{{(b.y - c.y) * inv_det, (c.x - b.x) * inv_det}},
                                                     {{(c.y - a.y) * inv_det, (a.x - c.x) * inv_det}},
                                                     {{(a.y - b.y) * inv_det, (b.x - a.x) * inv_det}}}};
    MixedStokesTri3LocalSystem result;
    result.area = area;
    const auto add = [&](std::size_t row, std::size_t col, double value) { result.matrix[row * 9U + col] += value; };
    for (std::size_t i = 0; i < 3U; ++i)
        for (std::size_t j = 0; j < 3U; ++j) {
            const double visc = viscosity * area * (grad[i][0] * grad[j][0] + grad[i][1] * grad[j][1]);
            add(2U * i, 2U * j, visc);
            add(2U * i + 1U, 2U * j + 1U, visc);
            for (std::size_t component = 0; component < 2U; ++component) {
                const double coupling = -area / 3.0 * grad[j][component];
                add(6U + i, 2U * j + component, coupling);
                add(2U * j + component, 6U + i, coupling);
            }
            add(6U + i, 6U + j, pressure_stabilization * area * (i == j ? 2.0 : 1.0) / 12.0);
        }
    return result;
}

} // namespace cfd::fem
