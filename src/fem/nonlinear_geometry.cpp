#include "cfd/fem/nonlinear_geometry.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::fem {
NonlinearTrussState green_lagrange_truss(double L, double l, double A, double E) {
    if (!(L > 0) || !(l > 0) || !(A > 0) || !(E > 0))
        throw std::invalid_argument("invalid nonlinear truss state");
    double lambda = l / L, eps = .5 * (lambda * lambda - 1.0), S = E * eps, force = A * S * lambda,
           tangent = A * E * (1.5 * lambda * lambda - .5) / L;
    return {lambda, eps, S, force, tangent};
}
NonlinearTriangleState total_lagrangian_triangle(const Mesh2D& mesh, const Tri3& triangle,
                                                 const std::array<Node2, 3>& current, double E, double nu) {
    if (!(E > 0.0) || !(nu > -1.0 && nu < 0.5))
        throw std::invalid_argument("invalid nonlinear triangle material");
    const auto a = mesh.nodes.at(triangle.node[0]), b = mesh.nodes.at(triangle.node[1]),
               c = mesh.nodes.at(triangle.node[2]);
    const double j11 = b.x - a.x, j12 = c.x - a.x, j21 = b.y - a.y, j22 = c.y - a.y, det = j11 * j22 - j12 * j21;
    if (!(det > 0.0))
        throw std::invalid_argument("nonlinear triangle requires positive reference orientation");
    const double area = 0.5 * det, inv = 1.0 / det;
    const std::array<std::array<double, 2>, 3> grad{{{{(b.y - c.y) * inv, (c.x - b.x) * inv}},
                                                     {{(c.y - a.y) * inv, (a.x - c.x) * inv}},
                                                     {{(a.y - b.y) * inv, (b.x - a.x) * inv}}}};
    const double d11 = current[1].x - current[0].x, d12 = current[2].x - current[0].x,
                 d21 = current[1].y - current[0].y, d22 = current[2].y - current[0].y;
    const double i11 = j22 * inv, i12 = -j12 * inv, i21 = -j21 * inv, i22 = j11 * inv;
    const double f11 = d11 * i11 + d12 * i21, f12 = d11 * i12 + d12 * i22, f21 = d21 * i11 + d22 * i21,
                 f22 = d21 * i12 + d22 * i22;
    const double exx = 0.5 * (f11 * f11 + f21 * f21 - 1.0), eyy = 0.5 * (f12 * f12 + f22 * f22 - 1.0),
                 exy = 0.5 * (f11 * f12 + f21 * f22);
    const double mu = E / (2.0 * (1.0 + nu)), lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
                 sxx = 2.0 * mu * exx + lambda * (exx + eyy), syy = 2.0 * mu * eyy + lambda * (exx + eyy),
                 sxy = 2.0 * mu * exy;
    const double p11 = f11 * sxx + f12 * sxy, p12 = f11 * sxy + f12 * syy, p21 = f21 * sxx + f22 * sxy,
                 p22 = f21 * sxy + f22 * syy;
    NonlinearTriangleState result;
    result.reference_area = area;
    result.determinant = f11 * f22 - f12 * f21;
    result.deformation_gradient = {f11, f12, f21, f22};
    result.green_lagrange_strain = {exx, eyy, exy};
    for (std::size_t i = 0; i < 3U; ++i) {
        result.internal_force[2U * i] = area * (p11 * grad[i][0] + p12 * grad[i][1]);
        result.internal_force[2U * i + 1U] = area * (p21 * grad[i][0] + p22 * grad[i][1]);
    }
    return result;
}
} // namespace cfd::fem
