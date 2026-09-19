#pragma once
#include "cfd/fem/mesh2d.hpp"
#include <array>
namespace cfd::fem {
struct NonlinearTrussState {
    double stretch{}, green_lagrange_strain{}, second_piola_stress{}, axial_force{}, tangent{};
};
[[nodiscard]] NonlinearTrussState green_lagrange_truss(double reference_length, double current_length, double area,
                                                       double youngs_modulus);
struct NonlinearTriangleState {
    double reference_area{};
    double determinant{};
    std::array<double, 4> deformation_gradient{};  // row-major F
    std::array<double, 3> green_lagrange_strain{}; // E_xx,E_yy,E_xy
    std::array<double, 6> internal_force{};        // nodal x/y residual
};
[[nodiscard]] NonlinearTriangleState total_lagrangian_triangle(const Mesh2D& mesh, const Tri3& triangle,
                                                               const std::array<Node2, 3>& current_nodes,
                                                               double youngs_modulus, double poisson_ratio);
} // namespace cfd::fem
