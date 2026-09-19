#include "cfd/fem/nonlinear_geometry.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto mesh = cfd::fem::make_rectangle_tri_mesh(1U, 1U);
    const auto triangle = mesh.triangles.front();
    std::array<cfd::fem::Node2, 3> current{mesh.nodes[triangle.node[0]], mesh.nodes[triangle.node[1]],
                                           mesh.nodes[triangle.node[2]]};
    const auto reference = cfd::fem::total_lagrangian_triangle(mesh, triangle, current, 100.0, 0.3);
    assert(std::abs(reference.determinant - 1.0) < 1.0e-12);
    for (double value : reference.internal_force)
        assert(std::abs(value) < 1.0e-12);
    current[1].x += 0.1;
    const auto deformed = cfd::fem::total_lagrangian_triangle(mesh, triangle, current, 100.0, 0.3);
    assert(deformed.determinant > 1.0 && std::abs(deformed.internal_force[2]) > 0.0);
    std::cout << "total-Lagrangian nonlinear triangle regression passed\n";
    return 0;
}
