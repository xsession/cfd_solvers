#include "cfd/solvers/fem/mixed_stokes.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto mesh = cfd::fem::make_rectangle_tri_mesh(1U, 1U);
    const auto local = cfd::fem::assemble_mixed_stokes_tri3(mesh, mesh.triangles.front(), 2.0);
    assert(local.area > 0.0);
    double coupling_sum = 0.0;
    for (std::size_t i = 0; i < 3U; ++i)
        for (std::size_t j = 0; j < 3U; ++j)
            for (std::size_t component = 0; component < 2U; ++component)
                coupling_sum += local.matrix[(6U + i) * 9U + 2U * j + component];
    assert(std::abs(coupling_sum) < 1.0e-12);
    for (double value : local.matrix)
        assert(std::isfinite(value));
    std::cout << "mixed velocity-pressure Stokes element regression passed\n";
    return 0;
}
