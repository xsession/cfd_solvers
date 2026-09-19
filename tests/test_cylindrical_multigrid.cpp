#include "cfd/solvers/fdtd/cylindrical_multigrid.hpp"

#include <cassert>
#include <iostream>

int main() {
    const auto levels = cfd::fdtd::make_cylindrical_multigrid_levels(16U, 8U, 1.0, 2.0, 4U);
    assert(levels.size() == 4U && levels[0].radial_cells == 16U && levels[1].radial_cells == 8U);
    assert(levels.back().radial_spacing > levels.front().radial_spacing);
    std::cout << "cylindrical multigrid hierarchy regression passed\n";
    return 0;
}
