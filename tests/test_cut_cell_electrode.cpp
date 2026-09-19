#include "cfd/electrochemistry/corrosion_models.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    using cfd::electrochemistry::cut_cell_electrode_interface;
    const auto cut = cut_cell_electrode_interface(-0.25, 0.75, 1.0);
    assert(std::abs(cut.electrode_fraction - 0.25) < 1.0e-12);
    assert(std::abs(cut.electrolyte_fraction - 0.75) < 1.0e-12);
    assert(cut.interface_area_fraction == 1.0);
    const auto fluid = cut_cell_electrode_interface(0.1, 0.9, 1.0);
    assert(fluid.electrolyte_fraction == 1.0 && fluid.interface_area_fraction == 0.0);
    const auto solid = cut_cell_electrode_interface(-0.1, -0.9, 1.0);
    assert(solid.electrode_fraction == 1.0 && solid.interface_area_fraction == 0.0);
    std::cout << "cut-cell electrode interface regression passed\n";
    return 0;
}
