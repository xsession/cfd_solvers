#include "cfd/solvers/fem/thermoelastic_block.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto block = cfd::fem::assemble_thermoelastic_block(2.0, 3.0, 4.0, 5.0, 0.5, 0.1);
    assert(std::abs(block.matrix[0] - 6.0) < 1.0e-12);
    assert(std::abs(block.matrix[8] - 28.0) < 1.0e-12);
    assert(std::abs(block.matrix[2] - block.matrix[6]) < 1.0e-12);
    std::cout << "thermoelastic monolithic block regression passed\n";
    return 0;
}
