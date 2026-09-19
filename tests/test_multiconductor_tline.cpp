#include "cfd/solvers/fdtd/multiconductor_tline.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <array>

int main() {
    cfd::fdtd::MulticonductorTlineConfig config{2U, {1.0e-6, 2.0e-6}, {1.0e-9, 2.0e-9}, 1.0e-6};
    cfd::fdtd::MulticonductorTlineState state{{1.0, 0.5}, {0.0, 0.0}};
    const std::array<double, 2> source{{0.0, 0.0}};
    cfd::fdtd::advance_multiconductor_tline(config, state, source);
    assert(state.current[0] > 0.0 && state.current[1] > 0.0);
    for (double value : state.voltage)
        assert(std::isfinite(value));
    std::cout << "multiconductor transmission-line regression passed\n";
    return 0;
}
