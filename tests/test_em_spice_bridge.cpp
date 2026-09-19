#include "cfd/multiphysics/em_spice_bridge.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto result = cfd::multiphysics::iterate_em_spice_port(
        {0.0, 0.0}, [](std::complex<double> current) { return std::complex<double>{1.0, 0.0} - 0.5 * current; },
        [](std::complex<double> voltage) { return 0.25 * voltage; },
        {.max_iterations = 80U, .relaxation = 0.8, .tolerance = 1.0e-12});
    assert(result.converged && std::abs(result.port_current - std::complex<double>{2.0 / 9.0, 0.0}) < 1.0e-9);
    assert(std::isfinite(result.port_voltage.real()));
    std::cout << "bidirectional EM/SPICE port bridge regression passed\n";
    return 0;
}
