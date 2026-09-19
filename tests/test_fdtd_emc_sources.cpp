#include "cfd/solvers/fdtd/emc_sources.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto axis = cfd::fdtd::make_graded_axis(1.0, 4U, 0.1, 1.5);
    assert(axis.size() == 5U && std::abs(axis.front()) < 1.0e-15 && std::abs(axis.back() - 1.0) < 1.0e-15);
    const auto weights = cfd::fdtd::deposit_thin_wire_subcell(0.45, 0.45, 0.45, 4U, 4U, 4U, 0.25, 0.25, 0.25);
    double sum = 0.0;
    for (const auto& entry : weights)
        sum += entry.weight;
    assert(std::abs(sum - 1.0) < 1.0e-12);
    const auto source = cfd::fdtd::huygens_equivalent_source({0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    assert(std::abs(source.electric_surface_current[0] + 1.0) < 1.0e-12);
    assert(std::abs(source.magnetic_surface_current[0] + 1.0) < 1.0e-12);
    assert(std::abs(cfd::fdtd::impedance_sheet_transmission({0.0, 0.0}, {50.0, 0.0}, {50.0, 0.0}) - 1.0) < 1.0e-12);
    const std::array<cfd::fdtd::EmcTransferSample, 2> transfer{{{1.0, {1.0, 0.0}}, {3.0, {3.0, 0.0}}}};
    assert(std::abs(cfd::fdtd::interpolate_emc_transfer(transfer, 2.0) - std::complex<double>{2.0, 0.0}) < 1.0e-12);
    std::cout << "FDTD EMC source/stencil regression passed\n";
    return 0;
}
