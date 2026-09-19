#include "cfd/electrochemistry/implicit_pnp.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    cfd::electrochemistry::ImplicitPnpElectrodeConfig config;
    config.boundary_concentration_mol_m3 = 2.0;
    config.exchange_current_a = 1.0e-9;
    const auto result = cfd::electrochemistry::solve_implicit_pnp_electrode_cell(1.0, 0.0, config);
    assert(result.converged && result.concentration_mol_m3 > 1.0 && std::isfinite(result.potential_v));
    std::cout << "implicit PNP/electrode cell regression passed\n";
    return 0;
}
