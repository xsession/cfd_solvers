#include "cfd/multiphysics/porous_electrothermal.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    cfd::multiphysics::PorousElectroThermalFlowState state{10.0, 300.0, 101325.0, 0.0};
    const cfd::multiphysics::PorousElectroThermalFlowConfig config{};
    const auto next = cfd::multiphysics::advance_porous_electrothermal_flow(state, config, 0.1, 1000.0);
    assert(next.species_concentration_mol_m3 < state.species_concentration_mol_m3);
    assert(next.temperature_k > state.temperature_k);
    assert(next.darcy_velocity_m_per_s < 0.0 && std::isfinite(next.pressure_pa));
    std::cout << "porous electrothermal/flow regression passed\n";
    return 0;
}
