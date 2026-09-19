#include "cfd/multiphysics/porous_electrothermal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {

PorousElectroThermalFlowState advance_porous_electrothermal_flow(PorousElectroThermalFlowState state,
                                                                 const PorousElectroThermalFlowConfig& config,
                                                                 double dt, double pressure_gradient) {
    if (!(config.porosity > 0.0 && config.porosity <= 1.0) || !(config.permeability_m2 >= 0.0) ||
        !(config.viscosity_pa_s > 0.0) || !(config.reaction_rate_per_s >= 0.0) ||
        !(config.reaction_heat_j_per_mol >= 0.0) || !(config.volumetric_heat_capacity_j_per_m3_k > 0.0) ||
        !(dt >= 0.0) || !(state.species_concentration_mol_m3 >= 0.0) || !(state.temperature_k > 0.0) ||
        !std::isfinite(pressure_gradient)) {
        throw std::invalid_argument("invalid porous electrothermal flow state");
    }
    const double rate = config.reaction_rate_per_s * state.species_concentration_mol_m3;
    state.species_concentration_mol_m3 = std::max(0.0, state.species_concentration_mol_m3 - dt * rate);
    state.temperature_k += dt * rate * config.reaction_heat_j_per_mol / config.volumetric_heat_capacity_j_per_m3_k;
    state.darcy_velocity_m_per_s = -config.permeability_m2 / config.viscosity_pa_s * pressure_gradient;
    state.pressure_pa += -pressure_gradient * state.darcy_velocity_m_per_s * dt;
    return state;
}

} // namespace cfd::multiphysics
