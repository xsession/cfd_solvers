#pragma once

namespace cfd::multiphysics {

struct PorousElectroThermalFlowConfig {
    double porosity{0.4};
    double permeability_m2{1.0e-12};
    double viscosity_pa_s{1.0e-3};
    double reaction_rate_per_s{1.0};
    double reaction_heat_j_per_mol{1.0e4};
    double volumetric_heat_capacity_j_per_m3_k{1.0e6};
};

struct PorousElectroThermalFlowState {
    double species_concentration_mol_m3{};
    double temperature_k{300.0};
    double pressure_pa{};
    double darcy_velocity_m_per_s{};
};

// Single-control-volume porous electrochemistry/thermal/flow baseline. It is
// intentionally explicit and conservative so it can serve as a didactic
// coupling kernel before a full battery/electrolyzer mesh is selected.
[[nodiscard]] PorousElectroThermalFlowState
advance_porous_electrothermal_flow(PorousElectroThermalFlowState state, const PorousElectroThermalFlowConfig& config,
                                   double dt_s, double pressure_gradient_pa_per_m);

} // namespace cfd::multiphysics
