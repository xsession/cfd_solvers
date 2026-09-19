#pragma once

#include <cstddef>

namespace cfd::electrochemistry {

struct ImplicitPnpElectrodeConfig {
    double dt_s{1.0e-3};
    double diffusivity_m2_s{1.0e-9};
    double control_volume_m3{1.0e-9};
    double boundary_concentration_mol_m3{1.0};
    double exchange_current_a{1.0e-6};
    double equilibrium_potential_v{};
    double double_layer_capacitance_f{1.0e-6};
    double temperature_k{298.15};
    double electrons{1.0};
};

struct ImplicitPnpElectrodeState {
    double concentration_mol_m3{};
    double potential_v{};
    std::size_t iterations{};
    bool converged{};
};

[[nodiscard]] ImplicitPnpElectrodeState solve_implicit_pnp_electrode_cell(double old_concentration_mol_m3,
                                                                          double old_potential_v,
                                                                          const ImplicitPnpElectrodeConfig& config = {},
                                                                          std::size_t max_iterations = 40U,
                                                                          double tolerance = 1.0e-10);

} // namespace cfd::electrochemistry
