#pragma once

#include <cstddef>

namespace cfd::electrochemistry {

struct CorrosionCell1DConfig {
    double electrolyte_length{1.0e-3};       // m
    double electrolyte_conductivity{5.0};    // S/m
    double temperature{298.15};              // K
    double metal_potential{0.0};              // V
    double bulk_electrolyte_potential{0.0};   // V
    double equilibrium_potential{0.0};        // V
    double exchange_current_density{1.0};     // A/m^2
    double electrons{2.0};
    double anodic_transfer{0.5};
    double cathodic_transfer{0.5};
    double metal_molar_mass{0.055845};        // kg/mol, iron default
    double metal_density{7874.0};             // kg/m^3, iron default
    std::size_t max_iterations{200};
    double current_tolerance{1.0e-12};
};

struct CorrosionCell1DResult {
    double current_density{};                 // A/m^2, anodic positive
    double surface_electrolyte_potential{};   // V
    double overpotential{};                   // V
    double molar_dissolution_flux{};          // mol/(m^2 s)
    double penetration_rate{};                // m/s
    std::size_t iterations{};
    bool converged{};
};

[[nodiscard]] CorrosionCell1DResult solve_corrosion_cell_1d(const CorrosionCell1DConfig& config);

} // namespace cfd::electrochemistry
