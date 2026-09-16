#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace cfd::electrochemistry {

struct ElectrodeReaction {
    std::string name;
    double equilibrium_potential{};       // V vs common reference
    double exchange_current_density{};    // A/m2
    double area{1.0};                     // m2
    double electrons{1.0};
    double alpha_anodic{0.5};
    double alpha_cathodic{0.5};
};

struct MixedPotentialResult {
    double potential{};
    double total_current{};
    std::vector<double> reaction_current_density;
    std::vector<double> reaction_current;
    std::size_t iterations{};
    bool converged{};
};

// Solve the mixed/corrosion potential where the algebraic sum of all
// electrode reaction currents equals target_current. Each reaction uses a
// Butler-Volmer polarization law against the same metal/electrode potential.
[[nodiscard]] MixedPotentialResult solve_mixed_potential(
    const std::vector<ElectrodeReaction>& reactions,
    double temperature = 298.15,
    double target_current = 0.0,
    std::size_t max_iterations = 200U,
    double relative_current_tolerance = 1.0e-10);

} // namespace cfd::electrochemistry
