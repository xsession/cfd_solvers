#pragma once

#include "cfd/chemistry/kinetics.hpp"
#include <span>

namespace cfd::chemistry {

struct ReactorConfig {
    double initial_step{1.0e-4};
    double minimum_step{1.0e-14};
    double relative_tolerance{1.0e-4};
    double absolute_tolerance{1.0e-10};
    std::size_t max_steps{100000U};
    std::size_t newton_iterations{30U};
};
struct ReactorResult {
    double time_advanced{};
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
};

// Isothermal, constant-volume batch reactor. Backward Euler with step doubling
// controls local error; damped Newton preserves non-negative concentrations.
// The caller state changes only after the entire requested interval succeeds.
[[nodiscard]] ReactorResult integrate_isothermal(const ReactionNetwork& network,
    std::span<double> concentrations,double temperature,double duration,const ReactorConfig& config = {});

} // namespace cfd::chemistry
