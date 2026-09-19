#pragma once

#include "cfd/solvers/fvm/particles.hpp"

#include <cstddef>
#include <vector>

namespace cfd::fvm {

// Deterministic cell-linked collision/coalescence/breakup controls. The
// interaction engine is deliberately independent from ParticleCloud so it
// can be used between Eulerian particle-transport stages.
struct ParticleInteractionConfig {
    double collision_distance_factor{1.0};
    double restitution{0.85};
    double fluid_density{1.2};
    double surface_tension{0.072};
    double coalescence_weber{6.0};
    double breakup_weber{12.0};
    double breakup_velocity_fraction{0.15};
    bool enable_coalescence{true};
    bool enable_breakup{true};
    std::size_t maximum_particles{1'000'000U};
};

struct ParticleInteractionStats {
    std::size_t candidate_pairs{};
    std::size_t collisions{};
    std::size_t coalescences{};
    std::size_t breakups{};
    std::size_t removed_particles{};
    double momentum_before{};
    double momentum_after{};
    double kinetic_energy_before{};
    double kinetic_energy_after{};
};

class ParticleInteractionSystem {
public:
    explicit ParticleInteractionSystem(ParticleInteractionConfig config = {});

    // Resolves each overlapping, approaching pair once. A uniform cell-linked
    // broad phase keeps the common case linear in the number of particles,
    // while the deterministic pair ordering makes regressions reproducible.
    void resolve(std::vector<Particle>& particles);

    [[nodiscard]] const ParticleInteractionConfig& config() const noexcept { return config_; }
    [[nodiscard]] const ParticleInteractionStats& last_stats() const noexcept { return stats_; }

private:
    ParticleInteractionConfig config_;
    ParticleInteractionStats stats_{};

    void validate_config() const;
    void validate_particles(const std::vector<Particle>& particles) const;
};

} // namespace cfd::fvm
