#include "cfd/solvers/fvm/dense_particle_rheology.hpp"
#include "cfd/solvers/fvm/particle_interactions.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void test_many_particle_coalescence() {
    using namespace cfd::fvm;
    std::vector<Particle> particles{{{0.0, 0.0, 0.0}, {0.1, 0.0, 0.0}, 0.01, 1000.0},
                                    {{0.009, 0.0, 0.0}, {-0.1, 0.0, 0.0}, 0.01, 1000.0},
                                    {{0.5, 0.5, 0.0}, {}, 0.01, 1000.0}};
    ParticleInteractionSystem system;
    system.resolve(particles);
    require(system.last_stats().candidate_pairs >= 1U, "linked-cell broad phase missed a pair");
    require(system.last_stats().coalescences == 1U, "low-Weber pair must coalesce");
    require(particles.size() == 2U, "coalescence must remove one parcel");
    const double expected_diameter = 0.01 * std::cbrt(2.0);
    require(std::abs(particles.front().diameter - expected_diameter) < 1.0e-12,
            "coalescence must conserve liquid volume");
}

void test_high_weber_breakup() {
    using namespace cfd::fvm;
    ParticleInteractionConfig config;
    config.enable_coalescence = false;
    config.maximum_particles = 8U;
    std::vector<Particle> particles{{{0.0, 0.0, 0.0}, {4.5, 0.0, 0.0}, 0.01, 1000.0},
                                    {{0.009, 0.0, 0.0}, {-4.5, 0.0, 0.0}, 0.01, 1000.0}};
    ParticleInteractionSystem system(config);
    system.resolve(particles);
    require(system.last_stats().breakups == 1U, "high-Weber pair must break up one parcel");
    require(particles.size() == 3U, "breakup must add one child parcel");
}

void test_dense_rheology_equal_opposite_sources() {
    using namespace cfd::fvm;
    DenseRheologyConfig config;
    config.grain_diameter = 0.2;
    DenseParticleRheology model(make_cartesian_hexa_mesh(2U, 1U, 1U), config);
    const std::vector<Particle> particles{{{0.25, 0.5, 0.5}, {0.4, 0.0, 0.0}, 0.2, 1000.0}};
    const std::vector<Vec3> carrier(model.mesh().cell_count(), Vec3{});
    const auto result = model.evaluate(particles, carrier);
    require(result.solid_fraction[0] > 0.0, "dense model must accumulate particle volume");
    require(result.granular_pressure[0] > 0.0, "slip must generate granular pressure");
    require(result.effective_viscosity[0] > 0.0, "slip must generate effective viscosity");
    require(std::abs(result.carrier_momentum_source[0].x + result.particle_momentum_source[0].x) < 1.0e-12,
            "dense carrier and particle sources must be equal and opposite");
}

} // namespace

int main() {
    try {
        test_many_particle_coalescence();
        test_high_weber_breakup();
        test_dense_rheology_equal_opposite_sources();
        std::cout << "particle interaction and dense rheology regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "particle interaction regression failed: " << error.what() << '\n';
        return 1;
    }
}
