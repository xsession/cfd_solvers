#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/dense_particle_rheology.hpp"
#include "cfd/solvers/fvm/particle_interactions.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t pair_count = options.quick ? 12U : 48U * options.scale;
    const std::size_t steps = options.quick ? 8U : 32U;
    Timer setup;
    std::vector<Particle> particles;
    particles.reserve(pair_count * 2U + 2U);
    for (std::size_t pair = 0; pair < pair_count; ++pair) {
        const double y = 0.02 + 0.015 * static_cast<double>(pair % 8U);
        const double x = 0.02 + 0.035 * static_cast<double>(pair / 8U);
        particles.push_back({{x - 0.004, y, 0.5}, {0.02, 0.0, 0.0}, 0.01, 1000.0});
        particles.push_back({{x + 0.004, y, 0.5}, {-0.02, 0.0, 0.0}, 0.01, 1000.0});
    }
    // A high-Weber pair demonstrates the breakup branch in the same run.
    particles.push_back({{0.72, 0.25, 0.5}, {4.5, 0.0, 0.0}, 0.01, 1000.0});
    particles.push_back({{0.724, 0.25, 0.5}, {-4.5, 0.0, 0.0}, 0.01, 1000.0});

    ParticleInteractionConfig interaction_config;
    interaction_config.maximum_particles = particles.size() + pair_count + 8U;
    ParticleInteractionSystem interactions(interaction_config);
    auto mesh = make_cartesian_hexa_mesh(12U, 5U, 1U, 1.0, 0.5, 1.0);
    DenseRheologyConfig rheology_config;
    rheology_config.grain_diameter = 0.01;
    rheology_config.particle_density = 1000.0;
    DenseParticleRheology rheology(std::move(mesh), rheology_config);
    const std::vector<Vec3> carrier_velocity(rheology.mesh().cell_count(), {0.0, 0.0, 0.0});
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    ParticleInteractionStats stats{};
    std::size_t total_collisions = 0U;
    std::size_t total_coalescences = 0U;
    std::size_t total_breakups = 0U;
    DenseRheologyResult dense;
    for (std::size_t step = 0; step < steps; ++step) {
        interactions.resolve(particles);
        stats = interactions.last_stats();
        total_collisions += stats.collisions;
        total_coalescences += stats.coalescences;
        total_breakups += stats.breakups;
        dense = rheology.evaluate(particles, carrier_velocity);
        for (auto& particle : particles)
            particle.position += particle.velocity * 2.0e-4;
    }
    const double simulation_ms = simulation.milliseconds();

    const std::filesystem::path path = options.output.empty() ? std::filesystem::path("particle_rheology_final.csv")
                                                              : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path);
    if (!file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    file << "x_m,y_m,solid_fraction,granular_pressure,effective_viscosity,carrier_source_x_n_m3\n";
    for (std::size_t cell = 0; cell < rheology.mesh().cell_count(); ++cell) {
        const auto& center = rheology.mesh().cells()[cell].center;
        file << center.x << ',' << center.y << ',' << dense.solid_fraction[cell] << ',' << dense.granular_pressure[cell]
             << ',' << dense.effective_viscosity[cell] << ',' << dense.carrier_momentum_source[cell].x << '\n';
    }
    std::cout << "particle_interactions particles=" << particles.size() << " collisions=" << total_collisions
              << " coalescences=" << total_coalescences << " breakups=" << total_breakups
              << " max_solid_fraction=" << *std::max_element(dense.solid_fraction.begin(), dense.solid_fraction.end())
              << " csv=" << path << '\n';
    emit({"phase03_fvm", "many_particle_collision_dense_rheology", "cpu", "particle_steps", particles.size(), steps,
          setup_ms, simulation_ms, static_cast<double>(particles.size() * steps),
          static_cast<double>(total_collisions + total_coalescences + total_breakups) +
              *std::max_element(dense.solid_fraction.begin(), dense.solid_fraction.end())});
    return 0;
}
