#include "cfd/solvers/fvm/dense_particle_rheology.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] std::size_t nearest_cell(const PolyMesh& mesh, Vec3 position) {
    std::size_t best = 0U;
    double best_distance = std::numeric_limits<double>::max();
    for (std::size_t cell = 0; cell < mesh.cell_count(); ++cell) {
        const Vec3 delta = mesh.cells()[cell].center - position;
        const double distance = dot(delta, delta);
        if (distance < best_distance) {
            best_distance = distance;
            best = cell;
        }
    }
    return best;
}

[[nodiscard]] double volume(const Particle& particle) {
    return std::numbers::pi * particle.diameter * particle.diameter * particle.diameter / 6.0;
}

} // namespace

DenseParticleRheology::DenseParticleRheology(PolyMesh mesh, DenseRheologyConfig config)
    : mesh_(std::move(mesh)), config_(config) {
    validate_config();
}

void DenseParticleRheology::validate_config() const {
    if (!(config_.fluid_dynamic_viscosity > 0.0) || !(config_.particle_density > 0.0) ||
        !(config_.grain_diameter > 0.0) ||
        !(config_.maximum_solid_fraction > 0.0 && config_.maximum_solid_fraction < 1.0) ||
        !(config_.static_friction >= 0.0) || !(config_.dynamic_friction >= config_.static_friction) ||
        !(config_.inertial_scale > 0.0) || !(config_.regularization > 0.0) ||
        !std::isfinite(config_.fluid_dynamic_viscosity) || !std::isfinite(config_.particle_density) ||
        !std::isfinite(config_.grain_diameter) || !std::isfinite(config_.maximum_solid_fraction) ||
        !std::isfinite(config_.static_friction) || !std::isfinite(config_.dynamic_friction) ||
        !std::isfinite(config_.inertial_scale) || !std::isfinite(config_.regularization)) {
        throw std::invalid_argument("invalid dense rheology controls");
    }
}

DenseRheologyResult DenseParticleRheology::evaluate(std::span<const Particle> particles,
                                                    std::span<const Vec3> carrier_velocity) const {
    if (carrier_velocity.size() != mesh_.cell_count()) {
        throw std::invalid_argument("dense rheology carrier field size mismatch");
    }
    DenseRheologyResult result;
    const std::size_t cells = mesh_.cell_count();
    result.solid_fraction.assign(cells, 0.0);
    result.granular_pressure.assign(cells, 0.0);
    result.effective_viscosity.assign(cells, 0.0);
    result.carrier_momentum_source.assign(cells, {});
    result.particle_momentum_source.assign(cells, {});
    std::vector<Vec3> mean_particle_velocity(cells);
    std::vector<double> weighted_volume(cells, 0.0);
    for (const auto& particle : particles) {
        if (!(particle.diameter > 0.0) || !(particle.density > 0.0) || !std::isfinite(particle.position.x) ||
            !std::isfinite(particle.position.y) || !std::isfinite(particle.position.z) ||
            !std::isfinite(particle.velocity.x) || !std::isfinite(particle.velocity.y) ||
            !std::isfinite(particle.velocity.z)) {
            throw std::invalid_argument("invalid dense particle state");
        }
        const std::size_t cell = nearest_cell(mesh_, particle.position);
        const double particle_volume = volume(particle);
        const double weighted = particle_volume / mesh_.cells()[cell].volume;
        result.solid_fraction[cell] += weighted;
        mean_particle_velocity[cell] += particle.velocity * weighted;
        weighted_volume[cell] += weighted;
    }

    for (std::size_t cell = 0; cell < cells; ++cell) {
        const double raw_fraction = result.solid_fraction[cell];
        if (raw_fraction <= 0.0)
            continue;
        const double phi = std::min(raw_fraction, config_.maximum_solid_fraction * (1.0 - 1.0e-10));
        result.solid_fraction[cell] = phi;
        const Vec3 mean_velocity =
            mean_particle_velocity[cell] / std::max(weighted_volume[cell], config_.regularization);
        const Vec3 slip = mean_velocity - carrier_velocity[cell];
        const double slip_speed = magnitude(slip);
        const double shear_rate = slip_speed / std::max(config_.grain_diameter, config_.regularization);
        const double packing_gap = std::max(config_.maximum_solid_fraction - phi, config_.regularization);
        const double pressure = config_.particle_density * config_.grain_diameter * config_.grain_diameter *
                                shear_rate * shear_rate * phi / packing_gap;
        const double inertial_number = config_.grain_diameter * shear_rate /
                                       std::sqrt(std::max(pressure / config_.particle_density, config_.regularization));
        const double friction = config_.static_friction + (config_.dynamic_friction - config_.static_friction) *
                                                              config_.inertial_scale /
                                                              (config_.inertial_scale + inertial_number);
        const double viscosity = pressure * friction / std::max(shear_rate, config_.regularization);
        const double hindered = std::pow(std::max(1.0 - phi / config_.maximum_solid_fraction, 1.0e-6), -2.0);
        const double beta = 18.0 * config_.fluid_dynamic_viscosity * phi * (1.0 + 2.65 * phi) * hindered /
                            (config_.grain_diameter * config_.grain_diameter);
        const Vec3 force_density = slip * beta;
        result.granular_pressure[cell] = pressure;
        result.effective_viscosity[cell] = viscosity;
        result.particle_momentum_source[cell] = force_density;
        result.carrier_momentum_source[cell] = force_density * -1.0;
    }
    return result;
}

} // namespace cfd::fvm
