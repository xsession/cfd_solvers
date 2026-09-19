#include "cfd/solvers/fvm/particle_interactions.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace cfd::fvm {
namespace {

struct BucketKey {
    std::int64_t x{};
    std::int64_t y{};
    std::int64_t z{};

    friend bool operator==(const BucketKey&, const BucketKey&) = default;
};

struct BucketHash {
    std::size_t operator()(BucketKey key) const noexcept {
        const auto mix = [](std::uint64_t value) {
            value += 0x9e3779b97f4a7c15ULL;
            value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
            value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
            return value ^ (value >> 31U);
        };
        const auto ux = static_cast<std::uint64_t>(key.x);
        const auto uy = static_cast<std::uint64_t>(key.y);
        const auto uz = static_cast<std::uint64_t>(key.z);
        return static_cast<std::size_t>(mix(ux) ^ (mix(uy) << 1U) ^ (mix(uz) >> 1U));
    }
};

[[nodiscard]] double particle_volume(const Particle& particle) {
    return std::numbers::pi * particle.diameter * particle.diameter * particle.diameter / 6.0;
}

[[nodiscard]] double particle_mass(const Particle& particle) {
    return particle.density * particle_volume(particle);
}

[[nodiscard]] BucketKey bucket_for(Vec3 position, double cell_size) {
    return {static_cast<std::int64_t>(std::floor(position.x / cell_size)),
            static_cast<std::int64_t>(std::floor(position.y / cell_size)),
            static_cast<std::int64_t>(std::floor(position.z / cell_size))};
}

[[nodiscard]] Vec3 unit_or_x(Vec3 vector) {
    const double length = magnitude(vector);
    return length > 1.0e-14 ? vector / length : Vec3{1.0, 0.0, 0.0};
}

[[nodiscard]] double momentum_magnitude(const std::vector<Particle>& particles) {
    Vec3 momentum{};
    for (const auto& particle : particles)
        momentum += particle.velocity * particle_mass(particle);
    return magnitude(momentum);
}

[[nodiscard]] double kinetic_energy(const std::vector<Particle>& particles) {
    double energy = 0.0;
    for (const auto& particle : particles) {
        const double speed = magnitude(particle.velocity);
        energy += 0.5 * particle_mass(particle) * speed * speed;
    }
    return energy;
}

} // namespace

ParticleInteractionSystem::ParticleInteractionSystem(ParticleInteractionConfig config) : config_(config) {
    validate_config();
}

void ParticleInteractionSystem::validate_config() const {
    if (!(config_.collision_distance_factor > 0.0) || !(config_.restitution >= 0.0 && config_.restitution <= 1.0) ||
        !(config_.fluid_density > 0.0) || !(config_.surface_tension > 0.0) || !(config_.coalescence_weber >= 0.0) ||
        !(config_.breakup_weber > config_.coalescence_weber) || !(config_.breakup_velocity_fraction >= 0.0) ||
        config_.maximum_particles == 0U || !std::isfinite(config_.collision_distance_factor) ||
        !std::isfinite(config_.restitution) || !std::isfinite(config_.fluid_density) ||
        !std::isfinite(config_.surface_tension) || !std::isfinite(config_.coalescence_weber) ||
        !std::isfinite(config_.breakup_weber) || !std::isfinite(config_.breakup_velocity_fraction)) {
        throw std::invalid_argument("invalid particle interaction controls");
    }
}

void ParticleInteractionSystem::validate_particles(const std::vector<Particle>& particles) const {
    if (particles.size() > config_.maximum_particles) {
        throw std::invalid_argument("particle population exceeds configured maximum");
    }
    for (const auto& particle : particles) {
        if (!std::isfinite(particle.position.x) || !std::isfinite(particle.position.y) ||
            !std::isfinite(particle.position.z) || !std::isfinite(particle.velocity.x) ||
            !std::isfinite(particle.velocity.y) || !std::isfinite(particle.velocity.z) || !(particle.diameter > 0.0) ||
            !std::isfinite(particle.diameter) || !(particle.density > 0.0) || !std::isfinite(particle.density)) {
            throw std::invalid_argument("invalid particle state");
        }
    }
}

void ParticleInteractionSystem::resolve(std::vector<Particle>& particles) {
    validate_particles(particles);
    stats_ = {};
    stats_.momentum_before = momentum_magnitude(particles);
    stats_.kinetic_energy_before = kinetic_energy(particles);
    if (particles.size() < 2U) {
        stats_.momentum_after = stats_.momentum_before;
        stats_.kinetic_energy_after = stats_.kinetic_energy_before;
        return;
    }

    double largest_diameter = 0.0;
    for (const auto& particle : particles)
        largest_diameter = std::max(largest_diameter, particle.diameter);
    const double cell_size = std::max(largest_diameter * config_.collision_distance_factor, 1.0e-14);
    std::unordered_map<BucketKey, std::vector<std::size_t>, BucketHash> buckets;
    buckets.reserve(particles.size() * 2U);
    for (std::size_t i = 0; i < particles.size(); ++i)
        buckets[bucket_for(particles[i].position, cell_size)].push_back(i);

    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    pairs.reserve(particles.size() * 2U);
    for (const auto& [key, members] : buckets) {
        for (std::int64_t dz = -1; dz <= 1; ++dz) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dx = -1; dx <= 1; ++dx) {
                    const BucketKey neighbour{key.x + dx, key.y + dy, key.z + dz};
                    const auto found = buckets.find(neighbour);
                    if (found == buckets.end())
                        continue;
                    for (const std::size_t i : members) {
                        for (const std::size_t j : found->second) {
                            if (i < j)
                                pairs.emplace_back(i, j);
                        }
                    }
                }
            }
        }
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());

    std::vector<bool> claimed(particles.size(), false);
    std::vector<bool> removed(particles.size(), false);
    std::vector<bool> breakup(particles.size(), false);
    for (const auto& [i, j] : pairs) {
        ++stats_.candidate_pairs;
        if (claimed[i] || claimed[j] || removed[i] || removed[j])
            continue;
        const Vec3 delta = particles[j].position - particles[i].position;
        const double distance = magnitude(delta);
        const double collision_distance =
            0.5 * (particles[i].diameter + particles[j].diameter) * config_.collision_distance_factor;
        if (distance > collision_distance)
            continue;
        const Vec3 normal = unit_or_x(delta);
        const double approach = dot(particles[j].velocity - particles[i].velocity, normal);
        if (!(approach < 0.0))
            continue;

        ++stats_.collisions;
        const double mean_diameter = 0.5 * (particles[i].diameter + particles[j].diameter);
        const double weber = config_.fluid_density * approach * approach * mean_diameter / config_.surface_tension;
        if (config_.enable_coalescence && weber <= config_.coalescence_weber) {
            particles[i] = coalesce_particles(particles[i], particles[j]);
            claimed[i] = true;
            claimed[j] = true;
            removed[j] = true;
            ++stats_.coalescences;
            continue;
        }

        collide_particles_elastic(particles[i], particles[j], config_.restitution);
        claimed[i] = true;
        claimed[j] = true;
        if (config_.enable_breakup && weber >= config_.breakup_weber) {
            const std::size_t target = particles[i].diameter >= particles[j].diameter ? i : j;
            breakup[target] = true;
        }
    }

    std::vector<Particle> result;
    result.reserve(particles.size() + stats_.breakups);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (removed[i]) {
            ++stats_.removed_particles;
            continue;
        }
        if (!breakup[i]) {
            result.push_back(particles[i]);
            continue;
        }
        if (result.size() + 2U > config_.maximum_particles) {
            throw std::runtime_error("particle breakup exceeds configured maximum");
        }
        const Particle parent = particles[i];
        const Vec3 direction = unit_or_x(parent.velocity);
        Particle child = parent;
        child.diameter = parent.diameter * std::cbrt(0.5);
        const double offset = 0.26 * parent.diameter;
        child.position = parent.position - direction * offset;
        result.push_back(child);
        child.position = parent.position + direction * offset;
        result.push_back(child);
        ++stats_.breakups;
    }
    particles.swap(result);
    stats_.momentum_after = momentum_magnitude(particles);
    stats_.kinetic_energy_after = kinetic_energy(particles);
}

} // namespace cfd::fvm
