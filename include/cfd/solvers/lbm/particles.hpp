#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::lbm {

struct ParticleVector3 {
    double x{}, y{}, z{};
};

struct ImmersedBoundaryParticle {
    ParticleVector3 position{};
    ParticleVector3 velocity{};
    double mass{1.0};
    double response_time{1.0};
};

struct ParticleCouplingResult {
    std::vector<float> acceleration_x;
    std::vector<float> acceleration_y;
    std::vector<float> acceleration_z;
    ParticleVector3 total_reaction_force{};
};

// Trilinear regularized-delta interpolation on a periodic unit-spacing lattice.
[[nodiscard]] ParticleVector3 interpolate_particle_velocity(
    std::size_t nx, std::size_t ny, std::size_t nz,
    std::span<const float> ux, std::span<const float> uy, std::span<const float> uz,
    ParticleVector3 position);

// Clean-room immersed-boundary point-particle baseline. Fluid velocity is
// interpolated to each Lagrangian particle with the same trilinear kernel used
// to spread the equal-and-opposite drag force back to Eulerian lattice cells.
// Unit lattice spacing, cell volume and reference fluid density are assumed.
[[nodiscard]] ParticleCouplingResult advance_immersed_boundary_particles(
    std::size_t nx, std::size_t ny, std::size_t nz,
    std::span<const float> ux, std::span<const float> uy, std::span<const float> uz,
    std::span<ImmersedBoundaryParticle> particles,
    double dt,
    bool two_way_coupling,
    ParticleVector3 gravity = {});

} // namespace cfd::lbm
