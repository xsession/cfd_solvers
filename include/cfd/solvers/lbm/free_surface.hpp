#pragma once

#include "cfd/solvers/lbm/particles.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cfd::lbm {

enum class FreeSurfacePhase : std::uint8_t { gas, interface_cell, fluid };

struct FreeSurfaceField3D {
    std::size_t nx{}, ny{}, nz{};
    std::vector<double> fill;
    std::vector<FreeSurfacePhase> phase;
};

[[nodiscard]] FreeSurfaceField3D make_free_surface_field(
    std::size_t nx, std::size_t ny, std::size_t nz,
    std::span<const double> fill,
    double interface_epsilon = 1.0e-8);

// First-order conservative periodic VOF advection. Velocity is cell-centered in
// lattice units; CFL must be chosen so no more than one cell is crossed per step.
void advect_free_surface_vof(
    FreeSurfaceField3D& surface,
    std::span<const float> ux, std::span<const float> uy, std::span<const float> uz,
    double dt,
    double interface_epsilon = 1.0e-8);

[[nodiscard]] std::vector<ParticleVector3> free_surface_normals(const FreeSurfaceField3D& surface);
[[nodiscard]] std::vector<double> free_surface_curvature(const FreeSurfaceField3D& surface);

// Continuum-surface-force baseline: a = sigma * kappa * grad(phi) / rho.
[[nodiscard]] std::vector<ParticleVector3> free_surface_surface_tension_acceleration(
    const FreeSurfaceField3D& surface,
    double surface_tension,
    double density = 1.0);

[[nodiscard]] double free_surface_volume(const FreeSurfaceField3D& surface) noexcept;

} // namespace cfd::lbm
