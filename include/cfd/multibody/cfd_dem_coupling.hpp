#pragma once

#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/multibody/rigid_body.hpp"
#include "cfd/multibody/contact.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::multibody {

// Symmetric viscous stress tensor in Cartesian components. Pressure is kept
// separate so resolved traction remains t = -p n + tau n.
struct SymmetricStress3 {
    double xx{}, yy{}, zz{}, xy{}, xz{}, yz{};
    [[nodiscard]] Vec3 traction(Vec3 normal) const noexcept;
};

struct RigidSurfaceSample {
    Vec3 position{};
    Vec3 normal{1.0,0.0,0.0};
    double area{};
};

struct ResolvedSurfaceLoad {
    Vec3 force{};
    Vec3 torque{};
};

// Low-order conservative sphere quadrature. sample_count==6 uses an octahedral
// rule that exactly integrates a linear pressure field's net force.
[[nodiscard]] std::vector<RigidSurfaceSample> sphere_surface_quadrature(
    const RigidBody& body, double radius, std::size_t sample_count = 6U);

// Samples pressure/stress from the nearest FVM cell and integrates traction.
// This is intentionally geometry-agnostic: callers may supply triangulated,
// immersed-boundary, or analytical surface quadrature points.
[[nodiscard]] ResolvedSurfaceLoad integrate_resolved_surface_load(
    const cfd::fvm::PolyMesh& mesh,
    std::span<const double> cell_pressure,
    std::span<const SymmetricStress3> cell_viscous_stress,
    std::span<const RigidSurfaceSample> surface,
    Vec3 moment_center);

void apply_resolved_surface_load(RigidBody& body, const ResolvedSurfaceLoad& load) noexcept;

struct UnresolvedCfdDemConfig {
    double fluid_density{1.0};
    double dynamic_viscosity{1.0e-3};
    double minimum_void_fraction{0.05};
    double void_drag_exponent{2.65};
    double saffman_coefficient{1.615};
    bool include_pressure_gradient{true};
    bool include_saffman_lift{true};
};

struct UnresolvedCfdDemResult {
    std::vector<double> solid_volume_fraction;
    std::vector<double> void_fraction;
    std::vector<Vec3> particle_force;
    // Reaction force density applied to the fluid [force / volume].
    std::vector<Vec3> fluid_momentum_source;
};

// Unresolved many-particle coupling. Sphere volume is deposited into the
// nearest Eulerian cell; drag/pressure/lift are evaluated from that cell and
// the equal-and-opposite reaction is returned as an Eulerian force density.
// The primitive is conservative by construction before optional alpha clipping.
[[nodiscard]] UnresolvedCfdDemResult unresolved_cfd_dem_coupling(
    const cfd::fvm::PolyMesh& mesh,
    std::span<const RigidBody> bodies,
    std::span<const SphereShape> spheres,
    std::span<const cfd::fvm::Vec3> fluid_velocity,
    std::span<const cfd::fvm::Vec3> pressure_gradient,
    std::span<const cfd::fvm::Vec3> vorticity,
    const UnresolvedCfdDemConfig& config = {});

void apply_unresolved_particle_forces(
    std::span<RigidBody> bodies,
    std::span<const SphereShape> spheres,
    std::span<const Vec3> sphere_forces);

} // namespace cfd::multibody
