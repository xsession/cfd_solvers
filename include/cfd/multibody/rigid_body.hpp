#pragma once

#include "cfd/multibody/math.hpp"
#include <cstddef>
#include <limits>
#include <vector>

namespace cfd::multibody {

inline constexpr std::size_t invalid_body=std::numeric_limits<std::size_t>::max();

struct RigidBodyState {
    Vec3 position{};
    Quaternion orientation{};
    Vec3 linear_velocity{};
    Vec3 angular_velocity{};
};

struct RigidBody {
    RigidBodyState state{};
    double mass{1.0};
    Vec3 inertia_diagonal{1.0,1.0,1.0};
    bool fixed{};
    Vec3 force{};
    Vec3 torque{};

    [[nodiscard]] double inverse_mass() const noexcept { return fixed ? 0.0 : 1.0/mass; }
    [[nodiscard]] Vec3 world_inverse_inertia_mul(Vec3 world_vector) const;
    [[nodiscard]] Vec3 point_velocity(Vec3 world_point) const noexcept;
    void add_force(Vec3 f) noexcept { force+=f; }
    void add_torque(Vec3 t) noexcept { torque+=t; }
    void add_force_at_point(Vec3 f,Vec3 world_point) noexcept;
    void apply_impulse(Vec3 impulse,Vec3 world_point) noexcept;
    void clear_accumulators() noexcept { force={};torque={}; }
    void validate() const;
};

enum class RigidBodyIntegrator { semi_implicit_euler, velocity_verlet };

void integrate_rigid_body(RigidBody& body,double dt,RigidBodyIntegrator integrator=RigidBodyIntegrator::semi_implicit_euler);

struct SphereShape {
    std::size_t body{invalid_body};
    double radius{0.5};
};

struct PlaneShape {
    Vec3 normal{0.0,1.0,0.0};
    double offset{}; // dot(normal,x)=offset
};

} // namespace cfd::multibody
