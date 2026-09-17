#pragma once

#include "cfd/multibody/rigid_body.hpp"
#include <variant>

namespace cfd::multibody {

struct ConstraintRow {
    std::size_t body_a{invalid_body};
    std::size_t body_b{invalid_body};
    Vec3 linear_a{};
    Vec3 angular_a{};
    Vec3 linear_b{};
    Vec3 angular_b{};
    double bias{};
    double lower{-std::numeric_limits<double>::infinity()};
    double upper{ std::numeric_limits<double>::infinity()};
    double impulse{};
};

struct DistanceConstraint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_anchor_a{},local_anchor_b{};
    double distance{};
};
struct SphericalJoint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_anchor_a{},local_anchor_b{};
};
struct RevoluteJoint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_anchor_a{},local_anchor_b{};
    Vec3 local_axis_a{0.0,0.0,1.0},local_axis_b{0.0,0.0,1.0};
};
struct PrismaticJoint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_anchor_a{},local_anchor_b{};
    Vec3 local_axis_a{1.0,0.0,0.0};
    Quaternion reference_relative_orientation{};
};
struct FixedJoint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_anchor_a{},local_anchor_b{};
    Quaternion reference_relative_orientation{};
};
struct GearConstraint {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_axis_a{0.0,0.0,1.0},local_axis_b{0.0,0.0,1.0};
    double ratio{1.0};
};

// Velocity motors are deliberately represented as bounded constraint rows.  The
// force/torque limits are converted to per-step impulse limits when rows are
// assembled, making them compatible with the existing projected solver.
struct LinearMotor {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_axis_a{1.0,0.0,0.0};
    double target_speed{};
    double max_force{std::numeric_limits<double>::infinity()};
};
struct AngularMotor {
    std::size_t body_a{invalid_body},body_b{invalid_body};
    Vec3 local_axis_a{0.0,0.0,1.0};
    double target_speed{};
    double max_torque{std::numeric_limits<double>::infinity()};
};

using Joint=std::variant<DistanceConstraint,SphericalJoint,RevoluteJoint,PrismaticJoint,FixedJoint,GearConstraint,LinearMotor,AngularMotor>;

[[nodiscard]] DistanceConstraint make_distance_constraint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_anchor_a,Vec3 world_anchor_b,double distance=-1.0);
[[nodiscard]] SphericalJoint make_spherical_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_anchor);
[[nodiscard]] RevoluteJoint make_revolute_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_anchor,Vec3 world_axis);
[[nodiscard]] PrismaticJoint make_prismatic_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_anchor,Vec3 world_axis);
[[nodiscard]] FixedJoint make_fixed_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_anchor);
[[nodiscard]] GearConstraint make_gear_constraint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_axis_a,Vec3 world_axis_b,double ratio);
[[nodiscard]] LinearMotor make_linear_motor(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_axis,double target_speed,double max_force=std::numeric_limits<double>::infinity());
[[nodiscard]] AngularMotor make_angular_motor(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 world_axis,double target_speed,double max_torque=std::numeric_limits<double>::infinity());

[[nodiscard]] std::vector<ConstraintRow> build_constraint_rows(const std::vector<RigidBody>& bodies,const Joint& joint,double dt,double baumgarte=0.2);
void solve_constraint_rows(std::vector<RigidBody>& bodies,std::vector<ConstraintRow>& rows,std::size_t iterations=12);

} // namespace cfd::multibody
