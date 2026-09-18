#pragma once

#include "cfd/multibody/rigid_body.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::multibody {

struct ImplicitSecondOrderConfig {
    double dt{1.0e-3};
    double beta{0.25};
    double gamma{0.5};
    // Generalized-alpha interpolation convention:
    // x_{n+1-alpha_f}=(1-alpha_f)x_{n+1}+alpha_f*x_n and
    // a_{n+1-alpha_m}=(1-alpha_m)a_{n+1}+alpha_m*a_n.
    double alpha_m{};
    double alpha_f{};
    std::size_t newton_iterations{16};
    double nonlinear_tolerance{1.0e-10};
    double finite_difference_epsilon{1.0e-7};

    [[nodiscard]] static ImplicitSecondOrderConfig newmark(double dt,double beta=0.25,double gamma=0.5);
    [[nodiscard]] static ImplicitSecondOrderConfig generalized_alpha(double dt,double spectral_radius_infinity);
    [[nodiscard]] static ImplicitSecondOrderConfig hht(double dt,double alpha);
};

using GeneralizedForceFunction=std::function<void(std::span<const double> q,
                                                  std::span<const double> v,
                                                  double time,
                                                  std::span<double> force)>;

// Dense implicit second-order integrator for generalized multibody coordinates.
// The mass matrix is constant over a step; nonlinear internal/external forces can
// depend on q, v and time and are resolved by Newton iterations.
class ImplicitGeneralizedIntegrator {
public:
    ImplicitGeneralizedIntegrator(std::vector<double> mass_matrix,
                                  std::size_t dofs,
                                  GeneralizedForceFunction force,
                                  ImplicitSecondOrderConfig config={});

    void initialize(std::vector<double> q,std::vector<double> v={},double time=0.0);
    void step();
    void run(std::size_t steps);

    [[nodiscard]] const std::vector<double>& q() const noexcept { return q_; }
    [[nodiscard]] const std::vector<double>& v() const noexcept { return v_; }
    [[nodiscard]] const std::vector<double>& a() const noexcept { return a_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t last_newton_iterations() const noexcept { return last_newton_iterations_; }
    [[nodiscard]] double last_residual_norm() const noexcept { return last_residual_norm_; }

private:
    std::vector<double> mass_{};
    std::size_t n_{};
    GeneralizedForceFunction force_{};
    ImplicitSecondOrderConfig config_{};
    std::vector<double> q_{},v_{},a_{};
    double time_{};
    std::size_t last_newton_iterations_{};
    double last_residual_norm_{};

    [[nodiscard]] std::vector<double> residual(std::span<const double> a1,
                                                std::span<const double> q0,
                                                std::span<const double> v0,
                                                std::span<const double> a0) const;
};

enum class ReducedJointType { revolute, prismatic };

struct ArticulatedLink {
    std::size_t parent{invalid_body};
    ReducedJointType joint_type{ReducedJointType::revolute};
    Vec3 parent_joint_offset{};      // parent-body local coordinates; world for root
    Vec3 joint_to_com_local{};       // child-frame vector from joint to child COM
    Vec3 axis_parent{0.0,0.0,1.0};   // joint axis in parent frame
    Quaternion reference_orientation{};
    RigidBody body{};
    double q{};
    double qd{};
    double qdd{};
    double applied_effort{};
    double damping{};
};

// Fixed-base tree reduced-coordinate baseline.  It assembles M(q)=sum J^T M_b J
// from body Jacobians and solves the dense generalized acceleration system.
// This is intentionally a readable O(n^3) baseline, not a Featherstone ABA.
class ArticulatedSystem {
public:
    explicit ArticulatedSystem(Vec3 gravity={0.0,-9.81,0.0});
    std::size_t add_link(ArticulatedLink link);
    void set_effort(std::size_t link,double effort);
    void forward_kinematics();
    void step(double dt);
    void run(std::size_t steps,double dt);

    [[nodiscard]] const std::vector<ArticulatedLink>& links() const noexcept { return links_; }
    [[nodiscard]] std::vector<ArticulatedLink>& links() noexcept { return links_; }
    [[nodiscard]] std::vector<double> mass_matrix() const;
    [[nodiscard]] std::vector<double> generalized_forces() const;
    [[nodiscard]] double kinetic_energy() const noexcept;

private:
    Vec3 gravity_{};
    std::vector<ArticulatedLink> links_{};
    std::vector<Vec3> joint_origin_world_{};
    std::vector<Vec3> joint_axis_world_{};

    [[nodiscard]] bool is_ancestor(std::size_t joint,std::size_t body) const;
    [[nodiscard]] Vec3 linear_jacobian(std::size_t joint,std::size_t body) const;
    [[nodiscard]] Vec3 angular_jacobian(std::size_t joint,std::size_t body) const;
};

struct FlexibleInterfaceNode {
    Vec3 reference_local{};
};

struct FlexibleMode {
    std::vector<Vec3> shape_local; // one displacement vector per interface node
    double modal_mass{1.0};
    double modal_stiffness{};
    double modal_damping{};
};

struct FlexibleInterfaceLoad {
    Vec3 body_force{};
    Vec3 body_torque{};
    std::vector<double> modal_force{};
};

// Conservative floating-frame interface between an FEM interface and a rigid
// multibody frame.  Body/modal state maps to nodal kinematics; nodal FEM
// reactions map back to a body wrench plus modal generalized forces.
class FlexibleBodyInterface {
public:
    FlexibleBodyInterface(std::vector<FlexibleInterfaceNode> nodes,
                          std::vector<FlexibleMode> modes={});
    void set_modal_state(std::vector<double> q,std::vector<double> qd={});

    [[nodiscard]] std::vector<Vec3> node_positions(const RigidBody& frame) const;
    [[nodiscard]] std::vector<Vec3> node_velocities(const RigidBody& frame) const;
    [[nodiscard]] FlexibleInterfaceLoad project_nodal_forces(const RigidBody& frame,
                                                              std::span<const Vec3> nodal_forces) const;
    [[nodiscard]] FlexibleInterfaceLoad apply_nodal_forces(RigidBody& frame,
                                                            std::span<const Vec3> nodal_forces) const;

    [[nodiscard]] const std::vector<double>& modal_q() const noexcept { return q_; }
    [[nodiscard]] const std::vector<double>& modal_qd() const noexcept { return qd_; }
    [[nodiscard]] const std::vector<FlexibleMode>& modes() const noexcept { return modes_; }
    [[nodiscard]] std::size_t node_count() const noexcept { return nodes_.size(); }

private:
    std::vector<FlexibleInterfaceNode> nodes_{};
    std::vector<FlexibleMode> modes_{};
    std::vector<double> q_{},qd_{};

    [[nodiscard]] Vec3 local_deformation(std::size_t node) const noexcept;
    [[nodiscard]] Vec3 local_deformation_velocity(std::size_t node) const noexcept;
};

} // namespace cfd::multibody
