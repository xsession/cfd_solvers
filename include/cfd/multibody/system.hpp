#pragma once

#include "cfd/multibody/constraints.hpp"
#include "cfd/multibody/contact.hpp"

namespace cfd::multibody {

enum class ContactMethod { smooth_penalty, nonsmooth_impulse };

struct MultibodyConfig {
    double dt{1.0e-3};
    Vec3 gravity{0.0,-9.81,0.0};
    RigidBodyIntegrator integrator{RigidBodyIntegrator::semi_implicit_euler};
    ContactMethod contact_method{ContactMethod::smooth_penalty};
    std::size_t constraint_iterations{16};
    double constraint_baumgarte{0.2};
    PenaltyContactModel penalty{};
    ImpulseContactModel impulse{};
};

class RigidBodySystem {
public:
    explicit RigidBodySystem(MultibodyConfig config={});
    std::size_t add_body(RigidBody body);
    void add_sphere(std::size_t body,double radius);
    void add_plane(PlaneShape plane);
    void add_triangle_mesh(TriangleMeshShape mesh);
    void add_joint(Joint joint);
    void step();
    void run(std::size_t steps);

    [[nodiscard]] std::vector<RigidBody>& bodies() noexcept { return bodies_; }
    [[nodiscard]] const std::vector<RigidBody>& bodies() const noexcept { return bodies_; }
    [[nodiscard]] const std::vector<SphereShape>& spheres() const noexcept { return spheres_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double total_kinetic_energy() const noexcept;
    [[nodiscard]] Vec3 linear_momentum() const noexcept;

private:
    MultibodyConfig config_{};
    std::vector<RigidBody> bodies_{};
    std::vector<SphereShape> spheres_{};
    std::vector<PlaneShape> planes_{};
    std::vector<TriangleMeshShape> triangle_meshes_{};
    std::vector<Joint> joints_{};
    double time_{};
};

struct DemConfig {
    double dt{1.0e-4};
    Vec3 gravity{0.0,-9.81,0.0};
    HertzMindlinContactModel contact{};
    bool history_dependent_tangential{true};
};

class ExplicitDemSystem {
public:
    explicit ExplicitDemSystem(DemConfig config={});
    std::size_t add_particle(Vec3 position,double radius,double density,Vec3 velocity={});
    void add_plane(PlaneShape plane);
    void add_triangle_mesh(TriangleMeshShape mesh);
    std::size_t add_bond(std::size_t body_a,std::size_t body_b,BondedParticleModel model={},double rest_length=0.0);
    void step();
    void run(std::size_t steps);
    [[nodiscard]] const std::vector<RigidBody>& bodies() const noexcept { return bodies_; }
    [[nodiscard]] const std::vector<SphereShape>& spheres() const noexcept { return spheres_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] const ContactManifoldCache& contact_history() const noexcept { return contact_history_; }
    [[nodiscard]] const std::vector<ParticleBond>& bonds() const noexcept { return bonds_; }
    [[nodiscard]] const BondUpdateStats& last_bond_stats() const noexcept { return last_bond_stats_; }

private:
    DemConfig config_{};
    std::vector<RigidBody> bodies_{};
    std::vector<SphereShape> spheres_{};
    std::vector<PlaneShape> planes_{};
    std::vector<TriangleMeshShape> triangle_meshes_{};
    ContactManifoldCache contact_history_{};
    std::vector<ParticleBond> bonds_{};
    BondUpdateStats last_bond_stats_{};
    double time_{};
};

// Simple two-way coupling primitive for spherical particles. The returned force
// is applied to the body; its negative is the reaction that a CFD cell/parcel
// coupling layer must receive to conserve momentum.
[[nodiscard]] Vec3 stokes_drag_force(const RigidBody& body,double radius,Vec3 fluid_velocity,double dynamic_viscosity);

} // namespace cfd::multibody
