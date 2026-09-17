#pragma once

#include "cfd/multibody/rigid_body.hpp"
#include <array>
#include <cstdint>
#include <utility>

namespace cfd::multibody {

struct Aabb { Vec3 minimum{},maximum{}; };
struct ContactPoint {
    std::size_t body_a{invalid_body};
    std::size_t body_b{invalid_body};
    Vec3 point{};
    Vec3 normal{1.0,0.0,0.0}; // outward force direction on body A
    double penetration{};
    std::uint64_t feature_id{};
};

// Convex vertices are stored in body-local coordinates.  This is intentionally
// a small clean-room shape primitive; mesh collision remains a separate item.
struct ConvexHullShape {
    std::size_t body{invalid_body};
    std::vector<Vec3> vertices{};
};

struct TriangleMeshShape {
    std::size_t body{invalid_body}; // invalid_body means a static world-space mesh
    std::vector<Vec3> vertices{};
    std::vector<std::array<std::size_t,3>> triangles{};
    std::uint32_t feature_namespace{};
};

struct GjkEpaResult {
    bool intersect{};
    ContactPoint contact{};
    std::size_t gjk_iterations{};
    std::size_t epa_iterations{};
};

struct PenaltyContactModel {
    double normal_stiffness{1.0e5};
    double normal_damping{100.0};
    double tangential_damping{50.0};
    double friction{0.5};
};
struct HertzMindlinContactModel {
    double normal_stiffness{1.0e5};
    double normal_damping{20.0};
    double tangential_stiffness{5.0e4};
    double tangential_damping{20.0};
    double friction{0.5};
    double rolling_resistance{0.01};
    double cohesion_force{};
};
struct ImpulseContactModel {
    double restitution{0.0};
    double friction{0.5};
    double baumgarte{0.15};
    std::size_t iterations{12};
};

// History-bearing bonded-particle law.  Strengths are force thresholds; this
// keeps the baseline independent of any assumed bond cross-section.
struct BondedParticleModel {
    double normal_stiffness{1.0e5};
    double shear_stiffness{5.0e4};
    double normal_damping{50.0};
    double shear_damping{25.0};
    double tensile_failure_force{1.0e3};
    double shear_failure_force{1.0e3};
    double damage_onset_ratio{0.7};
};

struct ParticleBond {
    std::size_t body_a{invalid_body};
    std::size_t body_b{invalid_body};
    double rest_length{};
    BondedParticleModel model{};
    Vec3 tangential_displacement{};
    double damage{};
    bool broken{};
};

struct BondUpdateStats {
    std::size_t active{};
    std::size_t damaged{};
    std::size_t broken_total{};
    std::size_t broken_this_step{};
};

struct PersistentContactState {
    ContactPoint contact{};
    Vec3 tangential_displacement{};
    std::size_t age{};
};

class ContactManifoldCache {
public:
    // Reconcile one contact set with persistent pair/feature state.  Tangential
    // history is projected onto the current tangent plane to remain frame-safe.
    void update(const std::vector<ContactPoint>& contacts);
    void clear() noexcept { states_.clear(); }
    [[nodiscard]] std::vector<PersistentContactState>& states() noexcept { return states_; }
    [[nodiscard]] const std::vector<PersistentContactState>& states() const noexcept { return states_; }
private:
    std::vector<PersistentContactState> states_{};
};

[[nodiscard]] Aabb sphere_aabb(const RigidBody& body,const SphereShape& sphere);
[[nodiscard]] Aabb convex_aabb(const RigidBody& body,const ConvexHullShape& shape);
[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> broad_phase_sweep_and_prune(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres);
[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> broad_phase_bvh(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres);
[[nodiscard]] std::vector<std::pair<std::size_t,std::size_t>> broad_phase_cell_linked(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,double cell_size);
[[nodiscard]] std::vector<ContactPoint> sphere_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const std::vector<std::pair<std::size_t,std::size_t>>& candidate_pairs);
[[nodiscard]] std::vector<ContactPoint> sphere_plane_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const std::vector<PlaneShape>& planes);
[[nodiscard]] std::vector<ContactPoint> sphere_triangle_mesh_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const TriangleMeshShape& mesh);
[[nodiscard]] GjkEpaResult convex_gjk_epa_contact(const std::vector<RigidBody>& bodies,const ConvexHullShape& a,const ConvexHullShape& b,std::size_t max_gjk_iterations=32,std::size_t max_epa_iterations=64,double tolerance=1.0e-9);

void apply_penalty_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const PenaltyContactModel& model);
void apply_hertz_mindlin_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const std::vector<SphereShape>& spheres,const HertzMindlinContactModel& model);
void apply_history_dependent_mindlin_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const std::vector<SphereShape>& spheres,ContactManifoldCache& cache,double dt,const HertzMindlinContactModel& model);
[[nodiscard]] ParticleBond make_particle_bond(const std::vector<RigidBody>& bodies,std::size_t body_a,std::size_t body_b,BondedParticleModel model={},double rest_length=0.0);
[[nodiscard]] BondUpdateStats apply_particle_bonds(std::vector<RigidBody>& bodies,std::vector<ParticleBond>& bonds,double dt);
void resolve_impulse_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,double dt,const ImpulseContactModel& model);

} // namespace cfd::multibody
