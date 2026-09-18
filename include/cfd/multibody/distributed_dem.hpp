#pragma once

#include "cfd/multibody/contact.hpp"
#include "cfd/multibody/rigid_body.hpp"

#include <cstdint>
#include <span>
#include <vector>

#if defined(CFD_HAS_MPI)
#include <mpi.h>
#endif

namespace cfd::multibody {

struct DistributedDemParticle {
    std::uint64_t global_id{};
    RigidBodyState state{};
    double mass{1.0};
    Vec3 inertia_diagonal{1.0,1.0,1.0};
    double radius{0.5};
    int owner_rank{};
};

struct DemSlabDecomposition {
    double x_min{};
    double x_max{1.0};
    int ranks{1};

    void validate() const;
    [[nodiscard]] double slab_width() const;
    [[nodiscard]] int owner_rank(double x) const;
    [[nodiscard]] double slab_min(int rank) const;
    [[nodiscard]] double slab_max(int rank) const;
};

struct DemExchangePlan {
    std::vector<std::vector<std::size_t>> migrate_indices_by_rank{};
    std::vector<std::vector<std::size_t>> ghost_indices_by_rank{};
};

// Persistent cross-rank contact state. Pair IDs are canonicalized such that
// particle_a < particle_b. The owner is therefore the rank owning particle_a.
struct DistributedDemContactHistory {
    std::uint64_t particle_a{};
    std::uint64_t particle_b{};
    Vec3 tangential_displacement{};
    std::size_t age{};
};

struct DistributedDemBond {
    std::uint64_t particle_a{};
    std::uint64_t particle_b{};
    double rest_length{};
    BondedParticleModel model{};
    Vec3 tangential_displacement{};
    double damage{};
    bool broken{};
};

struct DistributedDemForceContribution {
    std::uint64_t global_id{};
    int owner_rank{};
    Vec3 force{};
    Vec3 torque{};
};

struct DistributedDemInteractionResult {
    std::vector<DistributedDemForceContribution> local{};
    std::vector<DistributedDemForceContribution> remote{};
    std::size_t contacts{};
    std::size_t active_bonds{};
    std::size_t broken_bonds{};
};

struct DistributedDemStepStats {
    std::size_t owned{};
    std::size_t ghosts{};
    std::size_t contacts{};
    std::size_t active_bonds{};
    std::size_t broken_bonds{};
    std::size_t remote_force_messages{};
};

[[nodiscard]] DistributedDemParticle make_distributed_dem_particle(
    std::uint64_t global_id,
    const RigidBody& body,
    const SphereShape& sphere,
    int owner_rank);

[[nodiscard]] DistributedDemBond make_distributed_dem_bond(
    const DistributedDemParticle& a,
    const DistributedDemParticle& b,
    BondedParticleModel model={},
    double rest_length=0.0);

// Deterministic planning layer used by both tests and MPI exchange. Ownership
// is based on particle centers; ghosting expands the local sphere AABB by the
// requested halo width so cross-rank contacts are visible on adjacent slabs.
[[nodiscard]] DemExchangePlan plan_dem_slab_exchange(
    std::span<const DistributedDemParticle> particles,
    const DemSlabDecomposition& decomposition,
    int local_rank,
    double ghost_width);

// Evaluate local-local and local-ghost contacts. For a cross-rank pair only
// the rank owning the lower global particle ID evaluates the contact and owns
// its persistent Mindlin state. Remote reaction forces are returned explicitly
// for reverse exchange to the partner owner.
[[nodiscard]] DistributedDemInteractionResult evaluate_distributed_dem_contacts(
    std::span<const DistributedDemParticle> owned,
    std::span<const DistributedDemParticle> ghosts,
    std::vector<DistributedDemContactHistory>& history,
    int local_rank,
    double dt,
    const HertzMindlinContactModel& model);

// Bond records follow the same lower-global-ID ownership rule. Remote partner
// state must be available either as a local particle or a ghost copy.
[[nodiscard]] DistributedDemInteractionResult evaluate_distributed_dem_bonds(
    std::span<const DistributedDemParticle> owned,
    std::span<const DistributedDemParticle> ghosts,
    std::vector<DistributedDemBond>& bonds,
    int local_rank,
    double dt);

void accumulate_distributed_dem_forces(
    std::vector<DistributedDemForceContribution>& destination,
    std::span<const DistributedDemForceContribution> source);

// Apply accumulated forces and advance locally-owned spherical particles. This
// is intentionally transport-agnostic and is used by both emulated-rank tests
// and the MPI timestep driver.
void integrate_distributed_dem_particles(
    std::span<DistributedDemParticle> owned,
    std::span<const DistributedDemForceContribution> contributions,
    Vec3 gravity,
    double dt);

#if defined(CFD_HAS_MPI)
class MpiDemDomainExchange {
public:
    MpiDemDomainExchange(DemSlabDecomposition decomposition,MPI_Comm communicator=MPI_COMM_WORLD);
    ~MpiDemDomainExchange();
    MpiDemDomainExchange(const MpiDemDomainExchange&)=delete;
    MpiDemDomainExchange& operator=(const MpiDemDomainExchange&)=delete;

    // Migrates ownership first, then returns read-only ghost copies required
    // by neighboring slabs for contact detection. Global IDs are preserved.
    [[nodiscard]] std::vector<DistributedDemParticle> exchange(
        std::vector<DistributedDemParticle>& owned,
        double ghost_width);

    // State-aware exchange: persistent pair/bond history owned by a migrating
    // lower-ID particle follows that particle to its new owner rank.
    [[nodiscard]] std::vector<DistributedDemParticle> exchange(
        std::vector<DistributedDemParticle>& owned,
        std::vector<DistributedDemContactHistory>& contact_history,
        std::vector<DistributedDemBond>& bonds,
        double ghost_width);

    // Reverse-exchange equal-and-opposite forces generated by local-ghost
    // interactions. Returned contributions target particles owned by this rank.
    [[nodiscard]] std::vector<DistributedDemForceContribution> exchange_forces(
        std::span<const DistributedDemForceContribution> remote_contributions);

    // Complete distributed CPU DEM baseline. Migration/ghosting occurs before
    // each force evaluation; particles that cross slabs are migrated on the
    // next call. The MPI implementation is deliberately separate from the
    // resident SYCL solver so GPU-aware transport can use the same records.
    [[nodiscard]] DistributedDemStepStats step(
        std::vector<DistributedDemParticle>& owned,
        std::vector<DistributedDemContactHistory>& contact_history,
        std::vector<DistributedDemBond>& bonds,
        double ghost_width,
        double dt,
        Vec3 gravity,
        const HertzMindlinContactModel& contact_model);

    [[nodiscard]] int rank() const noexcept { return rank_; }
    [[nodiscard]] int size() const noexcept { return size_; }
    [[nodiscard]] MPI_Comm communicator() const noexcept { return communicator_; }
    [[nodiscard]] const DemSlabDecomposition& decomposition() const noexcept { return decomposition_; }

private:
    DemSlabDecomposition decomposition_{};
    MPI_Comm communicator_{MPI_COMM_NULL};
    int rank_{};
    int size_{1};
};
#endif

} // namespace cfd::multibody
