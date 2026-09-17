#pragma once

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

[[nodiscard]] DistributedDemParticle make_distributed_dem_particle(
    std::uint64_t global_id,
    const RigidBody& body,
    const SphereShape& sphere,
    int owner_rank);

// Deterministic planning layer used by both tests and MPI exchange.  Ownership
// is based on particle centers; ghosting expands the local sphere AABB by the
// requested halo width so cross-rank contacts are visible on adjacent slabs.
[[nodiscard]] DemExchangePlan plan_dem_slab_exchange(
    std::span<const DistributedDemParticle> particles,
    const DemSlabDecomposition& decomposition,
    int local_rank,
    double ghost_width);

#if defined(CFD_HAS_MPI)
class MpiDemDomainExchange {
public:
    MpiDemDomainExchange(DemSlabDecomposition decomposition,MPI_Comm communicator=MPI_COMM_WORLD);
    ~MpiDemDomainExchange();
    MpiDemDomainExchange(const MpiDemDomainExchange&)=delete;
    MpiDemDomainExchange& operator=(const MpiDemDomainExchange&)=delete;

    // Migrates ownership first, then returns read-only ghost copies required
    // by neighboring slabs for contact detection.  Global IDs are preserved.
    [[nodiscard]] std::vector<DistributedDemParticle> exchange(
        std::vector<DistributedDemParticle>& owned,
        double ghost_width);

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
