#pragma once

#include "cfd/core/device_residency.hpp"
#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/multibody/resident_dem_sycl.hpp"

#if defined(CFD_HAS_MPI) && defined(CFD_HAS_SYCL)
#include <mpi.h>
#include <sycl/sycl.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cfd::multibody {

enum class ResidentDemMpiMode {
    staged_host,
    direct_device
};

struct ResidentDemMpiExchangeConfig {
    DemSlabDecomposition decomposition{};
    ResidentDemMpiMode mode{ResidentDemMpiMode::staged_host};
    // Direct-device mode deliberately requires an explicit opt-in. The MPI
    // standard does not itself guarantee that an implementation accepts SYCL
    // USM pointers, so callers/CI must only enable this after runtime probing.
    bool assume_gpu_aware_mpi{};
    // Maximum resident bond records owned by this rank. Zero selects a
    // conservative default of 4 * particle capacity when storage is reserved.
    std::size_t bond_capacity{};
};

struct ResidentDemMpiExchangeStats {
    std::size_t migrated_out{};
    std::size_t migrated_in{};
    std::size_t ghosts_out{};
    std::size_t ghosts_in{};
    std::size_t staged_host_bytes{};
    std::size_t direct_device_bytes{};
    std::size_t scalar_synchronizations{};
    std::size_t reverse_forces_out{};
    std::size_t reverse_forces_in{};
    std::size_t nonblocking_exchanges{};
    std::size_t resident_bonds{};
};

// GPU-aware slab exchange for ResidentDemSycl. Particle/history payloads are
// packed on-device. staged_host mode transfers only packed boundary/migrant
// records through pinned-by-runtime host vectors; direct_device passes USM
// pointers directly to MPI and therefore requires a GPU-aware MPI stack.
//
// Received ghosts remain in a separate device buffer rather than polluting the
// owned particle array. step() uses that separation to overlap nonblocking halo
// exchange with interior contacts, then evaluates boundary/local-ghost contacts
// and reverse-exchanges compact force/torque records before integration.
class MpiResidentDemSyclExchange {
public:
    MpiResidentDemSyclExchange(
        sycl::queue& communication_queue,
        ResidentDemMpiExchangeConfig config,
        MPI_Comm communicator=MPI_COMM_WORLD);
    ~MpiResidentDemSyclExchange();
    MpiResidentDemSyclExchange(const MpiResidentDemSyclExchange&)=delete;
    MpiResidentDemSyclExchange& operator=(const MpiResidentDemSyclExchange&)=delete;

    // Migrate particles/history to their owning slab, compact the local
    // resident array on-device, append incoming migrants on-device, then pack
    // and exchange adjacent-slab ghosts. Returns transport statistics.
    [[nodiscard]] ResidentDemMpiExchangeStats exchange(
        ResidentDemSycl& dem,
        double ghost_width);

    // End-to-end distributed resident timestep. Migration is completed first,
    // then halo data moves nonblocking while interior owned-owned contacts run
    // on the device. Boundary/local-ghost contacts and bonds emit compact
    // reverse-force records, which are exchanged before final integration.
    [[nodiscard]] ResidentDemMpiExchangeStats step(
        ResidentDemSycl& dem,
        double ghost_width);

    // Bond records are owned by the rank owning the lower global particle ID.
    // When that particle migrates, the persistent damage/tangential state is
    // migrated with it on-device. Upload/download are setup/checkpoint seams,
    // not part of the ordinary distributed timestep.
    void upload_bonds(std::span<const ResidentDemPackedBond> bonds);
    [[nodiscard]] std::vector<ResidentDemPackedBond> download_bonds() const;
    [[nodiscard]] std::size_t bond_count() const noexcept { return bond_count_; }

    [[nodiscard]] const ResidentDemPackedParticle* ghost_particles_device() const noexcept { return ghost_receive_; }
    [[nodiscard]] std::size_t ghost_particle_count() const noexcept { return ghost_count_; }
    // One byte per currently-owned particle: 1 means the particle belongs to
    // the halo/boundary region and 0 means it can participate in interior work
    // while communication is in flight.
    [[nodiscard]] const std::uint8_t* boundary_mask_device() const noexcept { return boundary_mask_; }
    [[nodiscard]] std::size_t boundary_mask_count() const noexcept { return boundary_mask_count_; }

    [[nodiscard]] int rank() const noexcept { return rank_; }
    [[nodiscard]] int size() const noexcept { return size_; }
    [[nodiscard]] ResidentDemMpiMode mode() const noexcept { return config_.mode; }
    [[nodiscard]] const ResidentDemMpiExchangeStats& last_stats() const noexcept { return stats_; }

private:
    sycl::queue& queue_;
    ResidentDemMpiExchangeConfig config_{};
    MPI_Comm communicator_{MPI_COMM_NULL};
    int rank_{};
    int size_{1};
    std::size_t reserved_capacity_{};
    std::size_t reserved_slots_{};
    std::size_t ghost_count_{};
    std::size_t boundary_mask_count_{};
    std::size_t remote_force_capacity_{};
    std::size_t bond_capacity_{};
    std::size_t bond_count_{};
    ResidentDemMpiExchangeStats stats_{};

    std::uint32_t* migrate_counts_{};
    std::uint32_t* ghost_counts_{};
    std::uint32_t* receive_counts_{};
    std::uint32_t* migrate_indices_{};
    std::uint32_t* ghost_indices_{};
    std::uint8_t* boundary_mask_{};
    std::uint32_t* remote_force_count_{};
    std::uint32_t* reverse_counts_{};
    std::uint32_t* reverse_cursor_{};
    std::uint32_t* reverse_offsets_{};

    ResidentDemPackedParticle* migrate_send_{};
    ResidentDemPackedParticle* migrate_receive_{};
    ResidentDemPackedHistory* migrate_history_send_{};
    ResidentDemPackedHistory* migrate_history_receive_{};
    ResidentDemPackedParticle* ghost_send_{};
    ResidentDemPackedParticle* ghost_receive_{};
    ResidentDemPackedForce* remote_force_raw_{};
    ResidentDemPackedForce* reverse_force_send_{};
    ResidentDemPackedForce* reverse_force_receive_{};
    ResidentDemPackedBond* bonds_{};
    ResidentDemPackedBond* bond_send_{};
    ResidentDemPackedBond* bond_receive_{};
    std::uint32_t* bond_counts_{};
    std::uint32_t* bond_receive_counts_{};
    std::uint32_t* bond_compaction_count_{};

    void reserve(const ResidentDemSycl& dem);
    void reserve_bonds(std::size_t capacity);
    void release() noexcept;
    void clear_counts();
    void classify_migration(const ResidentDemSycl& dem);
    void classify_ghosts(const ResidentDemSycl& dem,double ghost_width);
    void pack_particle_segments(
        ResidentDemSycl& dem,
        const std::uint32_t* indices,
        const std::uint32_t* counts,
        ResidentDemPackedParticle* particles,
        ResidentDemPackedHistory* history);
    void exchange_migrants(ResidentDemSycl& dem);
    void exchange_ghosts(ResidentDemSycl& dem,double ghost_width);
    void exchange_bonds_for_migration(ResidentDemSycl& dem);
    void exchange_reverse_forces(ResidentDemSycl& dem);
};

} // namespace cfd::multibody
#endif
