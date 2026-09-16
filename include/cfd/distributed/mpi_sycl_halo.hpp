#pragma once

#include "cfd/distributed/mpi_runtime.hpp"
#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/lbm/distributed_sycl_pull.hpp"

#if defined(CFD_HAS_MPI)
#include <mpi.h>
#endif
#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#endif

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd::distributed {

#if defined(CFD_HAS_MPI) && defined(CFD_HAS_SYCL)

enum class GpuMpiMode {
    staged_host,
    direct_device
};

template<class Descriptor>
class MpiSyclSelectiveHaloExchange {
public:
    using Block = cfd::lbm::DistributedSyclPullBlock<Descriptor>;
    static constexpr std::size_t q = static_cast<std::size_t>(Descriptor::q);

    MpiSyclSelectiveHaloExchange(const MpiCartesianRuntime& runtime,
                                 Block& block,
                                 GpuMpiMode mode = GpuMpiMode::staged_host)
        : runtime_(runtime), block_(block), mode_(mode) {
        try {
            slots_.reserve(26);
            for (const auto offset : neighbor_offsets()) {
                std::size_t direction_count = 0;
                const auto directions = crossing_population_array<Descriptor>(offset, direction_count);
                if (direction_count == 0) continue;
                const int neighbor_rank = runtime_.neighbor(offset);
                if (neighbor_rank < 0) continue;
                slots_.emplace_back();
                auto& slot = slots_.back();
                slot.offset = offset;
                slot.neighbor = neighbor_rank;
                slot.directions = directions;
                slot.direction_count = direction_count;
                slot.cells_per_direction = boundary_cells(block_.brick().extent, offset);
                slot.count = direction_count * slot.cells_per_direction;
                allocate_slot(slot);
            }
            initialize_persistent_requests();
        } catch (...) {
            release_requests();
            for (auto& slot : slots_) free_slot(slot);
            throw;
        }
    }

    ~MpiSyclSelectiveHaloExchange() {
        try {
            if (active_) wait();
            block_.queue().wait_and_throw();
        } catch (...) {}
        release_requests();
        for (auto& slot : slots_) free_slot(slot);
    }

    MpiSyclSelectiveHaloExchange(const MpiSyclSelectiveHaloExchange&) = delete;
    MpiSyclSelectiveHaloExchange& operator=(const MpiSyclSelectiveHaloExchange&) = delete;

    void begin() {
        if (active_) throw std::logic_error("MPI SYCL halo exchange already active");
        for (auto& slot : slots_) enqueue_pack(slot);
        if (mode_ == GpuMpiMode::staged_host) {
            for (auto& slot : slots_) {
                block_.queue().memcpy(slot.host_send, slot.device_send, slot.count * sizeof(float));
            }
        }
        // MPI must not observe send buffers before the device pack/staging copy completes.
        block_.queue().wait_and_throw();
        if (!requests_.empty()) {
            check_mpi(MPI_Startall(mpi_count(requests_.size()), requests_.data()),
                      "MPI_Startall(SYCL selective halo)");
        }
        active_ = true;
    }

    void wait() {
        if (!active_) return;
        if (!requests_.empty()) {
            check_mpi(MPI_Waitall(mpi_count(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE),
                      "MPI_Waitall(SYCL selective halo)");
        }
        if (mode_ == GpuMpiMode::staged_host) {
            for (auto& slot : slots_) {
                block_.queue().memcpy(slot.device_recv, slot.host_recv, slot.count * sizeof(float));
            }
        }
        for (auto& slot : slots_) enqueue_unpack(slot);
        // No queue wait here: the in-order queue lets boundary kernels depend on unpack implicitly.
        active_ = false;
    }

    [[nodiscard]] GpuMpiMode mode() const noexcept { return mode_; }
    [[nodiscard]] std::size_t message_count() const noexcept { return slots_.size(); }
    [[nodiscard]] std::size_t values_per_step() const noexcept {
        std::size_t values = 0;
        for (const auto& slot : slots_) values += slot.count;
        return values;
    }
    [[nodiscard]] std::size_t bytes_per_step() const noexcept { return values_per_step() * sizeof(float); }

private:
    struct Slot {
        NeighborOffset offset{};
        int neighbor{-1};
        std::array<int, q> directions{};
        std::size_t direction_count{};
        std::size_t cells_per_direction{};
        std::size_t count{};
        float* device_send{nullptr};
        float* device_recv{nullptr};
        float* host_send{nullptr};
        float* host_recv{nullptr};
    };

    const MpiCartesianRuntime& runtime_;
    Block& block_;
    GpuMpiMode mode_{GpuMpiMode::staged_host};
    std::vector<Slot> slots_;
    std::vector<MPI_Request> requests_;
    bool active_{false};

    void allocate_slot(Slot& slot) {
        auto& queue = block_.queue();
        slot.device_send = sycl::malloc_device<float>(slot.count, queue);
        slot.device_recv = sycl::malloc_device<float>(slot.count, queue);
        if (!slot.device_send || !slot.device_recv) {
            if (slot.device_send) sycl::free(slot.device_send, queue);
            if (slot.device_recv) sycl::free(slot.device_recv, queue);
            slot.device_send = nullptr;
            slot.device_recv = nullptr;
            throw std::bad_alloc{};
        }
        if (mode_ == GpuMpiMode::staged_host) {
            slot.host_send = sycl::malloc_host<float>(slot.count, queue);
            slot.host_recv = sycl::malloc_host<float>(slot.count, queue);
            if (!slot.host_send || !slot.host_recv) {
                if (slot.host_send) sycl::free(slot.host_send, queue);
                if (slot.host_recv) sycl::free(slot.host_recv, queue);
                sycl::free(slot.device_send, queue);
                sycl::free(slot.device_recv, queue);
                slot = {};
                throw std::bad_alloc{};
            }
        }
    }

    void free_slot(Slot& slot) noexcept {
        auto& queue = block_.queue();
        try { if (slot.host_send) sycl::free(slot.host_send, queue); } catch (...) {}
        try { if (slot.host_recv) sycl::free(slot.host_recv, queue); } catch (...) {}
        try { if (slot.device_send) sycl::free(slot.device_send, queue); } catch (...) {}
        try { if (slot.device_recv) sycl::free(slot.device_recv, queue); } catch (...) {}
        slot = {};
    }

    void initialize_persistent_requests() {
        requests_.resize(slots_.size() * 2U, MPI_REQUEST_NULL);
        std::size_t request = 0;
        for (auto& slot : slots_) {
            float* recv = mode_ == GpuMpiMode::staged_host ? slot.host_recv : slot.device_recv;
            check_mpi(MPI_Recv_init(recv, mpi_count(slot.count), MPI_FLOAT,
                                    slot.neighbor, neighbor_tag(slot.offset.opposite()),
                                    runtime_.communicator(), &requests_[request++]),
                      "MPI_Recv_init(SYCL selective halo)");
        }
        for (auto& slot : slots_) {
            float* send = mode_ == GpuMpiMode::staged_host ? slot.host_send : slot.device_send;
            check_mpi(MPI_Send_init(send, mpi_count(slot.count), MPI_FLOAT,
                                    slot.neighbor, neighbor_tag(slot.offset),
                                    runtime_.communicator(), &requests_[request++]),
                      "MPI_Send_init(SYCL selective halo)");
        }
    }

    void release_requests() noexcept {
        for (auto& request : requests_) {
            if (request != MPI_REQUEST_NULL) MPI_Request_free(&request);
        }
    }

    void enqueue_pack(Slot& slot) {
        const auto e = block_.brick().extent;
        const auto total = block_.total_extent();
        const std::size_t total_cells = block_.total_cells();
        const auto offset = slot.offset;
        const float* current = block_.current_data();
        float* output = slot.device_send;
        block_.queue().parallel_for(sycl::range<1>(slot.count), [=](sycl::id<1> gid) {
            const std::size_t p = gid[0];
            const auto entry = selective_send_entry<Descriptor>(e, offset, p);
            const std::size_t cell = (entry.z * total.y + entry.y) * total.x + entry.x;
            output[p] = current[static_cast<std::size_t>(entry.direction) * total_cells + cell];
        });
    }

    void enqueue_unpack(Slot& slot) {
        const auto e = block_.brick().extent;
        const auto total = block_.total_extent();
        const std::size_t total_cells = block_.total_cells();
        const auto incoming = slot.offset;
        const float* input = slot.device_recv;
        float* current = block_.current_data();
        block_.queue().parallel_for(sycl::range<1>(slot.count), [=](sycl::id<1> gid) {
            const std::size_t p = gid[0];
            const auto entry = selective_receive_entry<Descriptor>(e, incoming, p);
            const std::size_t cell = (entry.z * total.y + entry.y) * total.x + entry.x;
            current[static_cast<std::size_t>(entry.direction) * total_cells + cell] = input[p];
        });
    }

    static int mpi_count(std::size_t n) {
        if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("MPI SYCL halo message exceeds INT_MAX elements");
        }
        return static_cast<int>(n);
    }
};

#endif

} // namespace cfd::distributed
