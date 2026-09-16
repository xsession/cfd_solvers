#pragma once

#include "cfd/distributed/mpi_runtime.hpp"
#include "cfd/distributed/selective_halo.hpp"

#if defined(CFD_HAS_MPI)
#include <mpi.h>
#endif

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfd::distributed {

#if defined(CFD_HAS_MPI)

template<class Descriptor>
class MpiSelectiveHaloExchange {
public:
    using Grid = HaloSoA3D<float, static_cast<std::size_t>(Descriptor::q)>;

    MpiSelectiveHaloExchange(const MpiCartesianRuntime& runtime, Grid& grid)
        : runtime_(runtime), grid_(grid) {
        try {
            slots_.reserve(26);
            for (const auto offset : neighbor_offsets()) {
                const std::size_t count = selective_halo_value_count<Descriptor>(grid_, offset);
                if (count == 0) continue;
                const int neighbor_rank = runtime_.neighbor(offset);
                if (neighbor_rank < 0) continue;
                Slot slot;
                slot.offset = offset;
                slot.neighbor = neighbor_rank;
                slot.send.resize(count);
                slot.recv.resize(count);
                slots_.push_back(std::move(slot));
            }
            requests_.resize(slots_.size() * 2U, MPI_REQUEST_NULL);
            std::size_t request = 0;
            for (auto& slot : slots_) {
                check_mpi(MPI_Recv_init(slot.recv.data(), mpi_count(slot.recv.size()), MPI_FLOAT,
                                        slot.neighbor, neighbor_tag(slot.offset.opposite()),
                                        runtime_.communicator(), &requests_[request++]),
                          "MPI_Recv_init(selective halo)");
            }
            for (auto& slot : slots_) {
                check_mpi(MPI_Send_init(slot.send.data(), mpi_count(slot.send.size()), MPI_FLOAT,
                                        slot.neighbor, neighbor_tag(slot.offset),
                                        runtime_.communicator(), &requests_[request++]),
                          "MPI_Send_init(selective halo)");
            }
        } catch (...) {
            release_requests();
            throw;
        }
    }

    ~MpiSelectiveHaloExchange() {
        try {
            if (active_) wait();
        } catch (...) {
            // Destructors must not throw. MPI errors will be surfaced by explicit wait().
        }
        release_requests();
    }

    MpiSelectiveHaloExchange(const MpiSelectiveHaloExchange&) = delete;
    MpiSelectiveHaloExchange& operator=(const MpiSelectiveHaloExchange&) = delete;

    void begin() {
        if (active_) throw std::logic_error("MPI selective halo exchange already active");
        for (auto& slot : slots_) {
            pack_selective_into<Descriptor>(grid_, slot.offset, slot.send.data(), slot.send.size());
        }
        if (!requests_.empty()) {
            check_mpi(MPI_Startall(mpi_count(requests_.size()), requests_.data()),
                      "MPI_Startall(selective halo)");
        }
        active_ = true;
    }

    void wait() {
        if (!active_) return;
        if (!requests_.empty()) {
            check_mpi(MPI_Waitall(mpi_count(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE),
                      "MPI_Waitall(selective halo)");
        }
        for (auto& slot : slots_) {
            unpack_selective<Descriptor>(grid_, slot.offset, slot.recv.data(), slot.recv.size());
        }
        active_ = false;
    }

    [[nodiscard]] std::size_t message_count() const noexcept { return slots_.size(); }
    [[nodiscard]] std::size_t values_per_step() const noexcept {
        std::size_t result = 0;
        for (const auto& slot : slots_) result += slot.send.size();
        return result;
    }
    [[nodiscard]] std::size_t bytes_per_step() const noexcept { return values_per_step() * sizeof(float); }

private:
    struct Slot {
        NeighborOffset offset{};
        int neighbor{-1};
        std::vector<float> send;
        std::vector<float> recv;
    };

    const MpiCartesianRuntime& runtime_;
    Grid& grid_;
    std::vector<Slot> slots_;
    std::vector<MPI_Request> requests_;
    bool active_{false};


    void release_requests() noexcept {
        for (auto& request : requests_) {
            if (request != MPI_REQUEST_NULL) MPI_Request_free(&request);
        }
    }

    static int mpi_count(std::size_t n) {
        if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("MPI selective halo message exceeds INT_MAX elements");
        }
        return static_cast<int>(n);
    }
};

#endif

} // namespace cfd::distributed
