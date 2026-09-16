#pragma once

#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/halo_grid.hpp"

#if defined(CFD_HAS_MPI)
#include <mpi.h>
#endif

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace cfd::distributed {

#if defined(CFD_HAS_MPI)

void check_mpi(int error_code, const char* operation);

class MpiEnvironment {
public:
    MpiEnvironment(int& argc, char**& argv, int requested_thread_level = MPI_THREAD_FUNNELED);
    ~MpiEnvironment();

    MpiEnvironment(const MpiEnvironment&) = delete;
    MpiEnvironment& operator=(const MpiEnvironment&) = delete;

    [[nodiscard]] int provided_thread_level() const noexcept { return provided_thread_level_; }

private:
    bool owns_mpi_{false};
    int provided_thread_level_{MPI_THREAD_SINGLE};
};

class MpiCartesianRuntime {
public:
    MpiCartesianRuntime(Extent3 global,
                        int dimensions = 3,
                        std::array<bool, 3> periodic = {true, true, true},
                        MPI_Comm parent = MPI_COMM_WORLD);
    ~MpiCartesianRuntime();

    MpiCartesianRuntime(const MpiCartesianRuntime&) = delete;
    MpiCartesianRuntime& operator=(const MpiCartesianRuntime&) = delete;

    [[nodiscard]] int rank() const noexcept { return rank_; }
    [[nodiscard]] int size() const noexcept { return size_; }
    [[nodiscard]] int local_rank() const noexcept { return local_rank_; }
    [[nodiscard]] int local_size() const noexcept { return local_size_; }
    [[nodiscard]] const std::string& processor_name() const noexcept { return processor_name_; }
    [[nodiscard]] const Brick& brick() const noexcept { return brick_; }
    [[nodiscard]] const ProcessGrid& process_grid() const noexcept { return grid_; }
    [[nodiscard]] const std::array<bool, 3>& periodic() const noexcept { return periodic_; }
    [[nodiscard]] MPI_Comm communicator() const noexcept { return cart_comm_; }

    [[nodiscard]] int neighbor(NeighborOffset offset) const noexcept;

private:
    MPI_Comm cart_comm_{MPI_COMM_NULL};
    MPI_Comm shared_comm_{MPI_COMM_NULL};
    int rank_{};
    int size_{1};
    int local_rank_{};
    int local_size_{1};
    std::string processor_name_;
    Extent3 global_{};
    ProcessGrid grid_{};
    Brick brick_{};
    std::array<bool, 3> periodic_{true, true, true};
    std::array<int, 27> neighbors_{};
};

template<class T, std::size_t Components>
class MpiHaloExchange {
public:
    MpiHaloExchange(const MpiCartesianRuntime& runtime, HaloSoA3D<T, Components>& grid)
        : runtime_(runtime), grid_(grid) {
        static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>, "MPI halo baseline supports float/double");
        slots_.reserve(26);
        for (const auto offset : neighbor_offsets()) {
            const int neighbor_rank = runtime_.neighbor(offset);
            if (neighbor_rank < 0) continue;
            Slot slot;
            slot.offset = offset;
            slot.neighbor = neighbor_rank;
            slot.send = grid_.pack(offset);
            slot.recv.resize(grid_.halo_value_count(offset));
            slots_.push_back(std::move(slot));
        }
    }

    void begin() {
        if (active_) throw std::logic_error("MPI halo exchange already active");
        requests_.clear();
        requests_.reserve(slots_.size() * 2U);
        for (auto& slot : slots_) {
            slot.send = grid_.pack(slot.offset);
            requests_.push_back(MPI_REQUEST_NULL);
            check_mpi(MPI_Irecv(slot.recv.data(), mpi_count(slot.recv.size()), mpi_type(), slot.neighbor,
                                neighbor_tag(slot.offset.opposite()), runtime_.communicator(), &requests_.back()),
                      "MPI_Irecv");
        }
        for (auto& slot : slots_) {
            requests_.push_back(MPI_REQUEST_NULL);
            check_mpi(MPI_Isend(slot.send.data(), mpi_count(slot.send.size()), mpi_type(), slot.neighbor,
                                neighbor_tag(slot.offset), runtime_.communicator(), &requests_.back()),
                      "MPI_Isend");
        }
        active_ = true;
    }

    void wait() {
        if (!active_) return;
        if (!requests_.empty()) check_mpi(MPI_Waitall(mpi_count(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE),
                                          "MPI_Waitall");
        for (auto& slot : slots_) grid_.unpack(slot.offset, slot.recv);
        active_ = false;
    }

    ~MpiHaloExchange() {
        if (active_) wait();
    }

private:
    struct Slot {
        NeighborOffset offset{};
        int neighbor{-1};
        std::vector<T> send;
        std::vector<T> recv;
    };

    const MpiCartesianRuntime& runtime_;
    HaloSoA3D<T, Components>& grid_;
    std::vector<Slot> slots_;
    std::vector<MPI_Request> requests_;
    bool active_{false};

    static int mpi_count(std::size_t n) {
        if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("MPI baseline halo message exceeds INT_MAX elements");
        }
        return static_cast<int>(n);
    }

    static MPI_Datatype mpi_type() noexcept {
        if constexpr (std::is_same_v<T, float>) return MPI_FLOAT;
        else return MPI_DOUBLE;
    }
};

#endif

} // namespace cfd::distributed
