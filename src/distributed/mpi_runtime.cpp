#include "cfd/distributed/mpi_runtime.hpp"

#if defined(CFD_HAS_MPI)

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace cfd::distributed {

void check_mpi(int error_code, const char* operation) {
    if (error_code == MPI_SUCCESS) return;
    std::array<char, MPI_MAX_ERROR_STRING> buffer{};
    int length = 0;
    MPI_Error_string(error_code, buffer.data(), &length);
    throw std::runtime_error(std::string(operation) + ": " +
                             std::string(buffer.data(), static_cast<std::size_t>(length)));
}

MpiEnvironment::MpiEnvironment(int& argc, char**& argv, int requested_thread_level) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        check_mpi(MPI_Init_thread(&argc, &argv, requested_thread_level, &provided_thread_level_), "MPI_Init_thread");
        owns_mpi_ = true;
    } else {
        check_mpi(MPI_Query_thread(&provided_thread_level_), "MPI_Query_thread");
    }
    if (provided_thread_level_ < requested_thread_level) {
        if (owns_mpi_) { MPI_Finalize(); owns_mpi_ = false; }
        throw std::runtime_error("MPI implementation does not provide requested thread level");
    }
}

MpiEnvironment::~MpiEnvironment() {
    if (!owns_mpi_) return;
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
}

MpiCartesianRuntime::MpiCartesianRuntime(Extent3 global, int dimensions,
                                         std::array<bool, 3> periodic, MPI_Comm parent)
    : global_(global), periodic_(periodic) {
    check_mpi(MPI_Comm_size(parent, &size_), "MPI_Comm_size");
    const std::size_t ranks = static_cast<std::size_t>(size_);
    grid_ = choose_process_grid(global_, ranks, dimensions);

    std::array<int, 3> dims{static_cast<int>(grid_.z), static_cast<int>(grid_.y), static_cast<int>(grid_.x)};
    std::array<int, 3> periods{periodic_[2] ? 1 : 0, periodic_[1] ? 1 : 0, periodic_[0] ? 1 : 0};
    check_mpi(MPI_Cart_create(parent, 3, dims.data(), periods.data(), 0, &cart_comm_), "MPI_Cart_create");
    if (cart_comm_ == MPI_COMM_NULL) throw std::runtime_error("MPI_Cart_create returned MPI_COMM_NULL");
    check_mpi(MPI_Comm_rank(cart_comm_, &rank_), "MPI_Comm_rank(cart)");
    std::array<int, 3> cart_coords{};
    check_mpi(MPI_Cart_coords(cart_comm_, rank_, 3, cart_coords.data()), "MPI_Cart_coords");
    const Coord3 logical_coords{cart_coords[2], cart_coords[1], cart_coords[0]};
    brick_ = make_brick(global_, grid_, logical_coords, static_cast<std::size_t>(rank_));

    check_mpi(MPI_Comm_split_type(parent, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shared_comm_), "MPI_Comm_split_type");
    check_mpi(MPI_Comm_rank(shared_comm_, &local_rank_), "MPI_Comm_rank(shared)");
    check_mpi(MPI_Comm_size(shared_comm_, &local_size_), "MPI_Comm_size(shared)");

    std::array<char, MPI_MAX_PROCESSOR_NAME> name{};
    int name_length = 0;
    check_mpi(MPI_Get_processor_name(name.data(), &name_length), "MPI_Get_processor_name");
    processor_name_.assign(name.data(), static_cast<std::size_t>(name_length));

    neighbors_.fill(-1);
    for (const auto offset : neighbor_offsets()) {
        Coord3 coord{brick_.coordinate.x + offset.x,
                     brick_.coordinate.y + offset.y,
                     brick_.coordinate.z + offset.z};
        const std::array<int, 3> logical_dims{static_cast<int>(grid_.x), static_cast<int>(grid_.y), static_cast<int>(grid_.z)};
        int* values[3]{&coord.x, &coord.y, &coord.z};
        bool valid = true;
        for (int axis = 0; axis < 3; ++axis) {
            if (*values[axis] < 0 || *values[axis] >= logical_dims[axis]) {
                if (!periodic_[static_cast<std::size_t>(axis)]) {
                    valid = false;
                    break;
                }
                int wrapped = *values[axis] % logical_dims[axis];
                if (wrapped < 0) wrapped += logical_dims[axis];
                *values[axis] = wrapped;
            }
        }
        if (!valid) continue;
        std::array<int, 3> target_coords{coord.z, coord.y, coord.x};
        int target_rank = -1;
        check_mpi(MPI_Cart_rank(cart_comm_, target_coords.data(), &target_rank), "MPI_Cart_rank");
        neighbors_[static_cast<std::size_t>(neighbor_tag(offset))] = target_rank;
    }
}

MpiCartesianRuntime::~MpiCartesianRuntime() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (finalized) return;
    if (shared_comm_ != MPI_COMM_NULL) MPI_Comm_free(&shared_comm_);
    if (cart_comm_ != MPI_COMM_NULL) MPI_Comm_free(&cart_comm_);
}

int MpiCartesianRuntime::neighbor(NeighborOffset offset) const noexcept {
    const int tag = neighbor_tag(offset);
    if (tag < 0 || tag >= static_cast<int>(neighbors_.size())) return -1;
    return neighbors_[static_cast<std::size_t>(tag)];
}

} // namespace cfd::distributed

#endif
