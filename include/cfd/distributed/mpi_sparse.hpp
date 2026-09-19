#pragma once

#include "cfd/core/distributed_csr.hpp"
#include "cfd/core/iterative_solvers.hpp"

#if defined(CFD_HAS_MPI)
#include <mpi.h>

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::distributed {

// Ownership-aware distributed CSR operator for contiguous row partitions.
// Setup exchanges requested global column IDs once; each multiply then moves
// only the vector values referenced by remote columns before local SpMV.
class MpiDistributedCsrOperator {
public:
    MpiDistributedCsrOperator(const cfd::core::CsrMatrix& global_matrix, std::size_t row_begin, std::size_t row_end,
                              MPI_Comm communicator = MPI_COMM_WORLD);
    ~MpiDistributedCsrOperator();

    MpiDistributedCsrOperator(const MpiDistributedCsrOperator&) = delete;
    MpiDistributedCsrOperator& operator=(const MpiDistributedCsrOperator&) = delete;

    [[nodiscard]] std::size_t local_rows() const noexcept { return partition_.local_rows(); }
    [[nodiscard]] std::size_t global_rows() const noexcept { return global_rows_; }
    [[nodiscard]] std::size_t row_begin() const noexcept { return partition_.row_begin(); }
    [[nodiscard]] MPI_Comm communicator() const noexcept { return communicator_; }
    [[nodiscard]] const cfd::core::DistributedCsrPartition& partition() const noexcept { return partition_; }
    [[nodiscard]] std::size_t halo_value_count() const noexcept { return partition_.halo_columns().size(); }
    [[nodiscard]] std::size_t halo_message_count() const noexcept;
    [[nodiscard]] std::size_t halo_bytes_per_step() const noexcept;
    [[nodiscard]] bool multiply_pending() const noexcept { return multiply_pending_; }

    // Start only the communication phase. The caller may execute independent
    // local work before finish_multiply() waits, unpacks the halo, and runs
    // the local sparse operator.
    void begin_multiply(std::span<const double> local_x);
    void finish_multiply(std::span<const double> local_x, std::span<double> local_y);
    void multiply(std::span<const double> local_x, std::span<double> local_y);
    [[nodiscard]] double global_sum(double local_value) const;

private:
    cfd::core::DistributedCsrPartition partition_;
    MPI_Comm communicator_{MPI_COMM_NULL};
    int rank_{};
    int size_{1};
    std::size_t global_rows_{};
    std::vector<std::size_t> ownership_begin_;
    std::vector<std::size_t> ownership_end_;

    // request_* describes halo columns this rank asks from each owner.
    std::vector<int> request_counts_;
    std::vector<int> request_displacements_;
    // response_* describes local values requested by peers.
    std::vector<int> response_counts_;
    std::vector<int> response_displacements_;
    std::vector<std::size_t> response_local_offsets_;
    std::vector<std::size_t> receive_to_halo_slot_;
    std::vector<double> send_values_;
    std::vector<double> receive_values_;
    std::vector<double> halo_values_;
    MPI_Request multiply_request_{MPI_REQUEST_NULL};
    bool multiply_pending_{false};

    [[nodiscard]] int owner_of(std::size_t global_index) const;
};

[[nodiscard]] cfd::core::IterativeSolverResult
mpi_distributed_conjugate_gradient(MpiDistributedCsrOperator& matrix, std::span<const double> local_rhs,
                                   std::span<double> local_x, cfd::core::KrylovWorkspace& workspace,
                                   std::size_t max_iterations, double relative_tolerance);

} // namespace cfd::distributed
#endif
