#include "cfd/distributed/mpi_sparse.hpp"

#if defined(CFD_HAS_MPI)

#include "cfd/core/parallel.hpp"
#include "cfd/distributed/mpi_runtime.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cfd::distributed {
namespace {

[[nodiscard]] int mpi_count(std::size_t n) {
    if (n > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("distributed sparse exchange exceeds MPI int count range");
    }
    return static_cast<int>(n);
}

[[nodiscard]] std::vector<int> displacements(const std::vector<int>& counts) {
    std::vector<int> result(counts.size(), 0);
    long long offset = 0;
    for (std::size_t i = 0U; i < counts.size(); ++i) {
        if (offset > static_cast<long long>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("distributed sparse displacement exceeds MPI int range");
        }
        result[i] = static_cast<int>(offset);
        offset += counts[i];
    }
    if (offset > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("distributed sparse payload exceeds MPI int range");
    }
    return result;
}

[[nodiscard]] std::size_t total_count(const std::vector<int>& counts) {
    std::size_t result = 0U;
    for (int count : counts) {
        if (count < 0) throw std::runtime_error("negative MPI exchange count");
        result += static_cast<std::size_t>(count);
    }
    return result;
}

} // namespace

MpiDistributedCsrOperator::MpiDistributedCsrOperator(const cfd::core::CsrMatrix& global_matrix,
                                                     std::size_t row_begin,
                                                     std::size_t row_end,
                                                     MPI_Comm communicator)
    : partition_(global_matrix, row_begin, row_end), global_rows_(global_matrix.rows()) {
    if (global_matrix.rows() != global_matrix.cols()) throw std::invalid_argument("distributed Krylov requires square CSR");
    check_mpi(MPI_Comm_dup(communicator, &communicator_), "MPI_Comm_dup(sparse)");
    try {
        check_mpi(MPI_Comm_rank(communicator_, &rank_), "MPI_Comm_rank(sparse)");
        check_mpi(MPI_Comm_size(communicator_, &size_), "MPI_Comm_size(sparse)");

        const unsigned long long local_range[2]{
            static_cast<unsigned long long>(row_begin),
            static_cast<unsigned long long>(row_end)};
        std::vector<unsigned long long> ranges(static_cast<std::size_t>(size_) * 2U);
        check_mpi(MPI_Allgather(local_range, 2, MPI_UNSIGNED_LONG_LONG,
                                ranges.data(), 2, MPI_UNSIGNED_LONG_LONG, communicator_),
                  "MPI_Allgather(sparse ownership)");
        ownership_begin_.resize(static_cast<std::size_t>(size_));
        ownership_end_.resize(static_cast<std::size_t>(size_));
        for (int r = 0; r < size_; ++r) {
            ownership_begin_[static_cast<std::size_t>(r)] = static_cast<std::size_t>(ranges[2U * static_cast<std::size_t>(r)]);
            ownership_end_[static_cast<std::size_t>(r)] = static_cast<std::size_t>(ranges[2U * static_cast<std::size_t>(r) + 1U]);
        }
        if (ownership_begin_.front() != 0U || ownership_end_.back() != global_rows_) {
            throw std::invalid_argument("distributed CSR ownership must cover all rows");
        }
        for (std::size_t r = 0U; r < ownership_begin_.size(); ++r) {
            if (ownership_begin_[r] >= ownership_end_[r]) throw std::invalid_argument("empty distributed CSR ownership range");
            if (r > 0U && ownership_begin_[r] != ownership_end_[r - 1U]) {
                throw std::invalid_argument("distributed CSR ownership must be contiguous and non-overlapping");
            }
        }

        request_counts_.assign(static_cast<std::size_t>(size_), 0);
        std::vector<std::vector<unsigned long long>> requested_by_owner(static_cast<std::size_t>(size_));
        std::vector<std::vector<std::size_t>> slots_by_owner(static_cast<std::size_t>(size_));
        const auto& halo_columns = partition_.halo_columns();
        for (std::size_t slot = 0U; slot < halo_columns.size(); ++slot) {
            const auto column = halo_columns[slot];
            const int owner = owner_of(column);
            if (owner == rank_) throw std::logic_error("local distributed CSR column classified as halo");
            requested_by_owner[static_cast<std::size_t>(owner)].push_back(static_cast<unsigned long long>(column));
            slots_by_owner[static_cast<std::size_t>(owner)].push_back(slot);
        }
        for (int r = 0; r < size_; ++r) {
            request_counts_[static_cast<std::size_t>(r)] = mpi_count(requested_by_owner[static_cast<std::size_t>(r)].size());
        }
        request_displacements_ = displacements(request_counts_);

        response_counts_.assign(static_cast<std::size_t>(size_), 0);
        check_mpi(MPI_Alltoall(request_counts_.data(), 1, MPI_INT,
                               response_counts_.data(), 1, MPI_INT, communicator_),
                  "MPI_Alltoall(sparse request counts)");
        response_displacements_ = displacements(response_counts_);

        std::vector<unsigned long long> outgoing_requests(total_count(request_counts_));
        receive_to_halo_slot_.resize(outgoing_requests.size());
        for (int r = 0; r < size_; ++r) {
            const auto rr = static_cast<std::size_t>(r);
            const auto base = static_cast<std::size_t>(request_displacements_[rr]);
            std::copy(requested_by_owner[rr].begin(), requested_by_owner[rr].end(), outgoing_requests.begin() + static_cast<std::ptrdiff_t>(base));
            std::copy(slots_by_owner[rr].begin(), slots_by_owner[rr].end(), receive_to_halo_slot_.begin() + static_cast<std::ptrdiff_t>(base));
        }

        std::vector<unsigned long long> incoming_requests(total_count(response_counts_));
        check_mpi(MPI_Alltoallv(outgoing_requests.data(), request_counts_.data(), request_displacements_.data(), MPI_UNSIGNED_LONG_LONG,
                                incoming_requests.data(), response_counts_.data(), response_displacements_.data(), MPI_UNSIGNED_LONG_LONG,
                                communicator_),
                  "MPI_Alltoallv(sparse column requests)");

        response_local_offsets_.resize(incoming_requests.size());
        for (std::size_t i = 0U; i < incoming_requests.size(); ++i) {
            const auto global_column = static_cast<std::size_t>(incoming_requests[i]);
            if (global_column < row_begin || global_column >= row_end) {
                throw std::runtime_error("peer requested sparse vector entry not owned by this rank");
            }
            response_local_offsets_[i] = global_column - row_begin;
        }

        send_values_.resize(response_local_offsets_.size());
        receive_values_.resize(outgoing_requests.size());
        halo_values_.resize(halo_columns.size());
    } catch (...) {
        if (communicator_ != MPI_COMM_NULL) MPI_Comm_free(&communicator_);
        communicator_ = MPI_COMM_NULL;
        throw;
    }
}

MpiDistributedCsrOperator::~MpiDistributedCsrOperator() {
    if (communicator_ == MPI_COMM_NULL) return;
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Comm_free(&communicator_);
}

int MpiDistributedCsrOperator::owner_of(std::size_t global_index) const {
    const auto it = std::upper_bound(ownership_end_.begin(), ownership_end_.end(), global_index);
    if (it == ownership_end_.end()) throw std::out_of_range("global sparse column has no owner");
    const auto owner = static_cast<std::size_t>(it - ownership_end_.begin());
    if (global_index < ownership_begin_[owner] || global_index >= ownership_end_[owner]) {
        throw std::out_of_range("global sparse column lies outside ownership partition");
    }
    return static_cast<int>(owner);
}

void MpiDistributedCsrOperator::multiply(std::span<const double> local_x,
                                         std::span<double> local_y) {
    if (local_x.size() != local_rows() || local_y.size() != local_rows()) {
        throw std::invalid_argument("distributed sparse multiply local vector size mismatch");
    }
    for (std::size_t i = 0U; i < response_local_offsets_.size(); ++i) {
        send_values_[i] = local_x[response_local_offsets_[i]];
    }

    MPI_Request request = MPI_REQUEST_NULL;
    check_mpi(MPI_Ialltoallv(send_values_.data(), response_counts_.data(), response_displacements_.data(), MPI_DOUBLE,
                             receive_values_.data(), request_counts_.data(), request_displacements_.data(), MPI_DOUBLE,
                             communicator_, &request),
              "MPI_Ialltoallv(sparse vector)");
    check_mpi(MPI_Wait(&request, MPI_STATUS_IGNORE), "MPI_Wait(sparse vector)");

    for (std::size_t i = 0U; i < receive_values_.size(); ++i) {
        halo_values_[receive_to_halo_slot_[i]] = receive_values_[i];
    }
    partition_.multiply(local_x, halo_values_, local_y);
}

double MpiDistributedCsrOperator::global_sum(double local_value) const {
    double global = 0.0;
    check_mpi(MPI_Allreduce(&local_value, &global, 1, MPI_DOUBLE, MPI_SUM, communicator_),
              "MPI_Allreduce(distributed Krylov)");
    return global;
}

cfd::core::IterativeSolverResult mpi_distributed_conjugate_gradient(
    MpiDistributedCsrOperator& matrix,
    std::span<const double> local_rhs,
    std::span<double> local_x,
    cfd::core::KrylovWorkspace& workspace,
    std::size_t max_iterations,
    double relative_tolerance) {
    if (local_rhs.empty() || local_rhs.size() != matrix.local_rows() || local_x.size() != matrix.local_rows()) {
        throw std::invalid_argument("distributed CG local vector size mismatch");
    }
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid distributed CG controls");
    }
    const std::size_t n = local_rhs.size();
    workspace.resize(n);
    auto& r = workspace.r;
    auto& p = workspace.p;
    auto& q = workspace.q;

    matrix.multiply(local_x, std::span<double>(q.data(), n));
    cfd::core::parallel_for(n, [&](std::size_t i) {
        r[i] = local_rhs[i] - q[i];
        p[i] = r[i];
    });
    const double local_rhs2 = cfd::core::parallel_sum(n, [&](std::size_t i) { return local_rhs[i] * local_rhs[i]; });
    const double rhs2 = matrix.global_sum(local_rhs2);
    double rr = matrix.global_sum(cfd::core::parallel_sum(n, [&](std::size_t i) { return r[i] * r[i]; }));
    const double scale = std::max(std::sqrt(rhs2 / static_cast<double>(matrix.global_rows())), 1.0);
    const double target = relative_tolerance * scale;
    double rms = std::sqrt(rr / static_cast<double>(matrix.global_rows()));
    if (rms <= target) return {0U, rms, true};

    for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
        matrix.multiply(std::span<const double>(p.data(), n), std::span<double>(q.data(), n));
        const double local_pq = cfd::core::parallel_sum(n, [&](std::size_t i) { return p[i] * q[i]; });
        const double pq = matrix.global_sum(local_pq);
        if (!(pq > 0.0) || !std::isfinite(pq)) return {iteration, rms, false};
        const double alpha = rr / pq;
        cfd::core::parallel_for(n, [&](std::size_t i) {
            local_x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        });
        const double rr_new = matrix.global_sum(cfd::core::parallel_sum(n, [&](std::size_t i) { return r[i] * r[i]; }));
        rms = std::sqrt(rr_new / static_cast<double>(matrix.global_rows()));
        const std::size_t completed = iteration + 1U;
        if (rms <= target) return {completed, rms, true};
        if (!std::isfinite(rr_new) || rr == 0.0) return {completed, rms, false};
        const double beta = rr_new / rr;
        cfd::core::parallel_for(n, [&](std::size_t i) { p[i] = r[i] + beta * p[i]; });
        rr = rr_new;
    }
    return {max_iterations, rms, false};
}

} // namespace cfd::distributed

#endif
