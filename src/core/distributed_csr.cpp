#include "cfd/core/distributed_csr.hpp"
#include <algorithm>
#include <stdexcept>
namespace cfd::core {
DistributedCsrPartition::DistributedCsrPartition(CsrMatrix A, std::size_t begin, std::size_t end)
    : matrix_(std::move(A)), row_begin_(begin), row_end_(end) {
    if (begin >= end || end > matrix_.rows() || matrix_.rows() != matrix_.cols()) {
        throw std::invalid_argument("invalid distributed CSR row partition");
    }

    const auto& offsets = matrix_.row_offsets();
    const auto& columns = matrix_.column_indices();
    column_refs_.reserve(offsets[end] - offsets[begin]);
    for (std::size_t r = begin; r < end; ++r) {
        for (std::size_t k = offsets[r]; k < offsets[r + 1U]; ++k) {
            const auto column = columns[k];
            if (column < begin || column >= end)
                halo_columns_.push_back(column);
        }
    }
    std::sort(halo_columns_.begin(), halo_columns_.end());
    halo_columns_.erase(std::unique(halo_columns_.begin(), halo_columns_.end()), halo_columns_.end());

    const auto local_count = local_rows();
    for (std::size_t r = begin; r < end; ++r) {
        for (std::size_t k = offsets[r]; k < offsets[r + 1U]; ++k) {
            const auto column = columns[k];
            if (column >= begin && column < end) {
                column_refs_.push_back(column - begin);
                continue;
            }
            const auto it = std::lower_bound(halo_columns_.begin(), halo_columns_.end(), column);
            if (it == halo_columns_.end() || *it != column) {
                throw std::runtime_error("distributed CSR missing halo column during setup");
            }
            column_refs_.push_back(local_count + static_cast<std::size_t>(it - halo_columns_.begin()));
        }
    }
}

void DistributedCsrPartition::multiply(std::span<const double> local, std::span<const double> halo,
                                       std::span<double> y) const {
    if (local.size() != local_rows() || halo.size() != halo_columns_.size() || y.size() != local_rows()) {
        throw std::invalid_argument("distributed CSR buffer mismatch");
    }

    const auto& offsets = matrix_.row_offsets();
    const auto& values = matrix_.values();
    const auto local_count = local_rows();
    parallel_for(local_count, [&](std::size_t local_row) {
        const std::size_t row = row_begin_ + local_row;
        const std::size_t ref_begin = offsets[row] - offsets[row_begin_];
        const std::size_t ref_end = offsets[row + 1U] - offsets[row_begin_];
        double sum = 0.0;
        for (std::size_t ref = ref_begin; ref < ref_end; ++ref) {
            const auto column_ref = column_refs_[ref];
            const double x = column_ref < local_count ? local[column_ref] : halo[column_ref - local_count];
            sum += values[offsets[row] + (ref - ref_begin)] * x;
        }
        y[local_row] = sum;
    });
}
} // namespace cfd::core
