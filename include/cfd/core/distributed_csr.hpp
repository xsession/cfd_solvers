#pragma once
#include "cfd/core/csr_matrix.hpp"
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::core {
struct HaloEntry {
    std::size_t global_column{};
    std::size_t halo_slot{};
};
class DistributedCsrPartition {
public:
    DistributedCsrPartition(CsrMatrix global, std::size_t row_begin, std::size_t row_end);
    [[nodiscard]] std::size_t local_rows() const noexcept { return row_end_ - row_begin_; }
    [[nodiscard]] std::size_t row_begin() const noexcept { return row_begin_; }
    [[nodiscard]] const std::vector<std::size_t>& halo_columns() const noexcept { return halo_columns_; }
    void multiply(std::span<const double> local_x, std::span<const double> halo_x, std::span<double> y) const;

private:
    CsrMatrix matrix_;
    std::size_t row_begin_{}, row_end_{};
    std::vector<std::size_t> halo_columns_;
    // One compact reference per stored coefficient. Local references are
    // [0, local_rows); halo references are local_rows + halo_slot. The
    // mapping is built once so SpMV never performs a lower_bound in its hot
    // inner loop.
    std::vector<std::size_t> column_refs_;
};
} // namespace cfd::core
