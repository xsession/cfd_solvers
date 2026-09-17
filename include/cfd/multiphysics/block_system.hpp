#pragma once
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::multiphysics {
class BlockCoupledLinearSystem {
public:
    explicit BlockCoupledLinearSystem(std::vector<std::size_t> block_sizes);
    void add(std::size_t row_block,std::size_t row,std::size_t col_block,std::size_t col,double value);
    void add_rhs(std::size_t block,std::size_t row,double value);
    [[nodiscard]] std::size_t blocks() const noexcept { return sizes_.size(); }
    [[nodiscard]] std::size_t total_size() const noexcept { return offsets_.back(); }
    [[nodiscard]] std::span<const double> rhs() const noexcept { return rhs_; }
    [[nodiscard]] cfd::core::CsrMatrix matrix() const;
    [[nodiscard]] cfd::core::IterativeSolverResult solve(std::span<double> solution,
        std::size_t max_iterations=1000,std::size_t restart=30,double relative_tolerance=1e-10) const;
private:
    std::vector<std::size_t> sizes_,offsets_;
    cfd::core::CsrBuilder builder_;
    std::vector<double> rhs_;
    [[nodiscard]] std::size_t index(std::size_t block,std::size_t local) const;
};
} // namespace cfd::multiphysics
