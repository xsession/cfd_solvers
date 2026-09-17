#pragma once
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::core {
struct RefinementResult {bool converged{};std::size_t outer_iterations{};double relative_residual{};};
[[nodiscard]] RefinementResult mixed_precision_iterative_refinement(const CsrMatrix& A,std::span<const double> b,std::span<double> x,std::size_t outer_iterations=12,std::size_t inner_jacobi_iterations=100,double relative_tolerance=1e-11);
} // namespace cfd::core
