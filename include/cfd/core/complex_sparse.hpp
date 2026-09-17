#pragma once

#include "cfd/core/iterative_solvers.hpp"

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::core {

using Complex = std::complex<double>;

class ComplexCsrBuilder {
public:
    explicit ComplexCsrBuilder(std::size_t dimension);
    void add(std::size_t row,std::size_t column,Complex value);
    [[nodiscard]] std::size_t dimension() const noexcept { return dimension_; }

    // Builds the equivalent real block matrix [Re(A) -Im(A); Im(A) Re(A)].
    [[nodiscard]] CsrMatrix build_real_block(double drop_tolerance=0.0) const;
private:
    std::size_t dimension_{};
    std::vector<std::vector<std::pair<std::size_t,Complex>>> entries_;
};

struct ComplexSparseSolveConfig {
    std::size_t max_iterations{4000U};
    std::size_t gmres_restart{60U};
    double relative_tolerance{1.0e-10};
    double ilu_diagonal_floor{1.0e-24};
};

struct ComplexSparseSolveResult {
    std::vector<Complex> solution;
    IterativeSolverResult linear_result;
};

[[nodiscard]] ComplexSparseSolveResult solve_complex_sparse(
    const ComplexCsrBuilder& matrix,
    std::span<const Complex> rhs,
    const ComplexSparseSolveConfig& config={});

} // namespace cfd::core
