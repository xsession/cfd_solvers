#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <cstddef>
#include <span>
#include <string>

namespace cfd::core {

// Persistent SYCL CSR operator. Matrix storage remains on the selected device;
// host-facing helpers only copy vectors at the API boundary.
class SyclCsrLinearAlgebra {
public:
    explicit SyclCsrLinearAlgebra(const CsrMatrix& matrix,
                                  sycl::device device = sycl::device{sycl::default_selector_v});
    ~SyclCsrLinearAlgebra() noexcept;

    SyclCsrLinearAlgebra(const SyclCsrLinearAlgebra&) = delete;
    SyclCsrLinearAlgebra& operator=(const SyclCsrLinearAlgebra&) = delete;
    SyclCsrLinearAlgebra(SyclCsrLinearAlgebra&&) = delete;
    SyclCsrLinearAlgebra& operator=(SyclCsrLinearAlgebra&&) = delete;

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return nonzeros_; }
    [[nodiscard]] std::string device_name() const;

    void multiply(std::span<const double> x, std::span<double> y);

    // Device-resident matrix-free-style CG loop for SPD CSR systems. SpMV,
    // vector updates and dot-product reductions all execute as SYCL kernels;
    // only scalar convergence values and the final solution cross to the host.
    [[nodiscard]] IterativeSolverResult conjugate_gradient(
        std::span<const double> rhs,
        std::span<double> x,
        std::size_t max_iterations,
        double relative_tolerance);

private:
    sycl::queue queue_;
    std::size_t rows_{};
    std::size_t cols_{};
    std::size_t nonzeros_{};
    std::size_t* row_offsets_{nullptr};
    std::size_t* column_indices_{nullptr};
    double* values_{nullptr};
    double* reduction_value_{nullptr};

    void spmv_device(const double* x, double* y);
    [[nodiscard]] double dot_device(const double* a, const double* b, std::size_t n);
    void release() noexcept;
};

} // namespace cfd::core
#endif
