#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/device_residency.hpp"
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
    // Reuse an existing queue/context so resident USM fields can be passed
    // directly into SpMV/CG without a context-crossing copy. The caller should
    // provide an in-order queue because the Krylov kernels rely on submission order.
    SyclCsrLinearAlgebra(const CsrMatrix& matrix, sycl::queue queue);
    ~SyclCsrLinearAlgebra() noexcept;

    SyclCsrLinearAlgebra(const SyclCsrLinearAlgebra&) = delete;
    SyclCsrLinearAlgebra& operator=(const SyclCsrLinearAlgebra&) = delete;
    SyclCsrLinearAlgebra(SyclCsrLinearAlgebra&&) = delete;
    SyclCsrLinearAlgebra& operator=(SyclCsrLinearAlgebra&&) = delete;

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] std::size_t nonzeros() const noexcept { return nonzeros_; }
    [[nodiscard]] std::string device_name() const;
    [[nodiscard]] sycl::queue& queue() noexcept { return queue_; }
    [[nodiscard]] const sycl::queue& queue() const noexcept { return queue_; }

    void multiply(std::span<const double> x, std::span<double> y);

    // Direct device-USM seam for resident solver fields. The caller owns the
    // allocations and must ensure they are accessible from this queue/context.
    // multiply_device() is synchronous at the API boundary but performs no
    // host field transfer.
    void multiply_device(const double* x_device, double* y_device);

    // Replace CSR numeric values without changing sparsity. This is the key
    // seam for resident nonlinear/FVM operators whose coefficients evolve
    // while connectivity stays fixed. The device variant performs only a
    // device-to-device copy and preserves row/column storage.
    void update_values(std::span<const double> values);
    void update_values_device(const double* values_device);

    [[nodiscard]] const std::size_t* row_offsets_device() const noexcept { return row_offsets_; }
    [[nodiscard]] const std::size_t* column_indices_device() const noexcept { return column_indices_; }

    // Device-resident matrix-free-style CG loop for SPD CSR systems. SpMV,
    // vector updates and dot-product reductions all execute as SYCL kernels.
    // This host-span compatibility entry point stages rhs/x at the API boundary.
    [[nodiscard]] IterativeSolverResult conjugate_gradient(
        std::span<const double> rhs,
        std::span<double> x,
        std::size_t max_iterations,
        double relative_tolerance);

    // Direct device-USM CG variant. rhs_device is read-only and x_device is
    // updated in place. Residual reductions currently synchronize a shared
    // scalar with the host for convergence decisions, but no volume vector is
    // staged through host memory.
    [[nodiscard]] IterativeSolverResult conjugate_gradient_device(
        const double* rhs_device,
        double* x_device,
        std::size_t max_iterations,
        double relative_tolerance);

    // Nonsymmetric resident Krylov solve for upwind FVM momentum/scalar
    // operators. Like the CG device seam, only scalar convergence reductions
    // synchronize with the host; no volume vector is staged through host RAM.
    [[nodiscard]] IterativeSolverResult bicgstab(
        std::span<const double> rhs,
        std::span<double> x,
        std::size_t max_iterations,
        double relative_tolerance);
    [[nodiscard]] IterativeSolverResult bicgstab_device(
        const double* rhs_device,
        double* x_device,
        std::size_t max_iterations,
        double relative_tolerance);

    [[nodiscard]] const DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

private:
    sycl::queue queue_;
    std::size_t rows_{};
    std::size_t cols_{};
    std::size_t nonzeros_{};
    std::size_t* row_offsets_{nullptr};
    std::size_t* column_indices_{nullptr};
    double* values_{nullptr};
    double* reduction_value_{nullptr};
    // Persistent Krylov/vector workspace: no per-solve device allocation churn.
    double* work_rhs_{nullptr};
    double* work_x_{nullptr};
    double* work_r_{nullptr};
    double* work_p_{nullptr};
    double* work_q_{nullptr};
    double* work_rhat_{nullptr};
    double* work_s_{nullptr};
    double* work_t_{nullptr};
    mutable DeviceTransferStats transfer_stats_{};

    void spmv_device(const double* x, double* y);
    [[nodiscard]] double dot_device(const double* a, const double* b, std::size_t n);
    void release() noexcept;
};

} // namespace cfd::core
#endif
