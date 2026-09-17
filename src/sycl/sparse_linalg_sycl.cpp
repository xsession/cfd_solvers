#include "cfd/core/sycl_sparse.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace cfd::core {
namespace {

template<class T>
void free_if(T*& ptr, sycl::queue& queue) noexcept {
    if (!ptr) return;
    try { sycl::free(ptr, queue); } catch (...) {}
    ptr = nullptr;
}

} // namespace

SyclCsrLinearAlgebra::SyclCsrLinearAlgebra(const CsrMatrix& matrix, sycl::device device)
    : queue_(device, sycl::property::queue::in_order{}),
      rows_(matrix.rows()), cols_(matrix.cols()), nonzeros_(matrix.nonzeros()) {
    if (rows_ == 0U || cols_ == 0U) throw std::invalid_argument("SYCL CSR matrix must be non-empty");
    if (rows_ != cols_) throw std::invalid_argument("SYCL Krylov baseline requires a square CSR matrix");
    if (!queue_.get_device().has(sycl::aspect::usm_device_allocations) ||
        !queue_.get_device().has(sycl::aspect::usm_shared_allocations)) {
        throw std::runtime_error("selected SYCL device lacks required USM capabilities");
    }

    row_offsets_ = sycl::malloc_device<std::size_t>(rows_ + 1U, queue_);
    column_indices_ = sycl::malloc_device<std::size_t>(nonzeros_, queue_);
    values_ = sycl::malloc_device<double>(nonzeros_, queue_);
    reduction_value_ = sycl::malloc_shared<double>(1U, queue_);
    if (!row_offsets_ || !column_indices_ || !values_ || !reduction_value_) {
        release();
        throw std::bad_alloc{};
    }

    try {
        queue_.memcpy(row_offsets_, matrix.row_offsets().data(), (rows_ + 1U) * sizeof(std::size_t));
        queue_.memcpy(column_indices_, matrix.column_indices().data(), nonzeros_ * sizeof(std::size_t));
        queue_.memcpy(values_, matrix.values().data(), nonzeros_ * sizeof(double));
        queue_.wait_and_throw();
    } catch (...) {
        release();
        throw;
    }
}

SyclCsrLinearAlgebra::~SyclCsrLinearAlgebra() noexcept { release(); }

void SyclCsrLinearAlgebra::release() noexcept {
    try { queue_.wait_and_throw(); } catch (...) {}
    free_if(row_offsets_, queue_);
    free_if(column_indices_, queue_);
    free_if(values_, queue_);
    free_if(reduction_value_, queue_);
}

std::string SyclCsrLinearAlgebra::device_name() const {
    return queue_.get_device().get_info<sycl::info::device::name>();
}

void SyclCsrLinearAlgebra::spmv_device(const double* x, double* y) {
    const auto* offsets = row_offsets_;
    const auto* columns = column_indices_;
    const auto* values = values_;
    queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
        const std::size_t row = id[0];
        double sum = 0.0;
        for (std::size_t k = offsets[row]; k < offsets[row + 1U]; ++k) {
            sum += values[k] * x[columns[k]];
        }
        y[row] = sum;
    });
}

double SyclCsrLinearAlgebra::dot_device(const double* a, const double* b, std::size_t n) {
    *reduction_value_ = 0.0;
    auto reduction = sycl::reduction(reduction_value_, sycl::plus<double>());
    queue_.parallel_for(sycl::range<1>(n), reduction,
                        [=](sycl::id<1> id, auto& sum) {
                            const std::size_t i = id[0];
                            sum.combine(a[i] * b[i]);
                        }).wait_and_throw();
    return *reduction_value_;
}

void SyclCsrLinearAlgebra::multiply(std::span<const double> x, std::span<double> y) {
    if (x.size() != cols_ || y.size() != rows_) throw std::invalid_argument("SYCL CSR multiply size mismatch");
    double* dx = sycl::malloc_device<double>(cols_, queue_);
    double* dy = sycl::malloc_device<double>(rows_, queue_);
    if (!dx || !dy) {
        free_if(dx, queue_);
        free_if(dy, queue_);
        throw std::bad_alloc{};
    }
    try {
        queue_.memcpy(dx, x.data(), cols_ * sizeof(double));
        spmv_device(dx, dy);
        queue_.memcpy(y.data(), dy, rows_ * sizeof(double)).wait_and_throw();
    } catch (...) {
        free_if(dx, queue_);
        free_if(dy, queue_);
        throw;
    }
    free_if(dx, queue_);
    free_if(dy, queue_);
}

IterativeSolverResult SyclCsrLinearAlgebra::conjugate_gradient(
    std::span<const double> rhs,
    std::span<double> x,
    std::size_t max_iterations,
    double relative_tolerance) {
    if (rhs.empty() || rhs.size() != rows_ || x.size() != rows_) {
        throw std::invalid_argument("SYCL CG vector size mismatch");
    }
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid SYCL CG controls");
    }

    double* drhs = sycl::malloc_device<double>(rows_, queue_);
    double* dx = sycl::malloc_device<double>(rows_, queue_);
    double* dr = sycl::malloc_device<double>(rows_, queue_);
    double* dp = sycl::malloc_device<double>(rows_, queue_);
    double* dq = sycl::malloc_device<double>(rows_, queue_);
    auto cleanup = [&]() noexcept {
        free_if(drhs, queue_); free_if(dx, queue_); free_if(dr, queue_); free_if(dp, queue_); free_if(dq, queue_);
    };
    if (!drhs || !dx || !dr || !dp || !dq) {
        cleanup();
        throw std::bad_alloc{};
    }

    try {
        queue_.memcpy(drhs, rhs.data(), rows_ * sizeof(double));
        queue_.memcpy(dx, x.data(), rows_ * sizeof(double));
        spmv_device(dx, dq);
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            dr[i] = drhs[i] - dq[i];
            dp[i] = dr[i];
        });
        queue_.wait_and_throw();

        const double rhs2 = dot_device(drhs, drhs, rows_);
        double rr = dot_device(dr, dr, rows_);
        const double scale = std::max(std::sqrt(rhs2 / static_cast<double>(rows_)), 1.0);
        const double target = relative_tolerance * scale;
        double rms = std::sqrt(rr / static_cast<double>(rows_));
        if (rms <= target) {
            queue_.memcpy(x.data(), dx, rows_ * sizeof(double)).wait_and_throw();
            cleanup();
            return {0U, rms, true};
        }

        for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
            spmv_device(dp, dq);
            queue_.wait_and_throw();
            const double dq_dot = dot_device(dp, dq, rows_);
            if (!(dq_dot > 0.0) || !std::isfinite(dq_dot)) {
                queue_.memcpy(x.data(), dx, rows_ * sizeof(double)).wait_and_throw();
                cleanup();
                return {iteration, rms, false};
            }
            const double alpha = rr / dq_dot;
            queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
                const std::size_t i = id[0];
                dx[i] += alpha * dp[i];
                dr[i] -= alpha * dq[i];
            });
            queue_.wait_and_throw();
            const double rr_new = dot_device(dr, dr, rows_);
            rms = std::sqrt(rr_new / static_cast<double>(rows_));
            const std::size_t completed = iteration + 1U;
            if (rms <= target) {
                queue_.memcpy(x.data(), dx, rows_ * sizeof(double)).wait_and_throw();
                cleanup();
                return {completed, rms, true};
            }
            if (!std::isfinite(rr_new) || rr == 0.0) {
                queue_.memcpy(x.data(), dx, rows_ * sizeof(double)).wait_and_throw();
                cleanup();
                return {completed, rms, false};
            }
            const double beta = rr_new / rr;
            queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
                const std::size_t i = id[0];
                dp[i] = dr[i] + beta * dp[i];
            });
            queue_.wait_and_throw();
            rr = rr_new;
        }
        queue_.memcpy(x.data(), dx, rows_ * sizeof(double)).wait_and_throw();
        cleanup();
        return {max_iterations, rms, false};
    } catch (...) {
        cleanup();
        throw;
    }
}

} // namespace cfd::core

#endif
