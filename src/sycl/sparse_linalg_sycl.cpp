#include "cfd/core/sycl_sparse.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

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
    : SyclCsrLinearAlgebra(matrix, sycl::queue(device, sycl::property::queue::in_order{})) {}

SyclCsrLinearAlgebra::SyclCsrLinearAlgebra(const CsrMatrix& matrix, sycl::queue queue)
    : queue_(std::move(queue)), rows_(matrix.rows()), cols_(matrix.cols()), nonzeros_(matrix.nonzeros()) {
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
    work_rhs_ = sycl::malloc_device<double>(rows_, queue_);
    work_x_ = sycl::malloc_device<double>(rows_, queue_);
    work_r_ = sycl::malloc_device<double>(rows_, queue_);
    work_p_ = sycl::malloc_device<double>(rows_, queue_);
    work_q_ = sycl::malloc_device<double>(rows_, queue_);
    work_rhat_ = sycl::malloc_device<double>(rows_, queue_);
    work_s_ = sycl::malloc_device<double>(rows_, queue_);
    work_t_ = sycl::malloc_device<double>(rows_, queue_);
    if (!row_offsets_ || !column_indices_ || !values_ || !reduction_value_ ||
        !work_rhs_ || !work_x_ || !work_r_ || !work_p_ || !work_q_ ||
        !work_rhat_ || !work_s_ || !work_t_) {
        release();
        throw std::bad_alloc{};
    }

    try {
        const std::size_t offset_bytes = (rows_ + 1U) * sizeof(std::size_t);
        const std::size_t column_bytes = nonzeros_ * sizeof(std::size_t);
        const std::size_t value_bytes = nonzeros_ * sizeof(double);
        queue_.memcpy(row_offsets_, matrix.row_offsets().data(), offset_bytes);
        queue_.memcpy(column_indices_, matrix.column_indices().data(), column_bytes);
        queue_.memcpy(values_, matrix.values().data(), value_bytes);
        queue_.wait_and_throw();
        transfer_stats_.record_host_to_device(offset_bytes + column_bytes + value_bytes);
        transfer_stats_.record_synchronization();
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
    free_if(work_rhs_, queue_);
    free_if(work_x_, queue_);
    free_if(work_r_, queue_);
    free_if(work_p_, queue_);
    free_if(work_q_, queue_);
    free_if(work_rhat_, queue_);
    free_if(work_s_, queue_);
    free_if(work_t_, queue_);
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
    transfer_stats_.record_synchronization();
    return *reduction_value_;
}

void SyclCsrLinearAlgebra::multiply_device(const double* x_device, double* y_device) {
    if (!x_device || !y_device) throw std::invalid_argument("SYCL CSR device multiply requires non-null pointers");
    spmv_device(x_device, y_device);
    queue_.wait_and_throw();
    transfer_stats_.record_synchronization();
}

void SyclCsrLinearAlgebra::update_values(std::span<const double> values) {
    if (values.size() != nonzeros_) throw std::invalid_argument("SYCL CSR value count mismatch");
    const std::size_t bytes = nonzeros_ * sizeof(double);
    queue_.memcpy(values_, values.data(), bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(bytes);
    transfer_stats_.record_synchronization();
}

void SyclCsrLinearAlgebra::update_values_device(const double* values_device) {
    if (!values_device) throw std::invalid_argument("SYCL CSR device values must be non-null");
    const std::size_t bytes = nonzeros_ * sizeof(double);
    queue_.memcpy(values_, values_device, bytes);
    transfer_stats_.record_device_to_device(bytes);
}

void SyclCsrLinearAlgebra::multiply(std::span<const double> x, std::span<double> y) {
    if (x.size() != cols_ || y.size() != rows_) throw std::invalid_argument("SYCL CSR multiply size mismatch");
    const std::size_t input_bytes = cols_ * sizeof(double);
    const std::size_t output_bytes = rows_ * sizeof(double);
    queue_.memcpy(work_x_, x.data(), input_bytes);
    transfer_stats_.record_host_to_device(input_bytes);
    spmv_device(work_x_, work_q_);
    queue_.memcpy(y.data(), work_q_, output_bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(output_bytes);
    transfer_stats_.record_synchronization();
}

IterativeSolverResult SyclCsrLinearAlgebra::conjugate_gradient_device(
    const double* rhs_device,
    double* x_device,
    std::size_t max_iterations,
    double relative_tolerance) {
    if (!rhs_device || !x_device) {
        throw std::invalid_argument("SYCL CG device vectors must be non-null");
    }
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid SYCL CG controls");
    }

    double* dr = work_r_;
    double* dp = work_p_;
    double* dq = work_q_;

    spmv_device(x_device, dq);
    queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
        const std::size_t i = id[0];
        dr[i] = rhs_device[i] - dq[i];
        dp[i] = dr[i];
    });
    queue_.wait_and_throw();
    transfer_stats_.record_synchronization();

    const double rhs2 = dot_device(rhs_device, rhs_device, rows_);
    double rr = dot_device(dr, dr, rows_);
    const double scale = std::max(std::sqrt(rhs2 / static_cast<double>(rows_)), 1.0);
    const double target = relative_tolerance * scale;
    double rms = std::sqrt(rr / static_cast<double>(rows_));
    if (rms <= target) return {0U, rms, true};

    for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
        spmv_device(dp, dq);
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
        const double dq_dot = dot_device(dp, dq, rows_);
        if (!(dq_dot > 0.0) || !std::isfinite(dq_dot)) {
            return {iteration, rms, false};
        }
        const double alpha = rr / dq_dot;
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            x_device[i] += alpha * dp[i];
            dr[i] -= alpha * dq[i];
        });
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
        const double rr_new = dot_device(dr, dr, rows_);
        rms = std::sqrt(rr_new / static_cast<double>(rows_));
        const std::size_t completed = iteration + 1U;
        if (rms <= target) return {completed, rms, true};
        if (!std::isfinite(rr_new) || rr == 0.0) return {completed, rms, false};

        const double beta = rr_new / rr;
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            dp[i] = dr[i] + beta * dp[i];
        });
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
        rr = rr_new;
    }
    return {max_iterations, rms, false};
}


IterativeSolverResult SyclCsrLinearAlgebra::bicgstab_device(
    const double* rhs_device,
    double* x_device,
    std::size_t max_iterations,
    double relative_tolerance) {
    if (!rhs_device || !x_device) {
        throw std::invalid_argument("SYCL BiCGStab device vectors must be non-null");
    }
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid SYCL BiCGStab controls");
    }

    double* r = work_r_;
    double* p = work_p_;
    double* v = work_q_;
    double* rhat = work_rhat_;
    double* s = work_s_;
    double* t = work_t_;

    spmv_device(x_device, v);
    queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
        const std::size_t i = id[0];
        r[i] = rhs_device[i] - v[i];
        rhat[i] = r[i];
        p[i] = 0.0;
        v[i] = 0.0;
    });
    queue_.wait_and_throw();
    transfer_stats_.record_synchronization();

    const double rhs2 = dot_device(rhs_device, rhs_device, rows_);
    double rr = dot_device(r, r, rows_);
    const double scale = std::max(std::sqrt(rhs2 / static_cast<double>(rows_)), 1.0);
    const double target = relative_tolerance * scale;
    double rms = std::sqrt(rr / static_cast<double>(rows_));
    if (rms <= target) return {0U, rms, true};

    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    constexpr double breakdown = 1.0e-300;

    for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
        const double rho = dot_device(rhat, r, rows_);
        if (!std::isfinite(rho) || std::abs(rho) <= breakdown ||
            !std::isfinite(omega) || std::abs(omega) <= breakdown) {
            return {iteration, rms, false};
        }
        const double beta = (rho / rho_old) * (alpha / omega);
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        });
        spmv_device(p, v);
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();

        const double rhat_v = dot_device(rhat, v, rows_);
        if (!std::isfinite(rhat_v) || std::abs(rhat_v) <= breakdown) {
            return {iteration, rms, false};
        }
        alpha = rho / rhat_v;
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            s[i] = r[i] - alpha * v[i];
        });
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();

        const double ss = dot_device(s, s, rows_);
        double s_rms = std::sqrt(ss / static_cast<double>(rows_));
        if (s_rms <= target) {
            queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
                const std::size_t i = id[0];
                x_device[i] += alpha * p[i];
            }).wait_and_throw();
            transfer_stats_.record_synchronization();
            return {iteration + 1U, s_rms, true};
        }

        spmv_device(s, t);
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
        const double tt = dot_device(t, t, rows_);
        const double ts = dot_device(t, s, rows_);
        if (!std::isfinite(tt) || tt <= breakdown || !std::isfinite(ts)) {
            return {iteration, rms, false};
        }
        omega = ts / tt;
        if (!std::isfinite(omega) || std::abs(omega) <= breakdown) {
            return {iteration, rms, false};
        }
        queue_.parallel_for(sycl::range<1>(rows_), [=](sycl::id<1> id) {
            const std::size_t i = id[0];
            x_device[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        });
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
        rr = dot_device(r, r, rows_);
        rms = std::sqrt(rr / static_cast<double>(rows_));
        if (rms <= target) return {iteration + 1U, rms, true};
        if (!std::isfinite(rr)) return {iteration + 1U, rms, false};
        rho_old = rho;
    }
    return {max_iterations, rms, false};
}

IterativeSolverResult SyclCsrLinearAlgebra::bicgstab(
    std::span<const double> rhs,
    std::span<double> x,
    std::size_t max_iterations,
    double relative_tolerance) {
    if (rhs.empty() || rhs.size() != rows_ || x.size() != rows_) {
        throw std::invalid_argument("SYCL BiCGStab vector size mismatch");
    }
    const std::size_t bytes = rows_ * sizeof(double);
    queue_.memcpy(work_rhs_, rhs.data(), bytes);
    queue_.memcpy(work_x_, x.data(), bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(2U * bytes);
    transfer_stats_.record_synchronization();
    const auto result = bicgstab_device(work_rhs_, work_x_, max_iterations, relative_tolerance);
    queue_.memcpy(x.data(), work_x_, bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(bytes);
    transfer_stats_.record_synchronization();
    return result;
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

    const std::size_t vector_bytes = rows_ * sizeof(double);
    queue_.memcpy(work_rhs_, rhs.data(), vector_bytes);
    queue_.memcpy(work_x_, x.data(), vector_bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(2U * vector_bytes);
    transfer_stats_.record_synchronization();

    const auto result = conjugate_gradient_device(
        work_rhs_, work_x_, max_iterations, relative_tolerance);

    queue_.memcpy(x.data(), work_x_, vector_bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(vector_bytes);
    transfer_stats_.record_synchronization();
    return result;
}

} // namespace cfd::core

#endif
