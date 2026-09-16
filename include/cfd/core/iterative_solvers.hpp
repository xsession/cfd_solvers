#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cfd::core {

struct IterativeSolverResult {
    std::size_t iterations{};
    double residual_rms{};
    bool converged{};
};

class JacobiPreconditioner {
public:
    explicit JacobiPreconditioner(std::span<const double> diagonal, double minimum_magnitude = 1.0e-30)
        : inverse_(diagonal.size()) {
        if (!(minimum_magnitude > 0.0)) throw std::invalid_argument("invalid Jacobi diagonal floor");
        parallel_for(diagonal.size(), [&](std::size_t i) {
            const double d = diagonal[i];
            inverse_[i] = std::abs(d) >= minimum_magnitude ? 1.0 / d : 1.0;
        });
    }

    void operator()(std::span<const double> rhs, std::span<double> result) const {
        if (rhs.size() != inverse_.size() || result.size() != inverse_.size()) {
            throw std::invalid_argument("Jacobi preconditioner size mismatch");
        }
        parallel_for(rhs.size(), [&](std::size_t i) { result[i] = inverse_[i] * rhs[i]; });
    }

private:
    std::vector<double> inverse_;
};

struct KrylovWorkspace {
    AlignedVector<double> r;
    AlignedVector<double> p;
    AlignedVector<double> q;
    AlignedVector<double> z;
    AlignedVector<double> r0;
    AlignedVector<double> v;
    AlignedVector<double> s;
    AlignedVector<double> t;

    void resize(std::size_t n) {
        r.resize(n); p.resize(n); q.resize(n); z.resize(n);
        r0.resize(n); v.resize(n); s.resize(n); t.resize(n);
    }
};

inline double vector_rms(std::span<const double> values) {
    if (values.empty()) return 0.0;
    const double sum2 = parallel_sum(values.size(), [&](std::size_t i) { return values[i] * values[i]; });
    return std::sqrt(sum2 / static_cast<double>(values.size()));
}

template<class ApplyOperator, class ApplyPreconditioner>
IterativeSolverResult preconditioned_conjugate_gradient(std::span<const double> rhs,
                                                         std::span<double> x,
                                                         ApplyOperator&& apply,
                                                         ApplyPreconditioner&& precondition,
                                                         KrylovWorkspace& workspace,
                                                         std::size_t max_iterations,
                                                         double relative_tolerance) {
    if (rhs.empty() || x.size() != rhs.size()) throw std::invalid_argument("PCG vector size mismatch");
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) throw std::invalid_argument("invalid PCG controls");
    const std::size_t n = rhs.size();
    workspace.resize(n);
    auto& r = workspace.r;
    auto& p = workspace.p;
    auto& q = workspace.q;
    auto& z = workspace.z;

    apply(std::span<const double>(x.data(), n), std::span<double>(q.data(), n));
    parallel_for(n, [&](std::size_t i) { r[i] = rhs[i] - q[i]; });
    const double rhs_rms = std::max(vector_rms(rhs), 1.0);
    const double target = relative_tolerance * rhs_rms;
    double rms = vector_rms(std::span<const double>(r.data(), n));
    if (rms <= target) return {0U, rms, true};

    precondition(std::span<const double>(r.data(), n), std::span<double>(z.data(), n));
    parallel_for(n, [&](std::size_t i) { p[i] = z[i]; });
    double rz = parallel_sum(n, [&](std::size_t i) { return r[i] * z[i]; });
    if (!std::isfinite(rz)) return {0U, rms, false};

    for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
        apply(std::span<const double>(p.data(), n), std::span<double>(q.data(), n));
        const double pq = parallel_sum(n, [&](std::size_t i) { return p[i] * q[i]; });
        if (!(pq > 0.0) || !std::isfinite(pq)) return {iteration, rms, false};
        const double alpha = rz / pq;
        parallel_for(n, [&](std::size_t i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
        });
        rms = vector_rms(std::span<const double>(r.data(), n));
        const std::size_t completed = iteration + 1U;
        if (rms <= target) return {completed, rms, true};
        precondition(std::span<const double>(r.data(), n), std::span<double>(z.data(), n));
        const double rz_new = parallel_sum(n, [&](std::size_t i) { return r[i] * z[i]; });
        if (!std::isfinite(rz_new) || rz == 0.0) return {completed, rms, false};
        const double beta = rz_new / rz;
        parallel_for(n, [&](std::size_t i) { p[i] = z[i] + beta * p[i]; });
        rz = rz_new;
    }
    return {max_iterations, rms, false};
}

template<class ApplyOperator, class ApplyPreconditioner>
IterativeSolverResult bicgstab(std::span<const double> rhs,
                               std::span<double> x,
                               ApplyOperator&& apply,
                               ApplyPreconditioner&& precondition,
                               KrylovWorkspace& workspace,
                               std::size_t max_iterations,
                               double relative_tolerance) {
    if (rhs.empty() || x.size() != rhs.size()) throw std::invalid_argument("BiCGStab vector size mismatch");
    if (max_iterations == 0U || !(relative_tolerance > 0.0)) throw std::invalid_argument("invalid BiCGStab controls");
    const std::size_t n = rhs.size();
    workspace.resize(n);
    auto& r = workspace.r;
    auto& r0 = workspace.r0;
    auto& p = workspace.p;
    auto& v = workspace.v;
    auto& z = workspace.z;
    auto& s = workspace.s;
    auto& t = workspace.t;
    auto& q = workspace.q;

    apply(std::span<const double>(x.data(), n), std::span<double>(q.data(), n));
    parallel_for(n, [&](std::size_t i) {
        r[i] = rhs[i] - q[i];
        r0[i] = r[i];
        p[i] = 0.0;
        v[i] = 0.0;
    });
    const double rhs_rms = std::max(vector_rms(rhs), 1.0);
    const double target = relative_tolerance * rhs_rms;
    double rms = vector_rms(std::span<const double>(r.data(), n));
    if (rms <= target) return {0U, rms, true};

    double rho_old = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    const double tiny = 64.0 * std::numeric_limits<double>::epsilon();

    for (std::size_t iteration = 0U; iteration < max_iterations; ++iteration) {
        const double rho = parallel_sum(n, [&](std::size_t i) { return r0[i] * r[i]; });
        if (std::abs(rho) <= tiny || !std::isfinite(rho)) return {iteration, rms, false};
        const double beta = (rho / rho_old) * (alpha / omega);
        parallel_for(n, [&](std::size_t i) { p[i] = r[i] + beta * (p[i] - omega * v[i]); });

        precondition(std::span<const double>(p.data(), n), std::span<double>(z.data(), n));
        apply(std::span<const double>(z.data(), n), std::span<double>(v.data(), n));
        const double denom = parallel_sum(n, [&](std::size_t i) { return r0[i] * v[i]; });
        if (std::abs(denom) <= tiny || !std::isfinite(denom)) return {iteration, rms, false};
        alpha = rho / denom;
        parallel_for(n, [&](std::size_t i) { s[i] = r[i] - alpha * v[i]; });
        const double s_rms = vector_rms(std::span<const double>(s.data(), n));
        if (s_rms <= target) {
            parallel_for(n, [&](std::size_t i) { x[i] += alpha * z[i]; });
            return {iteration + 1U, s_rms, true};
        }

        // Reuse q as preconditioned s, leaving z available for the p correction.
        precondition(std::span<const double>(s.data(), n), std::span<double>(q.data(), n));
        apply(std::span<const double>(q.data(), n), std::span<double>(t.data(), n));
        const double tt = parallel_sum(n, [&](std::size_t i) { return t[i] * t[i]; });
        if (tt <= tiny || !std::isfinite(tt)) return {iteration, rms, false};
        omega = parallel_sum(n, [&](std::size_t i) { return t[i] * s[i]; }) / tt;
        if (std::abs(omega) <= tiny || !std::isfinite(omega)) return {iteration, rms, false};
        parallel_for(n, [&](std::size_t i) {
            x[i] += alpha * z[i] + omega * q[i];
            r[i] = s[i] - omega * t[i];
        });
        rms = vector_rms(std::span<const double>(r.data(), n));
        const std::size_t completed = iteration + 1U;
        if (rms <= target) return {completed, rms, true};
        rho_old = rho;
    }
    return {max_iterations, rms, false};
}


class Ilu0Preconditioner {
public:
    explicit Ilu0Preconditioner(const CsrMatrix& matrix, double diagonal_floor = 1.0e-30)
        : n_(matrix.rows()),
          row_offsets_(matrix.row_offsets()),
          columns_(matrix.column_indices()),
          lu_(matrix.values()),
          diagonal_position_(n_, static_cast<std::size_t>(-1)),
          diagonal_floor_(diagonal_floor) {
        if (matrix.rows() != matrix.cols()) throw std::invalid_argument("ILU(0) requires a square matrix");
        if (!(diagonal_floor_ > 0.0)) throw std::invalid_argument("invalid ILU(0) diagonal floor");
        factorize();
    }

    void operator()(std::span<const double> rhs, std::span<double> result) const {
        if (rhs.size() != n_ || result.size() != n_) throw std::invalid_argument("ILU(0) size mismatch");
        // In-place forward/back substitution in result. This intentionally stays
        // serial: triangular dependencies dominate and the first implementation
        // favors correctness/readability over level-scheduling complexity.
        for (std::size_t i = 0; i < n_; ++i) {
            double sum = rhs[i];
            for (std::size_t k = row_offsets_[i]; k < diagonal_position_[i]; ++k) {
                sum -= lu_[k] * result[columns_[k]];
            }
            result[i] = sum;
        }
        for (std::size_t ii = n_; ii-- > 0U;) {
            double sum = result[ii];
            for (std::size_t k = diagonal_position_[ii] + 1U; k < row_offsets_[ii + 1U]; ++k) {
                sum -= lu_[k] * result[columns_[k]];
            }
            result[ii] = sum / lu_[diagonal_position_[ii]];
        }
    }

private:
    std::size_t n_{};
    std::vector<std::size_t> row_offsets_;
    std::vector<std::size_t> columns_;
    std::vector<double> lu_;
    std::vector<std::size_t> diagonal_position_;
    double diagonal_floor_{};

    [[nodiscard]] std::size_t find_in_row(std::size_t row, std::size_t column) const noexcept {
        const auto first = columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
        const auto last = columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1U]);
        const auto it = std::lower_bound(first, last, column);
        if (it == last || *it != column) return static_cast<std::size_t>(-1);
        return static_cast<std::size_t>(it - columns_.begin());
    }

    void factorize() {
        for (std::size_t i = 0; i < n_; ++i) {
            const std::size_t diag = find_in_row(i, i);
            if (diag == static_cast<std::size_t>(-1)) throw std::invalid_argument("ILU(0) matrix is missing a diagonal entry");
            diagonal_position_[i] = diag;

            for (std::size_t k = row_offsets_[i]; k < diag; ++k) {
                const std::size_t j = columns_[k];
                const std::size_t jdiag = diagonal_position_[j];
                if (jdiag == static_cast<std::size_t>(-1)) throw std::runtime_error("ILU(0) invalid row ordering");
                double pivot = lu_[jdiag];
                if (std::abs(pivot) < diagonal_floor_) {
                    pivot = std::copysign(diagonal_floor_, pivot == 0.0 ? 1.0 : pivot);
                    lu_[jdiag] = pivot;
                }
                const double lij = lu_[k] / pivot;
                lu_[k] = lij;
                for (std::size_t m = jdiag + 1U; m < row_offsets_[j + 1U]; ++m) {
                    const std::size_t target = find_in_row(i, columns_[m]);
                    if (target != static_cast<std::size_t>(-1)) lu_[target] -= lij * lu_[m];
                }
            }
            if (std::abs(lu_[diag]) < diagonal_floor_) {
                lu_[diag] = std::copysign(diagonal_floor_, lu_[diag] == 0.0 ? 1.0 : lu_[diag]);
            }
        }
    }
};

template<class ApplyOperator, class ApplyPreconditioner>
IterativeSolverResult restarted_gmres(std::span<const double> rhs,
                                      std::span<double> x,
                                      ApplyOperator&& apply,
                                      ApplyPreconditioner&& precondition,
                                      std::size_t max_iterations,
                                      std::size_t restart,
                                      double relative_tolerance) {
    if (rhs.empty() || x.size() != rhs.size()) throw std::invalid_argument("GMRES vector size mismatch");
    if (max_iterations == 0U || restart == 0U || !(relative_tolerance > 0.0)) {
        throw std::invalid_argument("invalid GMRES controls");
    }
    const std::size_t n = rhs.size();
    restart = std::min(restart, max_iterations);
    AlignedVector<double> ax(n), r(n), z(n), w(n);
    std::vector<AlignedVector<double>> v(restart + 1U, AlignedVector<double>(n));
    std::vector<double> h((restart + 1U) * restart, 0.0);
    std::vector<double> cs(restart, 0.0), sn(restart, 0.0), g(restart + 1U, 0.0), y(restart, 0.0);
    const double rhs_rms = std::max(vector_rms(rhs), 1.0);
    const double target = relative_tolerance * rhs_rms;
    std::size_t completed = 0U;
    double actual_rms = std::numeric_limits<double>::infinity();

    auto recompute_residual = [&]() {
        apply(std::span<const double>(x.data(), n), std::span<double>(ax.data(), n));
        parallel_for(n, [&](std::size_t i) { r[i] = rhs[i] - ax[i]; });
        actual_rms = vector_rms(std::span<const double>(r.data(), n));
        return actual_rms;
    };

    if (recompute_residual() <= target) return {0U, actual_rms, true};

    while (completed < max_iterations) {
        precondition(std::span<const double>(r.data(), n), std::span<double>(z.data(), n));
        const double beta2 = parallel_sum(n, [&](std::size_t i) { return z[i] * z[i]; });
        const double beta = std::sqrt(beta2);
        if (!(beta > 0.0) || !std::isfinite(beta)) return {completed, actual_rms, false};
        parallel_for(n, [&](std::size_t i) { v[0][i] = z[i] / beta; });
        std::fill(h.begin(), h.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        const std::size_t inner_limit = std::min(restart, max_iterations - completed);
        std::size_t used = 0U;
        for (std::size_t j = 0; j < inner_limit; ++j) {
            apply(std::span<const double>(v[j].data(), n), std::span<double>(ax.data(), n));
            precondition(std::span<const double>(ax.data(), n), std::span<double>(w.data(), n));

            for (std::size_t i = 0; i <= j; ++i) {
                const double hij = parallel_sum(n, [&](std::size_t k) { return w[k] * v[i][k]; });
                h[i * restart + j] = hij;
                parallel_for(n, [&](std::size_t k) { w[k] -= hij * v[i][k]; });
            }
            const double hn = std::sqrt(parallel_sum(n, [&](std::size_t k) { return w[k] * w[k]; }));
            h[(j + 1U) * restart + j] = hn;
            if (hn > std::numeric_limits<double>::epsilon()) {
                parallel_for(n, [&](std::size_t k) { v[j + 1U][k] = w[k] / hn; });
            }

            for (std::size_t i = 0; i < j; ++i) {
                const double a = h[i * restart + j];
                const double b = h[(i + 1U) * restart + j];
                h[i * restart + j] = cs[i] * a + sn[i] * b;
                h[(i + 1U) * restart + j] = -sn[i] * a + cs[i] * b;
            }
            const double a = h[j * restart + j];
            const double b = h[(j + 1U) * restart + j];
            const double denom = std::hypot(a, b);
            cs[j] = denom > 0.0 ? a / denom : 1.0;
            sn[j] = denom > 0.0 ? b / denom : 0.0;
            h[j * restart + j] = cs[j] * a + sn[j] * b;
            h[(j + 1U) * restart + j] = 0.0;
            const double gj = g[j];
            g[j] = cs[j] * gj;
            g[j + 1U] = -sn[j] * gj;
            used = j + 1U;
            ++completed;

            const double estimated_rms = std::abs(g[j + 1U]) / std::sqrt(static_cast<double>(n));
            if (estimated_rms <= target || hn <= std::numeric_limits<double>::epsilon() || completed == max_iterations) break;
        }

        for (std::size_t ii = used; ii-- > 0U;) {
            double sum = g[ii];
            for (std::size_t j = ii + 1U; j < used; ++j) sum -= h[ii * restart + j] * y[j];
            const double diag = h[ii * restart + ii];
            if (std::abs(diag) <= std::numeric_limits<double>::epsilon()) return {completed, actual_rms, false};
            y[ii] = sum / diag;
        }
        for (std::size_t j = 0; j < used; ++j) {
            const double alpha = y[j];
            parallel_for(n, [&](std::size_t i) { x[i] += alpha * v[j][i]; });
        }
        if (recompute_residual() <= target) return {completed, actual_rms, true};
    }
    return {completed, actual_rms, false};
}

} // namespace cfd::core
