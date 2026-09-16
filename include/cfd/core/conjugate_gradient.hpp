#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace cfd::core {

struct ConjugateGradientResult {
    std::size_t iterations{};
    double residual_rms{};
    bool converged{};
};

struct ConjugateGradientWorkspace {
    AlignedVector<double> residual;
    AlignedVector<double> direction;
    AlignedVector<double> applied;

    void resize(std::size_t n) {
        residual.resize(n);
        direction.resize(n);
        applied.resize(n);
    }
};

// Matrix-free CG for symmetric positive (semi-)definite operators. For a
// semidefinite operator the caller is responsible for supplying an RHS and
// initial guess in the compatible subspace (e.g. mean-zero periodic Poisson).
template<class ApplyOperator>
ConjugateGradientResult conjugate_gradient(std::span<const double> rhs,
                                           std::span<double> x,
                                           ApplyOperator&& apply,
                                           ConjugateGradientWorkspace& workspace,
                                           std::size_t max_iterations,
                                           double relative_tolerance) {
    if (rhs.empty() || x.size() != rhs.size()) throw std::invalid_argument("CG vector size mismatch");
    if (max_iterations == 0 || !(relative_tolerance > 0.0)) throw std::invalid_argument("invalid CG controls");
    const std::size_t n = rhs.size();
    workspace.resize(n);
    auto& r = workspace.residual;
    auto& d = workspace.direction;
    auto& q = workspace.applied;

    apply(std::span<const double>(x.data(), n), std::span<double>(q.data(), n));
    parallel_for(n, [&](std::size_t i) {
        r[i] = rhs[i] - q[i];
        d[i] = r[i];
    });
    const double rhs2 = parallel_sum(n, [&](std::size_t i) { return rhs[i] * rhs[i]; });
    double rr = parallel_sum(n, [&](std::size_t i) { return r[i] * r[i]; });
    const double scale = std::max(std::sqrt(rhs2 / static_cast<double>(n)), 1.0);
    const double target = relative_tolerance * scale;
    double rms = std::sqrt(rr / static_cast<double>(n));
    if (rms <= target) return {0, rms, true};

    std::size_t iteration = 0;
    for (; iteration < max_iterations; ++iteration) {
        apply(std::span<const double>(d.data(), n), std::span<double>(q.data(), n));
        const double dq = parallel_sum(n, [&](std::size_t i) { return d[i] * q[i]; });
        if (!(dq > 0.0) || !std::isfinite(dq)) break;
        const double alpha = rr / dq;
        parallel_for(n, [&](std::size_t i) {
            x[i] += alpha * d[i];
            r[i] -= alpha * q[i];
        });
        const double rr_new = parallel_sum(n, [&](std::size_t i) { return r[i] * r[i]; });
        rms = std::sqrt(rr_new / static_cast<double>(n));
        ++iteration;
        if (rms <= target) return {iteration, rms, true};
        const double beta = rr_new / rr;
        parallel_for(n, [&](std::size_t i) { d[i] = r[i] + beta * d[i]; });
        rr = rr_new;
    }
    return {iteration, rms, false};
}

} // namespace cfd::core
