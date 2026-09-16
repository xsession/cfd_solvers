#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cfd::core {

struct NewtonConfig {
    std::size_t max_iterations{30U};
    double absolute_tolerance{1.0e-10};
    double relative_tolerance{1.0e-9};
    std::size_t linear_max_iterations{1000U};
    std::size_t gmres_restart{40U};
    double linear_relative_tolerance{1.0e-11};
    std::size_t line_search_steps{12U};
    double minimum_step{1.0e-6};
};

struct NewtonResult {
    std::size_t iterations{};
    double residual_rms{};
    double accepted_step{1.0};
    bool converged{};
};

using NonlinearResidual = std::function<void(std::span<const double>, std::span<double>)>;
using NonlinearJacobian = std::function<CsrMatrix(std::span<const double>)>;

inline NewtonResult newton_solve(std::span<double> x,
                                 const NonlinearResidual& residual,
                                 const NonlinearJacobian& jacobian,
                                 const NewtonConfig& config = {}) {
    if (x.empty() || !residual || !jacobian) throw std::invalid_argument("Newton system is incomplete");
    if (config.max_iterations == 0U || !(config.absolute_tolerance > 0.0)
        || !(config.relative_tolerance > 0.0) || config.linear_max_iterations == 0U
        || config.gmres_restart == 0U || !(config.linear_relative_tolerance > 0.0)
        || config.line_search_steps == 0U || !(config.minimum_step > 0.0)
        || config.minimum_step > 1.0) {
        throw std::invalid_argument("invalid Newton controls");
    }

    const std::size_t n = x.size();
    std::vector<double> r(n, 0.0), trial_r(n, 0.0), rhs(n, 0.0), delta(n, 0.0), trial_x(n, 0.0);
    residual(std::span<const double>(x.data(), n), std::span<double>(r.data(), n));
    double rms = vector_rms(r);
    if (!std::isfinite(rms)) return {0U, rms, 0.0, false};
    const double initial_rms = std::max(rms, std::numeric_limits<double>::min());
    const double target = std::max(config.absolute_tolerance, config.relative_tolerance * initial_rms);
    if (rms <= target) return {0U, rms, 1.0, true};

    double accepted_step = 1.0;
    for (std::size_t iteration = 0U; iteration < config.max_iterations; ++iteration) {
        CsrMatrix j = jacobian(std::span<const double>(x.data(), n));
        if (j.rows() != n || j.cols() != n) throw std::runtime_error("Newton Jacobian size mismatch");
        parallel_for(n, [&](std::size_t i) {
            rhs[i] = -r[i];
            delta[i] = 0.0;
        });
        Ilu0Preconditioner ilu(j);
        const auto linear = restarted_gmres(
            std::span<const double>(rhs.data(), n), std::span<double>(delta.data(), n),
            [&](std::span<const double> v, std::span<double> out) { j.multiply(v, out); },
            [&](std::span<const double> v, std::span<double> out) { ilu(v, out); },
            config.linear_max_iterations, std::min(config.gmres_restart, n), config.linear_relative_tolerance);
        if (!linear.converged) return {iteration, rms, 0.0, false};

        double step = 1.0;
        bool accepted = false;
        double trial_rms = rms;
        for (std::size_t ls = 0U; ls < config.line_search_steps; ++ls) {
            parallel_for(n, [&](std::size_t i) { trial_x[i] = x[i] + step * delta[i]; });
            residual(std::span<const double>(trial_x.data(), n), std::span<double>(trial_r.data(), n));
            trial_rms = vector_rms(trial_r);
            if (std::isfinite(trial_rms) && trial_rms < rms) {
                accepted = true;
                break;
            }
            step *= 0.5;
            if (step < config.minimum_step) break;
        }
        if (!accepted) return {iteration + 1U, rms, step, false};
        std::copy(trial_x.begin(), trial_x.end(), x.begin());
        r.swap(trial_r);
        rms = trial_rms;
        accepted_step = step;
        if (rms <= target) return {iteration + 1U, rms, accepted_step, true};
    }
    return {config.max_iterations, rms, accepted_step, false};
}

} // namespace cfd::core
