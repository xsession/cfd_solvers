#include "cfd/multiphysics/partitioned.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {

PartitionedCouplerResult solve_partitioned_fixed_point(std::span<double> state,
                                                        const PartitionedUpdate& update,
                                                        const PartitionedCouplerConfig& config) {
    if (state.empty() || !update) throw std::invalid_argument("partitioned coupling requires state/update");
    if (config.max_iterations == 0U || !(config.absolute_tolerance > 0.0)
        || !(config.minimum_relaxation > 0.0)
        || !(config.maximum_relaxation >= config.minimum_relaxation)
        || !(config.initial_relaxation >= config.minimum_relaxation
             && config.initial_relaxation <= config.maximum_relaxation)) {
        throw std::invalid_argument("invalid partitioned coupling controls");
    }

    const std::size_t n = state.size();
    std::vector<double> candidate(n, 0.0);
    std::vector<double> residual(n, 0.0);
    std::vector<double> previous_residual(n, 0.0);
    double omega = config.initial_relaxation;
    double rms = 0.0;

    for (std::size_t iteration = 0U; ; ++iteration) {
        update(std::span<const double>(state.data(), n), std::span<double>(candidate.data(), n));
        double rr = 0.0;
        for (std::size_t i = 0U; i < n; ++i) {
            residual[i] = candidate[i] - state[i];
            rr += residual[i] * residual[i];
        }
        rms = std::sqrt(rr / static_cast<double>(n));
        if (!std::isfinite(rms)) return {iteration, rms, omega, false};
        if (rms <= config.absolute_tolerance) return {iteration, rms, omega, true};
        if (iteration == config.max_iterations) return {iteration, rms, omega, false};

        if (config.aitken && iteration > 0U) {
            double numerator = 0.0;
            double denominator = 0.0;
            for (std::size_t i = 0U; i < n; ++i) {
                const double delta = residual[i] - previous_residual[i];
                numerator += previous_residual[i] * delta;
                denominator += delta * delta;
            }
            if (denominator > 1.0e-30) {
                omega = std::clamp(-omega * numerator / denominator,
                                   config.minimum_relaxation,
                                   config.maximum_relaxation);
            }
        }
        for (std::size_t i = 0U; i < n; ++i) {
            state[i] += omega * residual[i];
            previous_residual[i] = residual[i];
        }
    }
}

} // namespace cfd::multiphysics
