#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::multiphysics {

struct PartitionedCouplerConfig {
    std::size_t max_iterations{100U};
    double absolute_tolerance{1.0e-10};
    double initial_relaxation{1.0};
    bool aitken{true};
    double minimum_relaxation{0.05};
    double maximum_relaxation{1.5};
};

struct PartitionedCouplerResult {
    std::size_t iterations{};
    double residual_rms{};
    double relaxation{};
    bool converged{};
};

using PartitionedUpdate = std::function<void(std::span<const double>, std::span<double>)>;

[[nodiscard]] PartitionedCouplerResult solve_partitioned_fixed_point(
    std::span<double> state,
    const PartitionedUpdate& update,
    const PartitionedCouplerConfig& config = {});

} // namespace cfd::multiphysics
