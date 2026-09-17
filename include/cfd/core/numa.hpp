#pragma once

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

#if defined(CFD_HAS_OPENMP)
#include <omp.h>
#endif

namespace cfd::core {

enum class NumaPlacementPolicy {
    inherit,
    compact,
    spread,
    explicit_node,
};

struct NumaNodeInfo {
    int id{-1};
    std::vector<unsigned> cpus;
};

struct NumaPlacementPlan {
    NumaPlacementPolicy policy{NumaPlacementPolicy::inherit};
    std::vector<unsigned> thread_cpus;
    std::vector<int> thread_nodes;

    [[nodiscard]] bool binds_threads() const noexcept { return !thread_cpus.empty(); }
    [[nodiscard]] std::size_t thread_count() const noexcept { return thread_cpus.size(); }
};

// Linux exposes CPU lists using strings such as "0-3,8,10-11". Keeping the
// parser public makes topology handling independently testable on every OS.
[[nodiscard]] std::vector<unsigned> parse_cpu_list(std::string_view text);

// Returns NUMA nodes intersected with the CPUs currently allowed to this
// process. On systems without discoverable NUMA topology, one synthetic node
// containing the available CPUs is returned.
[[nodiscard]] std::vector<NumaNodeInfo> discover_numa_nodes();

// Build a deterministic per-thread placement plan. A thread_count of zero uses
// the number of CPUs visible to the process. explicit_node requires node_id.
[[nodiscard]] NumaPlacementPlan plan_numa_thread_placement(
    std::size_t thread_count,
    NumaPlacementPolicy policy,
    std::optional<int> node_id = std::nullopt);

// Bind only the calling thread. False means the platform/runtime could not
// apply the requested CPU affinity; no exception is thrown for unsupported OSes.
[[nodiscard]] bool bind_current_thread_to_cpu(unsigned cpu) noexcept;

// Applies a plan to an OpenMP team (or to the calling thread without OpenMP).
// This is intentionally explicit and does not silently mutate OMP_* variables.
[[nodiscard]] bool apply_numa_thread_placement(const NumaPlacementPlan& plan) noexcept;

// Portable first-touch baseline used when no explicit placement is requested.
template<class T>
void numa_first_touch(std::span<T> values, const T& value = T{}) {
    parallel_for(values.size(), [&](std::size_t i) { values[i] = value; });
}

// Placement-aware first touch. Each worker binds before touching its contiguous
// range so page ownership follows the requested CPU/NUMA policy on OSes that
// implement thread affinity.
template<class T>
void numa_first_touch(std::span<T> values,
                      const T& value,
                      const NumaPlacementPlan& plan) {
    if (values.empty()) return;
    if (!plan.binds_threads()) {
        numa_first_touch(values, value);
        return;
    }
#if defined(CFD_HAS_OPENMP)
    const int requested = static_cast<int>(std::min<std::size_t>(
        plan.thread_count(), static_cast<std::size_t>(std::numeric_limits<int>::max())));
    #pragma omp parallel num_threads(requested)
    {
        const auto tid = static_cast<std::size_t>(omp_get_thread_num());
        const auto team = static_cast<std::size_t>(omp_get_num_threads());
        if (tid < plan.thread_cpus.size()) {
            (void)bind_current_thread_to_cpu(plan.thread_cpus[tid]);
        }
        const std::size_t begin = values.size() * tid / team;
        const std::size_t end = values.size() * (tid + 1U) / team;
        for (std::size_t i = begin; i < end; ++i) values[i] = value;
    }
#else
    (void)bind_current_thread_to_cpu(plan.thread_cpus.front());
    std::fill(values.begin(), values.end(), value);
#endif
}

} // namespace cfd::core
