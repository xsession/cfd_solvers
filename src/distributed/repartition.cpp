#include "cfd/distributed/repartition.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cfd::distributed {
namespace {
void validate_boundaries(std::span<const std::size_t> boundaries) {
    if (boundaries.size() < 2U || boundaries.front() != 0U) {
        throw std::invalid_argument("partition boundaries must begin at zero");
    }
    for (std::size_t i = 1U; i < boundaries.size(); ++i) {
        if (boundaries[i] < boundaries[i - 1U]) throw std::invalid_argument("partition boundaries must be monotone");
    }
}
} // namespace

ContiguousRepartitionPlan plan_contiguous_repartition(std::span<const double> weights,
                                                       std::size_t parts,
                                                       std::size_t minimum_cells_per_part) {
    if (weights.empty() || parts == 0U || minimum_cells_per_part == 0U ||
        weights.size() < parts * minimum_cells_per_part) {
        throw std::invalid_argument("invalid contiguous repartition dimensions");
    }
    std::vector<double> prefix(weights.size() + 1U, 0.0);
    for (std::size_t i = 0U; i < weights.size(); ++i) {
        if (!(weights[i] >= 0.0) || !std::isfinite(weights[i])) throw std::invalid_argument("invalid repartition weight");
        prefix[i + 1U] = prefix[i] + weights[i];
    }
    const double total = prefix.back();
    ContiguousRepartitionPlan plan;
    plan.boundaries.assign(parts + 1U, 0U);
    plan.boundaries.back() = weights.size();

    for (std::size_t part = 1U; part < parts; ++part) {
        const std::size_t minimum = plan.boundaries[part - 1U] + minimum_cells_per_part;
        const std::size_t maximum = weights.size() - (parts - part) * minimum_cells_per_part;
        if (!(total > 0.0)) {
            const std::size_t remaining_cells = weights.size() - plan.boundaries[part - 1U];
            const std::size_t remaining_parts = parts - part + 1U;
            plan.boundaries[part] = std::clamp(plan.boundaries[part - 1U] + remaining_cells / remaining_parts,
                                               minimum, maximum);
            continue;
        }
        const double target = total * static_cast<double>(part) / static_cast<double>(parts);
        auto it = std::lower_bound(prefix.begin() + static_cast<std::ptrdiff_t>(minimum),
                                   prefix.begin() + static_cast<std::ptrdiff_t>(maximum + 1U), target);
        std::size_t selected = static_cast<std::size_t>(it - prefix.begin());
        if (selected > minimum) {
            const std::size_t previous = selected - 1U;
            if (std::abs(prefix[previous] - target) <= std::abs(prefix[selected] - target)) selected = previous;
        }
        plan.boundaries[part] = std::clamp(selected, minimum, maximum);
    }

    plan.loads.resize(parts, 0.0);
    double maximum_load = 0.0;
    for (std::size_t part = 0U; part < parts; ++part) {
        plan.loads[part] = prefix[plan.boundaries[part + 1U]] - prefix[plan.boundaries[part]];
        maximum_load = std::max(maximum_load, plan.loads[part]);
    }
    plan.average_load = total / static_cast<double>(parts);
    plan.imbalance = plan.average_load > 0.0 ? maximum_load / plan.average_load : 1.0;
    return plan;
}

std::vector<MigrationSegment> build_migration_segments(std::span<const std::size_t> old_boundaries,
                                                        std::span<const std::size_t> new_boundaries) {
    validate_boundaries(old_boundaries);
    validate_boundaries(new_boundaries);
    if (old_boundaries.back() != new_boundaries.back()) {
        throw std::invalid_argument("old and new partitions must cover the same global range");
    }
    std::vector<MigrationSegment> segments;
    for (std::size_t source = 0U; source + 1U < old_boundaries.size(); ++source) {
        for (std::size_t destination = 0U; destination + 1U < new_boundaries.size(); ++destination) {
            const auto begin = std::max(old_boundaries[source], new_boundaries[destination]);
            const auto end = std::min(old_boundaries[source + 1U], new_boundaries[destination + 1U]);
            if (begin < end && source != destination) segments.push_back({source, destination, begin, end});
        }
    }
    return segments;
}

bool repartition_recommended(double current_imbalance,
                             double proposed_imbalance,
                             double trigger_imbalance,
                             double minimum_relative_improvement) {
    if (!(current_imbalance >= 1.0) || !(proposed_imbalance >= 1.0) || !(trigger_imbalance >= 1.0) ||
        !(minimum_relative_improvement >= 0.0 && minimum_relative_improvement < 1.0)) {
        throw std::invalid_argument("invalid repartition decision controls");
    }
    if (current_imbalance < trigger_imbalance || proposed_imbalance >= current_imbalance) return false;
    return (current_imbalance - proposed_imbalance) / current_imbalance >= minimum_relative_improvement;
}

} // namespace cfd::distributed
