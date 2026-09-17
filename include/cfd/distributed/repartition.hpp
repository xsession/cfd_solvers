#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::distributed {

struct ContiguousRepartitionPlan {
    std::vector<std::size_t> boundaries;
    std::vector<double> loads;
    double average_load{};
    double imbalance{};
};

struct MigrationSegment {
    std::size_t source_part{};
    std::size_t destination_part{};
    std::size_t global_begin{};
    std::size_t global_end{};
};

[[nodiscard]] ContiguousRepartitionPlan plan_contiguous_repartition(std::span<const double> weights,
                                                                     std::size_t parts,
                                                                     std::size_t minimum_cells_per_part = 1U);
[[nodiscard]] std::vector<MigrationSegment> build_migration_segments(std::span<const std::size_t> old_boundaries,
                                                                      std::span<const std::size_t> new_boundaries);
[[nodiscard]] bool repartition_recommended(double current_imbalance,
                                           double proposed_imbalance,
                                           double trigger_imbalance = 1.15,
                                           double minimum_relative_improvement = 0.05);

} // namespace cfd::distributed
