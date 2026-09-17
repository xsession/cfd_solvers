#include "cfd/core/autotune.hpp"
#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/repartition.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void test_device_autotune_database() {
    using namespace cfd::core;
    const DeviceTuningKey key{"ExampleVendor", "Example GPU", "1.2.3"};
    const std::vector<KernelTuningCandidate> candidates{
        {64U, 1U, 2U},
        {128U, 2U, 3U},
        {256U, 1U, 3U}
    };
    KernelTuningDatabase database;
    const auto best = database.autotune(key, "lbm.pull-collide-stream", candidates,
        [](const KernelTuningCandidate& candidate) {
            const double group_penalty = std::abs(static_cast<double>(candidate.work_group_size) - 128.0) / 1000.0;
            const double vector_penalty = candidate.vector_width == 2U ? 0.0 : 0.02;
            const double fusion_bonus = candidate.fused_stages == 3U ? 0.01 : 0.03;
            return 0.1 + group_penalty + vector_penalty + fusion_bonus;
        });
    require(best.candidate.work_group_size == 128U, "autotuner should select the best work-group size");
    require(best.candidate.vector_width == 2U, "autotuner should preserve the best vector width");
    require(best.candidate.fused_stages == 3U, "autotuner should preserve the fused stage count");

    const auto path = std::filesystem::temp_directory_path() / "cfd_solvers_v0142_tuning.tsv";
    database.save(path);
    const auto loaded = KernelTuningDatabase::load(path);
    const auto restored = loaded.lookup(key, "lbm.pull-collide-stream");
    require(restored.has_value(), "saved tuning entry should load by device key");
    require(restored->candidate.work_group_size == 128U, "loaded tuning entry should preserve parameters");
    require(device_tuning_key_string(key).find("Example GPU") != std::string::npos,
            "device key string should be diagnostic-friendly");
    std::filesystem::remove(path);
}

void test_multi_device_rank_assignment() {
    using namespace cfd::distributed;
    const auto rank0 = assign_device_group(0U, 2U, 8U);
    const auto rank1 = assign_device_group(1U, 2U, 8U);
    require(rank0.device_ordinals == std::vector<std::size_t>({0U, 1U, 2U, 3U}),
            "rank 0 should receive the first accelerator group");
    require(rank1.device_ordinals == std::vector<std::size_t>({4U, 5U, 6U, 7U}),
            "rank 1 should receive the second accelerator group");
    const auto capped = assign_device_group(0U, 2U, 8U, 2U);
    require(capped.device_ordinals == std::vector<std::size_t>({0U, 1U}),
            "per-rank device cap should be honored");
    const auto oversubscribed = assign_device_group(3U, 4U, 2U);
    require(oversubscribed.oversubscribed && oversubscribed.device_ordinals.size() == 1U &&
            oversubscribed.device_ordinals.front() == 1U,
            "rank-heavy configurations should fall back to deterministic oversubscription");
    const auto slices = split_rank_work(10U, capped);
    require(slices.size() == 2U && slices[0].begin == 0U && slices[0].end == 5U &&
            slices[1].begin == 5U && slices[1].end == 10U,
            "multi-device work should be divided into contiguous balanced slices");
}

void test_adaptive_contiguous_repartition() {
    using namespace cfd::distributed;
    const std::vector<double> weights{8.0, 8.0, 8.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const auto plan = plan_contiguous_repartition(weights, 2U);
    require(plan.boundaries.size() == 3U && plan.boundaries.front() == 0U && plan.boundaries.back() == weights.size(),
            "repartition plan should cover the entire domain");
    require(plan.imbalance < 1.25, "weighted repartition should reduce load imbalance");

    const std::vector<std::size_t> old_boundaries{0U, 4U, 8U};
    const auto segments = build_migration_segments(old_boundaries, plan.boundaries);
    require(!segments.empty(), "changed partition boundaries should generate migration segments");
    require(repartition_recommended(1.6, plan.imbalance, 1.15, 0.05),
            "meaningful imbalance reduction should trigger repartition");
    require(!repartition_recommended(1.05, 1.01, 1.15, 0.05),
            "balanced domains should not churn partitions");

    const std::vector<double> zero_weights(8U, 0.0);
    const auto even = plan_contiguous_repartition(zero_weights, 4U);
    require(even.boundaries == std::vector<std::size_t>({0U, 2U, 4U, 6U, 8U}),
            "zero-cost domains should fall back to deterministic even partitioning");
}
} // namespace

int main() {
    try {
        test_device_autotune_database();
        test_multi_device_rank_assignment();
        test_adaptive_contiguous_repartition();
        std::cout << "v0.14.2 portability/autotuning tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.14.2 portability/autotuning test failed: " << error.what() << '\n';
        return 1;
    }
}
