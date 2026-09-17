#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#endif

namespace cfd::core {

struct DeviceTuningKey {
    std::string vendor;
    std::string name;
    std::string driver;

    [[nodiscard]] bool operator==(const DeviceTuningKey&) const = default;
};

struct KernelTuningCandidate {
    std::size_t work_group_size{1};
    std::size_t vector_width{1};
    std::size_t fused_stages{1};
};

struct KernelTuningRecord {
    DeviceTuningKey device;
    std::string kernel_group;
    KernelTuningCandidate candidate;
    double score_seconds{};
};

class KernelTuningDatabase {
public:
    [[nodiscard]] static KernelTuningDatabase load(const std::filesystem::path& path);
    void save(const std::filesystem::path& path) const;

    void upsert(KernelTuningRecord record);
    [[nodiscard]] std::optional<KernelTuningRecord> lookup(const DeviceTuningKey& device,
                                                            std::string_view kernel_group) const;
    [[nodiscard]] std::span<const KernelTuningRecord> records() const noexcept { return records_; }

    template<class Benchmark>
    [[nodiscard]] KernelTuningRecord autotune(const DeviceTuningKey& device,
                                              std::string kernel_group,
                                              std::span<const KernelTuningCandidate> candidates,
                                              Benchmark&& benchmark) {
        if (kernel_group.empty() || candidates.empty()) {
            throw std::invalid_argument("autotune requires a kernel group and at least one candidate");
        }
        std::optional<KernelTuningRecord> best;
        for (const auto& candidate : candidates) {
            if (candidate.work_group_size == 0U || candidate.vector_width == 0U || candidate.fused_stages == 0U) {
                throw std::invalid_argument("autotune candidate values must be positive");
            }
            const double score = std::invoke(benchmark, candidate);
            if (!(score > 0.0) || !std::isfinite(score)) {
                throw std::runtime_error("autotune benchmark returned an invalid score");
            }
            KernelTuningRecord current{device, kernel_group, candidate, score};
            if (!best || current.score_seconds < best->score_seconds) best = std::move(current);
        }
        upsert(*best);
        return *best;
    }

private:
    std::vector<KernelTuningRecord> records_;
};

[[nodiscard]] std::string device_tuning_key_string(const DeviceTuningKey& key);

#if defined(CFD_HAS_SYCL)
[[nodiscard]] DeviceTuningKey sycl_device_tuning_key(const sycl::device& device);
#endif

} // namespace cfd::core
