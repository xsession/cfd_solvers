#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#endif

namespace cfd::distributed {

struct DeviceAssignment {
    std::size_t local_rank{};
    std::size_t local_size{1};
    std::size_t visible_devices{};
    std::size_t device_ordinal{};
    bool oversubscribed{};
};

struct RankDeviceGroup {
    std::size_t local_rank{};
    std::size_t local_size{1};
    std::size_t visible_devices{};
    std::vector<std::size_t> device_ordinals;
    bool oversubscribed{};
};

struct DeviceWorkSlice {
    std::size_t device_ordinal{};
    std::size_t begin{};
    std::size_t end{};
};

inline DeviceAssignment assign_device(std::size_t local_rank,
                                      std::size_t local_size,
                                      std::size_t visible_devices) {
    if (local_size == 0U || local_rank >= local_size) throw std::invalid_argument("invalid local rank topology");
    if (visible_devices == 0U) return {local_rank, local_size, 0U, 0U, false};
    return {local_rank, local_size, visible_devices, local_rank % visible_devices, local_size > visible_devices};
}

// Assign one or more accelerators to each local MPI rank. When accelerators
// outnumber local ranks, devices are partitioned into non-overlapping groups.
// When ranks outnumber devices, the single-device mapping is intentionally
// oversubscribed. max_devices_per_rank == 0 means no explicit cap.
inline RankDeviceGroup assign_device_group(std::size_t local_rank,
                                           std::size_t local_size,
                                           std::size_t visible_devices,
                                           std::size_t max_devices_per_rank = 0U) {
    if (local_size == 0U || local_rank >= local_size) throw std::invalid_argument("invalid local rank topology");
    RankDeviceGroup group{local_rank, local_size, visible_devices, {}, false};
    if (visible_devices == 0U) return group;
    if (visible_devices < local_size) {
        group.device_ordinals.push_back(local_rank % visible_devices);
        group.oversubscribed = true;
        return group;
    }

    const std::size_t base = visible_devices / local_size;
    const std::size_t remainder = visible_devices % local_size;
    const std::size_t natural_count = base + (local_rank < remainder ? 1U : 0U);
    const std::size_t count = max_devices_per_rank == 0U
        ? natural_count
        : std::min(natural_count, max_devices_per_rank);
    const std::size_t first = local_rank * base + std::min(local_rank, remainder);
    for (std::size_t i = 0U; i < count; ++i) group.device_ordinals.push_back(first + i);
    return group;
}

inline std::vector<DeviceWorkSlice> split_rank_work(std::size_t work_items,
                                                     const RankDeviceGroup& group) {
    if (group.device_ordinals.empty()) return {};
    std::vector<DeviceWorkSlice> slices;
    slices.reserve(group.device_ordinals.size());
    const std::size_t base = work_items / group.device_ordinals.size();
    const std::size_t remainder = work_items % group.device_ordinals.size();
    std::size_t begin = 0U;
    for (std::size_t i = 0U; i < group.device_ordinals.size(); ++i) {
        const std::size_t count = base + (i < remainder ? 1U : 0U);
        slices.push_back({group.device_ordinals[i], begin, begin + count});
        begin += count;
    }
    return slices;
}

#if defined(CFD_HAS_SYCL)
inline std::vector<sycl::device> visible_accelerators() {
    std::vector<sycl::device> devices;
    for (const auto& device : sycl::device::get_devices()) {
        if (device.is_gpu() || device.is_accelerator()) devices.push_back(device);
    }
    if (devices.empty()) {
        for (const auto& device : sycl::device::get_devices()) {
            if (!device.is_host()) devices.push_back(device);
        }
    }
    return devices;
}

inline sycl::device device_for_local_rank(std::size_t local_rank, std::size_t local_size) {
    auto devices = visible_accelerators();
    if (devices.empty()) return sycl::device{sycl::default_selector_v};
    const auto assignment = assign_device(local_rank, local_size, devices.size());
    return devices[assignment.device_ordinal];
}

inline std::vector<sycl::device> devices_for_local_rank(std::size_t local_rank,
                                                         std::size_t local_size,
                                                         std::size_t max_devices_per_rank = 0U) {
    auto devices = visible_accelerators();
    if (devices.empty()) return {sycl::device{sycl::default_selector_v}};
    const auto group = assign_device_group(local_rank, local_size, devices.size(), max_devices_per_rank);
    std::vector<sycl::device> selected;
    selected.reserve(group.device_ordinals.size());
    for (const auto ordinal : group.device_ordinals) selected.push_back(devices.at(ordinal));
    return selected;
}

inline std::vector<sycl::queue> queues_for_local_rank(std::size_t local_rank,
                                                       std::size_t local_size,
                                                       std::size_t max_devices_per_rank = 0U) {
    std::vector<sycl::queue> queues;
    for (const auto& device : devices_for_local_rank(local_rank, local_size, max_devices_per_rank)) {
        queues.emplace_back(device, sycl::property::queue::in_order{});
    }
    return queues;
}

struct SyclRankDeviceGroup {
    std::vector<sycl::queue> queues;
    std::vector<DeviceWorkSlice> slices;

    void wait() {
        for (auto& queue : queues) queue.wait_and_throw();
    }
};

inline SyclRankDeviceGroup make_rank_device_group(std::size_t local_rank,
                                                   std::size_t local_size,
                                                   std::size_t work_items,
                                                   std::size_t max_devices_per_rank = 0U) {
    auto devices = visible_accelerators();
    if (devices.empty()) devices.push_back(sycl::device{sycl::default_selector_v});
    const auto assignment = assign_device_group(local_rank, local_size, devices.size(), max_devices_per_rank);
    auto slices = split_rank_work(work_items, assignment);
    SyclRankDeviceGroup group;
    group.slices = std::move(slices);
    group.queues.reserve(assignment.device_ordinals.size());
    for (const auto ordinal : assignment.device_ordinals) {
        group.queues.emplace_back(devices.at(ordinal), sycl::property::queue::in_order{});
    }
    return group;
}
#endif

} // namespace cfd::distributed
