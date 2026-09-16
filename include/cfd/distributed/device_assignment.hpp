#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#include <vector>
#endif

namespace cfd::distributed {

struct DeviceAssignment {
    std::size_t local_rank{};
    std::size_t local_size{1};
    std::size_t visible_devices{};
    std::size_t device_ordinal{};
    bool oversubscribed{};
};

inline DeviceAssignment assign_device(std::size_t local_rank,
                                      std::size_t local_size,
                                      std::size_t visible_devices) {
    if (local_size == 0 || local_rank >= local_size) throw std::invalid_argument("invalid local rank topology");
    if (visible_devices == 0) return {local_rank, local_size, 0, 0, false};
    return {local_rank, local_size, visible_devices, local_rank % visible_devices, local_size > visible_devices};
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
#endif

} // namespace cfd::distributed
