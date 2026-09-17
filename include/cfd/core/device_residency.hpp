#pragma once

#include <cstddef>
#include <cstdint>

namespace cfd::core {

// Lightweight accounting for explicit host/device traffic. Hot accelerator
// loops should leave these counters unchanged; transfers belong at setup,
// checkpoint, visualization, and API boundaries.
struct DeviceTransferStats {
    std::uint64_t host_to_device_bytes{};
    std::uint64_t device_to_host_bytes{};
    std::uint64_t device_to_device_bytes{};
    std::uint64_t synchronization_points{};

    [[nodiscard]] constexpr std::uint64_t explicit_transfer_bytes() const noexcept {
        return host_to_device_bytes + device_to_host_bytes + device_to_device_bytes;
    }

    [[nodiscard]] constexpr std::uint64_t host_transfer_bytes() const noexcept {
        return host_to_device_bytes + device_to_host_bytes;
    }

    constexpr void record_host_to_device(std::size_t bytes) noexcept {
        host_to_device_bytes += static_cast<std::uint64_t>(bytes);
    }

    constexpr void record_device_to_host(std::size_t bytes) noexcept {
        device_to_host_bytes += static_cast<std::uint64_t>(bytes);
    }

    constexpr void record_device_to_device(std::size_t bytes) noexcept {
        device_to_device_bytes += static_cast<std::uint64_t>(bytes);
    }

    constexpr void record_synchronization(std::uint64_t count = 1U) noexcept {
        synchronization_points += count;
    }

    constexpr void reset() noexcept { *this = {}; }
};

} // namespace cfd::core
