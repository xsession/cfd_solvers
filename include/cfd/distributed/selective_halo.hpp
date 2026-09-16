#pragma once

#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/halo_grid.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfd::distributed {

template<class Descriptor>
[[nodiscard]] constexpr bool population_crosses(int direction, NeighborOffset offset) noexcept {
    if (offset.is_center()) return false;
    if (offset.x != 0 && Descriptor::cx(direction) != offset.x) return false;
    if (offset.y != 0 && Descriptor::cy(direction) != offset.y) return false;
    if (offset.z != 0 && Descriptor::cz(direction) != offset.z) return false;
    return true;
}

template<class Descriptor>
[[nodiscard]] inline std::vector<int> crossing_populations(NeighborOffset offset) {
    std::vector<int> result;
    result.reserve(static_cast<std::size_t>(Descriptor::q));
    for (int d = 0; d < Descriptor::q; ++d) {
        if (population_crosses<Descriptor>(d, offset)) result.push_back(d);
    }
    return result;
}

template<class Descriptor>
[[nodiscard]] constexpr std::array<int, static_cast<std::size_t>(Descriptor::q)>
crossing_population_array(NeighborOffset offset, std::size_t& count) noexcept {
    std::array<int, static_cast<std::size_t>(Descriptor::q)> result{};
    count = 0;
    for (int d = 0; d < Descriptor::q; ++d) {
        if (population_crosses<Descriptor>(d, offset)) result[count++] = d;
    }
    return result;
}

[[nodiscard]] constexpr std::size_t boundary_cells(Extent3 interior, NeighborOffset offset) noexcept {
    const std::size_t nx = offset.x == 0 ? interior.x : 1U;
    const std::size_t ny = offset.y == 0 ? interior.y : 1U;
    const std::size_t nz = offset.z == 0 ? interior.z : 1U;
    return nx * ny * nz;
}

template<class Descriptor>
[[nodiscard]] inline std::size_t selective_halo_value_count(Extent3 interior, NeighborOffset offset) {
    return crossing_populations<Descriptor>(offset).size() * boundary_cells(interior, offset);
}

template<class Descriptor, class T, std::size_t Components>
[[nodiscard]] inline std::size_t selective_halo_value_count(const HaloSoA3D<T, Components>& grid,
                                                             NeighborOffset offset) {
    static_assert(Components == static_cast<std::size_t>(Descriptor::q),
                  "descriptor population count must match halo-grid components");
    return selective_halo_value_count<Descriptor>(grid.interior_extent(), offset);
}

struct SelectiveHaloEntry {
    int direction{};
    std::size_t x{};
    std::size_t y{};
    std::size_t z{};
};

template<class Descriptor>
[[nodiscard]] constexpr SelectiveHaloEntry
selective_send_entry(Extent3 interior, NeighborOffset offset, std::size_t linear) {
    std::size_t direction_count = 0;
    const auto directions = crossing_population_array<Descriptor>(offset, direction_count);
    const std::size_t cells = boundary_cells(interior, offset);
    if (cells == 0 || linear >= direction_count * cells) return {};
    const std::size_t direction_index = linear / cells;
    const std::size_t local = linear - direction_index * cells;
    const std::size_t rx = offset.x == 0 ? interior.x : 1U;
    const std::size_t ry = offset.y == 0 ? interior.y : 1U;
    const std::size_t lx = local % rx;
    const std::size_t yz = local / rx;
    const std::size_t ly = yz % ry;
    const std::size_t lz = yz / ry;
    return {
        directions[direction_index],
        offset.x < 0 ? 1U : (offset.x > 0 ? interior.x : lx + 1U),
        offset.y < 0 ? 1U : (offset.y > 0 ? interior.y : ly + 1U),
        offset.z < 0 ? 1U : (offset.z > 0 ? interior.z : lz + 1U)};
}

template<class Descriptor>
[[nodiscard]] constexpr SelectiveHaloEntry
selective_receive_entry(Extent3 interior, NeighborOffset incoming_from, std::size_t linear) {
    std::size_t direction_count = 0;
    const auto directions = crossing_population_array<Descriptor>(incoming_from.opposite(), direction_count);
    const std::size_t cells = boundary_cells(interior, incoming_from);
    if (cells == 0 || linear >= direction_count * cells) return {};
    const std::size_t direction_index = linear / cells;
    const std::size_t local = linear - direction_index * cells;
    const std::size_t rx = incoming_from.x == 0 ? interior.x : 1U;
    const std::size_t ry = incoming_from.y == 0 ? interior.y : 1U;
    const std::size_t lx = local % rx;
    const std::size_t yz = local / rx;
    const std::size_t ly = yz % ry;
    const std::size_t lz = yz / ry;
    return {
        directions[direction_index],
        incoming_from.x < 0 ? 0U : (incoming_from.x > 0 ? interior.x + 1U : lx + 1U),
        incoming_from.y < 0 ? 0U : (incoming_from.y > 0 ? interior.y + 1U : ly + 1U),
        incoming_from.z < 0 ? 0U : (incoming_from.z > 0 ? interior.z + 1U : lz + 1U)};
}


template<class Descriptor, class T, std::size_t Components>
void pack_selective_into(const HaloSoA3D<T, Components>& grid,
                         NeighborOffset offset,
                         T* output,
                         std::size_t output_count) {
    static_assert(Components == static_cast<std::size_t>(Descriptor::q),
                  "descriptor population count must match halo-grid components");
    const std::size_t expected = selective_halo_value_count<Descriptor>(grid, offset);
    if (output_count != expected) throw std::invalid_argument("selective halo output size mismatch");
    if (expected != 0 && output == nullptr) throw std::invalid_argument("selective halo output is null");
    const auto e = grid.interior_extent();
    for (std::size_t p = 0; p < expected; ++p) {
        const auto entry = selective_send_entry<Descriptor>(e, offset, p);
        output[p] = grid(static_cast<std::size_t>(entry.direction), entry.x, entry.y, entry.z);
    }
}

template<class Descriptor, class T, std::size_t Components>
[[nodiscard]] std::vector<T> pack_selective(const HaloSoA3D<T, Components>& grid,
                                             NeighborOffset offset) {
    std::vector<T> result(selective_halo_value_count<Descriptor>(grid, offset));
    pack_selective_into<Descriptor>(grid, offset, result.data(), result.size());
    return result;
}

template<class Descriptor, class T, std::size_t Components>
void unpack_selective(HaloSoA3D<T, Components>& grid,
                      NeighborOffset incoming_from,
                      const T* values,
                      std::size_t value_count) {
    static_assert(Components == static_cast<std::size_t>(Descriptor::q),
                  "descriptor population count must match halo-grid components");
    const std::size_t expected = selective_halo_value_count<Descriptor>(grid, incoming_from);
    if (value_count != expected) throw std::invalid_argument("selective halo payload size mismatch");
    if (expected != 0 && values == nullptr) throw std::invalid_argument("selective halo payload is null");
    const auto e = grid.interior_extent();
    for (std::size_t p = 0; p < expected; ++p) {
        const auto entry = selective_receive_entry<Descriptor>(e, incoming_from, p);
        grid(static_cast<std::size_t>(entry.direction), entry.x, entry.y, entry.z) = values[p];
    }
}

template<class Descriptor, class T, std::size_t Components>
void unpack_selective(HaloSoA3D<T, Components>& grid,
                      NeighborOffset incoming_from,
                      const std::vector<T>& values) {
    unpack_selective<Descriptor>(grid, incoming_from, values.data(), values.size());
}

struct HaloTrafficStats {
    std::size_t messages{};
    std::size_t values{};
    std::size_t bytes{};
};

template<class Descriptor, class T = float>
[[nodiscard]] inline HaloTrafficStats selective_halo_traffic(Extent3 interior) {
    HaloTrafficStats stats{};
    for (const auto offset : neighbor_offsets()) {
        const auto values = selective_halo_value_count<Descriptor>(interior, offset);
        if (values == 0) continue;
        ++stats.messages;
        stats.values += values;
    }
    stats.bytes = stats.values * sizeof(T);
    return stats;
}

template<class Descriptor, class T = float>
[[nodiscard]] inline HaloTrafficStats full_halo_traffic(Extent3 interior) {
    HaloTrafficStats stats{};
    for (const auto offset : neighbor_offsets()) {
        const std::size_t values = static_cast<std::size_t>(Descriptor::q) * boundary_cells(interior, offset);
        ++stats.messages;
        stats.values += values;
    }
    stats.bytes = stats.values * sizeof(T);
    return stats;
}

} // namespace cfd::distributed
