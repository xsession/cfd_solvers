#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfd::distributed {

struct Extent3 {
    std::size_t x{1};
    std::size_t y{1};
    std::size_t z{1};

    [[nodiscard]] constexpr std::size_t cells() const noexcept { return x * y * z; }
};

struct Coord3 {
    int x{};
    int y{};
    int z{};
};

struct ProcessGrid {
    std::size_t x{1};
    std::size_t y{1};
    std::size_t z{1};

    [[nodiscard]] constexpr std::size_t size() const noexcept { return x * y * z; }
};

struct Brick {
    Extent3 global{};
    Extent3 begin{};
    Extent3 extent{};
    ProcessGrid processes{};
    Coord3 coordinate{};
    std::size_t rank{};
};

struct NeighborOffset {
    int x{};
    int y{};
    int z{};

    [[nodiscard]] constexpr NeighborOffset opposite() const noexcept { return {-x, -y, -z}; }
    [[nodiscard]] constexpr bool is_center() const noexcept { return x == 0 && y == 0 && z == 0; }
};

inline const std::array<NeighborOffset, 26>& neighbor_offsets() {
    static const std::array<NeighborOffset, 26> offsets = [] {
        std::array<NeighborOffset, 26> result{};
        std::size_t n = 0;
        for (int z = -1; z <= 1; ++z) {
            for (int y = -1; y <= 1; ++y) {
                for (int x = -1; x <= 1; ++x) {
                    if (x == 0 && y == 0 && z == 0) continue;
                    result[n++] = {x, y, z};
                }
            }
        }
        return result;
    }();
    return offsets;
}

[[nodiscard]] constexpr int neighbor_tag(NeighborOffset offset) noexcept {
    return (offset.z + 1) * 9 + (offset.y + 1) * 3 + (offset.x + 1);
}

inline ProcessGrid choose_process_grid(Extent3 global, std::size_t ranks, int dimensions = 3) {
    if (ranks == 0) throw std::invalid_argument("process count must be positive");
    if (global.x == 0 || global.y == 0 || global.z == 0) throw std::invalid_argument("global extents must be positive");
    if (dimensions != 2 && dimensions != 3) throw std::invalid_argument("dimensions must be 2 or 3");

    ProcessGrid best{};
    double best_score = std::numeric_limits<double>::infinity();
    bool found = false;

    for (std::size_t px = 1; px <= ranks; ++px) {
        if (ranks % px != 0 || px > global.x) continue;
        const std::size_t rem = ranks / px;
        for (std::size_t py = 1; py <= rem; ++py) {
            if (rem % py != 0 || py > global.y) continue;
            const std::size_t pz = rem / py;
            if (dimensions == 2 && pz != 1) continue;
            if (pz > global.z) continue;

            const double lx = static_cast<double>(global.x) / static_cast<double>(px);
            const double ly = static_cast<double>(global.y) / static_cast<double>(py);
            const double lz = static_cast<double>(global.z) / static_cast<double>(pz);
            const double surface = dimensions == 2 ? 2.0 * (lx + ly)
                                                   : 2.0 * (lx * ly + lx * lz + ly * lz);
            const double max_l = std::max({lx, ly, lz});
            const double min_l = std::min({lx, ly, lz});
            const double aspect_penalty = max_l / std::max(min_l, 1.0e-12);
            const double score = surface * (1.0 + 0.02 * aspect_penalty);
            if (score < best_score) {
                best_score = score;
                best = {px, py, pz};
                found = true;
            }
        }
    }

    if (!found) throw std::invalid_argument("cannot map ranks onto non-empty Cartesian bricks");
    return best;
}

[[nodiscard]] inline Coord3 coords_from_rank(std::size_t rank, ProcessGrid grid) {
    if (rank >= grid.size()) throw std::out_of_range("rank outside process grid");
    const std::size_t xy = grid.x * grid.y;
    return {
        static_cast<int>(rank % grid.x),
        static_cast<int>((rank / grid.x) % grid.y),
        static_cast<int>(rank / xy)};
}

[[nodiscard]] inline std::size_t rank_from_coords(Coord3 coord, ProcessGrid grid) {
    if (coord.x < 0 || coord.y < 0 || coord.z < 0 ||
        static_cast<std::size_t>(coord.x) >= grid.x ||
        static_cast<std::size_t>(coord.y) >= grid.y ||
        static_cast<std::size_t>(coord.z) >= grid.z) {
        throw std::out_of_range("coordinates outside process grid");
    }
    return (static_cast<std::size_t>(coord.z) * grid.y + static_cast<std::size_t>(coord.y)) * grid.x +
           static_cast<std::size_t>(coord.x);
}

[[nodiscard]] inline std::pair<std::size_t, std::size_t>
partition_axis(std::size_t global, std::size_t coordinate, std::size_t partitions) {
    if (partitions == 0 || coordinate >= partitions) throw std::invalid_argument("invalid axis partition");
    const std::size_t base = global / partitions;
    const std::size_t extra = global % partitions;
    const std::size_t begin = coordinate * base + std::min(coordinate, extra);
    const std::size_t extent = base + (coordinate < extra ? 1U : 0U);
    return {begin, extent};
}

[[nodiscard]] inline Brick make_brick(Extent3 global, ProcessGrid grid, Coord3 coord, std::size_t rank) {
    if (grid.size() == 0 || rank >= grid.size()) throw std::invalid_argument("invalid process grid/rank");
    if (coord.x < 0 || coord.y < 0 || coord.z < 0 ||
        static_cast<std::size_t>(coord.x) >= grid.x ||
        static_cast<std::size_t>(coord.y) >= grid.y ||
        static_cast<std::size_t>(coord.z) >= grid.z) {
        throw std::invalid_argument("invalid process-grid coordinates");
    }
    const auto [bx, nx] = partition_axis(global.x, static_cast<std::size_t>(coord.x), grid.x);
    const auto [by, ny] = partition_axis(global.y, static_cast<std::size_t>(coord.y), grid.y);
    const auto [bz, nz] = partition_axis(global.z, static_cast<std::size_t>(coord.z), grid.z);
    return {global, {bx, by, bz}, {nx, ny, nz}, grid, coord, rank};
}

[[nodiscard]] inline Brick make_brick(Extent3 global, ProcessGrid grid, std::size_t rank) {
    return make_brick(global, grid, coords_from_rank(rank, grid), rank);
}

[[nodiscard]] inline int neighbor_rank(const Brick& brick, NeighborOffset offset,
                                       std::array<bool, 3> periodic = {true, true, true}) {
    if (offset.is_center()) return static_cast<int>(brick.rank);
    Coord3 c{brick.coordinate.x + offset.x, brick.coordinate.y + offset.y, brick.coordinate.z + offset.z};
    const std::array<int, 3> dims{
        static_cast<int>(brick.processes.x), static_cast<int>(brick.processes.y), static_cast<int>(brick.processes.z)};
    int* values[3]{&c.x, &c.y, &c.z};
    for (int axis = 0; axis < 3; ++axis) {
        if (*values[axis] < 0 || *values[axis] >= dims[axis]) {
            if (!periodic[static_cast<std::size_t>(axis)]) return -1;
            int wrapped = *values[axis] % dims[axis];
            if (wrapped < 0) wrapped += dims[axis];
            *values[axis] = wrapped;
        }
    }
    return static_cast<int>(rank_from_coords(c, brick.processes));
}

} // namespace cfd::distributed
