#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace cfd::core {

struct SlabPartition {
    std::size_t begin{};
    std::size_t end{};
    std::size_t halo{1};

    [[nodiscard]] std::size_t size() const noexcept { return end - begin; }
};

inline SlabPartition partition_1d(std::size_t global_size,
                                  std::size_t rank,
                                  std::size_t ranks,
                                  std::size_t halo = 1) {
    if (ranks == 0 || rank >= ranks) throw std::invalid_argument("invalid partition rank/count");
    const std::size_t base = global_size / ranks;
    const std::size_t extra = global_size % ranks;
    const std::size_t begin = rank * base + std::min(rank, extra);
    const std::size_t local = base + (rank < extra ? 1U : 0U);
    return {begin, begin + local, halo};
}

} // namespace cfd::core
