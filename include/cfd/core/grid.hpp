#pragma once

#include <cstddef>
#include <stdexcept>

namespace cfd::core {

struct Grid2D {
    std::size_t nx{};
    std::size_t ny{};
    double dx{1.0};
    double dy{1.0};

    Grid2D(std::size_t nx_, std::size_t ny_, double dx_ = 1.0, double dy_ = 1.0)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_) {
        if (nx < 2 || ny < 2 || dx <= 0.0 || dy <= 0.0) {
            throw std::invalid_argument("Grid2D requires nx,ny >= 2 and positive spacing");
        }
    }

    [[nodiscard]] std::size_t cells() const noexcept { return nx * ny; }
    [[nodiscard]] std::size_t index(std::size_t x, std::size_t y) const noexcept { return y * nx + x; }
};

} // namespace cfd::core
