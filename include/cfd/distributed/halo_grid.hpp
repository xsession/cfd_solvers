#pragma once

#include "cfd/core/soa_field.hpp"
#include "cfd/distributed/cartesian.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfd::distributed {

template<class T, std::size_t Components>
class HaloSoA3D {
public:
    explicit HaloSoA3D(Extent3 interior, std::size_t halo = 1)
        : interior_(interior), halo_(halo),
          total_{interior.x + 2U * halo, interior.y + 2U * halo, interior.z + 2U * halo},
          data_(total_.cells()) {
        if (halo_ != 1) throw std::invalid_argument("Phase-3 halo grid currently supports halo width 1");
        if (interior_.x == 0 || interior_.y == 0 || interior_.z == 0) {
            throw std::invalid_argument("halo grid interior must be non-empty");
        }
    }

    [[nodiscard]] const Extent3& interior_extent() const noexcept { return interior_; }
    [[nodiscard]] const Extent3& total_extent() const noexcept { return total_; }
    [[nodiscard]] std::size_t halo() const noexcept { return halo_; }
    [[nodiscard]] std::size_t interior_cells() const noexcept { return interior_.cells(); }

    [[nodiscard]] std::size_t index(std::size_t x, std::size_t y, std::size_t z) const noexcept {
        return (z * total_.y + y) * total_.x + x;
    }

    [[nodiscard]] T& operator()(std::size_t component, std::size_t x, std::size_t y, std::size_t z) noexcept {
        return data_(component, index(x, y, z));
    }
    [[nodiscard]] const T& operator()(std::size_t component, std::size_t x, std::size_t y, std::size_t z) const noexcept {
        return data_(component, index(x, y, z));
    }

    [[nodiscard]] cfd::core::StaticSoA<T, Components>& storage() noexcept { return data_; }
    [[nodiscard]] const cfd::core::StaticSoA<T, Components>& storage() const noexcept { return data_; }

    [[nodiscard]] std::size_t halo_value_count(NeighborOffset offset) const noexcept {
        const std::size_t nx = offset.x == 0 ? interior_.x : 1U;
        const std::size_t ny = offset.y == 0 ? interior_.y : 1U;
        const std::size_t nz = offset.z == 0 ? interior_.z : 1U;
        return Components * nx * ny * nz;
    }

    [[nodiscard]] std::vector<T> pack(NeighborOffset offset) const {
        std::vector<T> values;
        values.reserve(halo_value_count(offset));
        const auto [xb, xe] = send_range(interior_.x, offset.x);
        const auto [yb, ye] = send_range(interior_.y, offset.y);
        const auto [zb, ze] = send_range(interior_.z, offset.z);
        for (std::size_t c = 0; c < Components; ++c) {
            for (std::size_t z = zb; z < ze; ++z) {
                for (std::size_t y = yb; y < ye; ++y) {
                    for (std::size_t x = xb; x < xe; ++x) values.push_back((*this)(c, x, y, z));
                }
            }
        }
        return values;
    }

    void unpack(NeighborOffset incoming_from, const std::vector<T>& values) {
        if (values.size() != halo_value_count(incoming_from)) throw std::invalid_argument("halo payload size mismatch");
        const auto [xb, xe] = receive_range(interior_.x, incoming_from.x);
        const auto [yb, ye] = receive_range(interior_.y, incoming_from.y);
        const auto [zb, ze] = receive_range(interior_.z, incoming_from.z);
        std::size_t p = 0;
        for (std::size_t c = 0; c < Components; ++c) {
            for (std::size_t z = zb; z < ze; ++z) {
                for (std::size_t y = yb; y < ye; ++y) {
                    for (std::size_t x = xb; x < xe; ++x) (*this)(c, x, y, z) = values[p++];
                }
            }
        }
    }

private:
    Extent3 interior_{};
    std::size_t halo_{1};
    Extent3 total_{};
    cfd::core::StaticSoA<T, Components> data_;

    [[nodiscard]] static std::pair<std::size_t, std::size_t> send_range(std::size_t n, int offset) noexcept {
        if (offset < 0) return {1U, 2U};
        if (offset > 0) return {n, n + 1U};
        return {1U, n + 1U};
    }

    [[nodiscard]] static std::pair<std::size_t, std::size_t> receive_range(std::size_t n, int offset) noexcept {
        if (offset < 0) return {0U, 1U};
        if (offset > 0) return {n + 1U, n + 2U};
        return {1U, n + 1U};
    }
};

} // namespace cfd::distributed
