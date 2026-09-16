#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/grid.hpp"
#include "cfd/core/soa_field.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace cfd::lbm {

struct D2Q9Config {
    std::size_t nx{128};
    std::size_t ny{128};
    float tau{0.6F};
    float acceleration_x{0.0F};
    float acceleration_y{0.0F};
};

struct D2Q9Descriptor {
    static constexpr int q = 9;

    static constexpr int cx(int i) noexcept {
        constexpr int values[q]{0, 1, 0, -1, 0, 1, -1, -1, 1};
        return values[i];
    }
    static constexpr int cy(int i) noexcept {
        constexpr int values[q]{0, 0, 1, 0, -1, 1, 1, -1, -1};
        return values[i];
    }
    static constexpr int opposite(int i) noexcept {
        constexpr int values[q]{0, 3, 4, 1, 2, 7, 8, 5, 6};
        return values[i];
    }
    static constexpr float weight(int i) noexcept {
        constexpr float values[q]{4.0F/9.0F, 1.0F/9.0F, 1.0F/9.0F, 1.0F/9.0F, 1.0F/9.0F,
                                  1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F, 1.0F/36.0F};
        return values[i];
    }

    static constexpr float equilibrium(int i, float rho, float ux, float uy) noexcept {
        const float cu = 3.0F * (static_cast<float>(cx(i)) * ux + static_cast<float>(cy(i)) * uy);
        const float uu = 1.5F * (ux * ux + uy * uy);
        return weight(i) * rho * (1.0F + cu + 0.5F * cu * cu - uu);
    }
};

class D2Q9Solver {
public:
    explicit D2Q9Solver(D2Q9Config config);

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F);
    void initialize_taylor_green(float amplitude = 0.03F);

    void set_body_acceleration(float ax, float ay) noexcept;
    void set_solid(std::size_t x, std::size_t y, bool solid = true);
    void set_wall_velocity(std::size_t x, std::size_t y, float ux, float uy);
    void set_velocity_inlet_left(float ux, float uy = 0.0F) noexcept;
    void set_pressure_outlet_right(float rho, float uy = 0.0F) noexcept;
    void clear_x_boundaries() noexcept;

    void step(std::size_t count = 1);
    void compute_macroscopic();

    [[nodiscard]] double mass() const;
    [[nodiscard]] double kinetic_energy() const;
    [[nodiscard]] float max_speed() const;
    [[nodiscard]] double average_ux(std::size_t x_begin, std::size_t x_end,
                                    std::size_t y_begin, std::size_t y_end) const;

    [[nodiscard]] const D2Q9Config& config() const noexcept { return config_; }
    [[nodiscard]] const cfd::core::Grid2D& grid() const noexcept { return grid_; }
    [[nodiscard]] const cfd::core::AlignedVector<float>& density() const noexcept { return rho_; }
    [[nodiscard]] const cfd::core::AlignedVector<float>& velocity_x() const noexcept { return ux_; }
    [[nodiscard]] const cfd::core::AlignedVector<float>& velocity_y() const noexcept { return uy_; }
    [[nodiscard]] const cfd::core::AlignedVector<std::uint8_t>& solid_mask() const noexcept { return solid_; }
    [[nodiscard]] const cfd::core::StaticSoA<float, D2Q9Descriptor::q>& populations() const noexcept { return f_; }

private:
    D2Q9Config config_;
    cfd::core::Grid2D grid_;
    cfd::core::StaticSoA<float, D2Q9Descriptor::q> f_;
    cfd::core::StaticSoA<float, D2Q9Descriptor::q> next_;
    cfd::core::AlignedVector<float> rho_;
    cfd::core::AlignedVector<float> ux_;
    cfd::core::AlignedVector<float> uy_;
    cfd::core::AlignedVector<std::uint8_t> solid_;
    cfd::core::AlignedVector<float> wall_ux_;
    cfd::core::AlignedVector<float> wall_uy_;

    bool left_velocity_inlet_{false};
    float inlet_ux_{0.0F};
    float inlet_uy_{0.0F};
    bool right_pressure_outlet_{false};
    float outlet_rho_{1.0F};
    float outlet_uy_{0.0F};

    void step_once();
    void reconstruct_left_velocity(std::array<float, D2Q9Descriptor::q>& fin) const noexcept;
    void reconstruct_right_pressure(std::array<float, D2Q9Descriptor::q>& fin) const noexcept;
};

} // namespace cfd::lbm
