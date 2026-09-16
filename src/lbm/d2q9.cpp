#include "cfd/solvers/lbm/d2q9.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::lbm {

namespace {

inline float guo_source(int q, float rho, float ux, float uy,
                        float ax, float ay, float omega) noexcept {
    if (ax == 0.0F && ay == 0.0F) return 0.0F;
    const float cx = static_cast<float>(D2Q9Descriptor::cx(q));
    const float cy = static_cast<float>(D2Q9Descriptor::cy(q));
    const float fx = rho * ax;
    const float fy = rho * ay;
    const float cu = cx * ux + cy * uy;
    const float cf = cx * fx + cy * fy;
    const float uf = ux * fx + uy * fy;
    return D2Q9Descriptor::weight(q) * (1.0F - 0.5F * omega) *
           (3.0F * (cf - uf) + 9.0F * cu * cf);
}

} // namespace

D2Q9Solver::D2Q9Solver(D2Q9Config config)
    : config_(config),
      grid_(config.nx, config.ny),
      f_(grid_.cells()),
      next_(grid_.cells()),
      rho_(grid_.cells(), 1.0F),
      ux_(grid_.cells(), 0.0F),
      uy_(grid_.cells(), 0.0F),
      solid_(grid_.cells(), 0U),
      wall_ux_(grid_.cells(), 0.0F),
      wall_uy_(grid_.cells(), 0.0F) {
    if (config_.tau <= 0.5F) {
        throw std::invalid_argument("D2Q9 tau must be > 0.5 for positive viscosity");
    }
    initialize_uniform();
}

void D2Q9Solver::initialize_uniform(float rho, float ux, float uy) {
    cfd::core::parallel_for(grid_.cells(), [&](std::size_t i) {
        rho_[i] = rho;
        ux_[i] = solid_[i] ? wall_ux_[i] : ux;
        uy_[i] = solid_[i] ? wall_uy_[i] : uy;
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            f_(static_cast<std::size_t>(q), i) = D2Q9Descriptor::equilibrium(q, rho, ux_[i], uy_[i]);
            next_(static_cast<std::size_t>(q), i) = f_(static_cast<std::size_t>(q), i);
        }
    });
}

void D2Q9Solver::initialize_taylor_green(float amplitude) {
    const float two_pi = 2.0F * std::numbers::pi_v<float>;
    cfd::core::parallel_for(grid_.cells(), [&](std::size_t i) {
        const std::size_t x = i % grid_.nx;
        const std::size_t y = i / grid_.nx;
        const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(grid_.nx);
        const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(grid_.ny);
        const float ux = amplitude * std::sin(two_pi * xf) * std::cos(two_pi * yf);
        const float uy = -amplitude * std::cos(two_pi * xf) * std::sin(two_pi * yf);
        rho_[i] = 1.0F;
        ux_[i] = ux;
        uy_[i] = uy;
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            f_(static_cast<std::size_t>(q), i) = D2Q9Descriptor::equilibrium(q, 1.0F, ux, uy);
            next_(static_cast<std::size_t>(q), i) = f_(static_cast<std::size_t>(q), i);
        }
    });
}

void D2Q9Solver::set_body_acceleration(float ax, float ay) noexcept {
    config_.acceleration_x = ax;
    config_.acceleration_y = ay;
}

void D2Q9Solver::set_solid(std::size_t x, std::size_t y, bool solid) {
    if (x >= grid_.nx || y >= grid_.ny) throw std::out_of_range("solid cell outside grid");
    const std::size_t i = grid_.index(x, y);
    solid_[i] = solid ? 1U : 0U;
    // set_solid() always denotes a stationary wall. Use set_wall_velocity()
    // explicitly for moving walls so a reused cell cannot retain stale velocity.
    wall_ux_[i] = 0.0F;
    wall_uy_[i] = 0.0F;
}

void D2Q9Solver::set_wall_velocity(std::size_t x, std::size_t y, float ux, float uy) {
    if (x >= grid_.nx || y >= grid_.ny) throw std::out_of_range("wall cell outside grid");
    const std::size_t i = grid_.index(x, y);
    solid_[i] = 1U;
    wall_ux_[i] = ux;
    wall_uy_[i] = uy;
}

void D2Q9Solver::set_velocity_inlet_left(float ux, float uy) noexcept {
    left_velocity_inlet_ = true;
    inlet_ux_ = ux;
    inlet_uy_ = uy;
}

void D2Q9Solver::set_pressure_outlet_right(float rho, float uy) noexcept {
    right_pressure_outlet_ = true;
    outlet_rho_ = rho;
    outlet_uy_ = uy;
}

void D2Q9Solver::clear_x_boundaries() noexcept {
    left_velocity_inlet_ = false;
    right_pressure_outlet_ = false;
}

void D2Q9Solver::step(std::size_t count) {
    for (std::size_t n = 0; n < count; ++n) step_once();
    compute_macroscopic();
}

void D2Q9Solver::reconstruct_left_velocity(std::array<float, D2Q9Descriptor::q>& f) const noexcept {
    const float ux = inlet_ux_;
    const float uy = inlet_uy_;
    const float rho = (f[0] + f[2] + f[4] + 2.0F * (f[3] + f[6] + f[7])) / (1.0F - ux);
    f[1] = f[3] + (2.0F / 3.0F) * rho * ux;
    f[5] = f[7] + 0.5F * (f[4] - f[2]) + (1.0F / 6.0F) * rho * ux + 0.5F * rho * uy;
    f[8] = f[6] + 0.5F * (f[2] - f[4]) + (1.0F / 6.0F) * rho * ux - 0.5F * rho * uy;
}

void D2Q9Solver::reconstruct_right_pressure(std::array<float, D2Q9Descriptor::q>& f) const noexcept {
    const float rho = outlet_rho_;
    const float uy = outlet_uy_;
    const float ux = -1.0F + (f[0] + f[2] + f[4] + 2.0F * (f[1] + f[5] + f[8])) / rho;
    f[3] = f[1] - (2.0F / 3.0F) * rho * ux;
    f[6] = f[8] + 0.5F * (f[4] - f[2]) - (1.0F / 6.0F) * rho * ux + 0.5F * rho * uy;
    f[7] = f[5] + 0.5F * (f[2] - f[4]) - (1.0F / 6.0F) * rho * ux - 0.5F * rho * uy;
}

void D2Q9Solver::step_once() {
    const float omega = 1.0F / config_.tau;
    const std::size_t nx = grid_.nx;
    const std::size_t ny = grid_.ny;
    const std::size_t cells = grid_.cells();

    cfd::core::parallel_for(cells, [&](std::size_t i) {
        const std::size_t x = i % nx;
        const std::size_t y = i / nx;

        if (solid_[i] != 0U) {
            for (int q = 0; q < D2Q9Descriptor::q; ++q) {
                next_(static_cast<std::size_t>(q), i) = f_(static_cast<std::size_t>(q), i);
            }
            return;
        }

        std::array<float, D2Q9Descriptor::q> fin{};
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            const int cx = D2Q9Descriptor::cx(q);
            const int cy = D2Q9Descriptor::cy(q);
            const int sx_raw = static_cast<int>(x) - cx;
            const int sy_raw = static_cast<int>(y) - cy;

            const bool missing_left = left_velocity_inlet_ && x == 0U && sx_raw < 0;
            const bool missing_right = right_pressure_outlet_ && x + 1U == nx && sx_raw >= static_cast<int>(nx);
            if (missing_left || missing_right) {
                fin[static_cast<std::size_t>(q)] = f_(static_cast<std::size_t>(q), i);
                continue;
            }

            const std::size_t sx = static_cast<std::size_t>((sx_raw + static_cast<int>(nx)) % static_cast<int>(nx));
            const std::size_t sy = static_cast<std::size_t>((sy_raw + static_cast<int>(ny)) % static_cast<int>(ny));
            const std::size_t source = sy * nx + sx;

            float value;
            if (solid_[source] != 0U) {
                const int qo = D2Q9Descriptor::opposite(q);
                const float wall_dot_c = static_cast<float>(cx) * wall_ux_[source] +
                                         static_cast<float>(cy) * wall_uy_[source];
                const float rho_wall = std::max(rho_[i], 1.0e-8F);
                value = f_(static_cast<std::size_t>(qo), i) +
                        6.0F * D2Q9Descriptor::weight(q) * rho_wall * wall_dot_c;
            } else {
                value = f_(static_cast<std::size_t>(q), source);
            }
            fin[static_cast<std::size_t>(q)] = value;
        }

        if (left_velocity_inlet_ && x == 0U) reconstruct_left_velocity(fin);
        if (right_pressure_outlet_ && x + 1U == nx) reconstruct_right_pressure(fin);

        float rho = 0.0F;
        float mx = 0.0F;
        float my = 0.0F;
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            const float value = fin[static_cast<std::size_t>(q)];
            rho += value;
            mx += value * static_cast<float>(D2Q9Descriptor::cx(q));
            my += value * static_cast<float>(D2Q9Descriptor::cy(q));
        }

        float ux = 0.0F;
        float uy = 0.0F;
        if (rho > 0.0F) {
            ux = mx / rho + 0.5F * config_.acceleration_x;
            uy = my / rho + 0.5F * config_.acceleration_y;
        }
        if (left_velocity_inlet_ && x == 0U) {
            ux = inlet_ux_;
            uy = inlet_uy_;
        }
        if (right_pressure_outlet_ && x + 1U == nx) {
            uy = outlet_uy_;
        }

        rho_[i] = rho;
        ux_[i] = ux;
        uy_[i] = uy;

        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            const float feq = D2Q9Descriptor::equilibrium(q, rho, ux, uy);
            next_(static_cast<std::size_t>(q), i) = fin[static_cast<std::size_t>(q)] -
                omega * (fin[static_cast<std::size_t>(q)] - feq) +
                guo_source(q, rho, ux, uy, config_.acceleration_x, config_.acceleration_y, omega);
        }
    });

    f_.raw().swap(next_.raw());
}

void D2Q9Solver::compute_macroscopic() {
    cfd::core::parallel_for(grid_.cells(), [&](std::size_t i) {
        if (solid_[i] != 0U) {
            rho_[i] = 1.0F;
            ux_[i] = wall_ux_[i];
            uy_[i] = wall_uy_[i];
            return;
        }
        float rho = 0.0F;
        float mx = 0.0F;
        float my = 0.0F;
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            const float value = f_(static_cast<std::size_t>(q), i);
            rho += value;
            mx += value * static_cast<float>(D2Q9Descriptor::cx(q));
            my += value * static_cast<float>(D2Q9Descriptor::cy(q));
        }
        rho_[i] = rho;
        if (rho > 0.0F) {
            ux_[i] = mx / rho + 0.5F * config_.acceleration_x;
            uy_[i] = my / rho + 0.5F * config_.acceleration_y;
        } else {
            ux_[i] = 0.0F;
            uy_[i] = 0.0F;
        }
        const std::size_t x = i % grid_.nx;
        if (left_velocity_inlet_ && x == 0U) {
            ux_[i] = inlet_ux_;
            uy_[i] = inlet_uy_;
        }
        if (right_pressure_outlet_ && x + 1U == grid_.nx) {
            uy_[i] = outlet_uy_;
        }
    });
}

double D2Q9Solver::mass() const {
    return cfd::core::parallel_sum(grid_.cells(), [&](std::size_t i) {
        return solid_[i] == 0U ? rho_[i] : 0.0F;
    });
}

double D2Q9Solver::kinetic_energy() const {
    return 0.5 * cfd::core::parallel_sum(grid_.cells(), [&](std::size_t i) {
        if (solid_[i] != 0U) return 0.0;
        return static_cast<double>(rho_[i]) *
               (static_cast<double>(ux_[i]) * static_cast<double>(ux_[i]) +
                static_cast<double>(uy_[i]) * static_cast<double>(uy_[i]));
    });
}

float D2Q9Solver::max_speed() const {
    return static_cast<float>(cfd::core::parallel_max(grid_.cells(), [&](std::size_t i) {
        if (solid_[i] != 0U) return 0.0;
        return std::sqrt(static_cast<double>(ux_[i]) * ux_[i] +
                         static_cast<double>(uy_[i]) * uy_[i]);
    }));
}

double D2Q9Solver::average_ux(std::size_t x_begin, std::size_t x_end,
                              std::size_t y_begin, std::size_t y_end) const {
    if (x_begin >= x_end || y_begin >= y_end || x_end > grid_.nx || y_end > grid_.ny) {
        throw std::out_of_range("average_ux region outside grid");
    }
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t y = y_begin; y < y_end; ++y) {
        for (std::size_t x = x_begin; x < x_end; ++x) {
            const std::size_t i = grid_.index(x, y);
            if (solid_[i] == 0U) {
                sum += ux_[i];
                ++count;
            }
        }
    }
    return count == 0U ? 0.0 : sum / static_cast<double>(count);
}

} // namespace cfd::lbm
