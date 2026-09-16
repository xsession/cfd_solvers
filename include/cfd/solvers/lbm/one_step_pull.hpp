#pragma once

#include "cfd/core/parallel.hpp"
#include "cfd/core/soa_field.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>

namespace cfd::lbm {

// Conventional two-lattice pull-stream/collide implementation. This is kept as a
// correctness oracle and performance baseline for the single-lattice solver.
template<class Descriptor>
class OneStepPullSolver {
public:
    static constexpr int q = Descriptor::q;

    explicit OneStepPullSolver(InPlaceLbmConfig config)
        : config_(normalized_config(config)),
          cells_(config_.nx * config_.ny * config_.nz),
          f_(cells_),
          next_(cells_) {
        if (config_.tau <= 0.5F) {
            throw std::invalid_argument("LBM tau must be > 0.5 for positive viscosity");
        }
        initialize_uniform();
    }

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F, float uz = 0.0F) {
        initialize([=](std::size_t, std::size_t, std::size_t) {
            return std::array<float, 4>{rho, ux, uy, uz};
        });
    }

    void initialize_taylor_green(float amplitude = 0.03F) {
        const float two_pi = 2.0F * std::numbers::pi_v<float>;
        initialize([=, this](std::size_t x, std::size_t y, std::size_t z) {
            const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(config_.nx);
            const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(config_.ny);
            const float sx = std::sin(two_pi * xf);
            const float cx = std::cos(two_pi * xf);
            const float sy = std::sin(two_pi * yf);
            const float cy = std::cos(two_pi * yf);
            float z_factor = 1.0F;
            if constexpr (Descriptor::dimensions == 3) {
                const float zf = (static_cast<float>(z) + 0.5F) / static_cast<float>(config_.nz);
                z_factor = std::cos(two_pi * zf);
            }
            return std::array<float, 4>{1.0F,
                                        amplitude * sx * cy * z_factor,
                                       -amplitude * cx * sy * z_factor,
                                        0.0F};
        });
    }

    void step(std::size_t count = 1) {
        for (std::size_t iteration = 0; iteration < count; ++iteration) step_once();
    }

    [[nodiscard]] MacroscopicFields compute_macroscopic() const {
        MacroscopicFields fields{
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F)
        };
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            for (int d = 0; d < q; ++d) {
                const float value = f_(static_cast<std::size_t>(d), n);
                rho += value;
                ux += value * static_cast<float>(Descriptor::cx(d));
                uy += value * static_cast<float>(Descriptor::cy(d));
                uz += value * static_cast<float>(Descriptor::cz(d));
            }
            if (rho > 0.0F) {
                ux /= rho;
                uy /= rho;
                uz /= rho;
            }
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;
            fields.rho[n] = rho;
            fields.ux[n] = ux;
            fields.uy[n] = uy;
            fields.uz[n] = uz;
        });
        return fields;
    }

    [[nodiscard]] double mass() const {
        return cfd::core::parallel_sum(cells_, [&](std::size_t n) {
            float rho = 0.0F;
            for (int d = 0; d < q; ++d) rho += f_(static_cast<std::size_t>(d), n);
            return rho;
        });
    }

    [[nodiscard]] std::size_t population_bytes() const noexcept {
        return 2U * static_cast<std::size_t>(q) * cells_ * sizeof(float);
    }

    [[nodiscard]] const InPlaceLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }

private:
    InPlaceLbmConfig config_;
    std::size_t cells_{};
    cfd::core::StaticSoA<float, static_cast<std::size_t>(q)> f_;
    cfd::core::StaticSoA<float, static_cast<std::size_t>(q)> next_;

    static InPlaceLbmConfig normalized_config(InPlaceLbmConfig config) {
        if (config.nx < 2 || config.ny < 2) {
            throw std::invalid_argument("LBM grid requires nx, ny >= 2");
        }
        if constexpr (Descriptor::dimensions == 2) {
            config.nz = 1;
        } else if (config.nz < 2) {
            throw std::invalid_argument("3-D LBM grid requires nz >= 2");
        }
        return config;
    }

    template<class Initializer>
    void initialize(Initializer&& initializer) {
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;
            const auto state = initializer(x, y, z);
            for (int d = 0; d < q; ++d) {
                f_(static_cast<std::size_t>(d), n) = detail::equilibrium<Descriptor>(
                    d, state[0], state[1], state[2], state[3]);
            }
        });
    }

    [[nodiscard]] std::size_t source_index(std::size_t x, std::size_t y, std::size_t z,
                                           int direction) const noexcept {
        const std::size_t xs = detail::periodic_shift(x, -Descriptor::cx(direction), config_.nx);
        const std::size_t ys = detail::periodic_shift(y, -Descriptor::cy(direction), config_.ny);
        const std::size_t zs = detail::periodic_shift(z, -Descriptor::cz(direction), config_.nz);
        return (zs * config_.ny + ys) * config_.nx + xs;
    }

    void step_once() {
        const float omega = 1.0F / config_.tau;
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;

            std::array<float, q> fin{};
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            for (int d = 0; d < q; ++d) {
                const std::size_t source = source_index(x, y, z, d);
                const float value = f_(static_cast<std::size_t>(d), source);
                fin[static_cast<std::size_t>(d)] = value;
                rho += value;
                ux += value * static_cast<float>(Descriptor::cx(d));
                uy += value * static_cast<float>(Descriptor::cy(d));
                uz += value * static_cast<float>(Descriptor::cz(d));
            }
            if (rho > 0.0F) {
                ux /= rho;
                uy /= rho;
                uz /= rho;
            }
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;
            for (int d = 0; d < q; ++d) {
                const std::size_t index = static_cast<std::size_t>(d);
                const float feq = detail::equilibrium<Descriptor>(d, rho, ux, uy, uz);
                next_(index, n) = fin[index] - omega * (fin[index] - feq) +
                    detail::guo_force<Descriptor>(d, rho, ux, uy, uz,
                        config_.acceleration_x, config_.acceleration_y, config_.acceleration_z, omega);
            }
        });
        std::swap(f_, next_);
    }
};

using D2Q9PullSolver = OneStepPullSolver<D2Q9InPlaceDescriptor>;
using D3Q19PullSolver = OneStepPullSolver<D3Q19Descriptor>;
using D3Q27PullSolver = OneStepPullSolver<D3Q27Descriptor>;

} // namespace cfd::lbm
