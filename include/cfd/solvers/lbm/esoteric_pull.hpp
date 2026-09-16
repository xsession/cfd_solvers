#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/core/soa_field.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace cfd::lbm {

struct InPlaceLbmConfig {
    std::size_t nx{128};
    std::size_t ny{128};
    std::size_t nz{128};
    float tau{0.6F};
    float acceleration_x{0.0F};
    float acceleration_y{0.0F};
    float acceleration_z{0.0F};
};

struct MacroscopicFields {
    cfd::core::AlignedVector<float> rho;
    cfd::core::AlignedVector<float> ux;
    cfd::core::AlignedVector<float> uy;
    cfd::core::AlignedVector<float> uz;
};

namespace detail {

inline std::size_t periodic_shift(std::size_t coordinate, int delta, std::size_t extent) noexcept {
    if (delta > 0) return coordinate + 1U == extent ? 0U : coordinate + 1U;
    if (delta < 0) return coordinate == 0U ? extent - 1U : coordinate - 1U;
    return coordinate;
}

template<class Descriptor>
inline float equilibrium(int q, float rho, float ux, float uy, float uz) noexcept {
    const float cu = 3.0F *
        (static_cast<float>(Descriptor::cx(q)) * ux +
         static_cast<float>(Descriptor::cy(q)) * uy +
         static_cast<float>(Descriptor::cz(q)) * uz);
    const float uu = 1.5F * (ux * ux + uy * uy + uz * uz);
    return Descriptor::weight(q) * rho * (1.0F + cu + 0.5F * cu * cu - uu);
}

template<class Descriptor>
inline float guo_force(int q, float rho, float ux, float uy, float uz,
                       float ax, float ay, float az, float omega) noexcept {
    if (ax == 0.0F && ay == 0.0F && az == 0.0F) return 0.0F;
    const float cx = static_cast<float>(Descriptor::cx(q));
    const float cy = static_cast<float>(Descriptor::cy(q));
    const float cz = static_cast<float>(Descriptor::cz(q));
    const float fx = rho * ax;
    const float fy = rho * ay;
    const float fz = rho * az;
    const float cu = cx * ux + cy * uy + cz * uz;
    const float cf = cx * fx + cy * fy + cz * fz;
    const float uf = ux * fx + uy * fy + uz * fz;
    return Descriptor::weight(q) * (1.0F - 0.5F * omega) *
           (3.0F * (cf - uf) + 9.0F * cu * cf);
}

} // namespace detail

// Thread-safe, single-population-grid implementation of the published Esoteric-Pull
// streaming pattern. No FluidX3D source is used here; the access pattern follows the
// public algorithm description in Lehmann, Computation 2022, 10(6), 92.
template<class Descriptor>
class EsotericPullSolver {
public:
    static constexpr int q = Descriptor::q;
    static constexpr std::size_t pair_count = static_cast<std::size_t>((q - 1) / 2);

    explicit EsotericPullSolver(InPlaceLbmConfig config)
        : config_(normalized_config(config)),
          cells_(config_.nx * config_.ny * config_.nz),
          f_(cells_) {
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
            std::array<float, q> fin{};
            load_cell(n, fin);
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            accumulate_macroscopic(fin, rho, ux, uy, uz);
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
            std::array<float, q> fin{};
            load_cell(n, fin);
            float rho = 0.0F;
            for (float value : fin) rho += value;
            return rho;
        });
    }

    [[nodiscard]] double kinetic_energy() const {
        return 0.5 * cfd::core::parallel_sum(cells_, [&](std::size_t n) {
            std::array<float, q> fin{};
            load_cell(n, fin);
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            accumulate_macroscopic(fin, rho, ux, uy, uz);
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;
            return static_cast<double>(rho) *
                   (static_cast<double>(ux) * ux +
                    static_cast<double>(uy) * uy +
                    static_cast<double>(uz) * uz);
        });
    }

    [[nodiscard]] float max_speed() const {
        return static_cast<float>(cfd::core::parallel_max(cells_, [&](std::size_t n) {
            std::array<float, q> fin{};
            load_cell(n, fin);
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            accumulate_macroscopic(fin, rho, ux, uy, uz);
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;
            return std::sqrt(ux * ux + uy * uy + uz * uz);
        }));
    }

    [[nodiscard]] std::vector<float> logical_populations() const {
        std::vector<float> result(static_cast<std::size_t>(q) * cells_);
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            std::array<float, q> fin{};
            load_cell(n, fin);
            for (int d = 0; d < q; ++d) {
                result[static_cast<std::size_t>(d) * cells_ + n] = fin[static_cast<std::size_t>(d)];
            }
        });
        return result;
    }

    [[nodiscard]] const InPlaceLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }
    [[nodiscard]] std::size_t population_bytes() const noexcept {
        return static_cast<std::size_t>(q) * cells_ * sizeof(float);
    }
    [[nodiscard]] std::size_t two_lattice_population_bytes() const noexcept {
        return 2U * population_bytes();
    }
    [[nodiscard]] const cfd::core::StaticSoA<float, static_cast<std::size_t>(q)>& raw_storage() const noexcept {
        return f_;
    }

private:
    InPlaceLbmConfig config_;
    std::size_t cells_{};
    cfd::core::StaticSoA<float, static_cast<std::size_t>(q)> f_;
    std::uint64_t time_step_{0};

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

    [[nodiscard]] std::size_t neighbor_index(std::size_t x, std::size_t y, std::size_t z,
                                             int direction) const noexcept {
        const std::size_t xn = detail::periodic_shift(x, Descriptor::cx(direction), config_.nx);
        const std::size_t yn = detail::periodic_shift(y, Descriptor::cy(direction), config_.ny);
        const std::size_t zn = detail::periodic_shift(z, Descriptor::cz(direction), config_.nz);
        return (zn * config_.ny + yn) * config_.nx + xn;
    }

    template<class Initializer>
    void initialize(Initializer&& initializer) {
        time_step_ = 0;
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;
            const auto state = initializer(x, y, z);
            const float rho = state[0];
            const float ux = state[1];
            const float uy = state[2];
            const float uz = state[3];

            f_(0, n) = detail::equilibrium<Descriptor>(0, rho, ux, uy, uz);
            for (int i = 1; i < q; i += 2) {
                const std::size_t neighbor = neighbor_index(x, y, z, i + 1);
                f_(static_cast<std::size_t>(i), neighbor) =
                    detail::equilibrium<Descriptor>(i, rho, ux, uy, uz);
                f_(static_cast<std::size_t>(i + 1), n) =
                    detail::equilibrium<Descriptor>(i + 1, rho, ux, uy, uz);
            }
        });
    }

    void load_cell(std::size_t n, std::array<float, q>& fin) const noexcept {
        const std::size_t x = n % config_.nx;
        const std::size_t yz = n / config_.nx;
        const std::size_t y = yz % config_.ny;
        const std::size_t z = yz / config_.ny;
        load_cell(n, x, y, z, fin);
    }

    void load_cell(std::size_t n, std::size_t x, std::size_t y, std::size_t z,
                   std::array<float, q>& fin) const noexcept {
        std::array<std::size_t, pair_count> neighbors{};
        load_cell(n, x, y, z, fin, neighbors);
    }

    void load_cell(std::size_t n, std::size_t x, std::size_t y, std::size_t z,
                   std::array<float, q>& fin,
                   std::array<std::size_t, pair_count>& neighbors) const noexcept {
        const bool odd = (time_step_ & 1U) != 0U;
        fin[0] = f_(0, n);
        std::size_t pair = 0;
        for (int i = 1; i < q; i += 2, ++pair) {
            const std::size_t neighbor = neighbor_index(x, y, z, i + 1);
            neighbors[pair] = neighbor;
            if (odd) {
                fin[static_cast<std::size_t>(i)] = f_(static_cast<std::size_t>(i + 1), neighbor);
                fin[static_cast<std::size_t>(i + 1)] = f_(static_cast<std::size_t>(i), n);
            } else {
                fin[static_cast<std::size_t>(i)] = f_(static_cast<std::size_t>(i), neighbor);
                fin[static_cast<std::size_t>(i + 1)] = f_(static_cast<std::size_t>(i + 1), n);
            }
        }
    }

    void store_cell(std::size_t n,
                    const std::array<std::size_t, pair_count>& neighbors,
                    const std::array<float, q>& fout) noexcept {
        const bool odd = (time_step_ & 1U) != 0U;
        f_(0, n) = fout[0];
        std::size_t pair = 0;
        for (int i = 1; i < q; i += 2, ++pair) {
            const std::size_t neighbor = neighbors[pair];
            if (odd) {
                f_(static_cast<std::size_t>(i), n) = fout[static_cast<std::size_t>(i)];
                f_(static_cast<std::size_t>(i + 1), neighbor) = fout[static_cast<std::size_t>(i + 1)];
            } else {
                f_(static_cast<std::size_t>(i + 1), n) = fout[static_cast<std::size_t>(i)];
                f_(static_cast<std::size_t>(i), neighbor) = fout[static_cast<std::size_t>(i + 1)];
            }
        }
    }

    static void accumulate_macroscopic(const std::array<float, q>& fin,
                                       float& rho, float& ux, float& uy, float& uz) noexcept {
        rho = 0.0F;
        ux = 0.0F;
        uy = 0.0F;
        uz = 0.0F;
        for (int d = 0; d < q; ++d) {
            const float value = fin[static_cast<std::size_t>(d)];
            rho += value;
            ux += value * static_cast<float>(Descriptor::cx(d));
            uy += value * static_cast<float>(Descriptor::cy(d));
            uz += value * static_cast<float>(Descriptor::cz(d));
        }
        if (rho > 0.0F) {
            ux /= rho;
            uy /= rho;
            uz /= rho;
        } else {
            ux = 0.0F;
            uy = 0.0F;
            uz = 0.0F;
        }
    }

    void step_once() {
        const float omega = 1.0F / config_.tau;
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;

            std::array<float, q> fout{};
            std::array<std::size_t, pair_count> neighbors{};
            load_cell(n, x, y, z, fout, neighbors);

            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
            accumulate_macroscopic(fout, rho, ux, uy, uz);
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;

            for (int d = 0; d < q; ++d) {
                const std::size_t index = static_cast<std::size_t>(d);
                const float feq = detail::equilibrium<Descriptor>(d, rho, ux, uy, uz);
                fout[index] -= omega * (fout[index] - feq);
                fout[index] += detail::guo_force<Descriptor>(
                    d, rho, ux, uy, uz,
                    config_.acceleration_x, config_.acceleration_y, config_.acceleration_z, omega);
            }
            store_cell(n, neighbors, fout);
        });
        ++time_step_;
    }
};

using D2Q9InPlaceSolver = EsotericPullSolver<D2Q9InPlaceDescriptor>;
using D3Q19Solver = EsotericPullSolver<D3Q19Descriptor>;
using D3Q27Solver = EsotericPullSolver<D3Q27Descriptor>;

extern template class EsotericPullSolver<D2Q9InPlaceDescriptor>;
extern template class EsotericPullSolver<D3Q19Descriptor>;
extern template class EsotericPullSolver<D3Q27Descriptor>;

} // namespace cfd::lbm
