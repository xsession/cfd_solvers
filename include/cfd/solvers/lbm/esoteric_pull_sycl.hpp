#pragma once

#include "cfd/solvers/lbm/esoteric_pull.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <string>
#include <vector>

namespace cfd::lbm {

template<class Descriptor>
class EsotericPullSyclSolver {
public:
    static constexpr int q = Descriptor::q;

    explicit EsotericPullSyclSolver(InPlaceLbmConfig config)
        : config_(normalized_config(config)),
          cells_(config_.nx * config_.ny * config_.nz),
          queue_(sycl::default_selector_v, sycl::property::queue::in_order{}) {
        if (config_.tau <= 0.5F) {
            throw std::invalid_argument("LBM tau must be > 0.5 for positive viscosity");
        }
        f_ = sycl::malloc_device<float>(static_cast<std::size_t>(q) * cells_, queue_);
        if (!f_) throw std::bad_alloc{};
        initialize_uniform();
    }

    ~EsotericPullSyclSolver() noexcept {
        try {
            queue_.wait_and_throw();
        } catch (...) {
            // Destructors must not propagate asynchronous device exceptions.
        }
        if (f_) sycl::free(f_, queue_);
    }

    EsotericPullSyclSolver(const EsotericPullSyclSolver&) = delete;
    EsotericPullSyclSolver& operator=(const EsotericPullSyclSolver&) = delete;

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F, float uz = 0.0F) {
        initialize_host([=](std::size_t, std::size_t, std::size_t) {
            return std::array<float, 4>{rho, ux, uy, uz};
        });
    }

    void initialize_taylor_green(float amplitude = 0.03F) {
        const float two_pi = 2.0F * std::numbers::pi_v<float>;
        initialize_host([=, this](std::size_t x, std::size_t y, std::size_t z) {
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
        const std::size_t nx = config_.nx;
        const std::size_t ny = config_.ny;
        const std::size_t nz = config_.nz;
        const std::size_t cells = cells_;
        const float omega = 1.0F / config_.tau;
        const float ax = config_.acceleration_x;
        const float ay = config_.acceleration_y;
        const float az = config_.acceleration_z;
        float* storage = f_;

        for (std::size_t iteration = 0; iteration < count; ++iteration) {
            const bool odd = (time_step_ & 1U) != 0U;
            queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
                const std::size_t n = gid[0];
                const std::size_t x = n % nx;
                const std::size_t yz = n / nx;
                const std::size_t y = yz % ny;
                const std::size_t z = yz / ny;

                float local[Descriptor::q];
                local[0] = storage[n];

                for (int i = 1; i < Descriptor::q; i += 2) {
                    const int d = i + 1;
                    std::size_t xn = x;
                    std::size_t yn = y;
                    std::size_t zn = z;
                    const int dx = Descriptor::cx(d);
                    const int dy = Descriptor::cy(d);
                    const int dz = Descriptor::cz(d);
                    if (dx > 0) xn = x + 1U == nx ? 0U : x + 1U;
                    else if (dx < 0) xn = x == 0U ? nx - 1U : x - 1U;
                    if (dy > 0) yn = y + 1U == ny ? 0U : y + 1U;
                    else if (dy < 0) yn = y == 0U ? ny - 1U : y - 1U;
                    if (dz > 0) zn = z + 1U == nz ? 0U : z + 1U;
                    else if (dz < 0) zn = z == 0U ? nz - 1U : z - 1U;
                    const std::size_t neighbor = (zn * ny + yn) * nx + xn;

                    if (odd) {
                        local[i] = storage[static_cast<std::size_t>(i + 1) * cells + neighbor];
                        local[i + 1] = storage[static_cast<std::size_t>(i) * cells + n];
                    } else {
                        local[i] = storage[static_cast<std::size_t>(i) * cells + neighbor];
                        local[i + 1] = storage[static_cast<std::size_t>(i + 1) * cells + n];
                    }
                }

                float rho = 0.0F;
                float ux = 0.0F;
                float uy = 0.0F;
                float uz = 0.0F;
                for (int d = 0; d < Descriptor::q; ++d) {
                    const float value = local[d];
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
                ux += 0.5F * ax;
                uy += 0.5F * ay;
                uz += 0.5F * az;
                const float uu = 1.5F * (ux * ux + uy * uy + uz * uz);
                const float fx = rho * ax;
                const float fy = rho * ay;
                const float fz = rho * az;
                const float uf = ux * fx + uy * fy + uz * fz;
                for (int d = 0; d < Descriptor::q; ++d) {
                    const float cx = static_cast<float>(Descriptor::cx(d));
                    const float cy = static_cast<float>(Descriptor::cy(d));
                    const float cz = static_cast<float>(Descriptor::cz(d));
                    const float cu_raw = cx * ux + cy * uy + cz * uz;
                    const float cu = 3.0F * cu_raw;
                    const float feq = Descriptor::weight(d) * rho * (1.0F + cu + 0.5F * cu * cu - uu);
                    const float cf = cx * fx + cy * fy + cz * fz;
                    const float force = Descriptor::weight(d) * (1.0F - 0.5F * omega) *
                        (3.0F * (cf - uf) + 9.0F * cu_raw * cf);
                    local[d] -= omega * (local[d] - feq);
                    local[d] += force;
                }

                storage[n] = local[0];
                for (int i = 1; i < Descriptor::q; i += 2) {
                    const int d = i + 1;
                    std::size_t xn = x;
                    std::size_t yn = y;
                    std::size_t zn = z;
                    const int dx = Descriptor::cx(d);
                    const int dy = Descriptor::cy(d);
                    const int dz = Descriptor::cz(d);
                    if (dx > 0) xn = x + 1U == nx ? 0U : x + 1U;
                    else if (dx < 0) xn = x == 0U ? nx - 1U : x - 1U;
                    if (dy > 0) yn = y + 1U == ny ? 0U : y + 1U;
                    else if (dy < 0) yn = y == 0U ? ny - 1U : y - 1U;
                    if (dz > 0) zn = z + 1U == nz ? 0U : z + 1U;
                    else if (dz < 0) zn = z == 0U ? nz - 1U : z - 1U;
                    const std::size_t neighbor = (zn * ny + yn) * nx + xn;

                    if (odd) {
                        storage[static_cast<std::size_t>(i) * cells + n] = local[i];
                        storage[static_cast<std::size_t>(i + 1) * cells + neighbor] = local[i + 1];
                    } else {
                        storage[static_cast<std::size_t>(i + 1) * cells + n] = local[i];
                        storage[static_cast<std::size_t>(i) * cells + neighbor] = local[i + 1];
                    }
                }
            });
            ++time_step_;
        }
    }

    void wait() { queue_.wait_and_throw(); }

    [[nodiscard]] MacroscopicFields download_macroscopic() const {
        queue_.wait_and_throw();
        std::vector<float> storage(static_cast<std::size_t>(q) * cells_);
        queue_.memcpy(storage.data(), f_, storage.size() * sizeof(float)).wait();
        MacroscopicFields fields{
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F)
        };
        const bool odd = (time_step_ & 1U) != 0U;
        for (std::size_t n = 0; n < cells_; ++n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;
            std::array<float, q> fin{};
            fin[0] = storage[n];
            for (int i = 1; i < q; i += 2) {
                const std::size_t neighbor = host_neighbor(x, y, z, i + 1);
                if (odd) {
                    fin[static_cast<std::size_t>(i)] =
                        storage[static_cast<std::size_t>(i + 1) * cells_ + neighbor];
                    fin[static_cast<std::size_t>(i + 1)] =
                        storage[static_cast<std::size_t>(i) * cells_ + n];
                } else {
                    fin[static_cast<std::size_t>(i)] =
                        storage[static_cast<std::size_t>(i) * cells_ + neighbor];
                    fin[static_cast<std::size_t>(i + 1)] =
                        storage[static_cast<std::size_t>(i + 1) * cells_ + n];
                }
            }
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            float uz = 0.0F;
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
            }
            ux += 0.5F * config_.acceleration_x;
            uy += 0.5F * config_.acceleration_y;
            uz += 0.5F * config_.acceleration_z;
            fields.rho[n] = rho;
            fields.ux[n] = ux;
            fields.uy[n] = uy;
            fields.uz[n] = uz;
        }
        return fields;
    }

    [[nodiscard]] std::string device_name() const {
        return queue_.get_device().get_info<sycl::info::device::name>();
    }
    [[nodiscard]] const InPlaceLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }
    [[nodiscard]] std::size_t population_bytes() const noexcept {
        return static_cast<std::size_t>(q) * cells_ * sizeof(float);
    }

private:
    InPlaceLbmConfig config_;
    std::size_t cells_{};
    mutable sycl::queue queue_;
    float* f_{nullptr};
    std::uint64_t time_step_{0};

    static InPlaceLbmConfig normalized_config(InPlaceLbmConfig config) {
        if (config.nx < 2 || config.ny < 2) throw std::invalid_argument("LBM grid requires nx, ny >= 2");
        if constexpr (Descriptor::dimensions == 2) config.nz = 1;
        else if (config.nz < 2) throw std::invalid_argument("3-D LBM grid requires nz >= 2");
        return config;
    }

    [[nodiscard]] std::size_t host_neighbor(std::size_t x,
                                            std::size_t y,
                                            std::size_t z,
                                            int direction) const noexcept {
        const std::size_t xn = detail::periodic_shift(x, Descriptor::cx(direction), config_.nx);
        const std::size_t yn = detail::periodic_shift(y, Descriptor::cy(direction), config_.ny);
        const std::size_t zn = detail::periodic_shift(z, Descriptor::cz(direction), config_.nz);
        return (zn * config_.ny + yn) * config_.nx + xn;
    }

    template<class Initializer>
    void initialize_host(Initializer&& initializer) {
        time_step_ = 0;
        std::vector<float> storage(static_cast<std::size_t>(q) * cells_, 0.0F);
        for (std::size_t n = 0; n < cells_; ++n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;
            const auto state = initializer(x, y, z);
            const float rho = state[0];
            const float ux = state[1];
            const float uy = state[2];
            const float uz = state[3];
            storage[n] = detail::equilibrium<Descriptor>(0, rho, ux, uy, uz);
            for (int i = 1; i < q; i += 2) {
                const std::size_t neighbor = host_neighbor(x, y, z, i + 1);
                storage[static_cast<std::size_t>(i) * cells_ + neighbor] =
                    detail::equilibrium<Descriptor>(i, rho, ux, uy, uz);
                storage[static_cast<std::size_t>(i + 1) * cells_ + n] =
                    detail::equilibrium<Descriptor>(i + 1, rho, ux, uy, uz);
            }
        }
        queue_.memcpy(f_, storage.data(), storage.size() * sizeof(float)).wait();
    }
};

using D2Q9InPlaceSyclSolver = EsotericPullSyclSolver<D2Q9InPlaceDescriptor>;
using D3Q19SyclSolver = EsotericPullSyclSolver<D3Q19Descriptor>;
using D3Q27SyclSolver = EsotericPullSyclSolver<D3Q27Descriptor>;

extern template class EsotericPullSyclSolver<D2Q9InPlaceDescriptor>;
extern template class EsotericPullSyclSolver<D3Q19Descriptor>;
extern template class EsotericPullSyclSolver<D3Q27Descriptor>;

} // namespace cfd::lbm
#endif
