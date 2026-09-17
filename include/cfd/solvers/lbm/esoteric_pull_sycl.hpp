#pragma once

#include "cfd/core/device_residency.hpp"
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
        if (macro_scratch_) sycl::free(macro_scratch_, queue_);
        if (reduction_scalar_) sycl::free(reduction_scalar_, queue_);
        if (local_acceleration_) sycl::free(local_acceleration_, queue_);
    }

    EsotericPullSyclSolver(const EsotericPullSyclSolver&) = delete;
    EsotericPullSyclSolver& operator=(const EsotericPullSyclSolver&) = delete;

    // Common initial conditions are generated on-device so startup does not
    // require a q*N host population staging array.
    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F, float uz = 0.0F) {
        if (!(rho > 0.0F) || !std::isfinite(rho) || !std::isfinite(ux) ||
            !std::isfinite(uy) || !std::isfinite(uz)) {
            throw std::invalid_argument("LBM uniform state must be finite with positive density");
        }
        initialize_device([=](std::size_t, std::size_t, std::size_t) {
            return std::array<float, 4>{rho, ux, uy, uz};
        });
    }

    void initialize_taylor_green(float amplitude = 0.03F) {
        if (!std::isfinite(amplitude)) throw std::invalid_argument("Taylor-Green amplitude must be finite");
        const float two_pi = 2.0F * std::numbers::pi_v<float>;
        const std::size_t nx = config_.nx;
        const std::size_t ny = config_.ny;
        const std::size_t nz = config_.nz;
        initialize_device([=](std::size_t x, std::size_t y, std::size_t z) {
            const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(nx);
            const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(ny);
            const float sx = sycl::sin(two_pi * xf);
            const float cx = sycl::cos(two_pi * xf);
            const float sy = sycl::sin(two_pi * yf);
            const float cy = sycl::cos(two_pi * yf);
            float z_factor = 1.0F;
            if constexpr (Descriptor::dimensions == 3) {
                const float zf = (static_cast<float>(z) + 0.5F) / static_cast<float>(nz);
                z_factor = sycl::cos(two_pi * zf);
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
        const float ax0 = config_.acceleration_x;
        const float ay0 = config_.acceleration_y;
        const float az0 = config_.acceleration_z;
        const float* local_ax = local_acceleration_;
        const float* local_ay = local_acceleration_ ? local_acceleration_ + cells_ : nullptr;
        const float* local_az = local_acceleration_ ? local_acceleration_ + 2U * cells_ : nullptr;
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
                const float ax = ax0 + (local_ax ? local_ax[n] : 0.0F);
                const float ay = ay0 + (local_ay ? local_ay[n] : 0.0F);
                const float az = az0 + (local_az ? local_az[n] : 0.0F);
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

    void wait() {
        queue_.wait_and_throw();
        transfer_stats_.record_synchronization();
    }

    struct DeviceMacroscopicView {
        float* rho{};
        float* ux{};
        float* uy{};
        float* uz{};
        std::size_t cells{};
    };

    // Materialize macroscopic fields in device scratch memory without copying
    // them to the host. Downstream visualization/coupling kernels can consume
    // this view directly while the simulation remains GPU-resident.
    [[nodiscard]] DeviceMacroscopicView refresh_device_macroscopic() const {
        ensure_macro_scratch();
        const bool odd = (time_step_ & 1U) != 0U;
        const std::size_t nx = config_.nx;
        const std::size_t ny = config_.ny;
        const std::size_t nz = config_.nz;
        const std::size_t cells = cells_;
        const float ax0 = config_.acceleration_x;
        const float ay0 = config_.acceleration_y;
        const float az0 = config_.acceleration_z;
        const float* local_ax = local_acceleration_;
        const float* local_ay = local_acceleration_ ? local_acceleration_ + cells_ : nullptr;
        const float* local_az = local_acceleration_ ? local_acceleration_ + 2U * cells_ : nullptr;
        const float* storage = f_;
        float* rho_out = macro_scratch_;
        float* ux_out = macro_scratch_ + cells;
        float* uy_out = macro_scratch_ + 2U * cells;
        float* uz_out = macro_scratch_ + 3U * cells;
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
            const float ax = ax0 + (local_ax ? local_ax[n] : 0.0F);
            const float ay = ay0 + (local_ay ? local_ay[n] : 0.0F);
            const float az = az0 + (local_az ? local_az[n] : 0.0F);
            rho_out[n] = rho;
            ux_out[n] = ux + 0.5F * ax;
            uy_out[n] = uy + 0.5F * ay;
            uz_out[n] = uz + 0.5F * az;
        });
        return {rho_out, ux_out, uy_out, uz_out, cells};
    }

    [[nodiscard]] MacroscopicFields download_macroscopic() const {
        const auto view = refresh_device_macroscopic();
        MacroscopicFields fields{
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F),
            cfd::core::AlignedVector<float>(cells_, 0.0F)
        };
        queue_.memcpy(fields.rho.data(), view.rho, cells_ * sizeof(float));
        queue_.memcpy(fields.ux.data(), view.ux, cells_ * sizeof(float));
        queue_.memcpy(fields.uy.data(), view.uy, cells_ * sizeof(float));
        queue_.memcpy(fields.uz.data(), view.uz, cells_ * sizeof(float)).wait_and_throw();
        const std::size_t bytes = 4U * cells_ * sizeof(float);
        transfer_stats_.record_device_to_host(bytes);
        transfer_stats_.record_synchronization();
        return fields;
    }

    [[nodiscard]] double total_mass() const {
        ensure_reduction_scalar();
        *reduction_scalar_ = 0.0F;
        const float* storage = f_;
        const std::size_t count = static_cast<std::size_t>(q) * cells_;
        auto reduction = sycl::reduction(reduction_scalar_, sycl::plus<float>());
        queue_.parallel_for(sycl::range<1>(count), reduction,
                            [=](sycl::id<1> gid, auto& sum) {
                                sum.combine(storage[gid[0]]);
                            }).wait_and_throw();
        transfer_stats_.record_synchronization();
        return static_cast<double>(*reduction_scalar_);
    }


    struct DeviceAccelerationView {
        float* ax{};
        float* ay{};
        float* az{};
        std::size_t cells{};
    };

    [[nodiscard]] DeviceAccelerationView device_local_acceleration(bool clear = false) {
        ensure_local_acceleration();
        if (clear) clear_local_body_acceleration();
        return {local_acceleration_, local_acceleration_ + cells_, local_acceleration_ + 2U * cells_, cells_};
    }

    void clear_local_body_acceleration() {
        if (!local_acceleration_) return;
        queue_.memset(local_acceleration_, 0, 3U * cells_ * sizeof(float));
    }

    void set_local_body_acceleration(std::span<const float> ax,
                                     std::span<const float> ay,
                                     std::span<const float> az) {
        if (ax.size() != cells_ || ay.size() != cells_ || az.size() != cells_) {
            throw std::invalid_argument("local LBM acceleration field size mismatch");
        }
        ensure_local_acceleration();
        queue_.memcpy(local_acceleration_, ax.data(), cells_ * sizeof(float));
        queue_.memcpy(local_acceleration_ + cells_, ay.data(), cells_ * sizeof(float));
        queue_.memcpy(local_acceleration_ + 2U * cells_, az.data(), cells_ * sizeof(float)).wait_and_throw();
        transfer_stats_.record_host_to_device(3U * cells_ * sizeof(float));
        transfer_stats_.record_synchronization();
    }

    [[nodiscard]] bool has_local_body_acceleration() const noexcept { return local_acceleration_ != nullptr; }
    [[nodiscard]] sycl::queue& queue() noexcept { return queue_; }
    [[nodiscard]] const sycl::queue& queue() const noexcept { return queue_; }

    [[nodiscard]] std::string device_name() const {
        return queue_.get_device().get_info<sycl::info::device::name>();
    }
    [[nodiscard]] const InPlaceLbmConfig& config() const noexcept { return config_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }
    [[nodiscard]] std::uint64_t time_step() const noexcept { return time_step_; }
    [[nodiscard]] std::size_t population_bytes() const noexcept {
        return static_cast<std::size_t>(q) * cells_ * sizeof(float);
    }
    [[nodiscard]] std::size_t resident_bytes() const noexcept {
        return population_bytes() + (macro_scratch_ ? 4U * cells_ * sizeof(float) : 0U) +
               (reduction_scalar_ ? sizeof(float) : 0U) +
               (local_acceleration_ ? 3U * cells_ * sizeof(float) : 0U);
    }
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

private:
    InPlaceLbmConfig config_;
    std::size_t cells_{};
    mutable sycl::queue queue_;
    float* f_{nullptr};
    mutable float* macro_scratch_{nullptr};
    mutable float* reduction_scalar_{nullptr};
    float* local_acceleration_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};
    std::uint64_t time_step_{0};

    static InPlaceLbmConfig normalized_config(InPlaceLbmConfig config) {
        if (config.nx < 2 || config.ny < 2) throw std::invalid_argument("LBM grid requires nx, ny >= 2");
        if constexpr (Descriptor::dimensions == 2) config.nz = 1;
        else if (config.nz < 2) throw std::invalid_argument("3-D LBM grid requires nz >= 2");
        return config;
    }

    void ensure_macro_scratch() const {
        if (macro_scratch_) return;
        macro_scratch_ = sycl::malloc_device<float>(4U * cells_, queue_);
        if (!macro_scratch_) throw std::bad_alloc{};
    }

    void ensure_reduction_scalar() const {
        if (reduction_scalar_) return;
        reduction_scalar_ = sycl::malloc_shared<float>(1U, queue_);
        if (!reduction_scalar_) throw std::bad_alloc{};
    }


    void ensure_local_acceleration() {
        if (local_acceleration_) return;
        local_acceleration_ = sycl::malloc_device<float>(3U * cells_, queue_);
        if (!local_acceleration_) throw std::bad_alloc{};
        queue_.memset(local_acceleration_, 0, 3U * cells_ * sizeof(float));
    }

    template<class Initializer>
    void initialize_device(Initializer initializer) {
        time_step_ = 0;
        const std::size_t nx = config_.nx;
        const std::size_t ny = config_.ny;
        const std::size_t nz = config_.nz;
        const std::size_t cells = cells_;
        float* storage = f_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
            const std::size_t n = gid[0];
            const std::size_t x = n % nx;
            const std::size_t yz = n / nx;
            const std::size_t y = yz % ny;
            const std::size_t z = yz / ny;
            const auto state = initializer(x, y, z);
            const float rho = state[0];
            const float ux = state[1];
            const float uy = state[2];
            const float uz = state[3];
            storage[n] = detail::equilibrium<Descriptor>(0, rho, ux, uy, uz);
            for (int i = 1; i < q; i += 2) {
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
                storage[static_cast<std::size_t>(i) * cells + neighbor] =
                    detail::equilibrium<Descriptor>(i, rho, ux, uy, uz);
                storage[static_cast<std::size_t>(i + 1) * cells + n] =
                    detail::equilibrium<Descriptor>(i + 1, rho, ux, uy, uz);
            }
        });
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
