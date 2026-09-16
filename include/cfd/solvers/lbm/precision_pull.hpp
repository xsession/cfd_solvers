#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/core/soa_field.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace cfd::lbm {

template<class Real>
struct TypedLbmConfig {
    std::size_t nx{128};
    std::size_t ny{128};
    std::size_t nz{128};
    Real tau{static_cast<Real>(0.6)};
};

template<class Real>
struct TypedMacroscopicFields {
    cfd::core::AlignedVector<Real> rho;
    cfd::core::AlignedVector<Real> ux;
    cfd::core::AlignedVector<Real> uy;
    cfd::core::AlignedVector<Real> uz;
};

namespace detail {

template<class Descriptor, class Real>
inline Real equilibrium_typed(int q, Real rho, Real ux, Real uy, Real uz) noexcept {
    const Real cx = static_cast<Real>(Descriptor::cx(q));
    const Real cy = static_cast<Real>(Descriptor::cy(q));
    const Real cz = static_cast<Real>(Descriptor::cz(q));
    const Real cu = static_cast<Real>(3) * (cx * ux + cy * uy + cz * uz);
    const Real uu = static_cast<Real>(1.5) * (ux * ux + uy * uy + uz * uz);
    return static_cast<Real>(Descriptor::weight(q)) * rho *
           (static_cast<Real>(1) + cu + static_cast<Real>(0.5) * cu * cu - uu);
}

} // namespace detail

template<class Descriptor, class Real>
class PrecisionPullSolver {
    static_assert(std::is_floating_point_v<Real>);
public:
    static constexpr int q = Descriptor::q;

    explicit PrecisionPullSolver(TypedLbmConfig<Real> config)
        : config_(normalized_config(config)),
          cells_(config_.nx * config_.ny * config_.nz),
          f_(cells_),
          next_(cells_) {
        if (config_.tau <= static_cast<Real>(0.5)) {
            throw std::invalid_argument("LBM tau must be > 0.5 for positive viscosity");
        }
        initialize_uniform();
    }

    void initialize_uniform(Real rho = static_cast<Real>(1),
                            Real ux = static_cast<Real>(0),
                            Real uy = static_cast<Real>(0),
                            Real uz = static_cast<Real>(0)) {
        initialize([=](std::size_t, std::size_t, std::size_t) {
            return std::array<Real, 4>{rho, ux, uy, uz};
        });
    }

    void initialize_taylor_green(Real amplitude = static_cast<Real>(0.03)) {
        const Real two_pi = static_cast<Real>(2) * std::numbers::pi_v<Real>;
        initialize([=, this](std::size_t x, std::size_t y, std::size_t z) {
            const Real xf = (static_cast<Real>(x) + static_cast<Real>(0.5)) / static_cast<Real>(config_.nx);
            const Real yf = (static_cast<Real>(y) + static_cast<Real>(0.5)) / static_cast<Real>(config_.ny);
            const Real sx = std::sin(two_pi * xf);
            const Real cx = std::cos(two_pi * xf);
            const Real sy = std::sin(two_pi * yf);
            const Real cy = std::cos(two_pi * yf);
            Real z_factor = static_cast<Real>(1);
            if constexpr (Descriptor::dimensions == 3) {
                const Real zf = (static_cast<Real>(z) + static_cast<Real>(0.5)) / static_cast<Real>(config_.nz);
                z_factor = std::cos(two_pi * zf);
            }
            return std::array<Real, 4>{
                static_cast<Real>(1),
                amplitude * sx * cy * z_factor,
                -amplitude * cx * sy * z_factor,
                static_cast<Real>(0)};
        });
    }

    void step(std::size_t count = 1) {
        for (std::size_t i = 0; i < count; ++i) step_once();
    }

    [[nodiscard]] TypedMacroscopicFields<Real> compute_macroscopic() const {
        TypedMacroscopicFields<Real> result{
            cfd::core::AlignedVector<Real>(cells_, static_cast<Real>(0)),
            cfd::core::AlignedVector<Real>(cells_, static_cast<Real>(0)),
            cfd::core::AlignedVector<Real>(cells_, static_cast<Real>(0)),
            cfd::core::AlignedVector<Real>(cells_, static_cast<Real>(0))};
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            Real rho = static_cast<Real>(0);
            Real mx = static_cast<Real>(0);
            Real my = static_cast<Real>(0);
            Real mz = static_cast<Real>(0);
            for (int d = 0; d < q; ++d) {
                const Real value = f_(static_cast<std::size_t>(d), n);
                rho += value;
                mx += value * static_cast<Real>(Descriptor::cx(d));
                my += value * static_cast<Real>(Descriptor::cy(d));
                mz += value * static_cast<Real>(Descriptor::cz(d));
            }
            result.rho[n] = rho;
            if (rho > static_cast<Real>(0)) {
                result.ux[n] = mx / rho;
                result.uy[n] = my / rho;
                result.uz[n] = mz / rho;
            }
        });
        return result;
    }

    [[nodiscard]] long double mass() const {
        long double sum = 0.0L;
        for (std::size_t n = 0; n < cells_; ++n) {
            for (int d = 0; d < q; ++d) {
                sum += static_cast<long double>(f_(static_cast<std::size_t>(d), n));
            }
        }
        return sum;
    }

private:
    TypedLbmConfig<Real> config_;
    std::size_t cells_{};
    cfd::core::StaticSoA<Real, static_cast<std::size_t>(q)> f_;
    cfd::core::StaticSoA<Real, static_cast<std::size_t>(q)> next_;

    static TypedLbmConfig<Real> normalized_config(TypedLbmConfig<Real> config) {
        if (config.nx < 2 || config.ny < 2) throw std::invalid_argument("LBM grid requires nx, ny >= 2");
        if constexpr (Descriptor::dimensions == 2) config.nz = 1;
        else if (config.nz < 2) throw std::invalid_argument("3-D LBM grid requires nz >= 2");
        return config;
    }

    [[nodiscard]] std::size_t source_index(std::size_t x, std::size_t y, std::size_t z,
                                           int direction) const noexcept {
        auto shift = [](std::size_t c, int delta, std::size_t extent) noexcept {
            if (delta > 0) return c + 1U == extent ? 0U : c + 1U;
            if (delta < 0) return c == 0U ? extent - 1U : c - 1U;
            return c;
        };
        const std::size_t xs = shift(x, -Descriptor::cx(direction), config_.nx);
        const std::size_t ys = shift(y, -Descriptor::cy(direction), config_.ny);
        const std::size_t zs = shift(z, -Descriptor::cz(direction), config_.nz);
        return (zs * config_.ny + ys) * config_.nx + xs;
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
                const Real value = detail::equilibrium_typed<Descriptor, Real>(
                    d, state[0], state[1], state[2], state[3]);
                f_(static_cast<std::size_t>(d), n) = value;
                next_(static_cast<std::size_t>(d), n) = value;
            }
        });
    }

    void step_once() {
        const Real omega = static_cast<Real>(1) / config_.tau;
        cfd::core::parallel_for(cells_, [&](std::size_t n) {
            const std::size_t x = n % config_.nx;
            const std::size_t yz = n / config_.nx;
            const std::size_t y = yz % config_.ny;
            const std::size_t z = yz / config_.ny;
            std::array<Real, q> fin{};
            Real rho = static_cast<Real>(0);
            Real mx = static_cast<Real>(0);
            Real my = static_cast<Real>(0);
            Real mz = static_cast<Real>(0);
            for (int d = 0; d < q; ++d) {
                const Real value = f_(static_cast<std::size_t>(d), source_index(x, y, z, d));
                fin[static_cast<std::size_t>(d)] = value;
                rho += value;
                mx += value * static_cast<Real>(Descriptor::cx(d));
                my += value * static_cast<Real>(Descriptor::cy(d));
                mz += value * static_cast<Real>(Descriptor::cz(d));
            }
            const Real ux = rho > static_cast<Real>(0) ? mx / rho : static_cast<Real>(0);
            const Real uy = rho > static_cast<Real>(0) ? my / rho : static_cast<Real>(0);
            const Real uz = rho > static_cast<Real>(0) ? mz / rho : static_cast<Real>(0);
            for (int d = 0; d < q; ++d) {
                const std::size_t index = static_cast<std::size_t>(d);
                const Real feq = detail::equilibrium_typed<Descriptor, Real>(d, rho, ux, uy, uz);
                next_(index, n) = fin[index] - omega * (fin[index] - feq);
            }
        });
        std::swap(f_, next_);
    }
};

} // namespace cfd::lbm
