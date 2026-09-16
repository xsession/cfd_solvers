#include "cfd/solvers/lbm/d2q9_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace cfd::lbm {

D2Q9SyclSolver::D2Q9SyclSolver(D2Q9Config config)
    : config_(config), cells_(config.nx * config.ny), queue_(sycl::default_selector_v, sycl::property::queue::in_order{}) {
    if (config_.nx < 2 || config_.ny < 2 || config_.tau <= 0.5F) {
        throw std::invalid_argument("invalid D2Q9 SYCL configuration");
    }
    f_ = sycl::malloc_device<float>(D2Q9Descriptor::q * cells_, queue_);
    next_ = sycl::malloc_device<float>(D2Q9Descriptor::q * cells_, queue_);
    solid_ = sycl::malloc_device<std::uint8_t>(cells_, queue_);
    if (!f_ || !next_ || !solid_) throw std::bad_alloc{};
    queue_.memset(solid_, 0, cells_ * sizeof(std::uint8_t)).wait();
    initialize_uniform();
}

D2Q9SyclSolver::~D2Q9SyclSolver() {
    queue_.wait_and_throw();
    if (f_) sycl::free(f_, queue_);
    if (next_) sycl::free(next_, queue_);
    if (solid_) sycl::free(solid_, queue_);
}

void D2Q9SyclSolver::upload_equilibrium(const std::vector<float>& rho,
                                        const std::vector<float>& ux,
                                        const std::vector<float>& uy) {
    std::vector<float> host(D2Q9Descriptor::q * cells_);
    for (std::size_t i = 0; i < cells_; ++i) {
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            host[static_cast<std::size_t>(q) * cells_ + i] =
                D2Q9Descriptor::equilibrium(q, rho[i], ux[i], uy[i]);
        }
    }
    queue_.memcpy(f_, host.data(), host.size() * sizeof(float)).wait();
}

void D2Q9SyclSolver::initialize_uniform(float rho, float ux, float uy) {
    upload_equilibrium(std::vector<float>(cells_, rho),
                       std::vector<float>(cells_, ux),
                       std::vector<float>(cells_, uy));
}

void D2Q9SyclSolver::initialize_taylor_green(float amplitude) {
    std::vector<float> rho(cells_, 1.0F);
    std::vector<float> ux(cells_, 0.0F);
    std::vector<float> uy(cells_, 0.0F);
    const float two_pi = 2.0F * std::numbers::pi_v<float>;
    for (std::size_t i = 0; i < cells_; ++i) {
        const std::size_t x = i % config_.nx;
        const std::size_t y = i / config_.nx;
        const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(config_.nx);
        const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(config_.ny);
        ux[i] = amplitude * std::sin(two_pi * xf) * std::cos(two_pi * yf);
        uy[i] = -amplitude * std::cos(two_pi * xf) * std::sin(two_pi * yf);
    }
    upload_equilibrium(rho, ux, uy);
}

void D2Q9SyclSolver::step(std::size_t count) {
    const std::size_t nx = config_.nx;
    const std::size_t ny = config_.ny;
    const std::size_t cells = cells_;
    const float omega = 1.0F / config_.tau;

    for (std::size_t iteration = 0; iteration < count; ++iteration) {
        float* in = f_;
        float* out = next_;
        std::uint8_t* solid = solid_;
        queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
            const std::size_t i = gid[0];
            const std::size_t x = i % nx;
            const std::size_t y = i / nx;

            if (solid[i] != 0U) {
                for (int q = 0; q < D2Q9Descriptor::q; ++q) {
                    const int qo = D2Q9Descriptor::opposite(q);
                    out[static_cast<std::size_t>(q) * cells + i] =
                        in[static_cast<std::size_t>(qo) * cells + i];
                }
                return;
            }

            float fin[D2Q9Descriptor::q]{};
            float rho = 0.0F;
            float ux = 0.0F;
            float uy = 0.0F;
            for (int q = 0; q < D2Q9Descriptor::q; ++q) {
                const int sx_raw = static_cast<int>(x) - D2Q9Descriptor::cx(q);
                const int sy_raw = static_cast<int>(y) - D2Q9Descriptor::cy(q);
                const std::size_t sx = static_cast<std::size_t>((sx_raw + static_cast<int>(nx)) % static_cast<int>(nx));
                const std::size_t sy = static_cast<std::size_t>((sy_raw + static_cast<int>(ny)) % static_cast<int>(ny));
                const std::size_t source = sy * nx + sx;
                const float value = solid[source] != 0U
                    ? in[static_cast<std::size_t>(D2Q9Descriptor::opposite(q)) * cells + i]
                    : in[static_cast<std::size_t>(q) * cells + source];
                fin[q] = value;
                rho += value;
                ux += value * static_cast<float>(D2Q9Descriptor::cx(q));
                uy += value * static_cast<float>(D2Q9Descriptor::cy(q));
            }
            if (rho > 0.0F) {
                ux /= rho;
                uy /= rho;
            } else {
                ux = 0.0F;
                uy = 0.0F;
            }
            for (int q = 0; q < D2Q9Descriptor::q; ++q) {
                const float feq = D2Q9Descriptor::equilibrium(q, rho, ux, uy);
                out[static_cast<std::size_t>(q) * cells + i] =
                    fin[q] - omega * (fin[q] - feq);
            }
        });
        std::swap(f_, next_);
    }
}

void D2Q9SyclSolver::wait() { queue_.wait_and_throw(); }

std::vector<float> D2Q9SyclSolver::download_density() const {
    std::vector<float> populations(D2Q9Descriptor::q * cells_);
    std::vector<float> rho(cells_, 0.0F);
    queue_.memcpy(populations.data(), f_, populations.size() * sizeof(float)).wait();
    for (std::size_t i = 0; i < cells_; ++i) {
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            rho[i] += populations[static_cast<std::size_t>(q) * cells_ + i];
        }
    }
    return rho;
}

std::string D2Q9SyclSolver::device_name() const {
    return queue_.get_device().get_info<sycl::info::device::name>();
}

} // namespace cfd::lbm
#endif
