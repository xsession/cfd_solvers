#include "cfd/solvers/lbm/d2q9_sycl.hpp"

#if defined(CFD_HAS_SYCL)

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

namespace cfd::lbm {
namespace {

inline float equilibrium_device(int q, float rho, float ux, float uy) {
    const float cx = static_cast<float>(D2Q9Descriptor::cx(q));
    const float cy = static_cast<float>(D2Q9Descriptor::cy(q));
    const float cu = 3.0F * (cx * ux + cy * uy);
    const float uu = 1.5F * (ux * ux + uy * uy);
    return D2Q9Descriptor::weight(q) * rho * (1.0F + cu + 0.5F * cu * cu - uu);
}

} // namespace

D2Q9SyclSolver::D2Q9SyclSolver(D2Q9Config config)
    : config_(config), cells_(config.nx * config.ny),
      queue_(sycl::default_selector_v, sycl::property::queue::in_order{}) {
    if (config_.nx < 2 || config_.ny < 2 || config_.tau <= 0.5F) {
        throw std::invalid_argument("invalid D2Q9 SYCL configuration");
    }
    f_ = sycl::malloc_device<float>(D2Q9Descriptor::q * cells_, queue_);
    next_ = sycl::malloc_device<float>(D2Q9Descriptor::q * cells_, queue_);
    solid_ = sycl::malloc_device<std::uint8_t>(cells_, queue_);
    rho_scratch_ = sycl::malloc_device<float>(cells_, queue_);
    reduction_scalar_ = sycl::malloc_shared<float>(1U, queue_);
    if (!f_ || !next_ || !solid_ || !rho_scratch_ || !reduction_scalar_) throw std::bad_alloc{};
    queue_.memset(solid_, 0, cells_ * sizeof(std::uint8_t)).wait();
    initialize_uniform();
}

D2Q9SyclSolver::~D2Q9SyclSolver() {
    try { queue_.wait_and_throw(); } catch (...) {}
    if (f_) sycl::free(f_, queue_);
    if (next_) sycl::free(next_, queue_);
    if (solid_) sycl::free(solid_, queue_);
    if (rho_scratch_) sycl::free(rho_scratch_, queue_);
    if (reduction_scalar_) sycl::free(reduction_scalar_, queue_);
}

void D2Q9SyclSolver::initialize_uniform(float rho, float ux, float uy) {
    if (!(rho > 0.0F) || !std::isfinite(rho) || !std::isfinite(ux) || !std::isfinite(uy)) {
        throw std::invalid_argument("D2Q9 uniform state must be finite with positive density");
    }
    const std::size_t cells = cells_;
    float* storage = f_;
    queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
        const std::size_t i = gid[0];
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            storage[static_cast<std::size_t>(q) * cells + i] = equilibrium_device(q, rho, ux, uy);
        }
    });
}

void D2Q9SyclSolver::initialize_taylor_green(float amplitude) {
    if (!std::isfinite(amplitude)) throw std::invalid_argument("Taylor-Green amplitude must be finite");
    const std::size_t nx = config_.nx;
    const std::size_t ny = config_.ny;
    const std::size_t cells = cells_;
    const float two_pi = 2.0F * std::numbers::pi_v<float>;
    float* storage = f_;
    queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
        const std::size_t i = gid[0];
        const std::size_t x = i % nx;
        const std::size_t y = i / nx;
        const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(nx);
        const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(ny);
        const float ux = amplitude * sycl::sin(two_pi * xf) * sycl::cos(two_pi * yf);
        const float uy = -amplitude * sycl::cos(two_pi * xf) * sycl::sin(two_pi * yf);
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            storage[static_cast<std::size_t>(q) * cells + i] = equilibrium_device(q, 1.0F, ux, uy);
        }
    });
}

void D2Q9SyclSolver::upload_solid_mask(std::span<const std::uint8_t> mask) {
    if (mask.size() != cells_) throw std::invalid_argument("D2Q9 solid mask size mismatch");
    const std::size_t bytes = cells_ * sizeof(std::uint8_t);
    queue_.memcpy(solid_, mask.data(), bytes).wait();
    transfer_stats_.record_host_to_device(bytes);
    transfer_stats_.record_synchronization();
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
                const float feq = equilibrium_device(q, rho, ux, uy);
                out[static_cast<std::size_t>(q) * cells + i] = fin[q] - omega * (fin[q] - feq);
            }
        });
        std::swap(f_, next_);
    }
}

void D2Q9SyclSolver::wait() {
    queue_.wait_and_throw();
    transfer_stats_.record_synchronization();
}

void D2Q9SyclSolver::enqueue_density(float* output) const {
    const std::size_t cells = cells_;
    const float* populations = f_;
    queue_.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> gid) {
        const std::size_t i = gid[0];
        float rho = 0.0F;
        for (int q = 0; q < D2Q9Descriptor::q; ++q) {
            rho += populations[static_cast<std::size_t>(q) * cells + i];
        }
        output[i] = rho;
    });
}

std::vector<float> D2Q9SyclSolver::download_density() const {
    std::vector<float> rho(cells_, 0.0F);
    enqueue_density(rho_scratch_);
    const std::size_t bytes = cells_ * sizeof(float);
    queue_.memcpy(rho.data(), rho_scratch_, bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(bytes);
    transfer_stats_.record_synchronization();
    return rho;
}

double D2Q9SyclSolver::total_mass() const {
    *reduction_scalar_ = 0.0F;
    const float* populations = f_;
    const std::size_t count = static_cast<std::size_t>(D2Q9Descriptor::q) * cells_;
    auto reduction = sycl::reduction(reduction_scalar_, sycl::plus<float>());
    queue_.parallel_for(sycl::range<1>(count), reduction,
                        [=](sycl::id<1> gid, auto& sum) {
                            sum.combine(populations[gid[0]]);
                        }).wait_and_throw();
    transfer_stats_.record_synchronization();
    return static_cast<double>(*reduction_scalar_);
}

std::string D2Q9SyclSolver::device_name() const {
    return queue_.get_device().get_info<sycl::info::device::name>();
}

std::size_t D2Q9SyclSolver::resident_bytes() const noexcept {
    return (2U * static_cast<std::size_t>(D2Q9Descriptor::q) + 1U) * cells_ * sizeof(float)
           + cells_ * sizeof(std::uint8_t) + sizeof(float);
}

} // namespace cfd::lbm
#endif
