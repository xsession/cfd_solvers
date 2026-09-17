#pragma once

#include "cfd/core/device_residency.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace cfd::lbm {

class D2Q9SyclSolver {
public:
    explicit D2Q9SyclSolver(D2Q9Config config);
    ~D2Q9SyclSolver();

    D2Q9SyclSolver(const D2Q9SyclSolver&) = delete;
    D2Q9SyclSolver& operator=(const D2Q9SyclSolver&) = delete;

    // Initialization executes directly on the selected device. No full lattice
    // is materialized on the host for these common initial conditions.
    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F);
    void initialize_taylor_green(float amplitude = 0.03F);
    void upload_solid_mask(std::span<const std::uint8_t> mask);
    void step(std::size_t count = 1);
    void wait();

    // Diagnostic extraction computes rho on-device and transfers only one
    // scalar per cell instead of all nine populations.
    [[nodiscard]] std::vector<float> download_density() const;
    [[nodiscard]] double total_mass() const;
    [[nodiscard]] std::string device_name() const;
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }
    [[nodiscard]] std::size_t resident_bytes() const noexcept;

private:
    D2Q9Config config_;
    std::size_t cells_{};
    mutable sycl::queue queue_;
    float* f_{nullptr};
    float* next_{nullptr};
    std::uint8_t* solid_{nullptr};
    mutable float* rho_scratch_{nullptr};
    mutable float* reduction_scalar_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};

    void enqueue_density(float* output) const;
};

} // namespace cfd::lbm
#endif
