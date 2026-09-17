#pragma once

#include "cfd/core/device_residency.hpp"
#include "cfd/multibody/math.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cfd::multibody {

struct ResidentDemParticle {
    Vec3 position{};
    Vec3 linear_velocity{};
    Vec3 angular_velocity{};
    double radius{0.5};
    double mass{1.0};
};

struct ResidentDemSyclConfig {
    std::size_t capacity{1024U};
    Vec3 domain_min{-1.0,-1.0,-1.0};
    Vec3 domain_max{ 1.0, 1.0, 1.0};
    double cell_size{0.1};
    double dt{1.0e-4};
    Vec3 gravity{0.0,-9.81,0.0};
    double normal_stiffness{1.0e5};
    double normal_damping{20.0};
    double tangential_stiffness{5.0e4};
    double tangential_damping{20.0};
    double friction{0.5};
    double rolling_resistance{0.01};
    double cohesion_force{};
    std::size_t history_slots_per_particle{24U};
};

struct ResidentDemDeviceView {
    std::size_t count{};
    float* x{}; float* y{}; float* z{};
    float* vx{}; float* vy{}; float* vz{};
    float* wx{}; float* wy{}; float* wz{};
    float* radius{}; float* inverse_mass{};
};

// Device-resident spherical DEM baseline. The ordinary step path performs no
// bulk host transfer: grid construction, neighbor bucketing, persistent
// Mindlin history, pair forces, rolling torque and integration all run on the
// supplied SYCL queue. Host traffic is restricted to explicit upload/download
// and scalar diagnostics.
class ResidentDemSycl {
public:
    ResidentDemSycl(sycl::queue& queue, ResidentDemSyclConfig config);
    ~ResidentDemSycl() noexcept;
    ResidentDemSycl(const ResidentDemSycl&) = delete;
    ResidentDemSycl& operator=(const ResidentDemSycl&) = delete;

    void upload(std::span<const ResidentDemParticle> particles);
    [[nodiscard]] std::vector<ResidentDemParticle> download() const;
    void build_neighbor_buckets();
    void step(std::size_t steps = 1U);

    [[nodiscard]] double total_kinetic_energy() const;
    [[nodiscard]] std::size_t particle_count() const noexcept { return count_; }
    [[nodiscard]] std::size_t cell_count() const noexcept { return cell_count_; }
    [[nodiscard]] std::size_t resident_bytes() const noexcept;
    [[nodiscard]] ResidentDemDeviceView device_view() noexcept;
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept { return transfer_stats_; }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

private:
    sycl::queue& queue_;
    ResidentDemSyclConfig config_{};
    std::size_t count_{};
    std::size_t nx_{},ny_{},nz_{},cell_count_{};
    std::uint32_t history_stamp_{1U};

    float *x_{},*y_{},*z_{},*vx_{},*vy_{},*vz_{},*wx_{},*wy_{},*wz_{};
    float *radius_{},*inverse_mass_{},*inverse_inertia_{};
    float *fx_{},*fy_{},*fz_{},*tx_{},*ty_{},*tz_{};
    std::uint32_t *particle_cell_{},*cell_counts_{},*cell_offsets_{},*cell_cursor_{},*sorted_indices_{};
    std::uint32_t *history_neighbor_{},*history_seen_{};
    float *history_tx_{},*history_ty_{},*history_tz_{};
    mutable float* reduction_scalar_{};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};

    void clear_forces_with_gravity();
    void apply_contacts();
    void integrate();
    void cleanup_history(std::uint32_t stamp);
    void allocate();
    void release() noexcept;
};

} // namespace cfd::multibody
#endif
