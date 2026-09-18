#pragma once

#if defined(CFD_HAS_SYCL)

#include "cfd/fvm/resident_sycl.hpp"
#include "cfd/multibody/resident_dem_sycl.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace cfd::multibody {

struct ResidentCfdDemSyclConfig {
    double fluid_density{1.0};
    double dynamic_viscosity{1.0e-3};
    double minimum_void_fraction{0.05};
    double void_drag_exponent{2.65};
    double saffman_coefficient{1.615};
    bool include_pressure_gradient{true};
    bool include_saffman_lift{true};
};

// Resident unresolved CFD<->DEM coupling baseline. Particle volume and reaction
// forces are projected to the nearest general-PolyMesh cell entirely on-device.
// The same particle force is added to ResidentDemSycl while its negative is
// accumulated as an acceleration source for the resident incompressible FVM
// momentum equation. No full field transfer is required in couple().
class ResidentCfdDemSyclCoupler {
public:
    ResidentCfdDemSyclCoupler(cfd::fvm::ResidentPolyMeshSycl& mesh,
                              ResidentDemSycl& dem,
                              ResidentCfdDemSyclConfig config = {});
    ~ResidentCfdDemSyclCoupler() noexcept;

    ResidentCfdDemSyclCoupler(const ResidentCfdDemSyclCoupler&) = delete;
    ResidentCfdDemSyclCoupler& operator=(const ResidentCfdDemSyclCoupler&) = delete;

    // dem.begin_distributed_contact_step() must be called before couple() if the
    // generated particle forces are to participate in the same DEM timestep.
    // fluid_acceleration_soa is optional; if present it is incremented, not reset.
    void couple(const double* fluid_velocity_soa,
                const double* pressure_gradient_soa = nullptr,
                const double* vorticity_soa = nullptr,
                double* fluid_acceleration_soa = nullptr);

    [[nodiscard]] const double* solid_volume_fraction_device() const noexcept { return solid_fraction_; }
    [[nodiscard]] const double* void_fraction_device() const noexcept { return void_fraction_; }
    [[nodiscard]] std::size_t cell_count() const noexcept { return mesh_.cell_count(); }
    [[nodiscard]] std::size_t resident_bytes() const noexcept;
    void download_solid_volume_fraction(std::vector<double>& out) const;
    void download_void_fraction(std::vector<double>& out) const;

private:
    cfd::fvm::ResidentPolyMeshSycl& mesh_;
    ResidentDemSycl& dem_;
    ResidentCfdDemSyclConfig config_{};
    float* alpha_accum_{};
    float* reaction_x_{};
    float* reaction_y_{};
    float* reaction_z_{};
    double* solid_fraction_{};
    double* void_fraction_{};

    void allocate();
    void release() noexcept;
};

} // namespace cfd::multibody
#endif
