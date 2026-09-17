#include "cfd/core/device_residency.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#include "cfd/solvers/lbm/resident_multiphysics_sycl.hpp"
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void test_residency_contract_marker() {
    cfd::core::DeviceTransferStats stats;
    require(stats.host_transfer_bytes() == 0U, "hot-loop residency contract starts without host traffic");
}

#if defined(CFD_HAS_SYCL)
void test_resident_multiphysics_chain() {
    constexpr std::size_t nx = 10U, ny = 8U, nz = 6U, cells = nx * ny * nz;
    cfd::lbm::D3Q19SyclSolver flow({nx, ny, nz, 0.72F});
    flow.initialize_uniform(1.0F, 0.0F, 0.0F, 0.0F);
    auto force = flow.device_local_acceleration(true);
    auto macro = flow.refresh_device_macroscopic();

    cfd::lbm::ResidentThermalD3Q7Sycl thermal(flow.queue(), {nx, ny, nz, 0.02F});
    thermal.initialize_uniform(300.0F);
    const double heat0 = thermal.total_scalar();
    thermal.reset_transfer_stats();
    thermal.step(macro, 3U);
    thermal.add_boussinesq_acceleration(force, 300.0F, 3.0e-3F, {0.0, -1.0, 0.0});
    const double heat1 = thermal.total_scalar();
    require(std::abs(heat1 - heat0) / heat0 < 2.0e-5, "resident thermal lattice conserves a uniform scalar");
    require(thermal.transfer_stats().host_transfer_bytes() == 0U,
            "resident thermal stepping does not transfer full fields through the host");

    cfd::lbm::ResidentFreeSurfaceSycl3D surface(flow.queue(), nx, ny, nz);
    surface.initialize_half_space_x(3.5F);
    const double volume0 = surface.volume();
    surface.advect(macro, 0.2F);
    const double volume1 = surface.volume();
    require(std::abs(volume1 - volume0) < 2.0e-4 * static_cast<double>(cells),
            "resident zero-flow VOF preserves volume");
    surface.add_surface_tension_acceleration(force, 0.072F, 1.0F);

    cfd::lbm::ImmersedBoundaryParticle seed{{2.25, 3.5, 2.75}, {}, 1.0, 0.5};
    cfd::lbm::ResidentParticlesSycl3D particles(flow.queue(), nx, ny, nz, 1U);
    particles.upload(std::span<const cfd::lbm::ImmersedBoundaryParticle>(&seed, 1U));
    particles.reset_transfer_stats();
    particles.advance(macro, force, 0.05F, true);
    require(particles.transfer_stats().host_transfer_bytes() == 0U,
            "resident particle push and two-way force spreading stay on-device");

    cfd::lbm::ResidentFlowDiagnosticsSycl3D diagnostics(flow.queue(), nx, ny, nz);
    require(diagnostics.q_criterion(macro) != nullptr, "resident Q-criterion returns a device field");

    flow.reset_transfer_stats();
    flow.step(2U);
    flow.wait();
    require(flow.transfer_stats().host_transfer_bytes() == 0U,
            "coupled local-force LBM step retains the zero-host-field-transfer contract");
}
#endif
} // namespace

int main() {
    try {
        test_residency_contract_marker();
#if defined(CFD_HAS_SYCL)
        test_resident_multiphysics_chain();
#endif
        std::cout << "v0.15.6 GPU LBM multiphysics tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.6 GPU LBM multiphysics test failure: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
