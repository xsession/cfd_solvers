#include "cfd/multiphysics/bioheat.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

cfd::multiphysics::VoxelTissueGrid3D sample_grid(double field) {
    using namespace cfd::multiphysics;
    VoxelTissueGrid3D grid;
    grid.nx = 4U;
    grid.ny = 4U;
    grid.nz = 2U;
    grid.dx_m = grid.dy_m = grid.dz_m = 0.01;
    grid.materials.push_back({1000.0, 3600.0, 0.5, 1.0, 0.01, 0.0});
    grid.material_index.assign(grid.nx * grid.ny * grid.nz, 0U);
    grid.electric_rms_v_per_m.assign(grid.material_index.size(), field);
    return grid;
}

void test_feedback() {
    using namespace cfd::multiphysics;
    const auto cold = temperature_feedback({}, 300.15);
    const auto hot = temperature_feedback({}, 330.15);
    require(hot.conductivity_s_per_m > cold.conductivity_s_per_m, "temperature feedback must change conductivity");
    require(hot.perfusion_per_s > cold.perfusion_per_s, "temperature feedback must change perfusion");
}

void test_exposure_scenarios() {
    using namespace cfd::multiphysics;
    BioelectromagneticExposureCase implant{"implant_fixture", ExposureScenario::implant, sample_grid(1.0), 0.01, 2.0};
    const auto implant_result = validate_exposure_case(implant);
    require(implant_result.maximum_mass_averaged_sar_w_per_kg > 0.0 && implant_result.within_limit,
            "implant SAR fixture should pass its limit");

    BioelectromagneticExposureCase wearable{"wearable_fixture", ExposureScenario::wearable, sample_grid(100.0), 0.01,
                                            2.0};
    require(!validate_exposure_case(wearable).within_limit, "wearable SAR fixture should expose an over-limit case");
}

} // namespace

int main() {
    try {
        test_feedback();
        test_exposure_scenarios();
        std::cout << "bioelectromagnetic feedback/exposure regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "bioelectromagnetic feedback regression failed: " << error.what() << '\n';
        return 1;
    }
}
