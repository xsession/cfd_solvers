#include "cfd/multiphysics/pipeline.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using cfd::multiphysics::CfdThermalStructuralOpticsPipeline;
    CfdThermalStructuralOpticsPipeline pipeline;
    pipeline.cfd_to_thermal = [](const auto& field) {
        CfdThermalStructuralOpticsPipeline::Field out;
        for (double value : field)
            out.push_back(value + 10.0);
        return out;
    };
    pipeline.thermal_to_structural = [](const auto& field) {
        CfdThermalStructuralOpticsPipeline::Field out;
        for (double value : field)
            out.push_back(2.0 * value);
        return out;
    };
    pipeline.structural_to_optics = [](const auto& field) {
        CfdThermalStructuralOpticsPipeline::Field out;
        for (double value : field)
            out.push_back(0.5 * value);
        return out;
    };

    const auto result = cfd::multiphysics::run_cfd_thermal_structural_optics(pipeline, {1.0, 2.0});
    assert(result.thermal.size() == 2 && result.structural.size() == 2 && result.optics.size() == 2);
    assert(std::abs(result.thermal[0] - 11.0) < 1e-12);
    assert(std::abs(result.structural[1] - 24.0) < 1e-12);
    assert(std::abs(result.optics[0] - 11.0) < 1e-12);

    bool rejected = false;
    try {
        (void)cfd::multiphysics::run_cfd_thermal_structural_optics(pipeline, {});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    pipeline.structural_to_optics = [](const auto&) { return CfdThermalStructuralOpticsPipeline::Field{NAN}; };
    rejected = false;
    try {
        (void)cfd::multiphysics::run_cfd_thermal_structural_optics(pipeline, {1.0});
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);

    cfd::multiphysics::CorrosionCfdStructuralPipeline corrosion;
    corrosion.recession_to_mesh = [](const auto& field) {
        return CfdThermalStructuralOpticsPipeline::Field{field[0] + 1.0};
    };
    corrosion.mesh_to_cfd = [](const auto& field) { return CfdThermalStructuralOpticsPipeline::Field{2.0 * field[0]}; };
    corrosion.cfd_to_structural = [](const auto& field) {
        return CfdThermalStructuralOpticsPipeline::Field{field[0] - 1.0};
    };
    const auto corrosion_result = cfd::multiphysics::run_corrosion_recession_mesh_cfd_structural(corrosion, {0.5});
    assert(std::abs(corrosion_result.structural[0] - 2.0) < 1.0e-12);

    std::cout << "CFD/thermal/structural/optics pipeline regression passed\n";
    return 0;
}
