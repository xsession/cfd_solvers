#include "cfd/multiphysics/pipeline.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {

PipelineResult run_cfd_thermal_structural_optics(const CfdThermalStructuralOpticsPipeline& pipeline,
                                                 const CfdThermalStructuralOpticsPipeline::Field& cfd_field) {
    if (!pipeline.cfd_to_thermal || !pipeline.thermal_to_structural || !pipeline.structural_to_optics ||
        cfd_field.empty()) {
        throw std::invalid_argument("incomplete multiphysics pipeline");
    }
    PipelineResult result;
    result.thermal = pipeline.cfd_to_thermal(cfd_field);
    result.structural = pipeline.thermal_to_structural(result.thermal);
    result.optics = pipeline.structural_to_optics(result.structural);
    for (const auto* field : {&result.thermal, &result.structural, &result.optics}) {
        if (field->empty())
            throw std::runtime_error("multiphysics pipeline stage returned an empty field");
        for (double value : *field)
            if (!std::isfinite(value))
                throw std::runtime_error("multiphysics pipeline produced non-finite data");
    }
    return result;
}

CorrosionCfdStructuralResult
run_corrosion_recession_mesh_cfd_structural(const CorrosionCfdStructuralPipeline& pipeline,
                                            const CorrosionCfdStructuralPipeline::Field& recession) {
    if (!pipeline.recession_to_mesh || !pipeline.mesh_to_cfd || !pipeline.cfd_to_structural || recession.empty()) {
        throw std::invalid_argument("incomplete corrosion/CFD/structural pipeline");
    }
    CorrosionCfdStructuralResult result;
    result.mesh = pipeline.recession_to_mesh(recession);
    result.cfd = pipeline.mesh_to_cfd(result.mesh);
    result.structural = pipeline.cfd_to_structural(result.cfd);
    for (const auto* field : {&result.mesh, &result.cfd, &result.structural}) {
        if (field->empty())
            throw std::runtime_error("corrosion pipeline stage returned an empty field");
        for (double value : *field) {
            if (!std::isfinite(value))
                throw std::runtime_error("corrosion pipeline produced non-finite data");
        }
    }
    return result;
}

} // namespace cfd::multiphysics
