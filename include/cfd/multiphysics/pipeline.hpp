#pragma once

#include <functional>
#include <vector>

namespace cfd::multiphysics {

struct CfdThermalStructuralOpticsPipeline {
    using Field = std::vector<double>;
    using Stage = std::function<Field(const Field&)>;

    Stage cfd_to_thermal;
    Stage thermal_to_structural;
    Stage structural_to_optics;
};

struct PipelineResult {
    CfdThermalStructuralOpticsPipeline::Field thermal;
    CfdThermalStructuralOpticsPipeline::Field structural;
    CfdThermalStructuralOpticsPipeline::Field optics;
};

struct CorrosionCfdStructuralPipeline {
    using Field = CfdThermalStructuralOpticsPipeline::Field;
    using Stage = CfdThermalStructuralOpticsPipeline::Stage;
    Stage recession_to_mesh;
    Stage mesh_to_cfd;
    Stage cfd_to_structural;
};

struct CorrosionCfdStructuralResult {
    CorrosionCfdStructuralPipeline::Field mesh;
    CorrosionCfdStructuralPipeline::Field cfd;
    CorrosionCfdStructuralPipeline::Field structural;
};

[[nodiscard]] PipelineResult
run_cfd_thermal_structural_optics(const CfdThermalStructuralOpticsPipeline& pipeline,
                                  const CfdThermalStructuralOpticsPipeline::Field& cfd_field);

[[nodiscard]] CorrosionCfdStructuralResult
run_corrosion_recession_mesh_cfd_structural(const CorrosionCfdStructuralPipeline& pipeline,
                                            const CorrosionCfdStructuralPipeline::Field& recession);

} // namespace cfd::multiphysics
