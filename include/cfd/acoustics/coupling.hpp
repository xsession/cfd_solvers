#pragma once

#include "cfd/fem/mesh2d.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/solvers/fem/dynamic_bar1d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::acoustics {

struct StructuralNodalForce2D { double x{}; double y{}; };

// Integrate pressure traction -p*n over selected FEM boundary edges and split
// each constant edge load equally between its two nodes.
[[nodiscard]] std::vector<StructuralNodalForce2D> pressure_to_structural_boundary_load(
    const cfd::fem::Mesh2D& mesh,
    std::span<const double> edge_pressure_pa,
    int patch = -1);

[[nodiscard]] double acoustic_pressure_to_bar_tip_force(double pressure_pa, double loaded_area_m2);
void drive_dynamic_bar_from_pressure_trace(cfd::fem::DynamicBar1D& bar,
                                           std::span<const double> pressure_pa,
                                           double loaded_area_m2);

void apply_acoustic_heating_to_pennes(cfd::multiphysics::PennesBioheat2D& solver,
                                      std::span<const double> intensity_w_m2,
                                      std::span<const double> absorption_np_per_m);

// Plane-wave absorption transfers momentum at approximately 2*alpha*I/c.
// Returned values are accelerations and can be passed directly to the FVM
// body-acceleration source.
[[nodiscard]] std::vector<cfd::fvm::Vec3> acoustic_streaming_acceleration(
    std::span<const double> intensity_w_m2,
    std::span<const double> absorption_np_per_m,
    std::span<const double> density_kg_m3,
    std::span<const double> sound_speed_m_s,
    cfd::fvm::Vec3 propagation_direction);

[[nodiscard]] double sensor_trace_l2_misfit(
    std::span<const std::vector<double>> predicted,
    std::span<const std::vector<double>> measured,
    bool normalize = true);

struct AcousticInverseResult {
    std::vector<double> parameters;
    double misfit{};
    std::size_t iterations{};
    bool converged{};
};

using AcousticTraceSimulator = std::function<std::vector<std::vector<double>>(std::span<const double>)>;

[[nodiscard]] AcousticInverseResult fit_acoustic_parameters(
    const AcousticTraceSimulator& simulator,
    std::span<const std::vector<double>> measured,
    std::vector<double> initial,
    std::vector<double> step,
    std::size_t max_iterations = 100,
    double tolerance = 1.0e-8);

} // namespace cfd::acoustics
