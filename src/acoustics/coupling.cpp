#include "cfd/acoustics/coupling.hpp"
#include "cfd/workflow/runner.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::acoustics {

std::vector<StructuralNodalForce2D> pressure_to_structural_boundary_load(
    const cfd::fem::Mesh2D& mesh,
    std::span<const double> edge_pressure_pa,
    int patch) {
    mesh.validate();
    if (edge_pressure_pa.size() != mesh.boundary_edges.size())
        throw std::invalid_argument("boundary pressure size mismatch");
    std::vector<StructuralNodalForce2D> out(mesh.node_count());
    for (std::size_t e = 0; e < mesh.boundary_edges.size(); ++e) {
        const auto& edge = mesh.boundary_edges[e];
        if (patch >= 0 && edge.patch != patch) continue;
        const double p = edge_pressure_pa[e];
        if (!std::isfinite(p)) throw std::invalid_argument("non-finite acoustic boundary pressure");
        const auto& a = mesh.nodes[edge.node[0]];
        const auto& b = mesh.nodes[edge.node[1]];
        const double dx = b.x - a.x, dy = b.y - a.y;
        const double length = std::hypot(dx, dy);
        if (!(length > 0.0)) throw std::runtime_error("degenerate FEM boundary edge");
        // Boundary-edge orientation is not guaranteed globally, so use the
        // right-hand normal of the stored edge. Conservation tests should use
        // a consistently oriented selected patch.
        const double nx = dy / length, ny = -dx / length;
        const double fx = -p * nx * length * 0.5;
        const double fy = -p * ny * length * 0.5;
        out[edge.node[0]].x += fx; out[edge.node[0]].y += fy;
        out[edge.node[1]].x += fx; out[edge.node[1]].y += fy;
    }
    return out;
}

double acoustic_pressure_to_bar_tip_force(double pressure_pa, double loaded_area_m2) {
    if (!std::isfinite(pressure_pa) || !(loaded_area_m2 >= 0.0) || !std::isfinite(loaded_area_m2))
        throw std::invalid_argument("invalid acoustic bar loading");
    return pressure_pa * loaded_area_m2;
}

void drive_dynamic_bar_from_pressure_trace(cfd::fem::DynamicBar1D& bar,
                                           std::span<const double> pressure_pa,
                                           double loaded_area_m2) {
    if (!(loaded_area_m2 >= 0.0) || !std::isfinite(loaded_area_m2)) throw std::invalid_argument("invalid bar loading area");
    for (double p:pressure_pa) bar.step(acoustic_pressure_to_bar_tip_force(p,loaded_area_m2));
}

void apply_acoustic_heating_to_pennes(cfd::multiphysics::PennesBioheat2D& solver,
                                      std::span<const double> intensity,
                                      std::span<const double> absorption) {
    if (intensity.size()!=absorption.size() || intensity.size()!=solver.temperature_k().size())
        throw std::invalid_argument("acoustic/Pennes grid size mismatch");
    std::vector<double> heating(intensity.size());
    for (std::size_t i=0;i<intensity.size();++i) {
        if (!(intensity[i]>=0.0)||!(absorption[i]>=0.0)||!std::isfinite(intensity[i])||!std::isfinite(absorption[i]))
            throw std::invalid_argument("invalid acoustic heating input");
        heating[i]=2.0*absorption[i]*intensity[i];
    }
    solver.set_volumetric_heating(heating);
}

std::vector<cfd::fvm::Vec3> acoustic_streaming_acceleration(
    std::span<const double> intensity,
    std::span<const double> absorption,
    std::span<const double> density,
    std::span<const double> sound_speed,
    cfd::fvm::Vec3 direction) {
    const std::size_t n = intensity.size();
    if (absorption.size()!=n || density.size()!=n || sound_speed.size()!=n)
        throw std::invalid_argument("acoustic streaming field size mismatch");
    const double mag = cfd::fvm::magnitude(direction);
    if (!(mag > 0.0) || !std::isfinite(mag)) throw std::invalid_argument("invalid acoustic propagation direction");
    direction = direction / mag;
    std::vector<cfd::fvm::Vec3> out(n);
    for (std::size_t i=0;i<n;++i) {
        if (!(intensity[i]>=0.0) || !(absorption[i]>=0.0) || !(density[i]>0.0) || !(sound_speed[i]>0.0) ||
            !std::isfinite(intensity[i]) || !std::isfinite(absorption[i]) || !std::isfinite(density[i]) || !std::isfinite(sound_speed[i]))
            throw std::invalid_argument("invalid acoustic streaming material/field value");
        const double acceleration = 2.0 * absorption[i] * intensity[i] / (density[i] * sound_speed[i]);
        out[i] = direction * acceleration;
    }
    return out;
}

double sensor_trace_l2_misfit(std::span<const std::vector<double>> predicted,
                              std::span<const std::vector<double>> measured,
                              bool normalize) {
    if (predicted.size()!=measured.size() || predicted.empty()) throw std::invalid_argument("sensor trace set mismatch");
    double error=0.0, reference=0.0;
    for (std::size_t s=0;s<predicted.size();++s) {
        if (predicted[s].size()!=measured[s].size() || predicted[s].empty()) throw std::invalid_argument("sensor trace length mismatch");
        for (std::size_t k=0;k<predicted[s].size();++k) {
            const double a=predicted[s][k], b=measured[s][k];
            if (!std::isfinite(a)||!std::isfinite(b)) throw std::invalid_argument("non-finite sensor trace");
            const double d=a-b; error += d*d; reference += b*b;
        }
    }
    if (normalize && reference>0.0) error/=reference;
    return std::sqrt(error);
}

AcousticInverseResult fit_acoustic_parameters(const AcousticTraceSimulator& simulator,
                                               std::span<const std::vector<double>> measured,
                                               std::vector<double> initial,
                                               std::vector<double> step,
                                               std::size_t max_iterations,
                                               double tolerance) {
    if (!simulator) throw std::invalid_argument("empty acoustic trace simulator");
    auto objective = [&](std::span<const double> x){ return sensor_trace_l2_misfit(simulator(x), measured, true); };
    auto result = cfd::workflow::inverse_coordinate_search(objective,std::move(initial),std::move(step),max_iterations,tolerance);
    return {std::move(result.parameters),result.objective,result.iterations,result.converged};
}

} // namespace cfd::acoustics
