#include "cfd/solvers/fvm/solid_heat.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {
ScalarTransportConfig scalar_config(const SolidHeatConfig& c) {
    if (!(c.dt > 0.0) || !(c.density > 0.0) || !(c.heat_capacity > 0.0) || !(c.conductivity >= 0.0) ||
        !std::isfinite(c.dt) || !std::isfinite(c.density) || !std::isfinite(c.heat_capacity) || !std::isfinite(c.conductivity)) {
        throw std::invalid_argument("invalid solid heat controls");
    }
    ScalarTransportConfig out;
    out.dt = c.dt;
    out.diffusivity = c.conductivity / (c.density * c.heat_capacity);
    out.linear_iterations = c.linear_iterations;
    out.gmres_restart = c.gmres_restart;
    out.linear_tolerance = c.linear_tolerance;
    out.temporal_scheme = c.temporal_scheme;
    out.crank_nicolson_off_centering = c.crank_nicolson_off_centering;
    return out;
}
}

SolidHeatConduction::SolidHeatConduction(PolyMesh mesh, SolidHeatConfig config)
    : config_(config), solver_(std::move(mesh), scalar_config(config)) {}

void SolidHeatConduction::initialize(double temperature) { solver_.initialize(temperature); }
void SolidHeatConduction::initialize(const std::function<double(Vec3)>& temperature) { solver_.initialize(temperature); }

void SolidHeatConduction::set_boundary(std::string_view patch, SolidThermalBoundaryType type, double temperature) {
    solver_.set_boundary(patch,
        type == SolidThermalBoundaryType::fixed_temperature ? ScalarBoundaryType::fixedValue : ScalarBoundaryType::zeroGradient,
        temperature);
}

void SolidHeatConduction::set_volumetric_heat_source(double q) {
    if (!std::isfinite(q)) throw std::invalid_argument("solid heat source must be finite");
    solver_.set_source(q / (config_.density * config_.heat_capacity));
}

void SolidHeatConduction::set_volumetric_heat_source(const std::function<double(Vec3)>& q) {
    if (!q) throw std::invalid_argument("solid heat source callback is empty");
    const double rho_cp = config_.density * config_.heat_capacity;
    solver_.set_source([q, rho_cp](Vec3 x) { const double value=q(x); if(!std::isfinite(value)) throw std::runtime_error("non-finite solid heat source"); return value/rho_cp; });
}

cfd::core::IterativeSolverResult SolidHeatConduction::step() { return solver_.step(); }
cfd::core::IterativeSolverResult SolidHeatConduction::run(std::size_t steps) { return solver_.run(steps); }

double SolidHeatConduction::total_sensible_energy() const {
    double total = 0.0;
    for (std::size_t c=0;c<solver_.mesh().cell_count();++c) {
        total += config_.density * config_.heat_capacity * solver_.values()[c] * solver_.mesh().cells()[c].volume;
    }
    return total;
}

} // namespace cfd::fvm
