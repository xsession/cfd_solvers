#include "cfd/solvers/fvm/soot_radiation.hpp"

#include "cfd/fvm/radiation.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {

SootRadiationModel::SootRadiationModel(SootRadiationConfig config) : config_(config) {
    if (!(config_.soot_yield >= 0.0) || !(config_.oxidation_rate >= 0.0) || !(config_.absorption_per_mass >= 0.0) ||
        !(config_.emissivity_floor >= 0.0) || !std::isfinite(config_.soot_yield) ||
        !std::isfinite(config_.oxidation_rate) || !std::isfinite(config_.absorption_per_mass) ||
        !std::isfinite(config_.emissivity_floor)) {
        throw std::invalid_argument("invalid soot-radiation controls");
    }
}

SootRadiationResult SootRadiationModel::evaluate(SootRadiationState state, double environment_temperature) const {
    if (!(state.soot_mass_fraction >= 0.0 && state.soot_mass_fraction <= 1.0) || !(state.temperature >= 0.0) ||
        !(environment_temperature >= 0.0) || !std::isfinite(state.soot_mass_fraction) ||
        !std::isfinite(state.temperature) || !std::isfinite(state.heat_release) ||
        !std::isfinite(environment_temperature)) {
        throw std::invalid_argument("invalid soot-radiation state");
    }
    const double production = config_.soot_yield * std::max(state.heat_release, 0.0);
    const double oxidation = config_.oxidation_rate * state.soot_mass_fraction;
    const double source = production - oxidation;
    const double absorption = config_.absorption_per_mass * state.soot_mass_fraction;
    const double radiation = optically_thin_radiation_source(state.temperature, environment_temperature, absorption);
    return {source, absorption, radiation};
}

} // namespace cfd::fvm
