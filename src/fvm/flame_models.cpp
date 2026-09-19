#include "cfd/solvers/fvm/flame_models.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {

void validate_state(FlameModelState state) {
    if (!(state.progress >= 0.0 && state.progress <= 1.0) ||
        !(state.mixture_fraction >= 0.0 && state.mixture_fraction <= 1.0) || !(state.temperature >= 0.0) ||
        !std::isfinite(state.progress) || !std::isfinite(state.mixture_fraction) || !std::isfinite(state.temperature)) {
        throw std::invalid_argument("invalid flame model state");
    }
}

} // namespace

PremixedFlameModel::PremixedFlameModel(PremixedFlameModelConfig config) : config_(config) {
    if (!(config_.reaction_rate >= 0.0) || !(config_.heat_release >= 0.0) || !(config_.ignition_temperature >= 0.0) ||
        !std::isfinite(config_.reaction_rate) || !std::isfinite(config_.heat_release) ||
        !std::isfinite(config_.ignition_temperature)) {
        throw std::invalid_argument("invalid premixed flame controls");
    }
}

FlameModelResult PremixedFlameModel::evaluate(FlameModelState state) const {
    validate_state(state);
    if (state.temperature < config_.ignition_temperature)
        return {};
    const double source = config_.reaction_rate * state.progress * (1.0 - state.progress);
    return {source, 0.0, config_.heat_release * source};
}

NonPremixedFlameModel::NonPremixedFlameModel(NonPremixedFlameModelConfig config) : config_(config) {
    if (!(config_.stoichiometric_mixture_fraction > 0.0 && config_.stoichiometric_mixture_fraction < 1.0) ||
        !(config_.mixing_width > 0.0) || !(config_.mixing_rate >= 0.0) || !(config_.heat_release >= 0.0) ||
        !(config_.ignition_temperature >= 0.0) || !std::isfinite(config_.stoichiometric_mixture_fraction) ||
        !std::isfinite(config_.mixing_width) || !std::isfinite(config_.mixing_rate) ||
        !std::isfinite(config_.heat_release) || !std::isfinite(config_.ignition_temperature)) {
        throw std::invalid_argument("invalid non-premixed flame controls");
    }
}

FlameModelResult NonPremixedFlameModel::evaluate(FlameModelState state) const {
    validate_state(state);
    if (state.temperature < config_.ignition_temperature)
        return {};
    const double normalized = (state.mixture_fraction - config_.stoichiometric_mixture_fraction) / config_.mixing_width;
    const double stoichiometric_window = std::exp(-normalized * normalized);
    const double progress_source = config_.mixing_rate * stoichiometric_window * (1.0 - state.progress);
    const double mixture_source =
        -config_.mixing_rate * (state.mixture_fraction - config_.stoichiometric_mixture_fraction);
    return {progress_source, mixture_source, config_.heat_release * progress_source};
}

} // namespace cfd::fvm
