#include "cfd/solvers/fdtd/multiconductor_tline.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fdtd {

void advance_multiconductor_tline(const MulticonductorTlineConfig& config, MulticonductorTlineState& state,
                                  std::span<const double> source) {
    if (config.conductors == 0U || config.inductance_h.size() != config.conductors ||
        config.capacitance_f.size() != config.conductors || state.voltage.size() != config.conductors ||
        state.current.size() != config.conductors || source.size() != config.conductors || !(config.dt_s > 0.0))
        throw std::invalid_argument("invalid multiconductor line controls");
    for (std::size_t i = 0; i < config.conductors; ++i) {
        if (!(config.inductance_h[i] > 0.0) || !(config.capacitance_f[i] > 0.0) || !std::isfinite(source[i]))
            throw std::invalid_argument("invalid multiconductor line state");
        const double current_old = state.current[i];
        state.current[i] += config.dt_s * (state.voltage[i] - source[i]) / config.inductance_h[i];
        state.voltage[i] -= config.dt_s * current_old / config.capacitance_f[i];
    }
}

} // namespace cfd::fdtd
