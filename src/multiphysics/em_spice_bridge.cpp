#include "cfd/multiphysics/em_spice_bridge.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace cfd::multiphysics {

EmSpiceCouplingResult iterate_em_spice_port(std::complex<double> current, const FieldPortSolve& field_solve,
                                            const CircuitPortSolve& circuit_solve,
                                            const EmSpiceCouplingConfig& config) {
    if (!field_solve || !circuit_solve || config.max_iterations == 0U ||
        !(config.relaxation > 0.0 && config.relaxation <= 1.0) || !(config.tolerance > 0.0) ||
        !std::isfinite(config.tolerance)) {
        throw std::invalid_argument("invalid EM/SPICE coupling controls");
    }
    EmSpiceCouplingResult result;
    for (std::size_t iteration = 0; iteration < config.max_iterations; ++iteration) {
        result.port_voltage = field_solve(current);
        const auto proposed = circuit_solve(result.port_voltage);
        const auto updated = (1.0 - config.relaxation) * current + config.relaxation * proposed;
        result.iterations = iteration + 1U;
        if (std::abs(updated - current) <= config.tolerance * std::max(1.0, std::abs(updated))) {
            result.port_current = updated;
            result.converged = true;
            return result;
        }
        current = updated;
    }
    result.port_current = current;
    return result;
}

} // namespace cfd::multiphysics
