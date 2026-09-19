#include "cfd/solvers/fem/transport_stabilization.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {

SupgParameters supg_parameters(double length, double speed, double diffusivity, double timestep) {
    if (!(length > 0.0) || !(diffusivity >= 0.0) || !(timestep > 0.0) || !std::isfinite(length) ||
        !std::isfinite(speed) || !std::isfinite(diffusivity) || !std::isfinite(timestep)) {
        throw std::invalid_argument("invalid SUPG controls");
    }
    const double advective = 2.0 * std::abs(speed) / length;
    const double diffusive = 4.0 * diffusivity / (length * length);
    const double temporal = 2.0 / timestep;
    const double inverse_tau = std::sqrt(temporal * temporal + advective * advective + diffusive * diffusive);
    const double peclet = std::abs(speed) * length / std::max(2.0 * diffusivity, 1.0e-30);
    return {1.0 / inverse_tau, peclet};
}

double dg_upwind_flux(double left, double right, double outward_speed, double penalty) {
    if (!std::isfinite(left) || !std::isfinite(right) || !std::isfinite(outward_speed) || !(penalty >= 0.0) ||
        !std::isfinite(penalty)) {
        throw std::invalid_argument("invalid DG flux controls");
    }
    const double upwind = outward_speed >= 0.0 ? left : right;
    return outward_speed * upwind + penalty * (left - right);
}

} // namespace cfd::fem
