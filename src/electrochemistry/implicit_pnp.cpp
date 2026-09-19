#include "cfd/electrochemistry/implicit_pnp.hpp"

#include "cfd/chemistry/kinetics.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::electrochemistry {

ImplicitPnpElectrodeState solve_implicit_pnp_electrode_cell(double old_c, double old_phi,
                                                            const ImplicitPnpElectrodeConfig& config,
                                                            std::size_t max_iterations, double tolerance) {
    if (!(old_c >= 0.0) || !std::isfinite(old_phi) || !(config.dt_s > 0.0) || !(config.diffusivity_m2_s >= 0.0) ||
        !(config.control_volume_m3 > 0.0) || !(config.boundary_concentration_mol_m3 >= 0.0) ||
        !(config.exchange_current_a >= 0.0) || !(config.double_layer_capacitance_f > 0.0) ||
        !(config.temperature_k > 0.0) || !(config.electrons > 0.0) || max_iterations == 0U || !(tolerance > 0.0)) {
        throw std::invalid_argument("invalid implicit PNP/electrode controls");
    }
    const double transport = config.diffusivity_m2_s / config.control_volume_m3;
    double c = old_c, phi = old_phi;
    ImplicitPnpElectrodeState result{c, phi, 0U, false};
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        const double argument = config.electrons * faraday_constant * (phi - config.equilibrium_potential_v) /
                                (2.0 * cfd::chemistry::gas_constant * config.temperature_k);
        const double current = config.exchange_current_a * (std::exp(argument) - std::exp(-argument));
        const double didphi = config.exchange_current_a * config.electrons * faraday_constant /
                              (2.0 * cfd::chemistry::gas_constant * config.temperature_k) *
                              (std::exp(argument) + std::exp(-argument));
        const double rc = c - old_c -
                          config.dt_s * (transport * (config.boundary_concentration_mol_m3 - c) -
                                         current / (config.electrons * faraday_constant * config.control_volume_m3));
        const double rp = config.double_layer_capacitance_f * (phi - old_phi) + config.dt_s * current;
        result.iterations = iteration + 1U;
        if (std::max(std::abs(rc), std::abs(rp)) <= tolerance) {
            result.concentration_mol_m3 = std::max(0.0, c);
            result.potential_v = phi;
            result.converged = true;
            return result;
        }
        const double j11 = 1.0 + config.dt_s * transport;
        const double j12 = config.dt_s * didphi / (config.electrons * faraday_constant * config.control_volume_m3);
        const double j21 = 0.0;
        const double j22 = config.double_layer_capacitance_f + config.dt_s * didphi;
        const double determinant = j11 * j22 - j12 * j21;
        if (std::abs(determinant) < 1.0e-30)
            throw std::runtime_error("implicit PNP Jacobian is singular");
        const double dc = (-rc * j22 + j12 * rp) / determinant;
        const double dphi = (-j11 * rp + j21 * rc) / determinant;
        c = std::max(0.0, c + dc);
        phi += dphi;
    }
    result.concentration_mol_m3 = c;
    result.potential_v = phi;
    return result;
}

} // namespace cfd::electrochemistry
