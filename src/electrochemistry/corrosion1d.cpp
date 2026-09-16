#include "cfd/solvers/electrochemistry/corrosion1d.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::electrochemistry {
namespace {

double residual(const CorrosionCell1DConfig& config, double current) {
    const double surface_phi = config.bulk_electrolyte_potential
                             + current * config.electrolyte_length / config.electrolyte_conductivity;
    const double eta = config.metal_potential - surface_phi - config.equilibrium_potential;
    return current - butler_volmer_current_density(config.exchange_current_density,
                                                   eta,
                                                   config.temperature,
                                                   config.electrons,
                                                   config.anodic_transfer,
                                                   config.cathodic_transfer);
}

} // namespace

CorrosionCell1DResult solve_corrosion_cell_1d(const CorrosionCell1DConfig& config) {
    if (!(config.electrolyte_length > 0.0) || !(config.electrolyte_conductivity > 0.0)
        || !(config.temperature > 0.0) || config.exchange_current_density < 0.0
        || !(config.electrons > 0.0) || config.max_iterations == 0U
        || !(config.current_tolerance > 0.0)) {
        throw std::invalid_argument("invalid 1-D corrosion cell configuration");
    }

    const double applied = config.metal_potential
                         - config.bulk_electrolyte_potential
                         - config.equilibrium_potential;
    if (std::abs(applied) <= 1.0e-18 || config.exchange_current_density == 0.0) {
        return {0.0, config.bulk_electrolyte_potential, applied, 0.0, 0.0, 0U, true};
    }

    double scale = std::max(config.exchange_current_density,
                            std::abs(applied) * config.electrolyte_conductivity / config.electrolyte_length);
    scale = std::max(scale, 1.0e-12);
    double lo = -scale;
    double hi = scale;
    double flo = residual(config, lo);
    double fhi = residual(config, hi);
    for (std::size_t expand = 0U; expand < 80U && !(flo <= 0.0 && fhi >= 0.0); ++expand) {
        lo *= 2.0;
        hi *= 2.0;
        flo = residual(config, lo);
        fhi = residual(config, hi);
        if (!std::isfinite(flo) || !std::isfinite(fhi)) break;
    }
    if (!(std::isfinite(flo) && std::isfinite(fhi) && flo <= 0.0 && fhi >= 0.0)) {
        throw std::runtime_error("failed to bracket corrosion-cell current density");
    }

    double current = 0.0;
    std::size_t iteration = 0U;
    bool converged = false;
    for (; iteration < config.max_iterations; ++iteration) {
        current = 0.5 * (lo + hi);
        const double f = residual(config, current);
        if (std::abs(f) <= config.current_tolerance * std::max(1.0, std::abs(current))
            || 0.5 * std::abs(hi - lo) <= config.current_tolerance * std::max(1.0, std::abs(current))) {
            converged = true;
            ++iteration;
            break;
        }
        if (f > 0.0) hi = current;
        else lo = current;
    }

    const double surface_phi = config.bulk_electrolyte_potential
                             + current * config.electrolyte_length / config.electrolyte_conductivity;
    const double eta = config.metal_potential - surface_phi - config.equilibrium_potential;
    const double molar_flux = faradaic_molar_flux(current, config.electrons);
    const double penetration = corrosion_penetration_rate(current,
                                                          config.metal_molar_mass,
                                                          config.metal_density,
                                                          config.electrons);
    return {current, surface_phi, eta, molar_flux, penetration, iteration, converged};
}

} // namespace cfd::electrochemistry
