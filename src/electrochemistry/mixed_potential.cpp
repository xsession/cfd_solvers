#include "cfd/solvers/electrochemistry/mixed_potential.hpp"

#include "cfd/electrochemistry/electrochemistry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::electrochemistry {
namespace {

double total_current_at(const std::vector<ElectrodeReaction>& reactions,
                        double potential,
                        double temperature,
                        double target_current) {
    double total = -target_current;
    for (const auto& reaction : reactions) {
        const double eta = potential - reaction.equilibrium_potential;
        total += reaction.area * butler_volmer_current_density(
            reaction.exchange_current_density, eta, temperature, reaction.electrons,
            reaction.alpha_anodic, reaction.alpha_cathodic);
    }
    return total;
}

} // namespace

MixedPotentialResult solve_mixed_potential(const std::vector<ElectrodeReaction>& reactions,
                                           double temperature,
                                           double target_current,
                                           std::size_t max_iterations,
                                           double relative_current_tolerance) {
    if (reactions.empty()) throw std::invalid_argument("mixed-potential solve requires reactions");
    if (!(temperature > 0.0) || max_iterations == 0U || !(relative_current_tolerance > 0.0)
        || !std::isfinite(target_current)) {
        throw std::invalid_argument("invalid mixed-potential controls");
    }

    double min_eq = std::numeric_limits<double>::infinity();
    double max_eq = -std::numeric_limits<double>::infinity();
    double current_scale = std::abs(target_current);
    for (const auto& reaction : reactions) {
        if (reaction.name.empty() || !(reaction.exchange_current_density >= 0.0)
            || !(reaction.area > 0.0) || !(reaction.electrons > 0.0)
            || !(reaction.alpha_anodic > 0.0) || !(reaction.alpha_cathodic > 0.0)
            || !std::isfinite(reaction.equilibrium_potential)) {
            throw std::invalid_argument("invalid electrode reaction");
        }
        min_eq = std::min(min_eq, reaction.equilibrium_potential);
        max_eq = std::max(max_eq, reaction.equilibrium_potential);
        current_scale += reaction.area * reaction.exchange_current_density;
    }
    current_scale = std::max(current_scale, 1.0e-30);
    const double tolerance = relative_current_tolerance * current_scale;

    double margin = 0.1;
    double lo = min_eq - margin;
    double hi = max_eq + margin;
    double flo = total_current_at(reactions, lo, temperature, target_current);
    double fhi = total_current_at(reactions, hi, temperature, target_current);
    for (std::size_t expansion = 0; flo * fhi > 0.0 && expansion < 12U; ++expansion) {
        margin *= 2.0;
        lo = min_eq - margin;
        hi = max_eq + margin;
        flo = total_current_at(reactions, lo, temperature, target_current);
        fhi = total_current_at(reactions, hi, temperature, target_current);
    }
    if (!(flo <= 0.0 && fhi >= 0.0)) {
        throw std::runtime_error("mixed-potential current balance could not be bracketed");
    }

    MixedPotentialResult result;
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        const double fm = total_current_at(reactions, mid, temperature, target_current);
        result.iterations = iteration + 1U;
        result.potential = mid;
        result.total_current = fm + target_current;
        if (std::abs(fm) <= tolerance || std::abs(hi - lo) <= 1.0e-14) {
            result.converged = true;
            break;
        }
        if (fm > 0.0) hi = mid;
        else lo = mid;
    }

    result.reaction_current_density.resize(reactions.size());
    result.reaction_current.resize(reactions.size());
    result.total_current = 0.0;
    for (std::size_t i = 0; i < reactions.size(); ++i) {
        const auto& reaction = reactions[i];
        const double j = butler_volmer_current_density(
            reaction.exchange_current_density,
            result.potential - reaction.equilibrium_potential,
            temperature,
            reaction.electrons,
            reaction.alpha_anodic,
            reaction.alpha_cathodic);
        result.reaction_current_density[i] = j;
        result.reaction_current[i] = reaction.area * j;
        result.total_current += result.reaction_current[i];
    }
    return result;
}

} // namespace cfd::electrochemistry
