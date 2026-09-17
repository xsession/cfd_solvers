#include "cfd/chemistry/kinetics.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace cfd::chemistry {

double EquilibriumConstant::value(double temperature) const {
    if (!(temperature>0.0) || !std::isfinite(temperature) || !(reference_temperature>0.0)
        || !std::isfinite(reference_temperature) || !(reference_value>0.0)
        || !std::isfinite(reference_value) || !std::isfinite(reaction_enthalpy))
        throw std::invalid_argument("invalid equilibrium constant parameters");
    const double result=reference_value*std::exp(-reaction_enthalpy/gas_constant
        *(1.0/temperature-1.0/reference_temperature));
    if (!(result>0.0) || !std::isfinite(result)) throw std::overflow_error("equilibrium constant out of range");
    return result;
}

double ArrheniusRate::rate_constant(double temperature) const {
    if (!(temperature > 0.0)||!std::isfinite(temperature)) throw std::invalid_argument("temperature must be finite and positive");
    if (!(pre_exponential >= 0.0)||!std::isfinite(pre_exponential)
        ||!std::isfinite(temperature_exponent)||!std::isfinite(activation_energy))
        throw std::invalid_argument("invalid Arrhenius parameters");
    return pre_exponential * std::pow(temperature, temperature_exponent)
         * std::exp(-activation_energy / (gas_constant * temperature));
}

ReactionNetwork::ReactionNetwork(std::vector<Species> species) : species_(std::move(species)) {}

std::size_t ReactionNetwork::add_species(Species species) {
    species_.push_back(std::move(species));
    return species_.size() - 1U;
}

void ReactionNetwork::validate_reaction(const ElementaryReaction& reaction) const {
    if (reaction.reactants.empty() && reaction.products.empty()) {
        throw std::invalid_argument("reaction must contain at least one stoichiometric term");
    }
    const auto check = [&](const StoichiometricTerm& term) {
        if (term.species >= species_.size()) throw std::out_of_range("reaction species index out of range");
        if (!(term.coefficient > 0.0)||!std::isfinite(term.coefficient)) throw std::invalid_argument("stoichiometric coefficient must be finite and positive");
    };
    for (const auto& term : reaction.reactants) check(term);
    for (const auto& term : reaction.products) check(term);
    for (const auto& order : reaction.forward_orders) {
        if (order.species >= species_.size()) throw std::out_of_range("reaction-order species index out of range");
        if (!std::isfinite(order.order)) throw std::invalid_argument("reaction order must be finite");
    }
    for (const auto& efficiency : reaction.third_body_efficiencies) {
        if (efficiency.species >= species_.size()) throw std::out_of_range("third-body species index out of range");
        if (!(efficiency.efficiency >= 0.0) || !std::isfinite(efficiency.efficiency))
            throw std::invalid_argument("third-body efficiency must be finite and non-negative");
    }
}

void ReactionNetwork::add_reaction(ElementaryReaction reaction) {
    validate_reaction(reaction);
    reactions_.push_back(std::move(reaction));
}

std::vector<double> ReactionNetwork::reaction_rates(std::span<const double> concentrations,
                                                    double temperature) const {
    std::vector<double> rates(reactions_.size(), 0.0);
    reaction_rates(concentrations,temperature,rates);
    return rates;
}

void ReactionNetwork::reaction_rates(std::span<const double> concentrations,double temperature,
                                     std::span<double> rates) const {
    if (concentrations.size()!=species_.size()||rates.size()!=reactions_.size())
        throw std::invalid_argument("reaction buffer size mismatch");
    if (!(temperature>0.0)||!std::isfinite(temperature)) throw std::invalid_argument("invalid reaction temperature");
    for(double c:concentrations) if(!(c>=0.0)||!std::isfinite(c))
        throw std::invalid_argument("invalid reaction concentration");
    for (std::size_t r = 0; r < reactions_.size(); ++r) {
        const auto& reaction = reactions_[r];
        const double k_inf = reaction.forward.rate_constant(temperature);
        const bool uses_third_body = reaction.third_body || reaction.low_pressure_limit.has_value();
        double collider = 1.0;
        if (uses_third_body) {
            collider = 0.0;
            for (double concentration : concentrations) collider += concentration;
            for (const auto& efficiency : reaction.third_body_efficiencies)
                collider += (efficiency.efficiency - 1.0) * concentrations[efficiency.species];
            if (!(collider >= 0.0) || !std::isfinite(collider))
                throw std::overflow_error("third-body concentration out of range");
        }
        double k_effective = k_inf;
        if (reaction.low_pressure_limit) {
            const double k0 = reaction.low_pressure_limit->rate_constant(temperature);
            if (k_inf == 0.0) k_effective = 0.0;
            else {
                const double reduced_pressure = k0 * collider / k_inf;
                k_effective = k_inf * reduced_pressure / (1.0 + reduced_pressure);
            }
        } else if (reaction.third_body) {
            k_effective *= collider;
        }

        double rate = k_effective;
        if (!reaction.forward_orders.empty()) {
            for (const auto& order : reaction.forward_orders)
                rate *= std::pow(concentrations[order.species], order.order);
        } else {
            for (const auto& term : reaction.reactants)
                rate *= std::pow(concentrations[term.species], term.coefficient);
        }
        if(reaction.equilibrium){
            double reverse=k_effective/reaction.equilibrium->value(temperature);
            for(const auto& term:reaction.products) reverse*=std::pow(concentrations[term.species],term.coefficient);
            rate-=reverse;
        }
        if(!std::isfinite(rate)) throw std::overflow_error("non-finite reaction rate");
        rates[r] = rate;
    }
}

std::vector<double> ReactionNetwork::source_terms(std::span<const double> concentrations,
                                                  double temperature) const {
    std::vector<double> source(species_.size(), 0.0);
    std::vector<double> rates(reactions_.size());
    source_terms(concentrations,temperature,source,rates);
    return source;
}

void ReactionNetwork::source_terms(std::span<const double> concentrations,double temperature,
                                  std::span<double> source,std::span<double> rates) const {
    if(source.size()!=species_.size()) throw std::invalid_argument("reaction source size mismatch");
    reaction_rates(concentrations,temperature,rates);
    std::fill(source.begin(),source.end(),0.0);
    for (std::size_t r = 0; r < reactions_.size(); ++r) {
        for (const auto& term : reactions_[r].reactants) source[term.species] -= term.coefficient * rates[r];
        for (const auto& term : reactions_[r].products) source[term.species] += term.coefficient * rates[r];
    }
}

} // namespace cfd::chemistry
