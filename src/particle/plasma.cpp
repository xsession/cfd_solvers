#include "cfd/particle/plasma.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfd::particle {
namespace {
constexpr double elementary_charge_c = 1.602176634e-19;

void validate_state(const PlasmaState& state){
    if(state.species.empty()||state.species.size()!=state.number_density_m3.size())throw std::invalid_argument("invalid plasma state size");
    for(std::size_t i=0;i<state.species.size();++i){
        if(!(state.species[i].mass_kg>0.0)||!std::isfinite(state.species[i].charge_c))throw std::invalid_argument("invalid plasma species");
        if(!(state.number_density_m3[i]>=0.0)||!std::isfinite(state.number_density_m3[i]))throw std::invalid_argument("invalid plasma density");
    }
}

double total_density(const PlasmaState& state){double total=0.0;for(double n:state.number_density_m3)total+=n;return total;}

double reaction_rate(const PlasmaState& state,const PlasmaReaction& reaction){
    if(!(reaction.rate_coefficient>=0.0)||!std::isfinite(reaction.rate_coefficient))throw std::invalid_argument("invalid plasma reaction rate coefficient");
    double rate=reaction.rate_coefficient;
    for(const auto& term:reaction.reactants){
        if(term.species>=state.number_density_m3.size()||!(term.coefficient>0.0)||!std::isfinite(term.coefficient))throw std::invalid_argument("invalid plasma reactant term");
        rate*=std::pow(state.number_density_m3[term.species],term.coefficient);
    }
    for(const auto& term:reaction.products){
        if(term.species>=state.number_density_m3.size()||!(term.coefficient>0.0)||!std::isfinite(term.coefficient))throw std::invalid_argument("invalid plasma product term");
    }
    return rate;
}
}

double plasma_charge_density_c_m3(const PlasmaState& state){
    validate_state(state);double charge=0.0;for(std::size_t i=0;i<state.species.size();++i)charge+=state.species[i].charge_c*state.number_density_m3[i];return charge;
}

PlasmaChemistryDiagnostics advance_plasma_chemistry(PlasmaState& state,std::span<const PlasmaReaction> reactions,double dt_s,std::size_t substeps){
    validate_state(state);if(reactions.empty()||!(dt_s>=0.0)||!std::isfinite(dt_s)||substeps==0U)throw std::invalid_argument("invalid plasma chemistry controls");
    PlasmaChemistryDiagnostics diagnostics;diagnostics.charge_density_before_c_m3=plasma_charge_density_c_m3(state);diagnostics.total_number_density_before_m3=total_density(state);diagnostics.reaction_rates_m3_s.assign(reactions.size(),0.0);
    const auto initial_density=state.number_density_m3;const double dt_sub=dt_s/static_cast<double>(substeps);
    for(std::size_t step=0;step<substeps;++step){
        std::vector<double> delta(state.number_density_m3.size(),0.0);
        for(std::size_t r=0;r<reactions.size();++r){
            const auto& reaction=reactions[r];const double rate=reaction_rate(state,reaction);diagnostics.reaction_rates_m3_s[r]=rate;
            for(const auto& term:reaction.reactants)delta[term.species]-=term.coefficient*rate*dt_sub;
            for(const auto& term:reaction.products)delta[term.species]+=term.coefficient*rate*dt_sub;
        }
        for(std::size_t i=0;i<state.number_density_m3.size();++i){
            const double updated=state.number_density_m3[i]+delta[i];
            if(updated<-1.0e-9*std::max(1.0,state.number_density_m3[i]))throw std::runtime_error("plasma chemistry step produced negative density");
            state.number_density_m3[i]=std::max(0.0,updated);
        }
    }
    diagnostics.charge_density_after_c_m3=plasma_charge_density_c_m3(state);diagnostics.total_number_density_after_m3=total_density(state);
    for(std::size_t i=0;i<state.number_density_m3.size();++i){
        const double scale=std::max(1.0,initial_density[i]);
        diagnostics.max_relative_density_change=std::max(diagnostics.max_relative_density_change,std::abs(state.number_density_m3[i]-initial_density[i])/scale);
    }
    return diagnostics;
}

double paschen_breakdown_voltage_v(const PaschenGas& gas,double pressure_pa,double gap_m){
    if(!(gas.townsend_a_per_m_pa>0.0)||!(gas.townsend_b_v_per_m_pa>0.0)||!(gas.secondary_emission_yield>0.0)||!(pressure_pa>0.0)||!(gap_m>0.0))throw std::invalid_argument("invalid Paschen controls");
    const double pd=pressure_pa*gap_m;const double argument=gas.townsend_a_per_m_pa*pd;const double secondary=std::log(1.0+1.0/gas.secondary_emission_yield);
    if(!(argument>1.0)||!(secondary>1.0))throw std::invalid_argument("Paschen controls outside breakdown-law domain");
    const double denominator=std::log(argument)-std::log(secondary);
    if(!(denominator>0.0))throw std::invalid_argument("Paschen denominator is not positive");
    return gas.townsend_b_v_per_m_pa*pd/denominator;
}

GasBreakdownEstimate estimate_gas_breakdown_threshold(const PaschenGas& gas,double pressure_pa,double gap_m,double applied_field_v_m){
    if(!(applied_field_v_m>=0.0)||!std::isfinite(applied_field_v_m))throw std::invalid_argument("invalid applied gas-breakdown field");
    GasBreakdownEstimate estimate;estimate.breakdown_voltage_v=paschen_breakdown_voltage_v(gas,pressure_pa,gap_m);estimate.breakdown_field_v_m=estimate.breakdown_voltage_v/gap_m;estimate.applied_field_exceeds_threshold=applied_field_v_m>=estimate.breakdown_field_v_m;return estimate;
}

MultipactorEstimate estimate_parallel_plate_multipactor(double frequency_hz,double gap_m,std::size_t order,double particle_mass_kg,double particle_charge_c,double secondary_yield_at_impact){
    if(!(frequency_hz>0.0)||!(gap_m>0.0)||order==0U||!(particle_mass_kg>0.0)||!(std::abs(particle_charge_c)>0.0)||!(secondary_yield_at_impact>=0.0))throw std::invalid_argument("invalid multipactor estimate controls");
    MultipactorEstimate estimate;estimate.order=order;estimate.transit_time_s=(static_cast<double>(2U*order-1U))/(2.0*frequency_hz);
    const double q=std::abs(particle_charge_c);estimate.resonant_field_v_m=2.0*particle_mass_kg*gap_m/(q*estimate.transit_time_s*estimate.transit_time_s);
    const double impact_velocity=2.0*gap_m/estimate.transit_time_s;estimate.impact_energy_ev=0.5*particle_mass_kg*impact_velocity*impact_velocity/elementary_charge_c;
    estimate.secondary_yield_sustains=secondary_yield_at_impact>1.0;return estimate;
}

} // namespace cfd::particle
