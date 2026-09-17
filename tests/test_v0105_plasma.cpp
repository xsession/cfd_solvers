#include "cfd/particle/plasma.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 

void plasma_chemistry_preserves_charge(){
    using namespace cfd::particle;
    constexpr double qe=1.602176634e-19;
    PlasmaState state;state.species={{"e",-qe,9.1093837015e-31},{"Ar",0.0,6.6335209e-26},{"Ar+",qe,6.6335209e-26}};
    state.number_density_m3={1.0e15,1.0e20,1.0e15};state.electron_temperature_ev=8.0;
    PlasmaReaction ionization;ionization.reactants={{0U,1.0},{1U,1.0}};ionization.products={{0U,2.0},{2U,1.0}};ionization.rate_coefficient=1.0e-20;
    const double charge0=plasma_charge_density_c_m3(state);
    const auto diagnostics=advance_plasma_chemistry(state,std::span<const PlasmaReaction>(&ionization,1U),1.0e-7,10U);
    check(state.number_density_m3[0]>1.0e15,"ionization increases electron density");
    check(state.number_density_m3[1]<1.0e20,"ionization consumes neutral density");
    check(state.number_density_m3[2]>1.0e15,"ionization increases ion density");
    check(std::abs(diagnostics.charge_density_after_c_m3-charge0)<1.0e-12,"ionization reaction preserves charge density");
    check(diagnostics.reaction_rates_m3_s.size()==1U&&diagnostics.reaction_rates_m3_s.front()>0.0,"reaction-rate diagnostic");
    check(diagnostics.max_relative_density_change>0.0,"density-change diagnostic");
}

void gas_breakdown_and_multipactor(){
    using namespace cfd::particle;
    PaschenGas air;air.townsend_a_per_m_pa=112.5;air.townsend_b_v_per_m_pa=2737.0;air.secondary_emission_yield=0.01;
    const auto below=estimate_gas_breakdown_threshold(air,101325.0,1.0e-3,1.0e6);
    const auto above=estimate_gas_breakdown_threshold(air,101325.0,1.0e-3,below.breakdown_field_v_m*1.1);
    check(below.breakdown_voltage_v>0.0&&below.breakdown_field_v_m>0.0,"positive Paschen threshold");
    check(!below.applied_field_exceeds_threshold&&above.applied_field_exceeds_threshold,"gas breakdown threshold comparison");
    const auto first=estimate_parallel_plate_multipactor(2.0e9,1.0e-3,1U,9.1093837015e-31,-1.602176634e-19,1.4);
    const auto third=estimate_parallel_plate_multipactor(2.0e9,1.0e-3,3U,9.1093837015e-31,-1.602176634e-19,0.8);
    check(first.transit_time_s<third.transit_time_s,"multipactor order transit time");
    check(first.resonant_field_v_m>third.resonant_field_v_m,"higher-order multipactor requires lower field");
    check(first.secondary_yield_sustains&&!third.secondary_yield_sustains,"multipactor secondary yield decision");
}

void rejects_invalid_inputs(){
    using namespace cfd::particle;
    bool caught=false;try{(void)paschen_breakdown_voltage_v({0.0,2737.0,0.01},101325.0,1.0e-3);}catch(const std::exception&){caught=true;}check(caught,"invalid Paschen gas rejected");
    caught=false;try{(void)estimate_parallel_plate_multipactor(0.0,1.0e-3,1U,1.0,-1.0,1.0);}catch(const std::exception&){caught=true;}check(caught,"invalid multipactor controls rejected");
}
}

int main(){
    try{
        plasma_chemistry_preserves_charge();
        gas_breakdown_and_multipactor();
        rejects_invalid_inputs();
    }catch(const std::exception& error){
        std::cerr<<"v0.10.5 plasma test failed: "<<error.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.10.5 plasma tests passed\n";
    return 0;
}
