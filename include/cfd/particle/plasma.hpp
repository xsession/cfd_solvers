#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace cfd::particle {

struct PlasmaSpecies {
    std::string name;
    double charge_c{};
    double mass_kg{1.0};
};

struct PlasmaState {
    std::vector<PlasmaSpecies> species;
    std::vector<double> number_density_m3;
    double electron_temperature_ev{};
};

struct PlasmaReactionTerm {
    std::size_t species{};
    double coefficient{1.0};
};

struct PlasmaReaction {
    std::vector<PlasmaReactionTerm> reactants;
    std::vector<PlasmaReactionTerm> products;
    double rate_coefficient{};
};

struct PlasmaChemistryDiagnostics {
    std::vector<double> reaction_rates_m3_s;
    double charge_density_before_c_m3{};
    double charge_density_after_c_m3{};
    double total_number_density_before_m3{};
    double total_number_density_after_m3{};
    double max_relative_density_change{};
};

[[nodiscard]] double plasma_charge_density_c_m3(const PlasmaState& state);
[[nodiscard]] PlasmaChemistryDiagnostics advance_plasma_chemistry(
    PlasmaState& state,std::span<const PlasmaReaction> reactions,double dt_s,std::size_t substeps = 1U);

struct PaschenGas {
    double townsend_a_per_m_pa{};
    double townsend_b_v_per_m_pa{};
    double secondary_emission_yield{0.01};
};

struct GasBreakdownEstimate {
    double breakdown_voltage_v{};
    double breakdown_field_v_m{};
    bool applied_field_exceeds_threshold{};
};

[[nodiscard]] double paschen_breakdown_voltage_v(const PaschenGas& gas,double pressure_pa,double gap_m);
[[nodiscard]] GasBreakdownEstimate estimate_gas_breakdown_threshold(
    const PaschenGas& gas,double pressure_pa,double gap_m,double applied_field_v_m);

struct MultipactorEstimate {
    std::size_t order{};
    double transit_time_s{};
    double resonant_field_v_m{};
    double impact_energy_ev{};
    bool secondary_yield_sustains{};
};

[[nodiscard]] MultipactorEstimate estimate_parallel_plate_multipactor(
    double frequency_hz,double gap_m,std::size_t order,double particle_mass_kg,double particle_charge_c,
    double secondary_yield_at_impact);

} // namespace cfd::particle
