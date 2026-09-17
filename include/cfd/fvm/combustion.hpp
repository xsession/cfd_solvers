#pragma once
#include <span>
#include <vector>
namespace cfd::fvm {
struct IdealGasSpeciesThermo { double molar_mass{}; double cp_mass{}; double formation_enthalpy_mass{}; };
struct CombustionThermoState { double temperature{},pressure{},density{},mixture_molar_mass{},gas_constant{},cp{},cv{},gamma{},specific_enthalpy{}; };
[[nodiscard]] CombustionThermoState ideal_gas_mixture_state(double temperature,double pressure,
    std::span<const double> mass_fractions,std::span<const IdealGasSpeciesThermo> species);
[[nodiscard]] double temperature_from_mixture_enthalpy(double specific_enthalpy,
    std::span<const double> mass_fractions,std::span<const IdealGasSpeciesThermo> species);
} // namespace cfd::fvm
