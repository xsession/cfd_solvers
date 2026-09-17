#pragma once
namespace cfd::fem {
struct PhaseChangeMaterial {double solidus{},liquidus{},heat_capacity_solid{},heat_capacity_liquid{},latent_heat{};[[nodiscard]] double liquid_fraction(double temperature)const;[[nodiscard]] double specific_enthalpy(double temperature)const;[[nodiscard]] double effective_heat_capacity(double temperature)const;[[nodiscard]] double temperature_from_enthalpy(double enthalpy)const;};
[[nodiscard]] double porous_effective_conductivity(double porosity,double solid_conductivity,double fluid_conductivity);
}
