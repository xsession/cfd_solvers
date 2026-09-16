#pragma once

#include <span>
#include <vector>

namespace cfd::chemistry {

// Ideal dilute acid/base families only. Concentrations and Ka use mol/L here;
// conversion to transport's mol/m^3 is explicit (multiply concentrations by 1000).
struct AcidFamily {
    double total_mol_per_litre{};
    int fully_protonated_charge{};
    std::vector<double> dissociation_constants;
};
struct AqueousEquilibriumResult {
    double ph{};
    double hydrogen_mol_per_litre{};
    double hydroxide_mol_per_litre{};
    double charge_residual_mol_per_litre{};
    std::vector<std::vector<double>> family_concentrations;
};

// strong_charge is the fixed signed sum of spectator-ion charge concentration.
// Kw and Ka are caller-supplied constants at the desired temperature.
[[nodiscard]] AqueousEquilibriumResult equilibrate_acids(std::span<const AcidFamily> families,
    double strong_charge_mol_per_litre=0.0,double water_ion_product=1.0e-14);

} // namespace cfd::chemistry
