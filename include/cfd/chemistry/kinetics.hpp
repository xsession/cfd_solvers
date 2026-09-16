#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace cfd::chemistry {

inline constexpr double gas_constant = 8.31446261815324; // J/(mol K)

struct Species {
    std::string name;
    double molar_mass{}; // kg/mol
    int charge{};
};

struct ArrheniusRate {
    double pre_exponential{};
    double temperature_exponent{};
    double activation_energy{}; // J/mol

    [[nodiscard]] double rate_constant(double temperature) const;
};

struct StoichiometricTerm {
    std::size_t species{};
    double coefficient{}; // positive magnitude
};

// Concentration-form Kc in the same concentration units as the reaction.
// Constant reaction enthalpy gives a van't Hoff temperature dependence.
struct EquilibriumConstant {
    double reference_value{1.0};
    double reference_temperature{298.15};
    double reaction_enthalpy{}; // J/mol
    [[nodiscard]] double value(double temperature) const;
};

struct ElementaryReaction {
    std::vector<StoichiometricTerm> reactants;
    std::vector<StoichiometricTerm> products;
    ArrheniusRate forward;
    std::optional<EquilibriumConstant> equilibrium;
};

class ReactionNetwork {
public:
    explicit ReactionNetwork(std::vector<Species> species = {});

    [[nodiscard]] const std::vector<Species>& species() const noexcept { return species_; }
    [[nodiscard]] const std::vector<ElementaryReaction>& reactions() const noexcept { return reactions_; }
    [[nodiscard]] std::size_t add_species(Species species);
    void add_reaction(ElementaryReaction reaction);

    [[nodiscard]] std::vector<double> reaction_rates(std::span<const double> concentrations,
                                                      double temperature) const;
    [[nodiscard]] std::vector<double> source_terms(std::span<const double> concentrations,
                                                   double temperature) const;
    // Allocation-free overloads; buffers must be separate from concentrations.
    void reaction_rates(std::span<const double> concentrations,double temperature,std::span<double> rates) const;
    void source_terms(std::span<const double> concentrations,double temperature,
                      std::span<double> source,std::span<double> rate_workspace) const;

private:
    std::vector<Species> species_;
    std::vector<ElementaryReaction> reactions_;
    void validate_reaction(const ElementaryReaction& reaction) const;
};

} // namespace cfd::chemistry
