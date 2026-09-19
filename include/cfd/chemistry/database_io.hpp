#pragma once

#include "cfd/chemistry/kinetics.hpp"

#include <functional>
#include <span>
#include <string_view>
#include <vector>

namespace cfd::chemistry {

// PHREEQC SOLUTION_MASTER_SPECIES subset. Molar masses are returned in kg/mol
// even when the source database uses g/mol.
[[nodiscard]] std::vector<Species> load_phreeqc_master_species(std::string_view database_text);

// Dependency-free Cantera YAML subset covering species names and elementary
// Arrhenius reactions. Vendor-specific YAML keys are ignored deliberately.
[[nodiscard]] ReactionNetwork load_cantera_yaml(std::string_view mechanism_text);

// Optional Reaktoro seam. The core stores no Reaktoro types or symbols; an
// application can adapt its own dynamic library through this evaluator.
class ReaktoroDynamicAdapter {
public:
    using Evaluator = std::function<std::vector<double>(std::span<const double>, double, double)>;

    explicit ReaktoroDynamicAdapter(Evaluator evaluator);
    [[nodiscard]] std::vector<double> equilibrate(std::span<const double> composition, double temperature,
                                                  double pressure) const;

private:
    Evaluator evaluator_;
};

} // namespace cfd::chemistry
