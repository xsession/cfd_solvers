#include "cfd/chemistry/database_io.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void test_phreeqc_master_species() {
    const auto species = cfd::chemistry::load_phreeqc_master_species("SOLUTION_MASTER_SPECIES\n"
                                                                     "H H+ 1.00794 1.0\n"
                                                                     "O H2O 18.01528 0.0\n"
                                                                     "SOLUTION_SPECIES\n");
    require(species.size() == 2U, "PHREEQC master species count");
    require(species[0].name == "H+" && std::abs(species[0].molar_mass - 1.00794e-3) < 1.0e-12,
            "PHREEQC molar mass conversion");
}

void test_cantera_yaml() {
    const auto network = cfd::chemistry::load_cantera_yaml("species:\n"
                                                           "- name: H2\n"
                                                           "- name: O2\n"
                                                           "- name: H2O\n"
                                                           "reactions:\n"
                                                           "- equation: 2 H2 + O2 => 2 H2O\n"
                                                           "  rate-constant: {A: 2.0e3, b: 0.0, Ea: 1.0e4}\n");
    require(network.species().size() == 3U && network.reactions().size() == 1U, "Cantera YAML mechanism shape");
    const std::vector<double> concentrations{2.0, 1.0, 0.0};
    const auto rates = network.reaction_rates(concentrations, 1000.0);
    require(rates.size() == 1U && rates.front() > 0.0, "Cantera Arrhenius reaction rate");
}

void test_reaktoro_seam() {
    cfd::chemistry::ReaktoroDynamicAdapter adapter([](std::span<const double> composition, double, double) {
        return std::vector<double>(composition.begin(), composition.end());
    });
    const std::vector<double> composition{0.25, 0.75};
    const auto result = adapter.equilibrate(composition, 298.15, 101325.0);
    require(result == composition, "Reaktoro callback seam preserves adapted result");
}

} // namespace

int main() {
    try {
        test_phreeqc_master_species();
        test_cantera_yaml();
        test_reaktoro_seam();
        std::cout << "chemistry database/adaptor regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "chemistry database regression failed: " << error.what() << '\n';
        return 1;
    }
}
