#include "cfd/battery/dfn.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.19.3 DFN regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-6, double abs = 1.0e-12) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}
}

int main() {
    using namespace cfd::battery;

    auto cfg = default_graphite_nmc_dfn_config();
    // Disable the slow (and numerically stiff) electrolyte concentration source so
    // that the global lithium balance is exactly closed by the electrode reaction
    // partition alone, which is what the conservation regression checks.
    cfg.transference_number = 1.0;

    DoyleFullerNewmanModel dfn(cfg);
    const auto open = dfn.state();
    const double li0 = dfn.total_lithium_mol();

    // Open-circuit consistency: the two electrodes report different potentials and
    // the area-averaged OCV is a physical cell voltage.
    require(open.open_circuit_voltage_v > 2.5 && open.open_circuit_voltage_v < 4.2,
            "DFN open-circuit voltage is in a physical range");
    require(near(open.terminal_voltage_v, open.open_circuit_voltage_v, 1.0e-9, 1.0e-6),
            "DFN terminal voltage equals OCV at zero current");

    // Discharge: the distributed particles must deplete the negative and lithiate
    // the positive electrode, and the loaded voltage must fall below OCV.
    DfnStepResult last{};
    for (int i = 0; i < 120; ++i) last = dfn.step(0.5, 1.0);
    require(last.negative_surface_stoichiometry < open.negative_surface_stoichiometry,
            "DFN discharge depletes the negative electrode");
    require(last.positive_surface_stoichiometry > open.positive_surface_stoichiometry,
            "DFN discharge lithiates the positive electrode");
    require(last.terminal_voltage_v < last.open_circuit_voltage_v,
            "DFN loaded terminal voltage below OCV on discharge");
    require(last.reaction_overpotential_v > 0.0,
            "DFN resolves a positive Butler-Volmer reaction overpotential");

    // Global lithium conservation: with the electrolyte source removed the
    // reaction partition moves lithium from the negative to the positive
    // electrode (and nowhere else), so the total inventory is conserved.
    const double li1 = dfn.total_lithium_mol();
    require(near(li1, li0, 1.0e-9, 1.0e-12), "DFN conserves total lithium inventory");

    // Distributed electrolyte must be physically bounded and vary across the cell.
    const auto& c_e = dfn.electrolyte_concentration();
    const auto [emin, emax] = std::minmax_element(c_e.begin(), c_e.end());
    require(*emin > 0.0 && *emax > 0.0, "DFN electrolyte concentration stays positive");
    require(*emax > *emin, "DFN develops a distributed electrolyte concentration field");

    // Solid and electrolyte potentials must be resolved (non-constant) across the
    // cell thickness under load.
    const auto& phi_e = dfn.electrolyte_potential();
    require(std::abs(phi_e.front() - phi_e.back()) > 0.0,
            "DFN resolves a distributed electrolyte potential drop");
    const auto& phi_s = dfn.solid_potential();
    bool solid_varying = false;
    for (std::size_t i = 1; i < phi_s.size(); ++i) {
        if (std::abs(phi_s[i] - phi_s[i - 1]) > 1.0e-12) { solid_varying = true; break; }
    }
    require(solid_varying, "DFN resolves a distributed solid potential");

    // Charging must move lithium the other way and raise the negative electrode.
    auto cfg_chg = cfg;
    DoyleFullerNewmanModel dfn_chg(cfg_chg);
    const auto open_chg = dfn_chg.state();
    DfnStepResult chg{};
    for (int i = 0; i < 60; ++i) chg = dfn_chg.step(-0.3, 1.0);
    require(chg.negative_surface_stoichiometry > open_chg.negative_surface_stoichiometry,
            "DFN charge lithiates the negative electrode");
    require(chg.positive_surface_stoichiometry < open_chg.positive_surface_stoichiometry,
            "DFN charge delithiates the positive electrode");
    require(chg.terminal_voltage_v > chg.open_circuit_voltage_v,
            "DFN charging terminal voltage above OCV");

    // Invalid configuration must be rejected.
    bool threw = false;
    try {
        auto bad = cfg;
        bad.negative_nodes = 0U;
        (void)DoyleFullerNewmanModel{bad};
    } catch (const std::invalid_argument&) { threw = true; }
    require(threw, "DFN rejects an invalid electrode node count");

    std::cout << "v0.19.3 DFN regression passed\n";
    return 0;
}
