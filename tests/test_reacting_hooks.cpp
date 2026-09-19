#include "cfd/solvers/fvm/flame_models.hpp"
#include "cfd/solvers/fvm/soot_radiation.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void test_flame_regimes() {
    using namespace cfd::fvm;
    PremixedFlameModel premixed({4.0, 1.0e6, 500.0});
    require(premixed.evaluate({0.4, 0.5, 400.0}).progress_source == 0.0,
            "premixed model must respect ignition threshold");
    const auto premixed_hot = premixed.evaluate({0.4, 0.5, 900.0});
    require(premixed_hot.progress_source > 0.0 && premixed_hot.heat_release > 0.0,
            "premixed model must release heat in the flame zone");

    NonPremixedFlameModel non_premixed({0.4, 0.05, 10.0, 1.0e6, 300.0});
    const auto at_stoich = non_premixed.evaluate({0.2, 0.4, 900.0});
    const auto off_stoich = non_premixed.evaluate({0.2, 0.1, 900.0});
    require(at_stoich.progress_source > off_stoich.progress_source, "non-premixed source must peak near stoichiometry");
}

void test_soot_radiation_hook() {
    using namespace cfd::fvm;
    SootRadiationModel model({0.02, 0.5, 2.0e4, 0.0});
    const auto result = model.evaluate({0.01, 1200.0, 1.0e6}, 300.0);
    require(result.soot_source > 0.0, "soot source must include production");
    require(result.absorption_coefficient > 0.0, "soot must produce absorption");
    require(result.radiative_source < 0.0, "hot soot must radiate energy from the cell");
}

} // namespace

int main() {
    try {
        test_flame_regimes();
        test_soot_radiation_hook();
        std::cout << "reacting flame and soot/radiation regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "reacting hooks regression failed: " << error.what() << '\n';
        return 1;
    }
}
