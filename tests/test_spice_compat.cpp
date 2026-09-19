#include "cfd/circuit/spice_compat.hpp"

#include <cassert>
#include <iostream>

int main() {
    for (const auto mode : {cfd::circuit::SpiceCompatibilityMode::pspice, cfd::circuit::SpiceCompatibilityMode::ltspice,
                            cfd::circuit::SpiceCompatibilityMode::hspice}) {
        const auto profile = cfd::circuit::spice_compatibility_profile(mode);
        assert(profile.engineering_suffixes && profile.common_model_cards && profile.subcircuits_and_parameters);
        assert(profile.supported_directives.size() >= 10U);
    }
    std::cout << "SPICE compatibility profile regression passed\n";
    return 0;
}
