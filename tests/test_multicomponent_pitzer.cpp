#include "cfd/chemistry/advanced_thermo.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    cfd::chemistry::MulticomponentPitzer model({{0U, 1U, 0.08, 0.02, 0.001}, {1U, 2U, 0.04, 0.01, 0.0005}});
    const int charge[] = {1, -1, 2};
    const double molality[] = {0.1, 0.1, 0.02};
    const auto gamma = model.activity_coefficients(charge, molality, 298.15);
    assert(gamma.size() == 3U);
    for (double value : gamma)
        assert(value > 0.0 && std::isfinite(value));
    std::cout << "multicomponent Pitzer regression passed\n";
    return 0;
}
