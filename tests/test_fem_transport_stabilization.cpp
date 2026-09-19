#include "cfd/solvers/fem/transport_stabilization.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const auto supg = cfd::fem::supg_parameters(0.1, 2.0, 1.0e-3, 1.0e-2);
        if (!(supg.tau_s > 0.0) || !(supg.element_peclet > 1.0))
            throw std::runtime_error("SUPG parameter baseline");
        const double positive = cfd::fem::dg_upwind_flux(2.0, 1.0, 3.0, 0.5);
        const double negative = cfd::fem::dg_upwind_flux(2.0, 1.0, -3.0, 0.5);
        if (!(positive > 0.0) || !(negative < 0.0) || !std::isfinite(positive + negative))
            throw std::runtime_error("DG upwind flux baseline");
        std::cout << "SUPG/DG FEM transport stabilization regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SUPG/DG FEM stabilization regression failed: " << error.what() << '\n';
        return 1;
    }
}
