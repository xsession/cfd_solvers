#include "cfd/electrochemistry/corrosion_models.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const std::vector<double> mesh{0.0, 0.25, 0.5, 0.75, 1.0};
    const auto moved = cfd::electrochemistry::recession_remesh_1d(mesh, 0.1);
    assert(moved.size() == mesh.size() && moved.back() < mesh.back());
    const auto update =
        cfd::electrochemistry::coupled_stress_corrosion_step(1.0, 1.0e7, 1.0e-6, 300.0, 10.0, 0.055, 7800.0);
    assert(update.stress_amplified_current_density > 1.0 && update.recession_distance > 0.0);
    const auto flux =
        cfd::electrochemistry::modified_nernst_planck_flux(1.0, 0.5, 1.0e-9, 10.0, 1, 300.0, 1.0e-3, 0.01, 10.0);
    assert(flux.activity > 0.0 && std::isfinite(flux.molar_flux));
    cfd::electrochemistry::PhaseFieldCorrosion2D phase(8U, 8U, 0.1, 0.1, 0.1, 0.01);
    phase.initialize_interface(0.35, 0.35, 0.2);
    phase.step(1.0e-3, std::vector<double>(64U, 0.0));
    for (double value : phase.order_parameter())
        assert(std::isfinite(value) && value >= -1.0 && value <= 1.0);
    std::cout << "advanced corrosion coupling regression passed\n";
    return 0;
}
