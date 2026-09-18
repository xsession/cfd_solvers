#include "cfd/tcad/semiconductor1d.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.18.1 TCAD analysis regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-7, double abs = 1.0e-15) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}

} // namespace

int main() {
    using namespace cfd::tcad;

    const double eg300 = silicon_bandgap_ev(300.0);
    const double eg400 = silicon_bandgap_ev(400.0);
    require(eg400 < eg300, "Varshni bandgap must decrease with temperature");
    const double ni400 = temperature_scaled_intrinsic_density(1.0e16, 300.0, 400.0, eg300, eg400);
    require(ni400 > 1.0e16, "intrinsic density must increase with temperature");
    require(temperature_scaled_mobility(0.135, 300.0, 400.0, 2.42) < 0.135,
            "phonon-limited mobility must decrease with temperature");

    Semiconductor1DConfig cfg;
    cfg.length_m = 1.0e-6;
    cfg.nodes = 41U;
    cfg.area_m2 = 2.0e-12;
    cfg.material.intrinsic_density_m3 = 1.0e16;
    cfg.max_gummel_iterations = 600U;
    cfg.relative_tolerance = 1.0e-9;
    cfg.under_relaxation = 0.4;
    constexpr double doping = 1.0e21;

    SemiconductorDevice1D device(cfg);
    device.set_net_doping([](double) { return doping; });
    device.initialize_charge_neutral();
    const auto zero = device.solve_dc();
    require(zero.converged, "zero-bias DC convergence");

    const double eps = vacuum_permittivity_f_per_m * cfg.material.relative_permittivity;
    const double expected_c = eps * cfg.area_m2 / cfg.length_m;
    const double expected_g = elementary_charge_c * cfg.material.electron_mobility_m2_per_vs
                            * doping * cfg.area_m2 / cfg.length_m;
    const auto ac = device.small_signal(1.0e6, 1.0e-5);
    require(ac.converged, "small-signal perturbation solves");
    require(near(ac.conductance_s, expected_g, 4.0e-4), "small-signal conductance");
    require(near(ac.capacitance_f, expected_c, 2.0e-5), "small-signal geometric capacitance");
    require(near(ac.admittance_s.imag(), 2.0 * std::acos(-1.0) * 1.0e6 * expected_c, 2.0e-5),
            "small-signal capacitive admittance");
    require(near(std::abs(ac.admittance_s * ac.impedance_ohm), 1.0, 1.0e-10),
            "admittance/impedance reciprocity");

    const auto cv = device.sweep_cv(-1.0e-3, 1.0e-3, 3U, 1.0e-5);
    require(cv.size() == 3U && std::all_of(cv.begin(), cv.end(), [](const CvPoint& p) { return p.converged; }),
            "C-V sweep convergence");
    for (const auto& point : cv) require(near(point.capacitance_f, expected_c, 3.0e-5), "uniform-bar C-V capacitance");

    constexpr double step_v = 1.0e-3;
    constexpr double dt = 1.0e-9;
    device.set_contact_voltages(step_v, 0.0);
    const auto transient = device.step_transient(dt);
    require(transient.converged, "backward-Euler transient convergence");
    require(near(transient.time_s, dt), "transient time advance");
    require(near(transient.displacement_current_a, expected_c * step_v / dt, 2.0e-4),
            "terminal displacement current");
    require(near(transient.conduction_current_a, expected_g * step_v, 5.0e-4),
            "transient conduction current");
    require(near(transient.terminal_current_a,
                 transient.conduction_current_a + transient.displacement_current_a, 1.0e-12),
            "transient total terminal current");

    SemiconductorDevice1D hot(cfg);
    hot.set_net_doping([](double) { return doping; });
    const double mu300 = hot.config().material.electron_mobility_m2_per_vs;
    const double ni300 = hot.config().material.intrinsic_density_m3;
    hot.set_temperature(400.0);
    require(hot.config().material.electron_mobility_m2_per_vs < mu300,
            "integrated temperature-dependent electron mobility");
    require(hot.config().material.intrinsic_density_m3 > ni300,
            "integrated temperature-dependent intrinsic density");
    require(hot.config().material.conduction_band_density_m3 > cfg.material.conduction_band_density_m3,
            "temperature-dependent density of states");

    std::cout << "v0.18.1 TCAD analyses passed: G=" << ac.conductance_s
              << " S, C=" << ac.capacitance_f
              << " F, Idisp=" << transient.displacement_current_a << " A\n";
    return 0;
}
