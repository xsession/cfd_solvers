#include "cfd/tcad/semiconductor1d.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <span>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.18.0 TCAD regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-8, double abs = 1.0e-12) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}

double trapezoid(std::span<const double> values, double dx) {
    double sum = 0.0;
    for (std::size_t i = 0U; i < values.size(); ++i) {
        sum += (i == 0U || i + 1U == values.size()) ? 0.5 * values[i] : values[i];
    }
    return sum * dx;
}

} // namespace

int main() {
    using namespace cfd::tcad;

    require(near(bernoulli_function(0.0), 1.0), "Bernoulli B(0)");
    const double xb = 0.73;
    require(near(bernoulli_function(-xb), std::exp(xb) * bernoulli_function(xb), 1.0e-12),
            "Bernoulli detailed-balance identity");

    RecombinationModel recomb;
    recomb.tau_n_s = 1.0e-6;
    recomb.tau_p_s = 2.0e-6;
    recomb.auger_n_m6_per_s = 2.0e-43;
    recomb.auger_p_m6_per_s = 1.0e-43;
    recomb.radiative_m3_per_s = 1.0e-16;
    const double ni = 1.0e16;
    const double req = total_recombination(ni, ni, ni, recomb);
    require(std::abs(req) < 1.0e-12,
            "equilibrium recombination must vanish");
    require(total_recombination(5.0e20, 2.0e20, ni, recomb) > 0.0,
            "excess carriers must recombine");

    const double mu_low_doping = caughey_thomas_mobility(1.0e20, 0.02, 0.14, 1.0e23, 0.7);
    const double mu_high_doping = caughey_thomas_mobility(1.0e25, 0.02, 0.14, 1.0e23, 0.7);
    require(mu_low_doping > mu_high_doping, "doping mobility degradation");
    require(high_field_mobility(0.14, 1.0e7, 1.0e5) < high_field_mobility(0.14, 1.0e3, 1.0e5),
            "high-field velocity saturation");

    Semiconductor1DConfig pn_cfg;
    pn_cfg.length_m = 2.0e-6;
    pn_cfg.nodes = 101U;
    pn_cfg.material.intrinsic_density_m3 = ni;
    pn_cfg.max_gummel_iterations = 1000U;
    pn_cfg.relative_tolerance = 2.0e-7;
    pn_cfg.under_relaxation = 0.05;
    SemiconductorDevice1D pn(pn_cfg);
    constexpr double doping = 1.0e21;
    pn.set_net_doping(pn_junction_doping(pn.x(), 1.0e-6, doping, doping));
    const auto equilibrium = pn.solve_equilibrium();
    require(equilibrium.converged, "PN equilibrium convergence");
    const double vt = thermal_voltage(pn_cfg.material.temperature_k);
    const double expected_vbi = 2.0 * vt * std::log(doping / ni);
    const double actual_vbi = pn.potential().back() - pn.potential().front();
    require(near(actual_vbi, expected_vbi, 2.0e-5), "PN built-in voltage");
    require(std::abs(equilibrium.terminal_current_a) < 1.0e-12, "PN equilibrium current");
    require(pn.electron_density().front() < pn.hole_density().front(), "P-side majority holes");
    require(pn.electron_density().back() > pn.hole_density().back(), "N-side majority electrons");

    const auto pin_profile = pin_diode_doping(pn.x(), 0.6e-6, 1.4e-6, doping, doping);
    require(pin_profile.front() < 0.0 && std::abs(pin_profile[50U]) < 1.0 && pin_profile.back() > 0.0,
            "PIN doping regions");

    Semiconductor1DConfig resistor_cfg;
    resistor_cfg.length_m = 1.0e-6;
    resistor_cfg.nodes = 51U;
    resistor_cfg.area_m2 = 1.0e-12;
    resistor_cfg.material.intrinsic_density_m3 = ni;
    resistor_cfg.right.voltage_v = -0.01;
    resistor_cfg.max_gummel_iterations = 300U;
    resistor_cfg.relative_tolerance = 1.0e-8;
    resistor_cfg.under_relaxation = 0.25;
    SemiconductorDevice1D resistor(resistor_cfg);
    resistor.set_net_doping([](double) { return doping; });
    resistor.initialize_charge_neutral();
    const auto dc = resistor.solve_dc();
    require(dc.converged, "uniform doped-bar DC convergence");
    const double electric_field = 0.01 / resistor_cfg.length_m;
    const double expected_j = elementary_charge_c * resistor_cfg.material.electron_mobility_m2_per_vs
                            * doping * electric_field;
    const auto jfaces = resistor.total_current_density_faces();
    const auto [jmin, jmax] = std::minmax_element(jfaces.begin(), jfaces.end());
    require(near(jfaces.front(), expected_j, 2.0e-4), "Scharfetter-Gummel drift current");
    require((*jmax - *jmin) / expected_j < 2.0e-8, "steady current continuity");
    require(near(dc.terminal_current_a, expected_j * resistor_cfg.area_m2, 2.0e-4), "terminal current");

    const auto heat = resistor.joule_heating_cells();
    const double dx = resistor_cfg.length_m / static_cast<double>(resistor_cfg.nodes - 1U);
    const double heat_power = trapezoid(heat, dx) * resistor_cfg.area_m2;
    require(near(heat_power, std::abs(dc.terminal_current_a * resistor_cfg.right.voltage_v), 3.0e-4),
            "electrothermal Joule power conservation");

    const auto sweep = resistor.sweep_right_contact(-0.005, -0.015, 3U);
    require(sweep.size() == 3U && std::all_of(sweep.begin(), sweep.end(), [](const IvPoint& p) { return p.converged; }),
            "DC voltage sweep convergence");
    require(std::abs(sweep[2].current_a) > std::abs(sweep[0].current_a), "ohmic IV monotonicity");

    Semiconductor1DConfig schottky_cfg;
    schottky_cfg.nodes = 11U;
    schottky_cfg.left.type = ContactType::schottky;
    schottky_cfg.left.electron_barrier_ev = 0.72;
    SemiconductorDevice1D schottky(schottky_cfg);
    schottky.set_net_doping([](double) { return doping; });
    schottky.initialize_charge_neutral();
    const double expected_schottky_n = schottky_cfg.material.conduction_band_density_m3
        * std::exp(-schottky_cfg.left.electron_barrier_ev / thermal_voltage(schottky_cfg.material.temperature_k));
    require(near(schottky.electron_density().front(), expected_schottky_n, 1.0e-12),
            "Schottky Boltzmann boundary density");

    std::cout << "v0.18.0 TCAD baseline regression passed: Vbi=" << actual_vbi
              << " V, J=" << jfaces.front() << " A/m^2, P=" << heat_power << " W\n";
    return 0;
}
