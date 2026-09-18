#include "cfd/battery/lithium_ion.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.19.0 battery regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-6, double abs = 1.0e-12) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}
}

int main() {
    using namespace cfd::battery;

    // Spherical finite-volume diffusion must conserve lithium against the
    // prescribed surface molar flux.
    ParticleConfig pc;
    pc.radius_m = 5.0e-6;
    pc.diffusivity_m2_per_s = 1.0e-14;
    pc.maximum_concentration_mol_per_m3 = 3.0e4;
    pc.initial_stoichiometry = 0.6;
    pc.radial_cells = 30U;
    SphericalDiffusionParticle particle(pc);
    const double c0 = particle.average_concentration();
    const double flux = 2.0e-6;
    const double dt = 0.02;
    particle.step(flux, dt);
    const double expected_delta = -3.0 * flux * dt / pc.radius_m;
    require(near(particle.average_concentration() - c0, expected_delta, 3.0e-10, 1.0e-8),
            "spherical particle lithium conservation");
    require(particle.surface_concentration(flux) < particle.average_concentration(),
            "outward flux lowers particle surface concentration");

    auto cfg = default_graphite_nmc_config();
    cfg.thermal.heat_transfer_coefficient_w_per_m2_k = 0.0;
    SingleParticleModel spm(cfg);
    const auto initial = spm.state();
    BatteryStepResult last{};
    for (int i = 0; i < 120; ++i) last = spm.step(0.5, 1.0);
    require(last.negative_surface_stoichiometry < initial.negative_surface_stoichiometry,
            "SPM discharge depletes negative electrode");
    require(last.positive_surface_stoichiometry > initial.positive_surface_stoichiometry,
            "SPM discharge lithiates positive electrode");
    require(last.terminal_voltage_v < last.open_circuit_voltage_v,
            "SPM loaded terminal voltage below OCV on discharge");
    require(last.temperature_k >= cfg.initial_temperature_k,
            "SPM irreversible losses heat the lumped cell");

    SingleParticleElectrolyteModel spme(cfg);
    BatteryStepResult spme_last{};
    for (int i = 0; i < 120; ++i) spme_last = spme.step(0.5, 1.0);
    const auto [emin, emax] = std::minmax_element(spme.electrolyte_concentration().begin(), spme.electrolyte_concentration().end());
    require(*emin > 0.0 && *emax > *emin,
            "SPMe develops a positive electrolyte concentration gradient");
    require(std::abs(spme_last.electrolyte_overpotential_v) > 0.0,
            "SPMe has electrolyte polarization");
    require(spme.electrolyte_potential().front() != spme.electrolyte_potential().back(),
            "SPMe reconstructs electrolyte potential drop");

    ThermalSlabConfig tc;
    tc.cells = 20U;
    tc.conductivity_w_per_m_k = 0.0;
    tc.convection_w_per_m2_k = 0.0;
    ThermalSlab1D thermal(tc);
    thermal.reset(300.0);
    std::vector<double> q(tc.cells, 2.0e5);
    thermal.step(2.0, q);
    const double expected_t = 300.0 + 2.0 * 2.0e5 / (tc.density_kg_per_m3 * tc.heat_capacity_j_per_kg_k);
    require(near(thermal.average_temperature(), expected_t, 1.0e-11),
            "1D thermal slab uniform heating energy balance");

    auto degradation_cfg = cfg;
    degradation_cfg.degradation.sei_rate_m2_per_s = 1.0e-17;
    degradation_cfg.degradation.plating_efficiency = 0.1;
    degradation_cfg.degradation.active_material_loss_per_coulomb = 1.0e-6;
    degradation_cfg.degradation.cracking_rate_per_gradient_second = 1.0e-2;
    SingleParticleModel ageing(degradation_cfg);
    const auto d0 = ageing.degradation_state();
    for (int i = 0; i < 30; ++i) (void)ageing.step(-0.5, 1.0);
    const auto d1 = ageing.degradation_state();
    require(d1.sei_thickness_m > d0.sei_thickness_m, "SEI thickness grows");
    require(d1.plated_lithium_mol > 0.0, "charging can accumulate plated lithium");
    require(d1.active_material_fraction < d0.active_material_fraction, "active-material fraction fades");
    require(d1.crack_damage >= d0.crack_damage, "particle crack damage is monotonic");

    const auto z_low = battery_impedance(0.01);
    const auto z_high = battery_impedance(1000.0);
    require(std::real(z_low) > std::real(z_high), "battery EIS resolves low-frequency polarization");
    require(std::imag(z_low) < 0.0 && std::imag(z_high) < 0.0,
            "battery EIS is capacitive/diffusive over test frequencies");

    SingleParticleModel cycle_model(cfg);
    const std::vector<DriveCyclePoint> profile{{60.0, 0.5}, {30.0, 0.0}, {60.0, -0.25}};
    const auto cycle = simulate_drive_cycle(cycle_model, profile, 2.0);
    require(cycle.samples.size() == 75U, "drive-cycle fixed maximum step sampling");
    require(near(cycle.discharged_capacity_ah, 0.5 * 60.0 / 3600.0, 1.0e-12),
            "drive-cycle discharge capacity accounting");
    require(std::isfinite(cycle.electrical_energy_wh), "drive-cycle energy accounting finite");

    std::cout << "v0.19.0 battery regression passed\n";
    return 0;
}
