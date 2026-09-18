#include "cfd/tcad/bipolar_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.18.3 TCAD bipolar/power regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-6, double abs = 1.0e-15) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}

} // namespace

int main() {
    using namespace cfd::tcad;

    BipolarJunctionConfig bjt_cfg;
    bjt_cfg.area_m2 = 1.0e-10;
    bjt_cfg.emitter_width_m = 0.5e-6;
    bjt_cfg.base_width_m = 0.25e-6;
    bjt_cfg.collector_width_m = 2.0e-6;
    bjt_cfg.emitter_doping_m3 = 1.0e24;
    bjt_cfg.base_doping_m3 = 2.0e22;
    bjt_cfg.collector_doping_m3 = 2.0e21;
    bjt_cfg.minority_electron_lifetime_base_s = 2.0e-7;
    bjt_cfg.minority_hole_lifetime_emitter_s = 5.0e-8;
    bjt_cfg.minority_hole_lifetime_collector_s = 2.0e-7;

    BipolarJunctionTransistor1D npn(bjt_cfg);
    require(npn.net_doping_m3(0.1e-6) > 0.0,
            "NPN emitter is n-type");
    require(npn.net_doping_m3(0.6e-6) < 0.0,
            "NPN base is p-type");
    require(npn.net_doping_m3(1.5e-6) > 0.0,
            "NPN collector is n-type");
    const auto profile = npn.doping_profile(101U);
    require(profile.size() == 101U && profile.front() > 0.0 && profile.back() > 0.0,
            "NPN profile generation");

    const auto derived = npn.derived();
    require(derived.forward_alpha > 0.9 && derived.forward_alpha < 1.0,
            "BJT forward transport factor");
    require(derived.forward_beta > 50.0,
            "BJT geometry produces useful forward gain");
    require(derived.reverse_beta > 0.1 && derived.reverse_beta < derived.forward_beta,
            "BJT asymmetric reverse gain");
    require(derived.emitter_base_built_in_v > derived.collector_base_built_in_v,
            "BJT emitter/base junction reflects heavier emitter doping");
    require(near(derived.forward_alpha * derived.emitter_saturation_current_a,
                 derived.reverse_alpha * derived.collector_saturation_current_a,
                 2.0e-12, 1.0e-30),
            "BJT reciprocity alphaF*IES=alphaR*ICS");

    const auto forward = npn.operating_point(5.0, 0.65, 0.0);
    require(forward.collector_current_a > 0.0 && forward.base_current_a > 0.0
            && forward.emitter_current_a < 0.0,
            "NPN forward-active current signs");
    require(std::abs(forward.collector_current_a / forward.base_current_a
                     - derived.forward_beta) < 0.03 * derived.forward_beta,
            "NPN forward-active beta follows derived diffusion gain");
    require(std::abs(forward.collector_current_a + forward.base_current_a
                     + forward.emitter_current_a) < 1.0e-14,
            "BJT terminal KCL");
    require(forward.dissipated_power_w > 0.0,
            "BJT dissipates positive power in forward active operation");

    auto pnp_cfg = bjt_cfg;
    pnp_cfg.polarity = BipolarPolarity::pnp;
    BipolarJunctionTransistor1D pnp(pnp_cfg);
    const auto pnp_forward = pnp.operating_point(0.0, 4.35, 5.0);
    require(pnp_forward.collector_current_a < 0.0 && pnp_forward.base_current_a < 0.0
            && pnp_forward.emitter_current_a > 0.0,
            "PNP polarity symmetry");

    const auto bjt_evaluator = make_bipolar_tcad_spice_evaluator(bjt_cfg);
    const double bjt_v[3]{5.0, 0.65, 0.0};
    const auto bjt_eval = bjt_evaluator(std::span<const double>(bjt_v, 3U));
    require(bjt_eval.terminal_current.size() == 3U && bjt_eval.jacobian.size() == 9U,
            "BJT SPICE evaluator dimensions");
    require(std::abs(std::accumulate(bjt_eval.terminal_current.begin(), bjt_eval.terminal_current.end(), 0.0)) < 1.0e-14,
            "BJT SPICE evaluator KCL");
    for (std::size_t col = 0U; col < 3U; ++col) {
        double sum = 0.0;
        for (std::size_t row = 0U; row < 3U; ++row) sum += bjt_eval.jacobian[row * 3U + col];
        require(std::abs(sum) < 1.0e-8, "BJT SPICE Jacobian KCL");
    }

    cfd::circuit::Circuit bjt_circuit;
    const auto bjt_vcc = bjt_circuit.node("vcc");
    const auto bjt_base = bjt_circuit.node("b");
    const auto bjt_collector = bjt_circuit.node("c");
    bjt_circuit.add_voltage_source("VCC", bjt_vcc, 0U, 5.0);
    bjt_circuit.add_voltage_source("VB", bjt_base, 0U, 0.62);
    bjt_circuit.add_resistor("RC", bjt_vcc, bjt_collector, 1000.0);
    bjt_circuit.add_static_device("QTCAD", {bjt_collector, bjt_base, 0U}, bjt_evaluator);
    const auto bjt_dc = bjt_circuit.dc_operating_point_homotopy();
    require(bjt_dc.converged && bjt_dc.node_voltage[bjt_collector] > 0.0
            && bjt_dc.node_voltage[bjt_collector] < 5.0,
            "TCAD-derived BJT converges inside SPICE");

    IgbtConfig igbt_cfg;
    igbt_cfg.channel.width_m = 100.0e-6;
    igbt_cfg.channel.length_m = 1.0e-6;
    igbt_cfg.channel.oxide_capacitance_f_per_m2 = 3.45e-3;
    igbt_cfg.channel.electron_mobility_m2_per_vs = 0.05;
    igbt_cfg.channel.threshold_voltage_v = 0.7;
    igbt_cfg.channel.channel_length_modulation_per_v = 0.01;
    igbt_cfg.bipolar = pnp_cfg;
    igbt_cfg.maximum_bipolar_gain = 5.0;
    igbt_cfg.drift_resistance_ohm = 0.08;
    igbt_cfg.conductivity_modulation_current_a = 0.02;
    igbt_cfg.junction_drop_v = 0.65;
    igbt_cfg.tail_lifetime_s = 2.0e-6;
    igbt_cfg.ambient_temperature_k = 300.0;
    igbt_cfg.thermal_resistance_k_per_w = 0.8;
    igbt_cfg.thermal_capacitance_j_per_k = 0.01;

    IgbtPowerDevice igbt(igbt_cfg);
    const auto off = igbt.operating_point(0.0, 3.0);
    const auto on = igbt.operating_point(4.0, 3.0);
    require(off.collector_current_a == 0.0, "IGBT gate cutoff");
    require(on.collector_current_a > 0.0 && on.channel_current_a > 0.0
            && on.bipolar_current_a > 0.0,
            "IGBT MOS-gated bipolar conduction");
    require(on.collector_current_a > on.channel_current_a,
            "IGBT bipolar conductivity modulation increases collector current");
    require(on.effective_drift_resistance_ohm < igbt_cfg.drift_resistance_ohm,
            "IGBT conductivity modulation lowers drift resistance");
    const auto hot = igbt.operating_point(4.0, 3.0, 400.0);
    require(hot.channel_current_a < on.channel_current_a,
            "IGBT channel mobility decreases with temperature");

    double last_on_tail = 0.0;
    for (std::size_t i = 0U; i < 200U; ++i) {
        const auto step = igbt.step(4.0, 3.0, 1.0e-7);
        last_on_tail = step.bipolar_current_a;
    }
    const double heated_temperature = igbt.state().junction_temperature_k;
    require(heated_temperature > igbt_cfg.ambient_temperature_k,
            "IGBT self-heating raises junction temperature");
    const auto turn_off = igbt.step(0.0, 3.0, 1.0e-7);
    require(turn_off.collector_current_a > 0.0 && turn_off.collector_current_a < last_on_tail,
            "IGBT turn-off tail current persists and decays");

    const auto igbt_evaluator = make_igbt_tcad_spice_evaluator(igbt_cfg);
    const double igbt_v[3]{3.0, 4.0, 0.0};
    const auto igbt_eval = igbt_evaluator(std::span<const double>(igbt_v, 3U));
    require(igbt_eval.terminal_current.size() == 3U && igbt_eval.jacobian.size() == 9U,
            "IGBT SPICE evaluator dimensions");
    require(std::abs(std::accumulate(igbt_eval.terminal_current.begin(), igbt_eval.terminal_current.end(), 0.0)) < 1.0e-12,
            "IGBT SPICE evaluator KCL");

    cfd::circuit::Circuit igbt_circuit;
    const auto supply = igbt_circuit.node("vdc");
    const auto collector = igbt_circuit.node("c");
    const auto gate = igbt_circuit.node("g");
    igbt_circuit.add_voltage_source("VDC", supply, 0U, 5.0);
    igbt_circuit.add_voltage_source("VG", gate, 0U, 4.0);
    igbt_circuit.add_resistor("RL", supply, collector, 20.0);
    igbt_circuit.add_static_device("IGBTTCAD", {collector, gate, 0U}, igbt_evaluator);
    const auto igbt_dc = igbt_circuit.dc_operating_point_homotopy();
    require(igbt_dc.converged && igbt_dc.node_voltage[collector] > 0.0
            && igbt_dc.node_voltage[collector] < 5.0,
            "TCAD-derived IGBT converges inside SPICE");

    std::cout << "v0.18.3 TCAD bipolar/power passed: betaF=" << derived.forward_beta
              << ", Ic(BJT)=" << forward.collector_current_a
              << " A, Ic(IGBT)=" << on.collector_current_a
              << " A, Tj=" << heated_temperature
              << " K, Vc(spice)=" << igbt_dc.node_voltage[collector] << " V\n";
    return 0;
}
