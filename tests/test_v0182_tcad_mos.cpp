#include "cfd/tcad/mos.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>

namespace {

void require(bool condition, const char* message) {
    if (!condition) {
        std::cerr << "v0.18.2 TCAD MOS regression failed: " << message << '\n';
        std::exit(1);
    }
}

bool near(double a, double b, double rel = 1.0e-6, double abs = 1.0e-15) {
    return std::abs(a - b) <= std::max(abs, rel * std::max(std::abs(a), std::abs(b)));
}

} // namespace

int main() {
    using namespace cfd::tcad;

    require(near(fermi_dirac_half(-10.0), std::exp(-10.0), 1.0e-8),
            "Fermi-Dirac nondegenerate limit");
    require(near(fermi_dirac_half(0.0), 0.7651470246, 1.0e-4),
            "Fermi-Dirac half integral at eta=0");
    for (double eta : {-5.0, -1.0, 0.0, 2.0, 8.0}) {
        require(near(inverse_fermi_dirac_half(fermi_dirac_half(eta)), eta, 2.0e-6, 2.0e-6),
                "inverse Fermi-Dirac round trip");
    }
    const double donors_low_ef = ionized_donor_density(1.0e22, -6.0, 0.045, 2.0, 300.0);
    const double donors_high_ef = ionized_donor_density(1.0e22, 2.0, 0.045, 2.0, 300.0);
    require(donors_low_ef > donors_high_ef && donors_low_ef > 0.9e22,
            "donor incomplete-ionization trend");
    const double acceptors_low_hole_ef = ionized_acceptor_density(1.0e22, -6.0, 0.045, 4.0, 300.0);
    const double acceptors_high_hole_ef = ionized_acceptor_density(1.0e22, 2.0, 0.045, 4.0, 300.0);
    require(acceptors_low_hole_ef > acceptors_high_hole_ef,
            "acceptor incomplete-ionization trend");

    Semiconductor1DConfig blocking_cfg;
    blocking_cfg.length_m = 1.0e-6;
    blocking_cfg.nodes = 81U;
    blocking_cfg.max_gummel_iterations = 1200U;
    blocking_cfg.relative_tolerance = 1.0e-8;
    blocking_cfg.under_relaxation = 0.10;
    blocking_cfg.left.type = ContactType::gate;
    blocking_cfg.left.voltage_v = 0.01;
    blocking_cfg.left.oxide_capacitance_f_per_m2 = 3.45e-3;
    blocking_cfg.right.type = ContactType::ohmic;
    SemiconductorDevice1D gated(blocking_cfg);
    gated.set_net_doping([](double) { return 1.0e20; });
    gated.initialize_charge_neutral();
    const auto gated_dc = gated.solve_dc();
    require(gated_dc.converged, "gate-boundary drift-diffusion convergence");
    require(std::abs(gated.electron_current_density_faces().front()) < 1.0e-20
            && std::abs(gated.hole_current_density_faces().front()) < 1.0e-20,
            "gate boundary blocks carrier current");
    require(std::isfinite(gated.left_terminal_charge_c()), "gate oxide charge is finite");

    auto insulated_cfg = blocking_cfg;
    insulated_cfg.left.type = ContactType::insulating;
    insulated_cfg.left.voltage_v = 0.0;
    insulated_cfg.left.oxide_capacitance_f_per_m2 = 0.0;
    SemiconductorDevice1D insulated(insulated_cfg);
    insulated.set_net_doping([](double) { return 1.0e20; });
    insulated.initialize_charge_neutral();
    const auto insulated_dc = insulated.solve_dc();
    require(insulated_dc.converged, "insulating-boundary convergence");
    require(std::abs(insulated.electron_current_density_faces().front()) < 1.0e-20
            && std::abs(insulated.hole_current_density_faces().front()) < 1.0e-20,
            "insulating boundary blocks carrier current");

    MosCapacitorConfig moscap_cfg;
    moscap_cfg.semiconductor_nodes = 161U;
    moscap_cfg.substrate_acceptor_m3 = 1.0e22;
    moscap_cfg.oxide_thickness_m = 10.0e-9;
    moscap_cfg.area_m2 = 1.0e-10;
    MosCapacitor1D moscap(moscap_cfg);
    const auto accumulation = moscap.solve(-1.0, 1.0e-3);
    const auto depletion = moscap.solve(0.5, 1.0e-3);
    const auto inversion = moscap.solve(1.0, 1.0e-3);
    const double cox_total = moscap.oxide_capacitance_per_area() * moscap_cfg.area_m2;
    require(accumulation.converged && depletion.converged && inversion.converged,
            "MOS capacitor operating points");
    require(accumulation.capacitance_f > 0.85 * cox_total,
            "MOS accumulation capacitance approaches Cox");
    require(depletion.capacitance_f < 0.5 * accumulation.capacitance_f,
            "MOS depletion capacitance decreases");
    require(inversion.capacitance_f > depletion.capacitance_f,
            "quasi-static inversion capacitance recovery");
    require(accumulation.surface_potential_v < 0.0 && inversion.surface_potential_v > 0.0,
            "MOS surface-potential polarity");
    const auto cv = moscap.sweep_cv(-1.0, 1.0, 5U, 1.0e-3);
    require(cv.size() == 5U && std::all_of(cv.begin(), cv.end(), [](const MosCapacitorPoint& p) {
        return p.converged && p.capacitance_f > 0.0;
    }), "MOS capacitor C-V sweep");

    LongChannelMosfetConfig mos_cfg;
    mos_cfg.width_m = 10.0e-6;
    mos_cfg.length_m = 1.0e-6;
    mos_cfg.oxide_capacitance_f_per_m2 = 3.45e-3;
    mos_cfg.electron_mobility_m2_per_vs = 0.05;
    mos_cfg.threshold_voltage_v = 0.5;
    mos_cfg.channel_length_modulation_per_v = 0.02;
    LongChannelMosfet mosfet(mos_cfg);
    require(mosfet.drain_current(0.4, 1.0) == 0.0, "MOSFET cutoff");
    const double beta = mos_cfg.electron_mobility_m2_per_vs * mos_cfg.oxide_capacitance_f_per_m2
                      * mos_cfg.width_m / mos_cfg.length_m;
    const double expected_linear = beta * (0.7 * 0.05 - 0.5 * 0.05 * 0.05);
    require(near(mosfet.drain_current(1.2, 0.05), expected_linear, 1.0e-12),
            "MOSFET integrated channel drift current");
    const double expected_sat = 0.5 * beta * 0.7 * 0.7 * (1.0 + 0.02 * (1.0 - 0.7));
    require(near(mosfet.drain_current(1.2, 1.0), expected_sat, 1.0e-12),
            "MOSFET saturation current");
    require(near(mosfet.drain_current(1.2, -0.1),
                 -mosfet.drain_current(1.3, 0.1, 0.1), 1.0e-12),
            "MOSFET drain/source reversal symmetry");
    const auto op = mosfet.operating_point(1.2, 0.1);
    require(op.transconductance_s > 0.0 && op.output_conductance_s > 0.0,
            "MOSFET differential conductances");

    const auto evaluator = make_long_channel_mosfet_spice_evaluator(mos_cfg);
    const double terminal_voltage_raw[4]{1.0, 1.2, 0.0, 0.0};
    const auto spice_eval = evaluator(std::span<const double>(terminal_voltage_raw, 4U));
    require(spice_eval.terminal_current.size() == 4U && spice_eval.jacobian.size() == 16U,
            "TCAD-SPICE evaluator dimensions");
    require(std::abs(std::accumulate(spice_eval.terminal_current.begin(), spice_eval.terminal_current.end(), 0.0)) < 1.0e-15,
            "TCAD-SPICE evaluator KCL");
    for (std::size_t col = 0U; col < 4U; ++col) {
        double column_sum = 0.0;
        for (std::size_t row = 0U; row < 4U; ++row) column_sum += spice_eval.jacobian[row * 4U + col];
        require(std::abs(column_sum) < 1.0e-9, "TCAD-SPICE Jacobian KCL");
    }

    cfd::circuit::Circuit circuit;
    const auto drain = circuit.node("d");
    const auto gate = circuit.node("g");
    const auto supply = circuit.node("vdd");
    circuit.add_voltage_source("VDD", supply, 0U, 1.0);
    circuit.add_voltage_source("VG", gate, 0U, 1.2);
    circuit.add_resistor("RD", supply, drain, 1000.0);
    circuit.add_static_device("MTCAD", {drain, gate, 0U, 0U}, evaluator);
    const auto circuit_dc = circuit.dc_operating_point_homotopy();
    require(circuit_dc.converged, "TCAD compact evaluator solves inside SPICE");
    require(circuit_dc.node_voltage[drain] > 0.0 && circuit_dc.node_voltage[drain] < 1.0,
            "TCAD MOSFET loads SPICE drain node");

    std::cout << "v0.18.2 TCAD MOS passed: Cacc/Cox="
              << accumulation.capacitance_f / cox_total
              << ", Cdep=" << depletion.capacitance_f
              << " F, Idsat=" << mosfet.drain_current(1.2, 1.0)
              << " A, Vd(spice)=" << circuit_dc.node_voltage[drain] << " V\n";
    return 0;
}
