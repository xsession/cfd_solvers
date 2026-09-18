#pragma once

#include "cfd/circuit/spice.hpp"
#include "cfd/tcad/semiconductor1d.hpp"

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::tcad {

// Complete Fermi-Dirac integral F_{1/2}(eta), normalized so that
// F_{1/2}(eta) -> exp(eta) in the nondegenerate limit.
[[nodiscard]] double fermi_dirac_half(double eta);
[[nodiscard]] double inverse_fermi_dirac_half(double value);
[[nodiscard]] double fermi_electron_density(double conduction_band_density_m3,
                                            double reduced_fermi_to_conduction);
[[nodiscard]] double fermi_hole_density(double valence_band_density_m3,
                                        double reduced_fermi_from_valence);

// Incomplete-ionization helpers. reduced_fermi_to_conduction is
// (E_F-E_C)/(kT); reduced_fermi_from_valence is (E_V-E_F)/(kT).
[[nodiscard]] double ionized_donor_density(double total_donor_m3,
                                           double reduced_fermi_to_conduction,
                                           double donor_binding_ev,
                                           double degeneracy,
                                           double temperature_k);
[[nodiscard]] double ionized_acceptor_density(double total_acceptor_m3,
                                              double reduced_fermi_from_valence,
                                              double acceptor_binding_ev,
                                              double degeneracy,
                                              double temperature_k);

struct MosCapacitorConfig {
    double semiconductor_thickness_m{1.0e-6};
    std::size_t semiconductor_nodes{161U};
    double area_m2{1.0e-10};
    double substrate_acceptor_m3{1.0e22};
    double oxide_thickness_m{10.0e-9};
    double oxide_relative_permittivity{3.9};
    double flatband_voltage_v{0.0};
    SemiconductorMaterial semiconductor{};
    std::size_t max_iterations{1200U};
    double tolerance{1.0e-8};
    double relaxation{0.18};
};

struct MosCapacitorPoint {
    double gate_voltage_v{};
    double surface_potential_v{};
    double gate_charge_c{};
    double capacitance_f{};
    bool converged{};
};

class MosCapacitor1D {
public:
    explicit MosCapacitor1D(MosCapacitorConfig config = {});

    [[nodiscard]] MosCapacitorPoint solve(double gate_voltage_v,
                                          double capacitance_perturbation_v = 1.0e-4) const;
    [[nodiscard]] std::vector<MosCapacitorPoint> sweep_cv(double start_v,
                                                          double stop_v,
                                                          std::size_t points,
                                                          double perturbation_v = 1.0e-4) const;
    [[nodiscard]] double oxide_capacitance_per_area() const;
    [[nodiscard]] const MosCapacitorConfig& config() const noexcept { return config_; }

private:
    [[nodiscard]] SemiconductorDevice1D make_device(double gate_voltage_v) const;
    MosCapacitorConfig config_{};
};

// Long-channel charge-sheet MOSFET. The channel current is the analytic
// integral of the 1-D drift equation I = W*mu*Qinv(V)*dV/dx using
// Qinv=-Cox*(Vgs-Vth-V). It is a compact drift-diffusion device baseline,
// not a 2-D semiconductor field solve.
struct LongChannelMosfetConfig {
    double width_m{10.0e-6};
    double length_m{1.0e-6};
    double oxide_capacitance_f_per_m2{3.45e-3};
    double electron_mobility_m2_per_vs{0.05};
    double threshold_voltage_v{0.5};
    double channel_length_modulation_per_v{0.02};
    double body_effect_coefficient_sqrt_v{0.0};
    double surface_potential_2phi_f_v{0.7};
};

struct MosfetOperatingPoint {
    double vgs_v{};
    double vds_v{};
    double vbs_v{};
    double drain_current_a{};
    double transconductance_s{};
    double output_conductance_s{};
    double threshold_voltage_v{};
};

struct MosfetIvSample {
    double gate_voltage_v{};
    double drain_voltage_v{};
    double drain_current_a{};
};

class LongChannelMosfet {
public:
    explicit LongChannelMosfet(LongChannelMosfetConfig config = {});
    [[nodiscard]] double threshold_voltage(double vbs_v = 0.0) const;
    [[nodiscard]] double drain_current(double vgs_v, double vds_v, double vbs_v = 0.0) const;
    [[nodiscard]] MosfetOperatingPoint operating_point(double vgs_v, double vds_v,
                                                       double vbs_v = 0.0) const;
    [[nodiscard]] std::vector<MosfetIvSample> sweep(double gate_start_v,
                                                   double gate_stop_v,
                                                   std::size_t gate_points,
                                                   double drain_start_v,
                                                   double drain_stop_v,
                                                   std::size_t drain_points,
                                                   double vbs_v = 0.0) const;
    [[nodiscard]] const LongChannelMosfetConfig& config() const noexcept { return config_; }

private:
    LongChannelMosfetConfig config_{};
};

// Four-terminal [D,G,S,B] static compact-device evaluator that can be passed
// directly to Circuit::add_static_device(). The Jacobian is evaluated by a
// centered finite difference of the TCAD-derived long-channel model.
[[nodiscard]] cfd::circuit::StaticDeviceEvaluator
make_long_channel_mosfet_spice_evaluator(LongChannelMosfetConfig config,
                                         double jacobian_step_v = 1.0e-6);

} // namespace cfd::tcad
