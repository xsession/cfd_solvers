#pragma once

#include "cfd/circuit/spice.hpp"
#include "cfd/tcad/mos.hpp"
#include "cfd/tcad/semiconductor1d.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::tcad {

enum class BipolarPolarity { npn, pnp };

struct BipolarJunctionConfig {
    BipolarPolarity polarity{BipolarPolarity::npn};
    double area_m2{1.0e-10};
    double emitter_width_m{0.5e-6};
    double base_width_m{0.25e-6};
    double collector_width_m{2.0e-6};
    double emitter_doping_m3{1.0e24};
    double base_doping_m3{2.0e22};
    double collector_doping_m3{2.0e21};
    double minority_electron_lifetime_base_s{2.0e-7};
    double minority_hole_lifetime_emitter_s{5.0e-8};
    double minority_hole_lifetime_collector_s{2.0e-7};
    SemiconductorMaterial material{};
};

struct BipolarDerivedParameters {
    double thermal_voltage_v{};
    double electron_diffusion_length_base_m{};
    double hole_diffusion_length_emitter_m{};
    double hole_diffusion_length_collector_m{};
    double base_transport_factor{};
    double emitter_injection_efficiency{};
    double collector_injection_efficiency{};
    double forward_alpha{};
    double reverse_alpha{};
    double forward_beta{};
    double reverse_beta{};
    double emitter_saturation_current_a{};
    double collector_saturation_current_a{};
    double transport_saturation_current_a{};
    double emitter_base_built_in_v{};
    double collector_base_built_in_v{};
};

struct BipolarOperatingPoint {
    double collector_voltage_v{};
    double base_voltage_v{};
    double emitter_voltage_v{};
    double collector_current_a{};
    double base_current_a{};
    double emitter_current_a{};
    double dissipated_power_w{};
    BipolarDerivedParameters parameters{};
};

// Physics-parameterized 1-D BJT charge-control baseline. The spatial helper
// exposes an abrupt emitter/base/collector profile while terminal transport is
// represented by an Ebers-Moll relation whose alpha and saturation currents
// are derived from minority-carrier diffusion in the three regions.
class BipolarJunctionTransistor1D {
public:
    explicit BipolarJunctionTransistor1D(BipolarJunctionConfig config = {});

    [[nodiscard]] double total_length_m() const noexcept;
    [[nodiscard]] double net_doping_m3(double x_m) const;
    [[nodiscard]] std::vector<double> doping_profile(std::size_t samples) const;
    [[nodiscard]] BipolarDerivedParameters derived(double temperature_k = 0.0) const;
    [[nodiscard]] BipolarOperatingPoint operating_point(double collector_v,
                                                        double base_v,
                                                        double emitter_v,
                                                        double temperature_k = 0.0) const;
    [[nodiscard]] const BipolarJunctionConfig& config() const noexcept { return config_; }

private:
    BipolarJunctionConfig config_{};
};

// Three-terminal [C,B,E] evaluator generated from BipolarJunctionTransistor1D.
[[nodiscard]] cfd::circuit::StaticDeviceEvaluator
make_bipolar_tcad_spice_evaluator(BipolarJunctionConfig config,
                                  double jacobian_step_v = 1.0e-6);

struct IgbtConfig {
    LongChannelMosfetConfig channel{};
    // Low-gain PNP path representing conductivity modulation of the drift
    // region. The PNP's forward beta is derived from semiconductor geometry.
    BipolarJunctionConfig bipolar{};
    double drift_resistance_ohm{0.08};
    double junction_drop_v{0.65};
    double conductivity_modulation_current_a{0.05};
    double maximum_bipolar_gain{12.0};
    double tail_lifetime_s{2.0e-6};
    double ambient_temperature_k{300.0};
    double thermal_resistance_k_per_w{0.6};
    double thermal_capacitance_j_per_k{0.02};
    double drift_resistance_temperature_coefficient_per_k{3.5e-3};
};

struct IgbtOperatingPoint {
    double gate_emitter_voltage_v{};
    double collector_emitter_voltage_v{};
    double collector_current_a{};
    double channel_current_a{};
    double bipolar_current_a{};
    double effective_drift_resistance_ohm{};
    double junction_temperature_k{};
    double dissipated_power_w{};
};

struct IgbtTransientState {
    double junction_temperature_k{300.0};
    double tail_current_a{};
    double time_s{};
};

// MOS-gated bipolar power-device foundation. It composes the TCAD-derived
// long-channel MOS channel and BJT transport baseline with drift-region
// resistance/conductivity modulation, turn-off tail storage and a lumped
// thermal RC. It is intentionally not a full Hefner/field-resolved IGBT.
class IgbtPowerDevice {
public:
    explicit IgbtPowerDevice(IgbtConfig config = {});

    [[nodiscard]] IgbtOperatingPoint operating_point(double vge_v,
                                                     double vce_v,
                                                     double temperature_k = 0.0) const;
    [[nodiscard]] IgbtOperatingPoint step(double vge_v, double vce_v, double dt_s);
    void reset(double temperature_k = 0.0);
    [[nodiscard]] const IgbtTransientState& state() const noexcept { return state_; }
    [[nodiscard]] const IgbtConfig& config() const noexcept { return config_; }

private:
    [[nodiscard]] double static_current(double vge_v, double vce_v, double temperature_k,
                                        double& channel_current_a,
                                        double& bipolar_current_a,
                                        double& effective_drift_resistance_ohm) const;

    IgbtConfig config_{};
    IgbtTransientState state_{};
};

// Three-terminal [C,G,E] static evaluator for the MOS-gated bipolar baseline.
[[nodiscard]] cfd::circuit::StaticDeviceEvaluator
make_igbt_tcad_spice_evaluator(IgbtConfig config,
                               double jacobian_step_v = 1.0e-6);

} // namespace cfd::tcad
