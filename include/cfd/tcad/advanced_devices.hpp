#pragma once

#include <complex>

namespace cfd::tcad {

struct GummelPoonConfig {
    double saturation_current_a{1.0e-15};
    double forward_beta{120.0};
    double reverse_beta{2.0};
    double early_voltage_forward_v{80.0};
    double early_voltage_reverse_v{20.0};
    double high_injection_current_a{0.02};
    double temperature_k{300.0};
};

struct GummelPoonOperatingPoint {
    double collector_current_a{};
    double base_current_a{};
    double emitter_current_a{};
    double forward_transport_current_a{};
    double reverse_transport_current_a{};
};

// Charge-control Gummel-Poon compact baseline with Early effect and smooth
// high-injection roll-off. It is suitable for circuit regressions, not a
// field-resolved TCAD replacement.
class GummelPoonBjt {
public:
    explicit GummelPoonBjt(GummelPoonConfig config = {});
    [[nodiscard]] GummelPoonOperatingPoint operating_point(double vbe_v, double vbc_v) const;

private:
    GummelPoonConfig config_;
};

struct MosLevel23Config {
    double oxide_capacitance_f_per_m2{3.45e-3};
    double width_m{10.0e-6};
    double length_m{1.0e-6};
    double mobility_m2_per_vs{0.05};
    double threshold_voltage_v{0.5};
    double channel_length_modulation_per_v{0.02};
    double mobility_degradation_per_v{0.08};
    double subthreshold_factor{1.5};
    double subthreshold_current_a{1.0e-12};
    double temperature_k{300.0};
};

struct MosLevel23OperatingPoint {
    double drain_current_a{};
    double transconductance_s{};
    double output_conductance_s{};
    bool subthreshold{};
};

class MosLevel23 {
public:
    explicit MosLevel23(MosLevel23Config config = {});
    [[nodiscard]] MosLevel23OperatingPoint operating_point(double vgs_v, double vds_v) const;

private:
    MosLevel23Config config_;
};

struct VdmosConfig {
    MosLevel23Config channel{};
    double drift_resistance_ohm{0.08};
    double body_diode_saturation_a{1.0e-12};
    double body_diode_ideality{1.4};
};

struct VdmosOperatingPoint {
    double drain_current_a{};
    double channel_current_a{};
    double body_diode_current_a{};
    double drift_drop_v{};
};

class VdmosPowerDevice {
public:
    explicit VdmosPowerDevice(VdmosConfig config = {});
    [[nodiscard]] VdmosOperatingPoint operating_point(double vgs_v, double vds_v) const;

private:
    VdmosConfig config_;
    MosLevel23 channel_;
};

struct HystereticCoreConfig {
    double saturation_induction_t{1.5};
    double coercive_field_a_per_m{40.0};
    double relaxation_time_s{1.0e-4};
    double temperature_coefficient_per_k{0.0};
    double reference_temperature_k{300.0};
};

struct HystereticCoreState {
    double induction_t{};
    double loss_density_j_per_m3{};
    double previous_field_a_per_m{};
};

// Rate-dependent tanh hysteresis baseline for nonlinear transformer/core
// coupling. The state exposes accumulated magnetic loss for circuit/thermal
// co-simulation.
class HystereticCoreModel {
public:
    explicit HystereticCoreModel(HystereticCoreConfig config = {});
    void reset() noexcept { state_ = {}; }
    void step(double field_a_per_m, double temperature_k, double dt_s);
    [[nodiscard]] const HystereticCoreState& state() const noexcept { return state_; }

private:
    HystereticCoreConfig config_;
    HystereticCoreState state_{};
};

} // namespace cfd::tcad
