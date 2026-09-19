#include "cfd/tcad/advanced_devices.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::tcad {
namespace {
constexpr double boltzmann_over_electron_charge = 8.617333262145e-5;

void require_finite_positive(double value, const char* message) {
    if (!(value > 0.0) || !std::isfinite(value))
        throw std::invalid_argument(message);
}

[[nodiscard]] double diode_current(double saturation, double voltage, double ideality, double temperature) {
    const double vt = boltzmann_over_electron_charge * temperature * ideality;
    return saturation * std::expm1(std::clamp(voltage / vt, -80.0, 80.0));
}

} // namespace

GummelPoonBjt::GummelPoonBjt(GummelPoonConfig config) : config_(config) {
    require_finite_positive(config_.saturation_current_a, "invalid Gummel-Poon saturation current");
    require_finite_positive(config_.forward_beta, "invalid Gummel-Poon forward beta");
    require_finite_positive(config_.reverse_beta, "invalid Gummel-Poon reverse beta");
    require_finite_positive(config_.early_voltage_forward_v, "invalid Gummel-Poon Early voltage");
    require_finite_positive(config_.early_voltage_reverse_v, "invalid reverse Early voltage");
    if (!(config_.high_injection_current_a > 0.0) || !(config_.temperature_k > 0.0) ||
        !std::isfinite(config_.high_injection_current_a) || !std::isfinite(config_.temperature_k))
        throw std::invalid_argument("invalid Gummel-Poon controls");
}

GummelPoonOperatingPoint GummelPoonBjt::operating_point(double vbe_v, double vbc_v) const {
    if (!std::isfinite(vbe_v) || !std::isfinite(vbc_v))
        throw std::invalid_argument("invalid BJT terminal voltage");
    const double vt = boltzmann_over_electron_charge * config_.temperature_k;
    const double forward_raw = config_.saturation_current_a * std::expm1(std::clamp(vbe_v / vt, -80.0, 80.0));
    const double reverse_raw = config_.saturation_current_a * std::expm1(std::clamp(vbc_v / vt, -80.0, 80.0));
    const double forward = forward_raw / (1.0 + std::abs(forward_raw) / config_.high_injection_current_a);
    const double reverse = reverse_raw / (1.0 + std::abs(reverse_raw) / config_.high_injection_current_a);
    const double alpha_forward = config_.forward_beta / (config_.forward_beta + 1.0);
    const double alpha_reverse = config_.reverse_beta / (config_.reverse_beta + 1.0);
    const double early_forward = 1.0 + (vbe_v - vbc_v) / config_.early_voltage_forward_v;
    const double early_reverse = 1.0 + (vbc_v - vbe_v) / config_.early_voltage_reverse_v;
    const double forward_transport = alpha_forward * forward * std::max(early_forward, 0.05);
    const double reverse_transport = alpha_reverse * reverse * std::max(early_reverse, 0.05);
    const double collector = forward_transport - reverse;
    const double base = (1.0 - alpha_forward) * forward + (1.0 - alpha_reverse) * reverse;
    const double emitter = -collector - base;
    return {collector, base, emitter, forward_transport, reverse_transport};
}

MosLevel23::MosLevel23(MosLevel23Config config) : config_(config) {
    require_finite_positive(config_.oxide_capacitance_f_per_m2, "invalid MOS oxide capacitance");
    require_finite_positive(config_.width_m, "invalid MOS width");
    require_finite_positive(config_.length_m, "invalid MOS length");
    require_finite_positive(config_.mobility_m2_per_vs, "invalid MOS mobility");
    if (!(config_.subthreshold_factor > 1.0) || !(config_.subthreshold_current_a > 0.0) ||
        !(config_.temperature_k > 0.0) || !std::isfinite(config_.threshold_voltage_v) ||
        !std::isfinite(config_.channel_length_modulation_per_v) || !std::isfinite(config_.mobility_degradation_per_v)) {
        throw std::invalid_argument("invalid MOS level 2/3 controls");
    }
}

MosLevel23OperatingPoint MosLevel23::operating_point(double vgs_v, double vds_v) const {
    if (!std::isfinite(vgs_v) || !std::isfinite(vds_v))
        throw std::invalid_argument("invalid MOS terminal voltage");
    const double vt = boltzmann_over_electron_charge * config_.temperature_k;
    const double overdrive = vgs_v - config_.threshold_voltage_v;
    if (overdrive <= 0.0) {
        const double exponent = std::clamp(overdrive / (config_.subthreshold_factor * vt), -80.0, 80.0);
        const double current = config_.subthreshold_current_a * std::exp(exponent) *
                               (1.0 - std::exp(std::clamp(-vds_v / vt, -80.0, 80.0)));
        return {current, current / (config_.subthreshold_factor * vt),
                config_.subthreshold_current_a * std::exp(exponent) / vt, true};
    }
    const double mobility = config_.mobility_m2_per_vs / (1.0 + config_.mobility_degradation_per_v * overdrive);
    const double beta = mobility * config_.oxide_capacitance_f_per_m2 * config_.width_m / config_.length_m;
    const double vds_abs = std::max(vds_v, 0.0);
    const double saturation = std::min(vds_abs, overdrive);
    const double current = beta * (overdrive * saturation - 0.5 * saturation * saturation) *
                           (1.0 + config_.channel_length_modulation_per_v * vds_abs);
    const double gm = beta * std::max(overdrive - saturation, 0.0) +
                      beta * saturation * (1.0 + config_.channel_length_modulation_per_v * vds_abs);
    const double go = current * config_.channel_length_modulation_per_v /
                      std::max(1.0 + config_.channel_length_modulation_per_v * vds_abs, 1.0e-30);
    return {current, gm, go, false};
}

VdmosPowerDevice::VdmosPowerDevice(VdmosConfig config) : config_(config), channel_(config.channel) {
    require_finite_positive(config_.drift_resistance_ohm, "invalid VDMOS drift resistance");
    require_finite_positive(config_.body_diode_saturation_a, "invalid VDMOS body diode");
    require_finite_positive(config_.body_diode_ideality, "invalid VDMOS diode ideality");
}

VdmosOperatingPoint VdmosPowerDevice::operating_point(double vgs_v, double vds_v) const {
    if (!std::isfinite(vgs_v) || !std::isfinite(vds_v))
        throw std::invalid_argument("invalid VDMOS terminal voltage");
    if (vds_v >= 0.0) {
        const auto channel = channel_.operating_point(vgs_v, vds_v);
        const double drop = channel.drain_current_a * config_.drift_resistance_ohm;
        return {channel.drain_current_a, channel.drain_current_a, 0.0, drop};
    }
    const double diode = diode_current(config_.body_diode_saturation_a, -vds_v, config_.body_diode_ideality,
                                       config_.channel.temperature_k);
    return {-diode, 0.0, -diode, 0.0};
}

HystereticCoreModel::HystereticCoreModel(HystereticCoreConfig config) : config_(config) {
    require_finite_positive(config_.saturation_induction_t, "invalid core saturation induction");
    require_finite_positive(config_.coercive_field_a_per_m, "invalid core coercive field");
    require_finite_positive(config_.relaxation_time_s, "invalid core relaxation time");
    require_finite_positive(config_.reference_temperature_k, "invalid core reference temperature");
}

void HystereticCoreModel::step(double field_a_per_m, double temperature_k, double dt_s) {
    if (!std::isfinite(field_a_per_m) || !(temperature_k > 0.0) || !(dt_s > 0.0) || !std::isfinite(temperature_k) ||
        !std::isfinite(dt_s)) {
        throw std::invalid_argument("invalid hysteretic core state");
    }
    const double temperature_factor =
        std::max(0.05, 1.0 - config_.temperature_coefficient_per_k * (temperature_k - config_.reference_temperature_k));
    const double target = config_.saturation_induction_t * temperature_factor *
                          std::tanh((field_a_per_m - std::copysign(config_.coercive_field_a_per_m,
                                                                   field_a_per_m - state_.previous_field_a_per_m)) /
                                    config_.coercive_field_a_per_m);
    const double blend = 1.0 - std::exp(-dt_s / config_.relaxation_time_s);
    const double old = state_.induction_t;
    state_.induction_t += blend * (target - state_.induction_t);
    state_.loss_density_j_per_m3 += std::abs((state_.induction_t - old) * field_a_per_m);
    state_.previous_field_a_per_m = field_a_per_m;
}

} // namespace cfd::tcad
