#include "cfd/tcad/bipolar_power.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::tcad {
namespace {

[[nodiscard]] double safe_expm1(double x) {
    return std::expm1(std::clamp(x, -80.0, 80.0));
}

[[nodiscard]] double positive_coth(double x) {
    if (!(x > 0.0)) throw std::invalid_argument("coth argument must be positive");
    if (x < 1.0e-4) {
        return 1.0 / x + x / 3.0;
    }
    if (x > 30.0) return 1.0;
    return 1.0 / std::tanh(x);
}

[[nodiscard]] double bounded_cosh_inverse(double x) {
    if (x > 40.0) return 0.0;
    return 1.0 / std::cosh(x);
}

[[nodiscard]] SemiconductorMaterial material_at_temperature(const SemiconductorMaterial& reference,
                                                            double temperature_k) {
    if (!(temperature_k > 0.0) || !std::isfinite(temperature_k)) {
        throw std::invalid_argument("semiconductor temperature must be positive and finite");
    }
    SemiconductorMaterial m = reference;
    if (std::abs(temperature_k - reference.temperature_k) < 1.0e-12) return m;
    const double eg_ref = silicon_bandgap_ev(reference.temperature_k);
    const double eg = silicon_bandgap_ev(temperature_k);
    m.intrinsic_density_m3 = temperature_scaled_intrinsic_density(
        reference.intrinsic_density_m3, reference.temperature_k, temperature_k, eg_ref, eg);
    m.electron_mobility_m2_per_vs = temperature_scaled_mobility(
        reference.electron_mobility_m2_per_vs, reference.temperature_k, temperature_k, 2.42);
    m.hole_mobility_m2_per_vs = temperature_scaled_mobility(
        reference.hole_mobility_m2_per_vs, reference.temperature_k, temperature_k, 2.20);
    m.conduction_band_density_m3 = reference.conduction_band_density_m3
        * std::pow(temperature_k / reference.temperature_k, 1.5);
    m.valence_band_density_m3 = reference.valence_band_density_m3
        * std::pow(temperature_k / reference.temperature_k, 1.5);
    m.temperature_k = temperature_k;
    return m;
}

[[nodiscard]] double built_in_voltage(double majority_a_m3, double majority_b_m3,
                                      double ni_m3, double vt_v) {
    const double ratio = majority_a_m3 * majority_b_m3 / (ni_m3 * ni_m3);
    return vt_v * std::log(std::max(ratio, 1.0));
}

[[nodiscard]] double finite_difference_step(double x, double base_step) {
    return std::max(base_step, std::abs(x) * 1.0e-7);
}

} // namespace

BipolarJunctionTransistor1D::BipolarJunctionTransistor1D(BipolarJunctionConfig config)
    : config_(config) {
    if (!(config_.area_m2 > 0.0)
        || !(config_.emitter_width_m > 0.0)
        || !(config_.base_width_m > 0.0)
        || !(config_.collector_width_m > 0.0)
        || !(config_.emitter_doping_m3 > 0.0)
        || !(config_.base_doping_m3 > 0.0)
        || !(config_.collector_doping_m3 > 0.0)
        || !(config_.minority_electron_lifetime_base_s > 0.0)
        || !(config_.minority_hole_lifetime_emitter_s > 0.0)
        || !(config_.minority_hole_lifetime_collector_s > 0.0)
        || !(config_.material.temperature_k > 0.0)
        || !(config_.material.intrinsic_density_m3 > 0.0)
        || !(config_.material.electron_mobility_m2_per_vs > 0.0)
        || !(config_.material.hole_mobility_m2_per_vs > 0.0)) {
        throw std::invalid_argument("invalid bipolar-junction configuration");
    }
}

double BipolarJunctionTransistor1D::total_length_m() const noexcept {
    return config_.emitter_width_m + config_.base_width_m + config_.collector_width_m;
}

double BipolarJunctionTransistor1D::net_doping_m3(double x_m) const {
    const double length = total_length_m();
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() * std::max(length, 1.0);
    if (!std::isfinite(x_m) || x_m < -tolerance || x_m > length + tolerance) {
        throw std::out_of_range("BJT profile coordinate outside device");
    }
    x_m = std::clamp(x_m, 0.0, length);
    const bool npn = config_.polarity == BipolarPolarity::npn;
    if (x_m < config_.emitter_width_m) {
        return (npn ? 1.0 : -1.0) * config_.emitter_doping_m3;
    }
    if (x_m < config_.emitter_width_m + config_.base_width_m) {
        return (npn ? -1.0 : 1.0) * config_.base_doping_m3;
    }
    return (npn ? 1.0 : -1.0) * config_.collector_doping_m3;
}

std::vector<double> BipolarJunctionTransistor1D::doping_profile(std::size_t samples) const {
    if (samples < 3U) throw std::invalid_argument("BJT doping profile requires at least three samples");
    std::vector<double> result(samples, 0.0);
    const double length = total_length_m();
    for (std::size_t i = 0U; i < samples; ++i) {
        const double x = length * static_cast<double>(i) / static_cast<double>(samples - 1U);
        result[i] = net_doping_m3(x);
    }
    return result;
}

BipolarDerivedParameters BipolarJunctionTransistor1D::derived(double temperature_k) const {
    if (temperature_k == 0.0) temperature_k = config_.material.temperature_k;
    const auto material = material_at_temperature(config_.material, temperature_k);
    const double vt = thermal_voltage(temperature_k);
    const double dn_base = material.electron_mobility_m2_per_vs * vt;
    const double dp_emitter = material.hole_mobility_m2_per_vs * vt;
    const double dp_collector = material.hole_mobility_m2_per_vs * vt;
    const double ln_base = std::sqrt(dn_base * config_.minority_electron_lifetime_base_s);
    const double lp_emitter = std::sqrt(dp_emitter * config_.minority_hole_lifetime_emitter_s);
    const double lp_collector = std::sqrt(dp_collector * config_.minority_hole_lifetime_collector_s);
    const double wb_over_ln = config_.base_width_m / ln_base;
    const double alpha_t = bounded_cosh_inverse(wb_over_ln);
    const double ni2 = material.intrinsic_density_m3 * material.intrinsic_density_m3;

    const double jn_base = elementary_charge_c * dn_base * ni2
        / (config_.base_doping_m3 * ln_base) * positive_coth(wb_over_ln);
    const double jp_emitter = elementary_charge_c * dp_emitter * ni2
        / (config_.emitter_doping_m3 * lp_emitter)
        * positive_coth(config_.emitter_width_m / lp_emitter);
    const double jp_collector = elementary_charge_c * dp_collector * ni2
        / (config_.collector_doping_m3 * lp_collector)
        * positive_coth(config_.collector_width_m / lp_collector);

    const double gamma_e = jn_base / (jn_base + jp_emitter);
    const double gamma_c = jn_base / (jn_base + jp_collector);
    const double alpha_f = std::clamp(alpha_t * gamma_e, 1.0e-12, 1.0 - 1.0e-12);
    const double alpha_r = std::clamp(alpha_t * gamma_c, 1.0e-12, 1.0 - 1.0e-12);
    const double ies = config_.area_m2 * (jn_base + jp_emitter);
    const double ics = config_.area_m2 * (jn_base + jp_collector);
    // Reciprocity follows from alpha_F*I_ES = alpha_R*I_CS = A*alpha_T*J_n.
    const double transport = config_.area_m2 * alpha_t * jn_base;

    BipolarDerivedParameters p;
    p.thermal_voltage_v = vt;
    p.electron_diffusion_length_base_m = ln_base;
    p.hole_diffusion_length_emitter_m = lp_emitter;
    p.hole_diffusion_length_collector_m = lp_collector;
    p.base_transport_factor = alpha_t;
    p.emitter_injection_efficiency = gamma_e;
    p.collector_injection_efficiency = gamma_c;
    p.forward_alpha = alpha_f;
    p.reverse_alpha = alpha_r;
    p.forward_beta = alpha_f / (1.0 - alpha_f);
    p.reverse_beta = alpha_r / (1.0 - alpha_r);
    p.emitter_saturation_current_a = ies;
    p.collector_saturation_current_a = ics;
    p.transport_saturation_current_a = transport;
    p.emitter_base_built_in_v = built_in_voltage(config_.emitter_doping_m3,
                                                  config_.base_doping_m3,
                                                  material.intrinsic_density_m3, vt);
    p.collector_base_built_in_v = built_in_voltage(config_.collector_doping_m3,
                                                    config_.base_doping_m3,
                                                    material.intrinsic_density_m3, vt);
    return p;
}

BipolarOperatingPoint BipolarJunctionTransistor1D::operating_point(double vc, double vb, double ve,
                                                                    double temperature_k) const {
    if (!std::isfinite(vc) || !std::isfinite(vb) || !std::isfinite(ve)) {
        throw std::invalid_argument("BJT terminal voltages must be finite");
    }
    const auto p = derived(temperature_k);
    const double polarity = config_.polarity == BipolarPolarity::npn ? 1.0 : -1.0;
    const double vbe = polarity * (vb - ve);
    const double vbc = polarity * (vb - vc);
    const double forward = p.emitter_saturation_current_a * safe_expm1(vbe / p.thermal_voltage_v);
    const double reverse = p.collector_saturation_current_a * safe_expm1(vbc / p.thermal_voltage_v);
    const double ic_eff = p.forward_alpha * forward - reverse;
    const double ie_eff = -forward + p.reverse_alpha * reverse;
    const double ib_eff = -(ic_eff + ie_eff);
    const double ic = polarity * ic_eff;
    const double ib = polarity * ib_eff;
    const double ie = polarity * ie_eff;
    const double power = vc * ic + vb * ib + ve * ie;
    return {vc, vb, ve, ic, ib, ie, power, p};
}

cfd::circuit::StaticDeviceEvaluator
make_bipolar_tcad_spice_evaluator(BipolarJunctionConfig config, double h) {
    if (!(h > 0.0) || !std::isfinite(h)) throw std::invalid_argument("BJT Jacobian step must be positive");
    const BipolarJunctionTransistor1D model(config);
    return [model, h](std::span<const double> v) -> cfd::circuit::StaticDeviceEvaluation {
        if (v.size() != 3U) throw std::invalid_argument("TCAD BJT evaluator expects [C,B,E]");
        auto currents = [&model](std::span<const double> x) {
            const auto op = model.operating_point(x[0U], x[1U], x[2U]);
            return std::vector<double>{op.collector_current_a, op.base_current_a, op.emitter_current_a};
        };
        const auto current = currents(v);
        std::vector<double> jac(9U, 0.0);
        std::vector<double> vp(v.begin(), v.end());
        std::vector<double> vm(v.begin(), v.end());
        for (std::size_t col = 0U; col < 3U; ++col) {
            vp.assign(v.begin(), v.end());
            vm.assign(v.begin(), v.end());
            const double step = finite_difference_step(v[col], h);
            vp[col] += step;
            vm[col] -= step;
            const auto ip = currents(vp);
            const auto im = currents(vm);
            for (std::size_t row = 0U; row < 3U; ++row) {
                jac[row * 3U + col] = (ip[row] - im[row]) / (2.0 * step);
            }
        }
        return {current, jac};
    };
}

IgbtPowerDevice::IgbtPowerDevice(IgbtConfig config) : config_(config) {
    if (!(config_.drift_resistance_ohm > 0.0)
        || !(config_.junction_drop_v >= 0.0)
        || !(config_.conductivity_modulation_current_a > 0.0)
        || !(config_.maximum_bipolar_gain >= 0.0)
        || !(config_.tail_lifetime_s > 0.0)
        || !(config_.ambient_temperature_k > 0.0)
        || !(config_.thermal_resistance_k_per_w > 0.0)
        || !(config_.thermal_capacitance_j_per_k > 0.0)
        || !(config_.drift_resistance_temperature_coefficient_per_k >= 0.0)) {
        throw std::invalid_argument("invalid IGBT configuration");
    }
    config_.bipolar.polarity = BipolarPolarity::pnp;
    state_.junction_temperature_k = config_.ambient_temperature_k;
}

double IgbtPowerDevice::static_current(double vge, double vce, double temperature_k,
                                       double& channel_current, double& bipolar_current,
                                       double& effective_drift_resistance) const {
    channel_current = 0.0;
    bipolar_current = 0.0;
    if (!(vce > config_.junction_drop_v)) {
        effective_drift_resistance = config_.drift_resistance_ohm;
        return 0.0;
    }
    const auto bjt = BipolarJunctionTransistor1D(config_.bipolar).derived(temperature_k);
    const double beta = std::min(config_.maximum_bipolar_gain, bjt.forward_beta);
    const double mobility_scale = temperature_scaled_mobility(1.0, config_.ambient_temperature_k,
                                                               temperature_k, 2.42);
    LongChannelMosfetConfig channel_cfg = config_.channel;
    channel_cfg.electron_mobility_m2_per_vs *= mobility_scale;
    const LongChannelMosfet mos(channel_cfg);
    const double r_temp = config_.drift_resistance_ohm
        * (1.0 + config_.drift_resistance_temperature_coefficient_per_k
                 * (temperature_k - config_.ambient_temperature_k));

    double total = 0.0;
    for (std::size_t iteration = 0U; iteration < 100U; ++iteration) {
        const double modulation = 1.0 + std::abs(channel_current) / config_.conductivity_modulation_current_a;
        effective_drift_resistance = std::max(r_temp / modulation, 1.0e-6);
        const double channel_vds = std::max(vce - config_.junction_drop_v - total * effective_drift_resistance, 0.0);
        const double new_channel = mos.drain_current(vge, channel_vds, 0.0);
        const double new_bipolar = beta * new_channel;
        const double target = new_channel + new_bipolar;
        const double next = 0.5 * total + 0.5 * target;
        channel_current = new_channel;
        bipolar_current = new_bipolar;
        if (std::abs(next - total) <= 1.0e-10 * std::max(1.0, std::abs(next))) {
            total = next;
            break;
        }
        total = next;
    }
    return std::max(total, 0.0);
}

IgbtOperatingPoint IgbtPowerDevice::operating_point(double vge, double vce, double temperature_k) const {
    if (!std::isfinite(vge) || !std::isfinite(vce)) throw std::invalid_argument("IGBT voltages must be finite");
    if (temperature_k == 0.0) temperature_k = config_.ambient_temperature_k;
    double channel = 0.0;
    double bipolar = 0.0;
    double resistance = config_.drift_resistance_ohm;
    const double current = static_current(vge, vce, temperature_k, channel, bipolar, resistance);
    return {vge, vce, current, channel, bipolar, resistance, temperature_k, current * std::max(vce, 0.0)};
}

IgbtOperatingPoint IgbtPowerDevice::step(double vge, double vce, double dt_s) {
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) throw std::invalid_argument("IGBT timestep must be positive");
    auto static_op = operating_point(vge, vce, state_.junction_temperature_k);
    const double target_tail = static_op.bipolar_current_a;
    const double decay = std::exp(-dt_s / config_.tail_lifetime_s);
    if (static_op.channel_current_a > 0.0) {
        state_.tail_current_a = target_tail + (state_.tail_current_a - target_tail) * decay;
    } else {
        state_.tail_current_a *= decay;
    }
    const double conduction = static_op.channel_current_a + state_.tail_current_a;
    const double power = conduction * std::max(vce, 0.0);
    const double cooling = (state_.junction_temperature_k - config_.ambient_temperature_k)
                         / config_.thermal_resistance_k_per_w;
    state_.junction_temperature_k += dt_s
        * (power - cooling) / config_.thermal_capacitance_j_per_k;
    state_.junction_temperature_k = std::max(state_.junction_temperature_k, 1.0);
    state_.time_s += dt_s;
    static_op.collector_current_a = conduction;
    static_op.bipolar_current_a = state_.tail_current_a;
    static_op.junction_temperature_k = state_.junction_temperature_k;
    static_op.dissipated_power_w = power;
    return static_op;
}

void IgbtPowerDevice::reset(double temperature_k) {
    if (temperature_k == 0.0) temperature_k = config_.ambient_temperature_k;
    if (!(temperature_k > 0.0) || !std::isfinite(temperature_k)) {
        throw std::invalid_argument("IGBT reset temperature must be positive");
    }
    state_ = {};
    state_.junction_temperature_k = temperature_k;
}

cfd::circuit::StaticDeviceEvaluator
make_igbt_tcad_spice_evaluator(IgbtConfig config, double h) {
    if (!(h > 0.0) || !std::isfinite(h)) throw std::invalid_argument("IGBT Jacobian step must be positive");
    const IgbtPowerDevice model(config);
    return [model, h](std::span<const double> v) -> cfd::circuit::StaticDeviceEvaluation {
        if (v.size() != 3U) throw std::invalid_argument("TCAD IGBT evaluator expects [C,G,E]");
        auto currents = [&model](std::span<const double> x) {
            const auto op = model.operating_point(x[1U] - x[2U], x[0U] - x[2U]);
            return std::vector<double>{op.collector_current_a, 0.0, -op.collector_current_a};
        };
        const auto current = currents(v);
        std::vector<double> jac(9U, 0.0);
        std::vector<double> vp(v.begin(), v.end());
        std::vector<double> vm(v.begin(), v.end());
        for (std::size_t col = 0U; col < 3U; ++col) {
            vp.assign(v.begin(), v.end());
            vm.assign(v.begin(), v.end());
            const double step = finite_difference_step(v[col], h);
            vp[col] += step;
            vm[col] -= step;
            const auto ip = currents(vp);
            const auto im = currents(vm);
            for (std::size_t row = 0U; row < 3U; ++row) {
                jac[row * 3U + col] = (ip[row] - im[row]) / (2.0 * step);
            }
        }
        return {current, jac};
    };
}

} // namespace cfd::tcad
