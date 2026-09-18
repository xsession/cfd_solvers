#include "cfd/tcad/mos.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::tcad {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

[[nodiscard]] double simpson_fermi_half(double eta) {
    const double upper = std::max(40.0, eta + 40.0);
    constexpr std::size_t intervals = 4096U;
    const double h = upper / static_cast<double>(intervals);
    auto integrand = [eta](double energy) {
        if (energy <= 0.0) return 0.0;
        const double z = energy - eta;
        double occupation{};
        if (z > 50.0) occupation = std::exp(-z);
        else if (z < -50.0) occupation = 1.0;
        else occupation = 1.0 / (1.0 + std::exp(z));
        return std::sqrt(energy) * occupation;
    };
    double sum = integrand(0.0) + integrand(upper);
    for (std::size_t i = 1U; i < intervals; ++i) {
        sum += (i % 2U == 0U ? 2.0 : 4.0) * integrand(static_cast<double>(i) * h);
    }
    return (2.0 / std::sqrt(pi)) * (h / 3.0) * sum;
}

[[nodiscard]] double safe_exp(double x) {
    return std::exp(std::clamp(x, -100.0, 100.0));
}

[[nodiscard]] double central_derivative(const std::function<double(double)>& f,
                                        double x,
                                        double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

} // namespace

double fermi_dirac_half(double eta) {
    if (!std::isfinite(eta)) throw std::invalid_argument("Fermi-Dirac argument must be finite");
    if (eta < -8.0) return std::exp(eta);
    if (eta > 40.0) {
        // Sommerfeld leading term; this branch is mainly to keep extreme
        // degeneracy finite without spending time in quadrature.
        return (4.0 / (3.0 * std::sqrt(pi))) * std::pow(eta, 1.5);
    }
    return simpson_fermi_half(eta);
}

double inverse_fermi_dirac_half(double value) {
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::invalid_argument("inverse Fermi-Dirac value must be positive and finite");
    }
    double lo = -100.0;
    double hi = 100.0;
    for (std::size_t i = 0U; i < 160U; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (fermi_dirac_half(mid) < value) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

double fermi_electron_density(double nc, double eta) {
    if (!(nc > 0.0)) throw std::invalid_argument("conduction-band density must be positive");
    return nc * fermi_dirac_half(eta);
}

double fermi_hole_density(double nv, double eta) {
    if (!(nv > 0.0)) throw std::invalid_argument("valence-band density must be positive");
    return nv * fermi_dirac_half(eta);
}

double ionized_donor_density(double total, double eta_n, double binding_ev,
                             double degeneracy, double temperature_k) {
    if (!(total >= 0.0) || !(binding_ev >= 0.0) || !(degeneracy > 0.0)) {
        throw std::invalid_argument("invalid donor-ionization parameters");
    }
    const double exponent = eta_n + binding_ev / thermal_voltage(temperature_k);
    return total / (1.0 + degeneracy * safe_exp(exponent));
}

double ionized_acceptor_density(double total, double eta_p, double binding_ev,
                                double degeneracy, double temperature_k) {
    if (!(total >= 0.0) || !(binding_ev >= 0.0) || !(degeneracy > 0.0)) {
        throw std::invalid_argument("invalid acceptor-ionization parameters");
    }
    const double exponent = eta_p + binding_ev / thermal_voltage(temperature_k);
    return total / (1.0 + degeneracy * safe_exp(exponent));
}

MosCapacitor1D::MosCapacitor1D(MosCapacitorConfig config) : config_(config) {
    if (!(config_.semiconductor_thickness_m > 0.0) || config_.semiconductor_nodes < 5U
        || !(config_.area_m2 > 0.0) || !(config_.substrate_acceptor_m3 > 0.0)
        || !(config_.oxide_thickness_m > 0.0) || !(config_.oxide_relative_permittivity > 0.0)
        || config_.max_iterations == 0U || !(config_.tolerance > 0.0)
        || !(config_.relaxation > 0.0 && config_.relaxation <= 1.0)) {
        throw std::invalid_argument("invalid MOS capacitor configuration");
    }
}

double MosCapacitor1D::oxide_capacitance_per_area() const {
    return vacuum_permittivity_f_per_m * config_.oxide_relative_permittivity / config_.oxide_thickness_m;
}

SemiconductorDevice1D MosCapacitor1D::make_device(double gate_voltage_v) const {
    Semiconductor1DConfig cfg;
    cfg.length_m = config_.semiconductor_thickness_m;
    cfg.nodes = config_.semiconductor_nodes;
    cfg.area_m2 = config_.area_m2;
    cfg.material = config_.semiconductor;
    cfg.max_gummel_iterations = config_.max_iterations;
    cfg.relative_tolerance = config_.tolerance;
    cfg.under_relaxation = config_.relaxation;

    const double ni = cfg.material.intrinsic_density_m3;
    const double vt = thermal_voltage(cfg.material.temperature_k);
    const double na = config_.substrate_acceptor_m3;
    const double root = std::sqrt(na * na + 4.0 * ni * ni);
    const double n_bulk = 0.5 * (-na + root);
    const double bulk_potential = vt * std::log(std::max(n_bulk, 1.0) / ni);

    cfg.left.type = ContactType::gate;
    cfg.left.voltage_v = gate_voltage_v;
    cfg.left.oxide_capacitance_f_per_m2 = oxide_capacitance_per_area();
    // SemiconductorDevice1D stores absolute electrostatic potential. Shift the
    // externally specified flatband voltage so the Robin condition is expressed
    // in surface potential relative to the neutral bulk.
    cfg.left.flatband_voltage_v = config_.flatband_voltage_v - bulk_potential;
    cfg.right.type = ContactType::ohmic;
    cfg.right.voltage_v = 0.0;

    SemiconductorDevice1D device(cfg);
    device.set_net_doping([na](double) { return -na; });
    device.initialize_charge_neutral();
    return device;
}

MosCapacitorPoint MosCapacitor1D::solve(double gate_voltage_v, double dv) const {
    if (!std::isfinite(gate_voltage_v) || !(dv > 0.0) || !std::isfinite(dv)) {
        throw std::invalid_argument("invalid MOS capacitor bias/perturbation");
    }
    const auto solve_one = [this](double vg) {
        const double ni = config_.semiconductor.intrinsic_density_m3;
        const double na = config_.substrate_acceptor_m3;
        const double vt = thermal_voltage(config_.semiconductor.temperature_k);
        const double eps = vacuum_permittivity_f_per_m * config_.semiconductor.relative_permittivity;
        const double cox = oxide_capacitance_per_area();
        const double root = std::sqrt(na * na + 4.0 * ni * ni);
        const double p0 = 0.5 * (na + root);
        const double n0 = ni * ni / p0;
        auto semiconductor_charge = [&](double psi) {
            const double u = psi / vt;
            const double ep = safe_exp(-u);
            const double en = safe_exp(u);
            const double energy_density = elementary_charge_c * vt
                * (p0 * (ep + u - 1.0) + n0 * (en - u - 1.0));
            const double magnitude = std::sqrt(std::max(0.0, 2.0 * eps * energy_density));
            if (std::abs(psi) < 1.0e-18) return 0.0;
            return psi > 0.0 ? -magnitude : magnitude;
        };
        auto residual = [&](double psi) {
            return config_.flatband_voltage_v + psi - semiconductor_charge(psi) / cox - vg;
        };
        double lo = -5.0;
        double hi = 5.0;
        double flo = residual(lo);
        double fhi = residual(hi);
        if (!(flo <= 0.0 && fhi >= 0.0)) {
            throw std::runtime_error("MOS capacitor surface-potential bracket failed");
        }
        for (std::size_t i = 0U; i < 180U; ++i) {
            const double mid = 0.5 * (lo + hi);
            const double fm = residual(mid);
            if (fm > 0.0) hi = mid;
            else lo = mid;
        }
        const double psi = 0.5 * (lo + hi);
        const double qs = semiconductor_charge(psi);
        const double qg = -qs * config_.area_m2;
        return std::pair<double, double>{psi, qg};
    };
    const auto [psi, q] = solve_one(gate_voltage_v);
    const auto [psi_plus, q_plus] = solve_one(gate_voltage_v + dv);
    const auto [psi_minus, q_minus] = solve_one(gate_voltage_v - dv);
    (void)psi_plus;
    (void)psi_minus;
    const double capacitance = (q_plus - q_minus) / (2.0 * dv);
    return {gate_voltage_v, psi, q, capacitance, true};
}

std::vector<MosCapacitorPoint> MosCapacitor1D::sweep_cv(double start_v, double stop_v,
                                                       std::size_t points, double dv) const {
    if (points < 2U || !std::isfinite(start_v) || !std::isfinite(stop_v)) {
        throw std::invalid_argument("invalid MOS capacitor C-V sweep");
    }
    std::vector<MosCapacitorPoint> out;
    out.reserve(points);
    for (std::size_t i = 0U; i < points; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(points - 1U);
        out.push_back(solve(start_v + t * (stop_v - start_v), dv));
    }
    return out;
}

LongChannelMosfet::LongChannelMosfet(LongChannelMosfetConfig config) : config_(config) {
    if (!(config_.width_m > 0.0) || !(config_.length_m > 0.0)
        || !(config_.oxide_capacitance_f_per_m2 > 0.0)
        || !(config_.electron_mobility_m2_per_vs > 0.0)
        || !(config_.channel_length_modulation_per_v >= 0.0)
        || !(config_.body_effect_coefficient_sqrt_v >= 0.0)
        || !(config_.surface_potential_2phi_f_v > 0.0)) {
        throw std::invalid_argument("invalid long-channel MOSFET configuration");
    }
}

double LongChannelMosfet::threshold_voltage(double vbs) const {
    // Vsb=-Vbs for an nMOS body-effect convention.
    const double vsb = std::max(-vbs, 0.0);
    return config_.threshold_voltage_v
         + config_.body_effect_coefficient_sqrt_v
         * (std::sqrt(config_.surface_potential_2phi_f_v + vsb)
            - std::sqrt(config_.surface_potential_2phi_f_v));
}

double LongChannelMosfet::drain_current(double vgs, double vds, double vbs) const {
    if (!std::isfinite(vgs) || !std::isfinite(vds) || !std::isfinite(vbs)) {
        throw std::invalid_argument("MOSFET voltages must be finite");
    }
    if (vds < 0.0) {
        // Swap drain/source while preserving the same gate/body absolute potentials.
        return -drain_current(vgs - vds, -vds, vbs - vds);
    }
    const double vth = threshold_voltage(vbs);
    const double overdrive = vgs - vth;
    if (!(overdrive > 0.0)) return 0.0;
    const double k = config_.electron_mobility_m2_per_vs * config_.oxide_capacitance_f_per_m2
                   * config_.width_m / config_.length_m;
    const double veff = std::min(vds, overdrive);
    double current = k * (overdrive * veff - 0.5 * veff * veff);
    if (vds > overdrive) {
        current *= 1.0 + config_.channel_length_modulation_per_v * (vds - overdrive);
    }
    return current;
}

MosfetOperatingPoint LongChannelMosfet::operating_point(double vgs, double vds, double vbs) const {
    constexpr double h = 1.0e-6;
    const auto id_vgs = [this, vds, vbs](double value) { return drain_current(value, vds, vbs); };
    const auto id_vds = [this, vgs, vbs](double value) { return drain_current(vgs, value, vbs); };
    return {vgs, vds, vbs,
            drain_current(vgs, vds, vbs),
            central_derivative(id_vgs, vgs, h),
            central_derivative(id_vds, vds, h),
            threshold_voltage(vbs)};
}

std::vector<MosfetIvSample> LongChannelMosfet::sweep(double vg0, double vg1, std::size_t ng,
                                                     double vd0, double vd1, std::size_t nd,
                                                     double vbs) const {
    if (ng < 2U || nd < 2U) throw std::invalid_argument("MOSFET sweep requires at least two points per axis");
    std::vector<MosfetIvSample> out;
    out.reserve(ng * nd);
    for (std::size_t ig = 0U; ig < ng; ++ig) {
        const double tg = static_cast<double>(ig) / static_cast<double>(ng - 1U);
        const double vg = vg0 + tg * (vg1 - vg0);
        for (std::size_t id = 0U; id < nd; ++id) {
            const double td = static_cast<double>(id) / static_cast<double>(nd - 1U);
            const double vd = vd0 + td * (vd1 - vd0);
            out.push_back({vg, vd, drain_current(vg, vd, vbs)});
        }
    }
    return out;
}

cfd::circuit::StaticDeviceEvaluator
make_long_channel_mosfet_spice_evaluator(LongChannelMosfetConfig config, double h) {
    if (!(h > 0.0) || !std::isfinite(h)) throw std::invalid_argument("SPICE MOSFET Jacobian step must be positive");
    const LongChannelMosfet model(config);
    return [model, h](std::span<const double> v) -> cfd::circuit::StaticDeviceEvaluation {
        if (v.size() != 4U) throw std::invalid_argument("TCAD MOSFET SPICE evaluator expects [D,G,S,B]");
        auto evaluate = [&model](std::span<const double> x) {
            const double id = model.drain_current(x[1U] - x[2U], x[0U] - x[2U], x[3U] - x[2U]);
            return std::vector<double>{id, 0.0, -id, 0.0};
        };
        const auto current = evaluate(v);
        std::vector<double> jac(16U, 0.0);
        std::vector<double> vp(v.begin(), v.end());
        std::vector<double> vm(v.begin(), v.end());
        for (std::size_t col = 0U; col < 4U; ++col) {
            vp.assign(v.begin(), v.end());
            vm.assign(v.begin(), v.end());
            vp[col] += h;
            vm[col] -= h;
            const auto ip = evaluate(vp);
            const auto im = evaluate(vm);
            for (std::size_t row = 0U; row < 4U; ++row) {
                jac[row * 4U + col] = (ip[row] - im[row]) / (2.0 * h);
            }
        }
        return {current, jac};
    };
}

} // namespace cfd::tcad
