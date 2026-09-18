#include "cfd/tcad/semiconductor1d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::tcad {
namespace {

void validate_material(const SemiconductorMaterial& m) {
    if (!(m.relative_permittivity > 0.0) || !(m.intrinsic_density_m3 > 0.0)
        || !(m.electron_mobility_m2_per_vs > 0.0) || !(m.hole_mobility_m2_per_vs > 0.0)
        || !(m.conduction_band_density_m3 > 0.0) || !(m.valence_band_density_m3 > 0.0)
        || !(m.temperature_k > 0.0)) {
        throw std::invalid_argument("invalid semiconductor material");
    }
}

std::vector<double> solve_tridiagonal(std::vector<double> lower,
                                      std::vector<double> diagonal,
                                      std::vector<double> upper,
                                      std::vector<double> rhs) {
    const std::size_t n = diagonal.size();
    if (n == 0U || lower.size() != n || upper.size() != n || rhs.size() != n) {
        throw std::invalid_argument("invalid tridiagonal system");
    }
    for (std::size_t i = 1U; i < n; ++i) {
        if (!(std::abs(diagonal[i - 1U]) > 1.0e-300)) throw std::runtime_error("singular tridiagonal system");
        const double factor = lower[i] / diagonal[i - 1U];
        diagonal[i] -= factor * upper[i - 1U];
        rhs[i] -= factor * rhs[i - 1U];
    }
    std::vector<double> x(n, 0.0);
    if (!(std::abs(diagonal[n - 1U]) > 1.0e-300)) throw std::runtime_error("singular tridiagonal system");
    x[n - 1U] = rhs[n - 1U] / diagonal[n - 1U];
    for (std::size_t ii = n - 1U; ii-- > 0U;) {
        if (!(std::abs(diagonal[ii]) > 1.0e-300)) throw std::runtime_error("singular tridiagonal system");
        x[ii] = (rhs[ii] - upper[ii] * x[ii + 1U]) / diagonal[ii];
    }
    return x;
}

std::pair<double, double> neutral_carriers(double net_doping_m3, double ni_m3) {
    const double root = std::sqrt(net_doping_m3 * net_doping_m3 + 4.0 * ni_m3 * ni_m3);
    const double n = 0.5 * (net_doping_m3 + root);
    const double p = ni_m3 * ni_m3 / std::max(n, std::numeric_limits<double>::min());
    return {n, p};
}

double safe_exp(double x) {
    return std::exp(std::clamp(x, -100.0, 100.0));
}

} // namespace

double thermal_voltage(double temperature_k) {
    if (!(temperature_k > 0.0)) throw std::invalid_argument("temperature must be positive");
    return boltzmann_j_per_k * temperature_k / elementary_charge_c;
}

double bernoulli_function(double x) {
    if (!std::isfinite(x)) throw std::invalid_argument("Bernoulli argument must be finite");
    const double ax = std::abs(x);
    if (ax < 1.0e-5) {
        const double x2 = x * x;
        return 1.0 - 0.5 * x + x2 / 12.0 - x2 * x2 / 720.0;
    }
    if (x > 50.0) return x * std::exp(-x);
    if (x < -50.0) return -x;
    return x / std::expm1(x);
}


double silicon_bandgap_ev(double temperature_k, double eg0_ev, double alpha_ev_per_k, double beta_k) {
    if (!(temperature_k > 0.0) || !(eg0_ev > 0.0) || !(alpha_ev_per_k >= 0.0) || !(beta_k > 0.0)) {
        throw std::invalid_argument("invalid Varshni parameters");
    }
    return eg0_ev - alpha_ev_per_k * temperature_k * temperature_k / (temperature_k + beta_k);
}

double temperature_scaled_intrinsic_density(double ni_reference_m3,
                                            double reference_temperature_k,
                                            double temperature_k,
                                            double bandgap_reference_ev,
                                            double bandgap_ev) {
    if (!(ni_reference_m3 > 0.0) || !(reference_temperature_k > 0.0) || !(temperature_k > 0.0)
        || !(bandgap_reference_ev > 0.0) || !(bandgap_ev > 0.0)) {
        throw std::invalid_argument("invalid intrinsic-density temperature scaling parameters");
    }
    const double vt_ref = thermal_voltage(reference_temperature_k);
    const double vt = thermal_voltage(temperature_k);
    const double density_scale = std::pow(temperature_k / reference_temperature_k, 1.5);
    const double exponent = bandgap_reference_ev / (2.0 * vt_ref) - bandgap_ev / (2.0 * vt);
    return ni_reference_m3 * density_scale * std::exp(std::clamp(exponent, -100.0, 100.0));
}

double temperature_scaled_mobility(double mobility_reference_m2_per_vs,
                                   double reference_temperature_k,
                                   double temperature_k,
                                   double exponent) {
    if (!(mobility_reference_m2_per_vs > 0.0) || !(reference_temperature_k > 0.0)
        || !(temperature_k > 0.0) || !(exponent >= 0.0)) {
        throw std::invalid_argument("invalid mobility temperature scaling parameters");
    }
    return mobility_reference_m2_per_vs * std::pow(temperature_k / reference_temperature_k, -exponent);
}

double srh_recombination(double n, double p, double ni, const RecombinationModel& m) {
    if (!(n >= 0.0) || !(p >= 0.0) || !(ni > 0.0)) throw std::invalid_argument("invalid SRH state");
    if (!(m.tau_n_s > 0.0) || !(m.tau_p_s > 0.0)) return 0.0;
    const double n1 = m.trap_n1_m3 > 0.0 ? m.trap_n1_m3 : ni;
    const double p1 = m.trap_p1_m3 > 0.0 ? m.trap_p1_m3 : ni;
    const double denominator = m.tau_p_s * (n + n1) + m.tau_n_s * (p + p1);
    const double excess = (n - ni) * p + ni * (p - ni);
    return excess / std::max(denominator, std::numeric_limits<double>::min());
}

double auger_recombination(double n, double p, double ni, const RecombinationModel& m) {
    if (!(n >= 0.0) || !(p >= 0.0) || !(ni > 0.0)) throw std::invalid_argument("invalid Auger state");
    const double excess = (n - ni) * p + ni * (p - ni);
    return (m.auger_n_m6_per_s * n + m.auger_p_m6_per_s * p) * excess;
}

double radiative_recombination(double n, double p, double ni, const RecombinationModel& m) {
    if (!(n >= 0.0) || !(p >= 0.0) || !(ni > 0.0)) throw std::invalid_argument("invalid radiative state");
    const double excess = (n - ni) * p + ni * (p - ni);
    return m.radiative_m3_per_s * excess;
}

double total_recombination(double n, double p, double ni, const RecombinationModel& m) {
    return srh_recombination(n, p, ni, m)
         + auger_recombination(n, p, ni, m)
         + radiative_recombination(n, p, ni, m);
}

double caughey_thomas_mobility(double doping_abs, double mu_min, double mu_max, double reference, double alpha) {
    if (!(doping_abs >= 0.0) || !(mu_min > 0.0) || !(mu_max >= mu_min) || !(reference > 0.0) || !(alpha > 0.0)) {
        throw std::invalid_argument("invalid Caughey-Thomas mobility parameters");
    }
    return mu_min + (mu_max - mu_min) / (1.0 + std::pow(doping_abs / reference, alpha));
}

double high_field_mobility(double mu0, double field, double vsat, double beta) {
    if (!(mu0 > 0.0) || !(vsat > 0.0) || !(beta > 0.0)) throw std::invalid_argument("invalid high-field mobility parameters");
    const double ratio = mu0 * std::abs(field) / vsat;
    return mu0 / std::pow(1.0 + std::pow(ratio, beta), 1.0 / beta);
}

SemiconductorDevice1D::SemiconductorDevice1D(Semiconductor1DConfig config) : config_(config), reference_material_(config.material) {
    validate_material(config_.material);
    if (!(config_.length_m > 0.0) || config_.nodes < 3U || !(config_.area_m2 > 0.0)
        || config_.max_gummel_iterations == 0U || !(config_.relative_tolerance > 0.0)
        || !(config_.under_relaxation > 0.0 && config_.under_relaxation <= 1.0)
        || !(config_.carrier_floor_m3 > 0.0)) {
        throw std::invalid_argument("invalid semiconductor device configuration");
    }
    const std::size_t n = config_.nodes;
    x_.resize(n);
    potential_v_.assign(n, 0.0);
    electron_m3_.assign(n, config_.material.intrinsic_density_m3);
    hole_m3_.assign(n, config_.material.intrinsic_density_m3);
    net_doping_m3_.assign(n, 0.0);
    const double dx = config_.length_m / static_cast<double>(n - 1U);
    for (std::size_t i = 0U; i < n; ++i) x_[i] = static_cast<double>(i) * dx;
}

void SemiconductorDevice1D::set_net_doping(std::span<const double> doping) {
    if (doping.size() != config_.nodes) throw std::invalid_argument("doping size mismatch");
    for (double value : doping) if (!std::isfinite(value)) throw std::invalid_argument("doping must be finite");
    net_doping_m3_.assign(doping.begin(), doping.end());
}

void SemiconductorDevice1D::set_net_doping(const std::function<double(double)>& doping) {
    if (!doping) throw std::invalid_argument("doping callback missing");
    for (std::size_t i = 0U; i < x_.size(); ++i) {
        net_doping_m3_[i] = doping(x_[i]);
        if (!std::isfinite(net_doping_m3_[i])) throw std::invalid_argument("doping callback produced non-finite value");
    }
}

void SemiconductorDevice1D::set_temperature(double temperature_k) {
    if (!(temperature_k > 0.0)) throw std::invalid_argument("temperature must be positive");
    config_.material.temperature_k = temperature_k;
    const auto& td = config_.temperature_dependence;
    if (!td.enabled) return;
    const double tref = reference_material_.temperature_k;
    const double eg_ref = silicon_bandgap_ev(tref, td.varshni_eg0_ev, td.varshni_alpha_ev_per_k, td.varshni_beta_k);
    const double eg = silicon_bandgap_ev(temperature_k, td.varshni_eg0_ev, td.varshni_alpha_ev_per_k, td.varshni_beta_k);
    config_.material.intrinsic_density_m3 = temperature_scaled_intrinsic_density(
        reference_material_.intrinsic_density_m3, tref, temperature_k, eg_ref, eg);
    config_.material.electron_mobility_m2_per_vs = temperature_scaled_mobility(
        reference_material_.electron_mobility_m2_per_vs, tref, temperature_k, td.electron_mobility_exponent);
    config_.material.hole_mobility_m2_per_vs = temperature_scaled_mobility(
        reference_material_.hole_mobility_m2_per_vs, tref, temperature_k, td.hole_mobility_exponent);
    const double dos_scale = std::pow(temperature_k / tref, td.density_of_states_exponent);
    config_.material.conduction_band_density_m3 = reference_material_.conduction_band_density_m3 * dos_scale;
    config_.material.valence_band_density_m3 = reference_material_.valence_band_density_m3 * dos_scale;
}

void SemiconductorDevice1D::set_contact_voltages(double left_v, double right_v) {
    if (!std::isfinite(left_v) || !std::isfinite(right_v)) throw std::invalid_argument("contact voltage must be finite");
    config_.left.voltage_v = left_v;
    config_.right.voltage_v = right_v;
}

double SemiconductorDevice1D::contact_electron_density(const Contact1D& c, double doping) const {
    const auto& m = config_.material;
    if (c.type == ContactType::ohmic) {
        return std::max(neutral_carriers(doping, m.intrinsic_density_m3).first, config_.carrier_floor_m3);
    }
    if (c.type == ContactType::schottky) {
        const double vt = thermal_voltage(m.temperature_k);
        return std::max(m.conduction_band_density_m3 * safe_exp(-c.electron_barrier_ev / vt),
                        config_.carrier_floor_m3);
    }
    throw std::logic_error("blocking contact has no prescribed electron density");
}

double SemiconductorDevice1D::contact_hole_density(const Contact1D& c, double doping) const {
    const auto& m = config_.material;
    if (c.type == ContactType::ohmic) {
        return std::max(neutral_carriers(doping, m.intrinsic_density_m3).second, config_.carrier_floor_m3);
    }
    if (c.type == ContactType::schottky) {
        const double n = contact_electron_density(c, doping);
        return std::max(m.intrinsic_density_m3 * m.intrinsic_density_m3 / n, config_.carrier_floor_m3);
    }
    throw std::logic_error("blocking contact has no prescribed hole density");
}

double SemiconductorDevice1D::contact_potential(const Contact1D& c, double doping) const {
    const double vt = thermal_voltage(config_.material.temperature_k);
    if (c.type == ContactType::ohmic) {
        const double n = contact_electron_density(c, doping);
        return vt * std::log(n / config_.material.intrinsic_density_m3) + c.voltage_v;
    }
    if (c.type == ContactType::schottky) return c.voltage_v;
    throw std::logic_error("non-Dirichlet contact has no prescribed electrostatic potential");
}

bool SemiconductorDevice1D::carrier_blocking(const Contact1D& c) const noexcept {
    return c.type == ContactType::insulating || c.type == ContactType::gate;
}

void SemiconductorDevice1D::enforce_contact_carriers() {
    if (!carrier_blocking(config_.left)) {
        electron_m3_.front() = contact_electron_density(config_.left, net_doping_m3_.front());
        hole_m3_.front() = contact_hole_density(config_.left, net_doping_m3_.front());
    }
    if (!carrier_blocking(config_.right)) {
        electron_m3_.back() = contact_electron_density(config_.right, net_doping_m3_.back());
        hole_m3_.back() = contact_hole_density(config_.right, net_doping_m3_.back());
    }
}

void SemiconductorDevice1D::apply_blocking_carrier_boundaries() {
    const double vt = thermal_voltage(config_.material.temperature_k);
    if (carrier_blocking(config_.left)) {
        const double a = (potential_v_[1U] - potential_v_[0U]) / vt;
        electron_m3_[0U] = std::max(electron_m3_[1U] * bernoulli_function(a) / bernoulli_function(-a),
                                    config_.carrier_floor_m3);
        hole_m3_[0U] = std::max(hole_m3_[1U] * bernoulli_function(-a) / bernoulli_function(a),
                                config_.carrier_floor_m3);
    }
    if (carrier_blocking(config_.right)) {
        const std::size_t n = config_.nodes;
        const double a = (potential_v_[n - 1U] - potential_v_[n - 2U]) / vt;
        electron_m3_[n - 1U] = std::max(electron_m3_[n - 2U] * bernoulli_function(-a) / bernoulli_function(a),
                                        config_.carrier_floor_m3);
        hole_m3_[n - 1U] = std::max(hole_m3_[n - 2U] * bernoulli_function(a) / bernoulli_function(-a),
                                    config_.carrier_floor_m3);
    }
}

void SemiconductorDevice1D::initialize_charge_neutral() {
    const double ni = config_.material.intrinsic_density_m3;
    const double vt = thermal_voltage(config_.material.temperature_k);
    for (std::size_t i = 0U; i < x_.size(); ++i) {
        const auto [n, p] = neutral_carriers(net_doping_m3_[i], ni);
        electron_m3_[i] = std::max(n, config_.carrier_floor_m3);
        hole_m3_[i] = std::max(p, config_.carrier_floor_m3);
        potential_v_[i] = vt * std::log(electron_m3_[i] / ni);
    }
    if (!carrier_blocking(config_.left)) {
        potential_v_.front() = contact_potential(config_.left, net_doping_m3_.front());
    }
    if (!carrier_blocking(config_.right)) {
        potential_v_.back() = contact_potential(config_.right, net_doping_m3_.back());
    }
    enforce_contact_carriers();
    apply_blocking_carrier_boundaries();
}

void SemiconductorDevice1D::solve_poisson_linear() {
    const std::size_t n = config_.nodes;
    const double dx = config_.length_m / static_cast<double>(n - 1U);
    const double eps = vacuum_permittivity_f_per_m * config_.material.relative_permittivity;
    if (config_.left.type == ContactType::insulating && config_.right.type == ContactType::insulating) {
        throw std::invalid_argument("Poisson problem with two insulating boundaries has no voltage reference");
    }
    std::vector<double> lo(n, -1.0), diag(n, 2.0), hi(n, -1.0), rhs(n, 0.0);
    for (std::size_t i = 1U; i + 1U < n; ++i) {
        const double rho = elementary_charge_c * (net_doping_m3_[i] + hole_m3_[i] - electron_m3_[i]);
        rhs[i] = rho * dx * dx / eps;
    }
    auto apply_left = [&](const Contact1D& c) {
        lo[0U] = 0.0;
        if (c.type == ContactType::ohmic || c.type == ContactType::schottky) {
            diag[0U] = 1.0;
            hi[0U] = 0.0;
            rhs[0U] = contact_potential(c, net_doping_m3_.front());
        } else if (c.type == ContactType::insulating) {
            diag[0U] = 1.0;
            hi[0U] = -1.0;
            rhs[0U] = 0.0;
        } else {
            if (!(c.oxide_capacitance_f_per_m2 > 0.0)) throw std::invalid_argument("gate contact requires positive oxide capacitance");
            const double gamma = c.oxide_capacitance_f_per_m2 * dx / eps;
            diag[0U] = 1.0 + gamma;
            hi[0U] = -1.0;
            rhs[0U] = gamma * (c.voltage_v - c.flatband_voltage_v);
        }
    };
    auto apply_right = [&](const Contact1D& c) {
        hi[n - 1U] = 0.0;
        if (c.type == ContactType::ohmic || c.type == ContactType::schottky) {
            lo[n - 1U] = 0.0;
            diag[n - 1U] = 1.0;
            rhs[n - 1U] = contact_potential(c, net_doping_m3_.back());
        } else if (c.type == ContactType::insulating) {
            lo[n - 1U] = -1.0;
            diag[n - 1U] = 1.0;
            rhs[n - 1U] = 0.0;
        } else {
            if (!(c.oxide_capacitance_f_per_m2 > 0.0)) throw std::invalid_argument("gate contact requires positive oxide capacitance");
            const double gamma = c.oxide_capacitance_f_per_m2 * dx / eps;
            lo[n - 1U] = -1.0;
            diag[n - 1U] = 1.0 + gamma;
            rhs[n - 1U] = gamma * (c.voltage_v - c.flatband_voltage_v);
        }
    };
    apply_left(config_.left);
    apply_right(config_.right);
    potential_v_ = solve_tridiagonal(std::move(lo), std::move(diag), std::move(hi), std::move(rhs));
}

void SemiconductorDevice1D::solve_electron_continuity() {
    const std::size_t n = config_.nodes;
    const std::size_t m = n - 2U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    std::vector<double> lo(m, 0.0), diag(m, 0.0), hi(m, 0.0), rhs(m, 0.0);
    for (std::size_t k = 0U; k < m; ++k) {
        const std::size_t i = k + 1U;
        const double am = (potential_v_[i] - potential_v_[i - 1U]) / vt;
        const double ap = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        lo[k] = -bernoulli_function(-am);
        diag[k] = bernoulli_function(am) + bernoulli_function(-ap);
        hi[k] = -bernoulli_function(ap);
        const double recomb = total_recombination(electron_m3_[i], hole_m3_[i],
                                                  config_.material.intrinsic_density_m3, config_.recombination);
        const double dx = config_.length_m / static_cast<double>(n - 1U);
        const double scale = dx * dx / (config_.material.electron_mobility_m2_per_vs * vt);
        rhs[k] = -recomb * scale;
    }
    if (carrier_blocking(config_.left)) {
        const double a = (potential_v_[1U] - potential_v_[0U]) / vt;
        diag.front() += lo.front() * bernoulli_function(a) / bernoulli_function(-a);
    } else {
        rhs.front() -= lo.front() * contact_electron_density(config_.left, net_doping_m3_.front());
    }
    lo.front() = 0.0;
    if (carrier_blocking(config_.right)) {
        const double a = (potential_v_[n - 1U] - potential_v_[n - 2U]) / vt;
        diag.back() += hi.back() * bernoulli_function(-a) / bernoulli_function(a);
    } else {
        rhs.back() -= hi.back() * contact_electron_density(config_.right, net_doping_m3_.back());
    }
    hi.back() = 0.0;
    auto interior = solve_tridiagonal(std::move(lo), std::move(diag), std::move(hi), std::move(rhs));
    for (std::size_t k = 0U; k < m; ++k) electron_m3_[k + 1U] = std::max(interior[k], config_.carrier_floor_m3);
    enforce_contact_carriers();
    apply_blocking_carrier_boundaries();
}

void SemiconductorDevice1D::solve_hole_continuity() {
    const std::size_t n = config_.nodes;
    const std::size_t m = n - 2U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    std::vector<double> lo(m, 0.0), diag(m, 0.0), hi(m, 0.0), rhs(m, 0.0);
    for (std::size_t k = 0U; k < m; ++k) {
        const std::size_t i = k + 1U;
        const double am = (potential_v_[i] - potential_v_[i - 1U]) / vt;
        const double ap = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        lo[k] = -bernoulli_function(am);
        diag[k] = bernoulli_function(-am) + bernoulli_function(ap);
        hi[k] = -bernoulli_function(-ap);
        const double recomb = total_recombination(electron_m3_[i], hole_m3_[i],
                                                  config_.material.intrinsic_density_m3, config_.recombination);
        const double dx = config_.length_m / static_cast<double>(n - 1U);
        const double scale = dx * dx / (config_.material.hole_mobility_m2_per_vs * vt);
        rhs[k] = -recomb * scale;
    }
    if (carrier_blocking(config_.left)) {
        const double a = (potential_v_[1U] - potential_v_[0U]) / vt;
        diag.front() += lo.front() * bernoulli_function(-a) / bernoulli_function(a);
    } else {
        rhs.front() -= lo.front() * contact_hole_density(config_.left, net_doping_m3_.front());
    }
    lo.front() = 0.0;
    if (carrier_blocking(config_.right)) {
        const double a = (potential_v_[n - 1U] - potential_v_[n - 2U]) / vt;
        diag.back() += hi.back() * bernoulli_function(a) / bernoulli_function(-a);
    } else {
        rhs.back() -= hi.back() * contact_hole_density(config_.right, net_doping_m3_.back());
    }
    hi.back() = 0.0;
    auto interior = solve_tridiagonal(std::move(lo), std::move(diag), std::move(hi), std::move(rhs));
    for (std::size_t k = 0U; k < m; ++k) hole_m3_[k + 1U] = std::max(interior[k], config_.carrier_floor_m3);
    enforce_contact_carriers();
    apply_blocking_carrier_boundaries();
}

void SemiconductorDevice1D::solve_electron_continuity_transient(std::span<const double> old_density, double dt_s) {
    if (old_density.size() != config_.nodes || !(dt_s > 0.0)) throw std::invalid_argument("invalid electron transient state");
    const std::size_t n = config_.nodes;
    const std::size_t m = n - 2U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    const double dx = config_.length_m / static_cast<double>(n - 1U);
    const double diffusion = config_.material.electron_mobility_m2_per_vs * vt;
    const double accumulation = dx * dx / (diffusion * dt_s);
    const double recombination_scale = dx * dx / diffusion;
    std::vector<double> lo(m, 0.0), diag(m, 0.0), hi(m, 0.0), rhs(m, 0.0);
    for (std::size_t k = 0U; k < m; ++k) {
        const std::size_t i = k + 1U;
        const double am = (potential_v_[i] - potential_v_[i - 1U]) / vt;
        const double ap = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        lo[k] = -bernoulli_function(-am);
        diag[k] = bernoulli_function(am) + bernoulli_function(-ap) + accumulation;
        hi[k] = -bernoulli_function(ap);
        const double recomb = total_recombination(electron_m3_[i], hole_m3_[i],
                                                  config_.material.intrinsic_density_m3,
                                                  config_.recombination);
        rhs[k] = accumulation * old_density[i] - recomb * recombination_scale;
    }
    if (carrier_blocking(config_.left)) {
        const double a = (potential_v_[1U] - potential_v_[0U]) / vt;
        diag.front() += lo.front() * bernoulli_function(a) / bernoulli_function(-a);
    } else {
        rhs.front() -= lo.front() * contact_electron_density(config_.left, net_doping_m3_.front());
    }
    lo.front() = 0.0;
    if (carrier_blocking(config_.right)) {
        const double a = (potential_v_[n - 1U] - potential_v_[n - 2U]) / vt;
        diag.back() += hi.back() * bernoulli_function(-a) / bernoulli_function(a);
    } else {
        rhs.back() -= hi.back() * contact_electron_density(config_.right, net_doping_m3_.back());
    }
    hi.back() = 0.0;
    const auto interior = solve_tridiagonal(std::move(lo), std::move(diag), std::move(hi), std::move(rhs));
    for (std::size_t k = 0U; k < m; ++k) electron_m3_[k + 1U] = std::max(interior[k], config_.carrier_floor_m3);
    enforce_contact_carriers();
    apply_blocking_carrier_boundaries();
}

void SemiconductorDevice1D::solve_hole_continuity_transient(std::span<const double> old_density, double dt_s) {
    if (old_density.size() != config_.nodes || !(dt_s > 0.0)) throw std::invalid_argument("invalid hole transient state");
    const std::size_t n = config_.nodes;
    const std::size_t m = n - 2U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    const double dx = config_.length_m / static_cast<double>(n - 1U);
    const double diffusion = config_.material.hole_mobility_m2_per_vs * vt;
    const double accumulation = dx * dx / (diffusion * dt_s);
    const double recombination_scale = dx * dx / diffusion;
    std::vector<double> lo(m, 0.0), diag(m, 0.0), hi(m, 0.0), rhs(m, 0.0);
    for (std::size_t k = 0U; k < m; ++k) {
        const std::size_t i = k + 1U;
        const double am = (potential_v_[i] - potential_v_[i - 1U]) / vt;
        const double ap = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        lo[k] = -bernoulli_function(am);
        diag[k] = bernoulli_function(-am) + bernoulli_function(ap) + accumulation;
        hi[k] = -bernoulli_function(-ap);
        const double recomb = total_recombination(electron_m3_[i], hole_m3_[i],
                                                  config_.material.intrinsic_density_m3,
                                                  config_.recombination);
        rhs[k] = accumulation * old_density[i] - recomb * recombination_scale;
    }
    if (carrier_blocking(config_.left)) {
        const double a = (potential_v_[1U] - potential_v_[0U]) / vt;
        diag.front() += lo.front() * bernoulli_function(-a) / bernoulli_function(a);
    } else {
        rhs.front() -= lo.front() * contact_hole_density(config_.left, net_doping_m3_.front());
    }
    lo.front() = 0.0;
    if (carrier_blocking(config_.right)) {
        const double a = (potential_v_[n - 1U] - potential_v_[n - 2U]) / vt;
        diag.back() += hi.back() * bernoulli_function(a) / bernoulli_function(-a);
    } else {
        rhs.back() -= hi.back() * contact_hole_density(config_.right, net_doping_m3_.back());
    }
    hi.back() = 0.0;
    const auto interior = solve_tridiagonal(std::move(lo), std::move(diag), std::move(hi), std::move(rhs));
    for (std::size_t k = 0U; k < m; ++k) hole_m3_[k + 1U] = std::max(interior[k], config_.carrier_floor_m3);
    enforce_contact_carriers();
    apply_blocking_carrier_boundaries();
}

double SemiconductorDevice1D::maximum_relative_change(std::span<const double> oldv, std::span<const double> newv) const {
    if (oldv.size() != newv.size()) throw std::invalid_argument("relative-change size mismatch");
    double result = 0.0;
    for (std::size_t i = 0U; i < oldv.size(); ++i) {
        const double scale = std::max({std::abs(oldv[i]), std::abs(newv[i]), config_.carrier_floor_m3});
        result = std::max(result, std::abs(newv[i] - oldv[i]) / scale);
    }
    return result;
}

DcResult SemiconductorDevice1D::solve_equilibrium() {
    initialize_charge_neutral();
    const double vt = thermal_voltage(config_.material.temperature_k);
    const double ni = config_.material.intrinsic_density_m3;
    double change = std::numeric_limits<double>::infinity();
    for (std::size_t it = 0U; it < config_.max_gummel_iterations; ++it) {
        const auto oldp = potential_v_;
        const auto oldn = electron_m3_;
        const auto oldh = hole_m3_;
        solve_poisson_linear();
        for (std::size_t i = 0U; i < x_.size(); ++i) {
            const double n_target = std::max(ni * safe_exp(potential_v_[i] / vt), config_.carrier_floor_m3);
            const double p_target = std::max(ni * safe_exp(-potential_v_[i] / vt), config_.carrier_floor_m3);
            electron_m3_[i] = (1.0 - config_.under_relaxation) * oldn[i] + config_.under_relaxation * n_target;
            hole_m3_[i] = (1.0 - config_.under_relaxation) * oldh[i] + config_.under_relaxation * p_target;
        }
        enforce_contact_carriers();
        apply_blocking_carrier_boundaries();
        const double cp = maximum_relative_change(oldp, potential_v_);
        const double cn = maximum_relative_change(oldn, electron_m3_);
        const double ch = maximum_relative_change(oldh, hole_m3_);
        change = std::max({cp, cn, ch});
        if (change < config_.relative_tolerance) return {true, it + 1U, change, terminal_current_a()};
    }
    return {false, config_.max_gummel_iterations, change, terminal_current_a()};
}

DcResult SemiconductorDevice1D::solve_dc() {
    if (electron_m3_.empty()) initialize_charge_neutral();
    enforce_contact_carriers();
    double change = std::numeric_limits<double>::infinity();
    for (std::size_t it = 0U; it < config_.max_gummel_iterations; ++it) {
        const auto oldp = potential_v_;
        const auto oldn = electron_m3_;
        const auto oldh = hole_m3_;
        solve_poisson_linear();
        solve_electron_continuity();
        solve_hole_continuity();
        for (std::size_t i = 1U; i + 1U < x_.size(); ++i) {
            electron_m3_[i] = (1.0 - config_.under_relaxation) * oldn[i] + config_.under_relaxation * electron_m3_[i];
            hole_m3_[i] = (1.0 - config_.under_relaxation) * oldh[i] + config_.under_relaxation * hole_m3_[i];
        }
        enforce_contact_carriers();
        apply_blocking_carrier_boundaries();
        const double cp = maximum_relative_change(oldp, potential_v_);
        const double cn = maximum_relative_change(oldn, electron_m3_);
        const double ch = maximum_relative_change(oldh, hole_m3_);
        change = std::max({cp, cn, ch});
        if (change < config_.relative_tolerance) return {true, it + 1U, change, terminal_current_a()};
    }
    return {false, config_.max_gummel_iterations, change, terminal_current_a()};
}

std::vector<IvPoint> SemiconductorDevice1D::sweep_right_contact(double start_v, double stop_v, std::size_t points) {
    if (points < 2U || !std::isfinite(start_v) || !std::isfinite(stop_v)) throw std::invalid_argument("invalid IV sweep");
    std::vector<IvPoint> out;
    out.reserve(points);
    const double left = config_.left.voltage_v;
    for (std::size_t k = 0U; k < points; ++k) {
        const double t = static_cast<double>(k) / static_cast<double>(points - 1U);
        config_.right.voltage_v = start_v + t * (stop_v - start_v);
        const auto result = solve_dc();
        out.push_back({config_.right.voltage_v - left, result.terminal_current_a, result.converged});
    }
    return out;
}


TransientStepResult SemiconductorDevice1D::step_transient(double dt_s) {
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) throw std::invalid_argument("transient timestep must be positive and finite");
    const auto old_potential = potential_v_;
    const auto old_electron = electron_m3_;
    const auto old_hole = hole_m3_;
    const auto old_field = electric_field_faces();
    enforce_contact_carriers();
    double change = std::numeric_limits<double>::infinity();
    std::size_t iterations = config_.max_gummel_iterations;
    bool converged = false;
    for (std::size_t it = 0U; it < config_.max_gummel_iterations; ++it) {
        const auto iter_potential = potential_v_;
        const auto iter_electron = electron_m3_;
        const auto iter_hole = hole_m3_;
        solve_poisson_linear();
        solve_electron_continuity_transient(old_electron, dt_s);
        solve_hole_continuity_transient(old_hole, dt_s);
        for (std::size_t i = 1U; i + 1U < x_.size(); ++i) {
            electron_m3_[i] = (1.0 - config_.under_relaxation) * iter_electron[i]
                              + config_.under_relaxation * electron_m3_[i];
            hole_m3_[i] = (1.0 - config_.under_relaxation) * iter_hole[i]
                          + config_.under_relaxation * hole_m3_[i];
        }
        enforce_contact_carriers();
        apply_blocking_carrier_boundaries();
        const double cp = maximum_relative_change(iter_potential, potential_v_);
        const double cn = maximum_relative_change(iter_electron, electron_m3_);
        const double ch = maximum_relative_change(iter_hole, hole_m3_);
        change = std::max({cp, cn, ch});
        if (change < config_.relative_tolerance) {
            iterations = it + 1U;
            converged = true;
            break;
        }
    }
    if (!converged) {
        potential_v_ = old_potential;
        electron_m3_ = old_electron;
        hole_m3_ = old_hole;
        return {false, iterations, change, time_s_, 0.0, 0.0, 0.0};
    }
    const auto new_field = electric_field_faces();
    const auto conduction_faces = total_current_density_faces();
    const double conduction_current = conduction_faces.front() * config_.area_m2;
    const double eps = vacuum_permittivity_f_per_m * config_.material.relative_permittivity;
    const double displacement_current = eps * config_.area_m2 * (new_field.front() - old_field.front()) / dt_s;
    time_s_ += dt_s;
    return {true, iterations, change, time_s_, conduction_current, displacement_current,
            conduction_current + displacement_current};
}

SmallSignalResult SemiconductorDevice1D::small_signal(double frequency_hz, double perturbation_v) const {
    if (!(frequency_hz >= 0.0) || !std::isfinite(frequency_hz) || !(perturbation_v > 0.0)
        || !std::isfinite(perturbation_v)) {
        throw std::invalid_argument("invalid small-signal parameters");
    }
    auto plus = *this;
    auto minus = *this;
    plus.config_.left.voltage_v += perturbation_v;
    minus.config_.left.voltage_v -= perturbation_v;
    const auto plus_dc = plus.solve_dc();
    const auto minus_dc = minus.solve_dc();
    const double dv = 2.0 * perturbation_v;
    const double conductance = (plus_dc.terminal_current_a - minus_dc.terminal_current_a) / dv;
    const double capacitance = (plus.left_terminal_charge_c() - minus.left_terminal_charge_c()) / dv;
    const double omega = 2.0 * std::acos(-1.0) * frequency_hz;
    const std::complex<double> admittance{conductance, omega * capacitance};
    std::complex<double> impedance{std::numeric_limits<double>::infinity(), 0.0};
    if (std::abs(admittance) > 0.0) impedance = 1.0 / admittance;
    return {config_.left.voltage_v - config_.right.voltage_v, frequency_hz, conductance, capacitance,
            admittance, impedance, plus_dc.converged && minus_dc.converged};
}

std::vector<CvPoint> SemiconductorDevice1D::sweep_cv(double start_left_v,
                                                     double stop_left_v,
                                                     std::size_t points,
                                                     double perturbation_v) const {
    if (points < 2U || !std::isfinite(start_left_v) || !std::isfinite(stop_left_v)
        || !(perturbation_v > 0.0) || !std::isfinite(perturbation_v)) {
        throw std::invalid_argument("invalid CV sweep parameters");
    }
    std::vector<CvPoint> out;
    out.reserve(points);
    for (std::size_t k = 0U; k < points; ++k) {
        const double t = static_cast<double>(k) / static_cast<double>(points - 1U);
        auto bias = *this;
        bias.config_.left.voltage_v = start_left_v + t * (stop_left_v - start_left_v);
        const auto dc = bias.solve_dc();
        const auto ac = bias.small_signal(0.0, perturbation_v);
        out.push_back({bias.config_.left.voltage_v - bias.config_.right.voltage_v,
                       ac.capacitance_f, bias.left_terminal_charge_c(), dc.converged && ac.converged});
    }
    return out;
}

std::vector<double> SemiconductorDevice1D::electron_current_density_faces() const {
    const std::size_t nf = config_.nodes - 1U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    const double dx = config_.length_m / static_cast<double>(nf);
    std::vector<double> out(nf, 0.0);
    const double prefactor = elementary_charge_c * config_.material.electron_mobility_m2_per_vs * vt / dx;
    for (std::size_t i = 0U; i < nf; ++i) {
        const double a = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        out[i] = prefactor * (electron_m3_[i + 1U] * bernoulli_function(a)
                            - electron_m3_[i] * bernoulli_function(-a));
    }
    if (carrier_blocking(config_.left)) out.front() = 0.0;
    if (carrier_blocking(config_.right)) out.back() = 0.0;
    return out;
}

std::vector<double> SemiconductorDevice1D::hole_current_density_faces() const {
    const std::size_t nf = config_.nodes - 1U;
    const double vt = thermal_voltage(config_.material.temperature_k);
    const double dx = config_.length_m / static_cast<double>(nf);
    std::vector<double> out(nf, 0.0);
    const double prefactor = elementary_charge_c * config_.material.hole_mobility_m2_per_vs * vt / dx;
    for (std::size_t i = 0U; i < nf; ++i) {
        const double a = (potential_v_[i + 1U] - potential_v_[i]) / vt;
        out[i] = prefactor * (hole_m3_[i] * bernoulli_function(a)
                            - hole_m3_[i + 1U] * bernoulli_function(-a));
    }
    if (carrier_blocking(config_.left)) out.front() = 0.0;
    if (carrier_blocking(config_.right)) out.back() = 0.0;
    return out;
}

std::vector<double> SemiconductorDevice1D::total_current_density_faces() const {
    auto out = electron_current_density_faces();
    const auto holes = hole_current_density_faces();
    for (std::size_t i = 0U; i < out.size(); ++i) out[i] += holes[i];
    return out;
}

std::vector<double> SemiconductorDevice1D::electric_field_faces() const {
    const std::size_t nf = config_.nodes - 1U;
    const double dx = config_.length_m / static_cast<double>(nf);
    std::vector<double> out(nf, 0.0);
    for (std::size_t i = 0U; i < nf; ++i) out[i] = -(potential_v_[i + 1U] - potential_v_[i]) / dx;
    return out;
}

std::vector<double> SemiconductorDevice1D::joule_heating_cells() const {
    const auto current = total_current_density_faces();
    const auto field = electric_field_faces();
    std::vector<double> face_heat(current.size(), 0.0);
    for (std::size_t i = 0U; i < current.size(); ++i) face_heat[i] = current[i] * field[i];
    std::vector<double> out(config_.nodes, 0.0);
    out.front() = face_heat.front();
    out.back() = face_heat.back();
    for (std::size_t i = 1U; i + 1U < config_.nodes; ++i) {
        out[i] = 0.5 * (face_heat[i - 1U] + face_heat[i]);
    }
    return out;
}

double SemiconductorDevice1D::terminal_current_a() const {
    const auto current = total_current_density_faces();
    if (current.empty()) return 0.0;
    return 0.5 * (current.front() + current.back()) * config_.area_m2;
}


double SemiconductorDevice1D::left_terminal_charge_c() const {
    if (config_.left.type == ContactType::gate) {
        return config_.left.oxide_capacitance_f_per_m2 * config_.area_m2
             * (config_.left.voltage_v - config_.left.flatband_voltage_v - potential_v_.front());
    }
    if (config_.left.type == ContactType::insulating) return 0.0;
    const auto field = electric_field_faces();
    if (field.empty()) return 0.0;
    const double eps = vacuum_permittivity_f_per_m * config_.material.relative_permittivity;
    return eps * config_.area_m2 * field.front();
}

double SemiconductorDevice1D::total_mobile_charge_c() const {
    const double dx = config_.length_m / static_cast<double>(config_.nodes - 1U);
    double integral = 0.0;
    for (std::size_t i = 0U; i < config_.nodes; ++i) {
        const double weight = (i == 0U || i + 1U == config_.nodes) ? 0.5 : 1.0;
        integral += weight * (hole_m3_[i] - electron_m3_[i]);
    }
    return elementary_charge_c * config_.area_m2 * dx * integral;
}

std::vector<double> pn_junction_doping(std::span<const double> x, double junction, double acceptor, double donor) {
    if (x.empty() || !(acceptor > 0.0) || !(donor > 0.0)) throw std::invalid_argument("invalid PN doping parameters");
    std::vector<double> out(x.size());
    for (std::size_t i = 0U; i < x.size(); ++i) out[i] = x[i] < junction ? -acceptor : donor;
    return out;
}

std::vector<double> pin_diode_doping(std::span<const double> x, double p_end, double n_begin,
                                     double acceptor, double donor, double intrinsic_net) {
    if (x.empty() || !(p_end >= 0.0) || !(n_begin > p_end) || !(acceptor > 0.0) || !(donor > 0.0)) {
        throw std::invalid_argument("invalid PIN doping parameters");
    }
    std::vector<double> out(x.size());
    for (std::size_t i = 0U; i < x.size(); ++i) {
        out[i] = x[i] < p_end ? -acceptor : (x[i] >= n_begin ? donor : intrinsic_net);
    }
    return out;
}

} // namespace cfd::tcad
