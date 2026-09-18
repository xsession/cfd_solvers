#include "cfd/battery/lithium_ion.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <stdexcept>

namespace cfd::battery {
namespace {

double clamp01(double x) noexcept { return std::clamp(x, 1.0e-6, 1.0 - 1.0e-6); }

void solve_tridiagonal(std::span<const double> lower,
                       std::span<const double> diagonal,
                       std::span<const double> upper,
                       std::span<double> rhs) {
    const std::size_t n = rhs.size();
    if (diagonal.size() != n || lower.size() + 1U != n || upper.size() + 1U != n) {
        throw std::invalid_argument("tridiagonal system size mismatch");
    }
    std::vector<double> cprime(n > 1U ? n - 1U : 0U, 0.0);
    double pivot = diagonal[0];
    if (std::abs(pivot) < 1.0e-30) throw std::runtime_error("singular tridiagonal system");
    if (n > 1U) cprime[0] = upper[0] / pivot;
    rhs[0] /= pivot;
    for (std::size_t i = 1U; i < n; ++i) {
        pivot = diagonal[i] - lower[i - 1U] * cprime[i - 1U];
        if (std::abs(pivot) < 1.0e-30) throw std::runtime_error("singular tridiagonal system");
        if (i + 1U < n) cprime[i] = upper[i] / pivot;
        rhs[i] = (rhs[i] - lower[i - 1U] * rhs[i - 1U]) / pivot;
    }
    for (std::size_t i = n - 1U; i > 0U; --i) rhs[i - 1U] -= cprime[i - 1U] * rhs[i];
}

double effective_transport(double bulk, double porosity, double exponent) {
    return bulk * std::pow(std::clamp(porosity, 1.0e-6, 1.0), exponent);
}

double region_porosity(const ElectrolyteConfig& e, double x_m) {
    if (x_m < e.negative_thickness_m) return e.negative_porosity;
    if (x_m < e.negative_thickness_m + e.separator_thickness_m) return e.separator_porosity;
    return e.positive_porosity;
}

} // namespace

SphericalDiffusionParticle::SphericalDiffusionParticle(ParticleConfig config) : config_(config) {
    if (!(config_.radius_m > 0.0) || !(config_.diffusivity_m2_per_s > 0.0)
        || !(config_.maximum_concentration_mol_per_m3 > 0.0) || config_.radial_cells < 3U) {
        throw std::invalid_argument("invalid spherical-particle configuration");
    }
    reset(config_.initial_stoichiometry);
}

void SphericalDiffusionParticle::reset(double stoichiometry) {
    const double c = clamp01(stoichiometry) * config_.maximum_concentration_mol_per_m3;
    concentration_.assign(config_.radial_cells, c);
}

void SphericalDiffusionParticle::step(double outward_flux, double dt) {
    if (!(dt > 0.0) || !std::isfinite(outward_flux)) throw std::invalid_argument("invalid particle step");
    const std::size_t n = concentration_.size();
    const double dr = config_.radius_m / static_cast<double>(n);
    const double d = config_.diffusivity_m2_per_s;
    std::vector<double> lower(n - 1U, 0.0), diag(n, 1.0), upper(n - 1U, 0.0), rhs = concentration_;
    for (std::size_t i = 0U; i + 1U < n; ++i) {
        const double rf = static_cast<double>(i + 1U) * dr;
        const double area = 4.0 * std::numbers::pi * rf * rf;
        const double rlo_i = static_cast<double>(i) * dr;
        const double rhi_i = static_cast<double>(i + 1U) * dr;
        const double vi = (4.0 / 3.0) * std::numbers::pi * (rhi_i * rhi_i * rhi_i - rlo_i * rlo_i * rlo_i);
        const double rlo_j = rhi_i;
        const double rhi_j = static_cast<double>(i + 2U) * dr;
        const double vj = (4.0 / 3.0) * std::numbers::pi * (rhi_j * rhi_j * rhi_j - rlo_j * rlo_j * rlo_j);
        const double ki = d * area / (dr * vi);
        const double kj = d * area / (dr * vj);
        diag[i] += dt * ki;
        upper[i] -= dt * ki;
        lower[i] -= dt * kj;
        diag[i + 1U] += dt * kj;
    }
    const double surface_area = 4.0 * std::numbers::pi * config_.radius_m * config_.radius_m;
    const double rlo = config_.radius_m - dr;
    const double v_last = (4.0 / 3.0) * std::numbers::pi
                        * (config_.radius_m * config_.radius_m * config_.radius_m - rlo * rlo * rlo);
    rhs.back() -= dt * surface_area * outward_flux / v_last;
    solve_tridiagonal(lower, diag, upper, rhs);
    const double max_c = config_.maximum_concentration_mol_per_m3;
    for (double& c : rhs) c = std::clamp(c, 1.0e-12 * max_c, (1.0 - 1.0e-12) * max_c);
    concentration_ = std::move(rhs);
}

double SphericalDiffusionParticle::average_concentration() const noexcept {
    const std::size_t n = concentration_.size();
    const double dr = config_.radius_m / static_cast<double>(n);
    double amount = 0.0, volume = 0.0;
    for (std::size_t i = 0U; i < n; ++i) {
        const double r0 = static_cast<double>(i) * dr;
        const double r1 = static_cast<double>(i + 1U) * dr;
        const double v = r1 * r1 * r1 - r0 * r0 * r0;
        amount += concentration_[i] * v;
        volume += v;
    }
    return amount / volume;
}

double SphericalDiffusionParticle::surface_concentration(double outward_flux) const noexcept {
    const double dr = config_.radius_m / static_cast<double>(concentration_.size());
    const double extrapolated = concentration_.back() - outward_flux * 0.5 * dr / config_.diffusivity_m2_per_s;
    return std::clamp(extrapolated, 1.0e-12 * config_.maximum_concentration_mol_per_m3,
                      (1.0 - 1.0e-12) * config_.maximum_concentration_mol_per_m3);
}

double SphericalDiffusionParticle::average_stoichiometry() const noexcept {
    return average_concentration() / config_.maximum_concentration_mol_per_m3;
}

double SphericalDiffusionParticle::surface_stoichiometry(double outward_flux) const noexcept {
    return surface_concentration(outward_flux) / config_.maximum_concentration_mol_per_m3;
}

double graphite_ocv_v(double x) {
    x = clamp01(x);
    return 0.08 + 0.75 * std::exp(-12.0 * x) + 0.12 * std::exp(-60.0 * (x - 0.12) * (x - 0.12))
         + 0.04 * (1.0 - x);
}

double nmc_ocv_v(double x) {
    x = clamp01(x);
    return 4.45 - 1.15 * x - 0.10 * std::tanh((x - 0.55) / 0.08) + 0.03 * std::sin(6.0 * x);
}

LithiumIonCellConfig default_graphite_nmc_config() {
    LithiumIonCellConfig c;
    c.area_m2 = 0.010;
    c.negative.particle = {5.5e-6, 3.0e-14, 3.1e4, 0.82, 24U};
    c.negative.thickness_m = 80.0e-6;
    c.negative.active_volume_fraction = 0.62;
    c.negative.solid_conductivity_s_per_m = 100.0;
    c.negative.exchange_current_density_a_per_m2 = 4.0;
    c.negative.open_circuit_voltage_v = graphite_ocv_v;
    c.positive.particle = {5.0e-6, 1.0e-13, 5.1e4, 0.42, 24U};
    c.positive.thickness_m = 75.0e-6;
    c.positive.active_volume_fraction = 0.58;
    c.positive.solid_conductivity_s_per_m = 10.0;
    c.positive.exchange_current_density_a_per_m2 = 4.0;
    c.positive.open_circuit_voltage_v = nmc_ocv_v;
    c.electrolyte.negative_thickness_m = c.negative.thickness_m;
    c.electrolyte.positive_thickness_m = c.positive.thickness_m;
    c.contact_resistance_ohm = 0.012;
    return c;
}

double symmetric_butler_volmer_overpotential_v(double j, double i0, double temperature) {
    if (!(i0 > 0.0) || !(temperature > 0.0)) throw std::invalid_argument("invalid Butler-Volmer parameters");
    return 2.0 * gas_constant * temperature / faraday_constant * std::asinh(j / (2.0 * i0));
}

ThermalSlab1D::ThermalSlab1D(ThermalSlabConfig config) : config_(config) {
    if (!(config_.thickness_m > 0.0) || config_.cells < 3U || !(config_.density_kg_per_m3 > 0.0)
        || !(config_.heat_capacity_j_per_kg_k > 0.0) || !(config_.conductivity_w_per_m_k >= 0.0)) {
        throw std::invalid_argument("invalid thermal-slab configuration");
    }
    reset(config_.ambient_temperature_k);
}

void ThermalSlab1D::reset(double temperature) { temperature_k_.assign(config_.cells, temperature); }

void ThermalSlab1D::step(double dt, std::span<const double> q) {
    if (!(dt > 0.0) || q.size() != temperature_k_.size()) throw std::invalid_argument("invalid thermal-slab step");
    const std::size_t n = temperature_k_.size();
    const double dx = config_.thickness_m / static_cast<double>(n);
    const double rho_cp = config_.density_kg_per_m3 * config_.heat_capacity_j_per_kg_k;
    const double a = config_.conductivity_w_per_m_k / (rho_cp * dx * dx);
    std::vector<double> lower(n - 1U, -dt * a), diag(n, 1.0 + 2.0 * dt * a), upper(n - 1U, -dt * a), rhs(n);
    for (std::size_t i = 0U; i < n; ++i) rhs[i] = temperature_k_[i] + dt * q[i] / rho_cp;
    const double beta = config_.convection_w_per_m2_k / (rho_cp * dx);
    diag.front() = 1.0 + dt * (a + beta);
    upper.front() = -dt * a;
    rhs.front() += dt * beta * config_.ambient_temperature_k;
    diag.back() = 1.0 + dt * (a + beta);
    lower.back() = -dt * a;
    rhs.back() += dt * beta * config_.ambient_temperature_k;
    solve_tridiagonal(lower, diag, upper, rhs);
    temperature_k_ = std::move(rhs);
}

double ThermalSlab1D::average_temperature() const noexcept {
    double s = 0.0; for (double t : temperature_k_) s += t; return s / static_cast<double>(temperature_k_.size());
}

SingleParticleModel::SingleParticleModel(LithiumIonCellConfig config)
    : config_(std::move(config)), negative_(config_.negative.particle), positive_(config_.positive.particle) {
    if (!(config_.area_m2 > 0.0)) throw std::invalid_argument("invalid battery area");
    if (!config_.negative.open_circuit_voltage_v) config_.negative.open_circuit_voltage_v = graphite_ocv_v;
    if (!config_.positive.open_circuit_voltage_v) config_.positive.open_circuit_voltage_v = nmc_ocv_v;
    reset();
}

void SingleParticleModel::set_external_temperature(double t) {
    if (!std::isfinite(t) || t <= 0.0) throw std::invalid_argument("invalid external cell temperature");
    temperature_k_ = t;
    external_thermal_ = true;
}

void SingleParticleModel::reset() {
    external_thermal_ = false;
    negative_.reset(config_.negative.particle.initial_stoichiometry);
    positive_.reset(config_.positive.particle.initial_stoichiometry);
    degradation_ = {config_.degradation.initial_sei_thickness_m, 0.0, 1.0, 0.0};
    temperature_k_ = config_.initial_temperature_k;
    time_s_ = 0.0;
}

double SingleParticleModel::negative_flux(double current) const noexcept {
    const double a = 3.0 * config_.negative.active_volume_fraction / config_.negative.particle.radius_m;
    return current / (faraday_constant * config_.area_m2 * config_.negative.thickness_m * a);
}

double SingleParticleModel::positive_flux(double current) const noexcept {
    const double a = 3.0 * config_.positive.active_volume_fraction / config_.positive.particle.radius_m;
    return -current / (faraday_constant * config_.area_m2 * config_.positive.thickness_m * a);
}

BatteryStepResult SingleParticleModel::evaluate(double current, double electrolyte_eta, double cmin, double cmax) const {
    const double jn_flux = negative_flux(current);
    const double jp_flux = positive_flux(current);
    const double xn = clamp01(negative_.surface_stoichiometry(jn_flux));
    const double xp = clamp01(positive_.surface_stoichiometry(jp_flux));
    const double un = config_.negative.open_circuit_voltage_v(xn);
    const double up = config_.positive.open_circuit_voltage_v(xp);
    const double ocv = up - un;
    const double an = 3.0 * config_.negative.active_volume_fraction / config_.negative.particle.radius_m;
    const double ap = 3.0 * config_.positive.active_volume_fraction / config_.positive.particle.radius_m;
    const double jn = current / (config_.area_m2 * config_.negative.thickness_m * an);
    const double jp = current / (config_.area_m2 * config_.positive.thickness_m * ap);
    const double i0n = config_.negative.exchange_current_density_a_per_m2 * 2.0 * std::sqrt(xn * (1.0 - xn));
    const double i0p = config_.positive.exchange_current_density_a_per_m2 * 2.0 * std::sqrt(xp * (1.0 - xp));
    const double reaction = std::abs(symmetric_butler_volmer_overpotential_v(jn, std::max(i0n, 1.0e-12), temperature_k_))
                          + std::abs(symmetric_butler_volmer_overpotential_v(jp, std::max(i0p, 1.0e-12), temperature_k_));
    const double solid_r = config_.negative.thickness_m / (config_.negative.solid_conductivity_s_per_m * config_.area_m2)
                         + config_.positive.thickness_m / (config_.positive.solid_conductivity_s_per_m * config_.area_m2);
    const double ohmic = current * (config_.contact_resistance_ohm + solid_r);
    const double voltage = ocv - std::copysign(reaction + std::abs(electrolyte_eta), current == 0.0 ? 1.0 : current) - ohmic;
    const double heat = std::abs(current) * std::max(0.0, std::abs(ocv - voltage));
    return {time_s_, current, voltage, ocv, reaction, electrolyte_eta, ohmic, xn, xp, cmin, cmax, temperature_k_, heat};
}

void SingleParticleModel::advance_degradation(double current, double dt, double gradient) {
    const auto& d = config_.degradation;
    if (d.enable_sei) {
        degradation_.sei_thickness_m += dt * d.sei_rate_m2_per_s / std::max(degradation_.sei_thickness_m, 1.0e-10);
    }
    if (d.enable_plating && current < 0.0) {
        const double charge = -current * dt;
        degradation_.plated_lithium_mol += d.plating_efficiency * charge / faraday_constant;
    }
    if (d.enable_active_material_loss) {
        degradation_.active_material_fraction = std::max(0.5,
            degradation_.active_material_fraction - d.active_material_loss_per_coulomb * std::abs(current) * dt);
    }
    if (d.enable_particle_cracking) {
        degradation_.crack_damage = std::clamp(degradation_.crack_damage + d.cracking_rate_per_gradient_second
                                               * std::abs(gradient) * dt, 0.0, 1.0);
    }
}

void SingleParticleModel::advance_thermal(double current, double voltage, double ocv, double dt) {
    if (external_thermal_) return;
    const auto& t = config_.thermal;
    const double q = std::abs(current) * std::abs(ocv - voltage);
    const double cooling = t.heat_transfer_coefficient_w_per_m2_k * t.cooling_area_m2
                         * (temperature_k_ - t.ambient_temperature_k);
    temperature_k_ += dt * (q - cooling) / (t.mass_kg * t.heat_capacity_j_per_kg_k);
}

BatteryStepResult SingleParticleModel::step(double current, double dt) {
    if (!(dt > 0.0) || !std::isfinite(current)) throw std::invalid_argument("invalid SPM step");
    const double jn = negative_flux(current), jp = positive_flux(current);
    negative_.step(jn, dt);
    positive_.step(jp, dt);
    time_s_ += dt;
    auto result = evaluate(current, 0.0, config_.electrolyte.initial_concentration_mol_per_m3,
                           config_.electrolyte.initial_concentration_mol_per_m3);
    const double gradient = negative_.surface_concentration(jn) - negative_.average_concentration();
    advance_degradation(current, dt, gradient / config_.negative.particle.maximum_concentration_mol_per_m3);
    advance_thermal(current, result.terminal_voltage_v, result.open_circuit_voltage_v, dt);
    result.temperature_k = temperature_k_;
    return result;
}

BatteryStepResult SingleParticleModel::state(double current) const {
    return evaluate(current, 0.0, config_.electrolyte.initial_concentration_mol_per_m3,
                    config_.electrolyte.initial_concentration_mol_per_m3);
}

SingleParticleElectrolyteModel::SingleParticleElectrolyteModel(LithiumIonCellConfig config)
    : SingleParticleModel(std::move(config)) {
    if (config_.electrolyte.cells < 3U) throw std::invalid_argument("SPMe needs at least three electrolyte cells");
    reset();
}

BatteryStepResult SingleParticleElectrolyteModel::state(double current) const {
    if (!std::isfinite(current)) throw std::invalid_argument("invalid SPMe current");
    const auto [lo, hi] = std::minmax_element(electrolyte_concentration_.begin(), electrolyte_concentration_.end());
    return evaluate(current, electrolyte_overpotential(current), *lo, *hi);
}

void SingleParticleElectrolyteModel::reset() {
    SingleParticleModel::reset();
    electrolyte_concentration_.assign(config_.electrolyte.cells, config_.electrolyte.initial_concentration_mol_per_m3);
    electrolyte_potential_v_.assign(config_.electrolyte.cells, 0.0);
    solid_potential_v_.assign(config_.electrolyte.cells, 0.0);
}

void SingleParticleElectrolyteModel::step_electrolyte(double current, double dt) {
    const auto& e = config_.electrolyte;
    const std::size_t n = electrolyte_concentration_.size();
    const double length = e.negative_thickness_m + e.separator_thickness_m + e.positive_thickness_m;
    const double dx = length / static_cast<double>(n);
    std::vector<double> lower(n - 1U, 0.0), diag(n, 1.0), upper(n - 1U, 0.0), rhs = electrolyte_concentration_;
    std::vector<double> porosity(n), deff(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        porosity[i] = region_porosity(e, x);
        deff[i] = effective_transport(e.diffusivity_m2_per_s, porosity[i], e.bruggeman_exponent);
    }
    for (std::size_t i = 0U; i + 1U < n; ++i) {
        const double df = 2.0 * deff[i] * deff[i + 1U] / std::max(deff[i] + deff[i + 1U], 1.0e-30);
        const double ki = df / (porosity[i] * dx * dx);
        const double kj = df / (porosity[i + 1U] * dx * dx);
        diag[i] += dt * ki; upper[i] -= dt * ki;
        lower[i] -= dt * kj; diag[i + 1U] += dt * kj;
    }
    const double source_n = (1.0 - e.transference_number) * current
                          / (faraday_constant * config_.area_m2 * e.negative_thickness_m);
    const double source_p = -(1.0 - e.transference_number) * current
                          / (faraday_constant * config_.area_m2 * e.positive_thickness_m);
    for (std::size_t i = 0U; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        double source = 0.0;
        if (x < e.negative_thickness_m) source = source_n;
        else if (x > e.negative_thickness_m + e.separator_thickness_m) source = source_p;
        rhs[i] += dt * source / porosity[i];
    }
    solve_tridiagonal(lower, diag, upper, rhs);
    for (double& c : rhs) c = std::max(c, 1.0e-6);
    electrolyte_concentration_ = std::move(rhs);
}

double SingleParticleElectrolyteModel::electrolyte_overpotential(double current) const {
    const auto& e = config_.electrolyte;
    const double length = e.negative_thickness_m + e.separator_thickness_m + e.positive_thickness_m;
    const double dx = length / static_cast<double>(electrolyte_concentration_.size());
    double r = 0.0;
    for (std::size_t i = 0U; i < electrolyte_concentration_.size(); ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        const double eps = region_porosity(e, x);
        const double kappa = effective_transport(e.conductivity_s_per_m, eps, e.bruggeman_exponent);
        r += dx / (kappa * config_.area_m2);
    }
    const std::size_t q = std::max<std::size_t>(1U, electrolyte_concentration_.size() / 8U);
    double cn = 0.0, cp = 0.0;
    for (std::size_t i = 0; i < q; ++i) cn += electrolyte_concentration_[i];
    for (std::size_t i = electrolyte_concentration_.size() - q; i < electrolyte_concentration_.size(); ++i) cp += electrolyte_concentration_[i];
    cn /= static_cast<double>(q); cp /= static_cast<double>(q);
    const double concentration = 2.0 * gas_constant * temperature_k_ * (1.0 - e.transference_number)
                               / faraday_constant * std::log(std::max(cp, 1.0e-9) / std::max(cn, 1.0e-9));
    return current * r - concentration;
}

void SingleParticleElectrolyteModel::update_potentials(double current) {
    const auto& e = config_.electrolyte;
    const std::size_t n = electrolyte_concentration_.size();
    const double length = e.negative_thickness_m + e.separator_thickness_m + e.positive_thickness_m;
    const double dx = length / static_cast<double>(n);
    electrolyte_potential_v_[0] = 0.0;
    solid_potential_v_[0] = config_.negative.open_circuit_voltage_v(negative_.surface_stoichiometry(negative_flux(current)));
    for (std::size_t i = 1U; i < n; ++i) {
        const double x = (static_cast<double>(i) - 0.5) * dx;
        const double eps = region_porosity(e, x);
        const double kappa = effective_transport(e.conductivity_s_per_m, eps, e.bruggeman_exponent);
        electrolyte_potential_v_[i] = electrolyte_potential_v_[i - 1U] - current * dx / (kappa * config_.area_m2);
        const double xpos = (static_cast<double>(i) + 0.5) * dx;
        double sigma = 1.0e30;
        if (xpos < e.negative_thickness_m) sigma = config_.negative.solid_conductivity_s_per_m;
        else if (xpos > e.negative_thickness_m + e.separator_thickness_m) sigma = config_.positive.solid_conductivity_s_per_m;
        solid_potential_v_[i] = solid_potential_v_[i - 1U] - (sigma > 1.0e20 ? 0.0 : current * dx / (sigma * config_.area_m2));
    }
}

BatteryStepResult SingleParticleElectrolyteModel::step(double current, double dt) {
    if (!(dt > 0.0) || !std::isfinite(current)) throw std::invalid_argument("invalid SPMe step");
    const double jn = negative_flux(current), jp = positive_flux(current);
    negative_.step(jn, dt); positive_.step(jp, dt); step_electrolyte(current, dt);
    time_s_ += dt;
    update_potentials(current);
    const auto [lo, hi] = std::minmax_element(electrolyte_concentration_.begin(), electrolyte_concentration_.end());
    const double eta_e = electrolyte_overpotential(current);
    auto result = evaluate(current, eta_e, *lo, *hi);
    const double gradient = negative_.surface_concentration(jn) - negative_.average_concentration();
    advance_degradation(current, dt, gradient / config_.negative.particle.maximum_concentration_mol_per_m3);
    advance_thermal(current, result.terminal_voltage_v, result.open_circuit_voltage_v, dt);
    result.temperature_k = temperature_k_;
    return result;
}

template<class Model>
DriveCycleResult drive_cycle_impl(Model& model, std::span<const DriveCyclePoint> segments, double max_step) {
    if (!(max_step > 0.0)) throw std::invalid_argument("drive-cycle maximum step must be positive");
    DriveCycleResult out;
    for (const auto& seg : segments) {
        if (!(seg.duration_s >= 0.0) || !std::isfinite(seg.current_a)) throw std::invalid_argument("invalid drive-cycle segment");
        double remaining = seg.duration_s;
        while (remaining > 1.0e-12) {
            const double dt = std::min(max_step, remaining);
            const auto s = model.step(seg.current_a, dt);
            out.samples.push_back(s);
            if (seg.current_a > 0.0) out.discharged_capacity_ah += seg.current_a * dt / 3600.0;
            out.electrical_energy_wh += s.terminal_voltage_v * seg.current_a * dt / 3600.0;
            remaining -= dt;
        }
    }
    return out;
}

DriveCycleResult simulate_drive_cycle(SingleParticleModel& model, std::span<const DriveCyclePoint> s, double dt) {
    return drive_cycle_impl(model, s, dt);
}
DriveCycleResult simulate_drive_cycle(SingleParticleElectrolyteModel& model, std::span<const DriveCyclePoint> s, double dt) {
    return drive_cycle_impl(model, s, dt);
}

std::complex<double> battery_impedance(double frequency, const BatteryEisConfig& c) {
    if (!(frequency > 0.0)) throw std::invalid_argument("battery EIS frequency must be positive");
    using C = std::complex<double>;
    const double omega = 2.0 * std::numbers::pi * frequency;
    const C j{0.0, 1.0};
    const auto parallel_rc = [&](double r, double cap) { return r / (C{1.0, 0.0} + j * omega * r * cap); };
    const C root = std::sqrt(j * omega * c.diffusion_time_s);
    const C warburg = c.diffusion_warburg_ohm_sqrt_s / std::sqrt(j * omega)
                    * (std::abs(root) > 1.0e-12 ? std::tanh(root) / root : C{1.0, 0.0});
    return C{c.series_resistance_ohm, 0.0}
         + parallel_rc(c.negative_charge_transfer_resistance_ohm, c.negative_double_layer_capacitance_f)
         + parallel_rc(c.positive_charge_transfer_resistance_ohm, c.positive_double_layer_capacitance_f)
         + warburg;
}

std::vector<std::complex<double>> battery_impedance_spectrum(std::span<const double> f, const BatteryEisConfig& c) {
    std::vector<std::complex<double>> z; z.reserve(f.size());
    for (double x : f) z.push_back(battery_impedance(x, c));
    return z;
}

} // namespace cfd::battery
