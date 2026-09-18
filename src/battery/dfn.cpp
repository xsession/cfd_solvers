#include "cfd/battery/dfn.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

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

} // namespace

DfnConfig default_graphite_nmc_dfn_config() {
    const LithiumIonCellConfig base = default_graphite_nmc_config();
    DfnConfig c;
    c.area_m2 = base.area_m2;
    c.negative_thickness_m = base.negative.thickness_m;
    c.positive_thickness_m = base.positive.thickness_m;
    c.negative_particle = base.negative.particle;
    c.positive_particle = base.positive.particle;
    c.negative_solid_conductivity_s_per_m = base.negative.solid_conductivity_s_per_m;
    c.positive_solid_conductivity_s_per_m = base.positive.solid_conductivity_s_per_m;
    c.negative_exchange_current_density_a_per_m2 = base.negative.exchange_current_density_a_per_m2;
    c.positive_exchange_current_density_a_per_m2 = base.positive.exchange_current_density_a_per_m2;
    c.negative_porosity = base.electrolyte.negative_porosity;
    c.positive_porosity = base.electrolyte.positive_porosity;
    c.separator_porosity = base.electrolyte.separator_porosity;
    c.bruggeman_exponent = base.electrolyte.bruggeman_exponent;
    c.negative_active_volume_fraction = base.negative.active_volume_fraction;
    c.positive_active_volume_fraction = base.positive.active_volume_fraction;
    c.electrolyte_initial_concentration_mol_per_m3 = base.electrolyte.initial_concentration_mol_per_m3;
    c.electrolyte_diffusivity_m2_per_s = base.electrolyte.diffusivity_m2_per_s;
    c.electrolyte_conductivity_s_per_m = base.electrolyte.conductivity_s_per_m;
    c.transference_number = base.electrolyte.transference_number;
    c.negative_ocv_v = graphite_ocv_v;
    c.positive_ocv_v = nmc_ocv_v;
    c.negative_nodes = 20U;
    c.separator_nodes = 9U;
    c.positive_nodes = 18U;
    c.contact_resistance_ohm = base.contact_resistance_ohm;
    c.initial_temperature_k = base.initial_temperature_k;
    return c;
}

DoyleFullerNewmanModel::DoyleFullerNewmanModel(DfnConfig config) : config_(std::move(config)) {
    if (!(config_.area_m2 > 0.0)) throw std::invalid_argument("invalid DFN area");
    if (!(config_.negative_thickness_m > 0.0) || !(config_.separator_thickness_m > 0.0)
        || !(config_.positive_thickness_m > 0.0)) throw std::invalid_argument("invalid DFN thickness");
    if (config_.negative_nodes < 2U || config_.separator_nodes < 1U || config_.positive_nodes < 2U) {
        throw std::invalid_argument("DFN needs at least two electrode nodes per side");
    }
    if (!config_.negative_ocv_v) config_.negative_ocv_v = graphite_ocv_v;
    if (!config_.positive_ocv_v) config_.positive_ocv_v = nmc_ocv_v;
    reset();
}

std::size_t DoyleFullerNewmanModel::total_nodes() const noexcept {
    return config_.negative_nodes + config_.separator_nodes + config_.positive_nodes;
}
double DoyleFullerNewmanModel::cell_length_m() const noexcept {
    return config_.negative_thickness_m + config_.separator_thickness_m + config_.positive_thickness_m;
}
double DoyleFullerNewmanModel::node_spacing_m() const noexcept {
    return cell_length_m() / static_cast<double>(total_nodes());
}

void DoyleFullerNewmanModel::reset() {
    negative_particles_.assign(config_.negative_nodes, SphericalDiffusionParticle{config_.negative_particle});
    positive_particles_.assign(config_.positive_nodes, SphericalDiffusionParticle{config_.positive_particle});
    electrolyte_concentration_.assign(total_nodes(), config_.electrolyte_initial_concentration_mol_per_m3);
    electrolyte_potential_v_.assign(total_nodes(), 0.0);
    solid_potential_v_.assign(total_nodes(), 0.0);
    local_reaction_current_.assign(total_nodes(), 0.0);
    temperature_k_ = config_.initial_temperature_k;
    time_s_ = 0.0;
}

double DoyleFullerNewmanModel::region_porosity(double x) const {
    if (x < config_.negative_thickness_m) return config_.negative_porosity;
    if (x < config_.negative_thickness_m + config_.separator_thickness_m) return config_.separator_porosity;
    return config_.positive_porosity;
}

double DoyleFullerNewmanModel::solid_conductivity(double x) const {
    if (x < config_.negative_thickness_m) return config_.negative_solid_conductivity_s_per_m;
    if (x < config_.negative_thickness_m + config_.separator_thickness_m) return 0.0;
    return config_.positive_solid_conductivity_s_per_m;
}

void DoyleFullerNewmanModel::update_reaction_partition(double current) {
    // Distribute the full electrode current across the active nodes of each
    // electrode, weighted by local exchange-current activity. Each electrode
    // carries the full terminal current (SPMe convention), which keeps the global
    // lithium balance closed.
    const std::size_t nn = config_.negative_nodes;
    const std::size_t ns = config_.separator_nodes;
    const double i_areal = current / config_.area_m2;
    for (std::size_t i = 0U; i < total_nodes(); ++i) local_reaction_current_[i] = 0.0;

    auto partition = [&](std::span<SphericalDiffusionParticle> particles,
                         double i0_base, double sign, std::size_t offset) {
        double weight_sum = 0.0;
        std::vector<double> w(particles.size(), 0.0);
        for (std::size_t i = 0U; i < particles.size(); ++i) {
            const double x = clamp01(particles[i].surface_stoichiometry());
            w[i] = i0_base * 2.0 * std::sqrt(x * (1.0 - x));
            weight_sum += w[i];
        }
        if (weight_sum <= 0.0) weight_sum = 1.0;
        for (std::size_t i = 0U; i < particles.size(); ++i) {
            local_reaction_current_[offset + i] = sign * i_areal * w[i] / weight_sum;
        }
    };
    partition(std::span<SphericalDiffusionParticle>{negative_particles_},
              config_.negative_exchange_current_density_a_per_m2, +1.0, 0U);
    partition(std::span<SphericalDiffusionParticle>{positive_particles_},
              config_.positive_exchange_current_density_a_per_m2, -1.0, nn + ns);
}

void DoyleFullerNewmanModel::advance_particles(double /*current*/, double dt) {
    const double dx = node_spacing_m();
    const double a_neg = 3.0 * config_.negative_active_volume_fraction / config_.negative_particle.radius_m;
    const double a_pos = 3.0 * config_.positive_active_volume_fraction / config_.positive_particle.radius_m;
    const std::size_t nn = config_.negative_nodes;
    const std::size_t ns = config_.separator_nodes;
    for (std::size_t i = 0U; i < config_.negative_nodes; ++i) {
        const double flux = local_reaction_current_[i] / (a_neg * dx * faraday_constant);
        negative_particles_[i].step(flux, dt);
    }
    for (std::size_t i = 0U; i < config_.positive_nodes; ++i) {
        const double flux = local_reaction_current_[nn + ns + i] / (a_pos * dx * faraday_constant);
        positive_particles_[i].step(flux, dt);
    }
}

void DoyleFullerNewmanModel::advance_electrolyte_concentration(double /*current*/, double dt) {
    const std::size_t n = total_nodes();
    const double dx = node_spacing_m();
    std::vector<double> porosity(n), deff(n);
    for (std::size_t i = 0U; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        porosity[i] = region_porosity(x);
        deff[i] = effective_transport(config_.electrolyte_diffusivity_m2_per_s, porosity[i], config_.bruggeman_exponent);
    }
    // Harmonic-mean face conductances make the scheme conservative for the
    // lithium mole count (porosity * concentration), not for concentration alone.
    std::vector<double> face(n + 1U, 0.0);
    for (std::size_t i = 1U; i < n; ++i) {
        const double df = 2.0 * deff[i - 1U] * deff[i] / std::max(deff[i - 1U] + deff[i], 1.0e-30);
        face[i] = df / (dx * dx);
    }
    std::vector<double> lower(n - 1U, 0.0), diag(n, 1.0), upper(n - 1U, 0.0), rhs(n);
    for (std::size_t i = 0U; i < n; ++i) {
        const double left = face[i];
        const double right = face[i + 1U];
        diag[i] = 1.0 + dt / porosity[i] * (left + right);
        if (i > 0U) lower[i - 1U] = -dt / porosity[i] * left;
        if (i + 1U < n) upper[i] = -dt / porosity[i] * right;
        rhs[i] = electrolyte_concentration_[i];
        const double source = (1.0 - config_.transference_number) * local_reaction_current_[i] / (dx * faraday_constant);
        rhs[i] += dt * source / porosity[i];
    }
    solve_tridiagonal(lower, diag, upper, rhs);
    for (double& c : rhs) c = std::max(c, 1.0e-6);
    electrolyte_concentration_ = std::move(rhs);
}

void DoyleFullerNewmanModel::update_potentials(double current) {
    const std::size_t n = total_nodes();
    const double dx = node_spacing_m();
    const double i_areal = current / config_.area_m2;

    // Electrolyte current density is the running sum of the local areal reaction
    // current; it is zero at both current collectors and peaks in the separator.
    std::vector<double> i_e(n, 0.0);
    double running = 0.0;
    for (std::size_t i = 0U; i < n; ++i) {
        i_e[i] = running;
        running += local_reaction_current_[i];
    }
    electrolyte_potential_v_[0] = 0.0;
    for (std::size_t i = 1U; i < n; ++i) {
        const double x = (static_cast<double>(i) - 0.5) * dx;
        const double eps = region_porosity(x);
        const double kappa = effective_transport(config_.electrolyte_conductivity_s_per_m, eps, config_.bruggeman_exponent);
        const double i_face = 0.5 * (i_e[i - 1U] + i_e[i]);
        electrolyte_potential_v_[i] = electrolyte_potential_v_[i - 1U] - i_face * dx / (kappa * config_.area_m2);
    }

    // Solid potential: solid current is i_areal at the negative collector and
    // decays by the local reaction; it is flat in the separator.
    solid_potential_v_[0] = config_.negative_ocv_v(negative_particles_[0].surface_stoichiometry());
    double i_s = i_areal;
    for (std::size_t i = 1U; i < n; ++i) {
        const double x = (static_cast<double>(i) - 0.5) * dx;
        const double sigma = solid_conductivity(x);
        if (sigma > 0.0) {
            solid_potential_v_[i] = solid_potential_v_[i - 1U] - i_s * dx / (sigma * config_.area_m2);
            i_s -= local_reaction_current_[i - 1U];
        } else {
            solid_potential_v_[i] = solid_potential_v_[i - 1U];
        }
    }
}

DfnStepResult DoyleFullerNewmanModel::evaluate(double current) const {
    const std::size_t nn = config_.negative_nodes;
    const std::size_t ns = config_.separator_nodes;
    const std::size_t np = config_.positive_nodes;
    const double dx = node_spacing_m();
    const double a_neg = 3.0 * config_.negative_active_volume_fraction / config_.negative_particle.radius_m;
    const double a_pos = 3.0 * config_.positive_active_volume_fraction / config_.positive_particle.radius_m;

    // Area-averaged OCV and reaction overpotential from the distributed particles.
    double ocv_neg = 0.0, ocv_pos = 0.0, reaction = 0.0;
    for (std::size_t i = 0U; i < nn; ++i) {
        const double x = clamp01(negative_particles_[i].surface_stoichiometry());
        ocv_neg += config_.negative_ocv_v(x);
        const double j_surf = local_reaction_current_[i] / (a_neg * dx);
        const double i0 = config_.negative_exchange_current_density_a_per_m2 * 2.0 * std::sqrt(x * (1.0 - x));
        reaction += std::abs(symmetric_butler_volmer_overpotential_v(j_surf, std::max(i0, 1.0e-12), temperature_k_));
    }
    for (std::size_t i = 0U; i < np; ++i) {
        const double x = clamp01(positive_particles_[i].surface_stoichiometry());
        ocv_pos += config_.positive_ocv_v(x);
        const double j_surf = local_reaction_current_[nn + ns + i] / (a_pos * dx);
        const double i0 = config_.positive_exchange_current_density_a_per_m2 * 2.0 * std::sqrt(x * (1.0 - x));
        reaction += std::abs(symmetric_butler_volmer_overpotential_v(j_surf, std::max(i0, 1.0e-12), temperature_k_));
    }
    ocv_neg /= static_cast<double>(nn);
    ocv_pos /= static_cast<double>(np);
    reaction /= static_cast<double>(nn + np);
    const double ocv = ocv_pos - ocv_neg;

    // Distributed electrolyte and solid ohmic resistances.
    double r_e = 0.0;
    const std::size_t total = total_nodes();
    for (std::size_t i = 0U; i < total; ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        const double eps = region_porosity(x);
        const double kappa = effective_transport(config_.electrolyte_conductivity_s_per_m, eps, config_.bruggeman_exponent);
        r_e += dx / (kappa * config_.area_m2);
    }
    const double r_s = config_.negative_thickness_m / (config_.negative_solid_conductivity_s_per_m * config_.area_m2)
                     + config_.positive_thickness_m / (config_.positive_solid_conductivity_s_per_m * config_.area_m2);

    // Nernst (concentration) polarization from the electrode-end concentrations.
    const std::size_t q = std::max<std::size_t>(1U, total / 8U);
    double cn = 0.0, cp = 0.0;
    for (std::size_t i = 0U; i < q; ++i) cn += electrolyte_concentration_[i];
    for (std::size_t i = total - q; i < total; ++i) cp += electrolyte_concentration_[i];
    cn /= static_cast<double>(q); cp /= static_cast<double>(q);
    const double concentration = 2.0 * gas_constant * temperature_k_ * (1.0 - config_.transference_number)
                                / faraday_constant * std::log(std::max(cp, 1.0e-9) / std::max(cn, 1.0e-9));
    const double eta_e = current * r_e - concentration;
    const double ohmic = current * (config_.contact_resistance_ohm + r_s);
    const double sign = (current == 0.0) ? 1.0 : current;
    const double voltage = ocv - std::copysign(reaction + std::abs(eta_e), sign) - ohmic;
    const double heat = std::abs(current) * std::max(0.0, std::abs(ocv - voltage));

    const double neg_surf = negative_particles_.front().surface_stoichiometry();
    const double pos_surf = positive_particles_.back().surface_stoichiometry();
    const auto [lo, hi] = std::minmax_element(electrolyte_concentration_.begin(), electrolyte_concentration_.end());
    return {time_s_, current, voltage, ocv, reaction, ohmic + eta_e, neg_surf, pos_surf, *lo, *hi, temperature_k_, heat};
}

DfnStepResult DoyleFullerNewmanModel::step(double current, double dt) {
    if (!(dt > 0.0) || !std::isfinite(current)) throw std::invalid_argument("invalid DFN step");
    update_reaction_partition(current);
    advance_particles(current, dt);
    advance_electrolyte_concentration(current, dt);
    time_s_ += dt;
    update_potentials(current);
    auto result = evaluate(current);
    result.temperature_k = temperature_k_;
    return result;
}

DfnStepResult DoyleFullerNewmanModel::state(double current) const {
    if (!std::isfinite(current)) throw std::invalid_argument("invalid DFN current");
    return evaluate(current);
}

double DoyleFullerNewmanModel::total_lithium_mol() const noexcept {
    const double dx = node_spacing_m();
    const double a = config_.area_m2;
    double total = 0.0;
    for (const auto& p : negative_particles_) total += p.average_concentration() * config_.negative_active_volume_fraction * a * dx;
    for (const auto& p : positive_particles_) total += p.average_concentration() * config_.positive_active_volume_fraction * a * dx;
    for (std::size_t i = 0U; i < electrolyte_concentration_.size(); ++i) {
        const double x = (static_cast<double>(i) + 0.5) * dx;
        total += electrolyte_concentration_[i] * region_porosity(x) * a * dx;
    }
    return total;
}

} // namespace cfd::battery
