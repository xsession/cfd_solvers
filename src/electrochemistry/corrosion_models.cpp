#include "cfd/electrochemistry/corrosion_models.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <numbers>
namespace cfd::electrochemistry {
double reaction_quotient(std::span<const RedoxActivityTerm> products, std::span<const RedoxActivityTerm> reactants) {
    double q = 1.0;
    auto apply = [&](auto terms, double sign) {
        for (auto t : terms) {
            if (!(t.activity > 0) || !std::isfinite(t.activity) || !std::isfinite(t.stoichiometric_power))
                throw std::invalid_argument("invalid redox activity");
            q *= std::pow(t.activity, sign * t.stoichiometric_power);
        }
    };
    apply(products, 1.0);
    apply(reactants, -1.0);
    if (!(q > 0) || !std::isfinite(q))
        throw std::overflow_error("reaction quotient out of range");
    return q;
}
double concentration_dependent_nernst(double e0, double t, double n, std::span<const RedoxActivityTerm> p,
                                      std::span<const RedoxActivityTerm> r) {
    return nernst_potential(e0, t, n, reaction_quotient(p, r));
}
namespace {
double cathodic_limit(double activation, double limit) {
    if (!(limit > 0) || !std::isfinite(limit))
        return activation;
    const double cath = std::max(0.0, -activation);
    const double limited = cath / (1.0 + cath / limit);
    return activation < 0 ? -limited : activation;
}
} // namespace
CathodicReactionResult oxygen_reduction_cathodic(double e, double e0, double oxygen, double proton, double i0,
                                                 double ilim, double t, double n) {
    const RedoxActivityTerm react[] = {{oxygen, 1.0}, {proton, 4.0}};
    const RedoxActivityTerm prod[] = {{1.0, 2.0}};
    const double eq = concentration_dependent_nernst(e0, t, n, prod, react);
    const double activation = butler_volmer_current_density(i0, e - eq, t, n);
    return {eq, activation, cathodic_limit(activation, ilim)};
}
CathodicReactionResult hydrogen_evolution_cathodic(double e, double h2, double proton, double i0, double ilim, double t,
                                                   double n) {
    const RedoxActivityTerm react[] = {{proton, 2.0}};
    const RedoxActivityTerm prod[] = {{h2, 1.0}};
    const double eq = concentration_dependent_nernst(0.0, t, n, prod, react);
    const double activation = butler_volmer_current_density(i0, e - eq, t, n);
    return {eq, activation, cathodic_limit(activation, ilim)};
}
double charge_transfer_resistance(double i0, double t, double n, double aa, double ac) {
    if (!(i0 > 0) || !(t > 0) || !(n > 0) || !(aa + ac > 0))
        throw std::invalid_argument("invalid charge-transfer controls");
    return cfd::chemistry::gas_constant * t / (i0 * n * faraday_constant * (aa + ac));
}
std::complex<double> randles_impedance(double f, double rs, double rct, double cdl, double sigma) {
    if (f < 0 || rs < 0 || !(rct > 0) || cdl < 0 || sigma < 0)
        throw std::invalid_argument("invalid EIS parameters");
    const double w = 2.0 * std::numbers::pi * f;
    std::complex<double> adm{1.0 / rct, w * cdl};
    if (sigma > 0 && w > 0) {
        const std::complex<double> zw = sigma / std::sqrt(w) * std::complex<double>{1, -1};
        adm += 1.0 / zw;
    }
    return rs + 1.0 / adm;
}
void DoubleLayerState::advance_current(double current, double dt) {
    if (!(capacitance_per_area > 0) || !(dt >= 0) || !std::isfinite(current))
        throw std::invalid_argument("invalid double-layer state");
    potential += current * dt / capacitance_per_area;
}

double passivation_current_density(double potential, const PassivationModel& model, double temperature) {
    if (!(model.exchange_current_density >= 0.0) || !(model.passive_current_density >= 0.0) ||
        !(model.transpassive_slope >= 0.0) || !(model.electrons > 0.0) || !(temperature > 0.0) ||
        !(model.transpassive_onset >= model.passive_onset))
        throw std::invalid_argument("invalid passivation model");
    if (potential < model.passive_onset)
        return butler_volmer_current_density(model.exchange_current_density,
                                             potential - model.active_equilibrium_potential, temperature,
                                             model.electrons);
    if (potential <= model.transpassive_onset)
        return model.passive_current_density;
    return model.passive_current_density + model.transpassive_slope * (potential - model.transpassive_onset);
}
void ProductLayerState::advance_from_anodic_current(double current, double dt, double molar_mass, double density,
                                                    double electrons) {
    if (thickness < 0.0 || resistivity < 0.0 || growth_fraction < 0.0 || growth_fraction > 1.0 || dt < 0.0)
        throw std::invalid_argument("invalid product-layer state");
    thickness +=
        growth_fraction * faradaic_recession_distance(std::max(0.0, current), dt, molar_mass, density, electrons);
}
double faradaic_recession_distance(double current, double dt, double molar_mass, double density, double electrons) {
    if (dt < 0.0 || !std::isfinite(dt))
        throw std::invalid_argument("invalid recession time");
    return corrosion_penetration_rate(std::abs(current), molar_mass, density, electrons) * dt;
}
double cathodic_protection_current_density(double protected_potential, double free_potential,
                                           double polarization_resistance) {
    if (!(polarization_resistance > 0.0) || !std::isfinite(protected_potential) || !std::isfinite(free_potential))
        throw std::invalid_argument("invalid cathodic-protection controls");
    return (protected_potential - free_potential) / polarization_resistance;
}

double bruggeman_effective_transport(double bulk, double porosity, double exponent) {
    if (!(bulk >= 0.0) || !(porosity >= 0.0 && porosity <= 1.0) || !(exponent >= 0.0))
        throw std::invalid_argument("invalid Bruggeman transport controls");
    return bulk * std::pow(porosity, exponent);
}
double solid_phase_ohmic_drop(double current, double length, double conductivity) {
    if (!(length >= 0.0) || !(conductivity > 0.0) || !std::isfinite(current))
        throw std::invalid_argument("invalid solid conduction controls");
    return current * length / conductivity;
}
double marcus_current_density(double i0, double eta, double lambda, double temperature, double electrons) {
    if (!(i0 >= 0.0) || !(lambda > 0.0) || !(temperature > 0.0) || !(electrons > 0.0))
        throw std::invalid_argument("invalid Marcus kinetics controls");
    const double work = electrons * faraday_constant * eta;
    const double denom = 4.0 * lambda * cfd::chemistry::gas_constant * temperature;
    const double forward = std::exp(-((lambda - work) * (lambda - work) - lambda * lambda) / denom);
    const double reverse = std::exp(-((lambda + work) * (lambda + work) - lambda * lambda) / denom);
    return i0 * (forward - reverse);
}
double bikerman_activity_correction(double concentration, double occupied_volume, double total) {
    if (!(concentration >= 0.0) || !(occupied_volume >= 0.0) || !(total >= 0.0))
        throw std::invalid_argument("invalid finite-size PNP state");
    const double free_fraction = 1.0 - occupied_volume * total;
    if (!(free_fraction > 0.0))
        throw std::domain_error("finite-size PNP packing fraction >= 1");
    return concentration / free_fraction;
}
ModifiedPnpFlux modified_nernst_planck_flux(double left, double right, double diffusivity, double potential_gradient,
                                            int charge, double temperature, double width, double occupied, double total,
                                            double velocity) {
    if (!(left >= 0.0) || !(right >= 0.0) || !(diffusivity >= 0.0) || !(temperature > 0.0) || !(width > 0.0) ||
        !(occupied >= 0.0) || !(total >= 0.0) || !std::isfinite(potential_gradient) || !std::isfinite(velocity))
        throw std::invalid_argument("invalid modified PNP flux controls");
    const double free_fraction = 1.0 - occupied * total;
    if (!(free_fraction > 0.0))
        throw std::domain_error("modified PNP packing fraction >= 1");
    const double concentration = 0.5 * (left + right), gradient = (right - left) / width;
    const double steric_gradient = gradient / free_fraction;
    const double electric_gradient = static_cast<double>(charge) * faraday_constant /
                                     (cfd::chemistry::gas_constant * temperature) * concentration * potential_gradient;
    return {velocity * concentration - diffusivity * (steric_gradient + electric_gradient),
            concentration / free_fraction};
}
ImmersedInterfaceWeight smooth_electrode_interface(double phi, double thickness) {
    if (!(thickness > 0.0))
        throw std::invalid_argument("invalid immersed-interface thickness");
    const double q = phi / thickness, t = std::tanh(q), c = std::cosh(q);
    return {0.5 * (1.0 + t), 0.5 / (thickness * c * c)};
}
CutCellInterface cut_cell_electrode_interface(double left, double right, double width) {
    if (!(width > 0.0) || !std::isfinite(left) || !std::isfinite(right))
        throw std::invalid_argument("invalid cut-cell electrode geometry");
    const bool left_electrode = left < 0.0, right_electrode = right < 0.0;
    if (left_electrode == right_electrode)
        return left_electrode ? CutCellInterface{0.0, 1.0, 0.0} : CutCellInterface{1.0, 0.0, 0.0};
    const double fraction = std::clamp(left / (left - right), 0.0, 1.0);
    const double electrode = left_electrode ? fraction : 1.0 - fraction;
    return {1.0 - electrode, electrode, 1.0};
}
CoupledElectrodeResult solve_anodic_bv_mass_transfer(double bulk, double km, double e0, double electrode, double i0,
                                                     double temperature, double electrons, double cref) {
    if (!(bulk > 0.0) || !(km > 0.0) || !(i0 >= 0.0) || !(temperature > 0.0) || !(electrons > 0.0) || !(cref > 0.0))
        throw std::invalid_argument("invalid coupled electrode controls");
    auto residual = [&](double cs) {
        const double eq =
            e0 + (cfd::chemistry::gas_constant * temperature / (electrons * faraday_constant)) * std::log(cs / cref);
        const double current = butler_volmer_current_density(i0, electrode - eq, temperature, electrons);
        return km * (cs - bulk) - current / (electrons * faraday_constant);
    };
    double lo = bulk, hi = std::max(2.0 * bulk, bulk + 1.0);
    double flo = residual(lo), fhi = residual(hi);
    for (int grow = 0; grow < 80 && flo * fhi > 0.0; ++grow) {
        hi *= 2.0;
        fhi = residual(hi);
    }
    if (flo * fhi > 0.0)
        throw std::runtime_error("could not bracket coupled electrode concentration");
    std::size_t it = 0;
    for (; it < 100; ++it) {
        const double mid = 0.5 * (lo + hi), fm = residual(mid);
        if (std::abs(fm) < 1e-12 * std::max(1.0, km * bulk)) {
            lo = hi = mid;
            break;
        }
        if (flo * fm <= 0.0) {
            hi = mid;
            fhi = fm;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    const double cs = 0.5 * (lo + hi);
    const double eq =
        e0 + (cfd::chemistry::gas_constant * temperature / (electrons * faraday_constant)) * std::log(cs / cref);
    const double current = butler_volmer_current_density(i0, electrode - eq, temperature, electrons);
    return {cs, current, current / (electrons * faraday_constant), it + 1};
}
double acid_speciation_redox_potential(double e0, double temperature, double electrons, double proton_power,
                                       std::span<const cfd::chemistry::AcidFamily> families, double strong_charge) {
    if (!(proton_power >= 0.0))
        throw std::invalid_argument("invalid proton stoichiometry");
    const auto eq = cfd::chemistry::equilibrate_acids(families, strong_charge);
    const double q = std::pow(eq.hydrogen_mol_per_litre, -proton_power);
    return nernst_potential(e0, temperature, electrons, q);
}
bool pitting_initiates(double chloride, double potential, const PittingCriterion& c) {
    if (!(c.critical_chloride >= 0.0) || !std::isfinite(c.critical_potential))
        throw std::invalid_argument("invalid pitting criterion");
    return chloride >= c.critical_chloride && potential >= c.critical_potential;
}
double stress_assisted_exchange_current(double base, double stress, double volume, double temperature) {
    if (!(base >= 0.0) || !(temperature > 0.0) || !std::isfinite(stress) || !std::isfinite(volume))
        throw std::invalid_argument("invalid stress-corrosion state");
    return base * std::exp(stress * volume / (cfd::chemistry::gas_constant * temperature));
}
double phase_field_corrosion_step(double phi, double laplacian, double driving, double mobility, double gradient,
                                  double dt) {
    if (!(mobility >= 0.0) || !(gradient >= 0.0) || !(dt >= 0.0))
        throw std::invalid_argument("invalid phase-field controls");
    const double derivative = phi * (phi * phi - 1.0) - gradient * laplacian - driving;
    return std::clamp(phi - mobility * dt * derivative, -1.0, 1.0);
}
double recessed_level_set(double phi, double distance) {
    if (!(distance >= 0.0))
        throw std::invalid_argument("invalid recession distance");
    return phi + distance;
}

double marcus_hush_chidsey_current_density(double i0, double eta, double lambda, double temperature, double electrons,
                                           std::size_t points) {
    if (!(i0 >= 0) || !(lambda > 0) || !(temperature > 0) || !(electrons > 0) || points < 100)
        throw std::invalid_argument("invalid MHC kinetics controls");
    const double lam = lambda / (cfd::chemistry::gas_constant * temperature),
                 v = electrons * faraday_constant * eta / (cfd::chemistry::gas_constant * temperature);
    const double lo = -20, hi = 20, h = (hi - lo) / static_cast<double>(points);
    auto integral = [&](double shift) {
        double q = 0;
        for (std::size_t k = 0; k < points; ++k) {
            double x = lo + (static_cast<double>(k) + 0.5) * h;
            double fermi = 1.0 / (1.0 + std::exp(x));
            q += std::exp(-(x - lam + shift) * (x - lam + shift) / (4 * lam)) * fermi;
        }
        return q * h;
    };
    double forward = integral(-v), reverse = integral(v);
    double norm = integral(0.0);
    return i0 * (forward - reverse) / std::max(norm, 1e-300);
}
std::vector<double> recession_remesh_1d(std::span<const double> coordinates, double recession) {
    if (coordinates.size() < 2U || !(recession >= 0.0) || !std::isfinite(recession))
        throw std::invalid_argument("invalid recession remesh controls");
    std::vector<double> result(coordinates.begin(), coordinates.end());
    for (std::size_t i = 1; i < result.size(); ++i) {
        if (!(coordinates[i] > coordinates[i - 1U]) || !std::isfinite(coordinates[i]))
            throw std::invalid_argument("remesh coordinates must be finite and increasing");
    }
    const double length = coordinates.back() - coordinates.front();
    if (!(recession < length))
        throw std::invalid_argument("recession consumes the mesh interval");
    const double scale = (length - recession) / length;
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = coordinates.front() + scale * (coordinates[i] - coordinates.front());
    return result;
}
StressCorrosionUpdate coupled_stress_corrosion_step(double base, double stress, double volume, double temperature,
                                                    double dt, double mass, double density, double electrons) {
    if (!(base >= 0.0) || !std::isfinite(stress) || !std::isfinite(volume) || !(temperature > 0.0) || !(dt >= 0.0) ||
        !(mass > 0.0) || !(density > 0.0) || !(electrons > 0.0))
        throw std::invalid_argument("invalid stress-corrosion coupling controls");
    const double current = stress_assisted_exchange_current(base, stress, volume, temperature);
    return {current, faradaic_recession_distance(current, dt, mass, density, electrons)};
}
PhaseFieldCorrosion2D::PhaseFieldCorrosion2D(std::size_t nx, std::size_t ny, double dx, double dy, double mobility,
                                             double gradient)
    : nx_(nx), ny_(ny), dx_(dx), dy_(dy), mobility_(mobility), gradient_(gradient), phi_(nx * ny), next_(nx * ny) {
    if (nx < 3U || ny < 3U || !(dx > 0.0) || !(dy > 0.0) || !(mobility >= 0.0) || !(gradient >= 0.0))
        throw std::invalid_argument("invalid 2-D phase-field grid");
}
void PhaseFieldCorrosion2D::initialize_interface(double x0, double y0, double width) {
    if (!(width > 0.0) || !std::isfinite(x0) || !std::isfinite(y0))
        throw std::invalid_argument("invalid 2-D phase-field interface");
    for (std::size_t j = 0; j < ny_; ++j)
        for (std::size_t i = 0; i < nx_; ++i) {
            const double x = static_cast<double>(i) * dx_, y = static_cast<double>(j) * dy_;
            const double d = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
            phi_[j * nx_ + i] = std::tanh((width - d) / width);
        }
}
void PhaseFieldCorrosion2D::step(double dt, std::span<const double> driving) {
    if (!(dt >= 0.0) || driving.size() != phi_.size())
        throw std::invalid_argument("2-D phase-field step mismatch");
    for (std::size_t j = 0; j < ny_; ++j)
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t k = j * nx_ + i, il = j * nx_ + (i ? i - 1 : i), ir = j * nx_ + (i + 1 < nx_ ? i + 1 : i),
                              ib = (j ? j - 1 : j) * nx_ + i, it = (j + 1 < ny_ ? j + 1 : j) * nx_ + i;
            const double lap =
                (phi_[il] - 2 * phi_[k] + phi_[ir]) / (dx_ * dx_) + (phi_[ib] - 2 * phi_[k] + phi_[it]) / (dy_ * dy_);
            next_[k] = phase_field_corrosion_step(phi_[k], lap, driving[k], mobility_, gradient_, dt);
        }
    phi_.swap(next_);
}
PhaseFieldCorrosion1D::PhaseFieldCorrosion1D(std::size_t n, double dx, double m, double g)
    : cells_(n), dx_(dx), mobility_(m), gradient_(g), phi_(n), next_(n) {
    if (n < 3 || !(dx > 0) || !(m >= 0) || !(g >= 0))
        throw std::invalid_argument("invalid phase-field grid");
}
void PhaseFieldCorrosion1D::initialize_interface(double x0, double w) {
    if (!(w > 0))
        throw std::invalid_argument("invalid phase-field width");
    for (std::size_t i = 0; i < cells_; ++i)
        phi_[i] = std::tanh((x0 - static_cast<double>(i) * dx_) / w);
}
void PhaseFieldCorrosion1D::step(double dt, std::span<const double> driving) {
    if (!(dt >= 0) || driving.size() != cells_)
        throw std::invalid_argument("phase-field step mismatch");
    for (std::size_t i = 0; i < cells_; ++i) {
        std::size_t l = i ? i - 1 : i, r = i + 1 < cells_ ? i + 1 : i;
        double lap = (phi_[l] - 2 * phi_[i] + phi_[r]) / (dx_ * dx_);
        next_[i] = phase_field_corrosion_step(phi_[i], lap, driving[i], mobility_, gradient_, dt);
    }
    phi_.swap(next_);
}

} // namespace cfd::electrochemistry
