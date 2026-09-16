#include "cfd/solvers/electrochemistry/nernst_planck1d.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace cfd::electrochemistry {

NernstPlanck1D::NernstPlanck1D(NernstPlanck1DConfig config) : config_(config) {
    if (config_.cells < 3U || !(config_.length > 0.0) || !(config_.temperature > 0.0)
        || !(config_.dt > 0.0)) {
        throw std::invalid_argument("invalid Nernst-Planck 1-D configuration");
    }
    dx_ = config_.length / static_cast<double>(config_.cells);
    potential_.assign(config_.cells, 0.0);
    flux_.assign(config_.cells + 1U, 0.0);
    next_.assign(config_.cells, 0.0);
}

void NernstPlanck1D::set_potential(const std::function<double(double)>& fn) {
    for (std::size_t i = 0; i < config_.cells; ++i) {
        potential_[i] = fn((static_cast<double>(i) + 0.5) * dx_);
    }
}

void NernstPlanck1D::set_potential(std::span<const double> values) {
    if (values.size() != config_.cells) throw std::invalid_argument("potential size mismatch");
    potential_.assign(values.begin(), values.end());
}

std::size_t NernstPlanck1D::add_species(std::string name,
                                        int charge,
                                        double diffusivity,
                                        const std::function<double(double)>& initial_concentration) {
    if (!(diffusivity > 0.0)) throw std::invalid_argument("species diffusivity must be positive");
    NernstPlanckSpecies1D species;
    species.name = std::move(name);
    species.charge = charge;
    species.diffusivity = diffusivity;
    species.concentration.resize(config_.cells);
    for (std::size_t i = 0; i < config_.cells; ++i) {
        const double value = initial_concentration((static_cast<double>(i) + 0.5) * dx_);
        if (value < 0.0 || !std::isfinite(value)) throw std::invalid_argument("invalid initial concentration");
        species.concentration[i] = value;
    }
    species_.push_back(std::move(species));
    return species_.size() - 1U;
}

std::size_t NernstPlanck1D::add_species(std::string name,
                                        int charge,
                                        double diffusivity,
                                        double uniform_concentration) {
    if (uniform_concentration < 0.0) throw std::invalid_argument("concentration must be non-negative");
    return add_species(std::move(name), charge, diffusivity,
                       [uniform_concentration](double) { return uniform_concentration; });
}

void NernstPlanck1D::step_species(NernstPlanckSpecies1D& species) {
    const std::size_t n = config_.cells;
    const double migration_scale = static_cast<double>(species.charge) * faraday_constant
                                 / (cfd::chemistry::gas_constant * config_.temperature);
    std::fill(flux_.begin(), flux_.end(), 0.0);

    const auto face_flux = [&](std::size_t left, std::size_t right) {
        const double c_left = species.concentration[left];
        const double c_right = species.concentration[right];
        const double c_face = 0.5 * (c_left + c_right);
        const double dc_dx = (c_right - c_left) / dx_;
        const double dphi_dx = (potential_[right] - potential_[left]) / dx_;
        const double advective = config_.advection_velocity >= 0.0
                               ? config_.advection_velocity * c_left
                               : config_.advection_velocity * c_right;
        return advective - species.diffusivity * (dc_dx + migration_scale * c_face * dphi_dx);
    };

    if (config_.boundary == TransportBoundary1D::periodic) {
        flux_[0] = face_flux(n - 1U, 0U);
        for (std::size_t face = 1U; face < n; ++face) {
            flux_[face] = face_flux(face - 1U, face);
        }
        flux_[n] = flux_[0];
    } else {
        flux_[0] = 0.0;
        for (std::size_t face = 1U; face < n; ++face) {
            flux_[face] = face_flux(face - 1U, face);
        }
        flux_[n] = 0.0;
    }

    for (std::size_t i = 0; i < n; ++i) {
        next_[i] = species.concentration[i] - (config_.dt / dx_) * (flux_[i + 1U] - flux_[i]);
        if (!std::isfinite(next_[i]) || next_[i] < -1.0e-12) {
            throw std::runtime_error("Nernst-Planck step became unstable or negative; reduce dt");
        }
        if (next_[i] < 0.0) next_[i] = 0.0;
    }
    species.concentration.swap(next_);
}

void NernstPlanck1D::step(std::size_t count) {
    if (species_.empty()) throw std::runtime_error("Nernst-Planck solver has no species");
    for (std::size_t step_index = 0U; step_index < count; ++step_index) {
        for (auto& species : species_) step_species(species);
        time_ += config_.dt;
    }
}

double NernstPlanck1D::total_amount(std::size_t species_index) const {
    if (species_index >= species_.size()) throw std::out_of_range("species index out of range");
    return std::accumulate(species_[species_index].concentration.begin(),
                           species_[species_index].concentration.end(), 0.0) * dx_;
}

double NernstPlanck1D::minimum_concentration(std::size_t species_index) const {
    if (species_index >= species_.size()) throw std::out_of_range("species index out of range");
    return *std::min_element(species_[species_index].concentration.begin(),
                             species_[species_index].concentration.end());
}

} // namespace cfd::electrochemistry
