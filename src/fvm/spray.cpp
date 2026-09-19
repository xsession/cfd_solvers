#include "cfd/solvers/fvm/spray.hpp"

#include "cfd/solvers/fvm/particles.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] bool finite(Vec3 value) noexcept {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

[[nodiscard]] std::size_t nearest_cell(const PolyMesh& mesh, Vec3 position) {
    std::size_t nearest = 0U;
    double best_distance = std::numeric_limits<double>::max();
    for (std::size_t cell = 0; cell < mesh.cell_count(); ++cell) {
        const Vec3 delta = mesh.cells()[cell].center - position;
        const double distance = dot(delta, delta);
        if (distance < best_distance) {
            best_distance = distance;
            nearest = cell;
        }
    }
    return nearest;
}

} // namespace

SprayInjectionEvaporation::SprayInjectionEvaporation(PolyMesh mesh, SprayConfig config)
    : mesh_(std::move(mesh)), config_(config), vapor_mass_source_(mesh_.cell_count(), 0.0),
      carrier_energy_source_(mesh_.cell_count(), 0.0), carrier_momentum_source_(mesh_.cell_count()) {
    validate_config();
}

void SprayInjectionEvaporation::validate_config() const {
    if (!(config_.dt > 0.0) || !std::isfinite(config_.dt) || !(config_.liquid_density > 0.0) ||
        !std::isfinite(config_.liquid_density) || !(config_.carrier_dynamic_viscosity > 0.0) ||
        !std::isfinite(config_.carrier_dynamic_viscosity) || !(config_.evaporation_constant >= 0.0) ||
        !std::isfinite(config_.evaporation_constant) || !(config_.latent_heat >= 0.0) ||
        !std::isfinite(config_.latent_heat) || !(config_.thermal_relaxation_time > 0.0) ||
        !std::isfinite(config_.thermal_relaxation_time) || !(config_.minimum_diameter > 0.0) ||
        !std::isfinite(config_.minimum_diameter) || config_.maximum_parcels == 0U) {
        throw std::invalid_argument("invalid spray controls");
    }
}

void SprayInjectionEvaporation::validate_vec(Vec3 value, const char* message) {
    if (!finite(value))
        throw std::invalid_argument(message);
}

void SprayInjectionEvaporation::validate_parcel(const SprayParcel& parcel) const {
    validate_vec(parcel.position, "spray parcel position must be finite");
    validate_vec(parcel.velocity, "spray parcel velocity must be finite");
    if (!(parcel.diameter >= config_.minimum_diameter) || !std::isfinite(parcel.diameter) || !(parcel.mass > 0.0) ||
        !std::isfinite(parcel.mass) || !std::isfinite(parcel.temperature)) {
        throw std::invalid_argument("invalid spray parcel");
    }
}

void SprayInjectionEvaporation::validate_injection(const SprayInjection& injection) const {
    validate_vec(injection.position, "spray injection position must be finite");
    validate_vec(injection.velocity, "spray injection velocity must be finite");
    if (!(injection.diameter >= config_.minimum_diameter) || !std::isfinite(injection.diameter) ||
        !(injection.temperature >= 0.0) || !std::isfinite(injection.temperature) ||
        !(injection.mass_flow_rate >= 0.0) || !std::isfinite(injection.mass_flow_rate) ||
        !(injection.start_time >= 0.0) || !std::isfinite(injection.start_time) || !(injection.duration > 0.0) ||
        !std::isfinite(injection.duration)) {
        throw std::invalid_argument("invalid spray injection");
    }
}

void SprayInjectionEvaporation::set_injection(SprayInjection injection) {
    validate_injection(injection);
    injection_ = injection;
}

void SprayInjectionEvaporation::add_parcel(SprayParcel parcel) {
    if (parcels_.size() >= config_.maximum_parcels)
        throw std::runtime_error("spray parcel capacity exceeded");
    validate_parcel(parcel);
    parcels_.push_back(parcel);
}

void SprayInjectionEvaporation::step(std::span<const Vec3> carrier_velocity,
                                     std::span<const double> carrier_temperature) {
    if (carrier_velocity.size() != mesh_.cell_count() || carrier_temperature.size() != mesh_.cell_count())
        throw std::invalid_argument("spray carrier field size mismatch");
    for (const auto& velocity : carrier_velocity)
        validate_vec(velocity, "spray carrier velocity must be finite");
    for (double temperature : carrier_temperature)
        if (!(temperature >= 0.0) || !std::isfinite(temperature))
            throw std::invalid_argument("spray carrier temperature must be non-negative and finite");

    std::fill(vapor_mass_source_.begin(), vapor_mass_source_.end(), 0.0);
    std::fill(carrier_energy_source_.begin(), carrier_energy_source_.end(), 0.0);
    std::fill(carrier_momentum_source_.begin(), carrier_momentum_source_.end(), Vec3{});

    if (injection_.has_value() && time_ >= injection_->start_time &&
        time_ < injection_->start_time + injection_->duration && injection_->mass_flow_rate > 0.0) {
        if (parcels_.size() >= config_.maximum_parcels)
            throw std::runtime_error("spray parcel capacity exceeded");
        add_parcel({injection_->position, injection_->velocity, injection_->diameter,
                    injection_->mass_flow_rate * config_.dt, injection_->temperature});
        injected_mass_ += injection_->mass_flow_rate * config_.dt;
    }

    std::vector<SprayParcel> active;
    active.reserve(parcels_.size());
    for (auto parcel : parcels_) {
        const std::size_t cell = nearest_cell(mesh_, parcel.position);
        const double old_mass = parcel.mass;
        const double tau =
            config_.liquid_density * parcel.diameter * parcel.diameter / (18.0 * config_.carrier_dynamic_viscosity);
        const double velocity_decay = std::exp(-config_.dt / tau);
        const Vec3 old_velocity = parcel.velocity;
        parcel.velocity = carrier_velocity[cell] + (parcel.velocity - carrier_velocity[cell]) * velocity_decay;
        const Vec3 drag_on_parcel = (parcel.velocity - old_velocity) * (old_mass / config_.dt);
        const double volume = mesh_.cells()[cell].volume;
        carrier_momentum_source_[cell] -= drag_on_parcel / volume;

        const double thermal_decay = std::exp(-config_.dt / config_.thermal_relaxation_time);
        parcel.temperature =
            carrier_temperature[cell] + (parcel.temperature - carrier_temperature[cell]) * thermal_decay;
        const double new_diameter = evaporate_d2_law(parcel.diameter, config_.evaporation_constant, config_.dt);
        const double diameter_ratio = new_diameter / parcel.diameter;
        parcel.mass = old_mass * diameter_ratio * diameter_ratio * diameter_ratio;
        const double vaporized = old_mass - parcel.mass;
        if (vaporized > 0.0) {
            vapor_mass_source_[cell] += vaporized / (config_.dt * volume);
            carrier_energy_source_[cell] -= config_.latent_heat * vaporized / (config_.dt * volume);
            evaporated_mass_ += vaporized;
        }
        parcel.diameter = new_diameter;
        parcel.position += parcel.velocity * config_.dt;
        if (parcel.diameter > config_.minimum_diameter && parcel.mass > 0.0)
            active.push_back(parcel);
    }
    parcels_.swap(active);
    time_ += config_.dt;
    ++steps_;
}

void SprayInjectionEvaporation::step(std::span<const Vec3> carrier_velocity, double carrier_temperature) {
    if (!(carrier_temperature >= 0.0) || !std::isfinite(carrier_temperature))
        throw std::invalid_argument("spray carrier temperature must be non-negative and finite");
    std::vector<double> temperature(mesh_.cell_count(), carrier_temperature);
    step(carrier_velocity, std::span<const double>(temperature));
}

double SprayInjectionEvaporation::total_liquid_mass() const noexcept {
    double total = 0.0;
    for (const auto& parcel : parcels_)
        total += parcel.mass;
    return total;
}

} // namespace cfd::fvm
