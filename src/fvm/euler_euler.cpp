#include "cfd/solvers/fvm/euler_euler.hpp"

#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] bool finite(Vec3 v) noexcept {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

[[nodiscard]] Vec3 componentwise_scale(Vec3 v, double scale) noexcept {
    return v * scale;
}

} // namespace

EulerEulerTransport::EulerEulerTransport(PolyMesh mesh, EulerEulerConfig config)
    : mesh_(std::move(mesh)), config_(config), dispersed_fraction_(mesh_.cell_count(), 0.0),
      primary_velocity_(mesh_.cell_count()), dispersed_velocity_(mesh_.cell_count()),
      primary_face_flux_(mesh_.face_count(), 0.0), dispersed_face_flux_(mesh_.face_count(), 0.0),
      drag_force_(mesh_.cell_count()) {
    validate_config();
    initialize(0.0);
}

void EulerEulerTransport::validate_config() const {
    if (!(config_.dt > 0.0) || !std::isfinite(config_.dt) || !(config_.primary_density > 0.0) ||
        !std::isfinite(config_.primary_density) || !(config_.dispersed_density > 0.0) ||
        !std::isfinite(config_.dispersed_density) || !(config_.drag_coefficient >= 0.0) ||
        !std::isfinite(config_.drag_coefficient) || !(config_.minimum_phase_fraction > 0.0) ||
        !(config_.minimum_phase_fraction < 0.5) || !std::isfinite(config_.minimum_phase_fraction) ||
        !finite(config_.primary_body_acceleration) || !finite(config_.dispersed_body_acceleration)) {
        throw std::invalid_argument("invalid Euler-Euler controls");
    }
}

void EulerEulerTransport::validate_velocity(Vec3 velocity) {
    if (!finite(velocity))
        throw std::invalid_argument("Euler-Euler velocity must be finite");
}

void EulerEulerTransport::initialize(double dispersed_fraction, Vec3 primary_velocity, Vec3 dispersed_velocity) {
    if (!std::isfinite(dispersed_fraction) || dispersed_fraction < 0.0 || dispersed_fraction > 1.0)
        throw std::invalid_argument("Euler-Euler volume fraction must lie in [0,1]");
    validate_velocity(primary_velocity);
    validate_velocity(dispersed_velocity);
    std::fill(dispersed_fraction_.begin(), dispersed_fraction_.end(), dispersed_fraction);
    std::fill(primary_velocity_.begin(), primary_velocity_.end(), primary_velocity);
    std::fill(dispersed_velocity_.begin(), dispersed_velocity_.end(), dispersed_velocity);
    std::fill(primary_face_flux_.begin(), primary_face_flux_.end(), 0.0);
    std::fill(dispersed_face_flux_.begin(), dispersed_face_flux_.end(), 0.0);
    std::fill(drag_force_.begin(), drag_force_.end(), Vec3{});
    time_ = 0.0;
    steps_ = 0U;
}

void EulerEulerTransport::initialize(std::span<const double> dispersed_fraction, std::span<const Vec3> primary_velocity,
                                     std::span<const Vec3> dispersed_velocity) {
    const auto n = mesh_.cell_count();
    if (dispersed_fraction.size() != n || primary_velocity.size() != n || dispersed_velocity.size() != n)
        throw std::invalid_argument("Euler-Euler initial field size mismatch");
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (!std::isfinite(dispersed_fraction[cell]) || dispersed_fraction[cell] < 0.0 ||
            dispersed_fraction[cell] > 1.0)
            throw std::invalid_argument("Euler-Euler volume fraction must lie in [0,1]");
        validate_velocity(primary_velocity[cell]);
        validate_velocity(dispersed_velocity[cell]);
    }
    std::copy(dispersed_fraction.begin(), dispersed_fraction.end(), dispersed_fraction_.begin());
    std::copy(primary_velocity.begin(), primary_velocity.end(), primary_velocity_.begin());
    std::copy(dispersed_velocity.begin(), dispersed_velocity.end(), dispersed_velocity_.begin());
    std::fill(primary_face_flux_.begin(), primary_face_flux_.end(), 0.0);
    std::fill(dispersed_face_flux_.begin(), dispersed_face_flux_.end(), 0.0);
    std::fill(drag_force_.begin(), drag_force_.end(), Vec3{});
    time_ = 0.0;
    steps_ = 0U;
}

void EulerEulerTransport::step() {
    const auto n = mesh_.cell_count();
    std::vector<Vec3> primary_momentum_flux(mesh_.face_count());
    std::vector<Vec3> dispersed_momentum_flux(mesh_.face_count());
    std::vector<Vec3> primary_momentum(n), dispersed_momentum(n);

    for (std::size_t cell = 0; cell < n; ++cell) {
        const double alpha = dispersed_fraction_[cell];
        const double beta = 1.0 - alpha;
        primary_momentum[cell] =
            componentwise_scale(primary_velocity_[cell], config_.primary_density * beta * mesh_.cells()[cell].volume);
        dispersed_momentum[cell] = componentwise_scale(dispersed_velocity_[cell],
                                                       config_.dispersed_density * alpha * mesh_.cells()[cell].volume);
        const Vec3 relative = dispersed_velocity_[cell] - primary_velocity_[cell];
        const double phase_weight = alpha * beta;
        double coefficient = config_.drag_coefficient * phase_weight;
        if (config_.quadratic_drag)
            coefficient *= magnitude(relative);
        drag_force_[cell] = relative * coefficient;
    }

    std::fill(primary_face_flux_.begin(), primary_face_flux_.end(), 0.0);
    std::fill(dispersed_face_flux_.begin(), dispersed_face_flux_.end(), 0.0);
    for (std::size_t face_index = 0; face_index < mesh_.face_count(); ++face_index) {
        const auto& face = mesh_.faces()[face_index];
        if (face.boundary())
            continue; // impermeable baseline boundary
        const auto owner = face.owner;
        const auto neighbour = face.neighbour;
        const Vec3 primary_face_velocity = (primary_velocity_[owner] + primary_velocity_[neighbour]) * 0.5;
        const Vec3 dispersed_face_velocity = (dispersed_velocity_[owner] + dispersed_velocity_[neighbour]) * 0.5;
        const double primary_normal_flux = dot(primary_face_velocity, face.area);
        const double dispersed_normal_flux = dot(dispersed_face_velocity, face.area);
        const bool primary_from_owner = primary_normal_flux >= 0.0;
        const bool dispersed_from_owner = dispersed_normal_flux >= 0.0;
        const double beta_owner = 1.0 - dispersed_fraction_[owner];
        const double beta_neighbour = 1.0 - dispersed_fraction_[neighbour];
        const double alpha_face = dispersed_from_owner ? dispersed_fraction_[owner] : dispersed_fraction_[neighbour];
        const double beta_face = primary_from_owner ? beta_owner : beta_neighbour;
        const double primary_flux = beta_face * primary_normal_flux;
        const double dispersed_flux = alpha_face * dispersed_normal_flux;
        primary_face_flux_[face_index] = primary_flux;
        dispersed_face_flux_[face_index] = dispersed_flux;
        const Vec3& primary_upwind = primary_from_owner ? primary_velocity_[owner] : primary_velocity_[neighbour];
        const Vec3& dispersed_upwind =
            dispersed_from_owner ? dispersed_velocity_[owner] : dispersed_velocity_[neighbour];
        primary_momentum_flux[face_index] = primary_upwind * (config_.primary_density * primary_flux);
        dispersed_momentum_flux[face_index] = dispersed_upwind * (config_.dispersed_density * dispersed_flux);
    }

    std::vector<double> next_fraction = dispersed_fraction_;
    std::vector<Vec3> next_primary = primary_velocity_;
    std::vector<Vec3> next_dispersed = dispersed_velocity_;
    for (std::size_t cell = 0; cell < n; ++cell) {
        double dispersed_net_flux = 0.0;
        Vec3 primary_net_momentum{};
        Vec3 dispersed_net_momentum{};
        for (const auto face_index : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[face_index];
            const double orientation = face.owner == cell ? 1.0 : -1.0;
            dispersed_net_flux += orientation * dispersed_face_flux_[face_index];
            primary_net_momentum += primary_momentum_flux[face_index] * orientation;
            dispersed_net_momentum += dispersed_momentum_flux[face_index] * orientation;
        }
        const double volume = mesh_.cells()[cell].volume;
        next_fraction[cell] =
            std::clamp(dispersed_fraction_[cell] - config_.dt * dispersed_net_flux / volume, 0.0, 1.0);
        const double alpha = dispersed_fraction_[cell];
        const double beta = 1.0 - alpha;
        const double primary_mass = config_.primary_density * beta * volume;
        const double dispersed_mass = config_.dispersed_density * alpha * volume;
        const Vec3 primary_force = config_.primary_body_acceleration * (primary_mass) + drag_force_[cell] * volume;
        const Vec3 dispersed_force = config_.dispersed_body_acceleration * (dispersed_mass)-drag_force_[cell] * volume;
        primary_momentum[cell] -= primary_net_momentum * config_.dt;
        dispersed_momentum[cell] -= dispersed_net_momentum * config_.dt;
        primary_momentum[cell] += primary_force * config_.dt;
        dispersed_momentum[cell] += dispersed_force * config_.dt;
        const double next_beta = 1.0 - next_fraction[cell];
        const double next_primary_mass = config_.primary_density * next_beta * volume;
        const double next_dispersed_mass = config_.dispersed_density * next_fraction[cell] * volume;
        if (next_primary_mass > config_.minimum_phase_fraction * config_.primary_density * volume)
            next_primary[cell] = primary_momentum[cell] / next_primary_mass;
        else
            next_primary[cell] = {};
        if (next_dispersed_mass > config_.minimum_phase_fraction * config_.dispersed_density * volume)
            next_dispersed[cell] = dispersed_momentum[cell] / next_dispersed_mass;
        else
            next_dispersed[cell] = {};
    }
    dispersed_fraction_.swap(next_fraction);
    primary_velocity_.swap(next_primary);
    dispersed_velocity_.swap(next_dispersed);
    time_ += config_.dt;
    ++steps_;
}

double EulerEulerTransport::volume_fraction_integral() const {
    return cfd::core::parallel_sum(
        mesh_.cell_count(), [&](std::size_t cell) { return dispersed_fraction_[cell] * mesh_.cells()[cell].volume; });
}

Vec3 EulerEulerTransport::total_momentum() const noexcept {
    Vec3 result{};
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell) {
        const double volume = mesh_.cells()[cell].volume;
        result += primary_velocity_[cell] * (config_.primary_density * (1.0 - dispersed_fraction_[cell]) * volume);
        result += dispersed_velocity_[cell] * (config_.dispersed_density * dispersed_fraction_[cell] * volume);
    }
    return result;
}

Vec3 EulerEulerTransport::integrated_drag_force() const noexcept {
    Vec3 result{};
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell)
        result += drag_force_[cell] * mesh_.cells()[cell].volume;
    return result;
}

double EulerEulerTransport::minimum_fraction() const {
    return *std::min_element(dispersed_fraction_.begin(), dispersed_fraction_.end());
}

double EulerEulerTransport::maximum_fraction() const {
    return *std::max_element(dispersed_fraction_.begin(), dispersed_fraction_.end());
}

} // namespace cfd::fvm
