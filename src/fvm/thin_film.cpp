#include "cfd/solvers/fvm/thin_film.hpp"

#include "cfd/core/parallel.hpp"
#include "cfd/fvm/advanced_models.hpp"
#include "cfd/fvm/operators.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::fvm {
namespace {

void require_finite(double value, const char* message) {
    if (!std::isfinite(value))
        throw std::invalid_argument(message);
}

double face_distance(const PolyMesh& mesh, const Face& face) {
    const Vec3 delta = face.boundary() ? face.center - mesh.cells()[face.owner].center
                                       : mesh.cells()[face.neighbour].center - mesh.cells()[face.owner].center;
    const double area = magnitude(face.area);
    if (!(area > 0.0))
        throw std::runtime_error("thin-film face has zero area");
    const Vec3 normal = face.area / area;
    const double projected = std::abs(dot(delta, normal));
    const double distance = projected > 1.0e-14 ? projected : magnitude(delta);
    if (!(distance > 0.0))
        throw std::runtime_error("thin-film face has zero center distance");
    return distance;
}

} // namespace

ThinFilmTransport::ThinFilmTransport(PolyMesh mesh, ThinFilmConfig config)
    : mesh_(std::move(mesh)), config_(config), thickness_(mesh_.cell_count()), pressure_(mesh_.cell_count()),
      face_flux_(mesh_.face_count()) {
    if (!(config_.dt > 0.0) || !(config_.density > 0.0) || !(config_.viscosity > 0.0) ||
        config_.surface_tension < 0.0 || !(config_.minimum_thickness >= 0.0) || !std::isfinite(config_.dt) ||
        !std::isfinite(config_.density) || !std::isfinite(config_.viscosity) ||
        !std::isfinite(config_.surface_tension) || !std::isfinite(config_.minimum_thickness)) {
        throw std::invalid_argument("invalid thin-film controls");
    }
    require_finite(config_.body_acceleration.x, "invalid thin-film body acceleration");
    require_finite(config_.body_acceleration.y, "invalid thin-film body acceleration");
    require_finite(config_.body_acceleration.z, "invalid thin-film body acceleration");
    initialize(config_.minimum_thickness);
}

void ThinFilmTransport::initialize(double thickness) {
    if (!std::isfinite(thickness) || thickness < 0.0)
        throw std::invalid_argument("invalid thin-film thickness");
    std::fill(thickness_.begin(), thickness_.end(), std::max(config_.minimum_thickness, thickness));
    std::fill(pressure_.begin(), pressure_.end(), 0.0);
    std::fill(face_flux_.begin(), face_flux_.end(), 0.0);
}

void ThinFilmTransport::set_boundary_pressure(std::vector<double> pressure) {
    if (pressure.size() != mesh_.face_count())
        throw std::invalid_argument("thin-film boundary pressure size mismatch");
    for (const double value : pressure)
        require_finite(value, "invalid thin-film boundary pressure");
    boundary_pressure_ = std::move(pressure);
}

void ThinFilmTransport::update_pressure() {
    const auto gradient = least_squares_gradient_scalar(mesh_, thickness_);
    std::vector<Vec3> normal(thickness_.size());
    for (std::size_t cell = 0; cell < thickness_.size(); ++cell) {
        const double magnitude_gradient = magnitude(gradient[cell]);
        if (magnitude_gradient > 1.0e-14)
            normal[cell] = gradient[cell] / magnitude_gradient;
    }
    const auto curvature = gauss_divergence_vector(mesh_, normal);
    for (std::size_t cell = 0; cell < pressure_.size(); ++cell)
        pressure_[cell] = -config_.surface_tension * curvature[cell];
}

void ThinFilmTransport::update_face_flux() {
    std::fill(face_flux_.begin(), face_flux_.end(), 0.0);
    for (std::size_t face_index = 0; face_index < mesh_.face_count(); ++face_index) {
        const auto& face = mesh_.faces()[face_index];
        if (face.boundary() && boundary_pressure_.empty())
            continue;

        const double area = magnitude(face.area);
        const double distance = face_distance(mesh_, face);
        const std::size_t owner = face.owner;
        const double owner_pressure = pressure_[owner];
        const double owner_thickness = thickness_[owner];
        double neighbour_pressure = owner_pressure;
        double face_thickness = owner_thickness;
        if (face.boundary()) {
            neighbour_pressure = boundary_pressure_[face_index];
        } else {
            neighbour_pressure = pressure_[face.neighbour];
            face_thickness = 0.5 * (owner_thickness + thickness_[face.neighbour]);
        }

        const Vec3 normal = face.area / area;
        const double pressure_gradient = (neighbour_pressure - owner_pressure) / distance;
        const double body_force = config_.density * dot(config_.body_acceleration, normal);
        const double pressure_gradient_with_body_force = pressure_gradient - body_force;
        const double flux_per_area =
            lubrication_thin_film_flux(face_thickness, pressure_gradient_with_body_force, config_.viscosity);
        face_flux_[face_index] = flux_per_area * area;
    }
}

void ThinFilmTransport::step() {
    update_pressure();
    update_face_flux();
    std::vector<double> next = thickness_;
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell) {
        double net_flux = 0.0;
        for (const std::size_t face_index : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[face_index];
            net_flux += face.owner == cell ? face_flux_[face_index] : -face_flux_[face_index];
        }
        next[cell] =
            std::max(config_.minimum_thickness, thickness_[cell] - config_.dt * net_flux / mesh_.cells()[cell].volume);
    }
    thickness_.swap(next);
}

double ThinFilmTransport::inventory() const {
    return cfd::core::parallel_sum(thickness_.size(),
                                   [&](std::size_t cell) { return thickness_[cell] * mesh_.cells()[cell].volume; });
}

} // namespace cfd::fvm
