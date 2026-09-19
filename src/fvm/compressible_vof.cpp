#include "cfd/solvers/fvm/compressible_vof.hpp"

#include "cfd/core/parallel.hpp"
#include "cfd/fvm/schemes.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::fvm {
namespace {

[[nodiscard]] bool finite(double value) noexcept {
    return std::isfinite(value);
}

} // namespace

CompressibleVofTransport::CompressibleVofTransport(PolyMesh mesh, CompressibleVofConfig config)
    : mesh_(std::move(mesh)), config_(config), liquid_fraction_(mesh_.cell_count(), 0.0),
      mass_(mesh_.cell_count(), 0.0), density_(mesh_.cell_count(), config_.gas_density),
      pressure_(mesh_.cell_count(), config_.reference_pressure),
      liquid_density_(mesh_.cell_count(), config_.liquid_density),
      gas_density_(mesh_.cell_count(), config_.gas_density), face_flux_(mesh_.face_count(), 0.0) {
    validate_config();
    initialize(0.0, config_.reference_pressure);
}

void CompressibleVofTransport::validate_config() const {
    if (!(config_.dt > 0.0) || !finite(config_.dt) || !(config_.liquid_density > 0.0) ||
        !finite(config_.liquid_density) || !(config_.gas_density > 0.0) || !finite(config_.gas_density) ||
        !(config_.liquid_bulk_modulus > 0.0) || !finite(config_.liquid_bulk_modulus) || !(config_.gas_gamma > 1.0) ||
        !finite(config_.gas_gamma) || !(config_.reference_pressure > 0.0) || !finite(config_.reference_pressure) ||
        !(config_.minimum_pressure > 0.0) || !finite(config_.minimum_pressure) || !(config_.minimum_density > 0.0) ||
        !finite(config_.minimum_density) || config_.minimum_pressure > config_.reference_pressure) {
        throw std::invalid_argument("invalid compressible VOF controls");
    }
}

double CompressibleVofTransport::mixture_density(double liquid_fraction, double pressure) const {
    const double liquid =
        config_.liquid_density * (1.0 + (pressure - config_.reference_pressure) / config_.liquid_bulk_modulus);
    const double gas = config_.gas_density * std::pow(pressure / config_.reference_pressure, 1.0 / config_.gas_gamma);
    return liquid_fraction * liquid + (1.0 - liquid_fraction) * gas;
}

double CompressibleVofTransport::pressure_for_density(double liquid_fraction, double density) const {
    if (!(liquid_fraction >= 0.0 && liquid_fraction <= 1.0) || !(density >= config_.minimum_density) ||
        !finite(density))
        throw std::runtime_error("invalid compressible VOF mixture state");
    double lower = config_.minimum_pressure;
    const double lower_density = mixture_density(liquid_fraction, lower);
    if (density <= lower_density)
        return lower;
    double upper = std::max(config_.reference_pressure * 2.0, lower * 2.0);
    for (std::size_t attempt = 0; attempt < 80U && mixture_density(liquid_fraction, upper) < density; ++attempt)
        upper *= 2.0;
    if (!(mixture_density(liquid_fraction, upper) >= density) || !finite(upper))
        throw std::runtime_error("compressible VOF pressure bracket failed");
    for (std::size_t iteration = 0; iteration < 80U; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        if (mixture_density(liquid_fraction, middle) < density)
            lower = middle;
        else
            upper = middle;
    }
    return 0.5 * (lower + upper);
}

void CompressibleVofTransport::initialize(double liquid_fraction, double pressure) {
    if (!finite(liquid_fraction) || liquid_fraction < 0.0 || liquid_fraction > 1.0)
        throw std::invalid_argument("compressible VOF fraction must lie in [0,1]");
    if (!(pressure >= config_.minimum_pressure) || !finite(pressure))
        throw std::invalid_argument("compressible VOF pressure must be positive");
    std::fill(liquid_fraction_.begin(), liquid_fraction_.end(), liquid_fraction);
    std::fill(pressure_.begin(), pressure_.end(), pressure);
    std::fill(face_flux_.begin(), face_flux_.end(), 0.0);
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell)
        mass_[cell] = mixture_density(liquid_fraction, pressure) * mesh_.cells()[cell].volume;
    update_properties();
    time_ = 0.0;
    steps_ = 0U;
}

void CompressibleVofTransport::initialize(std::span<const double> liquid_fraction, std::span<const double> pressure) {
    if (liquid_fraction.size() != mesh_.cell_count() || pressure.size() != mesh_.cell_count())
        throw std::invalid_argument("compressible VOF initial field size mismatch");
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell) {
        if (!finite(liquid_fraction[cell]) || liquid_fraction[cell] < 0.0 || liquid_fraction[cell] > 1.0)
            throw std::invalid_argument("compressible VOF fraction must lie in [0,1]");
        if (!(pressure[cell] >= config_.minimum_pressure) || !finite(pressure[cell]))
            throw std::invalid_argument("compressible VOF pressure must be positive");
        mass_[cell] = mixture_density(liquid_fraction[cell], pressure[cell]) * mesh_.cells()[cell].volume;
    }
    std::copy(liquid_fraction.begin(), liquid_fraction.end(), liquid_fraction_.begin());
    std::copy(pressure.begin(), pressure.end(), pressure_.begin());
    std::fill(face_flux_.begin(), face_flux_.end(), 0.0);
    update_properties();
    time_ = 0.0;
    steps_ = 0U;
}

void CompressibleVofTransport::set_face_flux(std::vector<double> volumetric_flux) {
    if (volumetric_flux.size() != mesh_.face_count())
        throw std::invalid_argument("compressible VOF face-flux size mismatch");
    for (double flux : volumetric_flux)
        if (!finite(flux))
            throw std::invalid_argument("compressible VOF face flux must be finite");
    face_flux_ = std::move(volumetric_flux);
}

void CompressibleVofTransport::update_properties() {
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell) {
        density_[cell] = mass_[cell] / mesh_.cells()[cell].volume;
        if (!(density_[cell] >= config_.minimum_density) || !finite(density_[cell]))
            throw std::runtime_error("compressible VOF update produced non-positive density");
        pressure_[cell] = pressure_for_density(liquid_fraction_[cell], density_[cell]);
        liquid_density_[cell] = config_.liquid_density *
                                (1.0 + (pressure_[cell] - config_.reference_pressure) / config_.liquid_bulk_modulus);
        gas_density_[cell] =
            config_.gas_density * std::pow(pressure_[cell] / config_.reference_pressure, 1.0 / config_.gas_gamma);
    }
}

void CompressibleVofTransport::step() {
    const auto face_fraction =
        interpolate_scalar_to_faces(mesh_, liquid_fraction_, FaceInterpolationScheme::upwind, face_flux_);
    const auto face_density = interpolate_scalar_to_faces(mesh_, density_, FaceInterpolationScheme::upwind, face_flux_);
    std::vector<double> next_fraction = liquid_fraction_;
    std::vector<double> next_mass = mass_;
    for (std::size_t cell = 0; cell < mesh_.cell_count(); ++cell) {
        double volume_flux = 0.0;
        double mass_flux = 0.0;
        for (const auto face_index : mesh_.cell_faces()[cell]) {
            const auto& face = mesh_.faces()[face_index];
            const double orientation = face.owner == cell ? 1.0 : -1.0;
            volume_flux += orientation * face_flux_[face_index] * face_fraction[face_index];
            mass_flux += orientation * face_flux_[face_index] * face_density[face_index];
        }
        const double volume = mesh_.cells()[cell].volume;
        next_fraction[cell] = std::clamp(liquid_fraction_[cell] - config_.dt * volume_flux / volume, 0.0, 1.0);
        next_mass[cell] = mass_[cell] - config_.dt * mass_flux;
        if (!(next_mass[cell] >= config_.minimum_density * volume) || !finite(next_mass[cell]))
            throw std::runtime_error("compressible VOF update produced non-positive mass");
    }
    liquid_fraction_.swap(next_fraction);
    mass_.swap(next_mass);
    update_properties();
    time_ += config_.dt;
    ++steps_;
}

double CompressibleVofTransport::liquid_volume() const {
    return cfd::core::parallel_sum(
        mesh_.cell_count(), [&](std::size_t cell) { return liquid_fraction_[cell] * mesh_.cells()[cell].volume; });
}

double CompressibleVofTransport::total_mass() const {
    return cfd::core::parallel_sum(mass_.size(), [&](std::size_t cell) { return mass_[cell]; });
}

double CompressibleVofTransport::minimum_pressure() const {
    return *std::min_element(pressure_.begin(), pressure_.end());
}

double CompressibleVofTransport::maximum_pressure() const {
    return *std::max_element(pressure_.begin(), pressure_.end());
}

} // namespace cfd::fvm
