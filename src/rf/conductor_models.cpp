#include "cfd/rf/conductor_models.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::rf {

GroundReflection finite_ground_reflection(double frequency, double conductivity, double relative_permittivity,
                                          double angle) {
    if (!(frequency > 0.0) || conductivity < 0.0 || !(relative_permittivity > 0.0) || !std::isfinite(frequency) ||
        !std::isfinite(conductivity) || !std::isfinite(relative_permittivity) || !std::isfinite(angle) || angle < 0.0 ||
        angle >= 0.5 * std::numbers::pi) {
        throw std::invalid_argument("invalid finite-ground controls");
    }
    constexpr double epsilon0 = 8.8541878128e-12;
    constexpr double mu0 = 4.0 * std::numbers::pi * 1.0e-7;
    const double omega = 2.0 * std::numbers::pi * frequency;
    const std::complex<double> eps_r{relative_permittivity, -conductivity / (omega * epsilon0)};
    const std::complex<double> root = std::sqrt(eps_r);
    const double sin_i = std::sin(angle);
    const double cos_i = std::cos(angle);
    const std::complex<double> cos_t = std::sqrt(std::complex<double>{1.0 - sin_i * sin_i, 0.0} / eps_r);
    const std::complex<double> eta_ratio =
        std::sqrt(std::complex<double>{mu0, 0.0} / (std::complex<double>{mu0, 0.0} * eps_r));
    (void)eta_ratio;
    const std::complex<double> te = (cos_i - root * cos_t) / (cos_i + root * cos_t);
    const std::complex<double> tm = (root * cos_i - cos_t) / (root * cos_i + cos_t);
    return {te, tm};
}

double round_wire_proximity_factor(double spacing, double radius, double frequency, double conductivity,
                                   double relative_permeability) {
    if (!(spacing > 2.0 * radius) || !(radius > 0.0) || !(frequency > 0.0) || !(conductivity > 0.0) ||
        !(relative_permeability > 0.0)) {
        throw std::invalid_argument("invalid round-wire proximity controls");
    }
    const double delta = skin_depth_m(conductivity, frequency, relative_permeability);
    const double gap = spacing - 2.0 * radius;
    const double crowding = std::min(10.0, radius / delta);
    return 1.0 + 0.5 * std::pow(radius / gap, 2.0) * crowding;
}

std::vector<FilamentSegment> subdivide_rectangular_peec(const FilamentSegment& centerline, double width, double height,
                                                        std::size_t ny, std::size_t nz) {
    if (!(filament_length(centerline) > 0.0) || !(width > 0.0) || !(height > 0.0) || ny == 0U || nz == 0U ||
        !(centerline.conductivity_s_per_m > 0.0)) {
        throw std::invalid_argument("invalid rectangular PEEC subdivision");
    }
    std::vector<FilamentSegment> result;
    result.reserve(ny * nz);
    const double dy = width / static_cast<double>(ny);
    const double dz = height / static_cast<double>(nz);
    for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t iz = 0; iz < nz; ++iz) {
            const double y = (static_cast<double>(iy) + 0.5) * dy - 0.5 * width;
            const double z = (static_cast<double>(iz) + 0.5) * dz - 0.5 * height;
            result.push_back({{centerline.start.x, centerline.start.y + y, centerline.start.z + z},
                              {centerline.end.x, centerline.end.y + y, centerline.end.z + z},
                              std::min(dy, dz) * 0.5,
                              centerline.conductivity_s_per_m});
        }
    }
    return result;
}

} // namespace cfd::rf
