#include "cfd/solvers/fdtd/emc_sources.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fdtd {
namespace {
std::array<double, 3> cross(std::array<double, 3> a, std::array<double, 3> b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
} // namespace

std::vector<double> make_graded_axis(double length, std::size_t cells, double first, double ratio) {
    if (!(length > 0.0) || cells == 0U || !(first > 0.0) || !(ratio > 0.0) || !std::isfinite(length) ||
        !std::isfinite(first) || !std::isfinite(ratio)) {
        throw std::invalid_argument("invalid graded-axis controls");
    }
    const double sum = std::abs(ratio - 1.0) < 1.0e-12
                           ? static_cast<double>(cells)
                           : (std::pow(ratio, static_cast<double>(cells)) - 1.0) / (ratio - 1.0);
    const double scale = length / (first * sum);
    std::vector<double> coordinates(cells + 1U);
    for (std::size_t i = 0; i < cells; ++i) {
        const double width = first * scale * std::pow(ratio, static_cast<double>(i));
        coordinates[i + 1U] = coordinates[i] + width;
    }
    coordinates.back() = length;
    return coordinates;
}

std::vector<SubcellWeight> deposit_thin_wire_subcell(double x, double y, double z, std::size_t nx, std::size_t ny,
                                                     std::size_t nz, double dx, double dy, double dz) {
    if (nx == 0U || ny == 0U || nz == 0U || !(dx > 0.0) || !(dy > 0.0) || !(dz > 0.0) || !std::isfinite(x) ||
        !std::isfinite(y) || !std::isfinite(z)) {
        throw std::invalid_argument("invalid thin-wire subcell controls");
    }
    const auto axis = [](double coordinate, double spacing, std::size_t count) {
        const double q = std::clamp(coordinate / spacing, 0.0, static_cast<double>(count) - 1.0);
        const std::size_t lo = static_cast<std::size_t>(std::floor(q));
        const std::size_t hi = std::min(lo + 1U, count - 1U);
        return std::array<std::size_t, 2>{lo, hi};
    };
    const auto ix = axis(x, dx, nx), iy = axis(y, dy, ny), iz = axis(z, dz, nz);
    const double fx = std::clamp(x / dx - static_cast<double>(ix[0]), 0.0, 1.0);
    const double fy = std::clamp(y / dy - static_cast<double>(iy[0]), 0.0, 1.0);
    const double fz = std::clamp(z / dz - static_cast<double>(iz[0]), 0.0, 1.0);
    std::vector<SubcellWeight> result;
    for (std::size_t a = 0; a < 2U; ++a)
        for (std::size_t b = 0; b < 2U; ++b)
            for (std::size_t c = 0; c < 2U; ++c) {
                const double weight = (a ? fx : 1.0 - fx) * (b ? fy : 1.0 - fy) * (c ? fz : 1.0 - fz);
                if (weight > 0.0)
                    result.push_back({ix[a], iy[b], iz[c], weight});
            }
    return result;
}

HuygensEquivalentSource huygens_equivalent_source(std::array<double, 3> normal, std::array<double, 3> electric,
                                                  std::array<double, 3> magnetic) {
    const double norm = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (!(norm > 0.0) || !std::isfinite(norm))
        throw std::invalid_argument("invalid Huygens normal");
    for (double value : electric)
        if (!std::isfinite(value))
            throw std::invalid_argument("invalid Huygens electric field");
    for (double value : magnetic)
        if (!std::isfinite(value))
            throw std::invalid_argument("invalid Huygens magnetic field");
    for (double& value : normal)
        value /= norm;
    const auto j = cross(normal, magnetic);
    auto m = cross(normal, electric);
    for (double& value : m)
        value = -value;
    return {j, m};
}

std::complex<double> impedance_sheet_transmission(std::complex<double> sheet, std::complex<double> medium1,
                                                  std::complex<double> medium2) {
    const auto denominator = medium1 + medium2 + sheet;
    if (std::abs(denominator) < 1.0e-30)
        throw std::invalid_argument("singular impedance-sheet transmission");
    return 2.0 * medium2 / denominator;
}

std::complex<double> interpolate_emc_transfer(std::span<const EmcTransferSample> samples, double frequency) {
    if (samples.empty() || !(frequency >= 0.0) || !std::isfinite(frequency))
        throw std::invalid_argument("invalid EMC transfer request");
    if (frequency <= samples.front().frequency_hz)
        return samples.front().gain;
    if (frequency >= samples.back().frequency_hz)
        return samples.back().gain;
    for (std::size_t i = 1U; i < samples.size(); ++i) {
        if (!(samples[i].frequency_hz > samples[i - 1U].frequency_hz))
            throw std::invalid_argument("EMC transfer frequencies must be strictly increasing");
        if (frequency <= samples[i].frequency_hz) {
            const double t =
                (frequency - samples[i - 1U].frequency_hz) / (samples[i].frequency_hz - samples[i - 1U].frequency_hz);
            return samples[i - 1U].gain + t * (samples[i].gain - samples[i - 1U].gain);
        }
    }
    return samples.back().gain;
}

} // namespace cfd::fdtd
