#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace cfd::rf {

struct ThinWireMomConfig {
    double length_m{};
    double radius_m{};
    double frequency_hz{};
    std::size_t segments{31U}; // odd: feed is the center segment
    double feed_voltage_v{1.0};
    std::size_t quadrature_order{8U};
};

struct ThinWireMomResult {
    std::vector<double> segment_center_z_m;
    std::vector<std::complex<double>> current_a;
    std::complex<double> feed_impedance_ohm{};
    double feed_current_a{};
};

// Thin straight PEC wire along z, pulse basis + point matching applied to
// Pocklington's thin-wire EFIE. A delta-gap electric field is imposed on the
// center segment. This is a readable reference solver, not a NEC replacement.
[[nodiscard]] ThinWireMomResult solve_center_fed_thin_wire(const ThinWireMomConfig& config);

} // namespace cfd::rf
