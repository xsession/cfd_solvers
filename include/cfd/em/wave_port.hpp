#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;

struct RectangularWaveguideMode {
    enum class Family { te, tm };
    Family family{Family::te};
    std::size_t m{};
    std::size_t n{};
    double frequency_hz{};
    double cutoff_frequency_hz{};
    Complex propagation_constant_rad_per_m{}; // beta for exp(-j beta z)
    Complex wave_impedance_ohm{};
    bool propagating{};
    // Scale applied to the transverse electric-field shape. Propagating modes
    // are normalized to exactly 1 W time-average forward power.
    double electric_scale{};
};

struct TransverseModeFields {
    std::array<Complex,2> electric_v_per_m{};
    std::array<Complex,2> magnetic_a_per_m{};
};

// Analytic rectangular PEC wave-port eigenmodes, sorted by cutoff frequency.
// TE modes allow one zero index (but not TE00); TM modes require m,n >= 1.
[[nodiscard]] std::vector<RectangularWaveguideMode> rectangular_waveguide_modes(
    double width_m,double height_m,double frequency_hz,
    double relative_permittivity=1.0,double relative_permeability=1.0,
    std::size_t mode_count=4U);

[[nodiscard]] TransverseModeFields rectangular_waveguide_transverse_fields(
    const RectangularWaveguideMode& mode,double width_m,double height_m,
    double x_m,double y_m);

// Numerical power check/integration utility for modal normalization.
[[nodiscard]] double rectangular_mode_power_w(const RectangularWaveguideMode& mode,
                                               double width_m,double height_m,
                                               std::size_t nx=64U,std::size_t ny=64U);

} // namespace cfd::em
