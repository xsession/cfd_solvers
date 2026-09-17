#pragma once

#include <complex>
#include <span>
#include <vector>

namespace cfd::particle {

using Complex = std::complex<double>;

// Causal single resonator longitudinal point-charge wake W(s) [V/C], sampled
// at distance s >= 0 behind the driving charge.
[[nodiscard]] double resonator_longitudinal_wake_v_per_c(
    double distance_m,double shunt_impedance_ohm,double resonant_frequency_hz,double quality_factor);

// Convolves a causal sampled point-charge wake with a bunch line-charge
// density [C/m]. Both arrays use the same distance spacing.
[[nodiscard]] std::vector<double> bunch_wake_potential_v(
    std::span<const double> wake_v_per_c,
    std::span<const double> line_charge_c_per_m,
    double spacing_m);

// Longitudinal impedance Z(omega)=integral W(t) exp(-j omega t) dt,
// with t=s/c and a causal sampled wake W(s).
[[nodiscard]] Complex wake_impedance_ohm(
    std::span<const double> wake_v_per_c,double spacing_m,double frequency_hz);

} // namespace cfd::particle
