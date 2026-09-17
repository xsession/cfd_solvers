#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;

struct FrequencyDomain1DConfig {
    std::size_t points{201U};
    double length_m{1.0};
    double frequency_hz{1.0e8};
    double relative_permittivity{1.0};
    double relative_permeability{1.0};
    double conductivity_s_per_m{};
};

struct FrequencyDomain1DResult {
    std::vector<Complex> electric_v_per_m;
    // Magnetic samples are staggered between adjacent electric samples.
    std::vector<Complex> magnetic_a_per_m;
    double spacing_m{};
};

// 1-D transverse driven Maxwell/Helmholtz baseline with PEC boundaries and
// exp(+j*omega*t) phasor convention. The impressed current density is sampled
// on the electric-field nodes and must have config.points entries.
[[nodiscard]] FrequencyDomain1DResult solve_pec_driven_maxwell_1d(
    const FrequencyDomain1DConfig& config,
    std::span<const Complex> impressed_current_density_a_per_m2);

struct CavityMode1D {
    double angular_frequency_rad_per_s{};
    double frequency_hz{};
    std::vector<double> electric_shape;
};

// Generalized finite-difference eigenmode baseline for a PEC 1-D cavity.
// relative_permittivity is nodal and must contain `points` positive values.
[[nodiscard]] std::vector<CavityMode1D> pec_cavity_eigenmodes_1d(
    std::size_t points,double length_m,
    std::span<const double> relative_permittivity,
    double relative_permeability,
    std::size_t mode_count);

[[nodiscard]] std::vector<CavityMode1D> pec_cavity_eigenmodes_1d(
    std::size_t points,double length_m,double relative_permittivity,
    double relative_permeability,std::size_t mode_count);

} // namespace cfd::em
