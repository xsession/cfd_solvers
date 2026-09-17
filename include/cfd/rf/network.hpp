#pragma once

#include <complex>
#include <cstddef>
#include <iosfwd>
#include <vector>

namespace cfd::rf {

using Complex = std::complex<double>;

struct Matrix2C {
    Complex a11{}, a12{}, a21{}, a22{};
};

[[nodiscard]] Matrix2C z_to_s(const Matrix2C& z, double reference_impedance = 50.0);
[[nodiscard]] Matrix2C s_to_z(const Matrix2C& s, double reference_impedance = 50.0);
[[nodiscard]] Matrix2C s_to_abcd(const Matrix2C& s, double reference_impedance = 50.0);
[[nodiscard]] Matrix2C abcd_to_s(const Matrix2C& abcd, double reference_impedance = 50.0);
[[nodiscard]] Matrix2C cascade_abcd(const Matrix2C& lhs, const Matrix2C& rhs) noexcept;
[[nodiscard]] Matrix2C cascade_s(const Matrix2C& lhs, const Matrix2C& rhs,
                                 double reference_impedance = 50.0);

struct TransmissionLine {
    double characteristic_impedance{50.0};
    Complex propagation_constant{}; // alpha + j*beta [1/m]
    double length{};                 // m

    [[nodiscard]] Matrix2C abcd() const;
    [[nodiscard]] Matrix2C s_parameters(double reference_impedance = 50.0) const;
};

struct TemTransmissionLine {
    double characteristic_impedance{50.0};
    double phase_velocity{299792458.0};
    double attenuation_nepers_per_m{};
    double length{};

    [[nodiscard]] TransmissionLine at_frequency(double frequency_hz) const;
    [[nodiscard]] Matrix2C s_parameters(double frequency_hz,double reference_impedance = 50.0) const;
};

struct MicrostripQuasiStatic {
    double characteristic_impedance{};
    double effective_permittivity{};
    double phase_velocity{};
    double guided_wavelength{};
};

struct CoaxQuasiStatic {
    double characteristic_impedance{};
    double phase_velocity{};
    double propagation_constant_rad_per_m{};
    double guided_wavelength{};
};

struct RectangularWaveguideTE10 {
    double cutoff_frequency_hz{};
    double propagation_constant_rad_per_m{};
    double guide_wavelength{};
    double wave_impedance{};
};

// Hammerstad/Jensen-style quasi-static baseline for an uncoated microstrip.
[[nodiscard]] MicrostripQuasiStatic microstrip_quasi_static(double width, double substrate_height,
                                                             double relative_permittivity,
                                                             double frequency_hz);

[[nodiscard]] CoaxQuasiStatic coax_quasi_static(double inner_radius,double outer_radius,
                                                   double relative_permittivity,double frequency_hz);
[[nodiscard]] RectangularWaveguideTE10 rectangular_waveguide_te10(double broad_wall,double narrow_wall,
                                                                   double relative_permittivity,
                                                                   double relative_permeability,
                                                                   double frequency_hz);



struct PlanarLineQuasiStatic {
    double characteristic_impedance{};
    double effective_permittivity{};
    double phase_velocity{};
    double guided_wavelength{};
};

// Symmetric homogeneous stripline baseline. width and ground_spacing are metres.
[[nodiscard]] PlanarLineQuasiStatic stripline_quasi_static(double width,double ground_spacing,
                                                            double relative_permittivity,double frequency_hz);

// Infinite-ground coplanar-waveguide baseline without conductor thickness.
// gap is the slot from centre strip edge to each ground plane.
[[nodiscard]] PlanarLineQuasiStatic coplanar_waveguide_quasi_static(double center_width,double gap,
                                                                     double relative_permittivity,double frequency_hz);

struct Touchstone2PortPoint {
    double frequency_hz{};
    Matrix2C s{};
};

// Touchstone v1-style S2P reader. Supports RI, MA and DB data with R reference impedance.
[[nodiscard]] std::vector<Touchstone2PortPoint> read_touchstone_s2p(std::istream& input,
                                                                    double* reference_impedance = nullptr);
void write_touchstone_s2p(std::ostream& output, const std::vector<Touchstone2PortPoint>& points,
                          double reference_impedance = 50.0);

} // namespace cfd::rf
