#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::rf {

struct HertzianFarField {
    std::complex<double> e_theta;
    std::complex<double> h_phi;
};

[[nodiscard]] HertzianFarField hertzian_dipole_far_field(double length_m,double current_a,
                                                          double frequency_hz,double theta_rad,
                                                          double distance_m);
[[nodiscard]] double shielding_effectiveness_db(double incident_magnitude,double transmitted_magnitude);

struct ArrayElement {
    double x_m{}, y_m{}, z_m{};
    std::complex<double> weight{1.0,0.0};
};

[[nodiscard]] std::complex<double> array_factor(std::span<const ArrayElement> elements,
                                                double frequency_hz,double theta_rad,double phi_rad);

struct PolarizationMetrics {
    double axial_ratio{};
    double tilt_angle_rad{};
    double ellipticity_angle_rad{};
    bool right_hand{};
};
[[nodiscard]] PolarizationMetrics polarization_metrics(std::complex<double> e_theta,
                                                        std::complex<double> e_phi);

struct PatternSample {
    double theta_rad{};
    double normalized_power{};
};

class SinusoidalDipoleSolver {
public:
    SinusoidalDipoleSolver(double length_m, double frequency_hz, double feed_current_a = 1.0);

    [[nodiscard]] double wavelength() const noexcept;
    [[nodiscard]] double electrical_length() const noexcept;
    [[nodiscard]] double pattern_factor(double theta_rad) const;
    [[nodiscard]] std::complex<double> far_e_theta(double theta_rad, double distance_m) const;
    [[nodiscard]] double radiated_power(std::size_t integration_points = 8192U) const;
    [[nodiscard]] double radiation_resistance(std::size_t integration_points = 8192U) const;
    [[nodiscard]] double directivity(std::size_t integration_points = 8192U) const;
    [[nodiscard]] std::vector<PatternSample> normalized_pattern(std::size_t samples = 181U) const;

private:
    double length_{};
    double frequency_{};
    double current_{};
};

} // namespace cfd::rf
