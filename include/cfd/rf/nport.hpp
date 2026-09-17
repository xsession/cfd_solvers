#pragma once

#include "cfd/rf/network.hpp"

#include <complex>
#include <cstddef>
#include <iosfwd>
#include <functional>
#include <span>
#include <vector>

namespace cfd::rf {

class ComplexMatrix {
public:
    ComplexMatrix() = default;
    explicit ComplexMatrix(std::size_t size, Complex value = {});

    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] Complex& operator()(std::size_t row, std::size_t column);
    [[nodiscard]] const Complex& operator()(std::size_t row, std::size_t column) const;
    [[nodiscard]] const std::vector<Complex>& data() const noexcept { return data_; }

    [[nodiscard]] static ComplexMatrix identity(std::size_t size);

private:
    std::size_t size_{};
    std::vector<Complex> data_;
};

[[nodiscard]] ComplexMatrix multiply(const ComplexMatrix& lhs, const ComplexMatrix& rhs);
[[nodiscard]] ComplexMatrix inverse(const ComplexMatrix& matrix);

// Power-wave N-port conversion for positive real per-port reference impedances.
[[nodiscard]] ComplexMatrix z_to_s(const ComplexMatrix& z, std::span<const double> reference_impedance);
[[nodiscard]] ComplexMatrix s_to_z(const ComplexMatrix& s, std::span<const double> reference_impedance);
[[nodiscard]] ComplexMatrix renormalize_s(const ComplexMatrix& s,
                                          std::span<const double> old_reference_impedance,
                                          std::span<const double> new_reference_impedance);

// Kurokawa power-wave conversion for complex per-port reference impedances.
// Each reference must have a strictly positive real part. These overloads are
// useful for waveguide, active-device and de-embedding workflows where the
// reference impedance is not purely real.
[[nodiscard]] ComplexMatrix z_to_s_power_wave(const ComplexMatrix& z,
                                              std::span<const Complex> reference_impedance);
[[nodiscard]] ComplexMatrix s_to_z_power_wave(const ComplexMatrix& s,
                                              std::span<const Complex> reference_impedance);
[[nodiscard]] ComplexMatrix renormalize_s_power_wave(const ComplexMatrix& s,
                                                     std::span<const Complex> old_reference_impedance,
                                                     std::span<const Complex> new_reference_impedance);

struct NPortPoint {
    double frequency_hz{};
    ComplexMatrix s;
    std::vector<double> reference_impedance;
};

// Touchstone SnP reader covering v1 option lines and the common v2 keyword form.
// Matrix ordering follows Touchstone column-major network-data order:
// S11,S21,...,SN1,S12,...,SNN.
[[nodiscard]] std::vector<NPortPoint> read_touchstone(std::istream& input, std::size_t port_count_hint = 0U);
void write_touchstone(std::ostream& output, std::span<const NPortPoint> points, bool version2 = true);

// Shift each reference plane toward the DUT by distance[i]. Positive distance
// removes that length of a uniform line. gamma is alpha+j*beta for each port.
[[nodiscard]] ComplexMatrix shift_reference_planes(const ComplexMatrix& s,
                                                    std::span<const Complex> gamma,
                                                    std::span<const double> distance);

[[nodiscard]] double return_loss_db(Complex reflection);
[[nodiscard]] double vswr(Complex reflection);

struct StabilityMetrics {
    double rollett_k{};
    double mu{};
    Complex determinant{};
};
[[nodiscard]] StabilityMetrics stability_metrics(const Matrix2C& s);

struct StabilityCircle {
    Complex center{};
    double radius{};
    bool valid{};
};
[[nodiscard]] StabilityCircle source_stability_circle(const Matrix2C& s);
[[nodiscard]] StabilityCircle load_stability_circle(const Matrix2C& s);
[[nodiscard]] double transducer_gain(const Matrix2C& s,Complex source_reflection = {},
                                     Complex load_reflection = {});

struct NoiseParameters {
    double minimum_noise_factor{1.0};
    Complex optimum_source_reflection{};
    double equivalent_noise_resistance_ohm{};
    double reference_impedance_ohm{50.0};
};
struct NoiseCircle {
    Complex center{};
    double radius{};
    double noise_factor{};
    bool valid{};
};
[[nodiscard]] double noise_factor(const NoiseParameters& parameters,Complex source_reflection);
[[nodiscard]] NoiseCircle noise_circle(const NoiseParameters& parameters,double target_noise_factor);

struct LoadPullSample {
    Complex load_reflection{};
    double transducer_gain{};
};
[[nodiscard]] std::vector<LoadPullSample> sample_load_pull(const Matrix2C& s,Complex source_reflection,
                                                           std::size_t radial_steps = 10U,
                                                           std::size_t angular_steps = 72U,
                                                           double maximum_radius = 0.95);

struct NetworkQuality {
    double maximum_singular_value{};
    double reciprocity_error{};
    bool passive{};
    bool reciprocal{};
};
// Passivity is tested from the largest singular value of S. Reciprocity uses
// the Frobenius norm of S-S^T. Tolerances are relative numerical guards.
[[nodiscard]] NetworkQuality network_quality(const ComplexMatrix& s,
                                             double passivity_tolerance = 1.0e-10,
                                             double reciprocity_tolerance = 1.0e-10);

// Conservative passivity projection baseline. If sigma_max(S) exceeds the
// requested bound the full matrix is scaled uniformly, preserving reciprocity
// and phase while guaranteeing the singular-value power bound.
[[nodiscard]] ComplexMatrix enforce_passivity(const ComplexMatrix& s,
                                              double maximum_singular_value = 1.0);

struct BroadbandCausality {
    double negative_time_energy_ratio{};
    bool causal{};
    std::size_t impulse_samples{};
};
// Discrete broadband causality gate for a uniformly sampled DC-to-Nyquist
// network sweep. A Hermitian extension is inverse-transformed and energy in
// wrapped negative-time samples is reported.
[[nodiscard]] BroadbandCausality check_broadband_causality(
    std::span<const NPortPoint> points, double negative_time_energy_tolerance = 1.0e-6,
    double frequency_uniformity_tolerance = 1.0e-8);

struct AdaptiveNetworkSweepConfig {
    std::size_t initial_points{5U};
    std::size_t maximum_points{257U};
    std::size_t maximum_refinements{12U};
    double relative_tolerance{1.0e-3};
    double absolute_tolerance{1.0e-6};
};
struct AdaptiveNetworkSweepResult {
    std::vector<NPortPoint> points;
    std::size_t evaluations{};
    std::size_t refinements{};
    bool converged{};
};
using NetworkSampler = std::function<ComplexMatrix(double frequency_hz)>;
// Log-frequency adaptive sweep using midpoint interpolation error as the
// refinement indicator. This is the sampling foundation for later rational MOR.
[[nodiscard]] AdaptiveNetworkSweepResult adaptive_network_sweep(
    double start_hz,double stop_hz,std::span<const double> reference_impedance,
    NetworkSampler sampler,const AdaptiveNetworkSweepConfig& config = {});

// 4-port single-ended -> [differential pair 1, differential pair 2,
// common pair 1, common pair 2] mixed-mode transform for port pairs (1,2),(3,4).
[[nodiscard]] ComplexMatrix single_ended_to_mixed_mode(const ComplexMatrix& four_port_s);

struct LMatch {
    // Series reactance followed by shunt susceptance at the load node.
    double series_reactance_ohm{};
    double shunt_susceptance_siemens{};
    bool shunt_first{};
};
// Match two positive real resistances at one frequency. Returns the low-pass
// solution; shunt_first reports whether the shunt element is placed at the source side.
[[nodiscard]] LMatch synthesize_l_match(double source_resistance, double load_resistance);

} // namespace cfd::rf
