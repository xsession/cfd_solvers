#pragma once

#include "cfd/rf/network.hpp"

#include <complex>
#include <cstddef>
#include <iosfwd>
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
