#pragma once

#include "cfd/circuit/spice.hpp"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::circuit {

struct AcPort {
    std::string positive;
    std::string negative{"0"};
};

// Extract a linear two-port Z matrix by independent 1-A AC current injection,
// then convert to S parameters. Existing independent AC sources in the circuit
// should be zero for a passive network characterization.
[[nodiscard]] cfd::rf::Matrix2C two_port_s_parameters(const Circuit& circuit,
                                                       const AcPort& first,
                                                       const AcPort& second,
                                                       double frequency_hz,
                                                       double reference_impedance = 50.0);

struct SensitivityResult {
    double baseline_output{};
    double derivative{}; // output volts per unit parameter
};

[[nodiscard]] SensitivityResult voltage_source_dc_sensitivity(Circuit circuit,
                                                               std::string_view source,
                                                               std::string_view output_node,
                                                               double perturbation,
                                                               const NewtonConfig& config = {});

struct FourierMeasurement {
    double frequency_hz{};
    Complex phasor{}; // peak-amplitude convention for a real signal
    double magnitude{};
    double phase_rad{};
};
[[nodiscard]] FourierMeasurement fourier_measurement(std::span<const TransientPoint> samples,
                                                      std::size_t node,
                                                      double frequency_hz);
[[nodiscard]] double total_harmonic_distortion(std::span<const TransientPoint> samples,
                                                std::size_t node,double fundamental_hz,
                                                std::size_t harmonics = 5U);

struct ResistorTolerance {
    std::string element;
    double relative_sigma{};
};
struct MonteCarloResult {
    std::size_t samples{};
    double mean{};
    double standard_deviation{};
    double minimum{};
    double maximum{};
};
[[nodiscard]] MonteCarloResult monte_carlo_dc_resistors(const Circuit& nominal,
                                                         std::string_view output_node,
                                                         std::span<const ResistorTolerance> tolerances,
                                                         std::size_t samples,
                                                         std::uint64_t seed = 1U,
                                                         const NewtonConfig& config = {});



struct ParameterSweepPoint {
    double parameter{};
    OperatingPoint operating_point;
};
using CircuitParameterSetter = std::function<void(Circuit&,double)>;
[[nodiscard]] std::vector<ParameterSweepPoint> parameter_sweep_dc(const Circuit& nominal,
    std::span<const double> values,const CircuitParameterSetter& setter,const NewtonConfig& config = {});
[[nodiscard]] std::vector<ParameterSweepPoint> temperature_sweep_dc(const Circuit& nominal,
    std::span<const double> temperature_k,const NewtonConfig& config = {});

struct DominantPoleEstimate {
    double dc_gain{};
    double pole_frequency_hz{};
    bool found{};
};
// Estimate a dominant low-pass pole from the first -3 dB crossing of the AC response.
[[nodiscard]] DominantPoleEstimate estimate_dominant_pole(const Circuit& circuit,
    std::string_view output_node,double start_hz,double stop_hz,std::size_t points,
    const NewtonConfig& config = {});

struct NoiseSpectrum {
    std::vector<NoiseResult> points;
    double integrated_rms_voltage{};
};
// Log-spaced output-noise sweep. Integrated RMS uses trapezoidal integration of PSD in Hz.
[[nodiscard]] NoiseSpectrum output_noise_spectrum(const Circuit& circuit,std::string_view output_node,
    double start_hz,double stop_hz,std::size_t points,double temperature_k = 300.0,
    const NewtonConfig& config = {});

// SPICE-style pnjlim core for PN-junction Newton stabilization.
[[nodiscard]] double limit_pn_junction_voltage(double proposed,double previous,
                                                double thermal_voltage,double saturation_current);

struct CircuitPort {
    std::size_t positive{};
    std::size_t negative{};
    double reference_impedance{50.0};
};
// Small-signal circuit N-port extraction. Each port is voltage-forced in turn;
// the MNA source currents form Y, then Y -> Z -> S. Internal AC sources should
// be disabled by the caller for passive network extraction.
[[nodiscard]] cfd::rf::ComplexMatrix extract_s_parameters(Circuit circuit,
    std::span<const CircuitPort> ports,double frequency_hz,const NewtonConfig& config = {});

} // namespace cfd::circuit
