#pragma once

#include "cfd/rf/nport.hpp"

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

namespace cfd::em {

using Complex = std::complex<double>;

struct EyeDiagramMetrics {
    std::size_t symbol_count{};
    std::size_t samples_per_symbol{};
    double zero_mean{};
    double one_mean{};
    double zero_sigma{};
    double one_sigma{};
    double decision_threshold{};
    double eye_height{};
    double eye_width_ui{};
    double estimated_ber{};
};

// Computes NRZ eye metrics from one uniformly sampled waveform and a known
// ideal bit stream. The decision sample defaults to the UI centre. Eye width is
// the fraction of sample phases where the worst-case one remains above the
// worst-case zero.
[[nodiscard]] EyeDiagramMetrics analyze_nrz_eye(
    std::span<const double> waveform_v,
    std::span<const int> symbols,
    std::size_t samples_per_symbol,
    std::size_t decision_sample = 0U);

struct DecouplingCapacitor {
    double capacitance_f{};
    double esr_ohm{};
    double esl_h{};
    std::size_t count{1U};
};

[[nodiscard]] Complex decoupling_impedance_ohm(const DecouplingCapacitor& capacitor,
                                               double frequency_hz);
[[nodiscard]] Complex pdn_parallel_impedance_ohm(
    std::span<const DecouplingCapacitor> capacitors,
    double frequency_hz,
    Complex plane_impedance_ohm = {});

struct DecouplingOptimizationResult {
    std::vector<DecouplingCapacitor> selected;
    double worst_ratio_before{};
    double worst_ratio_after{};
    std::size_t iterations{};
};
[[nodiscard]] DecouplingOptimizationResult optimize_decoupling_greedy(
    std::span<const double> frequencies_hz,
    std::span<const double> target_impedance_ohm,
    std::span<const DecouplingCapacitor> library,
    std::size_t maximum_parts,
    Complex plane_impedance_ohm = {});

struct PdnGridConfig {
    std::size_t nx{2U},ny{2U};
    double horizontal_resistance_ohm{1.0e-3};
    double vertical_resistance_ohm{1.0e-3};
    std::size_t max_iterations{20000U};
    double relative_tolerance{1.0e-12};
};
struct PdnVoltageSource { std::size_t node{}; double voltage_v{}; };
struct PdnIrDropResult {
    std::vector<double> voltage_v;
    double minimum_voltage_v{};
    double maximum_drop_v{};
    std::size_t iterations{};
    bool converged{};
};
// Positive load_current_a consumes current from the PDN node to return/ground.
[[nodiscard]] PdnIrDropResult solve_pdn_ir_drop_grid(
    const PdnGridConfig& config,
    std::span<const double> load_current_a,
    std::span<const PdnVoltageSource> voltage_sources);

[[nodiscard]] double normalized_double_exponential_pulse(double time_s,double amplitude,
                                                          double rise_tau_s,double fall_tau_s);
[[nodiscard]] std::vector<double> sample_double_exponential_pulse(double duration_s,double dt_s,
                                                                   double amplitude,double rise_tau_s,
                                                                   double fall_tau_s);
[[nodiscard]] std::vector<double> sample_damped_sine(double duration_s,double dt_s,double amplitude,
                                                     double frequency_hz,double decay_tau_s);
struct ProbeMetrics {
    double peak{};
    double rms{};
    double impulse{};
    double energy{};
};
[[nodiscard]] ProbeMetrics probe_waveform(std::span<const double> samples,double dt_s);

} // namespace cfd::em
