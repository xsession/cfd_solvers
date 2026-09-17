#pragma once

#include "cfd/rf/nport.hpp"

#include <complex>
#include <cstddef>
#include <functional>
#include <utility>
#include <istream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace cfd::circuit {

using Complex = std::complex<double>;

struct DiodeModel {
    double saturation_current{1.0e-14};
    double emission_coefficient{1.0};
    double temperature_k{300.0};
};

struct MosLevel1Model {
    bool pmos{};
    double threshold_voltage{0.7};
    double transconductance_parameter{1.0e-3}; // A/V^2
    double channel_length_modulation{};        // 1/V
};

struct BjtModel {
    bool pnp{};
    double saturation_current{1.0e-15};
    double forward_beta{100.0};
    double reverse_beta{1.0};
    double temperature_k{300.0};
};

struct JfetModel {
    bool pchannel{};
    double pinch_off_voltage{-2.0}; // V; n-channel convention is negative
    double beta{1.0e-3};            // A/V^2
    double channel_length_modulation{};
};

struct SwitchModel {
    double on_resistance{1.0};
    double off_resistance{1.0e12};
    double threshold_voltage{};
    double hysteresis_voltage{};
};

struct NewtonConfig {
    std::size_t max_iterations{80U};
    double voltage_tolerance{1.0e-10};
    double residual_tolerance{1.0e-10};
    double gmin{1.0e-12};
};

struct OperatingPoint {
    std::vector<double> node_voltage;   // includes node 0 = ground
    // Branch order: independent V, VCVS, CCVS, then inductors.
    std::vector<double> branch_current;
    std::size_t iterations{};
    bool converged{};
};

struct AcSolution {
    double frequency_hz{};
    std::vector<Complex> node_voltage;   // includes ground
    std::vector<Complex> branch_current;
};

struct TransientPoint {
    double time{};
    std::vector<double> node_voltage;
    std::vector<double> branch_current;
};

enum class TransientMethod {
    backward_euler,
    trapezoidal,
    bdf2
};

using SourceWaveform = std::function<double(double)>;

struct PulseWaveform {
    double initial{};
    double pulsed{};
    double delay{};
    double rise{};
    double fall{};
    double width{};
    double period{};
    [[nodiscard]] double operator()(double time) const;
};

struct SineWaveform {
    double offset{};
    double amplitude{};
    double frequency_hz{};
    double delay{};
    double damping_per_s{};
    double phase_deg{};
    [[nodiscard]] double operator()(double time) const;
};

struct PwlWaveform {
    std::vector<std::pair<double,double>> points;
    [[nodiscard]] double operator()(double time) const;
};

struct DcSweepPoint {
    double source_value{};
    OperatingPoint operating_point;
};

struct TemperatureSweepPoint {
    double temperature_k{};
    OperatingPoint operating_point;
};

struct NoiseResult {
    double frequency_hz{};
    double output_noise_v_per_sqrt_hz{};
    double output_noise_density_v2_per_hz{};
};

struct StaticDeviceEvaluation {
    // Terminal currents are positive flowing out of each listed terminal.
    std::vector<double> terminal_current;
    // Row-major dI_terminal[row]/dV_terminal[column].
    std::vector<double> jacobian;
};
using StaticDeviceEvaluator = std::function<StaticDeviceEvaluation(std::span<const double>)>;

class Circuit {
public:
    Circuit();

    [[nodiscard]] std::size_t node(std::string_view name);
    [[nodiscard]] std::size_t find_node(std::string_view name) const;
    [[nodiscard]] std::size_t node_count() const noexcept { return node_names_.size(); }
    [[nodiscard]] const std::vector<std::string>& node_names() const noexcept { return node_names_; }

    void add_resistor(std::string name,std::size_t positive,std::size_t negative,double resistance_ohm);
    void add_capacitor(std::string name,std::size_t positive,std::size_t negative,double capacitance_f);
    void add_inductor(std::string name,std::size_t positive,std::size_t negative,double inductance_h);
    void add_mutual_inductance(std::string name,std::string first_inductor,std::string second_inductor,double coupling);
    void add_voltage_source(std::string name,std::size_t positive,std::size_t negative,double dc_voltage,
                            Complex ac_voltage = {}, SourceWaveform transient_waveform = {});
    void add_current_source(std::string name,std::size_t positive,std::size_t negative,double dc_current,
                            Complex ac_current = {}, SourceWaveform transient_waveform = {});

    void add_vccs(std::string name,std::size_t positive,std::size_t negative,
                  std::size_t control_positive,std::size_t control_negative,double transconductance);
    void add_vcvs(std::string name,std::size_t positive,std::size_t negative,
                  std::size_t control_positive,std::size_t control_negative,double voltage_gain);
    void add_cccs(std::string name,std::size_t positive,std::size_t negative,
                  std::string control_voltage_source,double current_gain);
    void add_ccvs(std::string name,std::size_t positive,std::size_t negative,
                  std::string control_voltage_source,double transresistance);

    void add_diode_model(std::string name,DiodeModel model);
    void add_diode(std::string name,std::size_t anode,std::size_t cathode,std::string model);
    void add_mos_model(std::string name,MosLevel1Model model);
    void add_mosfet(std::string name,std::size_t drain,std::size_t gate,std::size_t source,
                    std::size_t bulk,std::string model);
    void add_bjt_model(std::string name,BjtModel model);
    void add_bjt(std::string name,std::size_t collector,std::size_t base,std::size_t emitter,std::string model);
    void add_jfet_model(std::string name,JfetModel model);
    void add_jfet(std::string name,std::size_t drain,std::size_t gate,std::size_t source,std::string model);
    void add_switch_model(std::string name,SwitchModel model);
    void add_switch(std::string name,std::size_t positive,std::size_t negative,
                    std::size_t control_positive,std::size_t control_negative,std::string model);
    // Frequency-domain N-port connected between each positive/negative node pair.
    // Sampled S-parameters are linearly interpolated in frequency.
    void add_nport(std::string name,std::vector<std::pair<std::size_t,std::size_t>> ports,
                   std::vector<cfd::rf::NPortPoint> samples);
    void add_tem_transmission_line(std::string name,
                                   std::pair<std::size_t,std::size_t> first_port,
                                   std::pair<std::size_t,std::size_t> second_port,
                                   const cfd::rf::TemTransmissionLine& line,
                                   std::span<const double> sample_frequencies_hz,
                                   double reference_impedance = 50.0);
    void add_microstrip_line(std::string name,
                             std::pair<std::size_t,std::size_t> first_port,
                             std::pair<std::size_t,std::size_t> second_port,
                             double width,double substrate_height,double relative_permittivity,
                             double length,std::span<const double> sample_frequencies_hz,
                             double reference_impedance = 50.0);
    void add_stripline_line(std::string name,
                            std::pair<std::size_t,std::size_t> first_port,
                            std::pair<std::size_t,std::size_t> second_port,
                            double width,double ground_spacing,double relative_permittivity,
                            double length,std::span<const double> sample_frequencies_hz,
                            double reference_impedance = 50.0);
    void add_coplanar_waveguide_line(std::string name,
                                     std::pair<std::size_t,std::size_t> first_port,
                                     std::pair<std::size_t,std::size_t> second_port,
                                     double center_width,double gap,double relative_permittivity,
                                     double length,std::span<const double> sample_frequencies_hz,
                                     double reference_impedance = 50.0);
    void add_rectangular_waveguide_te10_line(std::string name,
                                            std::pair<std::size_t,std::size_t> first_port,
                                            std::pair<std::size_t,std::size_t> second_port,
                                            double broad_wall,double narrow_wall,
                                            double relative_permittivity,double relative_permeability,
                                            double length,std::span<const double> sample_frequencies_hz,
                                            double reference_impedance = 50.0);
    // Generic memoryless nonlinear multi-terminal compact-device seam.
    // The evaluator returns terminal currents and their Jacobian. This is the
    // clean internal target for future Verilog-A/OSDI adapters.
    void add_static_device(std::string name,std::vector<std::size_t> terminals,StaticDeviceEvaluator evaluator);

    void set_voltage_source_dc(std::string_view name,double voltage);
    [[nodiscard]] double voltage_source_dc(std::string_view name) const;
    void set_resistance(std::string_view name,double resistance_ohm);
    [[nodiscard]] double resistance(std::string_view name) const;
    // Update temperature-dependent compact models (diode/BJT in the current baseline).
    void set_device_temperature(double temperature_k);

    [[nodiscard]] OperatingPoint dc_operating_point(const NewtonConfig& config = {}) const;
    [[nodiscard]] OperatingPoint dc_operating_point_homotopy(const NewtonConfig& config = {},
                                                              std::size_t source_steps = 8U,
                                                              std::size_t gmin_steps = 6U) const;
    [[nodiscard]] AcSolution ac(double frequency_hz,const OperatingPoint* operating_point = nullptr,
                                const NewtonConfig& config = {}) const;
    [[nodiscard]] std::vector<AcSolution> ac_log_sweep(double start_hz,double stop_hz,std::size_t points,
                                                       const NewtonConfig& config = {}) const;
    [[nodiscard]] std::vector<DcSweepPoint> dc_sweep_voltage_source(std::string_view source,
                                                                    double start,double stop,std::size_t points,
                                                                    const NewtonConfig& config = {});
    [[nodiscard]] std::vector<TemperatureSweepPoint> dc_temperature_sweep(double start_k,double stop_k,
                                                                            std::size_t points,
                                                                            const NewtonConfig& config = {}) const;
    [[nodiscard]] std::vector<TransientPoint> transient(double time_step,std::size_t steps,
                                                         const NewtonConfig& config = {}) const;
    [[nodiscard]] std::vector<TransientPoint> transient(double time_step,std::size_t steps,
                                                         TransientMethod method,
                                                         const NewtonConfig& config = {}) const;
    [[nodiscard]] NoiseResult output_noise(double frequency_hz,std::string_view output_node,
                                           double temperature_k = 300.0,
                                           const NewtonConfig& config = {}) const;

    [[nodiscard]] double voltage(const OperatingPoint& point,std::string_view node_name) const;
    [[nodiscard]] Complex voltage(const AcSolution& point,std::string_view node_name) const;
    [[nodiscard]] Complex branch_current(const AcSolution& point,std::string_view branch_name) const;

    [[nodiscard]] static Circuit parse_spice(std::istream& input);
    [[nodiscard]] static Circuit parse_spice(std::string_view text);
    // File-aware parser with recursive .INCLUDE and .LIB <file> <section> support.
    [[nodiscard]] static Circuit parse_spice_file(const std::string& path,
                                                  std::string_view library_section = {});

private:
    struct Resistor { std::string name; std::size_t p{},n{}; double value{}; };
    struct Capacitor { std::string name; std::size_t p{},n{}; double value{}; };
    struct Inductor { std::string name; std::size_t p{},n{}; double value{}; };
    struct MutualInductance { std::string name,first,second; double coupling{}; };
    struct VoltageSource { std::string name; std::size_t p{},n{}; double dc{}; Complex ac{}; SourceWaveform transient; };
    struct CurrentSource { std::string name; std::size_t p{},n{}; double dc{}; Complex ac{}; SourceWaveform transient; };
    struct Vccs { std::string name; std::size_t p{},n{},cp{},cn{}; double gain{}; };
    struct Vcvs { std::string name; std::size_t p{},n{},cp{},cn{}; double gain{}; };
    struct Cccs { std::string name; std::size_t p{},n{}; std::string control; double gain{}; };
    struct Ccvs { std::string name; std::size_t p{},n{}; std::string control; double gain{}; };
    struct Diode { std::string name; std::size_t a{},k{}; std::string model; };
    struct Mosfet { std::string name; std::size_t d{},g{},s{},b{}; std::string model; };
    struct Bjt { std::string name; std::size_t c{},b{},e{}; std::string model; };
    struct Jfet { std::string name; std::size_t d{},g{},s{}; std::string model; };
    struct Switch { std::string name; std::size_t p{},n{},cp{},cn{}; std::string model; };
    struct NPortDevice { std::string name; std::vector<std::pair<std::size_t,std::size_t>> ports; std::vector<cfd::rf::NPortPoint> samples; };
    struct StaticDevice { std::string name; std::vector<std::size_t> terminals; StaticDeviceEvaluator evaluator; };

    std::vector<std::string> node_names_;
    std::unordered_map<std::string,std::size_t> node_lookup_;
    std::vector<Resistor> resistors_;
    std::vector<Capacitor> capacitors_;
    std::vector<Inductor> inductors_;
    std::vector<MutualInductance> mutual_inductances_;
    std::vector<VoltageSource> voltage_sources_;
    std::vector<CurrentSource> current_sources_;
    std::vector<Vccs> vccs_;
    std::vector<Vcvs> vcvs_;
    std::vector<Cccs> cccs_;
    std::vector<Ccvs> ccvs_;
    std::vector<Diode> diodes_;
    std::vector<Mosfet> mosfets_;
    std::vector<Bjt> bjts_;
    std::vector<Jfet> jfets_;
    std::vector<Switch> switches_;
    std::vector<NPortDevice> nports_;
    std::vector<StaticDevice> static_devices_;
    std::unordered_map<std::string,DiodeModel> diode_models_;
    std::unordered_map<std::string,MosLevel1Model> mos_models_;
    std::unordered_map<std::string,BjtModel> bjt_models_;
    std::unordered_map<std::string,JfetModel> jfet_models_;
    std::unordered_map<std::string,SwitchModel> switch_models_;

    [[nodiscard]] std::size_t branch_count() const noexcept;
    [[nodiscard]] std::size_t branch_index(std::string_view name) const;
    [[nodiscard]] std::size_t inductor_branch(std::size_t inductor_index) const noexcept;
    [[nodiscard]] OperatingPoint dc_operating_point_scaled(const NewtonConfig& config,
                                                            double source_scale,
                                                            const OperatingPoint* initial) const;
};

[[nodiscard]] double parse_spice_number(std::string_view text);
[[nodiscard]] double evaluate_spice_expression(
    std::string_view expression,
    const std::unordered_map<std::string,double>& parameters = {});

} // namespace cfd::circuit
