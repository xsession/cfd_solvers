#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::tcad {

inline constexpr double elementary_charge_c = 1.602176634e-19;
inline constexpr double boltzmann_j_per_k = 1.380649e-23;
inline constexpr double vacuum_permittivity_f_per_m = 8.8541878128e-12;

[[nodiscard]] double thermal_voltage(double temperature_k);
[[nodiscard]] double bernoulli_function(double x);
[[nodiscard]] double silicon_bandgap_ev(double temperature_k,
                                        double eg0_ev = 1.17,
                                        double alpha_ev_per_k = 4.73e-4,
                                        double beta_k = 636.0);
[[nodiscard]] double temperature_scaled_intrinsic_density(double ni_reference_m3,
                                                          double reference_temperature_k,
                                                          double temperature_k,
                                                          double bandgap_reference_ev,
                                                          double bandgap_ev);
[[nodiscard]] double temperature_scaled_mobility(double mobility_reference_m2_per_vs,
                                                 double reference_temperature_k,
                                                 double temperature_k,
                                                 double exponent);

struct SemiconductorMaterial {
    double relative_permittivity{11.7};
    double intrinsic_density_m3{1.0e16};
    double electron_mobility_m2_per_vs{0.135};
    double hole_mobility_m2_per_vs{0.048};
    double conduction_band_density_m3{2.8e25};
    double valence_band_density_m3{1.04e25};
    double temperature_k{300.0};
};

struct TemperatureDependence {
    bool enabled{true};
    double varshni_eg0_ev{1.17};
    double varshni_alpha_ev_per_k{4.73e-4};
    double varshni_beta_k{636.0};
    double electron_mobility_exponent{2.42};
    double hole_mobility_exponent{2.20};
    double density_of_states_exponent{1.5};
};

struct RecombinationModel {
    double tau_n_s{0.0};
    double tau_p_s{0.0};
    double trap_n1_m3{0.0};
    double trap_p1_m3{0.0};
    double auger_n_m6_per_s{0.0};
    double auger_p_m6_per_s{0.0};
    double radiative_m3_per_s{0.0};
};

[[nodiscard]] double srh_recombination(double n_m3, double p_m3, double ni_m3,
                                       const RecombinationModel& model);
[[nodiscard]] double auger_recombination(double n_m3, double p_m3, double ni_m3,
                                         const RecombinationModel& model);
[[nodiscard]] double radiative_recombination(double n_m3, double p_m3, double ni_m3,
                                             const RecombinationModel& model);
[[nodiscard]] double total_recombination(double n_m3, double p_m3, double ni_m3,
                                         const RecombinationModel& model);

[[nodiscard]] double caughey_thomas_mobility(double doping_abs_m3,
                                             double mu_min,
                                             double mu_max,
                                             double reference_m3,
                                             double alpha);
[[nodiscard]] double high_field_mobility(double low_field_mobility,
                                         double electric_field_v_per_m,
                                         double saturation_velocity_m_per_s,
                                         double beta = 2.0);

enum class ContactType { ohmic, schottky, insulating, gate };

struct Contact1D {
    ContactType type{ContactType::ohmic};
    double voltage_v{0.0};
    double electron_barrier_ev{0.7};
    // For gate contacts, this is the oxide capacitance per unit area used by
    // the electrostatic Robin boundary. Carrier flux is zero through gate and
    // insulating contacts.
    double oxide_capacitance_f_per_m2{0.0};
    double flatband_voltage_v{0.0};
};

struct Semiconductor1DConfig {
    double length_m{1.0e-6};
    std::size_t nodes{101U};
    double area_m2{1.0e-12};
    SemiconductorMaterial material{};
    TemperatureDependence temperature_dependence{};
    RecombinationModel recombination{};
    Contact1D left{};
    Contact1D right{};
    std::size_t max_gummel_iterations{200U};
    double relative_tolerance{1.0e-8};
    double under_relaxation{0.35};
    double carrier_floor_m3{1.0};
};

struct DcResult {
    bool converged{};
    std::size_t iterations{};
    double relative_change{};
    double terminal_current_a{};
};

struct IvPoint {
    double voltage_v{};
    double current_a{};
    bool converged{};
};

struct TransientStepResult {
    bool converged{};
    std::size_t iterations{};
    double relative_change{};
    double time_s{};
    double conduction_current_a{};
    double displacement_current_a{};
    double terminal_current_a{};
};

struct SmallSignalResult {
    double bias_voltage_v{};
    double frequency_hz{};
    double conductance_s{};
    double capacitance_f{};
    std::complex<double> admittance_s{};
    std::complex<double> impedance_ohm{};
    bool converged{};
};

struct CvPoint {
    double voltage_v{};
    double capacitance_f{};
    double electrode_charge_c{};
    bool converged{};
};

class SemiconductorDevice1D {
public:
    explicit SemiconductorDevice1D(Semiconductor1DConfig config = {});

    void set_net_doping(std::span<const double> net_doping_m3);
    void set_net_doping(const std::function<double(double)>& net_doping_m3);
    void set_temperature(double temperature_k);
    void set_contact_voltages(double left_v, double right_v);

    void initialize_charge_neutral();
    [[nodiscard]] DcResult solve_equilibrium();
    [[nodiscard]] DcResult solve_dc();
    [[nodiscard]] std::vector<IvPoint> sweep_right_contact(double start_v, double stop_v, std::size_t points);
    [[nodiscard]] TransientStepResult step_transient(double dt_s);
    [[nodiscard]] SmallSignalResult small_signal(double frequency_hz,
                                                 double perturbation_v = 1.0e-4) const;
    [[nodiscard]] std::vector<CvPoint> sweep_cv(double start_left_v,
                                                double stop_left_v,
                                                std::size_t points,
                                                double perturbation_v = 1.0e-4) const;

    [[nodiscard]] const Semiconductor1DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& x() const noexcept { return x_; }
    [[nodiscard]] const std::vector<double>& potential() const noexcept { return potential_v_; }
    [[nodiscard]] const std::vector<double>& electron_density() const noexcept { return electron_m3_; }
    [[nodiscard]] const std::vector<double>& hole_density() const noexcept { return hole_m3_; }
    [[nodiscard]] const std::vector<double>& net_doping() const noexcept { return net_doping_m3_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }

    [[nodiscard]] std::vector<double> electron_current_density_faces() const;
    [[nodiscard]] std::vector<double> hole_current_density_faces() const;
    [[nodiscard]] std::vector<double> total_current_density_faces() const;
    [[nodiscard]] std::vector<double> electric_field_faces() const;
    [[nodiscard]] std::vector<double> joule_heating_cells() const;
    [[nodiscard]] double terminal_current_a() const;
    [[nodiscard]] double left_terminal_charge_c() const;
    [[nodiscard]] double total_mobile_charge_c() const;

private:
    [[nodiscard]] double contact_electron_density(const Contact1D& contact, double doping_m3) const;
    [[nodiscard]] double contact_hole_density(const Contact1D& contact, double doping_m3) const;
    [[nodiscard]] double contact_potential(const Contact1D& contact, double doping_m3) const;
    [[nodiscard]] bool carrier_blocking(const Contact1D& contact) const noexcept;
    void enforce_contact_carriers();
    void apply_blocking_carrier_boundaries();
    void solve_poisson_linear();
    void solve_electron_continuity();
    void solve_hole_continuity();
    void solve_electron_continuity_transient(std::span<const double> old_density_m3, double dt_s);
    void solve_hole_continuity_transient(std::span<const double> old_density_m3, double dt_s);
    [[nodiscard]] double maximum_relative_change(std::span<const double> old_values,
                                                 std::span<const double> new_values) const;

    Semiconductor1DConfig config_{};
    SemiconductorMaterial reference_material_{};
    std::vector<double> x_{};
    std::vector<double> potential_v_{};
    std::vector<double> electron_m3_{};
    std::vector<double> hole_m3_{};
    std::vector<double> net_doping_m3_{};
    double time_s_{};
};

[[nodiscard]] std::vector<double> pn_junction_doping(std::span<const double> x,
                                                     double junction_position_m,
                                                     double acceptor_m3,
                                                     double donor_m3);
[[nodiscard]] std::vector<double> pin_diode_doping(std::span<const double> x,
                                                   double p_end_m,
                                                   double n_begin_m,
                                                   double acceptor_m3,
                                                   double donor_m3,
                                                   double intrinsic_net_m3 = 0.0);

} // namespace cfd::tcad
