#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::battery {

inline constexpr double faraday_constant = 96485.33212;
inline constexpr double gas_constant = 8.31446261815324;

struct ParticleConfig {
    double radius_m{5.0e-6};
    double diffusivity_m2_per_s{1.0e-14};
    double maximum_concentration_mol_per_m3{3.0e4};
    double initial_stoichiometry{0.5};
    std::size_t radial_cells{24U};
};

class SphericalDiffusionParticle {
public:
    explicit SphericalDiffusionParticle(ParticleConfig config = {});
    void reset(double stoichiometry);
    void step(double outward_molar_flux_mol_per_m2_s, double dt_s);
    [[nodiscard]] double average_concentration() const noexcept;
    [[nodiscard]] double surface_concentration(double outward_molar_flux_mol_per_m2_s = 0.0) const noexcept;
    [[nodiscard]] double average_stoichiometry() const noexcept;
    [[nodiscard]] double surface_stoichiometry(double outward_molar_flux_mol_per_m2_s = 0.0) const noexcept;
    [[nodiscard]] const std::vector<double>& concentrations() const noexcept { return concentration_; }
    [[nodiscard]] const ParticleConfig& config() const noexcept { return config_; }
private:
    ParticleConfig config_{};
    std::vector<double> concentration_{};
};

struct ElectrodeConfig {
    ParticleConfig particle{};
    double thickness_m{80.0e-6};
    double active_volume_fraction{0.60};
    double solid_conductivity_s_per_m{100.0};
    double exchange_current_density_a_per_m2{3.0};
    std::function<double(double)> open_circuit_voltage_v{};
};

struct ElectrolyteConfig {
    double negative_thickness_m{80.0e-6};
    double separator_thickness_m{25.0e-6};
    double positive_thickness_m{75.0e-6};
    double initial_concentration_mol_per_m3{1000.0};
    double diffusivity_m2_per_s{2.0e-10};
    double conductivity_s_per_m{1.0};
    double transference_number{0.38};
    double negative_porosity{0.30};
    double separator_porosity{0.50};
    double positive_porosity{0.30};
    double bruggeman_exponent{1.5};
    std::size_t cells{45U};
};

struct LumpedThermalConfig {
    double mass_kg{0.045};
    double heat_capacity_j_per_kg_k{900.0};
    double heat_transfer_coefficient_w_per_m2_k{8.0};
    double cooling_area_m2{0.012};
    double ambient_temperature_k{298.15};
};

struct ThermalSlabConfig {
    double thickness_m{180.0e-6};
    std::size_t cells{32U};
    double density_kg_per_m3{2200.0};
    double heat_capacity_j_per_kg_k{1000.0};
    double conductivity_w_per_m_k{1.0};
    double convection_w_per_m2_k{8.0};
    double ambient_temperature_k{298.15};
};

class ThermalSlab1D {
public:
    explicit ThermalSlab1D(ThermalSlabConfig config = {});
    void reset(double temperature_k);
    void step(double dt_s, std::span<const double> volumetric_heat_w_per_m3);
    [[nodiscard]] const std::vector<double>& temperature() const noexcept { return temperature_k_; }
    [[nodiscard]] double average_temperature() const noexcept;
private:
    ThermalSlabConfig config_{};
    std::vector<double> temperature_k_{};
};

struct DegradationConfig {
    bool enable_sei{true};
    bool enable_plating{true};
    bool enable_active_material_loss{true};
    bool enable_particle_cracking{true};
    double initial_sei_thickness_m{5.0e-9};
    double sei_rate_m2_per_s{2.0e-20};
    double plating_efficiency{0.02};
    double active_material_loss_per_coulomb{2.0e-10};
    double cracking_rate_per_gradient_second{5.0e-10};
};

struct DegradationState {
    double sei_thickness_m{};
    double plated_lithium_mol{};
    double active_material_fraction{1.0};
    double crack_damage{};
};

struct LithiumIonCellConfig {
    double area_m2{0.010};
    ElectrodeConfig negative{};
    ElectrodeConfig positive{};
    ElectrolyteConfig electrolyte{};
    LumpedThermalConfig thermal{};
    DegradationConfig degradation{};
    double contact_resistance_ohm{0.010};
    double initial_temperature_k{298.15};
};

[[nodiscard]] LithiumIonCellConfig default_graphite_nmc_config();
[[nodiscard]] double graphite_ocv_v(double stoichiometry);
[[nodiscard]] double nmc_ocv_v(double stoichiometry);
[[nodiscard]] double symmetric_butler_volmer_overpotential_v(double interfacial_current_density_a_per_m2,
                                                              double exchange_current_density_a_per_m2,
                                                              double temperature_k);

struct BatteryStepResult {
    double time_s{};
    double current_a{}; // positive = discharge
    double terminal_voltage_v{};
    double open_circuit_voltage_v{};
    double reaction_overpotential_v{};
    double electrolyte_overpotential_v{};
    double ohmic_drop_v{};
    double negative_surface_stoichiometry{};
    double positive_surface_stoichiometry{};
    double electrolyte_min_concentration_mol_per_m3{};
    double electrolyte_max_concentration_mol_per_m3{};
    double temperature_k{};
    double heat_generation_w{};
};

class SingleParticleModel {
public:
    explicit SingleParticleModel(LithiumIonCellConfig config = {});
    void reset();
    [[nodiscard]] BatteryStepResult step(double current_a, double dt_s);
    [[nodiscard]] BatteryStepResult state(double current_a = 0.0) const;
    [[nodiscard]] const SphericalDiffusionParticle& negative_particle() const noexcept { return negative_; }
    [[nodiscard]] const SphericalDiffusionParticle& positive_particle() const noexcept { return positive_; }
    [[nodiscard]] const DegradationState& degradation_state() const noexcept { return degradation_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }
    [[nodiscard]] double temperature_k() const noexcept { return temperature_k_; }
    // Suppress internal lumped heating until reset; an external solver owns T.
    void set_external_temperature(double temperature_k);
    [[nodiscard]] const LithiumIonCellConfig& config() const noexcept { return config_; }
protected:
    [[nodiscard]] double negative_flux(double current_a) const noexcept;
    [[nodiscard]] double positive_flux(double current_a) const noexcept;
    [[nodiscard]] BatteryStepResult evaluate(double current_a,
                                             double electrolyte_overpotential_v,
                                             double electrolyte_min,
                                             double electrolyte_max) const;
    void advance_degradation(double current_a, double dt_s, double negative_surface_gradient);
    void advance_thermal(double current_a, double terminal_voltage_v, double ocv_v, double dt_s);

    LithiumIonCellConfig config_{};
    SphericalDiffusionParticle negative_{};
    SphericalDiffusionParticle positive_{};
    DegradationState degradation_{};
    double temperature_k_{};
    double time_s_{};
    bool external_thermal_{false};
};

class SingleParticleElectrolyteModel : public SingleParticleModel {
public:
    explicit SingleParticleElectrolyteModel(LithiumIonCellConfig config = {});
    void reset();
    [[nodiscard]] BatteryStepResult state(double current_a = 0.0) const;
    [[nodiscard]] BatteryStepResult step(double current_a, double dt_s);
    [[nodiscard]] const std::vector<double>& electrolyte_concentration() const noexcept { return electrolyte_concentration_; }
    [[nodiscard]] const std::vector<double>& electrolyte_potential() const noexcept { return electrolyte_potential_v_; }
    [[nodiscard]] const std::vector<double>& solid_potential() const noexcept { return solid_potential_v_; }
private:
    void step_electrolyte(double current_a, double dt_s);
    [[nodiscard]] double electrolyte_overpotential(double current_a) const;
    void update_potentials(double current_a);

    std::vector<double> electrolyte_concentration_{};
    std::vector<double> electrolyte_potential_v_{};
    std::vector<double> solid_potential_v_{};
};

struct DriveCyclePoint {
    double duration_s{};
    double current_a{};
};

struct DriveCycleResult {
    std::vector<BatteryStepResult> samples;
    double discharged_capacity_ah{};
    double electrical_energy_wh{};
};

[[nodiscard]] DriveCycleResult simulate_drive_cycle(SingleParticleModel& model,
                                                     std::span<const DriveCyclePoint> segments,
                                                     double maximum_step_s = 1.0);
[[nodiscard]] DriveCycleResult simulate_drive_cycle(SingleParticleElectrolyteModel& model,
                                                     std::span<const DriveCyclePoint> segments,
                                                     double maximum_step_s = 1.0);

struct BatteryEisConfig {
    double series_resistance_ohm{0.010};
    double negative_charge_transfer_resistance_ohm{0.020};
    double positive_charge_transfer_resistance_ohm{0.025};
    double negative_double_layer_capacitance_f{2.0};
    double positive_double_layer_capacitance_f{2.0};
    double diffusion_warburg_ohm_sqrt_s{0.015};
    double diffusion_time_s{200.0};
};

[[nodiscard]] std::complex<double> battery_impedance(double frequency_hz, const BatteryEisConfig& config = {});
[[nodiscard]] std::vector<std::complex<double>> battery_impedance_spectrum(std::span<const double> frequency_hz,
                                                                           const BatteryEisConfig& config = {});

} // namespace cfd::battery
