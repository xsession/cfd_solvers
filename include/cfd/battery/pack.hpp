#pragma once
#include "cfd/battery/lithium_ion.hpp"
namespace cfd::battery {
struct PackControl {
    double minimum_cell_voltage_v{2.5};
    double maximum_cell_voltage_v{4.3};
    double maximum_temperature_k{333.15};
    // Ideal regulated dissipative shunts, sampled at each step's start.
    double balancing_current_a{0.0};
    double balancing_voltage_threshold_v{0.01};
    double balancing_minimum_voltage_v{3.5};
};
struct PackStepResult {
    std::vector<BatteryStepResult> cells;
    std::vector<double> balancing_currents_a;
    double time_s{}, terminal_voltage_v{}, terminal_power_w{}, balancing_heat_w{};
    double net_terminal_charge_ah{}, terminal_energy_wh{}, balancing_energy_wh{};
};
// Series SPMe string; positive current = discharge. Shunt heat is external.
class SeriesBatteryPack {
public:
    explicit SeriesBatteryPack(std::span<const LithiumIonCellConfig> cells, PackControl control = {});
    // Checks loaded start/end limits. Failure leaves all state unchanged.
    // No within-step event localization; use appropriately short timesteps.
    [[nodiscard]] PackStepResult step(double terminal_current_a, double dt_s);
    [[nodiscard]] const std::vector<SingleParticleElectrolyteModel>& cells() const noexcept { return cells_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }
    void reset();
    // Enables external thermal ownership; reset restores lumped heating.
    void set_external_temperatures(std::span<const double> temperature_k);
private:
    void check_limits(const BatteryStepResult& state) const;
    PackControl control_;
    std::vector<SingleParticleElectrolyteModel> cells_;
    double time_s_{}, charge_ah_{}, energy_wh_{}, balance_wh_{};
};
}
