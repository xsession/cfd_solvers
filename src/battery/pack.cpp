#include "cfd/battery/pack.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace cfd::battery {
namespace {
bool nonnegative(double x) { return std::isfinite(x) && x >= 0.0; }
}
SeriesBatteryPack::SeriesBatteryPack(std::span<const LithiumIonCellConfig> configs, PackControl control)
    : control_(control) {
    if (configs.empty() || !nonnegative(control.minimum_cell_voltage_v)
        || !std::isfinite(control.maximum_cell_voltage_v)
        || control.maximum_cell_voltage_v <= control.minimum_cell_voltage_v
        || !std::isfinite(control.maximum_temperature_k) || control.maximum_temperature_k <= 0.0
        || !nonnegative(control.balancing_current_a)
        || !nonnegative(control.balancing_voltage_threshold_v)
        || !nonnegative(control.balancing_minimum_voltage_v))
        throw std::invalid_argument("invalid series-pack configuration");
    for (const auto& c : configs) cells_.emplace_back(c);
}
void SeriesBatteryPack::check_limits(const BatteryStepResult& s) const {
    if (!std::isfinite(s.terminal_voltage_v) || !std::isfinite(s.temperature_k)
        || s.terminal_voltage_v < control_.minimum_cell_voltage_v
        || s.terminal_voltage_v > control_.maximum_cell_voltage_v
        || s.temperature_k <= 0.0 || s.temperature_k > control_.maximum_temperature_k)
        throw std::domain_error("series-pack cell voltage or temperature cutoff");
}
void SeriesBatteryPack::set_external_temperatures(std::span<const double> temperatures) {
    if (temperatures.size()!=cells_.size()) throw std::invalid_argument("pack temperature count mismatch");
    auto trial=cells_;
    for (std::size_t i=0;i<trial.size();++i) {
        trial[i].set_external_temperature(temperatures[i]);
        check_limits(trial[i].state());
    }
    cells_=std::move(trial);
}
void SeriesBatteryPack::reset() {
    for (auto& c : cells_) c.reset();
    time_s_ = charge_ah_ = energy_wh_ = balance_wh_ = 0.0;
}
PackStepResult SeriesBatteryPack::step(double current, double dt) {
    if (!std::isfinite(current) || !std::isfinite(dt) || dt <= 0.0 || !std::isfinite(time_s_ + dt))
        throw std::invalid_argument("invalid series-pack step");
    auto trial = cells_;
    PackStepResult out;
    double lowest = std::numeric_limits<double>::infinity();
    std::vector<double> voltages;
    for (const auto& c : cells_) {
        const auto s = c.state();
        check_limits(s);
        voltages.push_back(s.terminal_voltage_v);
        lowest = std::min(lowest, s.terminal_voltage_v);
    }
    for (std::size_t i = 0; i < trial.size(); ++i) {
        const double bleed = voltages[i] >= control_.balancing_minimum_voltage_v
            && voltages[i] - lowest > control_.balancing_voltage_threshold_v
            ? control_.balancing_current_a : 0.0;
        const double cell_current = current + bleed;
        if (!std::isfinite(cell_current)) throw std::invalid_argument("pack cell current overflow");
        check_limits(trial[i].state(cell_current));
        const auto s = trial[i].step(cell_current, dt);
        check_limits(s);
        out.cells.push_back(s);
        out.balancing_currents_a.push_back(bleed);
        out.terminal_voltage_v += s.terminal_voltage_v;
        out.balancing_heat_w += bleed * s.terminal_voltage_v;
    }
    out.time_s = time_s_ + dt;
    out.terminal_power_w = current * out.terminal_voltage_v;
    out.net_terminal_charge_ah = charge_ah_ + current * dt / 3600.0;
    // Backward rectangle energy quadrature, as in the cell drive-cycle API.
    out.terminal_energy_wh = energy_wh_ + out.terminal_power_w * dt / 3600.0;
    out.balancing_energy_wh = balance_wh_ + out.balancing_heat_w * dt / 3600.0;
    if (!std::isfinite(out.terminal_voltage_v) || !std::isfinite(out.balancing_heat_w)
        || !std::isfinite(out.net_terminal_charge_ah) || !std::isfinite(out.terminal_energy_wh)
        || !std::isfinite(out.balancing_energy_wh))
        throw std::overflow_error("series-pack accounting overflow");
    cells_ = std::move(trial);
    time_s_ = out.time_s;
    charge_ah_ = out.net_terminal_charge_ah;
    energy_wh_ = out.terminal_energy_wh;
    balance_wh_ = out.balancing_energy_wh;
    return out;
}
}
