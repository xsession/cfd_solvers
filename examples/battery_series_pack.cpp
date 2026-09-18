#include "cfd/battery/pack.hpp"
#include <iostream>
int main() {
    using namespace cfd::battery;
    std::vector<LithiumIonCellConfig> configs(3, default_graphite_nmc_config());
    configs[0].positive.particle.initial_stoichiometry = 0.30;
    PackControl control;
    control.balancing_current_a = 0.05;
    SeriesBatteryPack pack(configs, control);
    std::cout << "time_s,pack_voltage_v,terminal_power_w,shunt_heat_w,terminal_energy_wh,shunt_energy_wh\n";
    for (int n = 0; n < 60; ++n) {
        const auto r = pack.step(-0.1, 1.0);
        std::cout << r.time_s << ',' << r.terminal_voltage_v << ',' << r.terminal_power_w << ','
                  << r.balancing_heat_w << ',' << r.terminal_energy_wh << ',' << r.balancing_energy_wh << '\n';
    }
}
