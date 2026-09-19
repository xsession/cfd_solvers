#include "common/benchmark.hpp"
#include "cfd/battery/pack.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::battery;

    const auto options = parse_options(argc, argv);
    const std::size_t cell_count = options.quick ? 3U : 6U * options.scale;
    const std::size_t steps = options.quick ? 12U : 60U;

    Timer setup;
    std::vector<LithiumIonCellConfig> configs(cell_count, default_graphite_nmc_config());
    for (std::size_t i = 0; i < configs.size(); ++i) {
        // Deliberately spread the positive-electrode state so balancing has
        // something visible to do in the CSV profile.
        configs[i].positive.particle.initial_stoichiometry = 0.28 + 0.015 * static_cast<double>(i);
    }
    PackControl control;
    control.balancing_current_a = 0.05;
    control.balancing_voltage_threshold_v = 0.002;
    SeriesBatteryPack pack(configs, control);
    const double setup_ms = setup.milliseconds();

    std::unique_ptr<std::ofstream> file;
    std::ostream* csv = &std::cout;
    if (!options.output.empty()) {
        const std::filesystem::path path(options.output);
        if (path.has_parent_path())
            std::filesystem::create_directories(path.parent_path());
        file = std::make_unique<std::ofstream>(path);
        if (!*file) {
            std::cerr << "cannot open output file: " << path << '\n';
            return 2;
        }
        csv = file.get();
    }

    *csv << "time_s,pack_voltage_v,terminal_power_w,shunt_heat_w,terminal_energy_wh,"
            "shunt_energy_wh\n";
    Timer simulation;
    PackStepResult last{};
    for (std::size_t step = 0; step < steps; ++step) {
        last = pack.step(-0.1, 1.0);
        *csv << last.time_s << ',' << last.terminal_voltage_v << ',' << last.terminal_power_w << ','
             << last.balancing_heat_w << ',' << last.terminal_energy_wh << ',' << last.balancing_energy_wh << '\n';
    }
    const double simulation_ms = simulation.milliseconds();
    emit({"phase16a_battery", "series_pack_balance", "cpu", "cell_steps", cell_count, steps, setup_ms, simulation_ms,
          static_cast<double>(cell_count * steps),
          last.terminal_voltage_v + last.balancing_energy_wh + last.terminal_energy_wh});
    return 0;
}
