#include "common/benchmark.hpp"
#include "cfd/battery/dfn.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::battery;

    const auto options = parse_options(argc, argv);
    const std::size_t steps = options.quick ? 24U : 120U;
    const double dt_s = 1.0;

    Timer setup;
    auto config = default_graphite_nmc_dfn_config();
    // The reduced quick mesh keeps the walkthrough responsive while preserving
    // the same distributed particle/electrolyte/potential calculation.
    if (options.quick) {
        config.negative_nodes = 8U;
        config.separator_nodes = 3U;
        config.positive_nodes = 8U;
        config.negative_particle.radial_cells = 8U;
        config.positive_particle.radial_cells = 8U;
    } else {
        config.negative_nodes *= options.scale;
        config.separator_nodes *= options.scale;
        config.positive_nodes *= options.scale;
    }
    // Keep the example's inventory check deterministic: the reaction partition
    // moves lithium between electrodes without an electrolyte source term.
    config.transference_number = 1.0;
    DoyleFullerNewmanModel dfn(config);
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

    *csv << "time_s,current_a,terminal_voltage_v,open_circuit_voltage_v,reaction_overpotential_v,"
            "ohmic_drop_v,negative_surface_stoichiometry,positive_surface_stoichiometry,"
            "electrolyte_min_mol_m3,electrolyte_max_mol_m3,temperature_k,heat_generation_w\n";
    auto initial = dfn.state();
    *csv << initial.time_s << ',' << initial.current_a << ',' << initial.terminal_voltage_v << ','
         << initial.open_circuit_voltage_v << ',' << initial.reaction_overpotential_v << ',' << initial.ohmic_drop_v
         << ',' << initial.negative_surface_stoichiometry << ',' << initial.positive_surface_stoichiometry << ','
         << initial.electrolyte_min_concentration_mol_per_m3 << ',' << initial.electrolyte_max_concentration_mol_per_m3
         << ',' << initial.temperature_k << ',' << initial.heat_generation_w << '\n';

    Timer simulation;
    DfnStepResult last = initial;
    for (std::size_t step = 0; step < steps; ++step) {
        last = dfn.step(0.5, dt_s);
        *csv << last.time_s << ',' << last.current_a << ',' << last.terminal_voltage_v << ','
             << last.open_circuit_voltage_v << ',' << last.reaction_overpotential_v << ',' << last.ohmic_drop_v << ','
             << last.negative_surface_stoichiometry << ',' << last.positive_surface_stoichiometry << ','
             << last.electrolyte_min_concentration_mol_per_m3 << ',' << last.electrolyte_max_concentration_mol_per_m3
             << ',' << last.temperature_k << ',' << last.heat_generation_w << '\n';
    }
    const double simulation_ms = simulation.milliseconds();
    const std::size_t node_count = config.negative_nodes + config.separator_nodes + config.positive_nodes;
    emit({"phase16a_battery", "dfn_discharge_profile", "cpu", "distributed_node_steps", node_count, steps, setup_ms,
          simulation_ms, static_cast<double>(node_count * steps),
          last.terminal_voltage_v + last.negative_surface_stoichiometry + last.positive_surface_stoichiometry +
              dfn.total_lithium_mol()});
    return 0;
}
