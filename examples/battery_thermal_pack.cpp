#include "common/benchmark.hpp"
#include "cfd/battery/thermal3d.hpp"

#include <algorithm>
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
    const std::size_t nx = options.quick ? 6U : 12U * options.scale;
    const std::size_t ny = options.quick ? 2U : 4U * options.scale;
    const std::size_t nz = 2U;
    const std::size_t steps = options.quick ? 12U : 60U;
    const std::size_t voxel_count = nx * ny * nz;

    Timer setup;
    std::vector<LithiumIonCellConfig> cells(3U, default_graphite_nmc_config());
    cells[0].contact_resistance_ohm = 0.05;

    ThermalGridConfig grid;
    grid.cells = {nx, ny, nz};
    grid.spacing_m = {0.01, 0.01, 0.01};
    for (auto& boundary : grid.boundaries)
        boundary = {8.0, 298.15};
    const double voxel_volume = grid.spacing_m[0] * grid.spacing_m[1] * grid.spacing_m[2];
    const double cell_capacity = cells[0].thermal.mass_kg * cells[0].thermal.heat_capacity_j_per_kg_k;
    const double voxels_per_cell = static_cast<double>(voxel_count) / 3.0;
    ThermalVoxel material{cell_capacity / (voxels_per_cell * voxel_volume), {2.0, 0.5, 0.5}};

    std::vector<std::vector<std::size_t>> regions(3U);
    for (std::size_t voxel = 0; voxel < voxel_count; ++voxel) {
        const std::size_t x = voxel % nx;
        const std::size_t region = std::min<std::size_t>(2U, (3U * x) / nx);
        regions[region].push_back(voxel);
    }
    ElectrothermalBatteryPack pack(cells, {}, ThermalGrid3D(grid, {material}), regions);
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

    *csv << "time_s,voltage_v,cell_heat_w,cooling_w,peak_temperature_k,energy_residual_j\n";
    Timer simulation;
    ElectrothermalPackResult last{};
    for (std::size_t step = 0; step < steps; ++step) {
        last = pack.step(0.5, 1.0);
        *csv << last.electrical.time_s << ',' << last.electrical.terminal_voltage_v << ','
             << last.thermal.supplied_power_w << ',' << last.thermal.outward_cooling_power_w << ','
             << last.thermal.maximum_temperature_k << ',' << last.thermal.energy_balance_error_j << '\n';
    }
    const double simulation_ms = simulation.milliseconds();
    emit({"phase16a_battery", "electrothermal_pack", "cpu", "voxel_steps", voxel_count, steps, setup_ms, simulation_ms,
          static_cast<double>(voxel_count * steps),
          last.electrical.terminal_voltage_v + last.thermal.maximum_temperature_k +
              last.thermal.energy_balance_error_j});
    return 0;
}
