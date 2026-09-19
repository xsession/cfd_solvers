#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/spray.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t nx = options.quick ? 32U : 96U * options.scale;
    const std::size_t ny = options.quick ? 8U : 24U * options.scale;
    const std::size_t steps = options.quick ? 30U : 180U;

    Timer setup;
    auto mesh = make_cartesian_hexa_mesh(nx, ny, 1U, 1.0, 0.25, 1.0);
    SprayConfig config;
    config.dt = options.quick ? 1.0e-4 : 5.0e-5;
    config.evaporation_constant = 5.0e-8;
    config.latent_heat = 2.2e6;
    config.thermal_relaxation_time = 2.0e-3;
    SprayInjectionEvaporation spray(std::move(mesh), config);
    spray.set_injection(
        {{0.10, 0.125, 0.5}, {0.8, 0.0, 0.0}, 1.0e-4, 300.0, 1.0e-5, 0.0, options.quick ? 2.0e-3 : 5.0e-3});
    const std::vector<Vec3> carrier_velocity(spray.mesh().cell_count(), {0.2, 0.0, 0.0});
    const std::vector<double> carrier_temperature(spray.mesh().cell_count(), 450.0);
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    for (std::size_t step = 0; step < steps; ++step)
        spray.step(std::span<const Vec3>(carrier_velocity), std::span<const double>(carrier_temperature));
    const double simulation_ms = simulation.milliseconds();

    const std::filesystem::path path =
        options.output.empty() ? std::filesystem::path("spray_final.csv") : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path);
    if (!file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    file << "x_m,y_m,vapor_mass_source_kg_m3_s,carrier_energy_source_w_m3,carrier_momentum_x_n_m3\n";
    for (std::size_t cell = 0; cell < spray.mesh().cell_count(); ++cell) {
        const auto& center = spray.mesh().cells()[cell].center;
        file << center.x << ',' << center.y << ',' << spray.vapor_mass_source()[cell] << ','
             << spray.carrier_energy_source()[cell] << ',' << spray.carrier_momentum_source()[cell].x << '\n';
    }

    std::cout << "spray parcels=" << spray.parcels().size() << " injected_mass=" << spray.injected_mass()
              << " liquid_mass_remaining=" << spray.total_liquid_mass()
              << " evaporated_mass=" << spray.evaporated_mass() << " csv=" << path << '\n';
    emit({"phase03_fvm", "coupled_spray_injection_evaporation", "cpu", "parcel_steps", spray.mesh().cell_count(), steps,
          setup_ms, simulation_ms, static_cast<double>(spray.mesh().cell_count() * steps),
          spray.injected_mass() + spray.total_liquid_mass() + spray.evaporated_mass()});
    return 0;
}
