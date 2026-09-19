#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/compressible_vof.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t nx = options.quick ? 48U : 128U * options.scale;
    const std::size_t ny = options.quick ? 8U : 32U * options.scale;
    const std::size_t steps = options.quick ? 24U : 180U;

    Timer setup;
    auto mesh = make_cartesian_hexa_mesh(nx, ny, 1U, 1.0, 0.25, 1.0);
    CompressibleVofConfig config;
    config.dt = options.quick ? 2.0e-5 : 1.0e-5;
    config.liquid_density = 1000.0;
    config.gas_density = 1.2;
    config.liquid_bulk_modulus = 2.2e9;
    config.gas_gamma = 1.4;
    config.reference_pressure = 1.0e5;
    CompressibleVofTransport vof(std::move(mesh), config);
    vof.initialize(
        [](Vec3 p) {
            const double transition = (p.x - 0.46) / 0.025;
            return 0.5 * (1.0 - std::tanh(transition));
        },
        config.reference_pressure);
    std::vector<double> face_flux(vof.mesh().face_count(), 0.0);
    for (std::size_t face = 0; face < vof.mesh().face_count(); ++face) {
        const auto& f = vof.mesh().faces()[face];
        if (!f.boundary())
            face_flux[face] = 1.0e-3 * (0.5 - f.center.x) * f.area.x;
    }
    vof.set_face_flux(std::move(face_flux));
    const double initial_volume = vof.liquid_volume();
    const double initial_mass = vof.total_mass();
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    for (std::size_t step = 0; step < steps; ++step)
        vof.step();
    const double simulation_ms = simulation.milliseconds();

    const std::filesystem::path path = options.output.empty() ? std::filesystem::path("compressible_vof_final.csv")
                                                              : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path);
    if (!file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    file << "x_m,y_m,liquid_volume_fraction,mixture_density_kg_m3,pressure_pa,liquid_density_kg_m3,gas_density_kg_m3\n";
    for (std::size_t cell = 0; cell < vof.mesh().cell_count(); ++cell) {
        const auto& center = vof.mesh().cells()[cell].center;
        file << center.x << ',' << center.y << ',' << vof.liquid_fraction()[cell] << ',' << vof.density()[cell] << ','
             << vof.pressure()[cell] << ',' << vof.liquid_density()[cell] << ',' << vof.gas_density()[cell] << '\n';
    }

    std::cout << "compressible_vof liquid_volume_initial=" << initial_volume
              << " liquid_volume_final=" << vof.liquid_volume() << " mass_initial=" << initial_mass
              << " mass_final=" << vof.total_mass() << " pressure_min=" << vof.minimum_pressure()
              << " pressure_max=" << vof.maximum_pressure() << " csv=" << path << '\n';
    emit({"phase03_fvm", "compressible_vof_barotropic", "cpu", "cell_steps", vof.mesh().cell_count(), steps, setup_ms,
          simulation_ms, static_cast<double>(vof.mesh().cell_count() * steps),
          vof.liquid_volume() + vof.total_mass() + vof.maximum_pressure()});
    return 0;
}
