#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/thin_film.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t nx = options.quick ? 32U : 128U * options.scale;
    const std::size_t ny = options.quick ? 12U : 48U * options.scale;
    const std::size_t steps = options.quick ? 20U : 160U;
    const double length_x = 1.0;
    const double length_y = 0.4;

    Timer setup;
    auto mesh = make_cartesian_hexa_mesh(nx, ny, 1U, length_x, length_y, 1.0);
    ThinFilmConfig config;
    config.dt = options.quick ? 2.0e-5 : 1.0e-5;
    config.density = 1000.0;
    config.viscosity = 1.0e-3;
    config.surface_tension = 0.02;
    ThinFilmTransport film(std::move(mesh), config);
    film.initialize([](Vec3 p) {
        const double dx = (p.x - 0.48) / 0.16;
        const double dy = (p.y - 0.20) / 0.10;
        return 0.01 + 0.004 * std::exp(-(dx * dx + dy * dy));
    });
    const double initial_inventory = film.inventory();
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    for (std::size_t step = 0; step < steps; ++step)
        film.step();
    const double simulation_ms = simulation.milliseconds();

    std::unique_ptr<std::ofstream> file;
    std::ostream* csv = &std::cout;
    const std::filesystem::path path =
        options.output.empty() ? std::filesystem::path("thin_film_final.csv") : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    file = std::make_unique<std::ofstream>(path);
    if (!*file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    csv = file.get();
    *csv << "x_m,y_m,thickness_m,capillary_pressure_pa\n";
    for (std::size_t cell = 0; cell < film.mesh().cell_count(); ++cell) {
        const auto& center = film.mesh().cells()[cell].center;
        *csv << center.x << ',' << center.y << ',' << film.thickness()[cell] << ',' << film.pressure()[cell] << '\n';
    }

    const auto [minimum, maximum] = std::minmax_element(film.thickness().begin(), film.thickness().end());
    std::cout << "thin_film inventory_initial=" << initial_inventory << " inventory_final=" << film.inventory()
              << " thickness_min=" << *minimum << " thickness_max=" << *maximum << " csv=" << path << '\n';
    emit({"phase03_fvm", "thin_film_capillary_relaxation", "cpu", "cell_steps", film.mesh().cell_count(), steps,
          setup_ms, simulation_ms, static_cast<double>(film.mesh().cell_count() * steps),
          film.inventory() + *maximum + *minimum});
    return 0;
}
