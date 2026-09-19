#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/euler_euler.hpp"

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
    const std::size_t nx = options.quick ? 40U : 120U * options.scale;
    const std::size_t ny = options.quick ? 12U : 36U * options.scale;
    const std::size_t steps = options.quick ? 24U : 180U;

    Timer setup;
    EulerEulerConfig config;
    config.dt = options.quick ? 2.0e-4 : 1.0e-4;
    config.primary_density = 1000.0;
    config.dispersed_density = 8.0;
    config.drag_coefficient = 18.0;
    config.dispersed_body_acceleration = {0.0, -0.8, 0.0};
    auto mesh = make_cartesian_hexa_mesh(nx, ny, 1U, 1.0, 0.30, 1.0);
    EulerEulerTransport phases(std::move(mesh), config);
    phases.initialize(
        [](Vec3 p) {
            const double dx = (p.x - 0.30) / 0.16;
            const double dy = (p.y - 0.20) / 0.09;
            return 0.04 + 0.42 * std::exp(-(dx * dx + dy * dy));
        },
        {0.08, 0.0, 0.0}, {0.34, 0.0, 0.0});
    const double initial_volume = phases.volume_fraction_integral();
    const Vec3 initial_momentum = phases.total_momentum();
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    for (std::size_t step = 0; step < steps; ++step)
        phases.step();
    const double simulation_ms = simulation.milliseconds();

    const std::filesystem::path path =
        options.output.empty() ? std::filesystem::path("euler_euler_final.csv") : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path);
    if (!file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    file << "x_m,y_m,dispersed_volume_fraction,primary_u_m_s,dispersed_u_m_s,drag_x_n_m3,drag_y_n_m3\n";
    for (std::size_t cell = 0; cell < phases.mesh().cell_count(); ++cell) {
        const auto& center = phases.mesh().cells()[cell].center;
        const auto& primary = phases.primary_velocity()[cell];
        const auto& dispersed = phases.dispersed_velocity()[cell];
        const auto& drag = phases.drag_force()[cell];
        file << center.x << ',' << center.y << ',' << phases.dispersed_fraction()[cell] << ',' << primary.x << ','
             << dispersed.x << ',' << drag.x << ',' << drag.y << '\n';
    }

    const Vec3 final_momentum = phases.total_momentum();
    std::cout << "euler_euler dispersed_volume_initial=" << initial_volume
              << " dispersed_volume_final=" << phases.volume_fraction_integral()
              << " fraction_min=" << phases.minimum_fraction() << " fraction_max=" << phases.maximum_fraction()
              << " momentum_initial_x=" << initial_momentum.x << " momentum_final_x=" << final_momentum.x
              << " csv=" << path << '\n';
    emit({"phase03_fvm", "coupled_euler_euler_drag", "cpu", "cell_steps", phases.mesh().cell_count(), steps, setup_ms,
          simulation_ms, static_cast<double>(phases.mesh().cell_count() * steps),
          phases.volume_fraction_integral() + phases.maximum_fraction() + std::abs(final_momentum.x)});
    return 0;
}
