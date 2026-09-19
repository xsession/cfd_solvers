#include "common/benchmark.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/solvers/fvm/vof.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t nx = options.quick ? 24U : 96U * options.scale;
    const std::size_t ny = options.quick ? 12U : 48U * options.scale;
    const std::size_t steps = options.quick ? 16U : 120U;
    const double length_x = 1.0;
    const double length_y = 0.5;
    const double velocity_x = 0.2;
    const double dx = length_x / static_cast<double>(nx);
    const double dt_s = 0.35 * dx / velocity_x;

    Timer setup;
    auto mesh = make_cartesian_hexa_mesh(nx, ny, 1U, length_x, length_y, 0.1);
    VofTransport vof(std::move(mesh), {dt_s, 1.0});
    vof.initialize([](Vec3 p) {
        // A sharp liquid column is transported through a closed rectangular
        // channel. The boundary flux is explicitly zero, so the volume check
        // below is a useful first conservation diagnostic.
        return (p.x < 0.30 && p.y < 0.34) ? 1.0 : 0.0;
    });
    std::vector<Vec3> velocity(vof.alpha().size(), {velocity_x, 0.0, 0.0});
    std::vector<Vec3> boundary_velocity(vof.mesh().face_count(), Vec3{});
    const auto flux = face_flux_from_velocity(vof.mesh(), velocity, boundary_velocity);
    vof.set_face_flux(flux);
    const double initial_volume = vof.volume();
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    for (std::size_t step = 0; step < steps; ++step)
        vof.step();
    const double simulation_ms = simulation.milliseconds();

    const auto capillary_force = csf_surface_tension(vof.mesh(), vof.alpha(), 0.072);
    double force_norm = 0.0;
    for (const auto force : capillary_force)
        force_norm += magnitude(force);
    double alpha_sum = std::accumulate(vof.alpha().begin(), vof.alpha().end(), 0.0);

    std::unique_ptr<std::ofstream> file;
    std::ostream* csv = &std::cout;
    const std::filesystem::path path = options.output.empty() ? std::filesystem::path("vof_dam_break_final.csv")
                                                              : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    file = std::make_unique<std::ofstream>(path);
    if (!*file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    csv = file.get();
    *csv << "x_m,y_m,volume_fraction,csf_force_x,csf_force_y,csf_force_z\n";
    for (std::size_t cell = 0; cell < vof.alpha().size(); ++cell) {
        const auto& p = vof.mesh().cells()[cell].center;
        const auto& force = capillary_force[cell];
        *csv << p.x << ',' << p.y << ',' << vof.alpha()[cell] << ',' << force.x << ',' << force.y << ',' << force.z
             << '\n';
    }

    std::cout << "VOF volume_initial=" << initial_volume << " volume_final=" << vof.volume()
              << " alpha_min=" << *std::min_element(vof.alpha().begin(), vof.alpha().end())
              << " alpha_max=" << *std::max_element(vof.alpha().begin(), vof.alpha().end()) << " csv=" << path << '\n';
    emit({"phase03_fvm", "vof_closed_channel", "cpu", "cell_steps", vof.mesh().cell_count(), steps, setup_ms,
          simulation_ms, static_cast<double>(vof.mesh().cell_count() * steps), vof.volume() + force_norm + alpha_sum});
    return 0;
}
