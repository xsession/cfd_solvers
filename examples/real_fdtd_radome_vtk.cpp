#include "common/benchmark.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::filesystem::path output_dir(int argc, char** argv) {
    std::filesystem::path out = "examples/output/fdtd_radome";
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--output")
            out = argv[++i];
    }
    return out;
}

void write_vtk_line(const cfd::fdtd::Maxwell1D& solver, const std::filesystem::path& path) {
    const auto& e = solver.electric();
    const auto& h = solver.magnetic();
    std::ofstream o(path);
    if (!o)
        throw std::runtime_error("cannot create VTK file: " + path.string());
    o << "# vtk DataFile Version 3.0\n"
      << "cfd_solvers 1D dielectric radome pulse\nASCII\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << e.size() << " 1 1\n"
      << "ORIGIN 0 0 0\n"
      << std::setprecision(17) << "SPACING 1 1 1\n"
      << "POINT_DATA " << e.size() << "\n"
      << "SCALARS Ez double 1\nLOOKUP_TABLE default\n";
    for (double x : e)
        o << x << '\n';
    o << "SCALARS Hy double 1\nLOOKUP_TABLE default\n";
    for (double x : h)
        o << x << '\n';
    o << "SCALARS energy_density double 1\nLOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < e.size(); ++i)
        o << 0.5 * (e[i] * e[i] + h[i] * h[i]) << '\n';
}

} // namespace

int main(int argc, char** argv) {
    using namespace cfd::examples;
    const auto opt = parse_options(argc, argv);
    const std::size_t cells = opt.quick ? 256U : 1600U * opt.scale;
    const std::size_t steps = opt.quick ? 80U : 900U;
    const auto out = output_dir(argc, argv);
    std::filesystem::create_directories(out);

    Timer setup;
    cfd::fdtd::Maxwell1D solver({cells, 2.0e-4, 0.9, 1.0, 1.0, cfd::fdtd::Boundary1D::cpml, 32U, 3.0, 1e-9, 6.0, 0.0});
    solver.set_material(cells * 52U / 100U, cells * 64U / 100U, 2.9, 0.0005);
    solver.set_material(cells * 70U / 100U, cells * 76U / 100U, 4.4, 0.002);
    solver.initialize_gaussian(0.22, 0.025);
    const double setup_ms = setup.milliseconds();

    Timer sim;
    const std::size_t first_snapshot = steps / 3U;
    const std::size_t second_snapshot = 2U * steps / 3U;
    for (std::size_t n = 0; n < steps; ++n) {
        solver.step();
        if (n + 1U == first_snapshot)
            write_vtk_line(solver, out / "radome_mid.vtk");
        if (n + 1U == second_snapshot)
            write_vtk_line(solver, out / "radome_late.vtk");
    }
    write_vtk_line(solver, out / "radome_final.vtk");
    const double sim_ms = sim.milliseconds();

    std::cout << "Wrote ParaView VTK snapshots in: " << out << "\n";
    emit({"real_fdtd_radome", "dielectric_window_pulse_vtk", "cpu", "cell_steps", cells, steps, setup_ms, sim_ms,
          double(cells) * double(steps), solver.energy()});
    return 0;
}
