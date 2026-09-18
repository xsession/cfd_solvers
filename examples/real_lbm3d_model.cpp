#include "common/benchmark.hpp"
#include "cfd/io/mesh_io.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/visualization.hpp"
#include "cfd/solvers/lbm/voxelize.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Point = cfd::fem::Point3;
using Surface = cfd::io::TriangleSurface;

std::filesystem::path output_dir(int argc, char** argv) {
    std::filesystem::path out = "examples/output/lbm3d_model";
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--output")
            out = argv[++i];
    return out;
}

std::optional<std::filesystem::path> stl_path(int argc, char** argv) {
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--stl")
            return std::filesystem::path(argv[++i]);
    return std::nullopt;
}

Surface make_ahmed_body() {
    // A watertight Ahmed-style bluff body: the sloped rear roof is useful for
    // exercising the same STL-to-voxel pipeline as an external vehicle mesh.
    const std::vector<std::pair<double, double>> profile{{0.15, 0.12}, {0.85, 0.12}, {0.85, 0.27},
                                                         {0.73, 0.46}, {0.28, 0.46}, {0.15, 0.36}};
    constexpr double y_min = 0.22;
    constexpr double y_max = 0.78;
    std::vector<Point> low, high;
    low.reserve(profile.size());
    high.reserve(profile.size());
    for (const auto [x, z] : profile) {
        low.push_back({x, y_min, z});
        high.push_back({x, y_max, z});
    }

    Surface surface;
    const auto n = profile.size();
    for (std::size_t i = 1; i + 1 < n; ++i) {
        surface.triangles.push_back({low[0], low[i], low[i + 1]});
        surface.triangles.push_back({high[0], high[i + 1], high[i]});
    }
    for (std::size_t i = 0; i < n; ++i) {
        const auto next = (i + 1U) % n;
        surface.triangles.push_back({low[i], low[next], high[next]});
        surface.triangles.push_back({low[i], high[next], high[i]});
    }
    return surface;
}

void fit_surface_to_unit_box(Surface& surface) {
    if (surface.triangles.empty())
        throw std::invalid_argument("STL surface contains no triangles");
    Point lo{surface.triangles.front().a.x, surface.triangles.front().a.y, surface.triangles.front().a.z};
    Point hi = lo;
    const auto visit = [&](Point p) {
        lo.x = std::min(lo.x, p.x);
        lo.y = std::min(lo.y, p.y);
        lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x);
        hi.y = std::max(hi.y, p.y);
        hi.z = std::max(hi.z, p.z);
    };
    for (const auto& triangle : surface.triangles) {
        visit(triangle.a);
        visit(triangle.b);
        visit(triangle.c);
    }
    const Point span{hi.x - lo.x, hi.y - lo.y, hi.z - lo.z};
    if (!(span.x > 0.0 && span.y > 0.0 && span.z > 0.0))
        throw std::invalid_argument("STL surface has degenerate bounds");
    const auto map = [&](Point p) {
        return Point{0.15 + 0.70 * (p.x - lo.x) / span.x, 0.15 + 0.70 * (p.y - lo.y) / span.y,
                     0.15 + 0.70 * (p.z - lo.z) / span.z};
    };
    for (auto& triangle : surface.triangles) {
        triangle.a = map(triangle.a);
        triangle.b = map(triangle.b);
        triangle.c = map(triangle.c);
    }
}

struct Snapshot {
    std::vector<cfd::lbm::Vector3> velocity;
    std::vector<double> speed;
};

Snapshot snapshot(const cfd::lbm::D3Q19Solver& solver) {
    const auto fields = solver.compute_macroscopic();
    Snapshot result{std::vector<cfd::lbm::Vector3>(solver.cells()), std::vector<double>(solver.cells())};
    for (std::size_t i = 0; i < solver.cells(); ++i) {
        result.velocity[i] = {fields.ux[i], fields.uy[i], fields.uz[i]};
        result.speed[i] = std::sqrt(static_cast<double>(fields.ux[i]) * fields.ux[i] +
                                    static_cast<double>(fields.uy[i]) * fields.uy[i] +
                                    static_cast<double>(fields.uz[i]) * fields.uz[i]);
    }
    return result;
}

void write_frame(const std::filesystem::path& out, const cfd::lbm::D3Q19Solver& solver,
                 const std::vector<unsigned char>& solid, std::size_t frame) {
    const auto fields = snapshot(solver);
    std::ostringstream stem;
    stem << "frame_" << std::setw(6) << std::setfill('0') << frame;
    cfd::lbm::write_vtk_structured_3d(out / "vtk" / (stem.str() + ".vtk"), solver.config().nx, solver.config().ny,
                                      solver.config().nz, fields.speed, fields.velocity, solid, "speed");
    cfd::lbm::write_json_state_3d(out / "live" / "state.json", solver.config().nx, solver.config().ny,
                                  solver.config().nz, solver.time_step(), fields.speed, solid, "speed");
}

} // namespace

int main(int argc, char** argv) {
    using namespace cfd::examples;
    const auto opt = parse_options(argc, argv);
    const std::size_t nx = opt.quick ? 48U : 96U * opt.scale;
    const std::size_t ny = opt.quick ? 24U : 48U * opt.scale;
    const std::size_t nz = opt.quick ? 24U : 48U * opt.scale;
    const std::size_t steps = opt.quick ? 12U : 80U;
    const std::size_t snapshot_interval = opt.quick ? 3U : 10U;
    const auto out = output_dir(argc, argv);
    const auto external_stl = stl_path(argc, argv);

    Timer setup;
    Surface surface = external_stl.has_value() ? cfd::io::read_ascii_stl(*external_stl) : make_ahmed_body();
    if (external_stl.has_value())
        fit_surface_to_unit_box(surface);
    const cfd::lbm::VoxelGrid3D grid{nx, ny, nz, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
    const auto solid = cfd::lbm::voxelize_surface(surface, grid);
    cfd::lbm::InPlaceLbmConfig config{nx, ny, nz, 0.58F, 3.0e-6F, 0.0F, 0.0F};
    config.smagorinsky_les = true;
    cfd::lbm::D3Q19Solver solver(config);
    solver.set_solid_mask(solid);
    solver.initialize_uniform();
    const double setup_ms = setup.milliseconds();

    std::filesystem::create_directories(out / "vtk");
    std::filesystem::create_directories(out / "live");
    write_frame(out, solver, solid, 0U);

    Timer simulation;
    for (std::size_t step = 1; step <= steps; ++step) {
        solver.step();
        if (step % snapshot_interval == 0U || step == steps)
            write_frame(out, solver, solid,
                        step / snapshot_interval + (step == steps && step % snapshot_interval != 0U));
    }
    const double simulation_ms = simulation.milliseconds();

    const auto fields = snapshot(solver);
    const auto q = cfd::lbm::q_criterion_3d(nx, ny, nz, fields.velocity);
    cfd::lbm::write_vtk_structured_3d(out / "vtk" / "q_criterion_final.vtk", nx, ny, nz, q, {}, solid, "q_criterion");
    const auto solid_cells = static_cast<std::size_t>(std::count(solid.begin(), solid.end(), 1U));
    const double checksum = solver.mass() + solver.kinetic_energy() + static_cast<double>(solid_cells);
    std::cout << "Wrote 3-D VTK frames: " << (out / "vtk") << "\n";
    std::cout << "Live state for examples/live_3d_viewer.py: " << (out / "live" / "state.json") << "\n";
    std::cout << "Voxelized model cells: " << solid_cells << " / " << solid.size() << "\n";
    emit({"real_lbm3d_model", "voxelized_ahmed_body_live_vtk", "cpu", "cell_steps", nx * ny * nz, steps, setup_ms,
          simulation_ms, static_cast<double>(nx * ny * nz) * static_cast<double>(steps), checksum});
    return 0;
}
