#include "common/benchmark.hpp"
#include "cfd/io/openfoam_io.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::filesystem::path output_dir(int argc, char** argv) {
    std::filesystem::path out = "examples/output/openfoam_channel";
    for (int i = 1; i + 1 < argc; ++i) {
        if (std::string(argv[i]) == "--output")
            out = argv[++i];
    }
    return out;
}

void write_structured_vtk(const std::filesystem::path& path, std::size_t nx, std::size_t ny, double lx, double ly,
                          const std::vector<cfd::fvm::Vec3>& u, const std::vector<double>& p) {
    std::ofstream o(path);
    if (!o)
        throw std::runtime_error("cannot create VTK file: " + path.string());
    o << "# vtk DataFile Version 3.0\n"
      << "cfd_solvers ventilated channel\nASCII\n"
      << "DATASET STRUCTURED_POINTS\n"
      << "DIMENSIONS " << nx << ' ' << ny << " 1\n"
      << "ORIGIN " << (0.5 * lx / static_cast<double>(nx)) << ' ' << (0.5 * ly / static_cast<double>(ny)) << " 0\n"
      << std::setprecision(17) << "SPACING " << (lx / static_cast<double>(nx)) << ' ' << (ly / static_cast<double>(ny))
      << " 1\n"
      << "POINT_DATA " << u.size() << "\n"
      << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for (double x : p)
        o << x << '\n';
    o << "VECTORS velocity double\n";
    for (const auto& v : u)
        o << v.x << ' ' << v.y << ' ' << v.z << '\n';
}

struct OpenFoamFace {
    std::array<std::size_t, 4> points{};
    std::size_t owner{};
    std::size_t neighbour{};
    bool internal{};
};

void write_openfoam_case(const std::filesystem::path& case_dir, std::size_t nx, std::size_t ny, std::size_t nz,
                         double lx, double ly, double lz, const std::vector<cfd::fvm::Vec3>& velocity,
                         const std::vector<double>& pressure) {
    const auto cell_count = nx * ny * nz;
    if (nx == 0U || ny == 0U || nz == 0U || velocity.size() != cell_count || pressure.size() != cell_count)
        throw std::invalid_argument("OpenFOAM export field size does not match the channel mesh");
    const auto mesh_dir = case_dir / "constant" / "polyMesh";
    std::filesystem::create_directories(mesh_dir);
    std::filesystem::create_directories(case_dir / "system");

    const auto point_index = [=](std::size_t i, std::size_t j, std::size_t k) {
        return (k * (ny + 1U) + j) * (nx + 1U) + i;
    };
    const auto cell_index = [=](std::size_t i, std::size_t j, std::size_t k) { return (k * ny + j) * nx + i; };
    const auto point = [=](std::size_t i, std::size_t j, std::size_t k) {
        return cfd::fvm::Vec3{lx * static_cast<double>(i) / static_cast<double>(nx),
                              ly * static_cast<double>(j) / static_cast<double>(ny),
                              lz * static_cast<double>(k) / static_cast<double>(nz)};
    };

    std::vector<cfd::fvm::Vec3> points;
    points.reserve((nx + 1U) * (ny + 1U) * (nz + 1U));
    for (std::size_t k = 0; k <= nz; ++k)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t i = 0; i <= nx; ++i)
                points.push_back(point(i, j, k));

    std::vector<OpenFoamFace> faces;
    faces.reserve((nx + 1U) * ny * nz + nx * (ny + 1U) * nz + nx * ny * (nz + 1U));
    const auto add_internal = [&](std::array<std::size_t, 4> face, std::size_t owner, std::size_t neighbour) {
        faces.push_back({face, owner, neighbour, true});
    };
    const auto add_boundary = [&](std::array<std::size_t, 4> face, std::size_t owner) {
        faces.push_back({face, owner, 0U, false});
    };

    // OpenFOAM requires internal faces first, followed by contiguous boundary patches.
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 1; i < nx; ++i)
                add_internal({point_index(i, j, k), point_index(i, j + 1U, k), point_index(i, j + 1U, k + 1U),
                              point_index(i, j, k + 1U)},
                             cell_index(i - 1U, j, k), cell_index(i, j, k));
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t i = 0; i < nx; ++i)
            for (std::size_t j = 1; j < ny; ++j)
                add_internal({point_index(i, j, k), point_index(i, j, k + 1U), point_index(i + 1U, j, k + 1U),
                              point_index(i + 1U, j, k)},
                             cell_index(i, j - 1U, k), cell_index(i, j, k));
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            for (std::size_t k = 1; k < nz; ++k)
                add_internal({point_index(i, j, k), point_index(i + 1U, j, k), point_index(i + 1U, j + 1U, k),
                              point_index(i, j + 1U, k)},
                             cell_index(i, j, k - 1U), cell_index(i, j, k));

    std::array<std::size_t, 6> starts{};
    std::array<std::size_t, 6> counts{};
    const auto begin_patch = [&](std::size_t patch) { starts[patch] = faces.size(); };
    const auto finish_patch = [&](std::size_t patch) { counts[patch] = faces.size() - starts[patch]; };

    begin_patch(0U);
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            add_boundary({point_index(0U, j, k), point_index(0U, j, k + 1U), point_index(0U, j + 1U, k + 1U),
                          point_index(0U, j + 1U, k)},
                         cell_index(0U, j, k));
    finish_patch(0U);

    begin_patch(1U);
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            add_boundary({point_index(nx, j, k), point_index(nx, j + 1U, k), point_index(nx, j + 1U, k + 1U),
                          point_index(nx, j, k + 1U)},
                         cell_index(nx - 1U, j, k));
    finish_patch(1U);

    begin_patch(2U);
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t i = 0; i < nx; ++i)
            add_boundary({point_index(i, 0U, k), point_index(i + 1U, 0U, k), point_index(i + 1U, 0U, k + 1U),
                          point_index(i, 0U, k + 1U)},
                         cell_index(i, 0U, k));
    finish_patch(2U);

    begin_patch(3U);
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t i = 0; i < nx; ++i)
            add_boundary({point_index(i, ny, k), point_index(i, ny, k + 1U), point_index(i + 1U, ny, k + 1U),
                          point_index(i + 1U, ny, k)},
                         cell_index(i, ny - 1U, k));
    finish_patch(3U);

    begin_patch(4U);
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            add_boundary({point_index(i, j, 0U), point_index(i, j + 1U, 0U), point_index(i + 1U, j + 1U, 0U),
                          point_index(i + 1U, j, 0U)},
                         cell_index(i, j, 0U));
    finish_patch(4U);

    begin_patch(5U);
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            add_boundary({point_index(i, j, nz), point_index(i + 1U, j, nz), point_index(i + 1U, j + 1U, nz),
                          point_index(i, j + 1U, nz)},
                         cell_index(i, j, nz - 1U));
    finish_patch(5U);

    const auto write_header = [](std::ofstream& out, const char* type, const char* object) {
        out << "FoamFile\n{\n    version 2.0;\n    format ascii;\n    class " << type << ";\n    object " << object
            << ";\n}\n";
    };
    {
        std::ofstream out(mesh_dir / "points");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM points");
        write_header(out, "vectorField", "points");
        out << std::setprecision(17) << points.size() << "\n(\n";
        for (const auto& p : points)
            out << '(' << p.x << ' ' << p.y << ' ' << p.z << ")\n";
        out << ")\n";
    }
    {
        std::ofstream out(mesh_dir / "faces");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM faces");
        write_header(out, "faceList", "faces");
        out << faces.size() << "\n(\n";
        for (const auto& face : faces)
            out << "4 (" << face.points[0] << ' ' << face.points[1] << ' ' << face.points[2] << ' ' << face.points[3]
                << ")\n";
        out << ")\n";
    }
    {
        std::ofstream out(mesh_dir / "owner");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM owner");
        write_header(out, "labelList", "owner");
        out << faces.size() << "\n(\n";
        for (const auto& face : faces)
            out << face.owner << '\n';
        out << ")\n";
    }
    {
        std::ofstream out(mesh_dir / "neighbour");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM neighbour");
        write_header(out, "labelList", "neighbour");
        const auto internal_count = static_cast<std::size_t>(
            std::count_if(faces.begin(), faces.end(), [](const OpenFoamFace& face) { return face.internal; }));
        out << internal_count << "\n(\n";
        for (const auto& face : faces)
            if (face.internal)
                out << face.neighbour << '\n';
        out << ")\n";
    }
    {
        std::ofstream out(mesh_dir / "boundary");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM boundary");
        write_header(out, "polyBoundaryMesh", "boundary");
        constexpr std::array<const char*, 6> names = {"left", "right", "bottom", "top", "front", "back"};
        out << names.size() << "\n(\n";
        for (std::size_t patch = 0; patch < names.size(); ++patch)
            out << names[patch] << "\n{\n    type patch;\n    nFaces " << counts[patch] << ";\n    startFace "
                << starts[patch] << ";\n}\n";
        out << ")\n";
    }

    {
        std::ofstream out(case_dir / "system" / "controlDict");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM controlDict");
        write_header(out, "dictionary", "controlDict");
        out << "application cfd_solvers_export;\nstartFrom latestTime;\nstartTime 0;\nstopAt endTime;\nendTime 0;\n"
               "deltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n";
    }
    {
        std::ofstream out(case_dir / "constant" / "transportProperties");
        if (!out)
            throw std::runtime_error("cannot create OpenFOAM transportProperties");
        write_header(out, "dictionary", "transportProperties");
        out << "transportModel Newtonian;\nnu [0 2 -1 0 0 0 0] 0.015;\n";
    }

    const std::string velocity_boundary = "boundaryField\n{\n"
                                          "    left { type fixedValue; value uniform (1.5 0 0); }\n"
                                          "    right { type zeroGradient; }\n"
                                          "    bottom { type fixedValue; value uniform (0 0 0); }\n"
                                          "    top { type fixedValue; value uniform (0 0 0); }\n"
                                          "    front { type slip; }\n"
                                          "    back { type slip; }\n"
                                          "}\n";
    const std::string pressure_boundary = "boundaryField\n{\n"
                                          "    left { type zeroGradient; }\n"
                                          "    right { type fixedValue; value uniform 0; }\n"
                                          "    bottom { type zeroGradient; }\n"
                                          "    top { type zeroGradient; }\n"
                                          "    front { type zeroGradient; }\n"
                                          "    back { type zeroGradient; }\n"
                                          "}\n";
    cfd::io::write_openfoam_vol_vector_field(case_dir / "0" / "U", "U", velocity, "[0 1 -1 0 0 0 0]",
                                             velocity_boundary);
    cfd::io::write_openfoam_vol_scalar_field(case_dir / "0" / "p", "p", pressure, "[0 2 -2 0 0 0 0]",
                                             pressure_boundary);
}

} // namespace

int main(int argc, char** argv) {
    using namespace cfd::examples;
    const auto opt = parse_options(argc, argv);
    const std::size_t nx = opt.quick ? 16U : 64U * opt.scale;
    const std::size_t ny = opt.quick ? 6U : 20U * opt.scale;
    const std::size_t steps = opt.quick ? 1U : 12U;
    constexpr double length = 4.0;
    constexpr double height = 1.0;

    Timer setup;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(nx, ny, 1, length, height, 0.08);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.01;
    cfg.kinematic_viscosity = 0.015;
    cfg.include_convection = true;
    cfg.momentum_iterations = 160;
    cfg.pressure_iterations = 400;
    cfg.momentum_tolerance = 1e-8;
    cfg.pressure_tolerance = 1e-8;
    cfg.pressure_correctors = 2;
    cfg.nonorthogonal_correctors = 0;
    cfd::fvm::CollocatedIncompressible flow(std::move(mesh), cfg);
    flow.set_velocity_boundary("left", cfd::fvm::VelocityBoundaryType::fixedValue, {1.5, 0.0, 0.0});
    flow.set_velocity_boundary("right", cfd::fvm::VelocityBoundaryType::zeroGradient);
    flow.set_pressure_boundary("right", cfd::fvm::PressureBoundaryType::fixedValue, 0.0);
    flow.set_velocity_boundary("bottom", cfd::fvm::VelocityBoundaryType::fixedValue, {0.0, 0.0, 0.0});
    flow.set_velocity_boundary("top", cfd::fvm::VelocityBoundaryType::fixedValue, {0.0, 0.0, 0.0});
    flow.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    flow.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    flow.initialize_fields([](cfd::fvm::Vec3 x) { return cfd::fvm::Vec3{1.5 * 4.0 * x.y * (1.0 - x.y), 0.0, 0.0}; },
                           [](cfd::fvm::Vec3 x) { return 0.05 * (4.0 - x.x); });
    const double setup_ms = setup.milliseconds();

    Timer sim;
    for (std::size_t s = 0; s < steps; ++s)
        flow.step_piso();
    const double sim_ms = sim.milliseconds();

    const auto out = output_dir(argc, argv);
    std::filesystem::create_directories(out / "openfoam" / "0");
    std::filesystem::create_directories(out / "paraview");
    write_openfoam_case(out / "openfoam", nx, ny, 1U, length, height, 0.08, flow.velocity(), flow.pressure());
    write_structured_vtk(out / "paraview" / "channel_final.vtk", nx, ny, length, height, flow.velocity(),
                         flow.pressure());

    double checksum = flow.continuity_l2();
    for (const auto& v : flow.velocity())
        checksum += std::abs(v.x) + std::abs(v.y);
    std::cout << "Wrote ParaView VTK: " << (out / "paraview" / "channel_final.vtk") << "\n";
    std::cout << "Wrote OpenFOAM case: " << (out / "openfoam") << "\n";
    emit({"real_openfoam_channel", "ventilated_channel_piso", "cpu", "cell_steps", flow.mesh().cell_count(), steps,
          setup_ms, sim_ms, double(flow.mesh().cell_count()) * double(steps), checksum});
    return 0;
}
