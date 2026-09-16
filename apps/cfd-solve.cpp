#include "cfd/core/parallel.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fem/poisson1d.hpp"
#include "cfd/solvers/fvm/diffusion2d.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/optics/ray.hpp"
#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/d2q9_sycl.hpp"
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>

namespace {

void print_list() {
    std::cout
        << "Available solvers:\n"
        << "  lbm-d2q9-cpu           legacy two-grid periodic BGK D2Q9 baseline (OpenMP)\n"
        << "  lbm-d2q9-inplace-cpu   single-grid periodic BGK D2Q9 (OpenMP)\n"
        << "  lbm-d3q19-cpu          single-grid periodic BGK D3Q19 (OpenMP)\n"
        << "  lbm-d3q27-cpu          single-grid periodic BGK D3Q27 (OpenMP)\n"
        << "  lbm-poiseuille         force-driven channel with halfway bounce-back walls\n"
        << "  lbm-cavity             moving-wall lid-driven cavity\n"
        << "  lbm-channel-io         Zou-He velocity inlet + pressure outlet channel\n"
#if defined(CFD_HAS_SYCL)
        << "  lbm-d2q9-sycl          legacy two-grid D2Q9 (SYCL)\n"
        << "  lbm-d2q9-inplace-sycl  single-grid D2Q9 (SYCL)\n"
        << "  lbm-d3q19-sycl         single-grid D3Q19 (SYCL)\n"
        << "  lbm-d3q27-sycl         single-grid D3Q27 (SYCL)\n"
#endif
        << "  fvm-diffusion2d        cell-centered finite-volume diffusion/Poisson solver\n"
        << "  fvm-operators3d        polyhedral owner/neighbour Gauss-operator manufactured test\n"
        << "  fvm-projection2d       staggered periodic incompressible pressure projection\n"
        << "  fvm-taylor-green2d     transient periodic incompressible Navier-Stokes Taylor-Green case\n"
        << "  fvm-collocated-channel pressure-driven collocated SIMPLE Poiseuille channel\n"
        << "  fvm-collocated-cavity  collocated SIMPLE lid-driven cavity\n"
        << "  fvm-collocated-skew    sheared-mesh PISO/PIMPLE pressure-coupling regression\n"
        << "  fem-poisson1d          linear finite-element Poisson solver\n"
        << "  fdtd-maxwell1d         Yee-grid electromagnetic time-domain solver\n"
        << "  optics-snell           geometric-optics refraction/reflection demo\n";
}

int run_lbm_cpu_legacy() {
    cfd::lbm::D2Q9Solver solver({256, 256, 0.60F});
    solver.initialize_taylor_green();
    const double mass0 = solver.mass();
    const auto t0 = std::chrono::steady_clock::now();
    solver.step(200);
    const auto t1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    const double mlups = (static_cast<double>(256U * 256U) * 200.0) / (seconds * 1.0e6);
    std::cout << "backend=cpu streaming=two-grid threads=" << cfd::core::cpu_thread_capacity()
              << " mass_error=" << std::abs(solver.mass() - mass0)
              << " kinetic_energy=" << solver.kinetic_energy()
              << " max_speed=" << solver.max_speed()
              << " MLUPS=" << mlups << '\n';
    return 0;
}

template<class Descriptor>
int run_lbm_inplace_cpu(cfd::lbm::InPlaceLbmConfig config, std::size_t steps) {
    cfd::lbm::EsotericPullSolver<Descriptor> solver(config);
    solver.initialize_taylor_green(0.02F);
    const double mass0 = solver.mass();
    const auto t0 = std::chrono::steady_clock::now();
    solver.step(steps);
    const auto t1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    const double mlups = static_cast<double>(solver.cells()) * static_cast<double>(steps) /
                         (seconds * 1.0e6);
    const double mass_error = std::abs(solver.mass() - mass0);
    std::cout << "backend=cpu streaming=in-place threads=" << cfd::core::cpu_thread_capacity()
              << " cells=" << solver.cells()
              << " population_MiB=" << static_cast<double>(solver.population_bytes()) / (1024.0 * 1024.0)
              << " mass_error_abs=" << mass_error
              << " mass_error_rel=" << mass_error / std::abs(mass0)
              << " kinetic_energy=" << solver.kinetic_energy()
              << " max_speed=" << solver.max_speed()
              << " MLUPS=" << mlups << '\n';
    return 0;
}

#if defined(CFD_HAS_SYCL)
int run_lbm_sycl_legacy() {
    cfd::lbm::D2Q9SyclSolver solver({256, 256, 0.60F});
    solver.initialize_taylor_green();
    const auto t0 = std::chrono::steady_clock::now();
    solver.step(200);
    solver.wait();
    const auto t1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    const double mlups = (static_cast<double>(256U * 256U) * 200.0) / (seconds * 1.0e6);
    const auto rho = solver.download_density();
    double mass = 0.0;
    for (float r : rho) mass += r;
    std::cout << "backend=sycl streaming=two-grid device=\"" << solver.device_name() << "\""
              << " mass=" << mass << " MLUPS=" << mlups << '\n';
    return 0;
}

template<class Descriptor>
int run_lbm_inplace_sycl(cfd::lbm::InPlaceLbmConfig config, std::size_t steps) {
    cfd::lbm::EsotericPullSyclSolver<Descriptor> solver(config);
    solver.initialize_taylor_green(0.02F);
    const auto before = solver.download_macroscopic();
    double mass0 = 0.0;
    for (float rho : before.rho) mass0 += rho;
    const auto t0 = std::chrono::steady_clock::now();
    solver.step(steps);
    solver.wait();
    const auto t1 = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(t1 - t0).count();
    const auto after = solver.download_macroscopic();
    double mass1 = 0.0;
    for (float rho : after.rho) mass1 += rho;
    const double mlups = static_cast<double>(solver.cells()) * static_cast<double>(steps) /
                         (seconds * 1.0e6);
    std::cout << "backend=sycl streaming=in-place device=\"" << solver.device_name() << "\""
              << " cells=" << solver.cells()
              << " population_MiB=" << static_cast<double>(solver.population_bytes()) / (1024.0 * 1024.0)
              << " mass_error_abs=" << std::abs(mass1 - mass0)
              << " mass_error_rel=" << std::abs(mass1 - mass0) / std::abs(mass0)
              << " MLUPS=" << mlups << '\n';
    return 0;
}
#endif

int run_poiseuille() {
    constexpr std::size_t nx = 96;
    constexpr std::size_t ny = 32;
    constexpr float tau = 0.8F;
    constexpr float acceleration = 5.0e-6F;
    cfd::lbm::D2Q9Solver solver({nx, ny, tau, acceleration, 0.0F});
    for (std::size_t x = 0; x < nx; ++x) {
        solver.set_solid(x, 0);
        solver.set_solid(x, ny - 1);
    }
    solver.initialize_uniform();
    solver.step(12000);
    const double nu = (static_cast<double>(tau) - 0.5) / 3.0;
    const double height = static_cast<double>(ny - 2);
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t y = 1; y + 1 < ny; ++y) {
        double ux = 0.0;
        for (std::size_t x = 0; x < nx; ++x) ux += solver.velocity_x()[y * nx + x];
        ux /= static_cast<double>(nx);
        const double wall_y = static_cast<double>(y) - 0.5;
        const double exact = static_cast<double>(acceleration) * wall_y * (height - wall_y) / (2.0 * nu);
        error2 += (ux - exact) * (ux - exact);
        exact2 += exact * exact;
    }
    std::cout << "case=poiseuille relative_L2=" << std::sqrt(error2 / exact2)
              << " max_speed=" << solver.max_speed() << '\n';
    return 0;
}

int run_cavity() {
    constexpr std::size_t n = 96;
    cfd::lbm::D2Q9Solver solver({n, n, 0.8F});
    for (std::size_t x = 0; x < n; ++x) {
        solver.set_solid(x, 0);
        solver.set_wall_velocity(x, n - 1, 0.05F, 0.0F);
    }
    for (std::size_t y = 0; y < n; ++y) {
        solver.set_solid(0, y);
        solver.set_solid(n - 1, y);
    }
    solver.initialize_uniform();
    solver.step(12000);
    const std::size_t center = (n / 2) * n + n / 2;
    std::cout << "case=lid-driven-cavity center_ux=" << solver.velocity_x()[center]
              << " center_uy=" << solver.velocity_y()[center]
              << " max_speed=" << solver.max_speed() << '\n';
    return 0;
}

int run_channel_io() {
    constexpr std::size_t nx = 160;
    constexpr std::size_t ny = 48;
    cfd::lbm::D2Q9Solver solver({nx, ny, 0.8F});
    for (std::size_t x = 0; x < nx; ++x) {
        solver.set_solid(x, 0);
        solver.set_solid(x, ny - 1);
    }
    solver.set_velocity_inlet_left(0.02F);
    solver.set_pressure_outlet_right(1.0F);
    solver.initialize_uniform(1.0F, 0.02F, 0.0F);
    solver.step(8000);
    std::cout << "case=velocity-pressure-channel"
              << " inlet_ux=" << solver.average_ux(0, 1, 1, ny - 1)
              << " mid_ux=" << solver.average_ux(nx / 2, nx / 2 + 1, 1, ny - 1)
              << " outlet_ux=" << solver.average_ux(nx - 1, nx, 1, ny - 1)
              << " max_speed=" << solver.max_speed() << '\n';
    return 0;
}

int run_fvm_operators() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(24, 20, 16, 1.2, 1.0, 0.8);
    std::vector<double> phi(mesh.cell_count());
    std::vector<double> boundary(mesh.face_count(), 0.0);
    for (std::size_t c = 0; c < mesh.cell_count(); ++c) {
        const auto x = mesh.cells()[c].center;
        phi[c] = x.x + 2.0 * x.y + 3.0 * x.z;
    }
    for (std::size_t f = 0; f < mesh.face_count(); ++f) {
        if (!mesh.faces()[f].boundary()) continue;
        const auto x = mesh.faces()[f].center;
        boundary[f] = x.x + 2.0 * x.y + 3.0 * x.z;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto grad = cfd::fvm::gauss_gradient_scalar(mesh, phi, boundary);
    const auto t1 = std::chrono::steady_clock::now();
    double max_error = 0.0;
    for (const auto g : grad) {
        max_error = std::max(max_error, std::abs(g.x - 1.0));
        max_error = std::max(max_error, std::abs(g.y - 2.0));
        max_error = std::max(max_error, std::abs(g.z - 3.0));
    }
    std::cout << "case=fvm-operators3d cells=" << mesh.cell_count()
              << " faces=" << mesh.face_count()
              << " threads=" << cfd::core::cpu_thread_capacity()
              << " max_gradient_error=" << max_error
              << " seconds=" << std::chrono::duration<double>(t1 - t0).count() << '\n';
    return 0;
}

int run_fvm_projection() {
    cfd::fvm::Projection2D solver({96, 80, 1.0, 1.0, 1.0, 1.0e-3, 1200});
    solver.initialize_divergent(0.05);
    const double before = solver.divergence_l2();
    const auto t0 = std::chrono::steady_clock::now();
    solver.project();
    const auto t1 = std::chrono::steady_clock::now();
    const double after = solver.divergence_l2();
    const auto& pressure = solver.pressure_result();
    std::cout << "case=fvm-projection2d threads=" << cfd::core::cpu_thread_capacity()
              << " divergence_before=" << before
              << " divergence_after=" << after
              << " reduction=" << (after / before)
              << " pressure_iterations=" << pressure.iterations
              << " pressure_residual_rms=" << pressure.residual_rms
              << " pressure_converged=" << (pressure.converged ? 1 : 0)
              << " seconds=" << std::chrono::duration<double>(t1 - t0).count() << '\n';
    return 0;
}

int run_fvm_taylor_green() {
    cfd::fvm::Incompressible2D solver({64, 64, 1.0, 1.0, 1.0, 1.0e-2, 1.0e-4, 800, 1.0e-10});
    constexpr double amplitude = 0.05;
    solver.initialize_taylor_green(amplitude);
    const double energy0 = solver.kinetic_energy();
    const auto t0 = std::chrono::steady_clock::now();
    solver.run(100);
    const auto t1 = std::chrono::steady_clock::now();
    const auto& pressure = solver.pressure_result();
    std::cout << "case=fvm-taylor-green2d threads=" << cfd::core::cpu_thread_capacity()
              << " nx=" << solver.config().nx
              << " ny=" << solver.config().ny
              << " steps=" << solver.steps()
              << " time=" << solver.time()
              << " divergence=" << solver.divergence_l2()
              << " velocity_rms_error=" << solver.taylor_green_velocity_error(amplitude)
              << " energy_ratio=" << solver.kinetic_energy()/energy0
              << " pressure_iterations=" << pressure.iterations
              << " pressure_converged=" << (pressure.converged ? 1 : 0)
              << " seconds=" << std::chrono::duration<double>(t1-t0).count() << '\n';
    return pressure.converged ? 0 : 2;
}


int run_fvm_collocated_channel() {
    constexpr double length = 2.0;
    constexpr double height = 1.0;
    constexpr double pressure_drop = 0.01;
    constexpr double nu = 0.05;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(32, 16, 1, length, height, 1.0);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.2;
    cfg.kinematic_viscosity = nu;
    cfg.include_convection = false;
    cfg.momentum_sweeps = 4;
    cfg.nonorthogonal_correctors = 0;
    cfg.velocity_relaxation = 0.8;
    cfg.pressure_relaxation = 0.5;
    cfg.pressure_tolerance = 1.0e-10;
    cfd::fvm::CollocatedIncompressible solver(std::move(mesh), cfg);
    solver.set_velocity_boundary("bottom", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("top", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("left", cfd::fvm::VelocityBoundaryType::zeroGradient);
    solver.set_velocity_boundary("right", cfd::fvm::VelocityBoundaryType::zeroGradient);
    solver.set_pressure_boundary("left", cfd::fvm::PressureBoundaryType::fixedValue, pressure_drop);
    solver.set_pressure_boundary("right", cfd::fvm::PressureBoundaryType::fixedValue, 0.0);
    solver.initialize_fields(
        [](cfd::fvm::Vec3) { return cfd::fvm::Vec3{}; },
        [&](cfd::fvm::Vec3 x) { return pressure_drop * (1.0 - x.x / length); });
    const auto t0 = std::chrono::steady_clock::now();
    const auto info = solver.solve_simple(300, 2.0e-7, 1.0e-6);
    const auto t1 = std::chrono::steady_clock::now();
    double error2 = 0.0;
    double exact2 = 0.0;
    double max_u = 0.0;
    std::size_t count = 0;
    const double gradient = pressure_drop / length;
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const auto x = solver.mesh().cells()[c].center;
        max_u = std::max(max_u, solver.velocity()[c].x);
        if (x.x < 0.8 || x.x > 1.2) continue;
        const double exact = gradient * x.y * (height - x.y) / (2.0 * nu);
        const double e = solver.velocity()[c].x - exact;
        error2 += e * e;
        exact2 += exact * exact;
        ++count;
    }
    std::cout << "case=fvm-collocated-channel"
              << " threads=" << cfd::core::cpu_thread_capacity()
              << " profile_relative_L2=" << std::sqrt(error2 / exact2)
              << " max_u=" << max_u
              << " continuity_L2=" << info.continuity_l2
              << " continuity_max=" << info.continuity_max
              << " velocity_change=" << info.velocity_rms_change
              << " pressure_iterations=" << info.pressure.iterations
              << " pressure_converged=" << (info.pressure.converged ? 1 : 0)
              << " samples=" << count
              << " seconds=" << std::chrono::duration<double>(t1 - t0).count() << '\n';
    return info.pressure.converged ? 0 : 2;
}

int run_fvm_collocated_cavity() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(20, 20, 1, 1.0, 1.0, 1.0);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.02;
    cfg.kinematic_viscosity = 0.05;
    cfg.include_convection = true;
    cfg.momentum_sweeps = 3;
    cfg.nonorthogonal_correctors = 0;
    cfg.velocity_relaxation = 0.7;
    cfg.pressure_relaxation = 0.4;
    cfg.pressure_tolerance = 1.0e-9;
    cfd::fvm::CollocatedIncompressible solver(std::move(mesh), cfg);
    solver.set_velocity_boundary("left", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("right", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("bottom", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("top", cfd::fvm::VelocityBoundaryType::fixedValue, {1.0, 0.0, 0.0});
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    solver.initialize_uniform();
    const auto t0 = std::chrono::steady_clock::now();
    const auto info = solver.solve_simple(300, 2.0e-6, 2.0e-6);
    const auto t1 = std::chrono::steady_clock::now();
    double center_distance = 1.0e9;
    cfd::fvm::Vec3 center{};
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const auto x = solver.mesh().cells()[c].center;
        const double d2 = (x.x - 0.5) * (x.x - 0.5) + (x.y - 0.5) * (x.y - 0.5);
        if (d2 < center_distance) {
            center_distance = d2;
            center = solver.velocity()[c];
        }
    }
    std::cout << "case=fvm-collocated-cavity"
              << " center_ux=" << center.x
              << " center_uy=" << center.y
              << " continuity_L2=" << info.continuity_l2
              << " continuity_max=" << info.continuity_max
              << " velocity_change=" << info.velocity_rms_change
              << " pressure_converged=" << (info.pressure.converged ? 1 : 0)
              << " seconds=" << std::chrono::duration<double>(t1 - t0).count() << '\n';
    return info.pressure.converged ? 0 : 2;
}

int run_fvm_collocated_skew() {
    auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(12, 10, 1, 1.0, 1.0, 1.0, 0.35);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.02;
    cfg.kinematic_viscosity = 0.02;
    cfg.include_convection = false;
    cfg.momentum_sweeps = 3;
    cfg.pressure_correctors = 2;
    cfg.outer_correctors = 2;
    cfg.nonorthogonal_correctors = 3;
    cfg.pressure_tolerance = 1.0e-10;
    cfd::fvm::CollocatedIncompressible solver(std::move(mesh), cfg);
    for (const auto name : {"left", "right", "bottom", "top"}) {
        solver.set_velocity_boundary(name, cfd::fvm::VelocityBoundaryType::fixedValue);
    }
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    auto velocity = [](cfd::fvm::Vec3 x) {
        return cfd::fvm::Vec3{0.1 * (x.x - 0.5), -0.08 * (x.y - 0.5), 0.0};
    };
    solver.initialize_fields(velocity, [](cfd::fvm::Vec3) { return 0.0; });
    const double before = solver.continuity_l2();
    const auto piso = solver.step_piso();
    solver.initialize_fields(velocity, [](cfd::fvm::Vec3) { return 0.0; });
    const auto pimple = solver.step_pimple();
    std::cout << "case=fvm-collocated-skew"
              << " continuity_before=" << before
              << " piso_continuity=" << piso.continuity_l2
              << " pimple_continuity=" << pimple.continuity_l2
              << " piso_pressure_converged=" << (piso.pressure.converged ? 1 : 0)
              << " pimple_pressure_converged=" << (pimple.pressure.converged ? 1 : 0)
              << '\n';
    return (piso.pressure.converged && pimple.pressure.converged) ? 0 : 2;
}

int run_fvm() {
    cfd::fvm::Diffusion2D solver({128, 128, 1.0, 1.0, 1.0});
    solver.set_dirichlet(0.0, 0.0, 0.0, 1.0);
    double residual = 0.0;
    for (int block = 0; block < 100; ++block) residual = solver.iterate(100);
    const auto& f = solver.field();
    std::cout << "residual=" << residual << " center=" << f[solver.grid().index(64, 64)] << '\n';
    return 0;
}

int run_fem() {
    cfd::fem::Poisson1D solver(128);
    solver.solve([](double) { return 1.0; });
    const auto& u = solver.u();
    const auto& x = solver.x();
    double max_error = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        const double exact = 0.5 * x[i] * (1.0 - x[i]);
        max_error = std::max(max_error, std::abs(u[i] - exact));
    }
    std::cout << "max_error=" << max_error << '\n';
    return 0;
}

int run_fdtd() {
    cfd::fdtd::Maxwell1D solver({2048, 1.0e-3, 0.99, 1.0, 1.0});
    solver.initialize_gaussian();
    const double e0 = solver.energy();
    solver.step(500);
    std::cout << "dt=" << solver.dt() << " energy_initial=" << e0
              << " energy_final=" << solver.energy() << '\n';
    return 0;
}

int run_optics() {
    const cfd::optics::Vec3 incident = cfd::optics::normalized({0.3, 0.0, 1.0});
    const auto transmitted = cfd::optics::refract(incident, {0.0, 0.0, -1.0}, 1.0, 1.5);
    if (!transmitted) return 2;
    std::cout << "refracted=(" << transmitted->x << ',' << transmitted->y << ',' << transmitted->z << ")\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || std::string_view(argv[1]) == "--list") {
        print_list();
        return 0;
    }
    const std::string_view solver = argv[1];
    if (solver == "lbm-d2q9-cpu") return run_lbm_cpu_legacy();
    if (solver == "lbm-d2q9-inplace-cpu") {
        return run_lbm_inplace_cpu<cfd::lbm::D2Q9InPlaceDescriptor>({256, 256, 1, 0.60F}, 200);
    }
    if (solver == "lbm-d3q19-cpu") {
        return run_lbm_inplace_cpu<cfd::lbm::D3Q19Descriptor>({64, 64, 64, 0.60F}, 100);
    }
    if (solver == "lbm-d3q27-cpu") {
        return run_lbm_inplace_cpu<cfd::lbm::D3Q27Descriptor>({56, 56, 56, 0.60F}, 100);
    }
    if (solver == "lbm-poiseuille") return run_poiseuille();
    if (solver == "lbm-cavity") return run_cavity();
    if (solver == "lbm-channel-io") return run_channel_io();
#if defined(CFD_HAS_SYCL)
    if (solver == "lbm-d2q9-sycl") return run_lbm_sycl_legacy();
    if (solver == "lbm-d2q9-inplace-sycl") {
        return run_lbm_inplace_sycl<cfd::lbm::D2Q9InPlaceDescriptor>({256, 256, 1, 0.60F}, 200);
    }
    if (solver == "lbm-d3q19-sycl") {
        return run_lbm_inplace_sycl<cfd::lbm::D3Q19Descriptor>({64, 64, 64, 0.60F}, 100);
    }
    if (solver == "lbm-d3q27-sycl") {
        return run_lbm_inplace_sycl<cfd::lbm::D3Q27Descriptor>({56, 56, 56, 0.60F}, 100);
    }
#endif
    if (solver == "fvm-diffusion2d") return run_fvm();
    if (solver == "fvm-operators3d") return run_fvm_operators();
    if (solver == "fvm-projection2d") return run_fvm_projection();
    if (solver == "fvm-taylor-green2d") return run_fvm_taylor_green();
    if (solver == "fvm-collocated-channel") return run_fvm_collocated_channel();
    if (solver == "fvm-collocated-cavity") return run_fvm_collocated_cavity();
    if (solver == "fvm-collocated-skew") return run_fvm_collocated_skew();
    if (solver == "fem-poisson1d") return run_fem();
    if (solver == "fdtd-maxwell1d") return run_fdtd();
    if (solver == "optics-snell") return run_optics();
    std::cerr << "Unknown solver: " << solver << "\n";
    print_list();
    return 1;
}
