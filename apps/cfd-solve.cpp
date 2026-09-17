#include "continuation_cases.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/solvers/electrochemistry/corrosion1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck_poly.hpp"
#include "cfd/solvers/electrochemistry/mixed_potential.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fdtd/maxwell3d.hpp"
#include "cfd/solvers/fdtd/port.hpp"
#include "cfd/solvers/fdtd/vtk_io.hpp"
#include "cfd/solvers/fem/poisson1d.hpp"
#include "cfd/solvers/fem/poisson3d.hpp"
#include "cfd/solvers/fem/nonlinear_poisson2d.hpp"
#include "cfd/solvers/fem/darcy2d.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include "cfd/solvers/fem/modal_bar1d.hpp"
#include "cfd/solvers/fvm/diffusion2d.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/optics/ray.hpp"
#include "cfd/solvers/optics/sequential.hpp"
#include "cfd/solvers/optics/paraxial.hpp"
#include "cfd/solvers/optics/materials.hpp"
#include "cfd/solvers/optics/polarization.hpp"
#include "cfd/multiphysics/field_registry.hpp"
#include "cfd/multiphysics/transfer.hpp"
#include "cfd/multiphysics/partitioned.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"
#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/d2q9_sycl.hpp"
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <vector>
#include <iostream>
#include <limits>
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
        << "  fvm-scalar-transport   implicit PolyMesh scalar advection-diffusion demo\n"
        << "  fem-poisson1d          linear finite-element Poisson solver\n"
        << "  fem-poisson3d          Tet4 finite-element manufactured Poisson case\n"
        << "  fem-nonlinear-poisson  nonlinear Tri3 Newton/ILU-GMRES manufactured case\n"
        << "  fem-darcy              saturated porous Darcy-flow channel\n"
        << "  fem-magnetostatic      2-D magnetic vector-potential manufactured case\n"
        << "  fem-modal-bar          fixed-free bar generalized eigenfrequency case\n"
        << "  fdtd-maxwell1d         Yee-grid electromagnetic time-domain solver\n"
        << "  fdtd-mur1d             dielectric/lossy 1-D FDTD with Mur absorbing boundaries\n"
        << "  fdtd-dispersive1d      Debye/Drude/Lorentz ADE material smoke case\n"
        << "  fdtd-port-vtk          PMC + synthetic S-parameters + VTK export\n"
        << "  fdtd-cpml-tfsf         CPML attenuation + +x total/scattered-field injection\n"
        << "  fdtd-lumped-rlc        field-coupled parallel R/L/C element smoke case\n"
        << "  multiphysics-coupling field registry, conservative transfer, Aitken fixed point\n"
        << "  multiphysics-electrothermal DC conduction -> Joule heat -> transient FEM thermal\n"
        << "  optics-snell           geometric-optics refraction/reflection demo\n"
        << "  optics-lens            sequential/paraxial lens + coating/material demo\n"
        << "  electrochem-corrosion1d ohmic electrolyte + Butler-Volmer corrosion cell\n"
        << "  electrochem-pnp1d      conservative 1-D Nernst-Planck transport demo\n"
        << "  electrochem-pnp-poly   PolyMesh Poisson-Nernst-Planck charged-slab demo\n"
        << "  electrochem-galvanic   multi-reaction galvanic mixed-potential demo\n"
        << "  em-frequency1d         driven 1-D frequency-domain Maxwell + PEC cavity modes\n"
        << "  em-edge2d              Tri3 Nedelec edge-element driven Maxwell baseline\n"
        << "  particle-pic1d         Boris charged-particle + electrostatic PIC baseline\n"
        << "  particle-em-pic1d      periodic 1-D/3-V electromagnetic PIC baseline\n"
        << "  particle-pic2d         periodic 2-D electrostatic PIC foundation\n"
        << "  particle-em-pic2d      periodic 2-D/3-V electromagnetic PIC baseline\n"
        << "  particle-staggered-em-pic2d true 2-D Yee EM-PIC with wall/sponge boundaries\n"
        << "  particle-pic3d         periodic 3-D electrostatic PIC foundation\n"
        << "  particle-em-pic3d      periodic 3-D/3-V electromagnetic PIC baseline\n"
        << "  particle-staggered-em-pic3d true 3-D Yee EM-PIC with wall/sponge boundaries\n"
        << "  particle-local-current3d local finite-volume 3-D current deposition check\n"
        << "  particle-sort-halo3d particle sorting and guard-halo classification check\n"
        << "  particle-domain-exchange3d serial PIC domain decomposition and scalar guard exchange check\n"
        << "  particle-comm-exchange3d communicator-ready particle/guard message packing check\n"
        << "  particle-transport3d communicator-shaped in-memory PIC transport check\n"
        << "  particle-serialized-transport3d serialized PIC transport envelope check\n"
        << "  particle-distributed-round3d serialized rank-local distributed PIC exchange round\n"
        << "  particle-field-guards3d six-component EM field guard exchange check\n"
        << "  particle-distributed-step3d serialized distributed staggered PIC step bridge\n"
        << "  particle-geant4-transport Geant4-style track/process/step transport baseline\n"
        << "  particle-transport-dose-bvh stochastic transport, BVH lookup, hit scoring, dose/SAR projection\n"
        << "  particle-campaign-doe SU2/csauto-style DOE, adapter, registry, finite-difference gradient\n"
        << "  particle-campaign-execution persistent campaign folders, external-solver command plans, log parsers, optimizer\n"
        << "  particle-campaign-local-runner native process launcher, runtime doctor, output discovery, registry update\n"
        << "  particle-campaign-control live-control directives, scheduler polling parse, output refresh\n"
        << "  particle-multiserver-deploy multi-server placement and Docker/Compose deploy scaffold\n"
        << "  particle-multiserver-execution supervised remote command plans and job metadata\n"
        << "  particle-multiserver-supervision retry, access-probe, log-tail and dashboard metadata\n"
        << "  particle-multiserver-controller dependency-free controller API scaffold over supervision artifacts\n"
        << "  particle-multiserver-dashboard controller static dashboard and event snapshot scaffold\n"
        << "  particle-plasma-chemistry multi-species plasma mass-action chemistry check\n"
        << "  particle-breakdown-threshold Paschen/multipactor threshold utility check\n"
        << "  bioheat-sar            SAR -> implicit Pennes bioheat baseline\n"
        << "  rf-dipole              sinusoidal half-wave dipole radiation baseline\n"
        << "  rf-microstrip          quasi-static microstrip impedance baseline\n"
        << "  rf-multiwire           coupled parallel thin-wire MoM + series-load baseline\n"
        << "  spice-rc               MNA small-signal RC low-pass baseline\n"
        << "  spice-diode            nonlinear diode DC operating-point baseline\n"
        << "  spice-adaptive         adaptive LTE-controlled RC transient baseline\n"
        << "  spice-pss-pz           periodic steady-state + pole/zero RC baseline\n";
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

int run_fem_poisson3d() {
    constexpr double pi=3.1415926535897932384626433832795;
    cfd::fem::Poisson3D solver(cfd::fem::make_box_tet_mesh(10,10,10));
    solver.solve([](cfd::fem::Point3 p){return 3.0*pi*pi*std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z);},
        [](cfd::fem::Point3){return 0.0;});
    double e2=0.0,x2=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){const auto p=solver.mesh().nodes[i];const double exact=std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z);const double e=solver.solution()[i]-exact;e2+=e*e;x2+=exact*exact;}
    std::cout << "case=fem-poisson3d nodes=" << solver.mesh().node_count()
              << " elements=" << solver.mesh().element_count()
              << " relative_L2=" << std::sqrt(e2/x2)
              << " iterations=" << solver.linear_result().iterations
              << " converged=" << (solver.linear_result().converged?1:0) << '\n';
    return solver.linear_result().converged?0:2;
}

int run_fem_nonlinear_poisson() {
    const double pi=std::numbers::pi_v<double>;
    constexpr double beta=3.0;
    cfd::fem::NonlinearPoisson2DConfig cfg; cfg.cubic_coefficient=beta;
    cfd::fem::NonlinearPoisson2D solver(cfd::fem::make_rectangle_tri_mesh(32,32),cfg);
    const auto exact=[&](cfd::fem::Node2 p){return std::sin(pi*p.x)*std::sin(pi*p.y);};
    solver.solve([&](cfd::fem::Node2 p){const double u=exact(p);return 2.0*pi*pi*u+beta*u*u*u;},
                 [](cfd::fem::Node2){return 0.0;});
    double e2=0.0,x2=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){const double ex=exact(solver.mesh().nodes[i]);const double e=solver.solution()[i]-ex;e2+=e*e;x2+=ex*ex;}
    std::cout << "case=fem-nonlinear-poisson relative_L2=" << std::sqrt(e2/x2)
              << " newton_iterations=" << solver.nonlinear_result().iterations
              << " residual_rms=" << solver.nonlinear_result().residual_rms
              << " converged=" << (solver.nonlinear_result().converged?1:0) << '\n';
    return solver.nonlinear_result().converged?0:2;
}

int run_fem_darcy() {
    using Type=cfd::fem::ScalarBoundaryType;
    constexpr double permeability=2.0e-12,viscosity=1.0e-3,p_left=1.0e5,length=2.0;
    cfd::fem::Darcy2DConfig cfg;cfg.viscosity=viscosity;
    cfd::fem::Darcy2D solver(cfd::fem::make_rectangle_tri_mesh(40,20,length,1.0),cfg);
    solver.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return p_left;},{}});
    solver.set_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    solver.solve([](cfd::fem::Node2){return permeability;});
    const auto v=solver.element_velocity([](cfd::fem::Node2){return permeability;});
    double mean=0.0;for(const auto&a:v)mean+=a.x;mean/=static_cast<double>(v.size());
    const double exact=permeability/viscosity*p_left/length;
    std::cout << "case=fem-darcy mean_vx=" << mean << " exact_vx=" << exact
              << " relative_error=" << std::abs(mean-exact)/exact
              << " iterations=" << solver.linear_result().iterations << '\n';
    return 0;
}

int run_fem_magnetostatic() {
    using Type=cfd::fem::ScalarBoundaryType;
    const double pi=std::numbers::pi_v<double>,nu=2.5;
    cfd::fem::Magnetostatics2D solver(cfd::fem::make_rectangle_tri_mesh(32,32));
    for(int p=0;p<4;++p)solver.set_boundary(p,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    const auto exact=[&](cfd::fem::Node2 x){return std::sin(pi*x.x)*std::sin(pi*x.y);};
    solver.solve([&](cfd::fem::Node2){return nu;},[&](cfd::fem::Node2 x){return 2.0*pi*pi*nu*exact(x);});
    double e2=0.0,x2=0.0;for(std::size_t i=0;i<solver.mesh().node_count();++i){const double ex=exact(solver.mesh().nodes[i]);const double e=solver.vector_potential()[i]-ex;e2+=e*e;x2+=ex*ex;}
    std::cout << "case=fem-magnetostatic relative_L2=" << std::sqrt(e2/x2)
              << " iterations=" << solver.linear_result().iterations
              << " converged=" << (solver.linear_result().converged?1:0) << '\n';
    return solver.linear_result().converged?0:2;
}

int run_fem_modal_bar() {
    cfd::fem::BarModal1DConfig cfg;cfg.elements=100;cfg.length=1.7;cfg.young_modulus=70.0e9;cfg.density=2700.0;
    cfg.eigen.relative_tolerance=2.0e-8;
    const cfd::fem::BarModal1D solver(cfg);
    const auto modes=solver.solve(3U);
    std::cout << "case=fem-modal-bar";
    for(std::size_t i=0;i<modes.size();++i)std::cout << " f" << (i+1U) << "_Hz=" << modes[i].frequency_hz;
    std::cout << '\n';
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


int run_scalar_transport() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(80, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ScalarTransport solver(std::move(mesh), {0.02, 0.05, 300, 20, 1.0e-11});
    solver.set_boundary("left", cfd::fvm::ScalarBoundaryType::fixedValue, 1.0);
    solver.set_boundary("right", cfd::fvm::ScalarBoundaryType::fixedValue, 0.0);
    solver.initialize(0.5);
    solver.run(120);
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const double exact = 1.0 - solver.mesh().cells()[c].center.x;
        const double error = solver.values()[c] - exact;
        error2 += error * error;
        exact2 += exact * exact;
    }
    std::cout << "case=fvm-scalar-transport"
              << " time=" << solver.time()
              << " relative_L2=" << std::sqrt(error2 / exact2)
              << " min=" << solver.minimum()
              << " max=" << solver.maximum()
              << " linear_iterations=" << solver.linear_result().iterations
              << " converged=" << solver.linear_result().converged << '\n';
    return solver.linear_result().converged ? 0 : 1;
}

int run_fdtd_mur1d() {
    cfd::fdtd::Maxwell1D solver({800,1.0e-3,0.95,1.0,1.0,cfd::fdtd::Boundary1D::mur1});
    solver.set_material(460,620,4.0,0.01);
    solver.initialize_gaussian(0.22,0.025);
    const double e0=solver.energy();
    solver.step(1200);
    std::cout << "case=fdtd-mur1d energy_ratio=" << solver.energy()/e0
              << " time_s=" << solver.time()
              << " dielectric_wave_speed=" << solver.local_wave_speed(500U) << '\n';
    return std::isfinite(solver.energy())?0:2;
}


int run_fdtd_dispersive1d() {
    constexpr std::size_t cells = 320U;
    cfd::fdtd::Maxwell1D solver({cells, 1.0e-4, 0.45, 1.0, 1.0, cfd::fdtd::Boundary1D::mur1});
    const double inv_dt = 1.0 / solver.dt();
    solver.set_debye_material(90U, 150U, 2.0, 3.0, 30.0 * solver.dt());
    solver.set_drude_material(150U, 210U, 1.0, 0.06 * inv_dt, 0.03 * inv_dt);
    solver.set_lorentz_material(210U, 270U, 1.5, 1.5, 0.055 * inv_dt, 0.025 * inv_dt);
    solver.initialize_gaussian(0.18, 0.025);
    const double initial = solver.energy();
    solver.step(220U);
    double pmax = 0.0;
    for (double p : solver.polarization()) pmax = std::max(pmax, std::abs(p));
    std::cout << "case=fdtd-dispersive1d"
              << " energy_ratio=" << solver.energy() / initial
              << " polarization_max=" << pmax
              << " finite=" << (std::isfinite(solver.energy()) ? 1 : 0) << '\n';
    return std::isfinite(solver.energy()) && pmax > 0.0 ? 0 : 2;
}

int run_fdtd_port_vtk() {
    constexpr double frequency = 25.0;
    constexpr double impedance = 50.0;
    constexpr double reflection = 0.2;
    constexpr double transmission = 0.75;
    constexpr double sample_dt = 5.0e-4;
    cfd::fdtd::WavePort1D reference(frequency, impedance, +1);
    cfd::fdtd::WavePort1D transmitted(frequency, impedance, +1);
    for (std::size_t n = 0; n <= 1600U; ++n) {
        const double t = sample_dt * static_cast<double>(n);
        const double incident = std::sin(2.0 * std::numbers::pi * frequency * t);
        const double reflected = reflection * incident;
        reference.sample(t, incident + reflected, (-incident + reflected) / impedance);
        const double through = transmission * incident;
        transmitted.sample(t, through, -through / impedance);
    }
    const auto sp = cfd::fdtd::s_parameters(reference, transmitted);

    cfd::fdtd::Maxwell3D solver({10, 9, 8, 1.0e-3, 1.0e-3, 1.0e-3, 0.7, 1.0, 1.0,
                                 cfd::fdtd::Boundary3D::pmc});
    solver.initialize_gaussian_ez(0.1, 0.15);
    solver.step(12U);
    const auto path = std::filesystem::temp_directory_path() / "cfd_solvers_fdtd_port_vtk_cli.vtk";
    cfd::fdtd::write_maxwell3d_vtk_ascii(solver, path);
    const auto bytes = std::filesystem::file_size(path);
    std::filesystem::remove(path);
    std::cout << "case=fdtd-port-vtk"
              << " S11_abs=" << std::abs(sp.s11)
              << " S21_abs=" << std::abs(sp.s21)
              << " vtk_bytes=" << bytes
              << " energy=" << solver.energy() << '\n';
    return (std::abs(std::abs(sp.s11) - reflection) < 3.0e-3
            && std::abs(std::abs(sp.s21) - transmission) < 3.0e-3 && bytes > 0U) ? 0 : 2;
}


int run_fdtd_cpml_tfsf() {
    using cfd::fdtd::Boundary1D;
    using cfd::fdtd::Maxwell1D;
    Maxwell1D absorber({400,1.0e-3,0.95,1.0,1.0,Boundary1D::cpml,24U,3.0,1.0e-8,5.0,0.0});
    absorber.initialize_gaussian(0.35,0.025);
    const double initial=absorber.energy();
    absorber.step(700U);
    const double residual=absorber.energy()/initial;

    Maxwell1D tfsf({500,1.0e-3,0.95,1.0,1.0,Boundary1D::cpml,32U,3.0,1.0e-10,6.0,0.0});
    const double dt=tfsf.dt(),t0=35.0*dt,tau=10.0*dt;
    tfsf.set_tfsf_source(120U,[=](double time){const double q=(time-t0)/tau;return std::exp(-q*q);});
    double scattered=0.0,total=0.0;
    for(std::size_t n=0;n<120U;++n){
        tfsf.step();
        scattered=std::max(scattered,std::abs(tfsf.electric()[80U]));
        total=std::max(total,std::abs(tfsf.electric()[160U]));
    }
    const double leakage=scattered/std::max(total,1.0e-300);
    std::cout << "case=fdtd-cpml-tfsf"
              << " cpml_energy_ratio=" << residual
              << " total_peak=" << total
              << " scattered_peak=" << scattered
              << " leakage_ratio=" << leakage << '\n';
    return residual<1.0e-4 && total>0.5 && leakage<1.0e-4 ? 0 : 2;
}

int run_fdtd_lumped_rlc() {
    using cfd::fdtd::Boundary1D;
    using cfd::fdtd::Maxwell1D;
    constexpr std::size_t cell=64U;
    constexpr double length=1.0e-3,area=1.0e-6,inductance=1.0e-6,field=0.2;
    Maxwell1D inductor({128,1.0e-3,0.5,1.0,1.0,Boundary1D::mur1});
    inductor.set_parallel_lumped_rlc(cell,length,area,std::numeric_limits<double>::infinity(),inductance,0.0);
    const double drive=length/(inductance*area);
    double expected=0.0;
    for(std::size_t n=0;n<20U;++n){
        inductor.set_hard_source(cell,field);
        inductor.step();
        expected+=inductor.dt()*drive*field;
    }

    Maxwell1D resistor({256,1.0e-3,0.8,1.0,1.0,Boundary1D::mur1});
    Maxwell1D reference({256,1.0e-3,0.8,1.0,1.0,Boundary1D::mur1});
    resistor.set_parallel_lumped_rlc(128U,1.0e-3,1.0e-6,25.0);
    resistor.initialize_gaussian(0.5,0.025);
    reference.initialize_gaussian(0.5,0.025);
    resistor.step(120U); reference.step(120U);
    const double dissipative_ratio=resistor.energy()/reference.energy();
    const double current_error=std::abs(inductor.lumped_inductor_current_density(cell)-expected);
    std::cout << "case=fdtd-lumped-rlc"
              << " inductor_current_density=" << inductor.lumped_inductor_current_density(cell)
              << " current_error=" << current_error
              << " resistor_energy_ratio=" << dissipative_ratio << '\n';
    return current_error<1.0e-12*std::max(1.0,std::abs(expected)) && dissipative_ratio<0.9 ? 0 : 2;
}

int run_multiphysics_coupling() {
    cfd::multiphysics::FieldRegistry registry;
    cfd::multiphysics::FieldMetadata temperature;
    temperature.name = "temperature";
    temperature.entities = 4U;
    temperature.units = cfd::multiphysics::kelvin_units();
    temperature.producer = "thermal";
    registry.add(temperature, 300.0);

    const cfd::multiphysics::CellGrid1D source{{0.0, 0.25, 0.5, 1.0}};
    const cfd::multiphysics::CellGrid1D target{{0.0, 0.5, 0.75, 1.0}};
    const std::vector<double> value{1.0, 2.0, 4.0};
    const auto mapped = cfd::multiphysics::conservative_cell_average_transfer(source, value, target);

    std::vector<double> state{1.0};
    const auto result = cfd::multiphysics::solve_partitioned_fixed_point(
        state,
        [](std::span<const double> x, std::span<double> y) { y[0] = std::cos(x[0]); });
    std::cout << "case=multiphysics-coupling"
              << " fields=" << registry.names().size()
              << " mapped0=" << mapped.front()
              << " fixed_point=" << state.front()
              << " iterations=" << result.iterations
              << " relaxation=" << result.relaxation
              << " converged=" << (result.converged ? 1 : 0) << '\n';
    return result.converged && std::abs(mapped.front() - 1.5) < 1.0e-14 ? 0 : 2;
}

int run_multiphysics_electrothermal() {
    using Type=cfd::fem::ScalarBoundaryType;
    cfd::fem::Heat2DConfig thermal;
    thermal.conductivity=1.0;
    thermal.volumetric_heat_capacity=2.0;
    thermal.dt=2.0e-3;
    cfd::multiphysics::JouleHeatingCoupler2D coupled(
        cfd::fem::make_rectangle_tri_mesh(24,12,2.0,1.0),thermal);
    coupled.set_electrical_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    coupled.set_electrical_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 4.0;},{}});
    coupled.initialize_temperature([](cfd::fem::Node2){return 300.0;});
    coupled.set_thermal_dirichlet([](cfd::fem::Node2,double){return 300.0;});
    coupled.solve_electrical([](cfd::fem::Node2){return 5.0;});
    coupled.thermal_run(10U);
    double mean_q=0.0,max_t=0.0;
    for(double q:coupled.joule_heating_density())mean_q+=q;
    mean_q/=static_cast<double>(coupled.joule_heating_density().size());
    for(double t:coupled.thermal().temperature())max_t=std::max(max_t,t);
    std::cout << "case=multiphysics-electrothermal"
              << " mean_joule_W_m3=" << mean_q
              << " max_temperature_K=" << max_t
              << " time_s=" << coupled.thermal().time()
              << " electric_converged=" << (coupled.electrical().linear_result().converged?1:0)
              << " thermal_converged=" << (coupled.thermal().linear_result().converged?1:0) << '\n';
    return coupled.electrical().linear_result().converged&&coupled.thermal().linear_result().converged
        &&std::abs(mean_q-20.0)<1.0e-6&&max_t>300.0?0:2;
}


int run_corrosion1d() {
    cfd::electrochemistry::CorrosionCell1DConfig config;
    config.metal_potential = 0.050;
    config.bulk_electrolyte_potential = 0.0;
    config.equilibrium_potential = 0.0;
    config.exchange_current_density = 1.0;
    config.electrolyte_conductivity = 5.0;
    config.electrolyte_length = 1.0e-3;
    config.electrons = 2.0;
    const auto result = cfd::electrochemistry::solve_corrosion_cell_1d(config);
    constexpr double seconds_per_year = 365.25 * 24.0 * 3600.0;
    std::cout << "case=electrochem-corrosion1d"
              << " current_A_per_m2=" << result.current_density
              << " surface_phi_V=" << result.surface_electrolyte_potential
              << " overpotential_V=" << result.overpotential
              << " dissolution_mol_per_m2_s=" << result.molar_dissolution_flux
              << " penetration_mm_per_year=" << result.penetration_rate * 1000.0 * seconds_per_year
              << " iterations=" << result.iterations
              << " converged=" << (result.converged ? 1 : 0) << '\n';
    return result.converged ? 0 : 2;
}

int run_nernst_planck1d() {
    constexpr double pi = 3.1415926535897932384626433832795;
    cfd::electrochemistry::NernstPlanck1D solver({256, 1.0, 298.15, 1.0e-4, 0.0,
                                                  cfd::electrochemistry::TransportBoundary1D::periodic});
    solver.set_potential([](double x) { return 2.0e-3 * std::sin(2.0 * pi * x); });
    const auto cation = solver.add_species("cation", 1, 1.0e-3, 1.0);
    const double amount0 = solver.total_amount(cation);
    solver.step(1000);
    std::cout << "case=electrochem-pnp1d"
              << " time=" << solver.time()
              << " amount_error=" << std::abs(solver.total_amount(cation) - amount0)
              << " min_concentration=" << solver.minimum_concentration(cation)
              << '\n';
    return 0;
}

int run_nernst_planck_poly() {
    constexpr double eps0 = 8.8541878128e-12;
    constexpr double eps_r = 78.5;
    constexpr double concentration = 1.0e-15;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(128, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(
        std::move(mesh),
        {298.15, 1.0e-6, eps_r, cfd::electrochemistry::ElectromigrationFluxScheme::scharfetterGummel});
    solver.add_species("space-charge", +1, 0.0, concentration);
    solver.set_potential_boundary("left", cfd::electrochemistry::PotentialBoundaryType::fixedPotential, 0.0);
    solver.set_potential_boundary("right", cfd::electrochemistry::PotentialBoundaryType::fixedPotential, 0.0);
    const auto result = solver.solve_poisson_potential(2000, 1.0e-12);
    const std::size_t mid = solver.potential().size() / 2U;
    const double x = solver.mesh().cells()[mid].center.x;
    const double rho = cfd::electrochemistry::faraday_constant * concentration;
    const double exact = rho * x * (1.0 - x) / (2.0 * eps0 * eps_r);
    std::cout << "case=electrochem-pnp-poly"
              << " potential_mid_V=" << solver.potential()[mid]
              << " exact_mid_V=" << exact
              << " relative_error=" << std::abs(solver.potential()[mid] - exact) / exact
              << " charge_C=" << solver.total_charge()
              << " iterations=" << result.iterations
              << " converged=" << (result.converged ? 1 : 0) << '\n';
    return result.converged ? 0 : 2;
}

int run_galvanic() {
    std::vector<cfd::electrochemistry::ElectrodeReaction> reactions{
        {"metal-A", 0.0, 2.0, 1.0, 1.0, 0.5, 0.5},
        {"metal-B", 0.2, 2.0, 1.0, 1.0, 0.5, 0.5},
    };
    const auto result = cfd::electrochemistry::solve_mixed_potential(reactions, 298.15);
    std::cout << "case=electrochem-galvanic"
              << " mixed_potential_V=" << result.potential
              << " net_current_A=" << result.total_current
              << " anodic_current_A=" << result.reaction_current[0]
              << " cathodic_current_A=" << result.reaction_current[1]
              << " iterations=" << result.iterations
              << " converged=" << (result.converged ? 1 : 0) << '\n';
    return result.converged ? 0 : 2;
}

int run_optics_lens() {
    cfd::optics::MaterialCatalog materials;
    const double n=materials.at("N-BK7").refractive_index_nm(587.5618);
    cfd::optics::SequentialOpticalSystem lens;
    lens.add_surface({cfd::optics::SurfaceType::sphere,0.0,50.0,12.0,n});
    lens.add_surface({cfd::optics::SurfaceType::sphere,5.0,-50.0,12.0,1.0});
    const double bfd=cfd::optics::paraxial_back_focal_distance(lens);
    std::vector<cfd::optics::Ray> rays;
    for(int i=-8;i<=8;++i) rays.push_back({{0.2*static_cast<double>(i),0.0,-10.0},{0.0,0.0,1.0},587.5618});
    const auto traced=lens.trace_many(rays);
    const double rms=cfd::optics::rms_spot_radius_at_plane(traced,5.0+bfd);
    const double qn=std::sqrt(n);
    const double ar=cfd::optics::thin_film_reflectance_normal(1.0,n,587.5618,{{qn,587.5618/(4.0*qn)}});
    std::cout << "case=optics-lens n_d=" << n << " back_focal_distance_mm=" << bfd
              << " rms_spot_mm=" << rms << " quarter_wave_R=" << ar << '\n';
    return (std::isfinite(bfd)&&std::isfinite(rms))?0:2;
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
    if (solver == "fvm-scalar-transport") return run_scalar_transport();
    if (solver == "fem-poisson1d") return run_fem();
    if (solver == "fem-poisson3d") return run_fem_poisson3d();
    if (solver == "fem-nonlinear-poisson") return run_fem_nonlinear_poisson();
    if (solver == "fem-darcy") return run_fem_darcy();
    if (solver == "fem-magnetostatic") return run_fem_magnetostatic();
    if (solver == "fem-modal-bar") return run_fem_modal_bar();
    if (solver == "fdtd-maxwell1d") return run_fdtd();
    if (solver == "fdtd-mur1d") return run_fdtd_mur1d();
    if (solver == "fdtd-dispersive1d") return run_fdtd_dispersive1d();
    if (solver == "fdtd-port-vtk") return run_fdtd_port_vtk();
    if (solver == "fdtd-cpml-tfsf") return run_fdtd_cpml_tfsf();
    if (solver == "fdtd-lumped-rlc") return run_fdtd_lumped_rlc();
    if (solver == "multiphysics-coupling") return run_multiphysics_coupling();
    if (solver == "multiphysics-electrothermal") return run_multiphysics_electrothermal();
    if (solver == "optics-snell") return run_optics();
    if (solver == "optics-lens") return run_optics_lens();
    if (solver == "electrochem-corrosion1d") return run_corrosion1d();
    if (solver == "electrochem-pnp1d") return run_nernst_planck1d();
    if (solver == "electrochem-pnp-poly") return run_nernst_planck_poly();
    if (solver == "electrochem-galvanic") return run_galvanic();
    const int continuation = run_continuation_case(solver);
    if (continuation >= 0) return continuation;
    std::cerr << "Unknown solver: " << solver << "\n";
    print_list();
    return 1;
}
