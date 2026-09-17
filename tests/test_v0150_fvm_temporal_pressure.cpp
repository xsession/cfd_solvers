#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"
#include "cfd/fvm/poly_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> run_diffusion(cfd::fvm::TemporalScheme scheme, double dt, std::size_t steps) {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(48, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ScalarTransportConfig cfg;
    cfg.dt = dt;
    cfg.diffusivity = 0.05;
    cfg.linear_iterations = 300;
    cfg.gmres_restart = 30;
    cfg.linear_tolerance = 1.0e-12;
    cfg.temporal_scheme = scheme;
    cfg.crank_nicolson_off_centering = 1.0;
    cfd::fvm::ScalarTransport solver(std::move(mesh), cfg);
    solver.set_boundary("left", cfd::fvm::ScalarBoundaryType::fixedValue, 0.0);
    solver.set_boundary("right", cfd::fvm::ScalarBoundaryType::fixedValue, 0.0);
    solver.initialize([](cfd::fvm::Vec3 x) { return std::sin(3.14159265358979323846 * x.x); });
    solver.run(steps);
    return solver.values();
}

double rms_difference(const std::vector<double>& a, const std::vector<double>& b) {
    require(a.size() == b.size(), "RMS comparison size mismatch");
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = a[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum / static_cast<double>(a.size()));
}

void test_scalar_second_order_time_integration() {
    constexpr double final_time = 0.12;
    const auto reference = run_diffusion(cfd::fvm::TemporalScheme::crank_nicolson, 0.0005,
                                         static_cast<std::size_t>(final_time / 0.0005 + 0.5));
    const auto euler = run_diffusion(cfd::fvm::TemporalScheme::euler, 0.02, 6U);
    const auto bdf2 = run_diffusion(cfd::fvm::TemporalScheme::backward_bdf2, 0.02, 6U);
    const auto cn = run_diffusion(cfd::fvm::TemporalScheme::crank_nicolson, 0.02, 6U);
    const double euler_error = rms_difference(euler, reference);
    const double bdf2_error = rms_difference(bdf2, reference);
    const double cn_error = rms_difference(cn, reference);
    std::cout << "scalar temporal errors euler=" << euler_error
              << " bdf2=" << bdf2_error << " cn=" << cn_error << '\n';
    require(bdf2_error < euler_error * 0.35, "BDF2 should materially reduce temporal error");
    require(cn_error < euler_error * 0.12, "Crank-Nicolson should materially reduce temporal error");
}

cfd::fvm::CollocatedIterationInfo run_shear_step(cfd::fvm::TemporalScheme scheme, bool pimple) {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(4, 18, 1, 1.0, 1.0, 1.0);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.03;
    cfg.kinematic_viscosity = 0.04;
    cfg.include_convection = false;
    cfg.temporal_scheme = scheme;
    cfg.crank_nicolson_off_centering = 1.0;
    cfg.pressure_correctors = 2;
    cfg.outer_correctors = 2;
    cfg.nonorthogonal_correctors = 0;
    cfg.pressure_tolerance = 1.0e-11;
    cfd::fvm::CollocatedIncompressible solver(std::move(mesh), cfg);
    solver.set_velocity_boundary("bottom", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("top", cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("left", cfd::fvm::VelocityBoundaryType::zeroGradient);
    solver.set_velocity_boundary("right", cfd::fvm::VelocityBoundaryType::zeroGradient);
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    solver.initialize_fields(
        [](cfd::fvm::Vec3 x) { return cfd::fvm::Vec3{std::sin(3.14159265358979323846 * x.y), 0.0, 0.0}; },
        [](cfd::fvm::Vec3) { return 0.0; });
    auto info = pimple ? solver.step_pimple() : solver.step_piso();
    require(info.pressure.converged, "first second-order collocated pressure solve should converge");
    info = pimple ? solver.step_pimple() : solver.step_piso();
    require(info.pressure.converged, "history-enabled second-order collocated pressure solve should converge");
    require(solver.steps() == 2U && std::abs(solver.time() - 0.06) < 1.0e-14,
            "collocated physical time history should advance once per step");
    require(std::isfinite(solver.kinetic_energy()) && solver.kinetic_energy() > 0.0,
            "second-order collocated integration should retain finite kinetic energy");
    require(solver.continuity_l2() < 1.0e-8, "second-order collocated integration should preserve continuity");
    return info;
}

void test_collocated_second_order_time_integration() {
    (void)run_shear_step(cfd::fvm::TemporalScheme::backward_bdf2, false);
    (void)run_shear_step(cfd::fvm::TemporalScheme::crank_nicolson, true);
}

void test_pressure_amg_callback_path() {
    auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(8, 6, 1, 1.0, 1.0, 1.0, 0.2);
    cfd::fvm::CollocatedIncompressibleConfig cfg;
    cfg.dt = 0.02;
    cfg.kinematic_viscosity = 0.02;
    cfg.include_convection = false;
    cfg.pressure_correctors = 1;
    cfg.nonorthogonal_correctors = 1;
    cfg.pressure_tolerance = 1.0e-10;
    cfg.pressure_preconditioner = cfd::fvm::PressurePreconditionerKind::external_amg;
    cfd::fvm::CollocatedIncompressible solver(std::move(mesh), cfg);
    for (const auto name : {"left", "right", "bottom", "top"})
        solver.set_velocity_boundary(name, cfd::fvm::VelocityBoundaryType::fixedValue);
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_pressure_boundary("right", cfd::fvm::PressureBoundaryType::fixedValue, 0.0);
    std::size_t preconditioner_calls = 0U;
    solver.set_pressure_amg_cycle([&](std::span<const double> r, std::span<double> z) {
        ++preconditioner_calls;
        std::copy(r.begin(), r.end(), z.begin());
    });
    solver.initialize_fields(
        [](cfd::fvm::Vec3 x) { return cfd::fvm::Vec3{0.08 * (x.x - 0.5), -0.05 * (x.y - 0.5), 0.0}; },
        [](cfd::fvm::Vec3) { return 0.0; });
    const auto info = solver.step_piso();
    require(info.pressure.converged, "externally preconditioned pressure solve should converge");
    require(preconditioner_calls > 0U, "pressure PCG should invoke the configured AMG cycle");
}
} // namespace

int main() {
    try {
        test_scalar_second_order_time_integration();
        test_collocated_second_order_time_integration();
        test_pressure_amg_callback_path();
        std::cout << "v0.15.0 FVM temporal/pressure tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.0 FVM temporal/pressure test failed: " << error.what() << '\n';
        return 1;
    }
}
