#include "cfd/fvm/resident_sycl.hpp"
#include "cfd/fvm/pressure_velocity.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition) throw std::runtime_error(std::string(message));
}

double csr_value(const cfd::core::CsrMatrix& a, std::size_t row, std::size_t col) {
    for (std::size_t k = a.row_offsets()[row]; k < a.row_offsets()[row + 1U]; ++k) {
        if (a.column_indices()[k] == col) return a.values()[k];
    }
    return 0.0;
}

void test_fixed_pressure_host_matrix() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(5U, 3U, 1U, 1.0, 0.5, 0.2);
    const auto a = cfd::fvm::build_orthogonal_pressure_matrix(mesh, 0.7, true, cfd::fvm::invalid_cell);
    require(a.rows() == mesh.cell_count(), "fixed-pressure host matrix dimensions");
    for (std::size_t r = 0U; r < a.rows(); ++r) {
        require(csr_value(a, r, r) > 0.0, "fixed-pressure host matrix positive diagonal");
        for (std::size_t k = a.row_offsets()[r]; k < a.row_offsets()[r + 1U]; ++k) {
            const auto c = a.column_indices()[k];
            require(std::abs(a.values()[k] - csr_value(a, c, r)) < 1.0e-13,
                    "fixed-pressure host matrix remains symmetric");
        }
    }
}

#if defined(CFD_HAS_SYCL)
void test_bicgstab_nonsymmetric() {
    cfd::core::CsrBuilder builder(3U, 3U);
    builder.add(0U, 0U, 4.0); builder.add(0U, 1U, -1.0);
    builder.add(1U, 0U, -2.0); builder.add(1U, 1U, 5.0); builder.add(1U, 2U, 1.0);
    builder.add(2U, 1U, -1.0); builder.add(2U, 2U, 3.0);
    cfd::core::SyclCsrLinearAlgebra solver(builder.build());
    const std::vector<double> exact{1.0, -2.0, 0.5};
    const std::vector<double> rhs{6.0, -11.5, 3.5};
    std::vector<double> x(3U, 0.0);
    const auto result = solver.bicgstab(rhs, x, 100U, 1.0e-12);
    require(result.converged, "resident BiCGStab converges on nonsymmetric system");
    for (std::size_t i = 0U; i < x.size(); ++i) {
        require(std::abs(x[i] - exact[i]) < 1.0e-9, "resident BiCGStab solution accuracy");
    }
}

void test_fixed_pressure_and_resident_scalar() {
    const auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(6U, 4U, 2U, 1.0, 0.7, 0.4, 0.11);
    cfd::fvm::ResidentIncompressibleConfig cfg{};
    cfg.kinematic_viscosity = 0.02;
    cfg.dt = 0.01;
    cfg.momentum_iterations = 120U;
    cfg.momentum_tolerance = 1.0e-11;
    cfg.pressure_iterations = 300U;
    cfg.pressure_tolerance = 1.0e-11;
    cfg.nonorthogonal_correctors = 2U;
    cfd::fvm::ResidentIncompressibleSycl flow(mesh, cfg, 0U);
    for (const auto& patch : mesh.patches()) {
        flow.set_velocity_boundary(patch.name, cfd::fvm::VelocityBoundaryType::fixedValue, {});
    }
    require(!mesh.patches().empty(), "test mesh has boundary patches");
    flow.set_pressure_boundary(mesh.patches().back().name, cfd::fvm::PressureBoundaryType::fixedValue, 0.0);
    flow.initialize_uniform({}, 0.0);
    flow.reset_transfer_stats();
    const auto info = flow.iterate_simple();
    require(info.pressure.converged, "fixed-pressure resident SIMPLE converges");
    require(info.continuity_l2 < 1.0e-12, "fixed-pressure quiescent continuity");
    require(flow.hot_loop_host_transfer_bytes() == 0U, "fixed-pressure hot loop has no host field transfer");

    cfd::fvm::ResidentScalarTransportConfig scfg{};
    scfg.dt = cfg.dt;
    scfg.diffusivity = 0.05;
    scfg.linear_iterations = 200U;
    scfg.linear_tolerance = 1.0e-12;
    cfd::fvm::ResidentScalarTransportSycl scalar(mesh, scfg, flow.mesh().queue());
    for (const auto& patch : mesh.patches()) {
        scalar.set_boundary(patch.name, cfd::fvm::PressureBoundaryType::fixedValue, 3.0);
    }
    scalar.initialize_uniform(3.0);
    scalar.reset_transfer_stats();
    const auto scalar_result = scalar.step(flow.face_flux_device());
    require(scalar_result.converged, "resident scalar BiCGStab converges");
    require(scalar.hot_loop_host_transfer_bytes() == 0U, "resident scalar hot loop has no host field transfer");
    require(scalar.steps() == 1U && std::abs(scalar.time() - scfg.dt) < 1.0e-14,
            "resident scalar advances time");

    std::vector<double> values(mesh.cell_count(), 0.0);
    scalar.download(values);
    for (double v : values) require(std::abs(v - 3.0) < 1.0e-9, "constant resident scalar is preserved");
}
#endif

} // namespace

int main() {
    try {
        test_fixed_pressure_host_matrix();
#if defined(CFD_HAS_SYCL)
        test_bicgstab_nonsymmetric();
        test_fixed_pressure_and_resident_scalar();
#endif
        std::cout << "v0.15.10 resident Krylov/scalar tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "v0.15.10 test failure: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
