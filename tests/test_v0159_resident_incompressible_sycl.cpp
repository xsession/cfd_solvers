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

void test_variable_mobility_pressure_matrix() {
    const auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(5U, 4U, 2U, 1.0, 0.8, 0.5, 0.21);
    std::vector<double> mobility(mesh.cell_count(), 0.0);
    for (std::size_t c = 0U; c < mobility.size(); ++c) {
        mobility[c] = 0.2 + 0.01 * static_cast<double>(c);
    }
    const auto a = cfd::fvm::build_orthogonal_pressure_matrix(mesh, mobility, false, 0U);
    require(a.rows() == mesh.cell_count(), "variable-mobility pressure matrix dimensions");
    require(std::abs(csr_value(a, 0U, 0U) - 1.0) < 1.0e-14,
            "variable-mobility pressure matrix pins reference cell");
    for (std::size_t r = 0U; r < a.rows(); ++r) {
        require(csr_value(a, r, r) > 0.0, "variable-mobility pressure diagonal positive");
        for (std::size_t k = a.row_offsets()[r]; k < a.row_offsets()[r + 1U]; ++k) {
            const std::size_t c = a.column_indices()[k];
            require(std::abs(a.values()[k] - csr_value(a, c, r)) < 1.0e-12,
                    "variable-mobility pressure matrix symmetric");
        }
    }

    bool found_nonorthogonal = false;
    for (std::size_t f = 0U; f < mesh.face_count(); ++f) {
        const auto d = cfd::fvm::decompose_face_area(mesh, f);
        const auto reconstructed = d.orthogonal_area + d.nonorthogonal_area;
        const auto sf = mesh.faces()[f].area;
        require(std::abs(reconstructed.x - sf.x) < 1.0e-13 &&
                std::abs(reconstructed.y - sf.y) < 1.0e-13 &&
                std::abs(reconstructed.z - sf.z) < 1.0e-13,
                "orthogonal plus nonorthogonal area reconstructs face area");
        found_nonorthogonal = found_nonorthogonal || cfd::fvm::magnitude(d.nonorthogonal_area) > 1.0e-12;
    }
    require(found_nonorthogonal, "sheared mesh exercises nonorthogonal correction geometry");
}

#if defined(CFD_HAS_SYCL)
void test_resident_simple_piso_pimple_quiescent() {
    const auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(6U, 5U, 2U, 1.0, 0.8, 0.4, 0.13);
    cfd::fvm::ResidentIncompressibleConfig cfg{};
    cfg.kinematic_viscosity = 0.03;
    cfg.dt = 0.01;
    cfg.momentum_sweeps = 12U;
    cfg.pressure_iterations = 300U;
    cfg.pressure_tolerance = 1.0e-11;
    cfg.nonorthogonal_correctors = 2U;
    cfg.pressure_correctors = 2U;
    cfg.outer_correctors = 2U;
    cfd::fvm::ResidentIncompressibleSycl solver(mesh, cfg, 0U);
    for (const auto& patch : mesh.patches()) {
        solver.set_velocity_boundary(patch.name, cfd::fvm::VelocityBoundaryType::fixedValue, {});
    }
    solver.initialize_uniform({}, 0.0);
    solver.reset_transfer_stats();

    const auto simple = solver.iterate_simple();
    require(simple.pressure.converged, "resident SIMPLE pressure solve converges");
    require(std::isfinite(simple.continuity_l2) && simple.continuity_l2 < 1.0e-12,
            "resident SIMPLE preserves quiescent continuity");

    const auto piso = solver.step_piso();
    require(piso.pressure.converged && piso.pressure_solves == cfg.pressure_correctors,
            "resident PISO executes requested pressure correctors");
    const auto pimple = solver.step_pimple();
    require(pimple.pressure.converged, "resident PIMPLE pressure solve converges");
    require(solver.steps() == 2U && std::abs(solver.time() - 2.0 * cfg.dt) < 1.0e-14,
            "resident PISO/PIMPLE advance physical time");
    require(solver.hot_loop_host_transfer_bytes() == 0U,
            "resident SIMPLE/PISO/PIMPLE hot loop performs no host field transfer");
    require(std::abs(solver.kinetic_energy()) < 1.0e-14,
            "resident quiescent kinetic energy remains zero");
}
#endif

} // namespace

int main() {
    try {
        test_variable_mobility_pressure_matrix();
#if defined(CFD_HAS_SYCL)
        test_resident_simple_piso_pimple_quiescent();
#endif
        std::cout << "v0.15.9 resident incompressible SYCL tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "v0.15.9 test failure: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
