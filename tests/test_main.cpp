#include "cfd/core/conjugate_gradient.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/electrochemistry/electrochemistry.hpp"
#include "cfd/solvers/electrochemistry/corrosion1d.hpp"
#include "cfd/solvers/electrochemistry/mixed_potential.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck_poly.hpp"
#include "cfd/core/decomposition.hpp"
#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fdtd/monitor.hpp"
#include "cfd/solvers/fdtd/maxwell3d.hpp"
#include "cfd/solvers/fdtd/port.hpp"
#include "cfd/solvers/fdtd/vtk_io.hpp"
#include "cfd/solvers/fem/poisson1d.hpp"
#include "cfd/fem/reference_element.hpp"
#include "cfd/fem/assembly.hpp"
#include "cfd/solvers/fem/scalar_diffusion2d.hpp"
#include "cfd/solvers/fem/electromagnetics2d.hpp"
#include "cfd/solvers/fem/poisson2d.hpp"
#include "cfd/solvers/fem/poisson3d.hpp"
#include "cfd/solvers/fem/heat2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include "cfd/solvers/fem/axisymmetric_elasticity.hpp"
#include "cfd/solvers/fem/nonlinear_poisson2d.hpp"
#include "cfd/fem/adaptivity.hpp"
#include "cfd/solvers/fem/darcy2d.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include "cfd/solvers/fem/modal_bar1d.hpp"
#include "cfd/solvers/fvm/diffusion2d.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/fvm/pressure_velocity.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/one_step_pull.hpp"
#include "cfd/solvers/lbm/precision_pull.hpp"
#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#endif
#include "cfd/solvers/optics/ray.hpp"
#include "cfd/solvers/optics/sequential.hpp"
#include "cfd/solvers/optics/materials.hpp"
#include "cfd/solvers/optics/paraxial.hpp"
#include "cfd/solvers/optics/polarization.hpp"
#include "cfd/multiphysics/field_registry.hpp"
#include "cfd/multiphysics/transfer.hpp"
#include "cfd/multiphysics/partitioned.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition) throw std::runtime_error(std::string(message));
}

void test_partition() {
    const auto p0 = cfd::core::partition_1d(10, 0, 3);
    const auto p1 = cfd::core::partition_1d(10, 1, 3);
    const auto p2 = cfd::core::partition_1d(10, 2, 3);
    require(p0.size() == 4 && p1.size() == 3 && p2.size() == 3, "balanced slab partition");
    require(p0.end == p1.begin && p1.end == p2.begin && p2.end == 10, "contiguous slab partition");
}


void test_cartesian_decomposition() {
    const cfd::distributed::Extent3 global{17, 15, 13};
    const auto grid = cfd::distributed::choose_process_grid(global, 8, 3);
    require(grid.x == 2 && grid.y == 2 && grid.z == 2, "cubic domain maps eight ranks to 2x2x2");

    std::size_t cells = 0;
    for (std::size_t rank = 0; rank < grid.size(); ++rank) {
        const auto brick = cfd::distributed::make_brick(global, grid, rank);
        require(brick.extent.x > 0 && brick.extent.y > 0 && brick.extent.z > 0,
                "Cartesian bricks are non-empty");
        cells += brick.extent.cells();
        const auto round_trip = cfd::distributed::rank_from_coords(brick.coordinate, grid);
        require(round_trip == rank, "Cartesian rank/coordinate round trip");
    }
    require(cells == global.cells(), "Cartesian bricks cover global domain exactly");

    const auto b0 = cfd::distributed::make_brick(global, grid, 0);
    require(cfd::distributed::neighbor_rank(b0, {-1, 0, 0}) == 1,
            "periodic x neighbor wraps on two-rank axis");
    require(cfd::distributed::neighbor_rank(b0, {-1, 0, 0}, {false, true, true}) == -1,
            "nonperiodic Cartesian boundary reports no neighbor");

    const auto a0 = cfd::distributed::assign_device(0, 6, 2);
    const auto a3 = cfd::distributed::assign_device(3, 6, 2);
    require(a0.device_ordinal == 0 && a3.device_ordinal == 1 && a3.oversubscribed,
            "local ranks map deterministically across visible devices");
}

template<class Descriptor>
void test_virtual_distributed_lbm(std::string_view label) {
    const cfd::distributed::Extent3 global{17, 15, 13};
    const auto grid = cfd::distributed::choose_process_grid(global, 8, 3);
    const cfd::lbm::DistributedLbmConfig config{global, 0.73F};
    std::vector<cfd::lbm::DistributedPullBlock<Descriptor>> blocks;
    blocks.reserve(grid.size());
    for (std::size_t rank = 0; rank < grid.size(); ++rank) {
        blocks.emplace_back(config, cfd::distributed::make_brick(global, grid, rank));
        blocks.back().initialize_taylor_green(0.017F);
    }

    cfd::lbm::PrecisionPullSolver<Descriptor, float> reference(
        {global.x, global.y, global.z, 0.73F});
    reference.initialize_taylor_green(0.017F);

    constexpr std::size_t steps = 9;
    for (std::size_t step = 0; step < steps; ++step) cfd::lbm::virtual_step(blocks);
    reference.step(steps);

    const auto distributed = cfd::lbm::gather_virtual_macroscopic(blocks);
    const auto serial = reference.compute_macroscopic();
    double max_error = 0.0;
    for (std::size_t i = 0; i < serial.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(distributed.rho[i] - serial.rho[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(distributed.ux[i] - serial.ux[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(distributed.uy[i] - serial.uy[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(distributed.uz[i] - serial.uz[i])));
    }
    require(max_error < 2.0e-6, label);

    double distributed_mass = 0.0;
    for (const auto& block : blocks) distributed_mass += block.local_mass();
    require(std::abs(distributed_mass - static_cast<double>(reference.mass())) < 5.0e-3,
            "distributed LBM preserves serial-reference mass");
}

void test_selective_halo_plan() {
    using cfd::distributed::NeighborOffset;
    const cfd::distributed::Extent3 e{20, 18, 16};

    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q19Descriptor>({1, 0, 0}).size() == 5,
            "D3Q19 face crossing population count");
    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q19Descriptor>({1, 1, 0}).size() == 1,
            "D3Q19 edge crossing population count");
    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q19Descriptor>({1, 1, 1}).empty(),
            "D3Q19 corner has no crossing population");
    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q27Descriptor>({1, 0, 0}).size() == 9,
            "D3Q27 face crossing population count");
    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q27Descriptor>({1, 1, 0}).size() == 3,
            "D3Q27 edge crossing population count");
    require(cfd::distributed::crossing_populations<cfd::lbm::D3Q27Descriptor>({1, 1, 1}).size() == 1,
            "D3Q27 corner crossing population count");

    const auto d19_full = cfd::distributed::full_halo_traffic<cfd::lbm::D3Q19Descriptor>(e);
    const auto d19_selective = cfd::distributed::selective_halo_traffic<cfd::lbm::D3Q19Descriptor>(e);
    const auto d27_full = cfd::distributed::full_halo_traffic<cfd::lbm::D3Q27Descriptor>(e);
    const auto d27_selective = cfd::distributed::selective_halo_traffic<cfd::lbm::D3Q27Descriptor>(e);
    require(d19_selective.bytes < d19_full.bytes / 3U, "D3Q19 selective halo cuts traffic by more than 3x");
    require(d27_selective.bytes < d27_full.bytes / 3U, "D3Q27 selective halo cuts traffic by more than 3x");
    require(d19_selective.messages == 18, "D3Q19 omits eight empty corner messages");
    require(d27_selective.messages == 26, "D3Q27 retains face/edge/corner messages");
}

template<class Descriptor>
void test_virtual_selective_distributed_lbm(std::string_view label) {
    const cfd::distributed::Extent3 global{17, 15, 13};
    const cfd::distributed::ProcessGrid grid{2, 2, 2};
    const cfd::lbm::DistributedLbmConfig config{global, 0.73F};
    std::vector<cfd::lbm::DistributedPullBlock<Descriptor>> full;
    std::vector<cfd::lbm::DistributedPullBlock<Descriptor>> selective;
    full.reserve(grid.size());
    selective.reserve(grid.size());
    for (std::size_t rank = 0; rank < grid.size(); ++rank) {
        const auto brick = cfd::distributed::make_brick(global, grid, rank);
        full.emplace_back(config, brick);
        selective.emplace_back(config, brick);
        full.back().initialize_taylor_green(0.017F);
        selective.back().initialize_taylor_green(0.017F);
    }
    constexpr std::size_t steps = 7;
    for (std::size_t i = 0; i < steps; ++i) {
        cfd::lbm::virtual_step(full);
        cfd::lbm::virtual_selective_step(selective);
    }
    const auto a = cfd::lbm::gather_virtual_macroscopic(full);
    const auto b = cfd::lbm::gather_virtual_macroscopic(selective);
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i] - b.rho[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i] - b.ux[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i] - b.uy[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uz[i] - b.uz[i])));
    }
    require(max_error == 0.0, label);
}

void test_virtual_selective_distributed_parity() {
    test_virtual_selective_distributed_lbm<cfd::lbm::D3Q19Descriptor>(
        "D3Q19 selective halo is bitwise-identical to full halo exchange");
    test_virtual_selective_distributed_lbm<cfd::lbm::D3Q27Descriptor>(
        "D3Q27 selective halo is bitwise-identical to full halo exchange");
}

void test_virtual_distributed_lbm_parity() {
    test_virtual_distributed_lbm<cfd::lbm::D3Q19Descriptor>("D3Q19 2x2x2 virtual-distributed parity");
    test_virtual_distributed_lbm<cfd::lbm::D3Q27Descriptor>("D3Q27 2x2x2 virtual-distributed parity");
}

void test_distributed_checkpoint_restart() {
    const cfd::distributed::Extent3 global{12, 10, 8};
    const cfd::distributed::ProcessGrid grid{1, 1, 1};
    const cfd::lbm::DistributedLbmConfig config{global, 0.71F};
    std::vector<cfd::lbm::DistributedPullBlock<cfd::lbm::D3Q19Descriptor>> blocks;
    blocks.emplace_back(config, cfd::distributed::make_brick(global, grid, 0));
    blocks[0].initialize_taylor_green(0.02F);
    for (int i = 0; i < 5; ++i) cfd::lbm::virtual_step(blocks);

    const auto path = std::filesystem::temp_directory_path() / "cfd_solvers_phase3_checkpoint.bin";
    blocks[0].save_checkpoint(path.string());

    std::vector<cfd::lbm::DistributedPullBlock<cfd::lbm::D3Q19Descriptor>> restored;
    restored.emplace_back(config, cfd::distributed::make_brick(global, grid, 0));
    restored[0].load_checkpoint(path.string());
    require(restored[0].time_step() == blocks[0].time_step(), "checkpoint restores time step");

    const auto before_a = blocks[0].compute_macroscopic();
    const auto before_b = restored[0].compute_macroscopic();
    double before_error = 0.0;
    for (std::size_t i = 0; i < before_a.rho.size(); ++i) {
        before_error = std::max(before_error, std::abs(static_cast<double>(before_a.rho[i] - before_b.rho[i])));
        before_error = std::max(before_error, std::abs(static_cast<double>(before_a.ux[i] - before_b.ux[i])));
        before_error = std::max(before_error, std::abs(static_cast<double>(before_a.uy[i] - before_b.uy[i])));
    }
    require(before_error == 0.0, "checkpoint reproduces populations exactly");

    cfd::lbm::virtual_step(blocks);
    cfd::lbm::virtual_step(restored);
    const auto after_a = blocks[0].compute_macroscopic();
    const auto after_b = restored[0].compute_macroscopic();
    double after_error = 0.0;
    for (std::size_t i = 0; i < after_a.rho.size(); ++i) {
        after_error = std::max(after_error, std::abs(static_cast<double>(after_a.rho[i] - after_b.rho[i])));
        after_error = std::max(after_error, std::abs(static_cast<double>(after_a.ux[i] - after_b.ux[i])));
        after_error = std::max(after_error, std::abs(static_cast<double>(after_a.uy[i] - after_b.uy[i])));
        after_error = std::max(after_error, std::abs(static_cast<double>(after_a.uz[i] - after_b.uz[i])));
    }
    require(after_error == 0.0, "checkpoint restart remains bitwise identical after next step");
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

void test_lbm_uniform() {
    cfd::lbm::D2Q9Solver solver({48, 32, 0.7F});
    solver.initialize_uniform(1.0F, 0.02F, 0.0F);
    const double mass0 = solver.mass();
    solver.step(30);
    require(std::abs(solver.mass() - mass0) < 2.0e-3, "LBM mass conservation");
    require(std::abs(solver.max_speed() - 0.02F) < 5.0e-4F, "LBM uniform flow invariance");
}


template<class Descriptor>
void test_descriptor_moments(std::string_view name) {
    double sum_w = 0.0;
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    for (int q = 0; q < Descriptor::q; ++q) {
        const double w = Descriptor::weight(q);
        const double cx = Descriptor::cx(q);
        const double cy = Descriptor::cy(q);
        const double cz = Descriptor::cz(q);
        sum_w += w;
        xx += w * cx * cx;
        yy += w * cy * cy;
        zz += w * cz * cz;
        xy += w * cx * cy;
        require(Descriptor::opposite(Descriptor::opposite(q)) == q, "LBM descriptor opposite symmetry");
    }
    require(std::abs(sum_w - 1.0) < 2.0e-7, name);
    require(std::abs(xx - 1.0 / 3.0) < 2.0e-7, "LBM descriptor xx isotropy");
    require(std::abs(yy - 1.0 / 3.0) < 2.0e-7, "LBM descriptor yy isotropy");
    if constexpr (Descriptor::dimensions == 3) {
        require(std::abs(zz - 1.0 / 3.0) < 2.0e-7, "LBM descriptor zz isotropy");
    }
    require(std::abs(xy) < 2.0e-7, "LBM descriptor xy isotropy");
}

void test_lbm_descriptors() {
    test_descriptor_moments<cfd::lbm::D2Q9InPlaceDescriptor>("D2Q9 weights");
    test_descriptor_moments<cfd::lbm::D3Q19Descriptor>("D3Q19 weights");
    test_descriptor_moments<cfd::lbm::D3Q27Descriptor>("D3Q27 weights");
}

void test_inplace_d2q9_against_pingpong() {
    constexpr std::size_t nx = 40;
    constexpr std::size_t ny = 32;
    cfd::lbm::D2Q9Solver reference({nx, ny, 0.7F});
    cfd::lbm::D2Q9InPlaceSolver inplace({nx, ny, 1, 0.7F});
    reference.initialize_taylor_green(0.02F);
    inplace.initialize_taylor_green(0.02F);
    reference.step(24);
    inplace.step(24);

    const auto macro = inplace.compute_macroscopic();
    double max_rho_error = 0.0;
    double max_u_error = 0.0;
    for (std::size_t i = 0; i < nx * ny; ++i) {
        max_rho_error = std::max(max_rho_error,
                                 std::abs(static_cast<double>(macro.rho[i] - reference.density()[i])));
        max_u_error = std::max(max_u_error,
                               std::abs(static_cast<double>(macro.ux[i] - reference.velocity_x()[i])));
        max_u_error = std::max(max_u_error,
                               std::abs(static_cast<double>(macro.uy[i] - reference.velocity_y()[i])));
    }
    require(max_rho_error < 2.0e-6, "Esoteric-Pull D2Q9 density parity with two-grid pull");
    require(max_u_error < 2.0e-6, "Esoteric-Pull D2Q9 velocity parity with two-grid pull");
    require(inplace.population_bytes() * 2U == inplace.two_lattice_population_bytes(),
            "single-grid population memory accounting");
}


template<class Descriptor>
void test_inplace_against_generic_pull(std::string_view label) {
    const std::size_t nz = Descriptor::dimensions == 3 ? 10U : 1U;
    const cfd::lbm::InPlaceLbmConfig config{18, 14, nz, 0.73F};
    cfd::lbm::EsotericPullSolver<Descriptor> inplace(config);
    cfd::lbm::OneStepPullSolver<Descriptor> reference(config);
    inplace.initialize_taylor_green(0.017F);
    reference.initialize_taylor_green(0.017F);

    // Odd and even counts exercise both physical layouts of the in-place scheme.
    for (const std::size_t steps : {7U, 8U}) {
        inplace.initialize_taylor_green(0.017F);
        reference.initialize_taylor_green(0.017F);
        inplace.step(steps);
        reference.step(steps);
        const auto a = inplace.compute_macroscopic();
        const auto b = reference.compute_macroscopic();
        double max_error = 0.0;
        for (std::size_t i = 0; i < a.rho.size(); ++i) {
            max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i] - b.rho[i])));
            max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i] - b.ux[i])));
            max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i] - b.uy[i])));
            max_error = std::max(max_error, std::abs(static_cast<double>(a.uz[i] - b.uz[i])));
        }
        require(max_error < 5.0e-6, label);
    }
}

void test_generic_pull_parity() {
    test_inplace_against_generic_pull<cfd::lbm::D2Q9InPlaceDescriptor>("D2Q9 in-place/two-grid parity");
    test_inplace_against_generic_pull<cfd::lbm::D3Q19Descriptor>("D3Q19 in-place/two-grid parity");
    test_inplace_against_generic_pull<cfd::lbm::D3Q27Descriptor>("D3Q27 in-place/two-grid parity");
}

void test_d3q19_uniform() {
    cfd::lbm::D3Q19Solver solver({24, 16, 12, 0.72F});
    solver.initialize_uniform(1.0F, 0.015F, -0.01F, 0.005F);
    const double mass0 = solver.mass();
    solver.step(20);
    const double mass1 = solver.mass();
    const auto macro = solver.compute_macroscopic();
    require(std::abs(mass1 - mass0) < 1.0e-2, "D3Q19 mass conservation");
    for (std::size_t i = 0; i < macro.rho.size(); ++i) {
        require(std::abs(macro.rho[i] - 1.0F) < 2.0e-5F, "D3Q19 uniform density invariance");
        require(std::abs(macro.ux[i] - 0.015F) < 2.0e-5F, "D3Q19 uniform ux invariance");
        require(std::abs(macro.uy[i] + 0.01F) < 2.0e-5F, "D3Q19 uniform uy invariance");
        require(std::abs(macro.uz[i] - 0.005F) < 2.0e-5F, "D3Q19 uniform uz invariance");
    }
}

void test_d3q27_taylor_green_decay() {
    cfd::lbm::D3Q27Solver solver({20, 18, 16, 0.8F});
    solver.initialize_taylor_green(0.025F);
    const double mass0 = solver.mass();
    const double energy0 = solver.kinetic_energy();
    solver.step(30);
    const double mass1 = solver.mass();
    const double energy1 = solver.kinetic_energy();
    require(std::abs(mass1 - mass0) < 2.0e-2, "D3Q27 mass conservation");
    require(std::isfinite(energy1) && energy1 > 0.0 && energy1 < energy0,
            "D3Q27 viscous Taylor-Green energy decay");
}

template<class Descriptor>
void test_forced_inplace_parity_descriptor(std::string_view label) {
    const std::size_t nz = Descriptor::dimensions == 3 ? 10U : 1U;
    const float az = Descriptor::dimensions == 3 ? 0.5e-6F : 0.0F;
    const cfd::lbm::InPlaceLbmConfig config{18, 14, nz, 0.73F, 2.0e-6F, -1.0e-6F, az};
    cfd::lbm::EsotericPullSolver<Descriptor> inplace(config);
    cfd::lbm::OneStepPullSolver<Descriptor> reference(config);
    inplace.initialize_uniform();
    reference.initialize_uniform();
    inplace.step(60);
    reference.step(60);
    const auto a = inplace.compute_macroscopic();
    const auto b = reference.compute_macroscopic();
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i] - b.rho[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i] - b.ux[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i] - b.uy[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uz[i] - b.uz[i])));
    }
    require(max_error < 7.0e-6, label);
    require(a.ux[a.ux.size() / 2] > 0.0F, "positive body acceleration drives positive velocity");
}

void test_forced_inplace_parity() {
    test_forced_inplace_parity_descriptor<cfd::lbm::D2Q9InPlaceDescriptor>(
        "D2Q9 Guo forcing in-place/two-grid parity");
    test_forced_inplace_parity_descriptor<cfd::lbm::D3Q19Descriptor>(
        "D3Q19 Guo forcing in-place/two-grid parity");
    test_forced_inplace_parity_descriptor<cfd::lbm::D3Q27Descriptor>(
        "D3Q27 Guo forcing in-place/two-grid parity");
}

void test_stationary_wall_resets_velocity() {
    cfd::lbm::D2Q9Solver solver({12, 10, 0.8F});
    solver.set_wall_velocity(4, 5, 0.05F, -0.02F);
    solver.set_solid(4, 5);
    solver.initialize_uniform();
    solver.compute_macroscopic();
    const std::size_t i = 5U * 12U + 4U;
    require(std::abs(solver.velocity_x()[i]) < 1.0e-12F,
            "set_solid resets previous moving-wall x velocity");
    require(std::abs(solver.velocity_y()[i]) < 1.0e-12F,
            "set_solid resets previous moving-wall y velocity");
}

void test_poiseuille_channel() {
    constexpr std::size_t nx = 24;
    constexpr std::size_t ny = 14;
    constexpr float tau = 0.8F;
    constexpr float acceleration = 1.0e-5F;
    cfd::lbm::D2Q9Solver solver({nx, ny, tau, acceleration, 0.0F});
    for (std::size_t x = 0; x < nx; ++x) {
        solver.set_solid(x, 0);
        solver.set_solid(x, ny - 1);
    }
    solver.initialize_uniform();
    solver.step(5000);

    const double nu = (static_cast<double>(tau) - 0.5) / 3.0;
    const double height = static_cast<double>(ny - 2);
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t y = 1; y + 1 < ny; ++y) {
        double numerical = 0.0;
        for (std::size_t x = 0; x < nx; ++x) numerical += solver.velocity_x()[y * nx + x];
        numerical /= static_cast<double>(nx);
        const double wall_y = static_cast<double>(y) - 0.5;
        const double exact = static_cast<double>(acceleration) * wall_y * (height - wall_y) / (2.0 * nu);
        const double error = numerical - exact;
        error2 += error * error;
        exact2 += exact * exact;
    }
    require(std::sqrt(error2 / exact2) < 0.01, "Poiseuille profile matches analytical parabola");
}

void test_lid_driven_cavity() {
    constexpr std::size_t n = 48;
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
    solver.step(6000);
    const auto idx = [](std::size_t x, std::size_t y) { return y * n + x; };
    const float center_ux = solver.velocity_x()[idx(n / 2, n / 2)];
    const float near_lid_ux = solver.velocity_x()[idx(n / 2, n - 2)];
    require(near_lid_ux > 0.03F, "moving lid transfers positive momentum to adjacent fluid");
    require(center_ux < -0.003F, "lid-driven cavity develops clockwise primary vortex");
    require(std::isfinite(solver.max_speed()) && solver.max_speed() < 0.08F,
            "lid-driven cavity remains stable and subsonic");
}

void test_velocity_inlet_pressure_outlet() {
    constexpr std::size_t nx = 80;
    constexpr std::size_t ny = 24;
    cfd::lbm::D2Q9Solver solver({nx, ny, 0.8F});
    for (std::size_t x = 0; x < nx; ++x) {
        solver.set_solid(x, 0);
        solver.set_solid(x, ny - 1);
    }
    solver.set_velocity_inlet_left(0.02F);
    solver.set_pressure_outlet_right(1.0F);
    solver.initialize_uniform(1.0F, 0.02F, 0.0F);
    solver.step(3000);

    const double inlet = solver.average_ux(0, 1, 1, ny - 1);
    const double mid = solver.average_ux(nx / 2, nx / 2 + 1, 1, ny - 1);
    const double outlet = solver.average_ux(nx - 1, nx, 1, ny - 1);
    require(std::abs(inlet - 0.02) < 1.0e-5, "Zou-He velocity inlet enforces target velocity");
    require(std::abs(mid - inlet) < 0.004, "inlet flow propagates through channel");
    require(std::abs(outlet - inlet) < 0.004, "pressure outlet preserves channel throughput");
    for (std::size_t y = 1; y + 1 < ny; ++y) {
        require(std::abs(solver.density()[y * nx + nx - 1] - 1.0F) < 2.0e-5F,
                "Zou-He pressure outlet enforces target density");
    }
}

template<class Real>
double taylor_green_relative_error(std::size_t n, std::size_t steps) {
    constexpr double tau = 0.8;
    constexpr double amplitude = 0.001;
    cfd::lbm::PrecisionPullSolver<cfd::lbm::D2Q9InPlaceDescriptor, Real> solver(
        {n, n, 1, static_cast<Real>(tau)});
    solver.initialize_taylor_green(static_cast<Real>(amplitude));
    solver.step(steps);
    const auto macro = solver.compute_macroscopic();
    const double nu = (tau - 0.5) / 3.0;
    const double k = 2.0 * std::numbers::pi / static_cast<double>(n);
    const double decay = std::exp(-2.0 * nu * k * k * static_cast<double>(steps));
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t y = 0; y < n; ++y) {
        for (std::size_t x = 0; x < n; ++x) {
            const double xf = (static_cast<double>(x) + 0.5) / static_cast<double>(n);
            const double yf = (static_cast<double>(y) + 0.5) / static_cast<double>(n);
            const double exact_ux = amplitude * std::sin(2.0 * std::numbers::pi * xf) *
                                    std::cos(2.0 * std::numbers::pi * yf) * decay;
            const double exact_uy = -amplitude * std::cos(2.0 * std::numbers::pi * xf) *
                                     std::sin(2.0 * std::numbers::pi * yf) * decay;
            const std::size_t i = y * n + x;
            const double ex = static_cast<double>(macro.ux[i]) - exact_ux;
            const double ey = static_cast<double>(macro.uy[i]) - exact_uy;
            error2 += ex * ex + ey * ey;
            exact2 += exact_ux * exact_ux + exact_uy * exact_uy;
        }
    }
    return std::sqrt(error2 / exact2);
}

void test_taylor_green_convergence() {
    const double e16 = taylor_green_relative_error<double>(16, 20);
    const double e32 = taylor_green_relative_error<double>(32, 20);
    const double e64 = taylor_green_relative_error<double>(64, 20);
    require(e32 < 0.35 * e16, "Taylor-Green error decreases from 16 to 32 cells");
    require(e64 < 0.35 * e32, "Taylor-Green error decreases from 32 to 64 cells");
    require(e64 < 0.002, "Taylor-Green fine-grid error bound");
}

void test_fp32_fp64_parity() {
    constexpr std::size_t n = 40;
    cfd::lbm::PrecisionPullSolver<cfd::lbm::D2Q9InPlaceDescriptor, float> fp32({n, n, 1, 0.73F});
    cfd::lbm::PrecisionPullSolver<cfd::lbm::D2Q9InPlaceDescriptor, double> fp64({n, n, 1, 0.73});
    fp32.initialize_taylor_green(0.02F);
    fp64.initialize_taylor_green(0.02);
    fp32.step(120);
    fp64.step(120);
    const auto a = fp32.compute_macroscopic();
    const auto b = fp64.compute_macroscopic();
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i]) - b.rho[i]));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i]) - b.ux[i]));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i]) - b.uy[i]));
    }
    require(max_error < 1.0e-6, "FP32 agrees with FP64 on low-Mach Taylor-Green regression");
}


#if defined(CFD_HAS_SYCL)
template<class Descriptor>
void test_cpu_sycl_parity(std::string_view label) {
    const cfd::lbm::InPlaceLbmConfig config{
        16, 14, Descriptor::dimensions == 3 ? 12U : 1U, 0.72F, 1.0e-6F, -0.5e-6F, 0.25e-6F};
    cfd::lbm::EsotericPullSolver<Descriptor> cpu(config);
    cfd::lbm::EsotericPullSyclSolver<Descriptor> gpu(config);
    cpu.initialize_taylor_green(0.018F);
    gpu.initialize_taylor_green(0.018F);
    cpu.step(12);
    gpu.step(12);
    gpu.wait();
    const auto a = cpu.compute_macroscopic();
    const auto b = gpu.download_macroscopic();
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i] - b.rho[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i] - b.ux[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i] - b.uy[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uz[i] - b.uz[i])));
    }
    require(max_error < 2.0e-5, label);
}

void test_sycl_parity() {
    test_cpu_sycl_parity<cfd::lbm::D2Q9InPlaceDescriptor>("D2Q9 CPU/SYCL parity");
    test_cpu_sycl_parity<cfd::lbm::D3Q19Descriptor>("D3Q19 CPU/SYCL parity");
    test_cpu_sycl_parity<cfd::lbm::D3Q27Descriptor>("D3Q27 CPU/SYCL parity");
}
#endif



void test_core_conjugate_gradient() {
    std::vector<double> rhs{2.0, 4.0, 6.0, 8.0};
    std::vector<double> x(rhs.size(), 0.0);
    cfd::core::ConjugateGradientWorkspace workspace;
    const auto result = cfd::core::conjugate_gradient(
        std::span<const double>(rhs.data(), rhs.size()),
        std::span<double>(x.data(), x.size()),
        [](std::span<const double> in, std::span<double> out) {
            for (std::size_t i = 0; i < in.size(); ++i) out[i] = 2.0 * in[i];
        },
        workspace, 16, 1.0e-13);
    require(result.converged, "shared matrix-free CG converges on diagonal SPD operator");
    for (std::size_t i = 0; i < x.size(); ++i) {
        require(std::abs(x[i] - 0.5 * rhs[i]) < 1.0e-12, "shared matrix-free CG solution");
    }
}

void test_fvm_poly_mesh_operators() {
    constexpr std::size_t nx = 8;
    constexpr std::size_t ny = 7;
    constexpr std::size_t nz = 6;
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(nx, ny, nz, 1.2, 0.9, 1.1);
    const std::size_t expected_internal = (nx - 1U) * ny * nz + nx * (ny - 1U) * nz + nx * ny * (nz - 1U);
    const std::size_t expected_boundary = 2U * (ny * nz + nx * nz + nx * ny);
    require(mesh.cell_count() == nx * ny * nz, "poly mesh cell count");
    require(mesh.face_count() == expected_internal + expected_boundary, "poly mesh face count");
    require(mesh.boundary_face_count() == expected_boundary, "poly mesh boundary face count");

    std::vector<double> linear(mesh.cell_count());
    std::vector<cfd::fvm::Vec3> vector_field(mesh.cell_count());
    std::vector<double> quadratic(mesh.cell_count());
    for (std::size_t c = 0; c < mesh.cell_count(); ++c) {
        const auto x = mesh.cells()[c].center;
        linear[c] = x.x + 2.0 * x.y + 3.0 * x.z;
        vector_field[c] = {x.x, 2.0 * x.y, 3.0 * x.z};
        quadratic[c] = x.x * x.x + x.y * x.y + x.z * x.z;
    }
    std::vector<double> linear_boundary(mesh.face_count(), 0.0);
    std::vector<cfd::fvm::Vec3> vector_boundary(mesh.face_count());
    std::vector<double> quadratic_boundary(mesh.face_count(), 0.0);
    for (std::size_t f = 0; f < mesh.face_count(); ++f) {
        if (!mesh.faces()[f].boundary()) continue;
        const auto x = mesh.faces()[f].center;
        linear_boundary[f] = x.x + 2.0 * x.y + 3.0 * x.z;
        vector_boundary[f] = {x.x, 2.0 * x.y, 3.0 * x.z};
        quadratic_boundary[f] = x.x * x.x + x.y * x.y + x.z * x.z;
    }

    const auto grad = cfd::fvm::gauss_gradient_scalar(mesh, linear, linear_boundary);
    double grad_error = 0.0;
    for (const auto g : grad) {
        grad_error = std::max(grad_error, std::abs(g.x - 1.0));
        grad_error = std::max(grad_error, std::abs(g.y - 2.0));
        grad_error = std::max(grad_error, std::abs(g.z - 3.0));
    }
    require(grad_error < 2.0e-12, "Gauss gradient exact for linear field on Cartesian hexa mesh");

    const auto div = cfd::fvm::gauss_divergence_vector(mesh, vector_field, vector_boundary);
    double div_error = 0.0;
    for (const double value : div) div_error = std::max(div_error, std::abs(value - 6.0));
    require(div_error < 2.0e-12, "Gauss divergence exact for linear vector field");

    const auto lap = cfd::fvm::orthogonal_laplacian_scalar(mesh, quadratic, 1.0, quadratic_boundary);
    auto idx = [](std::size_t i, std::size_t j, std::size_t k) { return (k * ny + j) * nx + i; };
    double lap_error = 0.0;
    for (std::size_t k = 1; k + 1 < nz; ++k)
        for (std::size_t j = 1; j + 1 < ny; ++j)
            for (std::size_t i = 1; i + 1 < nx; ++i)
                lap_error = std::max(lap_error, std::abs(lap[idx(i, j, k)] - 6.0));
    require(lap_error < 2.0e-11, "orthogonal FVM Laplacian exact for interior quadratic field");
}

void test_fvm_face_schemes() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(10, 8, 6, 1.0, 0.8, 0.6);
    std::vector<cfd::fvm::Vec3> velocity(mesh.cell_count(), {1.0, 0.0, 0.0});
    const auto flux = cfd::fvm::face_flux_from_velocity(mesh, velocity);

    std::vector<double> constant(mesh.cell_count(), 3.5);
    const auto constant_div = cfd::fvm::convective_divergence_scalar(
        mesh, constant, flux, cfd::fvm::FaceInterpolationScheme::upwind);
    double constant_error = 0.0;
    for (double value : constant_div) constant_error = std::max(constant_error, std::abs(value));
    require(constant_error < 1.0e-12, "upwind convection preserves constant scalar field");

    std::vector<double> linear(mesh.cell_count());
    std::vector<double> boundary(mesh.face_count(), 0.0);
    for (std::size_t c=0;c<mesh.cell_count();++c) linear[c]=mesh.cells()[c].center.x;
    for (std::size_t f=0;f<mesh.face_count();++f) {
        if (mesh.faces()[f].boundary()) boundary[f]=mesh.faces()[f].center.x;
    }
    const auto linear_div = cfd::fvm::convective_divergence_scalar(
        mesh, linear, flux, cfd::fvm::FaceInterpolationScheme::linear, boundary);
    double linear_error = 0.0;
    for (double value : linear_div) linear_error = std::max(linear_error, std::abs(value-1.0));
    require(linear_error < 2.0e-12, "linear face interpolation gives exact uniform-advection divergence");

    const auto face_upwind = cfd::fvm::interpolate_scalar_to_faces(
        mesh, linear, cfd::fvm::FaceInterpolationScheme::upwind, flux, boundary);
    bool checked_positive_internal = false;
    for (std::size_t f=0;f<mesh.face_count();++f) {
        const auto& face=mesh.faces()[f];
        if (!face.boundary() && face.area.x>0.0 && std::abs(face.area.y)<1.0e-15 && std::abs(face.area.z)<1.0e-15) {
            require(std::abs(face_upwind[f]-linear[face.owner])<1.0e-14,
                    "positive face flux selects owner value in upwind interpolation");
            checked_positive_internal=true;
            break;
        }
    }
    require(checked_positive_internal, "found positive internal x face for upwind test");

    const auto vector_div = cfd::fvm::convective_divergence_vector(
        mesh, velocity, flux, cfd::fvm::FaceInterpolationScheme::upwind);
    double vector_error=0.0;
    for (auto value : vector_div) vector_error=std::max(vector_error,cfd::fvm::magnitude(value));
    require(vector_error<1.0e-12,"uniform vector momentum has zero conservative convective divergence");
}


void test_fvm_rhie_chow_checkerboard() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 4, 1, 1.0, 1.0, 1.0);
    std::vector<cfd::fvm::Vec3> h_by_a(mesh.cell_count());
    std::vector<double> mobility(mesh.cell_count(), 0.01);
    std::vector<double> pressure(mesh.cell_count());
    for (std::size_t c = 0; c < pressure.size(); ++c) {
        pressure[c] = ((c % 8U) & 1U) ? 1.0 : -1.0;
    }
    std::vector<cfd::fvm::VelocityBoundaryCondition> velocity_bc(mesh.patches().size());
    std::vector<cfd::fvm::PressureBoundaryCondition> pressure_bc(mesh.patches().size());
    const auto flux = cfd::fvm::rhie_chow_face_flux(
        mesh, h_by_a, mobility, pressure, velocity_bc, pressure_bc, false);
    double internal_response = 0.0;
    for (std::size_t f = 0; f < mesh.face_count(); ++f) {
        if (!mesh.faces()[f].boundary()) internal_response += std::abs(flux[f]);
    }
    require(internal_response > 0.2,
            "Rhie-Chow direct face pressure difference responds to checkerboard pressure mode");
}

void test_fvm_collocated_nonorthogonal_coupling() {
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
        solver.set_velocity_boundary(name, cfd::fvm::VelocityBoundaryType::fixedValue, {0.0, 0.0, 0.0});
    }
    solver.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    solver.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    auto velocity = [](cfd::fvm::Vec3 x) {
        return cfd::fvm::Vec3{0.1 * (x.x - 0.5), -0.08 * (x.y - 0.5), 0.0};
    };
    solver.initialize_fields(velocity, [](cfd::fvm::Vec3) { return 0.0; });
    const double before = solver.continuity_l2();
    const auto piso = solver.step_piso();
    require(piso.pressure.converged, "collocated PISO pressure solve converges on sheared mesh");
    require(piso.continuity_l2 < before * 1.0e-3,
            "non-orthogonal PISO correction strongly reduces flux divergence");

    solver.initialize_fields(velocity, [](cfd::fvm::Vec3) { return 0.0; });
    const auto pimple = solver.step_pimple();
    require(pimple.pressure.converged, "collocated PIMPLE pressure solve converges on sheared mesh");
    require(pimple.continuity_l2 < piso.continuity_l2 * 0.35,
            "PIMPLE outer corrector further reduces sheared-mesh continuity error");
}

void test_fvm_collocated_pressure_channel() {
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
    const auto info = solver.solve_simple(300, 2.0e-7, 1.0e-6);
    require(info.pressure.converged, "SIMPLE pressure solve converges for pressure-driven channel");
    for (const auto& result : solver.momentum_results()) {
        require(result.converged, "ILU(0)-GMRES momentum solve converges for pressure-driven channel");
    }
    require(solver.continuity_l2() < 1.0e-6, "SIMPLE pressure-driven channel satisfies continuity");

    double error2 = 0.0;
    double exact2 = 0.0;
    std::size_t count = 0;
    const double gradient = pressure_drop / length;
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const auto x = solver.mesh().cells()[c].center;
        if (x.x < 0.8 || x.x > 1.2) continue;
        const double exact = gradient * x.y * (height - x.y) / (2.0 * nu);
        const double error = solver.velocity()[c].x - exact;
        error2 += error * error;
        exact2 += exact * exact;
        ++count;
    }
    require(count > 0, "channel validation samples center section");
    require(std::sqrt(error2 / exact2) < 0.02,
            "collocated SIMPLE channel matches analytical Poiseuille profile");
}

void test_fvm_collocated_lid_cavity() {
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
    const auto info = solver.solve_simple(300, 2.0e-6, 2.0e-6);
    require(info.pressure.converged, "SIMPLE pressure solve converges for lid-driven cavity");
    require(solver.continuity_l2() < 2.0e-6, "collocated cavity maintains conservative face flux");

    double center_distance = 1.0e9;
    cfd::fvm::Vec3 center_velocity{};
    double upper_ux = -1.0e9;
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const auto x = solver.mesh().cells()[c].center;
        const double d2 = (x.x - 0.5) * (x.x - 0.5) + (x.y - 0.5) * (x.y - 0.5);
        if (d2 < center_distance) {
            center_distance = d2;
            center_velocity = solver.velocity()[c];
        }
        if (x.y > 0.85 && std::abs(x.x - 0.5) < 0.1) upper_ux = std::max(upper_ux, solver.velocity()[c].x);
    }
    require(center_velocity.x < -0.05, "lid cavity develops clockwise recirculation at the center");
    require(upper_ux > 0.1, "lid cavity carries positive x velocity below moving lid");
}

void test_fvm_transient_taylor_green() {
    constexpr double amplitude=0.05;
    cfd::fvm::Incompressible2D coarse({24,24,1.0,1.0,1.0,1.0e-2,1.0e-4,600,1.0e-10});
    coarse.initialize_taylor_green(amplitude);
    coarse.run(50);
    const double coarse_error=coarse.taylor_green_velocity_error(amplitude);
    require(coarse.divergence_l2()<1.0e-9,"transient FVM projection maintains incompressibility");
    require(coarse.pressure_result().converged,"transient FVM pressure CG converges");
    require(coarse_error<3.0e-4,"transient FVM Taylor-Green velocity remains accurate");

    cfd::fvm::Incompressible2D fine({48,48,1.0,1.0,1.0,1.0e-2,2.5e-5,900,1.0e-10});
    fine.initialize_taylor_green(amplitude);
    fine.run(200);
    const double fine_error=fine.taylor_green_velocity_error(amplitude);
    require(fine.divergence_l2()<1.0e-9,"refined transient FVM remains incompressible");
    require(fine_error<coarse_error*0.45,"Taylor-Green refinement reduces transient FVM velocity error");
}

void test_fvm_pressure_projection() {
    cfd::fvm::Projection2D solver({32, 28, 1.0, 1.0, 1.0, 1.0e-3, 900});
    solver.initialize_divergent(0.05);
    const double before = solver.divergence_l2();
    solver.project();
    const double after = solver.divergence_l2();
    require(before > 1.0e-2, "projection regression starts with nonzero divergence");
    require(after < before * 2.0e-4, "pressure projection strongly reduces face-flux divergence");

    solver.initialize_taylor_green(0.04);
    const double tg_before = solver.divergence_l2();
    const double energy_before = solver.kinetic_energy();
    solver.project();
    const double tg_after = solver.divergence_l2();
    require(tg_after < tg_before * 1.0e-6, "projection removes staggered Taylor-Green sampling divergence");
    require(std::abs(solver.kinetic_energy() - energy_before) / energy_before < 1.0e-6,
            "projection minimally perturbs near-divergence-free Taylor-Green velocity");
}

void test_fvm_scalar_transport() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(40, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ScalarTransport solver(std::move(mesh), {0.05, 0.1, 300, 20, 1.0e-12});
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
    require(std::sqrt(error2 / exact2) < 8.0e-4, "implicit PolyMesh scalar diffusion reaches linear Dirichlet solution");
    require(solver.linear_result().converged, "scalar transport ILU-GMRES converges");

    auto adv_mesh = cfd::fvm::make_cartesian_hexa_mesh(50, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ScalarTransport adv(std::move(adv_mesh), {0.005, 1.0e-4, 300, 20, 1.0e-11});
    adv.set_boundary("left", cfd::fvm::ScalarBoundaryType::fixedValue, 1.0);
    adv.set_boundary("right", cfd::fvm::ScalarBoundaryType::zeroGradient);
    std::vector<cfd::fvm::Vec3> velocity(adv.mesh().cell_count(), {0.2, 0.0, 0.0});
    std::vector<cfd::fvm::Vec3> boundary_velocity(adv.mesh().face_count(), {0.2, 0.0, 0.0});
    adv.set_face_flux(cfd::fvm::face_flux_from_velocity(adv.mesh(), velocity, boundary_velocity));
    adv.initialize(0.0);
    adv.run(120);
    require(adv.minimum() >= -1.0e-12 && adv.maximum() <= 1.0 + 1.0e-12,
            "implicit upwind scalar advection remains bounded");
    require(adv.values().front() > adv.values().back(), "scalar advection carries inlet concentration downstream");
}

void test_fvm_zero() {
    cfd::fvm::Diffusion2D solver({32, 24, 1.0, 1.0, 1.0});
    solver.set_dirichlet(0.0, 0.0, 0.0, 0.0);
    solver.set_source(0.0);
    const double residual = solver.iterate(50);
    require(residual == 0.0, "FVM zero solution");
}

void test_fem_poisson() {
    cfd::fem::Poisson1D solver(64);
    solver.solve([](double) { return 1.0; });
    double max_error = 0.0;
    for (std::size_t i = 0; i < solver.u().size(); ++i) {
        const double x = solver.x()[i];
        const double exact = 0.5 * x * (1.0 - x);
        max_error = std::max(max_error, std::abs(solver.u()[i] - exact));
    }
    require(max_error < 2.0e-4, "FEM Poisson manufactured solution");
}

void test_fem_reference_elements() {
    using cfd::fem::ElementType;
    using cfd::fem::Point3;
    using cfd::fem::ReferencePoint;
    struct Case { ElementType type; ReferencePoint point; double reference_measure; };
    const std::array<Case,7> cases{{
        {ElementType::line2,{0.2,0.0,0.0},2.0},
        {ElementType::tri3,{0.2,0.3,0.0},0.5},
        {ElementType::quad4,{0.2,-0.3,0.0},4.0},
        {ElementType::tet4,{0.1,0.2,0.15},1.0/6.0},
        {ElementType::hex8,{0.1,-0.2,0.3},8.0},
        {ElementType::prism6,{0.2,0.3,-0.1},1.0},
        {ElementType::pyramid5,{0.05,-0.08,0.2},4.0/3.0}
    }};
    for (const auto& c : cases) {
        const auto shape=cfd::fem::evaluate_shape(c.type,c.point);
        double sum=0.0,gx=0.0,gy=0.0,gz=0.0;
        for(std::size_t i=0;i<shape.value.size();++i){sum+=shape.value[i];gx+=shape.gradient_reference[i].x;gy+=shape.gradient_reference[i].y;gz+=shape.gradient_reference[i].z;}
        require(std::abs(sum-1.0)<2.0e-14,"FEM reference shape functions form partition of unity");
        require(std::abs(gx)+std::abs(gy)+std::abs(gz)<5.0e-14,"FEM reference shape gradients sum to zero");
        double weight_sum=0.0; for(const auto&q:cfd::fem::gaussian_quadrature(c.type,2U)) weight_sum+=q.weight;
        require(std::abs(weight_sum-c.reference_measure)<2.0e-13,"FEM Gaussian quadrature integrates reference constant exactly");
    }

    const std::array<Point3,2> line{{{0.0,0.0,0.0},{2.0,0.0,0.0}}};
    const std::array<Point3,3> tri{{{0.0,0.0,0.0},{2.0,0.0,0.0},{0.0,3.0,0.0}}};
    const std::array<Point3,4> quad{{{0.0,0.0,0.0},{2.0,0.0,0.0},{2.0,4.0,0.0},{0.0,4.0,0.0}}};
    const std::array<Point3,4> tet{{{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}};
    const std::array<Point3,8> hex{{{0.0,0.0,0.0},{2.0,0.0,0.0},{2.0,4.0,0.0},{0.0,4.0,0.0},{0.0,0.0,6.0},{2.0,0.0,6.0},{2.0,4.0,6.0},{0.0,4.0,6.0}}};
    const std::array<Point3,6> prism{{{0.0,0.0,-1.0},{1.0,0.0,-1.0},{0.0,1.0,-1.0},{0.0,0.0,1.0},{1.0,0.0,1.0},{0.0,1.0,1.0}}};
    const std::array<Point3,5> pyramid{{{-1.0,-1.0,0.0},{1.0,-1.0,0.0},{1.0,1.0,0.0},{-1.0,1.0,0.0},{0.0,0.0,1.0}}};
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::line2,line,{0.0,0.0,0.0}).measure-1.0)<1.0e-14,"Line2 isoparametric metric");
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::tri3,tri,{0.2,0.2,0.0}).measure-6.0)<1.0e-14,"Tri3 isoparametric Jacobian");
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::quad4,quad,{0.1,-0.1,0.0}).measure-2.0)<1.0e-14,"Quad4 isoparametric Jacobian");
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::tet4,tet,{0.1,0.2,0.1}).measure-1.0)<1.0e-14,"Tet4 isoparametric Jacobian");
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::hex8,hex,{0.1,-0.2,0.3}).measure-6.0)<1.0e-13,"Hex8 isoparametric Jacobian");
    require(std::abs(cfd::fem::evaluate_isoparametric(ElementType::prism6,prism,{0.2,0.2,0.0}).measure-1.0)<1.0e-14,"Prism6 isoparametric Jacobian");
    require(cfd::fem::evaluate_isoparametric(ElementType::pyramid5,pyramid,{0.05,-0.05,0.2}).measure>0.0,"Pyramid5 isoparametric Jacobian positive");
}

void test_fem_matrix_free_laplace() {
    const auto mesh=cfd::fem::make_rectangle_tri_mesh(8,7,1.3,0.9);
    const auto A=cfd::fem::assemble_tri3_laplace_matrix(mesh,[](cfd::fem::Node2 p){return 1.0+0.2*p.x+0.1*p.y;});
    std::vector<double>x(mesh.node_count()),ya(mesh.node_count()),ym(mesh.node_count());
    for(std::size_t i=0;i<x.size();++i)x[i]=std::sin(0.37*static_cast<double>(i))+0.1*std::cos(0.11*static_cast<double>(i));
    A.multiply(x,ya);
    cfd::fem::apply_tri3_laplace_matrix_free(mesh,[](cfd::fem::Node2 p){return 1.0+0.2*p.x+0.1*p.y;},x,ym);
    double maxe=0.0; for(std::size_t i=0;i<x.size();++i)maxe=std::max(maxe,std::abs(ya[i]-ym[i]));
    require(maxe<2.0e-12,"matrix-free Tri3 Laplace matches assembled CSR action");
}

void test_fem_mixed_boundary_diffusion() {
    cfd::fem::ScalarDiffusion2D solver(cfd::fem::make_rectangle_tri_mesh(18,16));
    using Type=cfd::fem::ScalarBoundaryType;
    solver.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2 p){return 2.0*p.y;},{}});
    solver.set_boundary(2,{Type::dirichlet,[](cfd::fem::Node2 p){return p.x;},{}});
    solver.set_boundary(1,{Type::neumann,[](cfd::fem::Node2){return 1.0;},{}});
    solver.set_boundary(3,{Type::robin,[](cfd::fem::Node2 p){return 2.0+3.0*(p.x+2.0);},[](cfd::fem::Node2){return 3.0;}});
    solver.solve([](cfd::fem::Node2){return 1.0;},[](cfd::fem::Node2){return 0.0;},[](cfd::fem::Node2){return 0.0;});
    require(solver.linear_result().converged,"mixed-BC scalar FEM PCG converges");
    double maxe=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){const auto p=solver.mesh().nodes[i];maxe=std::max(maxe,std::abs(solver.solution()[i]-(p.x+2.0*p.y)));}
    require(maxe<2.0e-10,"Dirichlet/Neumann/Robin FEM reproduces exact affine field");
}

void test_fem_electromagnetics2d() {
    using Type=cfd::fem::ScalarBoundaryType;
    auto mesh=cfd::fem::make_rectangle_tri_mesh(20,12,2.0,1.0);
    cfd::fem::Electrostatics2D electro(mesh);
    electro.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    electro.set_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 3.0;},{}});
    electro.solve([](cfd::fem::Node2){return 2.5;},[](cfd::fem::Node2){return 0.0;});
    require(electro.linear_result().converged,"FEM electrostatic capacitor solve converges");
    double maxe=0.0;
    for(std::size_t i=0;i<electro.mesh().node_count();++i){const auto p=electro.mesh().nodes[i];maxe=std::max(maxe,std::abs(electro.potential()[i]-1.5*p.x));}
    require(maxe<2.0e-10,"FEM electrostatics reproduces linear parallel-plate potential");
    const auto ef=electro.element_electric_field();
    for(const auto&e:ef){require(std::abs(e.x+1.5)<2.0e-10 && std::abs(e.y)<2.0e-10,"FEM electrostatic element field matches exact capacitor field");}

    cfd::fem::DCConduction2D dc(std::move(mesh));
    dc.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    dc.set_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 4.0;},{}});
    dc.solve([](cfd::fem::Node2){return 5.0;});
    require(dc.linear_result().converged,"FEM DC conduction solve converges");
    const auto j=dc.element_current_density([](cfd::fem::Node2){return 5.0;});
    for(const auto&v:j){require(std::abs(v.x+10.0)<2.0e-9 && std::abs(v.y)<2.0e-9,"FEM DC current density obeys Ohm law");}
}

void test_fem_poisson2d() {
    constexpr double pi = 3.1415926535897932384626433832795;
    cfd::fem::Poisson2D solver(cfd::fem::make_rectangle_tri_mesh(32, 32));
    solver.solve(
        [](cfd::fem::Node2 p) { return 2.0 * pi * pi * std::sin(pi*p.x) * std::sin(pi*p.y); },
        [](cfd::fem::Node2) { return 0.0; });
    require(solver.linear_result().converged, "2-D Tri3 FEM Poisson PCG converges");
    double error2=0.0, exact2=0.0;
    for (std::size_t i=0;i<solver.mesh().node_count();++i) {
        const auto p=solver.mesh().nodes[i];
        const double exact=std::sin(pi*p.x)*std::sin(pi*p.y);
        const double e=solver.solution()[i]-exact;
        error2+=e*e; exact2+=exact*exact;
    }
    require(std::sqrt(error2/exact2) < 2.0e-3,
            "2-D Tri3 FEM Poisson converges to manufactured sine solution");
}

void test_fem_poisson3d() {
    constexpr double pi=3.1415926535897932384626433832795;
    const auto run=[](std::size_t n){
        cfd::fem::Poisson3D solver(cfd::fem::make_box_tet_mesh(n,n,n));
        solver.solve([](cfd::fem::Point3 p){return 3.0*pi*pi*std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z);},
            [](cfd::fem::Point3){return 0.0;});
        require(solver.linear_result().converged,"3-D Tet4 Poisson PCG converges");
        double e2=0.0,x2=0.0;
        for(std::size_t i=0;i<solver.mesh().node_count();++i){const auto p=solver.mesh().nodes[i];const double exact=std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z);const double e=solver.solution()[i]-exact;e2+=e*e;x2+=exact*exact;}
        return std::sqrt(e2/x2);
    };
    const double coarse=run(5U),fine=run(10U);
    require(fine<2.0e-2 && fine<0.4*coarse,"3-D Tet4 Poisson shows expected mesh-convergence trend");
}

void test_fem_heat2d() {
    constexpr double pi=3.1415926535897932384626433832795;
    constexpr double alpha=0.2;
    cfd::fem::Heat2D solver(cfd::fem::make_rectangle_tri_mesh(28,28), {alpha,1.0,1.0e-4,2000,1.0e-10});
    solver.initialize([](cfd::fem::Node2 p){return std::sin(pi*p.x)*std::sin(pi*p.y);});
    solver.set_dirichlet([](cfd::fem::Node2,double){return 0.0;});
    solver.run(10);
    require(solver.linear_result().converged,"2-D transient FEM heat solve converges");
    const double decay=std::exp(-2.0*pi*pi*alpha*solver.time());
    double e2=0.0, exact2=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){
        const auto p=solver.mesh().nodes[i]; const double exact=decay*std::sin(pi*p.x)*std::sin(pi*p.y);
        const double e=solver.temperature()[i]-exact; e2+=e*e; exact2+=exact*exact;
    }
    require(std::sqrt(e2/exact2)<3.0e-3,"implicit P1 FEM heat matches analytical Fourier decay");
}

void test_fem_elasticity2d_patch() {
    constexpr double strain=1.0e-3, nu=0.3;
    cfd::fem::Elasticity2D solver(cfd::fem::make_rectangle_tri_mesh(12,10), {200.0e9,nu,cfd::fem::ElasticityMode2D::planeStress,4000,1.0e-11});
    solver.set_dirichlet([](cfd::fem::Node2 p){return p.x==0.0||p.x==1.0||p.y==0.0||p.y==1.0;},
        [](cfd::fem::Node2 p){return cfd::fem::Displacement2{strain*p.x,-nu*strain*p.y};});
    solver.solve(); require(solver.linear_result().converged,"Tri3 plane-stress elasticity patch solve converges");
    double maxe=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){ const auto p=solver.mesh().nodes[i]; const auto u=solver.displacement()[i];
        maxe=std::max(maxe,std::abs(u.x-strain*p.x)); maxe=std::max(maxe,std::abs(u.y+nu*strain*p.y)); }
    require(maxe<2.0e-12,"Tri3 linear elasticity reproduces affine plane-stress patch exactly");
}

void test_fem_axisymmetric_elasticity_patch() {
    constexpr double strain=2.0e-4;
    auto mesh=cfd::fem::make_rectangle_tri_mesh(10,9,1.0,1.0);
    cfd::fem::AxisymmetricElasticity solver(std::move(mesh),{150.0e9,0.27,5000,1.0e-11});
    solver.set_dirichlet([](cfd::fem::Node2 p){return p.x==0.0||p.x==1.0||p.y==0.0||p.y==1.0;},
        [](cfd::fem::Node2 p){return cfd::fem::Displacement2{strain*p.x,strain*p.y};});
    solver.solve();
    require(solver.linear_result().converged,"axisymmetric elasticity patch solve converges");
    double maxe=0.0;
    for(std::size_t i=0;i<solver.mesh().node_count();++i){const auto p=solver.mesh().nodes[i];const auto u=solver.displacement()[i];maxe=std::max(maxe,std::abs(u.x-strain*p.x));maxe=std::max(maxe,std::abs(u.y-strain*p.y));}
    require(maxe<5.0e-12,"axisymmetric elasticity reproduces uniform dilatation patch");
}

double fem_element_l2_error(const cfd::fem::Mesh2D& mesh,
                            const std::vector<double>& nodal,
                            const std::function<double(cfd::fem::Node2)>& exact) {
    double error2=0.0,measure=0.0;
    for(const auto&tri:mesh.triangles){
        const auto&a=mesh.nodes[tri.node[0]],&b=mesh.nodes[tri.node[1]],&c=mesh.nodes[tri.node[2]];
        const double area=0.5*std::abs((b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
        const cfd::fem::Node2 p{(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
        const double uh=(nodal[tri.node[0]]+nodal[tri.node[1]]+nodal[tri.node[2]])/3.0;
        const double e=uh-exact(p); error2+=area*e*e; measure+=area;
    }
    return std::sqrt(error2/measure);
}

void test_fem_nonlinear_poisson_newton() {
    const double pi=std::numbers::pi_v<double>;
    const double beta=3.0;
    cfd::fem::NonlinearPoisson2DConfig cfg;
    cfg.cubic_coefficient=beta;
    cfg.newton.absolute_tolerance=1.0e-10;
    cfg.newton.relative_tolerance=1.0e-9;
    cfg.newton.linear_relative_tolerance=1.0e-11;
    cfd::fem::NonlinearPoisson2D solver(cfd::fem::make_rectangle_tri_mesh(24,24),cfg);
    const auto exact=[&](cfd::fem::Node2 p){return std::sin(pi*p.x)*std::sin(pi*p.y);};
    solver.solve([&](cfd::fem::Node2 p){const double u=exact(p);return 2.0*pi*pi*u+beta*u*u*u;},
                 [](cfd::fem::Node2){return 0.0;});
    const double error=fem_element_l2_error(solver.mesh(),solver.solution(),exact);
    require(solver.nonlinear_result().converged,"nonlinear FEM Newton solve converges");
    require(solver.nonlinear_result().iterations<12U,"nonlinear FEM Newton convergence is rapid");
    require(error<4.0e-3,"nonlinear FEM manufactured solution accuracy");
}

void test_fem_adaptive_refinement() {
    const double pi=std::numbers::pi_v<double>;
    const auto exact=[&](cfd::fem::Node2 p){return std::sin(pi*p.x)*std::sin(pi*p.y);};
    const auto source=[&](cfd::fem::Node2 p){return 2.0*pi*pi*exact(p);};
    auto coarse_mesh=cfd::fem::make_rectangle_tri_mesh(8,8);
    cfd::fem::Poisson2D coarse(coarse_mesh); coarse.solve(source,[](cfd::fem::Node2){return 0.0;});
    const double coarse_error=fem_element_l2_error(coarse.mesh(),coarse.solution(),exact);
    const auto estimate=cfd::fem::estimate_poisson_error_tri3(coarse.mesh(),coarse.solution(),source);
    const auto marked=cfd::fem::mark_dorfler(estimate.element_indicator,0.55);
    const auto marked_count=static_cast<std::size_t>(std::count(marked.begin(),marked.end(),static_cast<unsigned char>(1U)));
    require(marked_count>0U&&marked_count<marked.size(),"Dorfler marking selects a strict element subset");
    auto refined_mesh=cfd::fem::refine_tri3_longest_edges(coarse.mesh(),marked);
    require(refined_mesh.element_count()>coarse.mesh().element_count(),"adaptive FEM refinement adds elements");
    require(refined_mesh.element_count()<4U*coarse.mesh().element_count(),"adaptive refinement remains below full red refinement");
    cfd::fem::Poisson2D refined(refined_mesh); refined.solve(source,[](cfd::fem::Node2){return 0.0;});
    const double refined_error=fem_element_l2_error(refined.mesh(),refined.solution(),exact);
    require(refined_error<coarse_error,"adaptive FEM refinement reduces manufactured-solution error");
}

void test_fem_darcy_flow() {
    using Type=cfd::fem::ScalarBoundaryType;
    constexpr double permeability=2.0e-12,viscosity=1.0e-3,p_left=1.0e5,p_right=0.0,length=2.0;
    cfd::fem::Darcy2DConfig cfg;cfg.viscosity=viscosity;cfg.density=1000.0;
    cfd::fem::Darcy2D solver(cfd::fem::make_rectangle_tri_mesh(24,12,length,1.0),cfg);
    solver.set_boundary(0,{Type::dirichlet,[&](cfd::fem::Node2){return p_left;},{}});
    solver.set_boundary(1,{Type::dirichlet,[&](cfd::fem::Node2){return p_right;},{}});
    solver.solve([](cfd::fem::Node2){return permeability;});
    const auto velocity=solver.element_velocity([](cfd::fem::Node2){return permeability;});
    const double expected=permeability/viscosity*(p_left-p_right)/length;
    double mean_x=0.0,max_y=0.0;for(const auto&v:velocity){mean_x+=v.x;max_y=std::max(max_y,std::abs(v.y));}mean_x/=static_cast<double>(velocity.size());
    require(std::abs(mean_x-expected)/expected<2.0e-6,"Darcy FEM reproduces analytical uniform flow");
    require(max_y<5.0e-6*expected,"Darcy FEM impermeable walls suppress transverse flow");
}

void test_fem_magnetostatics() {
    using Type=cfd::fem::ScalarBoundaryType;
    const double pi=std::numbers::pi_v<double>,nu=2.5;
    cfd::fem::Magnetostatics2D solver(cfd::fem::make_rectangle_tri_mesh(24,24));
    for(int patch=0;patch<4;++patch) solver.set_boundary(patch,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    const auto exact=[&](cfd::fem::Node2 p){return std::sin(pi*p.x)*std::sin(pi*p.y);};
    solver.solve([&](cfd::fem::Node2){return nu;},[&](cfd::fem::Node2 p){return 2.0*pi*pi*nu*exact(p);});
    const double error=fem_element_l2_error(solver.mesh(),solver.vector_potential(),exact);
    require(error<4.0e-3,"magnetostatic A_z manufactured solution accuracy");
    const auto b=solver.element_flux_density();
    double max_mag=0.0;for(const auto&v:b)max_mag=std::max(max_mag,std::hypot(v.x,v.y));
    require(max_mag>2.0,"magnetostatic FEM reconstructs nonzero magnetic flux density");
}

void test_fem_modal_bar() {
    cfd::fem::BarModal1DConfig cfg;cfg.elements=80U;cfg.length=1.7;cfg.young_modulus=70.0e9;cfg.density=2700.0;
    cfg.eigen.max_outer_iterations=120U;cfg.eigen.relative_tolerance=2.0e-8;cfg.eigen.linear_relative_tolerance=1.0e-12;
    cfd::fem::BarModal1D solver(cfg);
    const auto modes=solver.solve(2U);
    const double wave=std::sqrt(cfg.young_modulus/cfg.density);
    for(std::size_t m=0;m<modes.size();++m){
        const double exact=(2.0*static_cast<double>(m)+1.0)*std::numbers::pi_v<double>*wave/(2.0*cfg.length);
        require(std::abs(modes[m].angular_frequency-exact)/exact<1.5e-3,"bar FEM eigenfrequency matches analytical fixed-free mode");
        require(modes[m].eigenpair.converged,"generalized FEM eigenpair converges");
    }
}

void test_fdtd_finite() {
    cfd::fdtd::Maxwell1D solver({256, 1.0e-3, 0.9, 1.0, 1.0});
    solver.initialize_gaussian();
    solver.step(80);
    const double energy = solver.energy();
    require(std::isfinite(energy) && energy > 0.0, "FDTD stable finite energy");
}


void test_fdtd_material_mur_monitors() {
    cfd::fdtd::Maxwell1D material({128,1.0e-3,0.9,1.0,1.0,cfd::fdtd::Boundary1D::mur1});
    const double c0=material.local_wave_speed(10U);
    material.set_material(60U,100U,4.0,0.02);
    require(std::abs(material.local_wave_speed(70U)/c0-0.5)<2.0e-13,"FDTD material preprocessing changes dielectric wave speed");
    material.set_hard_source(40U,0.25);
    require(material.electric()[40U]==0.25,"FDTD hard source sets electric field exactly");
    material.add_soft_source(40U,0.1);
    require(std::abs(material.electric()[40U]-0.35)<1.0e-15,"FDTD soft source adds to electric field");
    material.step(20U);
    require(std::isfinite(material.energy())&&material.energy()>0.0,"lossy dielectric FDTD remains finite");

    cfd::fdtd::Maxwell1D pec({400,1.0e-3,0.95,1.0,1.0,cfd::fdtd::Boundary1D::pec});
    cfd::fdtd::Maxwell1D mur({400,1.0e-3,0.95,1.0,1.0,cfd::fdtd::Boundary1D::mur1});
    cfd::fdtd::Maxwell1D pml({400,1.0e-3,0.95,1.0,1.0,cfd::fdtd::Boundary1D::pml,24U,3.0,1.0e-8});
    pec.initialize_gaussian(0.35,0.025); mur.initialize_gaussian(0.35,0.025); pml.initialize_gaussian(0.35,0.025);
    pec.step(700U); mur.step(700U); pml.step(700U);
    require(mur.energy()<1.0e-4*pec.energy(),"first-order Mur boundary strongly reduces reflected pulse energy versus PEC");
    require(pml.energy()<mur.energy(),"1-D matched-loss PML absorbs pulse more strongly than Mur in the regression case");

    cfd::fdtd::TimeProbe probe;
    cfd::fdtd::DftMonitor dft(10.0);
    constexpr double dt=1.0e-3;
    for(std::size_t i=0;i<=1000U;++i){const double t=dt*static_cast<double>(i);const double v=std::sin(2.0*std::numbers::pi*10.0*t);probe.sample(t,v);dft.sample(t,v);}
    require(probe.time().size()==1001U&&probe.value().size()==1001U,"FDTD time probe records samples");
    require(std::abs(dft.amplitude(1.0)-1.0)<2.0e-3,"FDTD DFT monitor recovers sinusoid amplitude");
}


void test_fdtd_cpml_tfsf_lumped() {
    using cfd::fdtd::Boundary1D;
    using cfd::fdtd::Maxwell1D;

    Maxwell1D pec({400,1.0e-3,0.95,1.0,1.0,Boundary1D::pec});
    Maxwell1D cpml({400,1.0e-3,0.95,1.0,1.0,Boundary1D::cpml,24U,3.0,1.0e-8,5.0,0.0});
    pec.initialize_gaussian(0.35,0.025);
    cpml.initialize_gaussian(0.35,0.025);
    const double cpml_initial=cpml.energy();
    pec.step(700U);
    cpml.step(700U);
    require(std::isfinite(cpml.energy()) && cpml.energy() < 1.0e-4*cpml_initial,
            "1-D CPML strongly attenuates the outgoing Gaussian pulse");
    require(cpml.energy() < 1.0e-4*pec.energy(),
            "1-D CPML leaves far less residual field energy than PEC");

    Maxwell1D tfsf({500,1.0e-3,0.95,1.0,1.0,Boundary1D::cpml,32U,3.0,1.0e-10,6.0,0.0});
    const double dt=tfsf.dt();
    const double t0=35.0*dt, tau=10.0*dt;
    tfsf.set_tfsf_source(120U,[=](double time){
        const double q=(time-t0)/tau;
        return std::exp(-q*q);
    });
    double scattered_max=0.0,total_max=0.0;
    for(std::size_t n=0;n<120U;++n){
        tfsf.step();
        scattered_max=std::max(scattered_max,std::abs(tfsf.electric()[80U]));
        total_max=std::max(total_max,std::abs(tfsf.electric()[160U]));
    }
    require(total_max>0.5,"TFSF launches a finite +x incident field into the total-field region");
    require(scattered_max/total_max<1.0e-4,
            "homogeneous 1-D TFSF keeps the exterior scattered-field region quiet");

    Maxwell1D rlc({128,1.0e-3,0.5,1.0,1.0,Boundary1D::mur1});
    constexpr std::size_t cell=64U;
    constexpr double length=1.0e-3,area=1.0e-6,inductance=1.0e-6,field=0.2;
    rlc.set_parallel_lumped_rlc(cell,length,area,std::numeric_limits<double>::infinity(),inductance,0.0);
    const double drive=length/(inductance*area);
    double expected_current_density=0.0;
    for(std::size_t n=0;n<20U;++n){
        rlc.set_hard_source(cell,field);
        rlc.step();
        expected_current_density+=rlc.dt()*drive*field;
    }
    require(std::abs(rlc.lumped_inductor_current_density(cell)-expected_current_density)
                < 1.0e-12*std::max(1.0,std::abs(expected_current_density)),
            "parallel lumped inductor current follows the analytical constant-voltage ramp");

    Maxwell1D resistor({256,1.0e-3,0.8,1.0,1.0,Boundary1D::mur1});
    Maxwell1D reference({256,1.0e-3,0.8,1.0,1.0,Boundary1D::mur1});
    resistor.set_parallel_lumped_rlc(128U,1.0e-3,1.0e-6,25.0,
                                     std::numeric_limits<double>::infinity(),0.0);
    resistor.initialize_gaussian(0.5,0.025);
    reference.initialize_gaussian(0.5,0.025);
    resistor.step(120U);
    reference.step(120U);
    require(std::isfinite(resistor.energy()) && resistor.energy()<0.9*reference.energy(),
            "parallel lumped resistor dissipates field energy");

    Maxwell1D capacitor({256,1.0e-3,0.8,1.0,1.0,Boundary1D::mur1});
    capacitor.set_parallel_lumped_rlc(128U,1.0e-3,1.0e-6,
                                      std::numeric_limits<double>::infinity(),
                                      std::numeric_limits<double>::infinity(),2.0e-12);
    capacitor.initialize_gaussian(0.5,0.025);
    capacitor.step(120U);
    require(std::isfinite(capacitor.energy()) && std::isfinite(capacitor.electric()[128U]),
            "parallel lumped capacitor remains numerically stable");
}

void test_sparse_krylov_solvers() {
    cfd::core::CsrBuilder spd_builder(3, 3);
    spd_builder.add(0, 0, 4.0); spd_builder.add(0, 1, -1.0);
    spd_builder.add(1, 0, -1.0); spd_builder.add(1, 1, 4.0); spd_builder.add(1, 2, -1.0);
    spd_builder.add(2, 1, -1.0); spd_builder.add(2, 2, 3.0);
    const auto spd = spd_builder.build();
    const std::vector<double> exact{1.0, 2.0, 3.0};
    std::vector<double> rhs(3, 0.0);
    spd.multiply(exact, rhs);
    std::vector<double> x(3, 0.0);
    const auto diagonal = spd.diagonal();
    const cfd::core::JacobiPreconditioner jacobi(diagonal);
    cfd::core::KrylovWorkspace workspace;
    const auto pcg = cfd::core::preconditioned_conjugate_gradient(
        rhs, x,
        [&](std::span<const double> in, std::span<double> out) { spd.multiply(in, out); },
        [&](std::span<const double> in, std::span<double> out) { jacobi(in, out); },
        workspace, 50, 1.0e-12);
    require(pcg.converged, "Jacobi-PCG converges on SPD CSR system");
    for (std::size_t i = 0; i < x.size(); ++i) {
        require(std::abs(x[i] - exact[i]) < 1.0e-10, "Jacobi-PCG solution accuracy");
    }

    cfd::core::CsrBuilder nonsym_builder(3, 3);
    nonsym_builder.add(0, 0, 4.0); nonsym_builder.add(0, 1, 1.0);
    nonsym_builder.add(1, 0, -2.0); nonsym_builder.add(1, 1, 3.0); nonsym_builder.add(1, 2, 1.0);
    nonsym_builder.add(2, 1, -1.0); nonsym_builder.add(2, 2, 2.0);
    const auto nonsym = nonsym_builder.build();
    const std::vector<double> exact2{1.0, -2.0, 0.5};
    std::vector<double> rhs2(3, 0.0);
    nonsym.multiply(exact2, rhs2);
    std::vector<double> x2(3, 0.0);
    const cfd::core::JacobiPreconditioner jacobi2(nonsym.diagonal());
    const auto bicg = cfd::core::bicgstab(
        rhs2, x2,
        [&](std::span<const double> in, std::span<double> out) { nonsym.multiply(in, out); },
        [&](std::span<const double> in, std::span<double> out) { jacobi2(in, out); },
        workspace, 100, 1.0e-12);
    require(bicg.converged, "Jacobi-BiCGStab converges on nonsymmetric CSR system");
    for (std::size_t i = 0; i < x2.size(); ++i) {
        require(std::abs(x2[i] - exact2[i]) < 1.0e-9, "BiCGStab solution accuracy");
    }

    std::vector<double> x3(3, 0.0);
    const cfd::core::Ilu0Preconditioner ilu0(nonsym);
    const auto gmres = cfd::core::restarted_gmres(
        rhs2, x3,
        [&](std::span<const double> in, std::span<double> out) { nonsym.multiply(in, out); },
        [&](std::span<const double> in, std::span<double> out) { ilu0(in, out); },
        50, 4, 1.0e-12);
    require(gmres.converged, "ILU(0)-GMRES converges on nonsymmetric CSR system");
    for (std::size_t i = 0; i < x3.size(); ++i) {
        require(std::abs(x3[i] - exact2[i]) < 1.0e-9, "GMRES solution accuracy");
    }
}

void test_fvm_bounded_linear_reconstruction() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(16, 2, 2, 1.0, 1.0, 1.0);
    std::vector<double> linear(mesh.cell_count(), 0.0);
    for (std::size_t c = 0; c < mesh.cell_count(); ++c) linear[c] = mesh.cells()[c].center.x;
    std::vector<double> boundary(mesh.face_count(), 0.0);
    std::vector<cfd::fvm::Vec3> velocity(mesh.cell_count(), {1.0, 0.0, 0.0});
    std::vector<cfd::fvm::Vec3> boundary_velocity(mesh.face_count(), {1.0, 0.0, 0.0});
    for (std::size_t f = 0; f < mesh.face_count(); ++f) boundary[f] = mesh.faces()[f].center.x;
    const auto flux = cfd::fvm::face_flux_from_velocity(mesh, velocity, boundary_velocity);
    const auto faces = cfd::fvm::interpolate_scalar_to_faces(
        mesh, linear, cfd::fvm::FaceInterpolationScheme::bounded_linear, flux, boundary);
    double max_error = 0.0;
    for (std::size_t f = 0; f < mesh.face_count(); ++f) {
        max_error = std::max(max_error, std::abs(faces[f] - mesh.faces()[f].center.x));
    }
    require(max_error < 2.0e-13, "bounded linear reconstruction preserves linear field");

    std::vector<double> step(mesh.cell_count(), 0.0);
    for (std::size_t c = 0; c < mesh.cell_count(); ++c) step[c] = mesh.cells()[c].center.x < 0.5 ? 0.0 : 1.0;
    for (std::size_t f = 0; f < mesh.face_count(); ++f) boundary[f] = mesh.faces()[f].center.x < 0.5 ? 0.0 : 1.0;
    const auto bounded = cfd::fvm::interpolate_scalar_to_faces(
        mesh, step, cfd::fvm::FaceInterpolationScheme::bounded_linear, flux, boundary);
    for (double value : bounded) {
        require(value >= -1.0e-14 && value <= 1.0 + 1.0e-14,
                "bounded linear reconstruction does not overshoot a discontinuity");
    }
}

void test_chemistry_kinetics() {
    cfd::chemistry::ReactionNetwork network({{"A", 0.010, 0}, {"B", 0.010, 0}});
    cfd::chemistry::ElementaryReaction reaction;
    reaction.reactants.push_back({0U, 1.0});
    reaction.products.push_back({1U, 1.0});
    reaction.forward = {2.0, 0.0, 0.0};
    network.add_reaction(reaction);
    const std::vector<double> c{3.0, 0.0};
    const auto rates = network.reaction_rates(c, 300.0);
    const auto source = network.source_terms(c, 300.0);
    require(std::abs(rates[0] - 6.0) < 1.0e-12, "mass-action elementary reaction rate");
    require(std::abs(source[0] + 6.0) < 1.0e-12 && std::abs(source[1] - 6.0) < 1.0e-12,
            "reaction source stoichiometry");
    require(std::abs(source[0] + source[1]) < 1.0e-12, "reaction network conserves one-to-one species count");
}

void test_electrochemistry_relations() {
    const double e0 = 0.77;
    require(std::abs(cfd::electrochemistry::nernst_potential(e0, 298.15, 2.0, 1.0) - e0) < 1.0e-15,
            "Nernst potential equals standard potential at unit reaction quotient");
    const double jp = cfd::electrochemistry::butler_volmer_current_density(2.0, 0.01, 298.15, 2.0);
    const double jm = cfd::electrochemistry::butler_volmer_current_density(2.0, -0.01, 298.15, 2.0);
    require(std::abs(jp + jm) < 1.0e-12, "symmetric Butler-Volmer is odd in overpotential");
    require(std::abs(cfd::electrochemistry::butler_volmer_current_density(2.0, 0.0, 298.15, 2.0)) < 1.0e-15,
            "Butler-Volmer zero-overpotential current");

    cfd::electrochemistry::CorrosionCell1DConfig cfg;
    cfg.metal_potential = 1.0e-3;
    cfg.equilibrium_potential = 0.0;
    cfg.bulk_electrolyte_potential = 0.0;
    cfg.exchange_current_density = 2.0;
    cfg.electrolyte_length = 1.0e-3;
    cfg.electrolyte_conductivity = 5.0;
    cfg.electrons = 2.0;
    const auto result = cfd::electrochemistry::solve_corrosion_cell_1d(cfg);
    require(result.converged && result.current_density > 0.0, "1-D corrosion cell converges to anodic current");
    const double conductance = cfg.exchange_current_density * cfg.electrons
        * cfd::electrochemistry::faraday_constant
        / (cfd::chemistry::gas_constant * cfg.temperature)
        * (cfg.anodic_transfer + cfg.cathodic_transfer);
    const double linear_current = conductance * cfg.metal_potential
        / (1.0 + conductance * cfg.electrolyte_length / cfg.electrolyte_conductivity);
    require(std::abs(result.current_density - linear_current) / linear_current < 3.0e-3,
            "corrosion cell matches small-overpotential linearized Butler-Volmer/ohmic solution");
    require(result.penetration_rate > 0.0, "Faraday law returns positive anodic penetration rate");
}

void test_mixed_potential_solver() {
    std::vector<cfd::electrochemistry::ElectrodeReaction> reactions{
        {"anodic", 0.0, 2.0, 1.0, 1.0, 0.5, 0.5},
        {"cathodic", 0.2, 2.0, 1.0, 1.0, 0.5, 0.5},
    };
    const auto result = cfd::electrochemistry::solve_mixed_potential(reactions, 298.15);
    require(result.converged, "mixed-potential galvanic solve converges");
    require(std::abs(result.potential - 0.1) < 1.0e-12,
            "symmetric galvanic pair has midpoint mixed potential");
    require(result.reaction_current.size() == 2U
            && result.reaction_current[0] > 0.0 && result.reaction_current[1] < 0.0,
            "mixed potential balances anodic and cathodic branches");
    require(std::abs(result.reaction_current[0] + result.reaction_current[1]) < 1.0e-10,
            "open-circuit galvanic mixed potential has zero net current");
}

void test_nernst_planck_transport() {
    constexpr double pi = 3.1415926535897932384626433832795;
    cfd::electrochemistry::NernstPlanck1D solver({128, 1.0, 298.15, 5.0e-4, 0.0,
                                                  cfd::electrochemistry::TransportBoundary1D::periodic});
    solver.set_potential([](double) { return 0.0; });
    const std::size_t species = solver.add_species("neutral", 0, 1.0e-2, [](double x) {
        return 1.0 + 0.1 * std::cos(2.0 * pi * x);
    });
    const double amount0 = solver.total_amount(species);
    solver.step(200);
    const double amount1 = solver.total_amount(species);
    require(std::abs(amount1 - amount0) < 2.0e-12, "Nernst-Planck finite-volume transport conserves periodic species amount");
    double amplitude = 0.0;
    const auto& concentration = solver.species()[species].concentration;
    for (std::size_t i = 0; i < concentration.size(); ++i) {
        const double x = (static_cast<double>(i) + 0.5) * solver.dx();
        amplitude += (concentration[i] - 1.0) * std::cos(2.0 * pi * x);
    }
    amplitude *= 2.0 / static_cast<double>(concentration.size());
    const double exact_amplitude = 0.1 * std::exp(-1.0e-2 * 4.0 * pi * pi * solver.time());
    require(std::abs(amplitude - exact_amplitude) / exact_amplitude < 1.5e-3,
            "Nernst-Planck diffusion matches analytical Fourier decay");

    cfd::electrochemistry::NernstPlanck1D charged({96, 1.0, 298.15, 1.0e-6, 0.0,
                                                   cfd::electrochemistry::TransportBoundary1D::periodic});
    charged.set_potential([](double x) { return 1.0e-3 * std::sin(2.0 * pi * x); });
    const std::size_t ion = charged.add_species("cation", 1, 1.0e-3, 1.0);
    const double ion0 = charged.total_amount(ion);
    charged.step(50);
    require(std::abs(charged.total_amount(ion) - ion0) < 2.0e-12,
            "electromigration transport remains conservative");
    require(charged.minimum_concentration(ion) > 0.0, "electromigration test remains positive at stable timestep");
}

void test_nernst_planck_polymesh() {
    constexpr double pi = 3.1415926535897932384626433832795;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(64, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(std::move(mesh), {298.15, 2.0e-4});
    const auto neutral = solver.add_species("neutral", 0, 1.0e-2, [](cfd::fvm::Vec3 x) {
        return 1.0 + 0.1 * std::cos(pi * x.x);
    });
    const auto ion = solver.add_species("cation", 1, 1.0e-3, 1.0);
    solver.set_potential([](cfd::fvm::Vec3 x) { return 1.0e-4 * std::cos(pi * x.x); });
    const double neutral0 = solver.total_amount(neutral);
    const double ion0 = solver.total_amount(ion);
    solver.step(250);
    require(std::abs(solver.total_amount(neutral) - neutral0) < 3.0e-12,
            "PolyMesh Nernst-Planck conserves neutral species with no-flux boundaries");
    require(std::abs(solver.total_amount(ion) - ion0) < 3.0e-12,
            "PolyMesh Nernst-Planck conserves charged species with no-flux boundaries");
    double amplitude = 0.0;
    const auto& c = solver.species()[neutral].concentration;
    for (std::size_t i = 0; i < c.size(); ++i) {
        const double x = solver.mesh().cells()[i].center.x;
        amplitude += (c[i] - 1.0) * std::cos(pi * x);
    }
    amplitude *= 2.0 / static_cast<double>(c.size());
    const double exact = 0.1 * std::exp(-1.0e-2 * pi * pi * solver.time());
    require(std::abs(amplitude - exact) / exact < 1.5e-3,
            "PolyMesh Nernst-Planck matches no-flux Fourier diffusion decay");
    require(solver.minimum_concentration(ion) > 0.0,
            "PolyMesh electromigration remains positive at stable timestep");
}

void test_electroneutral_liquid_junction() {
    constexpr double R = cfd::chemistry::gas_constant;
    constexpr double F = cfd::electrochemistry::faraday_constant;
    const double temperature = 298.15;
    const double d_plus = 2.0e-3;
    const double d_minus = 1.0e-3;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(96, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(std::move(mesh), {temperature, 1.0e-5});
    auto profile = [](cfd::fvm::Vec3 x) { return 1.0 + 0.5 * x.x; };
    solver.add_species("cation", 1, d_plus, profile);
    solver.add_species("anion", -1, d_minus, profile);
    const auto result = solver.solve_electroneutral_potential(1000, 1.0e-12);
    require(result.converged, "electroneutral liquid-junction potential solve converges");

    const double factor = -(R * temperature / F) * (d_plus - d_minus) / (d_plus + d_minus);
    std::vector<double> exact(solver.mesh().cell_count(), 0.0);
    double exact_mean = 0.0;
    for (std::size_t c = 0; c < exact.size(); ++c) {
        exact[c] = factor * std::log(profile(solver.mesh().cells()[c].center));
        exact_mean += exact[c];
    }
    exact_mean /= static_cast<double>(exact.size());
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t c = 0; c < exact.size(); ++c) {
        exact[c] -= exact_mean;
        const double error = solver.potential()[c] - exact[c];
        error2 += error * error;
        exact2 += exact[c] * exact[c];
    }
    require(std::sqrt(error2 / exact2) < 2.0e-4,
            "electroneutral potential matches binary liquid-junction analytical profile");
    const auto current = solver.current_flux_faces();
    double max_current = 0.0;
    for (double value : current) max_current = std::max(max_current, std::abs(value));
    require(max_current < 5.0e-9, "electroneutral potential cancels diffusive ionic current");
}

void test_poisson_nernst_planck_potential() {
    constexpr double eps0 = 8.8541878128e-12;
    constexpr double concentration = 1.0e-15;
    constexpr double eps_r = 78.5;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(128, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(
        std::move(mesh), {298.15, 1.0e-6, eps_r});
    solver.add_species("space-charge", +1, 0.0, concentration);
    solver.set_potential_boundary("left", cfd::electrochemistry::PotentialBoundaryType::fixedPotential, 0.0);
    solver.set_potential_boundary("right", cfd::electrochemistry::PotentialBoundaryType::fixedPotential, 0.0);
    const auto result = solver.solve_poisson_potential(2000, 1.0e-12);
    require(result.converged, "Poisson-Nernst-Planck potential solve converges");

    const double rho = cfd::electrochemistry::faraday_constant * concentration;
    const double epsilon = eps0 * eps_r;
    double error2 = 0.0;
    double exact2 = 0.0;
    for (std::size_t c = 0; c < solver.mesh().cell_count(); ++c) {
        const double x = solver.mesh().cells()[c].center.x;
        const double exact = rho * x * (1.0 - x) / (2.0 * epsilon);
        const double error = solver.potential()[c] - exact;
        error2 += error * error;
        exact2 += exact * exact;
    }
    require(std::sqrt(error2 / exact2) < 1.0e-2,
            "Poisson potential matches charged-slab analytical profile");
    require(std::abs(solver.total_charge() - rho) < 1.0e-22,
            "Poisson-Nernst-Planck total charge diagnostic is conservative");

    auto neutral_mesh = cfd::fvm::make_cartesian_hexa_mesh(16, 1, 1);
    cfd::electrochemistry::NernstPlanckPolyMesh neutral(std::move(neutral_mesh), {298.15, 1.0e-6, eps_r});
    neutral.add_species("cation", +1, 1.0e-3, 1.0);
    neutral.add_species("anion", -1, 1.0e-3, 1.0);
    const auto neutral_result = neutral.solve_poisson_potential(500, 1.0e-12);
    require(neutral_result.converged, "neutral insulating PNP potential solve converges with mean-zero gauge");
    neutral.step_poisson_nernst_planck(2, 500, 1.0e-12);
    require(neutral.minimum_concentration(0) > 0.0 && neutral.minimum_concentration(1) > 0.0,
            "operator-split Poisson-Nernst-Planck step remains positive");
}

void test_scharfetter_gummel_flux() {
    constexpr double dt = 1.0e-3;
    constexpr double temperature = 298.15;
    constexpr double diffusivity = 1.0e-3;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(2, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(
        std::move(mesh),
        {temperature, dt, 78.5, cfd::electrochemistry::ElectromigrationFluxScheme::scharfetterGummel});
    const auto ion = solver.add_species("cation", +1, diffusivity, [](cfd::fvm::Vec3 x) {
        return x.x < 0.5 ? 1.0 : 2.0;
    });
    solver.set_potential([](cfd::fvm::Vec3 x) { return 0.2 * x.x; });

    const double c0 = 1.0;
    const double c1 = 2.0;
    const double delta_phi = 0.1;
    const double psi = cfd::electrochemistry::faraday_constant * delta_phi
        / (cfd::chemistry::gas_constant * temperature);
    const auto B = [](double x) {
        if (std::abs(x) < 1.0e-8) return 1.0 - 0.5 * x + x * x / 12.0;
        return x / std::expm1(x);
    };
    const double metric = 2.0; // unit-area face, half-unit cell-center spacing
    const double flux_owner_to_neighbour = diffusivity * metric * (c0 * B(psi) - c1 * B(-psi));
    const double volume = 0.5;
    solver.step();
    const auto& c = solver.species()[ion].concentration;
    require(std::abs(c[0] - (c0 - dt * flux_owner_to_neighbour / volume)) < 2.0e-12,
            "Scharfetter-Gummel owner update matches exponential-fit flux");
    require(std::abs(c[1] - (c1 + dt * flux_owner_to_neighbour / volume)) < 2.0e-12,
            "Scharfetter-Gummel neighbour update matches exponential-fit flux");
    require(std::abs(solver.total_amount(ion) - 1.5) < 2.0e-13,
            "Scharfetter-Gummel internal transport conserves species amount");
}

void test_butler_volmer_polymesh_boundary_flux() {
    constexpr double dt = 1.0e-3;
    constexpr double electrode_potential = 0.05;
    constexpr double exchange_current = 1.0;
    constexpr int electrons = 2;
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(16, 1, 1, 1.0, 1.0, 1.0);
    cfd::electrochemistry::NernstPlanckPolyMesh solver(std::move(mesh), {298.15, dt, 78.5});
    const auto metal_ion = solver.add_species("M2+", +2, 0.0, 1.0);
    solver.set_potential(0.0);
    solver.set_species_butler_volmer_boundary(
        metal_ion, "left", electrode_potential, 0.0, exchange_current, electrons);
    const double amount0 = solver.total_amount(metal_ion);
    const double scale = static_cast<double>(electrons) * cfd::electrochemistry::faraday_constant
        / (cfd::chemistry::gas_constant * 298.15);
    const double current = exchange_current
        * (std::exp(0.5 * scale * electrode_potential) - std::exp(-0.5 * scale * electrode_potential));
    const double expected_added = dt * current
        / (static_cast<double>(electrons) * cfd::electrochemistry::faraday_constant);
    solver.step();
    const double actual_added = solver.total_amount(metal_ion) - amount0;
    require(current > 0.0, "positive electrode overpotential produces anodic Butler-Volmer current");
    require(std::abs(actual_added - expected_added) < 1.0e-12,
            "Butler-Volmer PolyMesh boundary injects Faradaic product species conservatively");
}

void test_fdtd_maxwell3d() {
    cfd::fdtd::Maxwell3D solver({14,13,12,1.0e-3,1.0e-3,1.0e-3,0.8,1.0,1.0});
    solver.initialize_gaussian_ez(0.1,0.14); const double e0=solver.energy();
    solver.add_soft_ez_source(7,6,6,0.01); solver.step(40); const double e1=solver.energy();
    require(e0>0.0 && std::isfinite(e1) && e1>0.0,"3-D Yee FDTD energy remains finite and positive");
    require(std::isfinite(solver.max_field()) && solver.max_field()<10.0,"3-D Yee FDTD remains stable below CFL limit");
    const auto& ez=solver.ez();
    const auto idx=[](std::size_t i,std::size_t j,std::size_t k){return (k*13U+j)*14U+i;};
    require(ez[idx(0,6,6)]==0.0 && ez[idx(13,6,6)]==0.0 && ez[idx(7,0,6)]==0.0 && ez[idx(7,12,6)]==0.0,
            "3-D Yee PEC boundary clamps tangential electric field");
}

void test_fdtd_pmc_ports_and_vtk() {
    cfd::fdtd::Maxwell3D pmc({12,11,10,1.0e-3,1.0e-3,1.0e-3,0.75,1.0,1.0,cfd::fdtd::Boundary3D::pmc});
    pmc.initialize_gaussian_ez(0.1,0.15);
    pmc.step(25U);
    const auto idx=[](std::size_t i,std::size_t j,std::size_t k){return (k*11U+j)*12U+i;};
    for(const auto* h:{&pmc.hx(),&pmc.hy(),&pmc.hz()}){
        require((*h)[idx(0U,5U,5U)]==0.0&&(*h)[idx(11U,5U,5U)]==0.0
                &&(*h)[idx(6U,0U,5U)]==0.0&&(*h)[idx(6U,10U,5U)]==0.0,
                "3-D Yee PMC boundary clamps magnetic field");
    }
    require(std::isfinite(pmc.energy())&&pmc.energy()>0.0,"3-D PMC FDTD remains finite");

    constexpr double frequency=10.0,impedance=50.0,reflection=0.25,transmission=0.8,dt=1.0e-3;
    cfd::fdtd::WavePort1D reference(frequency,impedance,+1);
    cfd::fdtd::WavePort1D transmitted(frequency,impedance,+1);
    for(std::size_t n=0;n<=1000U;++n){
        const double time=dt*static_cast<double>(n);
        const double incident=std::sin(2.0*std::numbers::pi*frequency*time);
        const double reflected=reflection*incident;
        reference.sample(time,incident+reflected,(-incident+reflected)/impedance);
        const double through=transmission*incident;
        transmitted.sample(time,through,-through/impedance);
    }
    const auto s=cfd::fdtd::s_parameters(reference,transmitted);
    require(std::abs(s.s11-std::complex<double>(reflection,0.0))<2.0e-3,"1-D wave port recovers synthetic S11");
    require(std::abs(s.s21-std::complex<double>(transmission,0.0))<2.0e-3,"1-D wave port recovers synthetic S21");

    const auto path=std::filesystem::temp_directory_path()/"cfd_solvers_maxwell3d_test.vtk";
    cfd::fdtd::write_maxwell3d_vtk_ascii(pmc,path);
    std::ifstream in(path);
    const std::string content((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
    require(content.find("DATASET STRUCTURED_POINTS")!=std::string::npos
            &&content.find("VECTORS E double")!=std::string::npos
            &&content.find("VECTORS H double")!=std::string::npos,
            "Maxwell3D legacy VTK export contains geometry and E/H vectors");
    std::filesystem::remove(path);
}


void test_multiphysics_field_registry_and_transfer() {
    cfd::multiphysics::FieldRegistry registry;
    cfd::multiphysics::FieldMetadata temperature;
    temperature.name = "temperature";
    temperature.entities = 4U;
    temperature.units = cfd::multiphysics::kelvin_units();
    temperature.producer = "heat";
    registry.add(temperature, 300.0);
    require(registry.contains("temperature"), "multiphysics field registry stores named field");
    require(registry.data("temperature").size() == 4U, "multiphysics field registry allocates entity data");
    require(registry.metadata("temperature").units.dimensionally_compatible(cfd::multiphysics::kelvin_units()),
            "multiphysics field registry preserves unit signature");
    bool duplicate_rejected = false;
    try { registry.add(temperature); } catch (const std::invalid_argument&) { duplicate_rejected = true; }
    require(duplicate_rejected, "multiphysics field registry rejects duplicate names");

    const cfd::multiphysics::CellGrid1D source{{0.0, 0.25, 0.5, 1.0}};
    const cfd::multiphysics::CellGrid1D target{{0.0, 0.5, 0.75, 1.0}};
    const std::vector<double> source_average{1.0, 2.0, 4.0};
    const auto mapped = cfd::multiphysics::conservative_cell_average_transfer(source, source_average, target);
    require(mapped.size() == 3U && std::abs(mapped[0] - 1.5) < 1.0e-14
            && std::abs(mapped[1] - 4.0) < 1.0e-14 && std::abs(mapped[2] - 4.0) < 1.0e-14,
            "exact-overlap cell transfer reconstructs expected target averages");
    auto integral = [](const cfd::multiphysics::CellGrid1D& grid, std::span<const double> value) {
        double total = 0.0;
        for (std::size_t i = 0U; i < value.size(); ++i) total += value[i] * (grid.edge[i + 1U] - grid.edge[i]);
        return total;
    };
    require(std::abs(integral(source, source_average) - integral(target, mapped)) < 1.0e-14,
            "cell-average transfer conserves integrated scalar exactly");

    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(2, 1, 1, 1.0, 1.0, 1.0);
    std::vector<double> face_flux(mesh.face_count(), 0.0);
    std::size_t internal = mesh.face_count();
    for (std::size_t f = 0U; f < mesh.face_count(); ++f) {
        if (!mesh.faces()[f].boundary()) { internal = f; break; }
    }
    require(internal < mesh.face_count(), "two-cell mesh has an internal face");
    face_flux[internal] = 3.0;
    const auto rate = cfd::multiphysics::conservative_face_flux_to_cell_rate(mesh, face_flux);
    double global = 0.0;
    for (std::size_t c = 0U; c < rate.size(); ++c) global += rate[c] * mesh.cells()[c].volume;
    require(std::abs(global) < 1.0e-14, "face-to-cell transfer cancels internal flux globally");
}

void test_multiphysics_partitioned_aitken() {
    std::vector<double> state{1.0};
    cfd::multiphysics::PartitionedCouplerConfig config;
    config.max_iterations = 100U;
    config.absolute_tolerance = 1.0e-12;
    config.initial_relaxation = 1.0;
    config.aitken = true;
    const auto result = cfd::multiphysics::solve_partitioned_fixed_point(
        state,
        [](std::span<const double> current, std::span<double> candidate) {
            candidate[0] = std::cos(current[0]);
        },
        config);
    require(result.converged, "Aitken partitioned fixed-point iteration converges");
    require(std::abs(state[0] - 0.7390851332151607) < 1.0e-10,
            "Aitken partitioned solver reaches cosine fixed point");
    require(result.iterations < 30U, "Aitken relaxation accelerates representative fixed point");
}


void test_multiphysics_joule_heating() {
    using Type=cfd::fem::ScalarBoundaryType;
    auto mesh=cfd::fem::make_rectangle_tri_mesh(20,10,2.0,1.0);
    cfd::fem::Heat2DConfig thermal;
    thermal.conductivity=1.0;
    thermal.volumetric_heat_capacity=2.0;
    thermal.dt=2.0e-3;
    cfd::multiphysics::JouleHeatingCoupler2D coupled(std::move(mesh),thermal);
    coupled.set_electrical_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
    coupled.set_electrical_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 4.0;},{}});
    coupled.initialize_temperature([](cfd::fem::Node2){return 300.0;});
    coupled.set_thermal_dirichlet([](cfd::fem::Node2,double){return 300.0;});
    coupled.solve_electrical([](cfd::fem::Node2){return 5.0;});
    require(coupled.electrical().linear_result().converged,"electro-thermal DC conduction solve converges");
    double max_q_error=0.0;
    for(double q:coupled.joule_heating_density()) max_q_error=std::max(max_q_error,std::abs(q-20.0));
    require(max_q_error<5.0e-8,"Joule-heating transfer matches sigma*|E|^2 analytical source");
    coupled.thermal_run(5U);
    require(coupled.thermal().linear_result().converged,"electro-thermal heat solve converges");
    double max_t=300.0;
    for(double t:coupled.thermal().temperature()) max_t=std::max(max_t,t);
    require(max_t>300.0&&std::isfinite(max_t),"Joule heating raises interior temperature");
}

void test_fdtd_dispersive_materials() {
    constexpr std::size_t cells = 240U;
    const double dx = 1.0e-4;
    cfd::fdtd::Maxwell1D debye({cells, dx, 0.6, 1.0, 1.0, cfd::fdtd::Boundary1D::mur1});
    const double dt = debye.dt();
    debye.set_debye_material(80U, 180U, 2.0, 3.0, 40.0 * dt, 0.0);
    require(debye.dispersion_model(100U) == cfd::fdtd::DispersionModel1D::debye,
            "Debye material metadata is retained");
    debye.initialize_gaussian(0.25, 0.035);
    debye.step(180U);
    require(std::isfinite(debye.energy()) && debye.energy() > 0.0,
            "Debye ADE FDTD remains finite and positive-energy");
    double pmax = 0.0;
    for (double p : debye.polarization()) pmax = std::max(pmax, std::abs(p));
    require(pmax > 0.0, "Debye material develops polarization response");

    cfd::fdtd::Maxwell1D drude({cells, dx, 0.35, 1.0, 1.0, cfd::fdtd::Boundary1D::mur1});
    const double wscale = 1.0 / drude.dt();
    drude.set_drude_material(80U, 180U, 1.0, 0.08 * wscale, 0.04 * wscale, 0.0);
    require(drude.dispersion_model(100U) == cfd::fdtd::DispersionModel1D::drude,
            "Drude material metadata is retained");
    drude.initialize_gaussian(0.25, 0.035);
    drude.step(120U);
    require(std::isfinite(drude.energy()) && std::isfinite(drude.electric()[120U]),
            "Drude ADE FDTD remains numerically finite");

    cfd::fdtd::Maxwell1D lorentz({cells, dx, 0.35, 1.0, 1.0, cfd::fdtd::Boundary1D::mur1});
    const double wl = 1.0 / lorentz.dt();
    lorentz.set_lorentz_material(80U, 180U, 1.5, 2.0, 0.07 * wl, 0.03 * wl, 0.0);
    require(lorentz.dispersion_model(100U) == cfd::fdtd::DispersionModel1D::lorentz,
            "Lorentz material metadata is retained");
    lorentz.initialize_gaussian(0.25, 0.035);
    lorentz.step(120U);
    require(std::isfinite(lorentz.energy()) && std::isfinite(lorentz.electric()[120U]),
            "Lorentz ADE FDTD remains numerically finite");

    lorentz.set_material(80U, 180U, 2.0, 0.0);
    require(lorentz.dispersion_model(100U) == cfd::fdtd::DispersionModel1D::none,
            "nondispersive material replacement clears ADE state");
}

void test_optics() {
    const auto n = cfd::optics::refract({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, 1.0, 1.5);
    require(n.has_value(), "normal incidence should refract");
    require(std::abs(n->x) < 1e-12 && std::abs(n->y) < 1e-12 && std::abs(n->z - 1.0) < 1e-12,
            "normal incidence direction");
    const auto tir = cfd::optics::refract(cfd::optics::normalized({0.9, 0.0, 0.435889894}),
                                          {0.0, 0.0, -1.0}, 1.5, 1.0);
    require(!tir.has_value(), "total internal reflection detection");

    cfd::optics::MaterialCatalog catalog;
    const double nbk7=catalog.at("N-BK7").refractive_index_nm(587.5618);
    require(std::abs(nbk7-1.5168)<2.0e-4,"Sellmeier N-BK7 d-line refractive index");
    require(catalog.at("N-BK7").refractive_index_nm(486.1327)>nbk7,"Sellmeier normal dispersion decreases index with wavelength");

    const auto fresnel=cfd::optics::fresnel_dielectric(1.0,1.5,1.0);
    require(!fresnel.total_internal_reflection&&std::abs(fresnel.reflectance_s-0.04)<1.0e-13&&std::abs(fresnel.reflectance_p-0.04)<1.0e-13,
            "Fresnel normal-incidence reflectance");
    const double qn=std::sqrt(1.5);
    const double qd=550.0/(4.0*qn);
    require(cfd::optics::thin_film_reflectance_normal(1.0,1.5,550.0,{{qn,qd}})<1.0e-12,
            "quarter-wave single-layer AR coating cancels normal-incidence reflection");

    const double invsqrt2=1.0/std::sqrt(2.0);
    const auto stokes=cfd::optics::stokes_from_jones({{invsqrt2,0.0},{0.0,invsqrt2}});
    require(std::abs(stokes.i-1.0)<1.0e-14&&std::abs(stokes.q)<1.0e-14&&std::abs(stokes.u)<1.0e-14&&std::abs(stokes.v-1.0)<1.0e-14,
            "Jones-to-Stokes conversion identifies circular polarization");

    cfd::optics::SequentialOpticalSystem lens(1.0);
    lens.add_surface({cfd::optics::SurfaceType::sphere,0.0,50.0,12.0,1.5});
    lens.add_surface({cfd::optics::SurfaceType::sphere,5.0,-50.0,12.0,1.0});
    const double bfd=cfd::optics::paraxial_back_focal_distance(lens);
    require(bfd>45.0&&bfd<55.0,"paraxial biconvex-lens back focal distance is physically plausible");
    std::vector<cfd::optics::Ray> rays;
    for(int i=-4;i<=4;++i) rays.push_back({{0.25*static_cast<double>(i),0.0,-10.0},{0.0,0.0,1.0},550.0});
    const auto traced=lens.trace_many(rays);
    for(const auto&r:traced) require(r.valid&&r.surfaces_traced==2U,"sequential real-ray bundle traverses both lens surfaces");
    const double rms=cfd::optics::rms_spot_radius_at_plane(traced,5.0+bfd);
    require(rms<2.0e-3,"real-ray near-axis bundle focuses near paraxial image plane");
}

} // namespace

int main() {
    try {
        test_partition();
        test_cartesian_decomposition();
        test_selective_halo_plan();
        test_virtual_distributed_lbm_parity();
        test_virtual_selective_distributed_parity();
        test_distributed_checkpoint_restart();
        test_lbm_uniform();
        test_lbm_descriptors();
        test_inplace_d2q9_against_pingpong();
        test_generic_pull_parity();
        test_d3q19_uniform();
        test_d3q27_taylor_green_decay();
        test_forced_inplace_parity();
        test_stationary_wall_resets_velocity();
        test_poiseuille_channel();
        test_lid_driven_cavity();
        test_velocity_inlet_pressure_outlet();
        test_taylor_green_convergence();
        test_fp32_fp64_parity();
#if defined(CFD_HAS_SYCL)
        test_sycl_parity();
#endif
        test_core_conjugate_gradient();
        test_sparse_krylov_solvers();
        test_fvm_poly_mesh_operators();
        test_fvm_face_schemes();
        test_fvm_bounded_linear_reconstruction();
        test_fvm_rhie_chow_checkerboard();
        test_fvm_collocated_nonorthogonal_coupling();
        test_fvm_collocated_pressure_channel();
        test_fvm_collocated_lid_cavity();
        test_fvm_pressure_projection();
        test_fvm_transient_taylor_green();
        test_fvm_scalar_transport();
        test_fvm_zero();
        test_fem_poisson();
        test_fem_reference_elements();
        test_fem_matrix_free_laplace();
        test_fem_mixed_boundary_diffusion();
        test_fem_electromagnetics2d();
        test_fem_poisson2d();
        test_fem_poisson3d();
        test_fem_heat2d();
        test_fem_elasticity2d_patch();
        test_fem_axisymmetric_elasticity_patch();
        test_fem_nonlinear_poisson_newton();
        test_fem_adaptive_refinement();
        test_fem_darcy_flow();
        test_fem_magnetostatics();
        test_fem_modal_bar();
        test_fdtd_finite();
        test_fdtd_material_mur_monitors();
        test_fdtd_cpml_tfsf_lumped();
        test_fdtd_maxwell3d();
        test_fdtd_dispersive_materials();
        test_fdtd_pmc_ports_and_vtk();
        test_multiphysics_field_registry_and_transfer();
        test_multiphysics_partitioned_aitken();
        test_multiphysics_joule_heating();
        test_optics();
        test_chemistry_kinetics();
        test_electrochemistry_relations();
        test_mixed_potential_solver();
        test_nernst_planck_transport();
        test_nernst_planck_polymesh();
        test_electroneutral_liquid_junction();
        test_poisson_nernst_planck_potential();
        test_scharfetter_gummel_flux();
        test_butler_volmer_polymesh_boundary_flux();
        std::cout << "all tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "test failure: " << e.what() << '\n';
        return 1;
    }
}
