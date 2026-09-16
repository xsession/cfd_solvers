#include "cfd/core/conjugate_gradient.hpp"
#include "cfd/core/decomposition.hpp"
#include "cfd/distributed/cartesian.hpp"
#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fem/poisson1d.hpp"
#include "cfd/solvers/fvm/diffusion2d.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"
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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
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

void test_fdtd_finite() {
    cfd::fdtd::Maxwell1D solver({256, 1.0e-3, 0.9, 1.0, 1.0});
    solver.initialize_gaussian();
    solver.step(80);
    const double energy = solver.energy();
    require(std::isfinite(energy) && energy > 0.0, "FDTD stable finite energy");
}

void test_optics() {
    const auto n = cfd::optics::refract({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, 1.0, 1.5);
    require(n.has_value(), "normal incidence should refract");
    require(std::abs(n->x) < 1e-12 && std::abs(n->y) < 1e-12 && std::abs(n->z - 1.0) < 1e-12,
            "normal incidence direction");
    const auto tir = cfd::optics::refract(cfd::optics::normalized({0.9, 0.0, 0.435889894}),
                                          {0.0, 0.0, -1.0}, 1.5, 1.0);
    require(!tir.has_value(), "total internal reflection detection");
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
        test_fvm_poly_mesh_operators();
        test_fvm_face_schemes();
        test_fvm_rhie_chow_checkerboard();
        test_fvm_collocated_nonorthogonal_coupling();
        test_fvm_collocated_pressure_channel();
        test_fvm_collocated_lid_cavity();
        test_fvm_pressure_projection();
        test_fvm_transient_taylor_green();
        test_fvm_zero();
        test_fem_poisson();
        test_fdtd_finite();
        test_optics();
        std::cout << "all tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "test failure: " << e.what() << '\n';
        return 1;
    }
}
