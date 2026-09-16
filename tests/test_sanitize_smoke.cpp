#include "cfd/distributed/cartesian.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include "cfd/solvers/lbm/precision_pull.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fvm/schemes.hpp"
#include "cfd/solvers/fvm/projection2d.hpp"
#include "cfd/solvers/fvm/incompressible2d.hpp"
#include "cfd/solvers/fvm/collocated_incompressible.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    cfd::lbm::D2Q9Solver channel({12, 8, 0.8F, 1.0e-6F, 0.0F});
    for (std::size_t x = 0; x < 12; ++x) {
        channel.set_solid(x, 0);
        channel.set_solid(x, 7);
    }
    channel.set_wall_velocity(5, 7, 0.02F, 0.0F);
    channel.set_velocity_inlet_left(0.01F);
    channel.set_pressure_outlet_right(1.0F);
    channel.initialize_uniform(1.0F, 0.01F, 0.0F);
    channel.step(8);
    if (!std::isfinite(channel.max_speed())) return 1;

    cfd::lbm::PrecisionPullSolver<cfd::lbm::D2Q9InPlaceDescriptor, double> fp64({8, 8, 1, 0.75});
    fp64.initialize_taylor_green(0.01);
    fp64.step(4);
    const auto macro = fp64.compute_macroscopic();
    if (macro.rho.size() != 64U || !std::isfinite(static_cast<double>(fp64.mass()))) return 2;

    const cfd::distributed::Extent3 global{9, 7, 5};
    const auto grid = cfd::distributed::choose_process_grid(global, 4, 3);
    std::vector<cfd::lbm::DistributedPullBlock<cfd::lbm::D3Q27Descriptor>> blocks;
    for (std::size_t rank = 0; rank < grid.size(); ++rank) {
        blocks.emplace_back(cfd::lbm::DistributedLbmConfig{global, 0.72F},
                            cfd::distributed::make_brick(global, grid, rank));
        blocks.back().initialize_taylor_green(0.01F);
    }
    cfd::lbm::virtual_selective_step(blocks);
    const auto distributed = cfd::lbm::gather_virtual_macroscopic(blocks);
    if (distributed.rho.size() != global.cells()) return 3;
    for (const float rho : distributed.rho) {
        if (!std::isfinite(rho)) return 4;
    }

    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(4, 3, 2);
    std::vector<double> scalar(mesh.cell_count(), 1.0);
    const auto lap = cfd::fvm::orthogonal_laplacian_scalar(mesh, scalar);
    if (lap.size() != mesh.cell_count()) return 5;

    cfd::fvm::Projection2D projection({12, 10, 1.0, 1.0, 1.0, 1.0e-3, 100});
    projection.initialize_divergent(0.02);
    const double div0 = projection.divergence_l2();
    projection.project();
    if (!std::isfinite(projection.divergence_l2()) || projection.divergence_l2() >= div0) return 6;

    std::vector<cfd::fvm::Vec3> velocity(mesh.cell_count(), {0.1, 0.0, 0.0});
    const auto flux = cfd::fvm::face_flux_from_velocity(mesh, velocity);
    const auto convection = cfd::fvm::convective_divergence_scalar(
        mesh, scalar, flux, cfd::fvm::FaceInterpolationScheme::upwind);
    if (convection.size() != mesh.cell_count()) return 7;

    cfd::fvm::Incompressible2D navier_stokes({12, 12, 1.0, 1.0, 1.0, 1.0e-2, 1.0e-4, 100});
    navier_stokes.initialize_taylor_green(0.02);
    navier_stokes.run(2);
    if (!std::isfinite(navier_stokes.divergence_l2()) || !navier_stokes.pressure_result().converged) return 8;

    auto skew_mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(6, 5, 1, 1.0, 1.0, 1.0, 0.2);
    cfd::fvm::CollocatedIncompressibleConfig collocated_cfg;
    collocated_cfg.dt = 0.02;
    collocated_cfg.kinematic_viscosity = 0.03;
    collocated_cfg.include_convection = false;
    collocated_cfg.momentum_sweeps = 2;
    collocated_cfg.pressure_correctors = 2;
    collocated_cfg.outer_correctors = 2;
    collocated_cfg.nonorthogonal_correctors = 2;
    collocated_cfg.pressure_iterations = 200;
    cfd::fvm::CollocatedIncompressible collocated(std::move(skew_mesh), collocated_cfg);
    for (const auto name : {"left", "right", "bottom", "top"}) {
        collocated.set_velocity_boundary(name, cfd::fvm::VelocityBoundaryType::fixedValue);
    }
    collocated.set_velocity_boundary("front", cfd::fvm::VelocityBoundaryType::slip);
    collocated.set_velocity_boundary("back", cfd::fvm::VelocityBoundaryType::slip);
    collocated.initialize_fields(
        [](cfd::fvm::Vec3 x) { return cfd::fvm::Vec3{0.02 * x.x, -0.01 * x.y, 0.0}; },
        [](cfd::fvm::Vec3) { return 0.0; });
    const auto collocated_result = collocated.step_pimple();
    if (!collocated_result.pressure.converged || !std::isfinite(collocated.continuity_l2())) return 9;

    std::cout << "sanitizer smoke passed\n";
    return 0;
}
