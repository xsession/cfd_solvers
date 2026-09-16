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
#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/solvers/electrochemistry/corrosion1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck1d.hpp"
#include "cfd/solvers/electrochemistry/nernst_planck_poly.hpp"
#include "cfd/solvers/electrochemistry/mixed_potential.hpp"
#include "cfd/solvers/fem/scalar_diffusion2d.hpp"
#include "cfd/solvers/fem/nonlinear_poisson2d.hpp"
#include "cfd/fem/adaptivity.hpp"
#include "cfd/solvers/fem/darcy2d.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include "cfd/solvers/fem/modal_bar1d.hpp"
#include "cfd/solvers/fdtd/maxwell1d.hpp"
#include "cfd/solvers/fdtd/monitor.hpp"
#include "cfd/solvers/optics/sequential.hpp"
#include "cfd/solvers/optics/paraxial.hpp"
#include "cfd/solvers/optics/polarization.hpp"
#include "cfd/multiphysics/field_registry.hpp"
#include "cfd/multiphysics/transfer.hpp"
#include "cfd/multiphysics/partitioned.hpp"
#include "cfd/multiphysics/electro_thermal.hpp"

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

    auto scalar_mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 1, 1);
    cfd::fvm::ScalarTransport scalar_transport(std::move(scalar_mesh), {0.01, 0.02, 100, 10, 1.0e-9});
    scalar_transport.set_boundary("left", cfd::fvm::ScalarBoundaryType::fixedValue, 1.0);
    scalar_transport.set_boundary("right", cfd::fvm::ScalarBoundaryType::fixedValue, 0.0);
    scalar_transport.initialize(0.5);
    scalar_transport.step();
    if (!scalar_transport.linear_result().converged || !std::isfinite(scalar_transport.volume_integral())) return 10;

    cfd::electrochemistry::CorrosionCell1DConfig corrosion_cfg;
    corrosion_cfg.metal_potential = 0.02;
    corrosion_cfg.equilibrium_potential = 0.0;
    corrosion_cfg.exchange_current_density = 1.0;
    const auto corrosion = cfd::electrochemistry::solve_corrosion_cell_1d(corrosion_cfg);
    if (!corrosion.converged || !std::isfinite(corrosion.current_density)) return 11;

    cfd::electrochemistry::NernstPlanck1D transport({32, 1.0, 298.15, 1.0e-4, 0.0,
                                                      cfd::electrochemistry::TransportBoundary1D::periodic});
    transport.set_potential([](double) { return 0.0; });
    const auto species = transport.add_species("ion", 1, 1.0e-3, 1.0);
    const double amount0 = transport.total_amount(species);
    transport.step(4);
    if (!std::isfinite(transport.minimum_concentration(species)) ||
        std::abs(transport.total_amount(species) - amount0) > 1.0e-10) return 12;

    auto electro_mesh = cfd::fvm::make_cartesian_hexa_mesh(12, 1, 1);
    cfd::electrochemistry::NernstPlanckPolyMesh electro_poly(std::move(electro_mesh), {298.15, 1.0e-5});
    electro_poly.add_species("cation", +1, 2.0e-3, [](cfd::fvm::Vec3 x) { return 1.0 + 0.05 * x.x; });
    electro_poly.add_species("anion", -1, 1.0e-3, [](cfd::fvm::Vec3 x) { return 1.0 + 0.05 * x.x; });
    const auto potential_result = electro_poly.solve_electroneutral_potential(500, 1.0e-10);
    if (!potential_result.converged) return 13;
    const auto current_flux = electro_poly.current_flux_faces();
    for (const double current : current_flux) {
        if (!std::isfinite(current)) return 14;
    }

    const auto mixed = cfd::electrochemistry::solve_mixed_potential({
        {"a", 0.0, 1.0, 1.0, 1.0, 0.5, 0.5},
        {"b", 0.1, 1.0, 1.0, 1.0, 0.5, 0.5},
    });
    if (!mixed.converged || !std::isfinite(mixed.potential)) return 15;


    {
        using Type=cfd::fem::ScalarBoundaryType;
        cfd::fem::ScalarDiffusion2D fem(cfd::fem::make_rectangle_tri_mesh(6,5));
        fem.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
        fem.set_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 1.0;},{}});
        fem.solve([](cfd::fem::Node2){return 1.0;},[](cfd::fem::Node2){return 0.0;},[](cfd::fem::Node2){return 0.0;});
        if(!fem.linear_result().converged) return 16;
    }
    {
        const double pi=3.14159265358979323846;
        cfd::fem::NonlinearPoisson2D nonlinear(cfd::fem::make_rectangle_tri_mesh(6,6),{1.0,1.0});
        nonlinear.solve([&](cfd::fem::Node2 p){const double u=std::sin(pi*p.x)*std::sin(pi*p.y);return 2.0*pi*pi*u+u*u*u;},
                        [](cfd::fem::Node2){return 0.0;});
        if(!nonlinear.nonlinear_result().converged) return 21;
        const auto estimate=cfd::fem::estimate_poisson_error_tri3(
            nonlinear.mesh(),nonlinear.solution(),[&](cfd::fem::Node2 p){const double u=std::sin(pi*p.x)*std::sin(pi*p.y);return 2.0*pi*pi*u;});
        const auto marked=cfd::fem::mark_dorfler(estimate.element_indicator,0.5);
        const auto refined=cfd::fem::refine_tri3_longest_edges(nonlinear.mesh(),marked);
        if(refined.element_count()<=nonlinear.mesh().element_count()) return 22;
    }
    {
        using Type=cfd::fem::ScalarBoundaryType;
        cfd::fem::Darcy2D darcy(cfd::fem::make_rectangle_tri_mesh(6,3,2.0,1.0));
        darcy.set_boundary(0,{Type::dirichlet,[](cfd::fem::Node2){return 1.0;},{}});
        darcy.set_boundary(1,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
        darcy.solve([](cfd::fem::Node2){return 1.0e-12;});
        if(darcy.element_velocity([](cfd::fem::Node2){return 1.0e-12;}).empty()) return 23;
        cfd::fem::Magnetostatics2D mag(cfd::fem::make_rectangle_tri_mesh(6,6));
        for(int patch=0;patch<4;++patch)mag.set_boundary(patch,{Type::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
        mag.solve([](cfd::fem::Node2){return 1.0;},[](cfd::fem::Node2 p){return std::sin(3.14159265358979323846*p.x)*std::sin(3.14159265358979323846*p.y);});
        if(!mag.linear_result().converged) return 24;
        cfd::fem::BarModal1DConfig modal_cfg;modal_cfg.elements=12U;modal_cfg.eigen.relative_tolerance=1.0e-6;
        if(cfd::fem::BarModal1D(modal_cfg).solve(1U).empty()) return 25;
    }
    {
        cfd::multiphysics::FieldRegistry registry;
        cfd::multiphysics::FieldMetadata meta; meta.name="T"; meta.entities=2U; meta.units=cfd::multiphysics::kelvin_units();
        registry.add(meta,300.0);
        const auto mapped=cfd::multiphysics::conservative_cell_average_transfer({{0.0,0.5,1.0}},{registry.data("T").data(),2U},{{0.0,0.25,0.75,1.0}});
        if(mapped.size()!=3U) return 26;
        std::vector<double> state{1.0};
        const auto result=cfd::multiphysics::solve_partitioned_fixed_point(state,[](std::span<const double> x,std::span<double> y){y[0]=std::cos(x[0]);});
        if(!result.converged) return 27;
    }
    {
        using BType=cfd::fem::ScalarBoundaryType;
        cfd::fem::Heat2DConfig th; th.dt=1.0e-3;
        cfd::multiphysics::JouleHeatingCoupler2D et(cfd::fem::make_rectangle_tri_mesh(4,3),th);
        et.set_electrical_boundary(0,{BType::dirichlet,[](cfd::fem::Node2){return 0.0;},{}});
        et.set_electrical_boundary(1,{BType::dirichlet,[](cfd::fem::Node2){return 1.0;},{}});
        et.initialize_temperature([](cfd::fem::Node2){return 300.0;});
        et.set_thermal_dirichlet([](cfd::fem::Node2,double){return 300.0;});
        et.solve_electrical([](cfd::fem::Node2){return 2.0;});
        et.thermal_step();
        if(!et.thermal().linear_result().converged) return 29;
    }
    {
        cfd::fdtd::Maxwell1D dispersive({96,1.0e-4,0.35,1.0,1.0,cfd::fdtd::Boundary1D::mur1});
        const double w=1.0/dispersive.dt();
        dispersive.set_lorentz_material(32,72,1.5,1.0,0.05*w,0.02*w);
        dispersive.initialize_gaussian(0.2,0.04); dispersive.step(20);
        if(!std::isfinite(dispersive.energy())) return 28;
    }
    {
        cfd::fdtd::Maxwell1D em({96,1.0e-3,0.9,1.0,1.0,cfd::fdtd::Boundary1D::mur1});
        em.set_material(48,80,2.5,0.01); em.initialize_gaussian(0.3,0.05); em.step(20);
        if(!std::isfinite(em.energy())) return 17;
        cfd::fdtd::DftMonitor dft(1.0); dft.sample(0.0,0.0); dft.sample(0.25,1.0); dft.sample(0.5,0.0);
        if(!std::isfinite(std::abs(dft.integral()))) return 18;
    }
    {
        cfd::optics::SequentialOpticalSystem lens;
        lens.add_surface({cfd::optics::SurfaceType::sphere,0.0,50.0,10.0,1.5});
        lens.add_surface({cfd::optics::SurfaceType::sphere,5.0,-50.0,10.0,1.0});
        const auto ray=lens.trace({{0.5,0.0,-5.0},{0.0,0.0,1.0},550.0});
        if(!ray.valid||!std::isfinite(cfd::optics::paraxial_back_focal_distance(lens))) return 19;
        const auto fr=cfd::optics::fresnel_dielectric(1.0,1.5,1.0);
        if(!std::isfinite(fr.reflectance_s)) return 20;
    }

    std::cout << "sanitizer smoke passed\n";
    return 0;
}
