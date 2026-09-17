#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include "cfd/solvers/lbm/advanced_d2q9.hpp"
#include "cfd/solvers/lbm/visualization.hpp"
#include "cfd/solvers/lbm/particles.hpp"
#include "cfd/solvers/lbm/compressed_pull.hpp"
#include "cfd/solvers/lbm/free_surface.hpp"
#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/voxelize_sycl.hpp"
#endif

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}

void test_smagorinsky_relaxation() {
    using D = cfd::lbm::D3Q19Descriptor;
    std::array<float, D::q> populations{};
    for (int d = 0; d < D::q; ++d) {
        populations[static_cast<std::size_t>(d)] =
            cfd::lbm::detail::equilibrium<D>(d, 1.0F, 0.0F, 0.0F, 0.0F);
    }
    const float base_tau = 0.56F;
    const float equilibrium_tau = cfd::lbm::detail::smagorinsky_relaxation_time<D>(
        populations, 1.0F, 0.0F, 0.0F, 0.0F, base_tau);
    require(std::abs(equilibrium_tau - base_tau) < 1.0e-7F,
            "Smagorinsky equilibrium state should retain base relaxation time");

    // A mass- and momentum-neutral opposite-direction perturbation creates a
    // non-equilibrium normal stress and must increase the local relaxation time.
    constexpr float epsilon = 0.01F;
    populations[0] -= 2.0F * epsilon;
    populations[1] += epsilon;
    populations[2] += epsilon;
    const float turbulent_tau = cfd::lbm::detail::smagorinsky_relaxation_time<D>(
        populations, 1.0F, 0.0F, 0.0F, 0.0F, base_tau);
    require(turbulent_tau > base_tau, "Smagorinsky closure should add local eddy viscosity");

    cfd::lbm::InPlaceLbmConfig cfg{20, 18, 16, 0.56F};
    cfg.smagorinsky_les = true;
    cfd::lbm::D3Q19Solver solver(cfg);
    solver.initialize_taylor_green(0.05F);
    const double mass0 = solver.mass();
    const double energy0 = solver.kinetic_energy();
    solver.step(24);
    const double mass1 = solver.mass();
    const double energy1 = solver.kinetic_energy();
    require(std::abs(mass1 - mass0) / mass0 < 5.0e-6,
            "Smagorinsky D3Q19 should conserve mass");
    require(std::isfinite(energy1) && energy1 > 0.0 && energy1 < energy0,
            "Smagorinsky D3Q19 Taylor-Green energy should decay");
}


void test_bouzidi_curved_wall() {
    using cfd::lbm::bouzidi_interpolated_bounce_back;
    require(std::abs(bouzidi_interpolated_bounce_back(0.5F, 2.0F, 7.0F, 11.0F) - 2.0F) < 1.0e-7F,
            "Bouzidi q=0.5 must reduce to halfway bounce-back");
    require(std::abs(bouzidi_interpolated_bounce_back(0.25F, 2.0F, 4.0F, 11.0F) - 3.0F) < 1.0e-7F,
            "Bouzidi near-wall interpolation branch");
    require(std::abs(bouzidi_interpolated_bounce_back(0.75F, 2.0F, 4.0F, 8.0F) - 4.0F) < 1.0e-6F,
            "Bouzidi far-wall interpolation branch");

    constexpr std::size_t nx = 16, ny = 12;
    cfd::lbm::D2Q9Solver halfway({nx, ny, 0.72F});
    cfd::lbm::D2Q9Solver interpolated({nx, ny, 0.72F});
    for (std::size_t y = 0; y < ny; ++y) {
        halfway.set_solid(0, y);
        interpolated.set_solid(0, y);
        for (int d : {1, 5, 8}) interpolated.set_interpolated_wall_link(1, y, d, 0.5F);
        for (int d : {3, 6, 7}) interpolated.set_interpolated_wall_link(nx - 1, y, d, 0.5F);
    }
    halfway.initialize_uniform(1.0F, 0.02F, 0.0F);
    interpolated.initialize_uniform(1.0F, 0.02F, 0.0F);
    halfway.step(4);
    interpolated.step(4);
    const auto& a = halfway.velocity_x();
    const auto& b = interpolated.velocity_x();
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) max_error = std::max(max_error, std::abs(static_cast<double>(a[i] - b[i])));
    require(max_error < 2.0e-6, "Bouzidi half-link configuration should match halfway bounce-back");
}


void test_immersed_boundary_particles() {
    constexpr std::size_t nx = 8, ny = 8, nz = 8;
    const std::size_t cells = nx * ny * nz;
    std::vector<float> ux(cells, 0.1F), uy(cells, -0.02F), uz(cells, 0.03F);
    cfd::lbm::ParticleVector3 p{2.25, 3.5, 4.75};
    const auto sampled = cfd::lbm::interpolate_particle_velocity(nx, ny, nz, ux, uy, uz, p);
    require(std::abs(sampled.x - 0.1) < 1.0e-7 && std::abs(sampled.y + 0.02) < 1.0e-7 && std::abs(sampled.z - 0.03) < 1.0e-7,
            "immersed-boundary interpolation should reproduce a uniform velocity field");

    cfd::lbm::ImmersedBoundaryParticle passive{p, {}, 2.0, 1.0};
    auto passive_result = cfd::lbm::advance_immersed_boundary_particles(nx, ny, nz, ux, uy, uz,
        std::span<cfd::lbm::ImmersedBoundaryParticle>(&passive, 1), 0.1, false);
    require(passive.velocity.x > 0.0 && passive_result.acceleration_x.empty(),
            "passive immersed-boundary particle should follow flow without fluid reaction");

    cfd::lbm::ImmersedBoundaryParticle coupled{p, {}, 2.0, 1.0};
    auto coupled_result = cfd::lbm::advance_immersed_boundary_particles(nx, ny, nz, ux, uy, uz,
        std::span<cfd::lbm::ImmersedBoundaryParticle>(&coupled, 1), 0.1, true);
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    for (std::size_t i = 0; i < cells; ++i) {
        sum_x += coupled_result.acceleration_x[i];
        sum_y += coupled_result.acceleration_y[i];
        sum_z += coupled_result.acceleration_z[i];
    }
    require(std::abs(sum_x - coupled_result.total_reaction_force.x) < 1.0e-7 &&
            std::abs(sum_y - coupled_result.total_reaction_force.y) < 1.0e-7 &&
            std::abs(sum_z - coupled_result.total_reaction_force.z) < 1.0e-7,
            "two-way immersed-boundary spreading should conserve the equal-and-opposite reaction force");
    require(coupled_result.total_reaction_force.x < 0.0,
            "two-way particle reaction should oppose positive fluid-to-particle drag");

    cfd::lbm::D3Q19Solver solver({nx, ny, nz, 0.7F});
    solver.initialize_uniform();
    solver.set_local_body_acceleration(coupled_result.acceleration_x,
                                       coupled_result.acceleration_y,
                                       coupled_result.acceleration_z);
    require(solver.has_local_body_acceleration(), "LBM solver should accept per-cell particle reaction acceleration");
    solver.step();
    const auto macro = solver.compute_macroscopic();
    double mean_ux = 0.0;
    for (float value : macro.ux) mean_ux += value;
    mean_ux /= static_cast<double>(macro.ux.size());
    require(mean_ux < 0.0, "two-way particle reaction should feed back into the LBM velocity field");
    solver.clear_local_body_acceleration();
    require(!solver.has_local_body_acceleration(), "local LBM force field should be clearable");
}


template<class Codec>
void test_compressed_codec_solver(const char* label, double mass_gate, double field_gate) {
    using D = cfd::lbm::D2Q9InPlaceDescriptor;
    using Solver = cfd::lbm::CompressedPopulationPullSolver<D, Codec>;
    Solver compressed({28, 24, 1, 0.72F});
    cfd::lbm::D2Q9InPlaceSolver reference({28, 24, 1, 0.72F});
    compressed.initialize_taylor_green(0.025F);
    reference.initialize_taylor_green(0.025F);
    const double mass0 = compressed.mass();
    compressed.step(16);
    reference.step(16);
    const double mass_error = std::abs(compressed.mass() - mass0) / mass0;
    require(mass_error < mass_gate, std::string(label) + " compressed-population mass gate");
    const auto a = compressed.compute_macroscopic();
    const auto b = reference.compute_macroscopic();
    double max_error = 0.0;
    for (std::size_t i = 0; i < a.rho.size(); ++i) {
        max_error = std::max(max_error, std::abs(static_cast<double>(a.rho[i] - b.rho[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.ux[i] - b.ux[i])));
        max_error = std::max(max_error, std::abs(static_cast<double>(a.uy[i] - b.uy[i])));
    }
    require(max_error < field_gate, std::string(label) + " compressed-population macroscopic error gate");
    require(compressed.population_bytes() == 28U * 24U * static_cast<std::size_t>(D::q) * 2U,
            std::string(label) + " compressed storage must use 16 bits/population");
}

void test_population_compression() {
    const float probe = 0.123456F;
    const float half = cfd::lbm::IeeeHalfCodec::decode(cfd::lbm::IeeeHalfCodec::encode(probe));
    const float bf16 = cfd::lbm::BFloat16Codec::decode(cfd::lbm::BFloat16Codec::encode(probe));
    require(std::abs(half - probe) < 1.0e-4F, "IEEE FP16 codec round-trip");
    require(std::abs(bf16 - probe) < 5.0e-4F, "BF16 codec round-trip");
    test_compressed_codec_solver<cfd::lbm::ShiftedHalfCodec<cfd::lbm::D2Q9InPlaceDescriptor>>(
        "shifted FP16", 2.0e-4, 4.0e-4);
    test_compressed_codec_solver<cfd::lbm::BFloat16Codec>("BF16", 2.0e-3, 4.0e-3);
    test_compressed_codec_solver<cfd::lbm::PackedDeviation16Codec<cfd::lbm::D2Q9InPlaceDescriptor>>(
        "custom packed deviation", 8.0e-4, 1.5e-3);
}


#if defined(CFD_HAS_SYCL)
void test_sycl_voxelization() {
    using P = cfd::fem::Point3;
    cfd::io::TriangleSurface cube;
    const P p000{.25,.25,.25},p100{.75,.25,.25},p010{.25,.75,.25},p110{.75,.75,.25};
    const P p001{.25,.25,.75},p101{.75,.25,.75},p011{.25,.75,.75},p111{.75,.75,.75};
    auto tri=[&](P a,P b,P c){cube.triangles.push_back({a,b,c});};
    tri(p000,p110,p100);tri(p000,p010,p110);tri(p001,p101,p111);tri(p001,p111,p011);
    tri(p000,p100,p101);tri(p000,p101,p001);tri(p010,p011,p111);tri(p010,p111,p110);
    tri(p000,p001,p011);tri(p000,p011,p010);tri(p100,p110,p111);tri(p100,p111,p101);
    const cfd::lbm::VoxelGrid3D grid{8,8,8,{0,0,0},{1,1,1}};
    const auto cpu=cfd::lbm::voxelize_surface(cube,grid);
    sycl::queue queue{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
    const auto gpu=cfd::lbm::voxelize_surface_sycl(cube,grid,queue);
    require(gpu==cpu,"SYCL GPU voxelization should match CPU parity mask");
}
#endif


void test_lightweight_html_viewer() {
    constexpr std::size_t nx=6,ny=5;
    std::vector<double> field(nx*ny);
    for(std::size_t y=0;y<ny;++y)for(std::size_t x=0;x<nx;++x)field[y*nx+x]=static_cast<double>(x)+0.25*static_cast<double>(y);
    const auto path=std::filesystem::temp_directory_path()/"cfd_solvers_lbm_viewer_v0141.html";
    cfd::lbm::write_html_scalar_viewer(path,nx,ny,field);
    std::ifstream in(path);const std::string html((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
    require(html.find("wheel: zoom")!=std::string::npos&&html.find("pointermove")!=std::string::npos&&html.find("const nx=6,ny=5")!=std::string::npos,
            "standalone LBM HTML viewer should embed data and interaction controls");
    std::filesystem::remove(path);
}


void test_free_surface_vof() {
    constexpr std::size_t nx=12,ny=8,nz=6,cells=nx*ny*nz;
    std::vector<double> phi(cells,0.0);
    for(std::size_t z=0;z<nz;++z)for(std::size_t y=0;y<ny;++y)for(std::size_t x=0;x<nx;++x){
        phi[(z*ny+y)*nx+x]=x<4?1.0:(x==4?0.5:0.0);
    }
    auto surface=cfd::lbm::make_free_surface_field(nx,ny,nz,phi);
    require(surface.phase[(2*ny+2)*nx+2]==cfd::lbm::FreeSurfacePhase::fluid&&
            surface.phase[(2*ny+2)*nx+4]==cfd::lbm::FreeSurfacePhase::interface_cell&&
            surface.phase[(2*ny+2)*nx+8]==cfd::lbm::FreeSurfacePhase::gas,
            "free-surface phase classification");
    const double volume0=cfd::lbm::free_surface_volume(surface);
    std::vector<float> ux(cells,0.1F),uy(cells,0.0F),uz(cells,0.0F);
    cfd::lbm::advect_free_surface_vof(surface,ux,uy,uz,0.25);
    require(std::abs(cfd::lbm::free_surface_volume(surface)-volume0)<1.0e-10,
            "periodic VOF advection should conserve fill volume");

    std::vector<double> ramp(cells);
    for(std::size_t z=0;z<nz;++z)for(std::size_t y=0;y<ny;++y)for(std::size_t x=0;x<nx;++x)
        ramp[(z*ny+y)*nx+x]=static_cast<double>(x)/static_cast<double>(nx-1);
    auto plane=cfd::lbm::make_free_surface_field(nx,ny,nz,ramp);
    const auto normals=cfd::lbm::free_surface_normals(plane);
    const auto curvature=cfd::lbm::free_surface_curvature(plane);
    const auto center=(3U*ny+4U)*nx+5U;
    require(normals[center].x>0.999&&std::abs(normals[center].y)<1e-12&&std::abs(normals[center].z)<1e-12,
            "planar VOF interface normal");
    require(std::abs(curvature[center])<1e-12,"planar VOF interface curvature");
    const auto capillary=cfd::lbm::free_surface_surface_tension_acceleration(plane,0.072,1.0);
    require(std::abs(capillary[center].x)<1e-12,"planar interface should have zero capillary force");
}


void test_cumulant_collision() {
    cfd::lbm::CumulantD2Q9Solver solver({48,40,0.58F});
    solver.initialize_taylor_green(0.035F);
    const double mass0=solver.mass(),energy0=solver.kinetic_energy();
    solver.step(60);
    const double mass1=solver.mass(),energy1=solver.kinetic_energy();
    require(std::abs(mass1-mass0)/mass0<5.0e-6,"D2Q9 cumulant collision mass conservation");
    require(std::isfinite(energy1)&&energy1>0.0&&energy1<energy0,"D2Q9 cumulant collision Taylor-Green decay");
}

void test_q_criterion() {
    constexpr std::size_t nx = 7, ny = 7, nz = 7;
    std::vector<cfd::lbm::Vector3> velocity(nx * ny * nz);
    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                const double xf = static_cast<double>(x) - 3.0;
                const double yf = static_cast<double>(y) - 3.0;
                velocity[(z * ny + y) * nx + x] = {-yf, xf, 0.0};
            }
        }
    }
    const auto q = cfd::lbm::q_criterion_3d(nx, ny, nz, velocity);
    const auto center = (3U * ny + 3U) * nx + 3U;
    require(std::abs(q[center] - 1.0) < 1.0e-12,
            "solid-body rotation should have Q=1 for unit angular speed");

    for (std::size_t z = 0; z < nz; ++z) {
        for (std::size_t y = 0; y < ny; ++y) {
            for (std::size_t x = 0; x < nx; ++x) {
                velocity[(z * ny + y) * nx + x] = {static_cast<double>(y), 0.0, 0.0};
            }
        }
    }
    const auto shear_q = cfd::lbm::q_criterion_3d(nx, ny, nz, velocity);
    require(std::abs(shear_q[center]) < 1.0e-12,
            "simple shear should balance strain and rotation in Q-criterion");
}
} // namespace

int main() {
    try {
        test_smagorinsky_relaxation();
        test_bouzidi_curved_wall();
        test_immersed_boundary_particles();
        test_population_compression();
#if defined(CFD_HAS_SYCL)
        test_sycl_voxelization();
#endif
        test_q_criterion();
        test_lightweight_html_viewer();
        test_free_surface_vof();
        test_cumulant_collision();
        std::cout << "v0.14.1 LBM LES/Q-criterion tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "v0.14.1 LBM LES/Q-criterion test failed: " << e.what() << '\n';
        return 1;
    }
}
