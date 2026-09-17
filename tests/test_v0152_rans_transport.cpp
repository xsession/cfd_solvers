#include "cfd/solvers/fvm/rans_transport.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void require_finite_positive(const std::vector<double>& values, double floor, const char* message) {
    for (const double value : values) {
        if (!std::isfinite(value) || value < floor) throw std::runtime_error(message);
    }
}

void test_linear_shear_strain_rate() {
    const auto mesh = cfd::fvm::make_cartesian_hexa_mesh(5, 5, 1, 1.0, 1.0, 1.0);
    std::vector<cfd::fvm::Vec3> velocity(mesh.cell_count());
    constexpr double shear = 3.5;
    for (std::size_t i = 0; i < mesh.cell_count(); ++i) {
        velocity[i] = {shear * mesh.cells()[i].center.y, 0.0, 0.0};
    }
    const auto strain = cfd::fvm::turbulence_strain_rate_magnitude(mesh, velocity);
    for (const double value : strain) {
        require(std::abs(value - shear) < 1.0e-10, "linear shear should give exact strain-rate magnitude");
    }
}

void test_spalart_allmaras_transport() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::SpalartAllmarasConfig cfg;
    cfg.transport.dt = 1.0e-4;
    cfg.transport.temporal_scheme = cfd::fvm::TemporalScheme::backward_bdf2;
    cfd::fvm::SpalartAllmarasTransport model(std::move(mesh), cfg);
    model.initialize(1.0e-5);
    model.set_wall_distance(10.0);
    model.set_strain_rate(50.0);
    const auto initial_nut = model.kinematic_eddy_viscosity();
    model.run(20U);
    require(model.time() > 0.0, "SA model time should advance");
    require_finite_positive(model.nu_tilde(), cfg.minimum_nu_tilde, "SA nu_tilde must remain bounded positive");
    const auto final_nut = model.kinematic_eddy_viscosity();
    require(final_nut.front() > initial_nut.front(), "SA production should raise eddy viscosity in homogeneous shear");
}

void test_kepsilon_decay_and_production() {
    {
        auto mesh = cfd::fvm::make_cartesian_hexa_mesh(6, 1, 1, 1.0, 1.0, 1.0);
        cfd::fvm::KEpsilonConfig cfg;
        cfg.transport.dt = 2.0e-4;
        cfd::fvm::KEpsilonTransport model(std::move(mesh), cfg);
        model.initialize(0.2, 0.05);
        model.set_strain_rate(0.0);
        const double initial_k = model.k().front();
        const double initial_epsilon = model.epsilon().front();
        model.run(20U);
        require(model.k().front() < initial_k, "k-epsilon k should decay without production");
        require(model.epsilon().front() < initial_epsilon, "k-epsilon epsilon should decay without production");
        require_finite_positive(model.k(), cfg.minimum_k, "k-epsilon k must stay positive");
        require_finite_positive(model.epsilon(), cfg.minimum_epsilon, "k-epsilon epsilon must stay positive");
    }
    {
        auto mesh = cfd::fvm::make_cartesian_hexa_mesh(6, 1, 1, 1.0, 1.0, 1.0);
        cfd::fvm::KEpsilonConfig cfg;
        cfg.transport.dt = 1.0e-5;
        cfd::fvm::KEpsilonTransport model(std::move(mesh), cfg);
        model.initialize(0.1, 0.01);
        model.set_strain_rate(20.0);
        const double initial_k = model.k().front();
        model.run(10U);
        require(model.k().front() > initial_k, "k-epsilon production should raise k in homogeneous shear");
        const auto dynamic = model.dynamic_eddy_viscosity();
        require(dynamic.front() > 0.0 && std::isfinite(dynamic.front()), "k-epsilon dynamic eddy viscosity should be finite positive");
    }
}

void test_sst_transport_and_blending() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::KOmegaSSTConfig cfg;
    cfg.transport.dt = 1.0e-5;
    cfg.transport.temporal_scheme = cfd::fvm::TemporalScheme::crank_nicolson;
    cfg.transport.crank_nicolson_off_centering = 0.9;
    cfd::fvm::KOmegaSSTTransport model(std::move(mesh), cfg);
    model.initialize(0.1, 10.0);
    model.set_wall_distance(0.01);
    model.set_strain_rate(20.0);
    const double initial_k = model.k().front();
    model.run(10U);
    require(model.k().front() > initial_k, "SST production should raise k in homogeneous shear");
    require_finite_positive(model.k(), cfg.minimum_k, "SST k must stay positive");
    require_finite_positive(model.omega(), cfg.minimum_omega, "SST omega must stay positive");
    for (const double f1 : model.blending_f1()) require(f1 >= 0.0 && f1 <= 1.0, "SST F1 must be bounded");
    for (const double f2 : model.blending_f2()) require(f2 >= 0.0 && f2 <= 1.0, "SST F2 must be bounded");
    const auto nut = model.kinematic_eddy_viscosity();
    require(nut.front() > 0.0 && std::isfinite(nut.front()), "SST eddy viscosity should be finite positive");
}


void test_sst_des_hybrid_switching() {
    cfd::fvm::KOmegaSSTConfig rans_cfg;
    rans_cfg.transport.dt = 1.0e-5;
    auto mesh_rans = cfd::fvm::make_cartesian_hexa_mesh(6, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::KOmegaSSTTransport rans(std::move(mesh_rans), rans_cfg);
    rans.initialize(0.1, 10.0);
    rans.set_wall_distance(1.0);
    rans.set_strain_rate(0.0);

    auto des_cfg = rans_cfg;
    des_cfg.des_enabled = true;
    des_cfg.des_zonal_filter = cfd::fvm::SSTDESZonalFilter::none;
    auto mesh_des = cfd::fvm::make_cartesian_hexa_mesh(6, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::KOmegaSSTTransport des(std::move(mesh_des), des_cfg);
    des.initialize(0.1, 10.0);
    des.set_wall_distance(1.0);
    des.set_grid_scale(1.0e-3);
    des.set_strain_rate(0.0);

    for (const double factor : rans.hybrid_dissipation_factor()) require(std::abs(factor - 1.0) < 1.0e-14, "RANS SST hybrid factor should be one");
    for (const double factor : des.hybrid_dissipation_factor()) require(factor > 1.0, "DES should increase modeled k dissipation on a fine LES grid");
    rans.step();
    des.step();
    require(des.k().front() < rans.k().front(), "DES length-scale switch should dissipate k faster than pure RANS in the same zero-production state");
}

void test_reynolds_stress_transport_framework() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ReynoldsStressConfig cfg;
    cfg.transport.dt = 1.0e-5;
    cfd::fvm::ReynoldsStressTransport model(std::move(mesh), cfg);
    model.initialize_isotropic(0.15, 0.02);

    std::vector<cfd::fvm::VelocityGradient3> gradient(model.mesh().cell_count());
    for (auto& g : gradient) {
        g.fill(0.0);
        g[1] = 20.0; // du/dy
    }
    model.set_velocity_gradient(gradient);
    model.run(20U);

    const auto stress = model.reynolds_stress();
    const auto k = model.turbulent_kinetic_energy();
    require(stress.front().xy < 0.0, "simple shear should generate the expected negative Rxy covariance");
    require(k.front() > 0.0 && std::isfinite(k.front()), "RSM turbulent kinetic energy should stay finite positive");
    require_finite_positive(model.epsilon(), cfg.minimum_epsilon, "RSM epsilon should remain positive");
    for (const auto& r : stress) require(cfd::fvm::reynolds_stress_is_realizable(r, 1.0e-10), "RSM realizability projection must keep the stress tensor positive semidefinite");

    const auto nut = model.kinematic_eddy_viscosity();
    require(nut.front() > 0.0 && std::isfinite(nut.front()), "RSM diffusion closure eddy viscosity should be finite positive");
}

void test_fixed_boundary_diffusion_path() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(12, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::KEpsilonConfig cfg;
    cfg.transport.dt = 1.0e-3;
    cfg.molecular_viscosity = 1.0e-2;
    cfg.minimum_epsilon = 1.0e-10;
    cfd::fvm::KEpsilonTransport model(std::move(mesh), cfg);
    model.initialize(0.15, 1.0e-10);
    model.set_strain_rate(0.0);
    model.set_k_boundary("left", cfd::fvm::TurbulenceScalarBoundaryType::fixedValue, 0.25);
    model.set_k_boundary("right", cfd::fvm::TurbulenceScalarBoundaryType::fixedValue, 0.05);
    model.run(30U);
    require(model.k().front() > model.k().back(), "fixed turbulence boundaries should drive the expected diffusive gradient");
    const auto [minimum, maximum] = std::minmax_element(model.k().begin(), model.k().end());
    require(*minimum >= 0.05 - 1.0e-10 && *maximum <= 0.25 + 1.0e-10, "bounded implicit diffusion should stay within boundary extrema");
}

} // namespace

int main() {
    try {
        test_linear_shear_strain_rate();
        test_spalart_allmaras_transport();
        test_kepsilon_decay_and_production();
        test_sst_transport_and_blending();
        test_sst_des_hybrid_switching();
        test_reynolds_stress_transport_framework();
        test_fixed_boundary_diffusion_path();
        std::cout << "v0.15.2 RANS transport tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.2 RANS transport test failure: " << error.what() << '\n';
        return 1;
    }
}
