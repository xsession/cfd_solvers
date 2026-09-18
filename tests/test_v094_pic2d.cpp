#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
constexpr double eps0 = 8.8541878128e-12;
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

double rms_error(const std::vector<double>& actual, const std::vector<double>& expected) {
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        const double delta = actual[i] - expected[i];
        numerator += delta * delta;
        denominator += expected[i] * expected[i];
    }
    return std::sqrt(numerator / std::max(denominator, 1.0e-300));
}

void two_dimensional_periodic_poisson() {
    using namespace cfd::particle;
    constexpr std::size_t nx = 18U, ny = 14U;
    constexpr double lx = 1.0, ly = 0.7, epsr = 2.5;
    const double kx = 2.0 * std::numbers::pi / lx;
    const double ky = 2.0 * std::numbers::pi / ly;
    std::vector<double> rho(nx * ny), ex_exact(nx * ny), ey_exact(nx * ny);
    for (std::size_t iy = 0; iy < ny; ++iy) {
        const double y = ly * static_cast<double>(iy) / static_cast<double>(ny);
        for (std::size_t ix = 0; ix < nx; ++ix) {
            const double x = lx * static_cast<double>(ix) / static_cast<double>(nx);
            const double phi = std::sin(kx * x) * std::cos(ky * y);
            const std::size_t idx = iy * nx + ix;
            rho[idx] = eps0 * epsr * (kx * kx + ky * ky) * phi;
            ex_exact[idx] = -kx * std::cos(kx * x) * std::cos(ky * y);
            ey_exact[idx] = ky * std::sin(kx * x) * std::sin(ky * y);
        }
    }
    const auto field = periodic_electric_field_from_charge_density_2d(rho, nx, ny, lx, ly, epsr, true);
    require(rms_error(field.electric_x_v_per_m, ex_exact) < 1.0e-12,
            "2-D periodic Poisson Ex matches manufactured sinusoid");
    require(rms_error(field.electric_y_v_per_m, ey_exact) < 1.0e-12,
            "2-D periodic Poisson Ey matches manufactured sinusoid");
}

void charge_conserving_current_2d() {
    using namespace cfd::particle;
    std::vector<PicParticle2D> before(5U), after(5U);
    for (std::size_t i = 0; i < before.size(); ++i) {
        before[i].mass_kg = 1.0;
        before[i].charge_c = (i % 2U == 0U ? 1.0 : -0.7) * 1.0e-12;
        before[i].weight = 1.0;
        before[i].position_m = {0.09 + 0.16 * static_cast<double>(i), 0.11 + 0.13 * static_cast<double>(i)};
        after[i] = before[i];
        after[i].position_m.x += 0.017 + 0.002 * static_cast<double>(i);
        after[i].position_m.y -= 0.011 + 0.003 * static_cast<double>(i);
    }
    const auto current = deposit_charge_conserving_current_2d(before, after, 20U, 16U, 1.0, 0.8, 2.0e-9);
    require(current.current_x_a_per_m2.size() == 20U * 16U && current.current_y_a_per_m2.size() == 20U * 16U,
            "2-D current deposition returns grid-sized vectors");
    require(current.continuity_linf_residual < 1.0e-9, "2-D spectral current reconstruction satisfies continuity");
}

void electrostatic_pic2d_and_gather() {
    using namespace cfd::particle;
    ElectrostaticPic2DConfig config;
    config.nx = 16U;
    config.ny = 12U;
    config.length_x_m = 1.0;
    config.length_y_m = 0.75;
    config.dt_s = 1.0e-5;
    ElectrostaticPic2D solver(config);
    std::vector<PicParticle2D> particles(2U);
    particles[0].mass_kg = 1.0;
    particles[0].charge_c = 1.0e-15;
    particles[0].position_m = {0.23, 0.27};
    particles[1].mass_kg = 1.0;
    particles[1].charge_c = -1.0e-15;
    particles[1].position_m = {0.61, 0.52};
    solver.set_particles(particles);
    solver.deposit_and_solve();
    double peak_field = 0.0;
    for (double value : solver.electric_x())
        peak_field = std::max(peak_field, std::abs(value));
    for (double value : solver.electric_y())
        peak_field = std::max(peak_field, std::abs(value));
    require(peak_field > 0.0, "2-D electrostatic PIC produces nonzero local field for separated neutral pairs");
    const auto field = solver.gather_electric_field({0.33, 0.29});
    require(std::isfinite(field.x) && std::isfinite(field.y), "2-D field gather returns finite bilinear values");
    solver.step(2U);
    require(std::abs(solver.time_s() - 2.0e-5) < 1.0e-15, "2-D electrostatic PIC advances time");
}

void boundary_primitives() {
    using namespace cfd::particle;
    PicParticle2D reflected;
    reflected.mass_kg = 1.0;
    reflected.charge_c = -1.0;
    reflected.position_m = {-0.1, 0.5};
    reflected.velocity_m_per_s = {-2.0, 3.0, 0.0};
    const auto reflection =
        apply_particle_box_boundary_2d(reflected, {{0.0, 0.0}, {1.0, 1.0}}, ParticleWallMode::specular_reflect);
    require(reflection.impacted && reflection.alive && reflected.position_m.x > 0.0 &&
                reflected.velocity_m_per_s.x > 0.0,
            "2-D specular particle wall mirrors position and normal velocity");
    PicParticle2D absorbed;
    absorbed.mass_kg = 9.1093837139e-31;
    absorbed.charge_c = -1.602176634e-19;
    absorbed.weight = 2.0;
    absorbed.position_m = {1.1, 0.2};
    absorbed.velocity_m_per_s = {1.0e7, 0.0, 0.0};
    SecondaryEmissionModel secondary;
    secondary.threshold_energy_ev = 1.0;
    secondary.energy_at_maximum_ev = 20.0;
    secondary.maximum_yield = 1.2;
    const auto absorption =
        apply_particle_box_boundary_2d(absorbed, {{0.0, 0.0}, {1.0, 1.0}}, ParticleWallMode::absorb, &secondary);
    require(absorption.impacted && !absorption.alive && absorption.secondary_macro_weight > 0.0,
            "2-D absorbing wall reports secondary-emission macro-weight");
    std::vector<double> ex(5U * 4U, 2.0), ey(5U * 4U, 3.0);
    GridBoundary2DConfig wall;
    wall.mode = GridBoundaryMode2D::electric_wall;
    apply_electrostatic_field_boundary_2d(ex, ey, 5U, 4U, wall);
    require(ey[0] == 0.0 && ey[4U] == 0.0 && ex[0] == 0.0 && ex[15U] == 0.0,
            "electric-wall primitive zeros tangential boundary fields");
    ex.assign(7U * 7U, 2.0);
    ey.assign(7U * 7U, 3.0);
    GridBoundary2DConfig sponge;
    sponge.mode = GridBoundaryMode2D::absorbing_sponge;
    sponge.sponge_cells = 2U;
    sponge.sponge_strength = 4.0;
    apply_electrostatic_field_boundary_2d(ex, ey, 7U, 7U, sponge);
    require(ex[0] < 2.0 && ey[0] < 3.0 && std::abs(ex[3U + 3U * 7U] - 2.0) < 1.0e-12,
            "absorbing sponge damps boundary layer while retaining interior values");
}
} // namespace

int main() {
    try {
        two_dimensional_periodic_poisson();
        charge_conserving_current_2d();
        electrostatic_pic2d_and_gather();
        boundary_primitives();
        std::cout << "v0.9.4 2-D PIC foundation tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.9.4 PIC2D test failure: " << error.what() << '\n';
        return 1;
    }
}
