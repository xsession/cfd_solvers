#include "cfd/solvers/fvm/thin_film.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void require_close(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance * (1.0 + std::abs(expected)))
        throw std::runtime_error(message);
}

void test_uniform_film_is_stationary() {
    using namespace cfd::fvm;
    ThinFilmConfig config;
    config.dt = 1.0e-4;
    config.surface_tension = 0.05;
    ThinFilmTransport film(make_cartesian_hexa_mesh(12U, 5U, 1U, 1.2, 0.5, 1.0), config);
    film.initialize(0.01);
    const auto before = film.thickness();
    for (int step = 0; step < 5; ++step)
        film.step();
    require(film.face_flux().size() == film.mesh().face_count(), "thin-film face flux size");
    for (std::size_t cell = 0; cell < before.size(); ++cell)
        require_close(film.thickness()[cell], before[cell], 1.0e-13, "uniform film must remain stationary");
    require_close(film.inventory(), 0.01 * 1.2 * 0.5, 1.0e-12, "uniform film inventory");
}

void test_capillary_relaxation_conserves_inventory() {
    using namespace cfd::fvm;
    ThinFilmConfig config;
    config.dt = 2.0e-5;
    config.surface_tension = 0.02;
    config.viscosity = 1.0e-3;
    ThinFilmTransport film(make_cartesian_hexa_mesh(20U, 8U, 1U, 1.0, 0.4, 1.0), config);
    film.initialize([](Vec3 point) {
        const double dx = (point.x - 0.5) / 0.16;
        const double dy = (point.y - 0.2) / 0.10;
        return 0.01 + 0.003 * std::exp(-(dx * dx + dy * dy));
    });
    const double initial_inventory = film.inventory();
    const auto initial_thickness = film.thickness();
    film.step();
    require_close(film.inventory(), initial_inventory, 1.0e-12, "closed capillary film must conserve inventory");
    require(*std::max_element(film.pressure().begin(), film.pressure().end()) -
                    *std::min_element(film.pressure().begin(), film.pressure().end()) >
                0.0,
            "curved film must generate a pressure range");
    bool changed = false;
    for (std::size_t cell = 0; cell < initial_thickness.size(); ++cell)
        changed = changed || std::abs(film.thickness()[cell] - initial_thickness[cell]) > 1.0e-14;
    require(changed, "capillary pressure must drive a thickness update");
    for (const double value : film.thickness())
        require(value >= config.minimum_thickness, "thin-film positivity bound");
}

void test_boundary_pressure_drives_flux() {
    using namespace cfd::fvm;
    ThinFilmConfig config;
    config.dt = 1.0e-6;
    config.surface_tension = 0.0;
    ThinFilmTransport film(make_cartesian_hexa_mesh(8U, 1U, 1U, 1.0, 0.1, 1.0), config);
    film.initialize(0.01);
    std::vector<double> boundary(film.mesh().face_count(), 0.0);
    for (std::size_t face = 0; face < film.mesh().face_count(); ++face)
        if (film.mesh().faces()[face].boundary())
            boundary[face] = 1.0;
    film.set_boundary_pressure(std::move(boundary));
    const double before = film.inventory();
    film.step();
    bool nonzero_flux = false;
    for (const double value : film.face_flux())
        nonzero_flux = nonzero_flux || std::abs(value) > 0.0;
    require(nonzero_flux, "fixed boundary pressure must drive a film flux");
    require(film.inventory() > before, "positive boundary pressure must add film inventory");
}

} // namespace

int main() {
    try {
        test_uniform_film_is_stationary();
        test_capillary_relaxation_conserves_inventory();
        test_boundary_pressure_drives_flux();
        std::cout << "thin-film transport regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "thin-film transport regression failed: " << error.what() << '\n';
        return 1;
    }
}
