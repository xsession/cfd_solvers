#include "cfd/solvers/fvm/compressible_vof.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition)
        throw std::runtime_error(std::string(message));
}

void require_close(double actual, double expected, double tolerance, std::string_view message) {
    if (std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message));
}

void test_uniform_state() {
    using namespace cfd::fvm;
    CompressibleVofConfig config;
    config.dt = 1.0e-4;
    auto mesh = make_cartesian_hexa_mesh(12U, 2U, 1U, 1.0, 0.2, 1.0);
    CompressibleVofTransport vof(std::move(mesh), config);
    vof.initialize(0.35, 1.0e5);
    const auto fraction = vof.liquid_fraction();
    const auto density = vof.density();
    const auto pressure = vof.pressure();
    for (std::size_t step = 0; step < 8U; ++step)
        vof.step();
    for (std::size_t cell = 0; cell < vof.mesh().cell_count(); ++cell) {
        require_close(vof.liquid_fraction()[cell], fraction[cell], 1.0e-14,
                      "zero compressible VOF flux preserves fraction");
        require_close(vof.density()[cell], density[cell], 1.0e-12, "zero compressible VOF flux preserves density");
        require_close(vof.pressure()[cell], pressure[cell], 1.0e-8, "zero compressible VOF flux preserves pressure");
    }
}

void test_closed_mass_and_volume_conservation() {
    using namespace cfd::fvm;
    CompressibleVofConfig config;
    config.dt = 2.0e-5;
    auto mesh = make_cartesian_hexa_mesh(20U, 3U, 1U, 1.0, 0.3, 1.0);
    CompressibleVofTransport vof(std::move(mesh), config);
    vof.initialize([](Vec3 p) {
        return 0.1 + 0.75 * std::exp(-35.0 * ((p.x - 0.35) * (p.x - 0.35) + (p.y - 0.15) * (p.y - 0.15)));
    });
    std::vector<double> flux(vof.mesh().face_count(), 0.0);
    for (std::size_t face = 0; face < vof.mesh().face_count(); ++face) {
        const auto& f = vof.mesh().faces()[face];
        if (!f.boundary())
            flux[face] = 2.0e-4 * (0.5 - f.center.x) * f.area.x;
    }
    vof.set_face_flux(std::move(flux));
    const double initial_volume = vof.liquid_volume();
    const double initial_mass = vof.total_mass();
    for (std::size_t step = 0; step < 30U; ++step)
        vof.step();
    require_close(vof.liquid_volume(), initial_volume, 1.0e-12, "closed compressible VOF conserves liquid volume");
    require_close(vof.total_mass(), initial_mass, 1.0e-10, "closed compressible VOF conserves mixture mass");
    require(vof.minimum_pressure() >= config.minimum_pressure, "compressible VOF pressure remains positive");
    for (double alpha : vof.liquid_fraction())
        require(alpha >= 0.0 && alpha <= 1.0, "compressible VOF fraction remains bounded");
}

void test_barotropic_pressure_response() {
    using namespace cfd::fvm;
    CompressibleVofConfig config;
    config.dt = 1.0e-3;
    config.liquid_bulk_modulus = 1.0e6;
    config.gas_density = 1.0;
    auto mesh = make_cartesian_hexa_mesh(16U, 1U, 1U, 1.0, 1.0, 1.0);
    CompressibleVofTransport vof(std::move(mesh), config);
    vof.initialize([](Vec3 p) { return p.x < 0.45 ? 0.75 : 0.05; });
    std::vector<double> flux(vof.mesh().face_count(), 0.0);
    for (std::size_t face = 0; face < vof.mesh().face_count(); ++face) {
        const auto& f = vof.mesh().faces()[face];
        if (!f.boundary())
            flux[face] = 2.0e-2 * (0.5 - f.center.x) * f.area.x;
    }
    vof.set_face_flux(std::move(flux));
    for (std::size_t step = 0; step < 12U; ++step)
        vof.step();
    require(vof.maximum_pressure() > vof.minimum_pressure() + 1.0,
            "compressible VOF pressure responds to local compression");
    for (double value : vof.density())
        require(std::isfinite(value) && value > 0.0, "compressible VOF density remains physical");
}

} // namespace

int main() {
    try {
        test_uniform_state();
        test_closed_mass_and_volume_conservation();
        test_barotropic_pressure_response();
        std::cout << "compressible VOF regression passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "compressible VOF regression failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
