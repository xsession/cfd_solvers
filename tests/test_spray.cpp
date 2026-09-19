#include "cfd/solvers/fvm/spray.hpp"

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

void test_scheduled_injection() {
    using namespace cfd::fvm;
    SprayConfig config;
    config.dt = 1.0e-4;
    config.evaporation_constant = 0.0;
    auto mesh = make_cartesian_hexa_mesh(8U, 2U, 1U, 1.0, 0.25, 1.0);
    SprayInjectionEvaporation spray(std::move(mesh), config);
    spray.set_injection({{0.2, 0.1, 0.5}, {0.4, 0.0, 0.0}, 1.0e-4, 300.0, 2.0e-6, 0.0, 3.0e-4});
    const std::vector<Vec3> velocity(spray.mesh().cell_count());
    for (std::size_t step = 0; step < 5U; ++step)
        spray.step(velocity, 300.0);
    require(spray.steps() == 5U && spray.parcels().size() == 3U, "spray injection follows its active time window");
    require_close(spray.injected_mass(), 6.0e-10, 1.0e-20, "spray injection mass flow is integrated by timestep");
}

void test_evaporation_mass_and_energy_sources() {
    using namespace cfd::fvm;
    SprayConfig config;
    config.dt = 1.0e-4;
    config.evaporation_constant = 5.0e-6;
    config.latent_heat = 2.0e6;
    auto mesh = make_cartesian_hexa_mesh(6U, 2U, 1U, 1.0, 0.25, 1.0);
    SprayInjectionEvaporation spray(std::move(mesh), config);
    spray.add_parcel({{0.4, 0.1, 0.5}, {}, 1.0e-3, 1.0e-6, 300.0});
    const double initial_mass = spray.total_liquid_mass();
    const std::vector<Vec3> velocity(spray.mesh().cell_count());
    for (std::size_t step = 0; step < 8U; ++step)
        spray.step(velocity, 450.0);
    require(spray.evaporated_mass() > 0.0 && spray.total_liquid_mass() < initial_mass,
            "spray D2 evaporation decreases liquid mass");
    require_close(spray.evaporated_mass() + spray.total_liquid_mass(), initial_mass, 1.0e-18,
                  "spray liquid mass closes against evaporation inventory");
    double vapor_rate = 0.0;
    double energy_rate = 0.0;
    for (std::size_t cell = 0; cell < spray.mesh().cell_count(); ++cell) {
        vapor_rate += spray.vapor_mass_source()[cell] * spray.mesh().cells()[cell].volume;
        energy_rate += spray.carrier_energy_source()[cell] * spray.mesh().cells()[cell].volume;
    }
    require(vapor_rate > 0.0 && energy_rate < 0.0, "spray evaporation emits vapor and removes latent heat");
}

void test_drag_reaction_source() {
    using namespace cfd::fvm;
    SprayConfig config;
    config.dt = 1.0e-4;
    config.evaporation_constant = 0.0;
    auto mesh = make_cartesian_hexa_mesh(4U, 1U, 1U, 1.0, 1.0, 1.0);
    SprayInjectionEvaporation spray(std::move(mesh), config);
    spray.add_parcel({{0.4, 0.5, 0.5}, {2.0, 0.0, 0.0}, 1.0e-4, 1.0e-6, 300.0});
    const double initial_velocity = spray.parcels().front().velocity.x;
    const std::vector<Vec3> carrier(spray.mesh().cell_count());
    spray.step(carrier, 300.0);
    require(spray.parcels().front().velocity.x < initial_velocity, "spray Stokes drag reduces parcel slip");
    double carrier_impulse = 0.0;
    for (std::size_t cell = 0; cell < spray.mesh().cell_count(); ++cell)
        carrier_impulse += spray.carrier_momentum_source()[cell].x * spray.mesh().cells()[cell].volume * config.dt;
    const double parcel_impulse =
        spray.parcels().front().mass * (spray.parcels().front().velocity.x - initial_velocity);
    require_close(carrier_impulse + parcel_impulse, 0.0, 1.0e-15,
                  "spray drag source is equal and opposite to parcel momentum change");
}

} // namespace

int main() {
    try {
        test_scheduled_injection();
        test_evaporation_mass_and_energy_sources();
        test_drag_reaction_source();
        std::cout << "spray injection/evaporation regression passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "spray regression failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
