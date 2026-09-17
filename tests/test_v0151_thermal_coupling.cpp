#include "cfd/solvers/fvm/solid_heat.hpp"
#include "cfd/solvers/fvm/thermal_transport.hpp"
#include "cfd/fvm/radiation.hpp"
#include "cfd/fvm/radiation_network.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void test_solid_uniform_source_energy_balance() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(8, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::SolidHeatConfig cfg;
    cfg.dt = 0.01;
    cfg.density = 2.0;
    cfg.heat_capacity = 5.0;
    cfg.conductivity = 0.0;
    cfg.temporal_scheme = cfd::fvm::TemporalScheme::backward_bdf2;
    cfd::fvm::SolidHeatConduction solid(std::move(mesh), cfg);
    solid.initialize(300.0);
    const double initial_energy = solid.total_sensible_energy();
    solid.set_volumetric_heat_source(100.0);
    solid.run(10U);

    require(std::abs(solid.time() - 0.1) < 1.0e-13, "solid time should advance");
    for (const double value : solid.temperature())
        require(std::abs(value - 301.0) < 1.0e-9, "uniform source should produce the exact uniform temperature rise");
    require(std::abs((solid.total_sensible_energy() - initial_energy) - 10.0) < 1.0e-8,
            "solid source energy should equal integrated volumetric power");
}

void test_solid_fixed_temperature_conduction() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(24, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::SolidHeatConfig cfg;
    cfg.dt = 0.02;
    cfg.density = 1.0;
    cfg.heat_capacity = 1.0;
    cfg.conductivity = 0.1;
    cfg.temporal_scheme = cfd::fvm::TemporalScheme::crank_nicolson;
    cfd::fvm::SolidHeatConduction solid(std::move(mesh), cfg);
    solid.set_boundary("left", cfd::fvm::SolidThermalBoundaryType::fixed_temperature, 400.0);
    solid.set_boundary("right", cfd::fvm::SolidThermalBoundaryType::fixed_temperature, 300.0);
    solid.initialize(300.0);
    solid.run(100U);

    const auto& t = solid.temperature();
    require(t.front() > t.back(), "hot-wall conduction should establish a streamwise thermal gradient");
    require(*std::min_element(t.begin(), t.end()) >= 299.999999, "solid conduction should remain lower bounded by the cold wall");
    require(*std::max_element(t.begin(), t.end()) <= 400.000001, "solid conduction should remain upper bounded by the hot wall");
    for (std::size_t i = 1; i < t.size(); ++i)
        require(t[i - 1U] + 1.0e-8 >= t[i], "solid steady trend should be monotone from hot to cold wall");
}

void test_reactive_radiative_energy_source() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(1, 1, 1, 1.0, 1.0, 1.0);
    cfd::fvm::ThermalTransportConfig cfg;
    cfg.dt = 1.0e-4;
    cfg.density = 1.0;
    cfg.default_heat_capacity = 1.0;
    cfg.default_conductivity = 0.0;
    cfd::fvm::ThermalTransport thermal(std::move(mesh), cfg);
    thermal.initialize(500.0);

    auto radiation = std::make_shared<cfd::fvm::OpticallyThinRadiation>(0.2);
    constexpr double chemical_heat_release = 200.0;
    constexpr double environment_temperature = 300.0;
    const double radiative = radiation->volumetric_source(500.0, environment_temperature);
    thermal.set_reactive_radiative_source(
        [](cfd::fvm::Vec3, double, double) { return chemical_heat_release; },
        radiation, environment_temperature);
    thermal.step();

    const double expected = 500.0 + cfg.dt * (chemical_heat_release + radiative);
    require(std::abs(thermal.temperature().front() - expected) < 1.0e-10,
            "thermal transport should add chemistry and radiation in the same energy source contract");
    require(thermal.temperature().front() < 500.0,
            "radiative loss should dominate the selected chemistry source in this regression");
}

void test_boussinesq_helper() {
    auto mesh = cfd::fvm::make_cartesian_hexa_mesh(2, 1, 1, 1.0, 1.0, 1.0);
    const std::vector<double> temperature{300.0, 310.0};
    const auto acceleration = cfd::fvm::boussinesq_acceleration(
        mesh, temperature, 300.0, 3.0e-3, cfd::fvm::Vec3{0.0, -9.81, 0.0});
    require(std::abs(acceleration[0].y) < 1.0e-14, "reference-temperature cell should have no Boussinesq correction");
    require(acceleration[1].y > 0.0, "hotter cell should accelerate opposite gravity under Boussinesq buoyancy");
}
} // namespace

int main() {
    try {
        test_solid_uniform_source_energy_balance();
        test_solid_fixed_temperature_conduction();
        test_reactive_radiative_energy_source();
        test_boussinesq_helper();
        std::cout << "v0.15.1 thermal coupling tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.1 thermal coupling test failed: " << error.what() << '\n';
        return 1;
    }
}
