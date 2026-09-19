#include "cfd/solvers/fvm/euler_euler.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition)
        throw std::runtime_error(std::string(message));
}

void require_close(double actual, double expected, double tolerance, std::string_view message) {
    if (std::abs(actual - expected) > tolerance)
        throw std::runtime_error(std::string(message));
}

void test_closed_volume_fraction_conservation() {
    using namespace cfd::fvm;
    EulerEulerConfig config;
    config.dt = 1.0e-5;
    config.drag_coefficient = 0.0;
    auto mesh = make_cartesian_hexa_mesh(12U, 4U, 1U, 1.0, 0.4, 1.0);
    EulerEulerTransport phases(std::move(mesh), config);
    phases.initialize(
        [](Vec3 p) {
            return 0.15 + 0.30 * std::exp(-30.0 * ((p.x - 0.35) * (p.x - 0.35) + (p.y - 0.2) * (p.y - 0.2)));
        },
        {0.2, 0.0, 0.0}, {0.45, 0.0, 0.0});
    const double initial = phases.volume_fraction_integral();
    for (std::size_t step = 0; step < 20U; ++step)
        phases.step();
    require_close(phases.volume_fraction_integral(), initial, 1.0e-12,
                  "closed Euler-Euler volume fraction is conservative");
    require(phases.minimum_fraction() >= 0.0 && phases.maximum_fraction() <= 1.0,
            "Euler-Euler volume fraction remains bounded");
}

void test_drag_action_reaction() {
    using namespace cfd::fvm;
    EulerEulerConfig config;
    config.dt = 1.0e-4;
    config.primary_density = 1000.0;
    config.dispersed_density = 4.0;
    config.drag_coefficient = 80.0;
    auto mesh = make_cartesian_hexa_mesh(8U, 2U, 1U, 1.0, 0.25, 1.0);
    EulerEulerTransport phases(std::move(mesh), config);
    phases.initialize(0.25, {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0});
    const double initial_relative = phases.dispersed_velocity().front().x - phases.primary_velocity().front().x;
    const Vec3 initial_momentum = phases.total_momentum();
    for (std::size_t step = 0; step < 40U; ++step)
        phases.step();
    const Vec3 final_momentum = phases.total_momentum();
    require_close(final_momentum.x, initial_momentum.x, 1.0e-10, "Euler-Euler drag conserves total momentum");
    const double final_relative = phases.dispersed_velocity().front().x - phases.primary_velocity().front().x;
    require(final_relative >= 0.0 && final_relative < initial_relative, "Euler-Euler drag reduces phase slip");
    require(phases.integrated_drag_force().x > 0.0, "Euler-Euler drag force acts on primary phase");
}

void test_body_force_and_diagnostics() {
    using namespace cfd::fvm;
    EulerEulerConfig config;
    config.dt = 2.0e-4;
    config.drag_coefficient = 12.0;
    config.dispersed_body_acceleration = {0.0, -1.0, 0.0};
    auto mesh = make_cartesian_hexa_mesh(6U, 3U, 1U, 1.0, 0.5, 1.0);
    EulerEulerTransport phases(std::move(mesh), config);
    phases.initialize(0.2, {0.0, 0.0, 0.0}, {0.0, 0.2, 0.0});
    for (std::size_t step = 0; step < 10U; ++step)
        phases.step();
    require(phases.steps() == 10U && std::abs(phases.time() - 2.0e-3) < 1.0e-15, "Euler-Euler advances physical time");
    for (const auto& velocity : phases.dispersed_velocity())
        require(std::isfinite(velocity.x) && std::isfinite(velocity.y),
                "Euler-Euler body-force velocity remains finite");
}

} // namespace

int main() {
    try {
        test_closed_volume_fraction_conservation();
        test_drag_action_reaction();
        test_body_force_and_diagnostics();
        std::cout << "Euler-Euler phase transport regression passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Euler-Euler regression failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
