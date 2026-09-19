#include "cfd/rf/surface_mom.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<cfd::rf::SurfaceTriangle> square_patch() {
    using P = cfd::fem::Point3;
    return {
        {P{0.0, 0.0, 0.0}, P{1.0, 0.0, 0.0}, P{1.0, 1.0, 0.0}},
        {P{0.0, 0.0, 0.0}, P{1.0, 1.0, 0.0}, P{0.0, 1.0, 0.0}},
    };
}

void test_delta_gap_patch() {
    using namespace cfd::rf;
    SurfaceMomConfig config;
    config.frequency_hz = 3.0e8;
    config.excitation = SurfaceMomExcitation::delta_gap;
    config.feed_edge = 2U; // diagonal shared by the two triangles
    config.feed_voltage_v = {1.0, 0.0};
    const auto result = solve_pec_surface_mom(square_patch(), config);
    require(result.basis.size() == 5U, "square patch must produce five RWG/half-RWG functions");
    require(result.basis[2].minus_triangle != std::numeric_limits<std::size_t>::max(),
            "shared patch edge must have two RWG supports");
    require(result.current_a.size() == result.basis.size(), "surface current size");
    require(std::isfinite(result.relative_residual) && result.relative_residual < 1.0e-10,
            "surface MoM direct solve residual");
    require(std::isfinite(std::abs(result.feed_impedance_ohm)) && std::abs(result.feed_current_a) > 1.0e-14,
            "surface MoM feed impedance");
    const auto field = surface_mom_far_field(result, {0.0, 1.0, 1.0}, 10.0);
    require(std::isfinite(field.power_density_w_m2) && field.power_density_w_m2 >= 0.0, "surface MoM far-field power");

    config.feed_voltage_v = {2.0, 0.0};
    const auto doubled = solve_pec_surface_mom(square_patch(), config);
    require(std::abs(doubled.feed_current_a - 2.0 * result.feed_current_a) <
                1.0e-10 * std::max(1.0, std::abs(doubled.feed_current_a)),
            "surface MoM linear feed scaling");
    require(std::abs(doubled.feed_impedance_ohm - result.feed_impedance_ohm) < 1.0e-9,
            "surface MoM feed impedance scaling invariance");
}

void test_plane_wave() {
    using namespace cfd::rf;
    SurfaceMomConfig config;
    config.frequency_hz = 1.0e8;
    config.plane_wave_direction = {0.0, 0.0, 1.0};
    config.plane_wave_electric = {1.0, 0.0, 0.0};
    const auto result = solve_pec_surface_mom(square_patch(), config);
    require(result.feed_edge == std::numeric_limits<std::size_t>::max(), "plane-wave solve must not expose a feed");
    require(result.relative_residual < 1.0e-10, "plane-wave surface MoM residual");
    require(std::isfinite(std::abs(result.current_a.front())), "plane-wave induced current");

    cfd::io::TriangleSurface imported;
    imported.triangles = {
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}},
    };
    const auto imported_result = solve_pec_surface_mom(imported, config);
    require(imported_result.basis.size() == result.basis.size(), "OBJ/STL surface bridge topology");
}

} // namespace

int main() {
    try {
        test_delta_gap_patch();
        test_plane_wave();
        std::cout << "RWG PEC surface MoM regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RWG PEC surface MoM regression failed: " << error.what() << '\n';
        return 1;
    }
}
