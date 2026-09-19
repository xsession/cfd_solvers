#include "cfd/rf/conductor_models.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
    const auto ground = cfd::rf::finite_ground_reflection(1.0e8, 5.8e7, 1.0, 0.2);
    assert(std::abs(ground.TE) > 0.8 && std::abs(ground.TM) > 0.8);
    const double crowding = cfd::rf::round_wire_proximity_factor(4.0e-3, 5.0e-4, 1.0e9, 5.8e7);
    assert(crowding >= 1.0 && std::isfinite(crowding));
    const cfd::rf::FilamentSegment center{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, 1.0e-4, 5.8e7};
    const auto filaments = cfd::rf::subdivide_rectangular_peec(center, 2.0e-3, 1.0e-3, 2U, 3U);
    assert(filaments.size() == 6U);
    assert(std::abs(filaments.front().start.y + 5.0e-4) < 1.0e-12);
    std::cout << "RF conductor/ground model regression passed\n";
    return 0;
}
