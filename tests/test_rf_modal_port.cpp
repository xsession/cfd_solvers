#include "cfd/em/edge_fem3d.hpp"
#include "cfd/rf/modal_elements.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void test_modal_elements() {
    using namespace cfd::rf;
    const auto series = series_modal_discontinuity({10.0, 2.0});
    require(std::abs(series.a12 - Complex{10.0, 2.0}) < 1.0e-14, "series modal ABCD impedance");
    const auto shunt = shunt_modal_discontinuity({0.01, -0.002});
    require(std::abs(shunt.a21 - Complex{0.01, -0.002}) < 1.0e-14, "shunt modal ABCD admittance");
    const auto s = modal_discontinuity_s({10.0, 0.0}, {0.01, 0.0});
    require(std::isfinite(std::abs(s.a11)) && std::abs(s.a11) <= 1.0, "modal S conversion is finite");
}

void test_edge_thevenin() {
    using namespace cfd::em;
    EdgeMaxwell3DResult field;
    field.edge_voltage_v = {Complex{2.0, 1.0}};
    const auto point = edge_port_thevenin(field, 0U, Complex{0.5, 0.0});
    require(std::abs(point.open_circuit_voltage_v - Complex{2.0, 1.0}) < 1.0e-14, "edge Thevenin voltage");
    require(std::abs(point.port_impedance_ohm - Complex{4.0, 2.0}) < 1.0e-14, "edge Thevenin impedance");
}

} // namespace

int main() {
    try {
        test_modal_elements();
        test_edge_thevenin();
        std::cout << "RF modal/edge-port bridge regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RF modal/edge-port regression failed: " << error.what() << '\n';
        return 1;
    }
}
