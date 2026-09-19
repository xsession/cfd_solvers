#include "cfd/circuit/osdi_adapter.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    try {
        cfd::circuit::OsdiDescriptorAdapter adapter(
            2U, [](std::span<const double> voltage, std::span<double> current, std::span<double> jacobian) {
                current[0] = voltage[0] - voltage[1];
                current[1] = -current[0];
                jacobian[0] = 1.0;
                jacobian[1] = -1.0;
                jacobian[2] = -1.0;
                jacobian[3] = 1.0;
            });
        const std::vector<double> voltage{2.0, 0.5};
        std::vector<double> current(2U), jacobian(4U);
        adapter.stamp(voltage, current, jacobian);
        if (current[0] != 1.5 || current[1] != -1.5 || jacobian[1] != -1.0)
            throw std::runtime_error("OSDI stamping adapter result");
        const auto plan = cfd::circuit::make_openvaf_compile_plan("model.va", "model.so");
        if (plan.argv.size() != 4U || plan.argv[0] != "openvaf" || plan.argv[2] != "model.so")
            throw std::runtime_error("OpenVAF compile plan");
        std::cout << "OSDI descriptor/OpenVAF adapter regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "OSDI adapter regression failed: " << error.what() << '\n';
        return 1;
    }
}
