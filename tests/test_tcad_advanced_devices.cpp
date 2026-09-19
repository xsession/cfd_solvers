#include "cfd/tcad/advanced_devices.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void test_compact_devices() {
    using namespace cfd::tcad;
    GummelPoonBjt bjt;
    const auto bjt_point = bjt.operating_point(0.72, 0.0);
    require(bjt_point.collector_current_a > 0.0 && bjt_point.base_current_a > 0.0, "Gummel-Poon forward operation");
    require(std::abs(bjt_point.collector_current_a + bjt_point.base_current_a + bjt_point.emitter_current_a) < 1.0e-15,
            "Gummel-Poon KCL");

    MosLevel23 mos;
    const auto subthreshold = mos.operating_point(0.2, 1.0);
    const auto strong = mos.operating_point(1.2, 1.0);
    require(subthreshold.subthreshold && strong.drain_current_a > subthreshold.drain_current_a,
            "MOS level 2/3 subthreshold-to-strong inversion transition");

    VdmosPowerDevice vdmos;
    const auto forward = vdmos.operating_point(1.5, 0.5);
    const auto reverse = vdmos.operating_point(0.0, -0.7);
    require(forward.channel_current_a > 0.0 && reverse.body_diode_current_a < 0.0, "VDMOS channel/body-diode branches");
}

void test_hysteretic_core() {
    using namespace cfd::tcad;
    HystereticCoreModel core;
    for (int i = 0; i < 10; ++i)
        core.step(1000.0, 300.0, 1.0e-4);
    const double positive = core.state().induction_t;
    for (int i = 0; i < 20; ++i)
        core.step(-1000.0, 300.0, 1.0e-4);
    require(positive > 0.0 && core.state().induction_t < positive, "hysteretic core reverses induction");
    require(core.state().loss_density_j_per_m3 > 0.0, "hysteretic core accumulates loss");
}

} // namespace

int main() {
    try {
        test_compact_devices();
        test_hysteretic_core();
        std::cout << "advanced TCAD compact-device regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "advanced TCAD device regression failed: " << error.what() << '\n';
        return 1;
    }
}
