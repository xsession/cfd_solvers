#include "common/benchmark.hpp"
#include "cfd/solvers/fvm/flame_models.hpp"
#include "cfd/solvers/fvm/soot_radiation.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::fvm;

    const auto options = parse_options(argc, argv);
    const std::size_t samples = options.quick ? 48U : 192U * options.scale;
    Timer setup;
    PremixedFlameModel premixed({8.0, 1.1e6, 500.0});
    NonPremixedFlameModel non_premixed({0.42, 0.08, 12.0, 1.1e6, 500.0});
    SootRadiationModel soot({0.015, 2.0, 2.0e4, 0.0});
    const double setup_ms = setup.milliseconds();

    const std::filesystem::path path = options.output.empty() ? std::filesystem::path("reacting_hooks_final.csv")
                                                              : std::filesystem::path(options.output);
    if (path.has_parent_path())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream file(path);
    if (!file) {
        std::cerr << "cannot open output file: " << path << '\n';
        return 2;
    }
    file << "mixture_fraction,premixed_source,nonpremixed_source,heat_release,soot_source,absorption_m1,radiative_"
            "source_w_m3\n";
    Timer simulation;
    double checksum = 0.0;
    for (std::size_t sample = 0; sample < samples; ++sample) {
        const double mixture_fraction = static_cast<double>(sample) / static_cast<double>(samples - 1U);
        const double temperature = 450.0 + 900.0 * mixture_fraction;
        const auto premixed_result = premixed.evaluate({0.35, mixture_fraction, temperature});
        const auto non_premixed_result = non_premixed.evaluate({0.25, mixture_fraction, temperature});
        const auto soot_result = soot.evaluate(
            {0.002 * mixture_fraction, temperature, premixed_result.heat_release + non_premixed_result.heat_release},
            300.0);
        file << mixture_fraction << ',' << premixed_result.progress_source << ',' << non_premixed_result.progress_source
             << ',' << premixed_result.heat_release + non_premixed_result.heat_release << ',' << soot_result.soot_source
             << ',' << soot_result.absorption_coefficient << ',' << soot_result.radiative_source << '\n';
        checksum += soot_result.soot_source + soot_result.radiative_source;
    }
    const double simulation_ms = simulation.milliseconds();
    std::cout << "reacting_hooks samples=" << samples << " csv=" << path << " checksum=" << checksum << '\n';
    emit({"phase03_fvm", "premixed_nonpremixed_soot_radiation_hooks", "cpu", "flame_cells", samples, samples, setup_ms,
          simulation_ms, static_cast<double>(samples), checksum});
    return 0;
}
