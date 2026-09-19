#include "common/benchmark.hpp"
#include "cfd/rf/surface_mom.hpp"

#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    using namespace cfd::rf;
    const auto options = parse_options(argc, argv);
    using P = cfd::fem::Point3;
    const std::vector<SurfaceTriangle> patch = {
        {P{0.0, 0.0, 0.0}, P{0.10, 0.0, 0.0}, P{0.10, 0.10, 0.0}},
        {P{0.0, 0.0, 0.0}, P{0.10, 0.10, 0.0}, P{0.0, 0.10, 0.0}},
    };
    SurfaceMomConfig config;
    config.frequency_hz = options.quick ? 1.0e9 : 2.4e9;
    config.excitation = SurfaceMomExcitation::delta_gap;
    config.feed_edge = 2U;
    config.feed_voltage_v = {1.0, 0.0};

    Timer setup;
    const auto setup_ms = setup.milliseconds();
    Timer simulation;
    const auto result = solve_pec_surface_mom(patch, config);
    const auto far_field = surface_mom_far_field(result, {0.0, 1.0, 1.0}, 1.0);
    const auto simulation_ms = simulation.milliseconds();

    if (!options.output.empty()) {
        const std::filesystem::path output_path(options.output);
        if (!output_path.parent_path().empty())
            std::filesystem::create_directories(output_path.parent_path());
        std::ofstream output(output_path);
        if (!output)
            return 2;
        output << std::setprecision(17) << "edge,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,current_real_a,current_imag_a\n";
        for (std::size_t edge = 0; edge < result.basis.size(); ++edge) {
            const auto& basis = result.basis[edge];
            output << edge << ',' << basis.edge_start.x << ',' << basis.edge_start.y << ',' << basis.edge_start.z << ','
                   << basis.edge_end.x << ',' << basis.edge_end.y << ',' << basis.edge_end.z << ','
                   << result.current_a[edge].real() << ',' << result.current_a[edge].imag() << '\n';
        }
    }
    double checksum = std::abs(result.feed_impedance_ohm) + far_field.power_density_w_m2 + result.relative_residual;
    for (const auto current : result.current_a)
        checksum += std::abs(current);
    emit({"phase12_em", "pec_surface_rwg_mom", "cpu", "rwg_unknowns", result.basis.size(), 1U, setup_ms, simulation_ms,
          static_cast<double>(result.basis.size()), checksum});
    return 0;
}
