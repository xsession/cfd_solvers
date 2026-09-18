#include "common/benchmark.hpp"
#include "cfd/chemistry/kinetics.hpp"
#include "cfd/chemistry/implicit_reactor.hpp"
#include "cfd/electrochemistry/corrosion_models.hpp"

#include <numeric>
#include <utility>
#include <vector>

int main(int argc, char** argv) {
    using namespace cfd::examples;
    const auto options = parse_options(argc, argv);
    const std::size_t cells = options.quick ? 64U : 1024U * options.scale;
    const std::size_t steps = options.quick ? 10U : 1000U;

    Timer setup;
    cfd::chemistry::ReactionNetwork network({{"A", 1.0, 0}, {"B", 1.0, 0}});
    cfd::chemistry::ElementaryReaction reaction;
    reaction.reactants = {{0U, 1.0}};
    reaction.products = {{1U, 1.0}};
    reaction.forward = {20.0, 0.0, 4000.0};
    network.add_reaction(std::move(reaction));

    std::vector<double> concentration{1.0, 0.0};
    cfd::electrochemistry::PhaseFieldCorrosion1D phase(cells, 1.0e-6, 1.0e-3, 1.0e-10);
    phase.initialize_interface(0.35 * static_cast<double>(cells) * 1.0e-6, 5.0e-6);
    std::vector<double> drive(cells, 0.2);
    const double setup_ms = setup.milliseconds();

    Timer simulation;
    const auto reactor_result = cfd::chemistry::integrate_isothermal(network, concentration, 800.0, 0.5);
    for (std::size_t i = 0; i < steps; ++i) phase.step(1.0e-3, drive);
    const double simulation_ms = simulation.milliseconds();

    const double checksum = std::accumulate(
        phase.order_parameter().begin(), phase.order_parameter().end(), concentration[0] + concentration[1]);
    emit({"phase07_chemistry_corrosion",
          "reactor_plus_phase_field_corrosion",
          "cpu",
          "cell_steps",
          cells,
          steps,
          setup_ms,
          simulation_ms,
          static_cast<double>(cells) * static_cast<double>(steps) + static_cast<double>(reactor_result.accepted_steps),
          checksum});
    return 0;
}
