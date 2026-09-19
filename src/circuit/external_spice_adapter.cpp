#include "cfd/circuit/external_spice_adapter.hpp"

#include <stdexcept>

namespace cfd::circuit {

ExternalSpiceRunPlan make_external_spice_run_plan(ExternalSpiceEngine engine, std::string netlist,
                                                  std::string raw_output, std::vector<std::string> extra_arguments) {
    if (netlist.empty() || raw_output.empty()) {
        throw std::invalid_argument("external SPICE plan paths must not be empty");
    }
    ExternalSpiceRunPlan plan;
    plan.engine = engine;
    plan.netlist = std::move(netlist);
    plan.raw_output = std::move(raw_output);
    plan.argv.push_back(engine == ExternalSpiceEngine::ngspice ? "ngspice" : "Xyce");
    if (engine == ExternalSpiceEngine::ngspice) {
        plan.argv.push_back("-b");
        plan.argv.push_back("-r");
        plan.argv.push_back(plan.raw_output);
        plan.argv.push_back(plan.netlist);
    } else {
        plan.argv.insert(plan.argv.end(), extra_arguments.begin(), extra_arguments.end());
        plan.argv.push_back(plan.netlist);
    }
    return plan;
}

} // namespace cfd::circuit
