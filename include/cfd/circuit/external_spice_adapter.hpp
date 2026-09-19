#pragma once

#include <string>
#include <vector>

namespace cfd::circuit {

enum class ExternalSpiceEngine { ngspice, xyce };

// Reproducible command description for optional external cross-validation.
// The clean solver never invokes the process; callers decide whether the
// selected engine is installed and how its stdout/stderr are supervised.
struct ExternalSpiceRunPlan {
    ExternalSpiceEngine engine{};
    std::vector<std::string> argv;
    std::string netlist;
    std::string raw_output;
};

[[nodiscard]] ExternalSpiceRunPlan make_external_spice_run_plan(ExternalSpiceEngine engine, std::string netlist,
                                                                std::string raw_output,
                                                                std::vector<std::string> extra_arguments = {});

} // namespace cfd::circuit
