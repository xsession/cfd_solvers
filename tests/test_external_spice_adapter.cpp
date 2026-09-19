#include "cfd/circuit/external_spice_adapter.hpp"

#include <cassert>
#include <iostream>

int main() {
    using cfd::circuit::ExternalSpiceEngine;
    const auto ng = cfd::circuit::make_external_spice_run_plan(ExternalSpiceEngine::ngspice, "case.cir", "case.raw");
    assert(ng.argv.size() == 5U && ng.argv[0] == "ngspice" && ng.argv[1] == "-b" && ng.argv[2] == "-r" &&
           ng.argv[3] == "case.raw" && ng.argv[4] == "case.cir");
    const auto xyce = cfd::circuit::make_external_spice_run_plan(ExternalSpiceEngine::xyce, "case.cir", "case.raw",
                                                                 {"-o", "case.prn"});
    assert(xyce.argv.size() == 4U && xyce.argv[0] == "Xyce" && xyce.argv[1] == "-o" && xyce.argv[2] == "case.prn" &&
           xyce.argv[3] == "case.cir");
    std::cout << "external ngspice/Xyce adapter regression passed\n";
    return 0;
}
