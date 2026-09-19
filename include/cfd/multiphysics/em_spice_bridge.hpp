#pragma once

#include <complex>
#include <cstddef>
#include <functional>

namespace cfd::multiphysics {

struct EmSpiceCouplingConfig {
    std::size_t max_iterations{50U};
    double relaxation{0.7};
    double tolerance{1.0e-10};
};

struct EmSpiceCouplingResult {
    std::complex<double> port_voltage{};
    std::complex<double> port_current{};
    std::size_t iterations{};
    bool converged{};
};

using FieldPortSolve = std::function<std::complex<double>(std::complex<double>)>;
using CircuitPortSolve = std::function<std::complex<double>(std::complex<double>)>;

[[nodiscard]] EmSpiceCouplingResult iterate_em_spice_port(std::complex<double> initial_current,
                                                          const FieldPortSolve& field_solve,
                                                          const CircuitPortSolve& circuit_solve,
                                                          const EmSpiceCouplingConfig& config = {});

} // namespace cfd::multiphysics
