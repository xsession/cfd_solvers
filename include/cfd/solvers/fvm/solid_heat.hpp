#pragma once

#include "cfd/solvers/fvm/scalar_transport.hpp"

#include <functional>
#include <string_view>

namespace cfd::fvm {

struct SolidHeatConfig {
    double dt{1.0e-3};
    double density{1.0};
    double heat_capacity{1.0};
    double conductivity{1.0};
    TemporalScheme temporal_scheme{TemporalScheme::backward_bdf2};
    double crank_nicolson_off_centering{1.0};
    std::size_t linear_iterations{400};
    std::size_t gmres_restart{30};
    double linear_tolerance{1.0e-10};
};

enum class SolidThermalBoundaryType { adiabatic, fixed_temperature };

class SolidHeatConduction {
public:
    explicit SolidHeatConduction(PolyMesh mesh, SolidHeatConfig config = {});
    void initialize(double temperature);
    void initialize(const std::function<double(Vec3)>& temperature);
    void set_boundary(std::string_view patch, SolidThermalBoundaryType type, double temperature = 0.0);
    void set_volumetric_heat_source(double watts_per_cubic_metre);
    void set_volumetric_heat_source(const std::function<double(Vec3)>& watts_per_cubic_metre);
    cfd::core::IterativeSolverResult step();
    cfd::core::IterativeSolverResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return solver_.mesh(); }
    [[nodiscard]] const std::vector<double>& temperature() const noexcept { return solver_.values(); }
    [[nodiscard]] double time() const noexcept { return solver_.time(); }
    [[nodiscard]] double total_sensible_energy() const;
    [[nodiscard]] const SolidHeatConfig& config() const noexcept { return config_; }
private:
    SolidHeatConfig config_;
    ScalarTransport solver_;
};

} // namespace cfd::fvm
