#pragma once

#include "cfd/core/conjugate_gradient.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::multiphysics {

[[nodiscard]] double sar_from_rms_electric_field(double conductivity_s_per_m,
                                                  double tissue_density_kg_per_m3,
                                                  double electric_rms_v_per_m);

enum class BioheatBoundary { periodic, fixed_temperature };

struct PennesBioheat2DConfig {
    std::size_t nx{32U},ny{32U};
    double width_m{0.1},height_m{0.1};
    double tissue_density_kg_per_m3{1000.0};
    double tissue_specific_heat_j_per_kg_k{3600.0};
    double thermal_conductivity_w_per_m_k{0.5};
    double blood_density_kg_per_m3{1060.0};
    double blood_specific_heat_j_per_kg_k{3770.0};
    double blood_perfusion_per_s{0.0};
    double blood_temperature_k{310.15};
    double metabolic_heat_w_per_m3{};
    double dt_s{0.1};
    BioheatBoundary boundary{BioheatBoundary::periodic};
    double fixed_boundary_temperature_k{310.15};
    std::size_t max_iterations{4000U};
    double relative_tolerance{1.0e-11};
};

class PennesBioheat2D {
public:
    explicit PennesBioheat2D(PennesBioheat2DConfig config = {});
    void initialize(double temperature_k);
    void set_sar(double sar_w_per_kg);
    void set_sar(std::span<const double> sar_w_per_kg);
    void step(std::size_t steps = 1U);

    [[nodiscard]] const PennesBioheat2DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& temperature_k() const noexcept { return temperature_; }
    [[nodiscard]] const std::vector<double>& sar_w_per_kg() const noexcept { return sar_; }
    [[nodiscard]] double time_s() const noexcept { return time_s_; }
    [[nodiscard]] const cfd::core::ConjugateGradientResult& linear_result() const noexcept { return linear_result_; }

private:
    PennesBioheat2DConfig config_;
    std::vector<double> temperature_;
    std::vector<double> sar_;
    cfd::core::ConjugateGradientWorkspace workspace_;
    cfd::core::ConjugateGradientResult linear_result_{};
    double time_s_{};

    [[nodiscard]] std::size_t index(std::size_t x,std::size_t y) const noexcept { return y*config_.nx+x; }
    [[nodiscard]] bool boundary_node(std::size_t x,std::size_t y) const noexcept;
    void step_once();
};

} // namespace cfd::multiphysics
