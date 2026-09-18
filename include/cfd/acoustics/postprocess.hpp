#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::acoustics {

[[nodiscard]] double plane_wave_intensity_w_m2(double pressure_rms_pa,double density_kg_m3,double sound_speed_m_s);
[[nodiscard]] double radiation_pressure_pa(double intensity_w_m2,double sound_speed_m_s,double reflection_coefficient=0.0);
[[nodiscard]] double acoustic_heating_w_m3(double intensity_w_m2,double absorption_np_per_m);

struct BeamformSample {
    double x_m{};
    double y_m{};
    std::span<const double> pressure;
};

[[nodiscard]] std::vector<double> delay_and_sum_beamform(std::span<const BeamformSample> sensors,
                                                          double focus_x_m,double focus_y_m,
                                                          double sample_dt,double sound_speed_m_s,
                                                          std::size_t output_samples);

} // namespace cfd::acoustics
