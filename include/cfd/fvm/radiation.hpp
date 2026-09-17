#pragma once
namespace cfd::fvm {
inline constexpr double stefan_boltzmann=5.670374419e-8;
[[nodiscard]] double gray_surface_exchange(double temperature_a,double emissivity_a,double temperature_b,double emissivity_b,double view_factor=1.0);
[[nodiscard]] double optically_thin_radiation_source(double temperature,double environment_temperature,double absorption_coefficient);
}
