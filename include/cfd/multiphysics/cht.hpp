#pragma once
namespace cfd::multiphysics {
struct InterfaceHeatTransfer {double heat_flux{};double interface_temperature{};};
[[nodiscard]] InterfaceHeatTransfer conjugate_interface_heat_flux(double temperature_a,double conductivity_a,double distance_a,double temperature_b,double conductivity_b,double distance_b);
}
