#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <span>
#include <string_view>
#include <vector>
namespace cfd::fvm {
enum class TurbulenceModel {laminar,mixing_length,smagorinsky,wale};
[[nodiscard]] TurbulenceModel parse_turbulence_model(std::string_view name);
struct TurbulenceConfig {TurbulenceModel model{TurbulenceModel::laminar};double mixing_length{0.01};double smagorinsky_constant{0.17};double wale_constant{0.325};};
[[nodiscard]] std::vector<double> eddy_viscosity(const PolyMesh& mesh,std::span<const Vec3> velocity,const TurbulenceConfig& config={});
}
