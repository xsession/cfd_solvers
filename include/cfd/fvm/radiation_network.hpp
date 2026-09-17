#pragma once
#include <span>
#include <vector>
namespace cfd::fvm {
struct RadiosityResult { std::vector<double> radiosity; std::vector<double> net_heat; };
// Dense diffuse-gray enclosure baseline. view_factors is row-major NxN and
// must satisfy non-negative row sums close to one; net_heat is positive out.
[[nodiscard]] RadiosityResult solve_gray_radiosity(std::span<const double> temperature,
    std::span<const double> emissivity,std::span<const double> area,std::span<const double> view_factors);
class RadiationSourceModel {public:virtual ~RadiationSourceModel()=default;[[nodiscard]] virtual double volumetric_source(double temperature,double environment_temperature)const=0;};
class OpticallyThinRadiation final:public RadiationSourceModel{public:explicit OpticallyThinRadiation(double absorption):absorption_(absorption){}[[nodiscard]] double volumetric_source(double temperature,double environment_temperature)const override;private:double absorption_{};};
} // namespace cfd::fvm
