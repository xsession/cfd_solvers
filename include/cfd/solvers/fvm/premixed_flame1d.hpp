#pragma once
#include <cstddef>
#include <vector>
namespace cfd::fvm {
struct PremixedFlame1DConfig {std::size_t cells{400};double length{1.0};double diffusivity{1e-3};double reaction_rate{1.0};double dt{1e-4};double unburned_temperature{300};double burned_temperature{1800};};
// Fisher-KPP progress-variable laminar-flame reference. c=0 unburned, c=1 burned.
class PremixedFlame1D {
public:explicit PremixedFlame1D(PremixedFlame1DConfig config={});void initialize_front(double location,double thickness);void step(std::size_t count=1);[[nodiscard]] double front_location()const;[[nodiscard]] double theoretical_speed()const;[[nodiscard]] const std::vector<double>& progress()const noexcept{return c_;}[[nodiscard]] double temperature(std::size_t i)const;
private:PremixedFlame1DConfig cfg_;std::vector<double>c_,next_;
};
} // namespace cfd::fvm
