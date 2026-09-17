#pragma once
#include "cfd/core/soa_field.hpp"
#include "cfd/core/grid.hpp"
#include "cfd/solvers/lbm/d2q9.hpp"
#include <array>
#include <cstddef>
namespace cfd::lbm {
struct AdvancedD2Q9Config {std::size_t nx{64},ny{64};float tau{0.7F};};
class RegularizedD2Q9Solver {
public:explicit RegularizedD2Q9Solver(AdvancedD2Q9Config config={});void initialize_uniform(float rho=1,float ux=0,float uy=0);void initialize_taylor_green(float amplitude=.03F);void step(std::size_t count=1);[[nodiscard]] double mass()const;[[nodiscard]] double kinetic_energy()const;
private:AdvancedD2Q9Config c_;cfd::core::Grid2D grid_;cfd::core::StaticSoA<float,9>f_,next_;void step_once();};
class MrtD2Q9Solver {
public:explicit MrtD2Q9Solver(AdvancedD2Q9Config config={});void initialize_uniform(float rho=1,float ux=0,float uy=0);void initialize_taylor_green(float amplitude=.03F);void step(std::size_t count=1);[[nodiscard]] double mass()const;[[nodiscard]] double kinetic_energy()const;
private:AdvancedD2Q9Config c_;cfd::core::Grid2D grid_;cfd::core::StaticSoA<float,9>f_,next_;std::array<std::array<double,9>,9>inv_{};void step_once();};
}
