#pragma once
#include "cfd/fem/mesh2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include <functional>
#include <vector>
namespace cfd::fem {
struct PenaltyStokes2DConfig{double viscosity{1.0};double divergence_penalty{1000.0};std::size_t max_iterations{5000};double relative_tolerance{1e-10};};
class PenaltyStokes2D {
public:PenaltyStokes2D(Mesh2D mesh,PenaltyStokes2DConfig config={});void set_dirichlet(std::function<Displacement2(Node2)> velocity);void solve(std::function<Displacement2(Node2)> body_force);[[nodiscard]] const std::vector<Displacement2>& velocity()const noexcept{return velocity_;}[[nodiscard]] double divergence_l2()const;
private:Mesh2D mesh_;PenaltyStokes2DConfig config_;std::function<Displacement2(Node2)> boundary_;std::vector<Displacement2>velocity_;
};
} // namespace cfd::fem
