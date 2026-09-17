#pragma once
#include "cfd/fem/mesh2d.hpp"
#include "cfd/solvers/fem/elasticity2d.hpp"
#include <functional>
#include <vector>
namespace cfd::fem {
struct Brinkman2DConfig {double viscosity{1e-3};double effective_viscosity{1e-3};double permeability{1e-8};std::size_t max_iterations{3000};double tolerance{1e-10};};
class Brinkman2D {
public:
    explicit Brinkman2D(Mesh2D mesh,Brinkman2DConfig config={});
    void set_zero_velocity_boundary(bool all_boundary=true){zero_boundary_=all_boundary;}
    void solve(const std::function<Displacement2(Node2)>& body_force);
    [[nodiscard]] const std::vector<Displacement2>& velocity()const noexcept{return velocity_;}
private:Mesh2D mesh_;Brinkman2DConfig c_;std::vector<Displacement2>velocity_;bool zero_boundary_{true};
};
} // namespace cfd::fem
