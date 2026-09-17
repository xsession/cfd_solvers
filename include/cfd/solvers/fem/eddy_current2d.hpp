#pragma once
#include "cfd/fem/mesh2d.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include <complex>
#include <functional>
#include <vector>
namespace cfd::fem {
struct EddyCurrent2DConfig {double frequency_hz{50.0};std::size_t max_iterations{4000};std::size_t gmres_restart{50};double relative_tolerance{1e-10};};
class EddyCurrent2D {
public:EddyCurrent2D(Mesh2D mesh,EddyCurrent2DConfig config={});void solve(const std::function<double(Node2)>& reluctivity,const std::function<double(Node2)>& conductivity,const std::function<std::complex<double>(Node2)>& source_current);[[nodiscard]] const std::vector<std::complex<double>>& vector_potential()const noexcept{return a_;}[[nodiscard]] const cfd::core::IterativeSolverResult& linear_result()const noexcept{return result_;}[[nodiscard]] const Mesh2D& mesh()const noexcept{return mesh_;}
private:Mesh2D mesh_;EddyCurrent2DConfig config_;std::vector<std::complex<double>>a_;cfd::core::IterativeSolverResult result_{};
};
}
