#pragma once
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::core {
struct MultigridResult {bool converged{};std::size_t cycles{};double relative_residual{};};
class PoissonMultigrid2D {
public:
    explicit PoissonMultigrid2D(std::size_t interior_points_per_axis);
    [[nodiscard]] MultigridResult solve(std::span<const double> rhs,std::span<double> solution,std::size_t max_cycles=50,double tolerance=1e-10)const;
    [[nodiscard]] std::size_t size()const noexcept{return n_;}
private:std::size_t n_{};static void smooth(std::size_t n,const std::vector<double>&f,std::vector<double>&u,int sweeps);static void residual(std::size_t n,const std::vector<double>&f,const std::vector<double>&u,std::vector<double>&r);static void vcycle(std::size_t n,const std::vector<double>&f,std::vector<double>&u);
};
} // namespace cfd::core
