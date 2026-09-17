#pragma once
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include <cstddef>
#include <functional>
#include <vector>
namespace cfd::fem {
struct DynamicBar1DConfig {std::size_t elements{80};double length{1.0},area{1.0},young_modulus{200e9},density{7800},dt{1e-5},newmark_beta{0.25},newmark_gamma{0.5};std::size_t linear_iterations{2000};double linear_tolerance{1e-11};};
class DynamicBar1D {
public:
    explicit DynamicBar1D(DynamicBar1DConfig config={});
    void initialize(const std::function<double(double)>& displacement,const std::function<double(double)>& velocity={});
    void step(double tip_force=0.0);void run(std::size_t steps,double tip_force=0.0);
    [[nodiscard]] const std::vector<double>& displacement()const noexcept{return u_;}[[nodiscard]] const std::vector<double>& velocity()const noexcept{return v_;}[[nodiscard]] double time()const noexcept{return time_;}
    [[nodiscard]] double kinetic_energy()const;[[nodiscard]] double strain_energy()const;
private:DynamicBar1DConfig c_;cfd::core::CsrMatrix k_,m_,effective_;std::vector<double>u_,v_,a_,rhs_,candidate_;double time_{};void build_matrices();
};
} // namespace cfd::fem
