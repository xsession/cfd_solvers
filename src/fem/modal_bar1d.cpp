#include "cfd/solvers/fem/modal_bar1d.hpp"

#include "cfd/core/csr_matrix.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {

BarModal1D::BarModal1D(BarModal1DConfig config):config_(config){
    if(config_.elements<2U||!(config_.length>0.0)||!(config_.area>0.0)||!(config_.young_modulus>0.0)||!(config_.density>0.0)) {
        throw std::invalid_argument("invalid modal bar properties");
    }
}

std::vector<BarMode1D> BarModal1D::solve(std::size_t mode_count) const{
    const std::size_t n=config_.elements; // free nodes 1..elements
    const double h=config_.length/static_cast<double>(config_.elements);
    const double ke=config_.young_modulus*config_.area/h;
    const double me=config_.density*config_.area*h/6.0;
    cfd::core::CsrBuilder kb(n,n),mb(n,n);
    for(std::size_t e=0U;e<config_.elements;++e){
        const std::size_t global[2]{e,e+1U};
        const double k[2][2]{{ke,-ke},{-ke,ke}};
        const double m[2][2]{{2.0*me,me},{me,2.0*me}};
        for(std::size_t i=0;i<2U;++i){
            if(global[i]==0U) continue;
            const std::size_t row=global[i]-1U;
            for(std::size_t j=0;j<2U;++j){
                if(global[j]==0U) continue;
                const std::size_t col=global[j]-1U;
                kb.add(row,col,k[i][j]);mb.add(row,col,m[i][j]);
            }
        }
    }
    auto eigen_config=config_.eigen;eigen_config.mode_count=mode_count;
    const auto pairs=cfd::core::lowest_generalized_eigenpairs(kb.build(),mb.build(),eigen_config);
    if(pairs.size()!=mode_count) throw std::runtime_error("bar generalized eigen solve did not converge");
    std::vector<BarMode1D> out;out.reserve(mode_count);
    constexpr double two_pi=6.283185307179586476925286766559;
    for(const auto&pair:pairs){
        if(!(pair.eigenvalue>0.0)) throw std::runtime_error("non-positive bar eigenvalue");
        const double omega=std::sqrt(pair.eigenvalue);
        BarMode1D mode;mode.angular_frequency=omega;mode.frequency_hz=omega/two_pi;mode.displacement.assign(config_.elements+1U,0.0);
        for(std::size_t i=0U;i<n;++i) mode.displacement[i+1U]=pair.vector[i];
        mode.eigenpair=pair;out.push_back(std::move(mode));
    }
    return out;
}

} // namespace cfd::fem
