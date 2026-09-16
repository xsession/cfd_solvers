#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace cfd::core {

struct GeneralizedEigenpair {
    double eigenvalue{};
    std::vector<double> vector;
    std::size_t iterations{};
    double relative_residual{};
    bool converged{};
};

struct GeneralizedEigenConfig {
    std::size_t mode_count{1U};
    std::size_t max_outer_iterations{200U};
    double relative_tolerance{1.0e-9};
    std::size_t linear_max_iterations{4000U};
    double linear_relative_tolerance{1.0e-12};
};

inline double csr_bilinear(const CsrMatrix& matrix,
                           std::span<const double> a,
                           std::span<const double> b,
                           std::vector<double>& workspace) {
    if (a.size()!=b.size() || a.size()!=matrix.rows() || matrix.rows()!=matrix.cols()) {
        throw std::invalid_argument("bilinear form size mismatch");
    }
    workspace.resize(a.size());
    matrix.multiply(b,workspace);
    return parallel_sum(a.size(),[&](std::size_t i){return a[i]*workspace[i];});
}

// Sequential inverse iteration with M-orthogonal deflation for the lowest
// generalized symmetric eigenpairs K phi = lambda M phi. K and M must be SPD
// after essential constraints have been eliminated.
inline std::vector<GeneralizedEigenpair> lowest_generalized_eigenpairs(
    const CsrMatrix& stiffness,
    const CsrMatrix& mass,
    const GeneralizedEigenConfig& config = {}) {
    if(stiffness.rows()!=stiffness.cols()||mass.rows()!=mass.cols()||stiffness.rows()!=mass.rows()) {
        throw std::invalid_argument("generalized eigen matrices must be same-size square matrices");
    }
    const std::size_t n=stiffness.rows();
    if(n==0U||config.mode_count==0U||config.mode_count>n||config.max_outer_iterations==0U
       ||!(config.relative_tolerance>0.0)||config.linear_max_iterations==0U||!(config.linear_relative_tolerance>0.0)) {
        throw std::invalid_argument("invalid generalized eigen controls");
    }
    const auto diagonal=stiffness.diagonal();
    JacobiPreconditioner jacobi(diagonal);
    KrylovWorkspace krylov;
    std::vector<GeneralizedEigenpair> modes;
    modes.reserve(config.mode_count);
    std::vector<double> x(n),rhs(n),y(n),kx(n),mx(n),tmp(n),residual(n);

    auto m_inner=[&](std::span<const double>a,std::span<const double>b){return csr_bilinear(mass,a,b,tmp);};
    auto normalize_m=[&](std::vector<double>&v){
        const double norm2=m_inner(v,v);
        if(!(norm2>0.0)||!std::isfinite(norm2)) throw std::runtime_error("invalid M norm in eigen solve");
        const double inv=1.0/std::sqrt(norm2);
        parallel_for(n,[&](std::size_t i){v[i]*=inv;});
    };
    auto orthogonalize=[&](std::vector<double>&v){
        for(const auto&mode:modes){
            const double alpha=m_inner(mode.vector,v);
            parallel_for(n,[&](std::size_t i){v[i]-=alpha*mode.vector[i];});
        }
    };

    for(std::size_t mode_index=0U;mode_index<config.mode_count;++mode_index){
        for(std::size_t i=0U;i<n;++i){
            const double phase=static_cast<double>((i+1U)*(mode_index+1U));
            x[i]=std::sin(0.731*phase)+0.25*std::cos(0.371*phase);
        }
        orthogonalize(x);normalize_m(x);
        GeneralizedEigenpair pair;pair.vector=x;pair.relative_residual=std::numeric_limits<double>::infinity();
        for(std::size_t iteration=0U;iteration<config.max_outer_iterations;++iteration){
            mass.multiply(x,rhs);
            std::fill(y.begin(),y.end(),0.0);
            const auto linear=preconditioned_conjugate_gradient(
                rhs,y,
                [&](std::span<const double>v,std::span<double>out){stiffness.multiply(v,out);},
                [&](std::span<const double>r,std::span<double>z){jacobi(r,z);},
                krylov,config.linear_max_iterations,config.linear_relative_tolerance);
            if(!linear.converged){pair.iterations=iteration;pair.converged=false;break;}
            orthogonalize(y);normalize_m(y);x.swap(y);
            stiffness.multiply(x,kx);mass.multiply(x,mx);
            const double numerator=parallel_sum(n,[&](std::size_t i){return x[i]*kx[i];});
            const double denominator=parallel_sum(n,[&](std::size_t i){return x[i]*mx[i];});
            const double lambda=numerator/denominator;
            parallel_for(n,[&](std::size_t i){residual[i]=kx[i]-lambda*mx[i];});
            const double scale=std::max(std::abs(lambda)*vector_rms(mx),std::numeric_limits<double>::min());
            const double rel=vector_rms(residual)/scale;
            pair.eigenvalue=lambda;pair.vector=x;pair.iterations=iteration+1U;pair.relative_residual=rel;
            if(rel<=config.relative_tolerance){pair.converged=true;break;}
        }
        if(!pair.converged) return modes;
        modes.push_back(std::move(pair));
    }
    return modes;
}

} // namespace cfd::core
