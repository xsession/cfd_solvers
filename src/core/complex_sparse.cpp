#include "cfd/core/complex_sparse.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::core {

ComplexCsrBuilder::ComplexCsrBuilder(std::size_t dimension)
    :dimension_(dimension),entries_(dimension) {
    if(dimension_==0U)throw std::invalid_argument("complex CSR dimension must be non-zero");
}

void ComplexCsrBuilder::add(std::size_t row,std::size_t column,Complex value){
    if(row>=dimension_||column>=dimension_)throw std::out_of_range("complex CSR builder index out of range");
    if(!std::isfinite(value.real())||!std::isfinite(value.imag()))throw std::invalid_argument("complex CSR entry must be finite");
    entries_[row].emplace_back(column,value);
}

CsrMatrix ComplexCsrBuilder::build_real_block(double drop_tolerance) const {
    if(!(drop_tolerance>=0.0))throw std::invalid_argument("complex CSR drop tolerance must be non-negative");
    CsrBuilder builder(2U*dimension_,2U*dimension_);
    for(std::size_t row=0;row<dimension_;++row){
        for(const auto& [column,value]:entries_[row]){
            const double re=value.real(),im=value.imag();
            // Keep both real-block diagonals structurally present when possible;
            // ILU(0) requires a diagonal in every row.
            if(std::abs(re)>drop_tolerance||row==column){
                const double structural_re=(row==column&&std::abs(re)<=drop_tolerance)
                    ?std::numeric_limits<double>::min():re;
                builder.add(row,column,structural_re);
                builder.add(dimension_+row,dimension_+column,structural_re);
            }
            if(std::abs(im)>drop_tolerance){
                builder.add(row,dimension_+column,-im);
                builder.add(dimension_+row,column,im);
            }
        }
    }
    return builder.build(0.0);
}

ComplexSparseSolveResult solve_complex_sparse(const ComplexCsrBuilder& matrix,
                                               std::span<const Complex> rhs,
                                               const ComplexSparseSolveConfig& config){
    const std::size_t n=matrix.dimension();
    if(rhs.size()!=n)throw std::invalid_argument("complex sparse RHS size mismatch");
    if(config.max_iterations==0U||config.gmres_restart==0U||!(config.relative_tolerance>0.0)
       ||!(config.ilu_diagonal_floor>0.0))throw std::invalid_argument("invalid complex sparse solver controls");
    const auto real_matrix=matrix.build_real_block();
    Ilu0Preconditioner ilu(real_matrix,config.ilu_diagonal_floor);
    std::vector<double> real_rhs(2U*n),x(2U*n,0.0);
    for(std::size_t i=0;i<n;++i){real_rhs[i]=rhs[i].real();real_rhs[n+i]=rhs[i].imag();}
    auto linear=restarted_gmres(real_rhs,x,
        [&](std::span<const double> in,std::span<double> out){real_matrix.multiply(in,out);},
        [&](std::span<const double> in,std::span<double> out){ilu(in,out);},
        config.max_iterations,std::min(config.gmres_restart,2U*n),config.relative_tolerance);
    // Indefinite curl-curl systems can occasionally be a poor match for ILU(0).
    // Retrying with unpreconditioned GMRES keeps the solve sparse and is robust
    // for the moderate reference problems used by the shared EM baseline.
    if(!linear.converged){
        std::fill(x.begin(),x.end(),0.0);
        linear=restarted_gmres(real_rhs,x,
            [&](std::span<const double> in,std::span<double> out){real_matrix.multiply(in,out);},
            [](std::span<const double> in,std::span<double> out){std::copy(in.begin(),in.end(),out.begin());},
            config.max_iterations,std::min(std::max(config.gmres_restart,2U*n),config.max_iterations),config.relative_tolerance);
    }
    ComplexSparseSolveResult result;result.linear_result=linear;result.solution.resize(n);
    for(std::size_t i=0;i<n;++i)result.solution[i]={x[i],x[n+i]};
    return result;
}

} // namespace cfd::core
