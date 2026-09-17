#include "cfd/multiphysics/block_system.hpp"
#include <stdexcept>
namespace cfd::multiphysics {
namespace {
std::size_t total(const std::vector<std::size_t>& sizes){std::size_t n=0;for(auto s:sizes){if(s==0)throw std::invalid_argument("block size must be positive");n+=s;}if(n==0)throw std::invalid_argument("block system requires blocks");return n;}
}
BlockCoupledLinearSystem::BlockCoupledLinearSystem(std::vector<std::size_t> sizes)
    :sizes_(std::move(sizes)),offsets_(sizes_.size()+1,0),builder_(total(sizes_),total(sizes_)),rhs_(total(sizes_),0.0){
    for(std::size_t b=0;b<sizes_.size();++b)offsets_[b+1]=offsets_[b]+sizes_[b];
}
std::size_t BlockCoupledLinearSystem::index(std::size_t b,std::size_t i) const{if(b>=sizes_.size()||i>=sizes_[b])throw std::out_of_range("block system index");return offsets_[b]+i;}
void BlockCoupledLinearSystem::add(std::size_t rb,std::size_t r,std::size_t cb,std::size_t c,double v){builder_.add(index(rb,r),index(cb,c),v);}
void BlockCoupledLinearSystem::add_rhs(std::size_t b,std::size_t r,double v){rhs_[index(b,r)]+=v;}
cfd::core::CsrMatrix BlockCoupledLinearSystem::matrix() const{return builder_.build();}
cfd::core::IterativeSolverResult BlockCoupledLinearSystem::solve(std::span<double> x,std::size_t maxit,std::size_t restart,double tol) const{
    if(x.size()!=total_size())throw std::invalid_argument("block solution size mismatch");
    auto A=matrix();cfd::core::Ilu0Preconditioner ilu(A);
    return cfd::core::restarted_gmres(rhs_,x,[&](auto in,auto out){A.multiply(in,out);},[&](auto in,auto out){ilu(in,out);},maxit,restart,tol);
}
} // namespace cfd::multiphysics
