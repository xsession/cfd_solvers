#include "cfd/fvm/fv_matrix.hpp"
#include <stdexcept>
namespace cfd::fvm {
FvScalarMatrix::FvScalarMatrix(std::size_t n):builder_(n,n),source_(n){if(n==0)throw std::invalid_argument("FVM matrix requires cells");}
void FvScalarMatrix::add_diagonal(std::size_t c,double v){builder_.add(c,c,v);}void FvScalarMatrix::add_coefficient(std::size_t r,std::size_t c,double v){builder_.add(r,c,v);}void FvScalarMatrix::add_source(std::size_t c,double v){if(c>=source_.size())throw std::out_of_range("FVM source index");source_[c]+=v;}cfd::core::CsrMatrix FvScalarMatrix::matrix()const{return builder_.build();}
}
