#pragma once
#include "cfd/core/csr_matrix.hpp"
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::fvm {
class FvScalarMatrix {
public:
    explicit FvScalarMatrix(std::size_t cells);
    void add_diagonal(std::size_t cell,double value);
    void add_coefficient(std::size_t row,std::size_t column,double value);
    void add_source(std::size_t cell,double value);
    [[nodiscard]] cfd::core::CsrMatrix matrix()const;
    [[nodiscard]] std::span<const double> source()const noexcept{return source_;}
private:cfd::core::CsrBuilder builder_;std::vector<double>source_;
};
}
