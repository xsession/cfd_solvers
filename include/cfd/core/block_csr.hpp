#pragma once
#include "cfd/core/parallel.hpp"
#include <array>
#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
namespace cfd::core {
template<std::size_t B> class BlockCsrMatrix {
public:
    using Block=std::array<double,B*B>;
    BlockCsrMatrix()=default;
    BlockCsrMatrix(std::size_t rows,std::size_t cols,std::vector<std::size_t> offsets,std::vector<std::size_t> columns,std::vector<Block> values)
      :rows_(rows),cols_(cols),offsets_(std::move(offsets)),columns_(std::move(columns)),values_(std::move(values)){validate();}
    [[nodiscard]] std::size_t block_rows() const noexcept{return rows_;} [[nodiscard]] std::size_t block_cols() const noexcept{return cols_;}
    [[nodiscard]] std::size_t block_nonzeros() const noexcept{return values_.size();}
    void multiply(std::span<const double> x,std::span<double> y) const {if(x.size()!=cols_*B||y.size()!=rows_*B)throw std::invalid_argument("block CSR multiply size mismatch");parallel_for(rows_,[&](std::size_t r){std::array<double,B> sum{};for(std::size_t k=offsets_[r];k<offsets_[r+1];++k){const auto c=columns_[k];const auto& a=values_[k];for(std::size_t i=0;i<B;++i)for(std::size_t j=0;j<B;++j)sum[i]+=a[i*B+j]*x[c*B+j];}for(std::size_t i=0;i<B;++i)y[r*B+i]=sum[i];});}
private:
    std::size_t rows_{},cols_{};std::vector<std::size_t> offsets_,columns_;std::vector<Block> values_;
    void validate() const {if(offsets_.size()!=rows_+1||columns_.size()!=values_.size()||offsets_.empty()||offsets_.front()!=0||offsets_.back()!=values_.size())throw std::invalid_argument("invalid block CSR structure");for(auto c:columns_)if(c>=cols_)throw std::invalid_argument("block CSR column out of range");}
};
template<std::size_t B> class BlockCsrBuilder {
public:
    using Block=typename BlockCsrMatrix<B>::Block;
    BlockCsrBuilder(std::size_t rows,std::size_t cols):rows_(rows),cols_(cols),entries_(rows){}
    void add(std::size_t row,std::size_t col,const Block& block){if(row>=rows_||col>=cols_)throw std::out_of_range("block CSR index");entries_[row].push_back({col,block});}
    [[nodiscard]] BlockCsrMatrix<B> build() const {std::vector<std::size_t> off(rows_+1);std::vector<std::size_t> cols;std::vector<Block> vals;for(std::size_t r=0;r<rows_;++r){auto e=entries_[r];std::sort(e.begin(),e.end(),[](const auto&a,const auto&b){return a.first<b.first;});for(std::size_t i=0;i<e.size();){auto col=e[i].first;Block block{};do{for(std::size_t q=0;q<B*B;++q)block[q]+=e[i].second[q];++i;}while(i<e.size()&&e[i].first==col);cols.push_back(col);vals.push_back(block);}off[r+1]=vals.size();}return {rows_,cols_,std::move(off),std::move(cols),std::move(vals)};}
private:std::size_t rows_,cols_;std::vector<std::vector<std::pair<std::size_t,Block>>> entries_;
};
} // namespace cfd::core
