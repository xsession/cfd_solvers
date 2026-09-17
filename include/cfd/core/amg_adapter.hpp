#pragma once
#include <functional>
#include <span>
#include <stdexcept>
#include <utility>
namespace cfd::core {
class AmgPreconditionerAdapter {
public:
    using Apply=std::function<void(std::span<const double>,std::span<double>)>;
    explicit AmgPreconditionerAdapter(Apply apply):apply_(std::move(apply)){if(!apply_)throw std::invalid_argument("AMG adapter requires apply callback");}
    void operator()(std::span<const double> r,std::span<double> z)const{if(r.size()!=z.size())throw std::invalid_argument("AMG adapter size mismatch");apply_(r,z);}
private:Apply apply_;
};
} // namespace cfd::core
