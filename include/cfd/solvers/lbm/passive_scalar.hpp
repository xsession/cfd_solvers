#pragma once
#include "cfd/core/parallel.hpp"
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
namespace cfd::lbm {
struct PassiveScalarConfig { std::size_t nx{64},ny{64};double diffusivity{0.01};double ux{},uy{}; };
// Periodic D2Q5 BGK advection-diffusion distribution. Temperature is an alias
// of the same conserved scalar equation; no buoyancy feedback is implied.
class PassiveScalarD2Q5 {
public:
    explicit PassiveScalarD2Q5(PassiveScalarConfig config):config_(config),cells_(config.nx*config.ny),g_(5*cells_),next_(5*cells_){if(config.nx<2||config.ny<2||!(config.diffusivity>0))throw std::invalid_argument("invalid passive scalar lattice");tau_=0.5+3.0*config.diffusivity;initialize([](double,double){return 0.0;});}
    template<class F> void initialize(F value){cfd::core::parallel_for(cells_,[&](std::size_t n){const auto x=n%config_.nx,y=n/config_.nx;const double phi=value((x+.5)/config_.nx,(y+.5)/config_.ny);for(int q=0;q<5;++q)g(q,n)=equilibrium(q,phi);});}
    void step(std::size_t count=1){for(std::size_t it=0;it<count;++it){const double omega=1.0/tau_;cfd::core::parallel_for(cells_,[&](std::size_t n){const std::size_t x=n%config_.nx,y=n/config_.nx;std::array<double,5> fin{};double phi=0;for(int q=0;q<5;++q){const auto xs=shift(x,-cx[q],config_.nx),ys=shift(y,-cy[q],config_.ny);fin[q]=g(q,ys*config_.nx+xs);phi+=fin[q];}for(int q=0;q<5;++q)next_[q*cells_+n]=fin[q]-omega*(fin[q]-equilibrium(q,phi));});g_.swap(next_);}}
    [[nodiscard]] std::vector<double> scalar() const {std::vector<double> out(cells_);cfd::core::parallel_for(cells_,[&](std::size_t n){double s=0;for(int q=0;q<5;++q)s+=g(q,n);out[n]=s;});return out;}
    [[nodiscard]] double total() const {return cfd::core::parallel_sum(cells_,[&](std::size_t n){double s=0;for(int q=0;q<5;++q)s+=g(q,n);return s;});}
    [[nodiscard]] const PassiveScalarConfig& config()const noexcept{return config_;}
private:
    PassiveScalarConfig config_;std::size_t cells_{};double tau_{};std::vector<double> g_,next_;
    static constexpr int cx[5]{0,1,0,-1,0},cy[5]{0,0,1,0,-1};static constexpr double w[5]{1.0/3.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};
    static std::size_t shift(std::size_t v,int delta,std::size_t n)noexcept{return static_cast<std::size_t>((static_cast<long long>(v)+delta+static_cast<long long>(n))%static_cast<long long>(n));}
    double& g(int q,std::size_t n){return g_[static_cast<std::size_t>(q)*cells_+n];} double g(int q,std::size_t n)const{return g_[static_cast<std::size_t>(q)*cells_+n];}
    double equilibrium(int q,double phi)const noexcept{return w[q]*phi*(1.0+3.0*(cx[q]*config_.ux+cy[q]*config_.uy));}
};
using ThermalD2Q5=PassiveScalarD2Q5;
} // namespace cfd::lbm
