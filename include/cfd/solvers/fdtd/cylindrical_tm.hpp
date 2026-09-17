#pragma once
#include "cfd/core/aligned_allocator.hpp"
#include <cstddef>
namespace cfd::fdtd {
struct CylindricalTMConfig { std::size_t nr{64},nz{128}; double dr{1e-3},dz{1e-3},courant{0.45},epsilon_r{1.0},mu_r{1.0}; };
// Axisymmetric m=0 TM(r,z) Maxwell FDTD: Er, Ez, Hphi. PEC outer-r and z boundaries;
// the r=0 Ez update uses the analytic cylindrical-axis limit.
class CylindricalTM {
public:
    explicit CylindricalTM(CylindricalTMConfig config={});
    void initialize_gaussian_ez(double r0,double z0,double width,double amplitude=1.0);
    void add_soft_ez_source(std::size_t ir,std::size_t iz,double value);
    void step(std::size_t count=1);
    [[nodiscard]] double energy() const;
    [[nodiscard]] double max_field() const;
    [[nodiscard]] double dt() const noexcept{return dt_;}
    [[nodiscard]] const cfd::core::AlignedVector<double>& ez()const noexcept{return ez_;}
private:
    CylindricalTMConfig c_;double dt_{},eps_{},mu_{};
    cfd::core::AlignedVector<double> er_,ez_,hphi_;
    std::size_t er_idx(std::size_t i,std::size_t j)const noexcept{return j*(c_.nr-1)+i;}
    std::size_t ez_idx(std::size_t i,std::size_t j)const noexcept{return j*c_.nr+i;}
    std::size_t h_idx(std::size_t i,std::size_t j)const noexcept{return j*(c_.nr-1)+i;}
    void step_once();
};
} // namespace cfd::fdtd
