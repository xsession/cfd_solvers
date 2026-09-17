#include "cfd/core/complex_sparse.hpp"
#include "cfd/em/edge_fem3d.hpp"
#include "cfd/em/adaptivity3d.hpp"
#include "cfd/em/wave_port.hpp"
#include "cfd/em/cable.hpp"
#include "cfd/particle/electromagnetic.hpp"
#include "cfd/particle/wakefield.hpp"
#include "cfd/solvers/fem/magnetostatics2d.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>
int main(){
    using C=std::complex<double>;
    cfd::core::ComplexCsrBuilder a(1U);a.add(0,0,C{2,3});const std::array<C,1> b{{C{1,-1}}};auto x=cfd::core::solve_complex_sparse(a,b);if(!x.linear_result.converged)throw std::runtime_error("complex sparse");
    auto mesh3=cfd::fem::make_box_tet_mesh(2,2,2);cfd::em::EdgeMaxwell3DConfig ec;ec.frequency_hz=5e7;ec.linear.relative_tolerance=1e-8;ec.linear.gmres_restart=200;auto field=cfd::em::solve_driven_edge_maxwell_3d(mesh3,ec,[](auto p)->cfd::em::ComplexVec3{return {C{},C{0,std::sin(3.141592653589793*p.x)},C{}};});auto ind=cfd::em::maxwell_face_jump_indicators_3d(mesh3,field);auto mark=cfd::em::mark_maxwell_dorfler(ind,0.5);auto refined=cfd::em::refine_tet4_marked_longest_edges(mesh3,mark);refined.validate();
    auto modes=cfd::em::rectangular_waveguide_modes(22.86e-3,10.16e-3,10e9);if(std::abs(cfd::em::rectangular_mode_power_w(modes.front(),22.86e-3,10.16e-3)-1.0)>1e-9)throw std::runtime_error("wave port");
    cfd::em::MulticonductorRlcg line{1,{0},{250e-9},{0},{100e-12}};cfd::rf::ComplexMatrix load(1);load(0,0)=50.0;std::array<C,1> vin{{1.0}};auto cable=cfd::em::solve_multiconductor_cable(0.1,50,1e9,line,vin,load);if(cable.load_voltage_v.empty())throw std::runtime_error("cable");
    cfd::particle::ChargedParticle p;p.velocity_m_per_s={2e8,0,0};for(int i=0;i<100;++i)cfd::particle::relativistic_boris_push(p,{{0,0,0},{0,0,1}},1e-4);std::vector<double>w(256);for(std::size_t i=0;i<w.size();++i)w[i]=cfd::particle::resonator_longitudinal_wake_v_per_c(i*1e-4,100,1e9,5);(void)cfd::particle::wake_impedance_ohm(w,1e-4,1e9);
    auto mesh2=cfd::fem::make_rectangle_tri_mesh(4,4);cfd::fem::Magnetostatics2D mag(mesh2);for(int patch=0;patch<4;++patch)mag.set_boundary(patch,{cfd::fem::ScalarBoundaryType::dirichlet,[](auto){return 0.0;},{}});cfd::fem::PiecewiseLinearBHCurve bh({0,1},{0,1e5});auto nr=mag.solve_nonlinear([&](double B){return bh.reluctivity(B);},[](auto){return 1e3;},5,1e-6,1.0);if(!nr.converged)throw std::runtime_error("nonlinear magnetics");
    return 0;
}
