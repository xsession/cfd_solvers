#include "cfd/solvers/fem/elasticity2d.hpp"

#include "cfd/core/csr_matrix.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::fem {
namespace {
struct Geometry { double area; std::array<std::array<double,2>,3> grad; Node2 centroid; };
Geometry geometry(const Mesh2D& m,const Tri3& t){
    const auto&a=m.nodes[t.node[0]],&b=m.nodes[t.node[1]],&c=m.nodes[t.node[2]];
    const double det=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x); const double area=0.5*std::abs(det);
    if (!(area > 0.0)) throw std::runtime_error("degenerate elasticity triangle");
    const double inv = 1.0 / det;
    return {area,{{{(b.y-c.y)*inv,(c.x-b.x)*inv},{(c.y-a.y)*inv,(a.x-c.x)*inv},{(a.y-b.y)*inv,(b.x-a.x)*inv}}},
            {(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0}};
}
std::array<std::array<double,3>,3> constitutive(const Elasticity2DConfig& c){
    std::array<std::array<double,3>,3>D{}; const double E=c.youngs_modulus,nu=c.poisson_ratio;
    if(c.mode==ElasticityMode2D::planeStress){
        const double f=E/(1.0-nu*nu); D={{{f,f*nu,0.0},{f*nu,f,0.0},{0.0,0.0,f*(1.0-nu)/2.0}}};
    }else{
        const double f=E/((1.0+nu)*(1.0-2.0*nu));
        D={{{f*(1.0-nu),f*nu,0.0},{f*nu,f*(1.0-nu),0.0},{0.0,0.0,f*(1.0-2.0*nu)/2.0}}};
    } return D;
}
}

Elasticity2D::Elasticity2D(Mesh2D mesh,Elasticity2DConfig config)
 :mesh_(std::move(mesh)),config_(config),displacement_(mesh_.node_count()),
  prescribed_x_(mesh_.node_count(),std::numeric_limits<double>::quiet_NaN()),
  prescribed_y_(mesh_.node_count(),std::numeric_limits<double>::quiet_NaN()){
    mesh_.validate();
    if(!(config_.youngs_modulus>0.0)||!(config_.poisson_ratio>-1.0&&config_.poisson_ratio<0.5)||config_.max_iterations==0U||!(config_.relative_tolerance>0.0))
        throw std::invalid_argument("invalid Elasticity2D controls");
}
void Elasticity2D::set_dirichlet(const std::function<bool(Node2)>& predicate,
                                  const std::function<Displacement2(Node2)>& value,
                                  bool fix_x,bool fix_y){
    for(std::size_t i=0;i<mesh_.node_count();++i) if(predicate(mesh_.nodes[i])){
        const auto v=value(mesh_.nodes[i]); if(fix_x) prescribed_x_[i]=v.x; if(fix_y) prescribed_y_[i]=v.y;
    }
}
void Elasticity2D::solve(const std::function<Displacement2(Node2)>& body_force){
    solve_impl(body_force,{},0.0,0.0);
}
void Elasticity2D::solve_thermal(std::span<const double> temperature,double expansion,double reference,
                                const std::function<Displacement2(Node2)>& body_force){
    if(temperature.size()!=mesh_.node_count()||!std::isfinite(expansion)||!std::isfinite(reference))
        throw std::invalid_argument("invalid thermoelastic input");
    for(double t:temperature)if(!std::isfinite(t))throw std::invalid_argument("non-finite thermal field");
    solve_impl(body_force,temperature,expansion,reference);
}
void Elasticity2D::solve_impl(const std::function<Displacement2(Node2)>& body_force,
                             std::span<const double> temperature,double expansion,double reference){
    const std::size_t ndof=2U*mesh_.node_count(); std::vector<std::size_t> map(ndof,static_cast<std::size_t>(-1));
    std::size_t nf=0U;
    for(std::size_t n=0;n<mesh_.node_count();++n){
        if(std::isfinite(prescribed_x_[n])) displacement_[n].x=prescribed_x_[n]; else map[2U*n]=nf++;
        if(std::isfinite(prescribed_y_[n])) displacement_[n].y=prescribed_y_[n]; else map[2U*n+1U]=nf++;
    }
    if(nf==0U) throw std::runtime_error("elasticity system has no free DOFs");
    cfd::core::CsrBuilder builder(nf,nf); std::vector<double> rhs(nf,0.0); const auto D=constitutive(config_);
    for(const auto&t:mesh_.triangles){
        const auto g=geometry(mesh_,t); double B[3][6]{};
        for(std::size_t i=0;i<3U;++i){ B[0][2*i]=g.grad[i][0]; B[1][2*i+1]=g.grad[i][1]; B[2][2*i]=g.grad[i][1]; B[2][2*i+1]=g.grad[i][0]; }
        double ke[6][6]{};
        for(std::size_t i=0;i<6U;++i) for(std::size_t j=0;j<6U;++j)
            for(std::size_t a=0;a<3U;++a) for(std::size_t b=0;b<3U;++b) ke[i][j]+=g.area*B[a][i]*D[a][b]*B[b][j];
        const auto bf=body_force?body_force(g.centroid):Displacement2{};
        double thermal_strain=0.0;
        if(!temperature.empty()) thermal_strain=expansion*((temperature[t.node[0]]+temperature[t.node[1]]+temperature[t.node[2]])/3.0-reference);
        if(config_.mode==ElasticityMode2D::planeStrain)thermal_strain*=1.0+config_.poisson_ratio;
        std::array<double,3> thermal_stress{};
        for(std::size_t a=0;a<3;++a)thermal_stress[a]=(D[a][0]+D[a][1])*thermal_strain;
        for(std::size_t i=0;i<3U;++i) for(std::size_t comp=0;comp<2U;++comp){
            const std::size_t gi=2U*t.node[i]+comp; if(map[gi]==static_cast<std::size_t>(-1)) continue; const auto row=map[gi];
            rhs[row]+=(comp==0?bf.x:bf.y)*g.area/3.0;
            for(std::size_t a=0;a<3;++a)rhs[row]+=g.area*B[a][2U*i+comp]*thermal_stress[a];
            for(std::size_t j=0;j<3U;++j) for(std::size_t cj=0;cj<2U;++cj){
                const std::size_t gj=2U*t.node[j]+cj; const double kij=ke[2U*i+comp][2U*j+cj];
                if(map[gj]==static_cast<std::size_t>(-1)){
                    const double prescribed=cj==0?displacement_[t.node[j]].x:displacement_[t.node[j]].y; rhs[row]-=kij*prescribed;
                } else builder.add(row,map[gj],kij);
            }
        }
    }
    const auto A=builder.build(); cfd::core::JacobiPreconditioner jacobi(A.diagonal()); cfd::core::KrylovWorkspace ws; std::vector<double>x(nf,0.0);
    linear_result_=cfd::core::preconditioned_conjugate_gradient(rhs,x,[&](std::span<const double>v,std::span<double>o){A.multiply(v,o);},
        [&](std::span<const double>r,std::span<double>z){jacobi(r,z);},ws,config_.max_iterations,config_.relative_tolerance);
    if(!linear_result_.converged) throw std::runtime_error("Elasticity2D PCG did not converge");
    for(std::size_t n=0;n<mesh_.node_count();++n){if(map[2*n]!=static_cast<std::size_t>(-1))displacement_[n].x=x[map[2*n]];if(map[2*n+1]!=static_cast<std::size_t>(-1))displacement_[n].y=x[map[2*n+1]];}
}
} // namespace cfd::fem
