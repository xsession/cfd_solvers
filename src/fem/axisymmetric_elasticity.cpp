#include "cfd/solvers/fem/axisymmetric_elasticity.hpp"

#include "cfd/core/csr_matrix.hpp"
#include "cfd/fem/reference_element.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::fem {
namespace {

constexpr double pi = 3.1415926535897932384626433832795;

std::array<std::array<double,4>,4> axisymmetric_constitutive(double E,double nu) {
    const double lambda=E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double mu=E/(2.0*(1.0+nu));
    std::array<std::array<double,4>,4>D{};
    for(std::size_t i=0;i<3U;++i) for(std::size_t j=0;j<3U;++j) D[i][j]=lambda;
    for(std::size_t i=0;i<3U;++i) D[i][i]+=2.0*mu;
    D[3][3]=mu;
    return D;
}

} // namespace

AxisymmetricElasticity::AxisymmetricElasticity(Mesh2D mesh,AxisymmetricElasticityConfig config)
    :mesh_(std::move(mesh)),config_(config),displacement_(mesh_.node_count()),
     prescribed_r_(mesh_.node_count(),std::numeric_limits<double>::quiet_NaN()),
     prescribed_z_(mesh_.node_count(),std::numeric_limits<double>::quiet_NaN()) {
    mesh_.validate();
    if(!(config_.youngs_modulus>0.0)||!(config_.poisson_ratio>-1.0&&config_.poisson_ratio<0.5)||
       config_.max_iterations==0U||!(config_.relative_tolerance>0.0)) throw std::invalid_argument("invalid axisymmetric elasticity controls");
    for(const auto&p:mesh_.nodes) if(p.x<-1.0e-14) throw std::invalid_argument("axisymmetric FEM requires non-negative radial coordinate");
}

void AxisymmetricElasticity::set_dirichlet(const std::function<bool(Node2)>& predicate,
                                            const std::function<Displacement2(Node2)>& value,
                                            bool fix_radial,bool fix_axial) {
    for(std::size_t i=0;i<mesh_.node_count();++i) if(predicate(mesh_.nodes[i])) {
        const auto v=value(mesh_.nodes[i]);
        if(fix_radial) prescribed_r_[i]=v.x;
        if(fix_axial) prescribed_z_[i]=v.y;
    }
}

void AxisymmetricElasticity::solve(const std::function<Displacement2(Node2)>& body_force) {
    const std::size_t ndof=2U*mesh_.node_count();
    std::vector<std::size_t> map(ndof,static_cast<std::size_t>(-1));
    std::size_t nf=0U;
    for(std::size_t n=0;n<mesh_.node_count();++n) {
        if(std::isfinite(prescribed_r_[n])) displacement_[n].x=prescribed_r_[n]; else map[2U*n]=nf++;
        if(std::isfinite(prescribed_z_[n])) displacement_[n].y=prescribed_z_[n]; else map[2U*n+1U]=nf++;
    }
    if(nf==0U) throw std::runtime_error("axisymmetric elasticity system has no free DOFs");
    cfd::core::CsrBuilder builder(nf,nf);
    std::vector<double> rhs(nf,0.0);
    const auto D=axisymmetric_constitutive(config_.youngs_modulus,config_.poisson_ratio);
    const auto qrule=gaussian_quadrature(ElementType::tri3,2U);

    for(const auto&tri:mesh_.triangles) {
        std::array<Point3,3> nodes{{
            {mesh_.nodes[tri.node[0]].x,mesh_.nodes[tri.node[0]].y,0.0},
            {mesh_.nodes[tri.node[1]].x,mesh_.nodes[tri.node[1]].y,0.0},
            {mesh_.nodes[tri.node[2]].x,mesh_.nodes[tri.node[2]].y,0.0}}};
        double ke[6][6]{};
        double fe[6]{};
        for(const auto&qp:qrule) {
            const auto shape=evaluate_shape(ElementType::tri3,qp.point);
            const auto iso=evaluate_isoparametric(ElementType::tri3,nodes,qp.point);
            const double r=iso.point.x;
            if(!(r>1.0e-14)) continue;
            double B[4][6]{};
            for(std::size_t i=0;i<3U;++i) {
                const auto&g=iso.gradient_physical[i];
                B[0][2U*i]=g.x;
                B[1][2U*i+1U]=g.y;
                B[2][2U*i]=shape.value[i]/r;
                B[3][2U*i]=g.y;
                B[3][2U*i+1U]=g.x;
            }
            const double w=2.0*pi*r*qp.weight*iso.measure;
            for(std::size_t i=0;i<6U;++i) for(std::size_t j=0;j<6U;++j)
                for(std::size_t a=0;a<4U;++a) for(std::size_t b=0;b<4U;++b) ke[i][j]+=w*B[a][i]*D[a][b]*B[b][j];
            const auto bf=body_force?body_force({iso.point.x,iso.point.y}):Displacement2{};
            for(std::size_t i=0;i<3U;++i) {
                fe[2U*i]+=w*shape.value[i]*bf.x;
                fe[2U*i+1U]+=w*shape.value[i]*bf.y;
            }
        }
        for(std::size_t i=0;i<3U;++i) for(std::size_t comp=0;comp<2U;++comp) {
            const std::size_t gi=2U*tri.node[i]+comp;
            if(map[gi]==static_cast<std::size_t>(-1)) continue;
            const auto row=map[gi]; rhs[row]+=fe[2U*i+comp];
            for(std::size_t j=0;j<3U;++j) for(std::size_t cj=0;cj<2U;++cj) {
                const std::size_t gj=2U*tri.node[j]+cj;
                const double kij=ke[2U*i+comp][2U*j+cj];
                if(map[gj]==static_cast<std::size_t>(-1)) {
                    const double prescribed=cj==0U?displacement_[tri.node[j]].x:displacement_[tri.node[j]].y;
                    rhs[row]-=kij*prescribed;
                } else builder.add(row,map[gj],kij);
            }
        }
    }
    const auto A=builder.build();
    cfd::core::JacobiPreconditioner jacobi(A.diagonal());
    cfd::core::KrylovWorkspace ws;
    std::vector<double>x(nf,0.0);
    linear_result_=cfd::core::preconditioned_conjugate_gradient(rhs,x,
        [&](std::span<const double>v,std::span<double>o){A.multiply(v,o);},
        [&](std::span<const double>r,std::span<double>z){jacobi(r,z);},ws,
        config_.max_iterations,config_.relative_tolerance);
    if(!linear_result_.converged) throw std::runtime_error("AxisymmetricElasticity PCG did not converge");
    for(std::size_t n=0;n<mesh_.node_count();++n) {
        if(map[2U*n]!=static_cast<std::size_t>(-1)) displacement_[n].x=x[map[2U*n]];
        if(map[2U*n+1U]!=static_cast<std::size_t>(-1)) displacement_[n].y=x[map[2U*n+1U]];
    }
}

} // namespace cfd::fem
