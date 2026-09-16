#include "cfd/solvers/fem/darcy2d.hpp"

#include "cfd/fem/reference_element.hpp"

#include <array>
#include <stdexcept>

namespace cfd::fem {
namespace {

Point3 as3(Node2 p){return {p.x,p.y,0.0};}
Node2 centroid(const Mesh2D&mesh,const Tri3&t){const auto&a=mesh.nodes[t.node[0]],&b=mesh.nodes[t.node[1]],&c=mesh.nodes[t.node[2]];return{(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};}

} // namespace

Darcy2D::Darcy2D(Mesh2D mesh,Darcy2DConfig config)
    : config_(config), solver_(std::move(mesh),config.linear){
    if(!(config_.viscosity>0.0)||!(config_.density>=0.0)) throw std::invalid_argument("invalid Darcy properties");
}

void Darcy2D::set_boundary(int patch,ScalarBoundaryCondition2D condition){solver_.set_boundary(patch,std::move(condition));}

void Darcy2D::solve(const std::function<double(Node2)>& permeability,
                    const std::function<double(Node2)>& volumetric_source){
    if(!permeability) throw std::invalid_argument("Darcy permeability callback missing");
    const auto zero=[](Node2){return 0.0;};
    const auto source=volumetric_source?volumetric_source:std::function<double(Node2)>(zero);
    solver_.solve([&](Node2 p){const double k=permeability(p);if(!(k>0.0))throw std::runtime_error("non-positive Darcy permeability");return k/config_.viscosity;},
                  zero,source);
}

std::vector<DarcyVelocity2> Darcy2D::element_velocity(const std::function<double(Node2)>& permeability) const{
    if(!permeability) throw std::invalid_argument("Darcy permeability callback missing");
    std::vector<DarcyVelocity2> velocity;velocity.reserve(mesh().triangles.size());
    for(const auto&tri:mesh().triangles){
        std::array<Point3,3> node{{as3(mesh().nodes[tri.node[0]]),as3(mesh().nodes[tri.node[1]]),as3(mesh().nodes[tri.node[2]])}};
        const auto iso=evaluate_isoparametric(ElementType::tri3,node,{1.0/3.0,1.0/3.0,0.0});
        double gx=0.0,gy=0.0;
        for(std::size_t i=0;i<3U;++i){gx+=pressure()[tri.node[i]]*iso.gradient_physical[i].x;gy+=pressure()[tri.node[i]]*iso.gradient_physical[i].y;}
        const double k=permeability(centroid(mesh(),tri));if(!(k>0.0))throw std::runtime_error("non-positive Darcy permeability");
        const double mobility=k/config_.viscosity;
        velocity.push_back({-mobility*(gx-config_.density*config_.gravity.x),-mobility*(gy-config_.density*config_.gravity.y)});
    }
    return velocity;
}

} // namespace cfd::fem
