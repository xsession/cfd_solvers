#include "cfd/solvers/fem/magnetostatics2d.hpp"

#include "cfd/fem/reference_element.hpp"

#include <array>
#include <stdexcept>

namespace cfd::fem {

Magnetostatics2D::Magnetostatics2D(Mesh2D mesh,ScalarDiffusion2DConfig config):solver_(std::move(mesh),config){}
void Magnetostatics2D::set_boundary(int patch,ScalarBoundaryCondition2D condition){solver_.set_boundary(patch,std::move(condition));}
void Magnetostatics2D::solve(const std::function<double(Node2)>& reluctivity,const std::function<double(Node2)>& current_density_z){
    if(!reluctivity||!current_density_z)throw std::invalid_argument("magnetostatic callback missing");
    solver_.solve(reluctivity,[](Node2){return 0.0;},current_density_z);
}
std::vector<MagneticFluxDensity2> Magnetostatics2D::element_flux_density() const{
    std::vector<MagneticFluxDensity2> out;out.reserve(mesh().triangles.size());
    for(const auto&tri:mesh().triangles){
        std::array<Point3,3> node{{{mesh().nodes[tri.node[0]].x,mesh().nodes[tri.node[0]].y,0.0},
                                  {mesh().nodes[tri.node[1]].x,mesh().nodes[tri.node[1]].y,0.0},
                                  {mesh().nodes[tri.node[2]].x,mesh().nodes[tri.node[2]].y,0.0}}};
        const auto iso=evaluate_isoparametric(ElementType::tri3,node,{1.0/3.0,1.0/3.0,0.0});
        double gx=0.0,gy=0.0;
        for(std::size_t i=0;i<3U;++i){gx+=vector_potential()[tri.node[i]]*iso.gradient_physical[i].x;gy+=vector_potential()[tri.node[i]]*iso.gradient_physical[i].y;}
        out.push_back({gy,-gx});
    }
    return out;
}

} // namespace cfd::fem
