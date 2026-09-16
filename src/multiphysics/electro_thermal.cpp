#include "cfd/multiphysics/electro_thermal.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::multiphysics {
namespace {

cfd::fem::Node2 centroid(const cfd::fem::Mesh2D& mesh,const cfd::fem::Tri3& tri) {
    const auto& a=mesh.nodes[tri.node[0]];
    const auto& b=mesh.nodes[tri.node[1]];
    const auto& c=mesh.nodes[tri.node[2]];
    return {(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
}

} // namespace

JouleHeatingCoupler2D::JouleHeatingCoupler2D(cfd::fem::Mesh2D mesh,
                                             cfd::fem::Heat2DConfig thermal_config,
                                             cfd::fem::ScalarDiffusion2DConfig electrical_config)
    : electrical_(mesh,electrical_config),thermal_(std::move(mesh),thermal_config),
      joule_heating_density_(thermal_.mesh().element_count(),0.0) {}

void JouleHeatingCoupler2D::set_electrical_boundary(int patch,cfd::fem::ScalarBoundaryCondition2D condition) {
    electrical_.set_boundary(patch,std::move(condition));
    electrical_solved_=false;
}

void JouleHeatingCoupler2D::initialize_temperature(const std::function<double(cfd::fem::Node2)>& temperature) {
    thermal_.initialize(temperature);
}

void JouleHeatingCoupler2D::set_thermal_dirichlet(
    const std::function<double(cfd::fem::Node2,double)>& boundary_temperature) {
    thermal_.set_dirichlet(boundary_temperature);
}

void JouleHeatingCoupler2D::solve_electrical(
    const std::function<double(cfd::fem::Node2)>& conductivity,
    const std::function<double(cfd::fem::Node2)>& current_source) {
    if(!conductivity) throw std::invalid_argument("electro-thermal coupling requires conductivity");
    electrical_solved_=false;
    electrical_.solve(conductivity,current_source);
    const auto current=electrical_.element_current_density(conductivity);
    const auto& mesh=electrical_.mesh();
    if(current.size()!=mesh.element_count()) throw std::runtime_error("electrical element field size mismatch");
    joule_heating_density_.resize(current.size());
    for(std::size_t e=0;e<current.size();++e){
        const double sigma=conductivity(centroid(mesh,mesh.triangles[e]));
        if(!(sigma>0.0)||!std::isfinite(sigma)) throw std::runtime_error("invalid conductivity in Joule-heating transfer");
        joule_heating_density_[e]=(current[e].x*current[e].x+current[e].y*current[e].y)/sigma;
    }
    electrical_solved_=true;
}

void JouleHeatingCoupler2D::thermal_step(){
    if(!electrical_solved_||joule_heating_density_.size()!=thermal_.mesh().element_count())
        throw std::runtime_error("Joule-heating field is not initialized");
    thermal_.step_element_source(joule_heating_density_);
}

void JouleHeatingCoupler2D::thermal_run(std::size_t steps){
    for(std::size_t i=0;i<steps;++i) thermal_step();
}

} // namespace cfd::multiphysics
