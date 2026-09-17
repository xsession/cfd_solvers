#include "cfd/multibody/system.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfd::multibody;

namespace {
void require(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}
bool near(double a,double b,double tol){return std::abs(a-b)<=tol;}

void test_free_fall_and_quaternion(){
    MultibodyConfig cfg;cfg.dt=1e-3;cfg.gravity={0,-9.81,0};
    RigidBodySystem sys(cfg);RigidBody b;b.mass=2.0;b.inertia_diagonal={0.2,0.3,0.4};b.state.angular_velocity={0,0,4};
    const auto id=sys.add_body(b);sys.run(1000);
    const auto& s=sys.bodies().at(id).state;
    require(near(s.linear_velocity.y,-9.81,2e-10),"free-fall velocity");
    require(std::abs(norm_squared(s.orientation)-1.0)<1e-12,"quaternion normalization");
}

void test_distance_and_revolute_constraints(){
    std::vector<RigidBody> bodies(2);bodies[0].fixed=true;bodies[0].mass=1;bodies[1].mass=1;
    bodies[1].state.position={1.2,0,0};
    auto d=make_distance_constraint(bodies,0,1,{0,0,0},{1.2,0,0},1.0);
    bodies[1].state.linear_velocity={1,0,0};
    auto rows=build_constraint_rows(bodies,Joint{d},1e-2,0.2);solve_constraint_rows(bodies,rows,20);
    require(bodies[1].state.linear_velocity.x<0.0,"distance constraint correction direction");

    bodies[1].state.position={0,0,0};bodies[1].state.linear_velocity={0,0,0};bodies[1].state.angular_velocity={1,2,3};
    auto rev=make_revolute_joint(bodies,0,1,{0,0,0},{0,0,1});
    rows=build_constraint_rows(bodies,Joint{rev},1e-2,0.2);solve_constraint_rows(bodies,rows,30);
    require(std::abs(bodies[1].state.angular_velocity.x)<1e-10,"revolute x lock");
    require(std::abs(bodies[1].state.angular_velocity.y)<1e-10,"revolute y lock");
    require(std::abs(bodies[1].state.angular_velocity.z-3.0)<1e-10,"revolute free axis");
}

void test_collision_and_contact(){
    std::vector<RigidBody> bodies(2);bodies[0].mass=1;bodies[1].mass=1;bodies[0].state.position={-0.45,0,0};bodies[1].state.position={0.45,0,0};
    std::vector<SphereShape> spheres{{0,0.5},{1,0.5}};
    auto pairs=broad_phase_sweep_and_prune(bodies,spheres);require(pairs.size()==1,"broad phase pair");auto bvh=broad_phase_bvh(bodies,spheres);require(bvh==pairs,"BVH broad phase parity");
    auto contacts=sphere_contacts(bodies,spheres,pairs);require(contacts.size()==1,"narrow phase contact");require(near(contacts[0].penetration,0.1,1e-12),"penetration");
    PenaltyContactModel pm;pm.normal_stiffness=1000;pm.normal_damping=0;pm.tangential_damping=0;
    apply_penalty_contacts(bodies,contacts,pm);
    require(bodies[0].force.x<0 && bodies[1].force.x>0,"equal/opposite penalty force");
    require(norm(bodies[0].force+bodies[1].force)<1e-12,"contact momentum conservation");

    bodies[0].force={};bodies[1].force={};bodies[0].state.linear_velocity={1,0,0};bodies[1].state.linear_velocity={-1,0,0};
    ImpulseContactModel im;im.restitution=1.0;im.friction=0;im.iterations=4;im.baumgarte=0;
    resolve_impulse_contacts(bodies,contacts,1e-3,im);
    require(bodies[0].state.linear_velocity.x<0 && bodies[1].state.linear_velocity.x>0,"nonsmooth restitution impulse");
}

void test_dem_and_drag(){
    DemConfig cfg;cfg.dt=2e-5;cfg.gravity={0,-9.81,0};cfg.contact.normal_stiffness=2e5;cfg.contact.normal_damping=80;cfg.contact.friction=0.4;
    ExplicitDemSystem dem(cfg);dem.add_plane({{0,1,0},0});const auto id=dem.add_particle({0,0.12,0},0.1,1000.0);
    dem.run(2500);const auto& b=dem.bodies().at(id);
    require(std::isfinite(b.state.position.y)&&b.state.position.y>0.07,"DEM plane support remains finite");
    RigidBody rb;rb.state.linear_velocity={1,0,0};const auto f=stokes_drag_force(rb,0.01,{0,0,0},1e-3);
    require(f.x<0 && near(f.y,0,1e-15),"Stokes drag direction");
}
}

int main(){
    try{test_free_fall_and_quaternion();test_distance_and_revolute_constraints();test_collision_and_contact();test_dem_and_drag();}
    catch(const std::exception& e){std::cerr<<"v0.16.0 multibody/DEM regression failed: "<<e.what()<<'\n';return 1;}
    std::cout<<"v0.16.0 multibody/DEM regression passed\n";return 0;
}
