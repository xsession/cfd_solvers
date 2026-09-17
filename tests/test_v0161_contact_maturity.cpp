#include "cfd/multibody/system.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfd::multibody;

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
bool near(double a,double b,double tol){return std::abs(a-b)<=tol;}

std::vector<Vec3> box_vertices(double hx,double hy,double hz){
    std::vector<Vec3> vertices;
    for(double x:{-hx,hx}) for(double y:{-hy,hy}) for(double z:{-hz,hz}) vertices.push_back({x,y,z});
    return vertices;
}

void test_linear_and_angular_motors(){
    std::vector<RigidBody> bodies(2);
    bodies[0].fixed=true;
    bodies[1].mass=2.0;
    bodies[1].inertia_diagonal={2.0,2.0,2.0};

    auto linear=make_linear_motor(bodies,0,1,{1,0,0},2.0,1000.0);
    auto rows=build_constraint_rows(bodies,Joint{linear},0.01,0.0);
    solve_constraint_rows(bodies,rows,8);
    require(near(bodies[1].state.linear_velocity.x,2.0,1e-12),"linear motor target speed");

    bodies[1].state.angular_velocity={};
    auto angular=make_angular_motor(bodies,0,1,{0,0,1},3.0,1000.0);
    rows=build_constraint_rows(bodies,Joint{angular},0.01,0.0);
    solve_constraint_rows(bodies,rows,8);
    require(near(bodies[1].state.angular_velocity.z,3.0,1e-12),"angular motor target speed");

    bodies[1].state.linear_velocity={};
    linear=make_linear_motor(bodies,0,1,{1,0,0},100.0,2.0);
    rows=build_constraint_rows(bodies,Joint{linear},0.01,0.0);
    solve_constraint_rows(bodies,rows,1);
    // max impulse = F*dt = .02 N s, dv = J/m = .01 m/s.
    require(near(bodies[1].state.linear_velocity.x,0.01,1e-12),"linear motor force limit");
}

void test_convex_gjk_epa(){
    std::vector<RigidBody> bodies(2);
    bodies[0].state.position={0,0,0};
    bodies[1].state.position={0.8,0,0};
    ConvexHullShape a{0,box_vertices(0.5,0.5,0.5)};
    ConvexHullShape b{1,box_vertices(0.5,0.5,0.5)};

    const auto hit=convex_gjk_epa_contact(bodies,a,b,64,128,1e-10);
    require(hit.intersect,"GJK overlapping boxes");
    require(hit.contact.penetration>0.19 && hit.contact.penetration<0.21,"EPA penetration depth");
    require(hit.contact.normal.x < -0.9,"EPA normal force direction on body A");
    require(std::abs(hit.contact.normal.y)<1e-6 && std::abs(hit.contact.normal.z)<1e-6,"EPA axis-aligned normal");
    require(hit.gjk_iterations>0 && hit.epa_iterations>0,"GJK/EPA iteration diagnostics");

    bodies[1].state.position={1.2,0,0};
    const auto miss=convex_gjk_epa_contact(bodies,a,b);
    require(!miss.intersect,"GJK separated boxes");

    bodies[1].state.position={0.65,0.15,0};
    bodies[1].state.orientation=quaternion_from_axis_angle({0,0,1},0.25);
    const auto rotated=convex_gjk_epa_contact(bodies,a,b,64,128,1e-9);
    require(rotated.intersect && rotated.contact.penetration>0.0,"GJK/EPA rotated convex overlap");
    require(near(norm(rotated.contact.normal),1.0,1e-9),"EPA normalized contact normal");
}


void test_triangle_mesh_contact(){
    std::vector<RigidBody> bodies(1);
    bodies[0].state.position={0.2,0.08,0.2};
    std::vector<SphereShape> spheres{{0,0.1}};
    TriangleMeshShape mesh;
    mesh.body=invalid_body;
    mesh.feature_namespace=9U;
    mesh.vertices={{0,0,0},{1,0,0},{0,0,1}};
    mesh.triangles={{{0,1,2}}};
    const auto contacts=sphere_triangle_mesh_contacts(bodies,spheres,mesh);
    require(contacts.size()==1,"sphere-triangle contact");
    require(near(contacts[0].penetration,0.02,1e-12),"sphere-triangle penetration");
    require(contacts[0].normal.y>0.999,"sphere-triangle normal");
    require((contacts[0].feature_id>>32U)==9U,"triangle feature namespace");

    ContactManifoldCache cache;
    cache.update(contacts);
    require(cache.states().size()==1 && cache.states()[0].age==1,"triangle manifold insert");
    cache.states()[0].tangential_displacement={0.01,0,0};
    cache.update(contacts);
    require(cache.states()[0].age==2 && cache.states()[0].tangential_displacement.x>0.0099,"triangle manifold persists");
}

void test_persistent_manifold_and_mindlin_history(){
    std::vector<RigidBody> bodies(1);
    bodies[0].mass=1.0;
    bodies[0].inertia_diagonal={1,1,1};
    bodies[0].state.linear_velocity={1,0,0};
    std::vector<SphereShape> spheres{{0,0.5}};
    std::vector<ContactPoint> contacts{{0,invalid_body,{0,-0.5,0},{0,1,0},0.01,7U}};

    HertzMindlinContactModel model;
    model.normal_stiffness=1000.0;
    model.normal_damping=0.0;
    model.tangential_stiffness=100.0;
    model.tangential_damping=0.0;
    model.friction=100.0;
    model.rolling_resistance=0.0;

    ContactManifoldCache cache;
    apply_history_dependent_mindlin_contacts(bodies,contacts,spheres,cache,0.01,model);
    require(cache.states().size()==1,"persistent manifold insertion");
    require(cache.states()[0].age==1,"persistent manifold initial age");
    require(cache.states()[0].tangential_displacement.x>0.0099,"Mindlin tangential spring accumulation");
    require(bodies[0].force.x<0.0,"Mindlin restoring force opposes slip");

    bodies[0].clear_accumulators();
    apply_history_dependent_mindlin_contacts(bodies,contacts,spheres,cache,0.01,model);
    require(cache.states()[0].age==2,"persistent manifold age increments");
    require(cache.states()[0].tangential_displacement.x>0.0199,"Mindlin history persists across steps");
    require(std::abs(cache.states()[0].tangential_displacement.y)<1e-12,"Mindlin displacement remains tangent");

    cache.update({});
    require(cache.states().empty(),"stale manifold removal");
}
}

int main(){
    try{
        test_linear_and_angular_motors();
        test_convex_gjk_epa();
        test_triangle_mesh_contact();
        test_persistent_manifold_and_mindlin_history();
    }catch(const std::exception& e){
        std::cerr<<"v0.16.1 contact maturity regression failed: "<<e.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.16.1 contact maturity regression passed\n";
    return 0;
}
