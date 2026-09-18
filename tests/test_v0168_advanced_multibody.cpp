#include "cfd/multibody/advanced_dynamics.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfd::multibody;

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
bool near(double a,double b,double tol){return std::abs(a-b)<=tol;}

double oscillator_energy(const ImplicitGeneralizedIntegrator& solver,double stiffness){
    return 0.5*solver.v()[0]*solver.v()[0]+0.5*stiffness*solver.q()[0]*solver.q()[0];
}

void test_implicit_newmark(){
    constexpr double k=100.0;
    auto force=[](std::span<const double> q,std::span<const double>,double,std::span<double> f){f[0]=-k*q[0];};
    auto controls=ImplicitSecondOrderConfig::newmark(0.01);
    controls.nonlinear_tolerance=1e-12;
    ImplicitGeneralizedIntegrator solver({1.0},1,force,controls);
    solver.initialize({1.0},{0.0});
    const double e0=oscillator_energy(solver,k);
    solver.run(1000);
    const double e1=oscillator_energy(solver,k);
    require(std::abs(e1-e0)/e0<2e-9,"Newmark average-acceleration energy conservation");
    require(solver.last_residual_norm()<1e-9,"Newmark nonlinear residual");
}

void test_hht_algorithmic_damping(){
    constexpr double k=40000.0;
    auto force=[](std::span<const double> q,std::span<const double>,double,std::span<double> f){f[0]=-k*q[0];};
    auto nc=ImplicitSecondOrderConfig::newmark(0.01);nc.nonlinear_tolerance=1e-11;
    auto hc=ImplicitSecondOrderConfig::hht(0.01,-0.2);hc.nonlinear_tolerance=1e-11;
    ImplicitGeneralizedIntegrator newmark({1.0},1,force,nc),hht({1.0},1,force,hc);
    newmark.initialize({1.0},{0.0});hht.initialize({1.0},{0.0});
    newmark.run(100);hht.run(100);
    const double en=oscillator_energy(newmark,k),eh=oscillator_energy(hht,k);
    require(eh<0.15*en,"HHT high-frequency algorithmic damping");

    const auto ga=ImplicitSecondOrderConfig::generalized_alpha(0.01,0.5);
    require(ga.beta>0.0&&ga.gamma>=0.5,"generalized-alpha parameter construction");
}

void test_articulated_reduced_coordinates(){
    ArticulatedSystem slider({0,0,0});
    ArticulatedLink link;link.joint_type=ReducedJointType::prismatic;link.axis_parent={1,0,0};link.body.mass=2.0;link.body.inertia_diagonal={1,1,1};link.applied_effort=4.0;
    slider.add_link(link);
    const auto m=slider.mass_matrix();
    require(m.size()==1&&near(m[0],2.0,1e-12),"prismatic reduced mass");
    slider.step(0.01);
    require(near(slider.links()[0].qdd,2.0,1e-12),"prismatic generalized acceleration");
    require(near(slider.links()[0].qd,0.02,1e-12),"prismatic generalized velocity");
    require(near(slider.links()[0].body.state.position.x,0.0002,1e-12),"prismatic forward kinematics");

    ArticulatedSystem arm({0,-9.81,0});
    ArticulatedLink a;a.joint_type=ReducedJointType::revolute;a.axis_parent={0,0,1};a.joint_to_com_local={0.5,0,0};a.body.mass=1.5;a.body.inertia_diagonal={0.1,0.1,0.2};
    const auto root=arm.add_link(a);
    ArticulatedLink b;b.parent=root;b.joint_type=ReducedJointType::revolute;b.axis_parent={0,0,1};b.parent_joint_offset={0.5,0,0};b.joint_to_com_local={0.5,0,0};b.body.mass=1.0;b.body.inertia_diagonal={0.08,0.08,0.12};
    arm.add_link(b);
    const auto mm=arm.mass_matrix();
    require(mm.size()==4,"two-link mass matrix size");
    require(mm[0]>0&&mm[3]>0&&near(mm[1],mm[2],1e-12)&&std::abs(mm[1])>1e-6,"articulated symmetric coupled mass matrix");
    const auto gf=arm.generalized_forces();
    require(gf[0]<0.0&&gf[1]<0.0,"articulated gravity projection");
    arm.step(1e-3);
    require(std::isfinite(arm.links()[0].q)&&std::isfinite(arm.links()[1].q),"articulated finite integration");
}

void test_flexible_fem_multibody_interface(){
    std::vector<FlexibleInterfaceNode> nodes{{{-1,0,0}},{{1,0,0}}};
    FlexibleMode bending;bending.shape_local={{0,1,0},{0,1,0}};bending.modal_mass=2.0;bending.modal_stiffness=100.0;
    FlexibleBodyInterface interface(nodes,{bending});
    interface.set_modal_state({0.1},{0.2});

    RigidBody frame;frame.mass=5.0;frame.inertia_diagonal={2,2,2};frame.state.position={1,2,0};frame.state.linear_velocity={0.5,0,0};frame.state.angular_velocity={0,0,0.3};
    const auto positions=interface.node_positions(frame);
    require(near(positions[0].x,0.0,1e-12)&&near(positions[0].y,2.1,1e-12),"flexible displacement transfer");
    const auto velocities=interface.node_velocities(frame);
    require(near(velocities[0].x,0.47,1e-12)&&near(velocities[0].y,-0.1,1e-12),"flexible velocity transfer");

    const std::vector<Vec3> nodal_forces{{0,2,0},{0,-1,0}};
    const auto load=interface.project_nodal_forces(frame,nodal_forces);
    require(near(load.body_force.y,1.0,1e-12),"FEM reaction force conservation");
    require(near(load.body_torque.z,-3.0,1e-12),"FEM reaction torque conservation");
    require(load.modal_force.size()==1&&near(load.modal_force[0],1.0,1e-12),"modal force projection");

    const Vec3 virtual_translation{0.02,-0.01,0};
    const Vec3 virtual_rotation{0,0,0.03};
    const double virtual_mode=0.04;
    double nodal_work=0.0;
    for(std::size_t i=0;i<nodes.size();++i){
        const Vec3 r=positions[i]-frame.state.position;
        const Vec3 modal_world=rotate(frame.state.orientation,bending.shape_local[i])*virtual_mode;
        const Vec3 dx=virtual_translation+cross(virtual_rotation,r)+modal_world;
        nodal_work+=dot(nodal_forces[i],dx);
    }
    const double reduced_work=dot(load.body_force,virtual_translation)+dot(load.body_torque,virtual_rotation)+load.modal_force[0]*virtual_mode;
    require(near(nodal_work,reduced_work,1e-12),"FEM/multibody virtual-work conservation");

    const auto applied=interface.apply_nodal_forces(frame,nodal_forces);
    require(near(frame.force.y,applied.body_force.y,1e-12)&&near(frame.torque.z,applied.body_torque.z,1e-12),"FEM reaction applied to multibody frame");
}
}

int main(){
    try{
        test_implicit_newmark();
        test_hht_algorithmic_damping();
        test_articulated_reduced_coordinates();
        test_flexible_fem_multibody_interface();
    }catch(const std::exception& e){
        std::cerr<<"v0.16.8 advanced multibody regression failed: "<<e.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.16.8 advanced multibody regression passed\n";
    return 0;
}
