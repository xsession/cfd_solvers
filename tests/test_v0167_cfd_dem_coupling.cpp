#include "cfd/multibody/cfd_dem_coupling.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/multibody/resident_cfd_dem_sycl.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
void require(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
bool near(double a,double b,double tol=1.0e-11){return std::abs(a-b)<=tol;}
bool near_vec(cfd::multibody::Vec3 a,cfd::multibody::Vec3 b,double tol=1.0e-10){
    return near(a.x,b.x,tol)&&near(a.y,b.y,tol)&&near(a.z,b.z,tol);
}
}

int main(){
    using namespace cfd::multibody;
    using cfd::fvm::Vec3;

    const auto mesh=cfd::fvm::make_cartesian_hexa_mesh(4U,2U,2U,1.0,1.0,1.0);

    // Resolved traction primitive: a single quadrature element gives an exact
    // force and moment, independent of interpolation details.
    std::vector<double> pressure(mesh.cell_count(),3.0);
    std::vector<SymmetricStress3> stress(mesh.cell_count());
    RigidSurfaceSample sample{{0.2,1.0,0.0},{1.0,0.0,0.0},2.0};
    const auto load=integrate_resolved_surface_load(mesh,pressure,stress,
        std::span<const RigidSurfaceSample>(&sample,1U),{0.2,0.0,0.0});
    require(near_vec(load.force,{-6.0,0.0,0.0}),"resolved pressure traction force");
    require(near_vec(load.torque,{0.0,0.0,6.0}),"resolved pressure traction torque");

    RigidBody body; body.state.position={0.5,0.5,0.5}; body.mass=2.0;
    auto surface=sphere_surface_quadrature(body,0.1,6U);
    const auto closed=integrate_resolved_surface_load(mesh,pressure,stress,surface,body.state.position);
    require(near_vec(closed.force,{},1.0e-12),"uniform pressure on closed sphere has zero resultant");
    require(near_vec(closed.torque,{},1.0e-12),"uniform pressure on closed sphere has zero torque");
    apply_resolved_surface_load(body,load);
    require(near_vec(body.force,load.force)&&near_vec(body.torque,load.torque),"resolved load application");

    // Unresolved many-particle coupling conserves particle volume and momentum.
    std::vector<RigidBody> bodies(2);
    bodies[0].state.position={0.15,0.25,0.25}; bodies[0].state.linear_velocity={0.0,0.0,0.0}; bodies[0].mass=1.0;
    bodies[1].state.position={0.78,0.75,0.75}; bodies[1].state.linear_velocity={0.2,0.0,0.0}; bodies[1].mass=1.0;
    std::vector<SphereShape> spheres{{0U,0.04},{1U,0.05}};
    std::vector<Vec3> velocity(mesh.cell_count(),{1.0,0.0,0.0});
    std::vector<Vec3> gradp(mesh.cell_count(),{0.2,0.0,0.0});
    std::vector<Vec3> vort(mesh.cell_count(),{0.0,0.0,2.0});
    UnresolvedCfdDemConfig cfg; cfg.fluid_density=1.2; cfg.dynamic_viscosity=1.8e-5;
    auto coupled=unresolved_cfd_dem_coupling(mesh,bodies,spheres,velocity,gradp,vort,cfg);
    require(coupled.particle_force.size()==2U,"two particle coupling forces");

    double deposited=0.0;
    cfd::multibody::Vec3 fluid_reaction{},particle_total{};
    for(std::size_t c=0;c<mesh.cell_count();++c){
        deposited+=coupled.solid_volume_fraction[c]*mesh.cells()[c].volume;
        fluid_reaction+=coupled.fluid_momentum_source[c]*mesh.cells()[c].volume;
        require(coupled.void_fraction[c]>=cfg.minimum_void_fraction&&coupled.void_fraction[c]<=1.0,"physical void fraction");
    }
    for(const auto f:coupled.particle_force) particle_total+=f;
    const double exact_volume=(4.0/3.0)*std::numbers::pi*(0.04*0.04*0.04+0.05*0.05*0.05);
    require(near(deposited,exact_volume,1.0e-14),"conservative particle-volume projection");
    require(near_vec(fluid_reaction+particle_total,{},1.0e-10),"two-way CFD/DEM momentum conservation");
    require(coupled.particle_force[0].x>0.0&&coupled.particle_force[1].x>0.0,"fluid accelerates slower particles");

    apply_unresolved_particle_forces(bodies,spheres,coupled.particle_force);
    require(near_vec(bodies[0].force,coupled.particle_force[0])&&near_vec(bodies[1].force,coupled.particle_force[1]),
            "unresolved forces apply to rigid bodies");


#if defined(CFD_HAS_SYCL)
    {
        sycl::queue q{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
        const auto gpu_mesh=cfd::fvm::make_cartesian_hexa_mesh(4U,1U,1U,1.0,0.2,0.2);
        cfd::fvm::ResidentPolyMeshSycl resident_mesh(gpu_mesh,q);
        ResidentDemSyclConfig dcfg;dcfg.capacity=4U;dcfg.domain_min={0.0,0.0,0.0};dcfg.domain_max={1.0,0.2,0.2};
        dcfg.cell_size=0.25;dcfg.dt=1.0e-3;dcfg.gravity={};
        ResidentDemSycl dem(q,dcfg);
        ResidentDemParticle rp;rp.global_id=7U;rp.position={0.37,0.1,0.1};rp.radius=0.02;rp.mass=1.0;
        dem.upload(std::span<const ResidentDemParticle>(&rp,1U));

        double* u=sycl::malloc_device<double>(3U*gpu_mesh.cell_count(),q);
        double* grad=sycl::malloc_device<double>(3U*gpu_mesh.cell_count(),q);
        double* vortd=sycl::malloc_device<double>(3U*gpu_mesh.cell_count(),q);
        double* accel=sycl::malloc_device<double>(3U*gpu_mesh.cell_count(),q);
        require(u&&grad&&vortd&&accel,"resident CFD/DEM test allocation");
        std::vector<cfd::fvm::Vec3> uh(gpu_mesh.cell_count(),{1.0,0.0,0.0});
        resident_mesh.upload_cell_vector(uh,u);resident_mesh.fill(grad,3U*gpu_mesh.cell_count(),0.0);
        resident_mesh.fill(vortd,3U*gpu_mesh.cell_count(),0.0);resident_mesh.fill(accel,3U*gpu_mesh.cell_count(),0.0);resident_mesh.wait();

        ResidentCfdDemSyclConfig ccfg;ccfg.fluid_density=1.2;ccfg.dynamic_viscosity=1.8e-5;
        ccfg.include_pressure_gradient=false;ccfg.include_saffman_lift=false;
        ResidentCfdDemSyclCoupler rc(resident_mesh,dem,ccfg);
        resident_mesh.reset_transfer_stats();dem.reset_transfer_stats();
        dem.begin_distributed_contact_step();
        rc.couple(u,grad,vortd,accel);
        dem.finish_distributed_contact_step();
        require(resident_mesh.transfer_stats().host_to_device_bytes==0U&&resident_mesh.transfer_stats().device_to_host_bytes==0U,
                "resident CFD/DEM hot loop has no FVM bulk host transfer");
        require(dem.transfer_stats().host_to_device_bytes==0U&&dem.transfer_stats().device_to_host_bytes==0U,
                "resident CFD/DEM hot loop has no DEM bulk host transfer");

        std::vector<double> alpha;rc.download_solid_volume_fraction(alpha);
        double vdep=0.0;for(std::size_t c=0;c<alpha.size();++c)vdep+=alpha[c]*gpu_mesh.cells()[c].volume;
        const double vp=(4.0/3.0)*std::numbers::pi*0.02*0.02*0.02;
        require(near(vdep,vp,2.0e-10),"resident conservative void-fraction projection");
        auto after=dem.download();require(after.size()==1U&&after[0].linear_velocity.x>0.0,"resident fluid drag accelerates DEM particle");
        std::vector<double> ah(3U*gpu_mesh.cell_count());q.memcpy(ah.data(),accel,ah.size()*sizeof(double)).wait_and_throw();
        double fluid_impulse_rate=0.0;for(std::size_t c=0;c<gpu_mesh.cell_count();++c)fluid_impulse_rate+=ccfg.fluid_density*ah[c]*gpu_mesh.cells()[c].volume;
        const double particle_force_x=rp.mass*(after[0].linear_velocity.x-rp.linear_velocity.x)/dcfg.dt;
        require(near(fluid_impulse_rate+particle_force_x,0.0,2.0e-4),"resident CFD/DEM action-reaction conservation");
        sycl::free(u,q);sycl::free(grad,q);sycl::free(vortd,q);sycl::free(accel,q);
    }
#endif

    std::cout << "v0.16.7 CFD/DEM coupling tests passed\n";
    return 0;
}
