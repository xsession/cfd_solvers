#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/multibody/mpi_resident_dem_sycl.hpp"
#include "cfd/multibody/resident_dem_sycl.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace cfd::multibody;

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void test_cpu_reference_contract(){
    DemSlabDecomposition d;d.x_min=0.0;d.x_max=2.0;d.ranks=2;d.validate();
    require(d.owner_rank(0.5)==0&&d.owner_rank(1.5)==1,"distributed resident DEM slab reference");
#if defined(CFD_HAS_SYCL)
    static_assert(std::is_trivially_copyable_v<ResidentDemPackedForce>);
    static_assert(std::is_trivially_copyable_v<ResidentDemPackedBond>);
#endif
}

#if defined(CFD_HAS_SYCL)
ResidentDemSyclConfig base_config(){
    ResidentDemSyclConfig cfg;cfg.capacity=8U;cfg.domain_min={0,0,0};cfg.domain_max={2,1,1};cfg.cell_size=0.25;cfg.dt=1.0e-3;cfg.gravity={0,0,0};cfg.normal_stiffness=2.0e4;cfg.normal_damping=0.0;cfg.tangential_stiffness=5.0e3;cfg.tangential_damping=0.0;cfg.friction=1.0;cfg.rolling_resistance=0.0;cfg.history_slots_per_particle=8U;return cfg;
}

void test_local_ghost_action_reaction_and_history(){
    sycl::queue q0{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
    sycl::queue q1{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
    auto cfg=base_config();ResidentDemSycl left(q0,cfg),right(q1,cfg);
    ResidentDemParticle a;a.global_id=10;a.owner_rank=0;a.position={0.95,0.5,0.5};a.linear_velocity={0,0.10,0};a.radius=0.10;a.mass=1.0;
    ResidentDemParticle b;b.global_id=20;b.owner_rank=1;b.position={1.05,0.5,0.5};b.radius=0.10;b.mass=1.0;
    left.upload(std::span<const ResidentDemParticle>(&a,1));right.upload(std::span<const ResidentDemParticle>(&b,1));

    auto* boundary0=sycl::malloc_device<std::uint8_t>(1U,q0);auto* boundary1=sycl::malloc_device<std::uint8_t>(1U,q1);
    const std::uint8_t one=1U;q0.memcpy(boundary0,&one,1U);q1.memcpy(boundary1,&one,1U);
    auto* ghost0=sycl::malloc_device<ResidentDemPackedParticle>(1U,q0);auto* ghost1=sycl::malloc_device<ResidentDemPackedParticle>(1U,q1);
    auto* remote0=sycl::malloc_device<ResidentDemPackedForce>(8U,q0);auto* remote1=sycl::malloc_device<ResidentDemPackedForce>(8U,q1);
    auto* count0=sycl::malloc_shared<std::uint32_t>(1U,q0);auto* count1=sycl::malloc_shared<std::uint32_t>(1U,q1);
    ResidentDemPackedParticle gb;gb.global_id=20;gb.owner_rank=1;gb.x=1.05F;gb.y=0.5F;gb.z=0.5F;gb.vx=0;gb.vy=0;gb.vz=0;gb.wx=gb.wy=gb.wz=0;gb.radius=0.10F;gb.inverse_mass=1.0F;gb.inverse_inertia=250.0F;
    ResidentDemPackedParticle ga=gb;ga.global_id=10;ga.owner_rank=0;ga.x=0.95F;ga.vy=0.10F;
    q0.memcpy(ghost0,&gb,sizeof(gb));q1.memcpy(ghost1,&ga,sizeof(ga));q0.wait_and_throw();q1.wait_and_throw();

    *count0=0U;*count1=0U;left.begin_distributed_contact_step();right.begin_distributed_contact_step();
    left.apply_owned_contacts_interior_device(boundary0);right.apply_owned_contacts_interior_device(boundary1);
    left.apply_owned_contacts_boundary_device(boundary0);right.apply_owned_contacts_boundary_device(boundary1);
    left.apply_ghost_contacts_device(ghost0,1U,boundary0,0,remote0,count0,8U);
    right.apply_ghost_contacts_device(ghost1,1U,boundary1,1,remote1,count1,8U);
    q0.wait_and_throw();q1.wait_and_throw();
    require(*count0==1U,"lower-ID rank must own exactly one cross-rank contact");
    require(*count1==0U,"higher-ID rank must not duplicate cross-rank contact");
    ResidentDemPackedForce reaction{};q0.memcpy(&reaction,remote0,sizeof(reaction)).wait_and_throw();
    require(reaction.global_id==20U&&reaction.owner_rank==1,"remote reaction target");
    q1.memcpy(remote1,&reaction,sizeof(reaction)).wait_and_throw();right.apply_remote_forces_device(remote1,1U);
    left.finish_distributed_contact_step();right.finish_distributed_contact_step();q0.wait_and_throw();q1.wait_and_throw();
    const auto la=left.download().front(),rb=right.download().front();
    require(std::abs(la.linear_velocity.x+rb.linear_velocity.x)<2.0e-4,"cross-rank normal momentum conservation");
    require(std::abs(la.linear_velocity.y+rb.linear_velocity.y-0.10)<2.0e-4,"cross-rank tangential momentum conservation");

    // Repeat with refreshed ghost state. Tangential spring history remains on
    // the lower-ID owner and should continue producing a tangential reaction.
    ResidentDemPackedParticle gb2=gb;gb2.x=static_cast<float>(rb.position.x);gb2.y=static_cast<float>(rb.position.y);gb2.z=static_cast<float>(rb.position.z);gb2.vx=static_cast<float>(rb.linear_velocity.x);gb2.vy=static_cast<float>(rb.linear_velocity.y);gb2.vz=static_cast<float>(rb.linear_velocity.z);
    q0.memcpy(ghost0,&gb2,sizeof(gb2)).wait_and_throw();*count0=0U;left.begin_distributed_contact_step();left.apply_owned_contacts_boundary_device(boundary0);left.apply_ghost_contacts_device(ghost0,1U,boundary0,0,remote0,count0,8U);q0.wait_and_throw();require(*count0==1U,"persistent ghost contact second step");ResidentDemPackedForce reaction2{};q0.memcpy(&reaction2,remote0,sizeof(reaction2)).wait_and_throw();require(std::abs(reaction2.fy)>1.0e-5F,"persistent Mindlin history/tangential force");left.finish_distributed_contact_step();q0.wait_and_throw();

    sycl::free(boundary0,q0);sycl::free(boundary1,q1);sycl::free(ghost0,q0);sycl::free(ghost1,q1);sycl::free(remote0,q0);sycl::free(remote1,q1);sycl::free(count0,q0);sycl::free(count1,q1);
}

void test_cross_rank_bond_device_history(){
    sycl::queue queue{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};auto cfg=base_config();ResidentDemSycl dem(queue,cfg);
    ResidentDemParticle a;a.global_id=100;a.owner_rank=0;a.position={0.9,0.5,0.5};a.radius=0.05;a.mass=1.0;dem.upload(std::span<const ResidentDemParticle>(&a,1));
    auto* ghost=sycl::malloc_device<ResidentDemPackedParticle>(1U,queue);ResidentDemPackedParticle b{};b.global_id=200;b.owner_rank=1;b.x=1.2F;b.y=0.5F;b.z=0.5F;b.radius=0.05F;b.inverse_mass=1.0F;b.inverse_inertia=1000.0F;queue.memcpy(ghost,&b,sizeof(b));
    auto* remote=sycl::malloc_device<ResidentDemPackedForce>(4U,queue);auto* count=sycl::malloc_shared<std::uint32_t>(1U,queue);auto* bond=sycl::malloc_device<ResidentDemPackedBond>(1U,queue);
    ResidentDemPackedBond host_bond{};host_bond.particle_a=100;host_bond.particle_b=200;host_bond.rest_length=0.2F;host_bond.normal_stiffness=1000.0F;host_bond.shear_stiffness=500.0F;host_bond.normal_damping=0;host_bond.shear_damping=0;host_bond.tensile_failure_force=10000.0F;host_bond.shear_failure_force=10000.0F;host_bond.damage_onset_ratio=0.7F;queue.memcpy(bond,&host_bond,sizeof(host_bond)).wait_and_throw();
    *count=0U;dem.begin_distributed_contact_step();dem.apply_bonds_device(bond,1U,ghost,1U,0,remote,count,4U);queue.wait_and_throw();require(*count==1U,"cross-rank bond must emit one reverse force");ResidentDemPackedForce rf{};ResidentDemPackedBond bout{};queue.memcpy(&rf,remote,sizeof(rf));queue.memcpy(&bout,bond,sizeof(bout)).wait_and_throw();require(rf.global_id==200U&&rf.fx<0.0F,"bond reverse force direction");require(bout.broken==0U&&bout.damage<1.0F,"bond history remains intact below strength");dem.finish_distributed_contact_step();queue.wait_and_throw();
    sycl::free(ghost,queue);sycl::free(remote,queue);sycl::free(count,queue);sycl::free(bond,queue);
}
#endif

#if defined(CFD_HAS_SYCL) && defined(CFD_HAS_MPI)
void test_one_rank_mpi_resident_step(){
    sycl::queue queue{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};auto cfg=base_config();ResidentDemSycl dem(queue,cfg);
    std::vector<ResidentDemParticle> p(2);p[0].global_id=1;p[0].owner_rank=0;p[0].position={0.45,0.5,0.5};p[0].radius=0.1;p[0].mass=1.0;p[1].global_id=2;p[1].owner_rank=0;p[1].position={0.55,0.5,0.5};p[1].radius=0.1;p[1].mass=1.0;dem.upload(p);
    ResidentDemMpiExchangeConfig ec;ec.decomposition.x_min=0.0;ec.decomposition.x_max=2.0;ec.mode=ResidentDemMpiMode::staged_host;ec.bond_capacity=4U;MpiResidentDemSyclExchange exchange(queue,ec,MPI_COMM_WORLD);
    ResidentDemPackedBond bond{};bond.particle_a=1;bond.particle_b=2;bond.rest_length=0.1F;bond.tensile_failure_force=1.0e6F;bond.shear_failure_force=1.0e6F;exchange.upload_bonds(std::span<const ResidentDemPackedBond>(&bond,1));
    const auto stats=exchange.step(dem,0.15);require(stats.nonblocking_exchanges>=1U,"resident MPI step must use nonblocking halo exchange");require(stats.resident_bonds==1U,"resident bond retained");require(exchange.bond_count()==1U,"resident bond count");require(dem.particle_count()==2U,"one-rank resident MPI step particle count");
}
#endif
}

int main(){
    try{
        test_cpu_reference_contract();
#if defined(CFD_HAS_SYCL)
        test_local_ghost_action_reaction_and_history();
        test_cross_rank_bond_device_history();
#endif
#if defined(CFD_HAS_SYCL) && defined(CFD_HAS_MPI)
        test_one_rank_mpi_resident_step();
#endif
    }catch(const std::exception& e){std::cerr<<"v0.16.6 distributed resident DEM regression failed: "<<e.what()<<'\n';return 1;}
    std::cout<<"v0.16.6 distributed resident DEM regression passed\n";return 0;
}
