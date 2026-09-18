#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/multibody/resident_dem_sycl.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace cfd::multibody;

namespace {
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}

void test_cpu_slab_reference(){
    DemSlabDecomposition d;d.x_min=0.0;d.x_max=4.0;d.ranks=4;d.validate();
    require(d.owner_rank(0.1)==0 && d.owner_rank(1.1)==1 && d.owner_rank(3.9)==3,"slab owner mapping");
    require(std::abs(d.slab_min(2)-2.0)<1.0e-12 && std::abs(d.slab_max(2)-3.0)<1.0e-12,"slab bounds");
}

#if defined(CFD_HAS_SYCL)
void test_device_pack_restore_stable_ids(){
    static_assert(std::is_trivially_copyable_v<ResidentDemPackedParticle>);
    static_assert(std::is_trivially_copyable_v<ResidentDemPackedHistory>);
    sycl::queue queue{sycl::device{sycl::default_selector_v},sycl::property::queue::in_order{}};
    ResidentDemSyclConfig cfg;cfg.capacity=8U;cfg.domain_min={0,0,0};cfg.domain_max={4,1,1};cfg.cell_size=0.25;cfg.gravity={0,0,0};cfg.history_slots_per_particle=4U;
    ResidentDemSycl dem(queue,cfg);
    std::vector<ResidentDemParticle> p(3);
    p[0].global_id=101;p[0].owner_rank=0;p[0].position={0.25,0.5,0.5};p[0].radius=0.05;p[0].mass=1.0;
    p[1].global_id=202;p[1].owner_rank=0;p[1].position={1.25,0.5,0.5};p[1].radius=0.05;p[1].mass=2.0;
    p[2].global_id=303;p[2].owner_rank=0;p[2].position={0.75,0.5,0.5};p[2].radius=0.05;p[2].mass=3.0;
    dem.upload(p);

    auto* index=sycl::malloc_device<std::uint32_t>(1U,queue);
    auto* packed=sycl::malloc_device<ResidentDemPackedParticle>(1U,queue);
    auto* history=sycl::malloc_device<ResidentDemPackedHistory>(cfg.history_slots_per_particle,queue);
    const std::uint32_t host_index=1U;queue.memcpy(index,&host_index,sizeof(host_index)).wait();
    dem.pack_indices_device(index,1U,packed,history);queue.wait_and_throw();
    dem.retain_x_slab_device(0.0,1.0,0);
    require(dem.particle_count()==2U,"device slab compaction should remove migrant");
    dem.append_packed_device(packed,history,1U,1);queue.wait_and_throw();
    const auto out=dem.download();
    require(out.size()==3U,"packed migrant append size");
    bool found=false;
    for(const auto& q:out)if(q.global_id==202U){found=true;require(q.owner_rank==1,"owner override after device migration");require(std::abs(q.mass-2.0)<1.0e-5,"migrant mass preserved");}
    require(found,"stable global ID survives device pack/compact/append");
    sycl::free(index,queue);sycl::free(packed,queue);sycl::free(history,queue);
}
#endif
}

int main(){
    try{
        test_cpu_slab_reference();
#if defined(CFD_HAS_SYCL)
        test_device_pack_restore_stable_ids();
#endif
    }catch(const std::exception& e){std::cerr<<"v0.16.5 GPU-aware DEM exchange regression failed: "<<e.what()<<'\n';return 1;}
    std::cout<<"v0.16.5 GPU-aware DEM exchange regression passed\n";return 0;
}
