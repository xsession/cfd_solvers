#include "cfd/core/device_residency.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/fdtd/resident_sycl.hpp"
#endif

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}

void test_residency_marker(){
    cfd::core::DeviceTransferStats stats;
    require(stats.host_transfer_bytes()==0U,"FDTD residency contract begins without host traffic");
}

#if defined(CFD_HAS_SYCL)
void test_resident_fdtd_chain(){
    cfd::fdtd::ResidentMaxwell3DConfig c;
    c.nx=18U;c.ny=16U;c.nz=14U;c.courant=0.55F;c.cpml_cells=3U;c.boundary=cfd::fdtd::ResidentBoundary3D::cpml;
    cfd::fdtd::ResidentMaxwell3DSycl solver(c);
    solver.initialize_gaussian_ez(0.2F,0.14F);
    solver.set_material_box(6U,12U,5U,11U,4U,10U,2.0F,2.5F,3.0F,1.0F,1.0e-5F);
    auto* probe=sycl::malloc_device<float>(8U,solver.queue());
    auto* halo=sycl::malloc_device<float>(solver.halo_value_count(cfd::fdtd::HaloFace3D::x_minus),solver.queue());
    require(probe!=nullptr&&halo!=nullptr,"resident FDTD scratch allocations succeed");
    solver.reset_transfer_stats();
    for(std::size_t s=0;s<4U;++s){
        solver.add_soft_ez_source(c.nx/2U,c.ny/2U,c.nz/2U,0.01F);
        solver.step();
        solver.capture_ez_to_device(c.nx/2U,c.ny/2U,c.nz/2U,probe,s);
    }
    solver.pack_halo(cfd::fdtd::HaloFace3D::x_minus,halo);
    solver.unpack_halo(cfd::fdtd::HaloFace3D::x_minus,halo);
    solver.wait();
    require(solver.transfer_stats().host_transfer_bytes()==0U,
            "resident FDTD timestep/probe/halo chain performs no bulk host transfer");
    require(solver.resident_bytes()>6U*solver.cells()*sizeof(float),"resident FDTD accounts for materials and CPML state");
    sycl::free(halo,solver.queue());sycl::free(probe,solver.queue());
}
#endif
}

int main(){
    try{
        test_residency_marker();
#if defined(CFD_HAS_SYCL)
        test_resident_fdtd_chain();
#endif
        std::cout<<"v0.15.7 GPU FDTD residency tests passed\n";
        return EXIT_SUCCESS;
    }catch(const std::exception& error){
        std::cerr<<"v0.15.7 GPU FDTD residency test failure: "<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
}
