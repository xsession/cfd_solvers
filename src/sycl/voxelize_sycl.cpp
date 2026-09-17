#include "cfd/solvers/lbm/voxelize_sycl.hpp"

#if defined(CFD_HAS_SYCL)
#include <array>
#include <stdexcept>
#include <vector>

namespace cfd::lbm {
namespace {
struct TriangleF { float ax,ay,az,bx,by,bz,cx,cy,cz; };
}

std::vector<unsigned char> voxelize_surface_sycl(const cfd::io::TriangleSurface& surface,
                                                  const VoxelGrid3D& grid,
                                                  sycl::queue& queue) {
    if (surface.triangles.empty() || grid.nx == 0 || grid.ny == 0 || grid.nz == 0 ||
        !(grid.maximum.x > grid.minimum.x) || !(grid.maximum.y > grid.minimum.y) ||
        !(grid.maximum.z > grid.minimum.z)) {
        throw std::invalid_argument("invalid SYCL voxelization inputs");
    }
    std::vector<TriangleF> host_triangles;
    host_triangles.reserve(surface.triangles.size());
    for (const auto& t : surface.triangles) {
        host_triangles.push_back({static_cast<float>(t.a.x),static_cast<float>(t.a.y),static_cast<float>(t.a.z),
                                  static_cast<float>(t.b.x),static_cast<float>(t.b.y),static_cast<float>(t.b.z),
                                  static_cast<float>(t.c.x),static_cast<float>(t.c.y),static_cast<float>(t.c.z)});
    }
    const std::size_t cells = grid.nx * grid.ny * grid.nz;
    auto* triangles = sycl::malloc_device<TriangleF>(host_triangles.size(), queue);
    auto* mask = sycl::malloc_shared<unsigned char>(cells, queue);
    if (!triangles || !mask) {
        if (triangles) sycl::free(triangles, queue);
        if (mask) sycl::free(mask, queue);
        throw std::runtime_error("SYCL voxelization allocation failed");
    }
    queue.memcpy(triangles, host_triangles.data(), host_triangles.size()*sizeof(TriangleF)).wait();
    const float minx=static_cast<float>(grid.minimum.x),miny=static_cast<float>(grid.minimum.y),minz=static_cast<float>(grid.minimum.z);
    const float dx=static_cast<float>((grid.maximum.x-grid.minimum.x)/grid.nx);
    const float dy=static_cast<float>((grid.maximum.y-grid.minimum.y)/grid.ny);
    const float dz=static_cast<float>((grid.maximum.z-grid.minimum.z)/grid.nz);
    const std::size_t nx=grid.nx,ny=grid.ny,nz=grid.nz,triangle_count=host_triangles.size();
    queue.parallel_for(sycl::range<1>(cells), [=](sycl::id<1> item) {
        const std::size_t n=item[0], i=n%nx, j=(n/nx)%ny, k=n/(nx*ny);
        const float ox=minx+(static_cast<float>(i)+0.5F)*dx;
        const float oy=miny+(static_cast<float>(j)+0.5F)*dy+1.0e-6F*dx;
        const float oz=minz+(static_cast<float>(k)+0.5F)*dz+2.0e-6F*dx;
        unsigned hits=0;
        for(std::size_t ti=0;ti<triangle_count;++ti){
            const TriangleF t=triangles[ti];
            const float e1x=t.bx-t.ax,e1y=t.by-t.ay,e1z=t.bz-t.az;
            const float e2x=t.cx-t.ax,e2y=t.cy-t.ay,e2z=t.cz-t.az;
            // ray d=(1,0,0); h=d x e2=(0,-e2z,e2y)
            const float hx=0.0F,hy=-e2z,hz=e2y;
            const float det=e1x*hx+e1y*hy+e1z*hz;
            const float adet=det<0.0F?-det:det;
            if(adet<1.0e-8F)continue;
            const float inv=1.0F/det,sx=ox-t.ax,sy=oy-t.ay,sz=oz-t.az;
            const float u=inv*(sx*hx+sy*hy+sz*hz);if(u<0.0F||u>1.0F)continue;
            const float qx=sy*e1z-sz*e1y,qy=sz*e1x-sx*e1z,qz=sx*e1y-sy*e1x;
            const float v=inv*qx;if(v<0.0F||u+v>1.0F)continue;
            const float dist=inv*(e2x*qx+e2y*qy+e2z*qz);if(dist>1.0e-7F)++hits;
        }
        mask[n]=static_cast<unsigned char>(hits&1U);
    }).wait();
    std::vector<unsigned char> result(mask, mask+cells);
    sycl::free(mask, queue); sycl::free(triangles, queue);
    return result;
}
} // namespace cfd::lbm
#endif
