#include "cfd/fem/mesh3d.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fem {
namespace {

double tet_det(const Point3&a,const Point3&b,const Point3&c,const Point3&d) {
    const double ax=b.x-a.x, ay=b.y-a.y, az=b.z-a.z;
    const double bx=c.x-a.x, by=c.y-a.y, bz=c.z-a.z;
    const double cx=d.x-a.x, cy=d.y-a.y, cz=d.z-a.z;
    return ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx);
}

} // namespace

void Mesh3D::validate() const {
    if(nodes.empty()||tetrahedra.empty()) throw std::invalid_argument("FEM 3-D mesh must contain nodes/elements");
    if(boundary_node.size()!=nodes.size()) throw std::invalid_argument("FEM 3-D boundary-node mask size mismatch");
    for(const auto&t:tetrahedra) {
        for(const auto n:t.node) if(n>=nodes.size()) throw std::out_of_range("FEM Tet4 node index out of range");
        const double d=tet_det(nodes[t.node[0]],nodes[t.node[1]],nodes[t.node[2]],nodes[t.node[3]]);
        if(!(std::abs(d)>1.0e-30)) throw std::invalid_argument("degenerate FEM Tet4");
    }
}

Mesh3D make_box_tet_mesh(std::size_t nx,std::size_t ny,std::size_t nz,double width,double height,double depth) {
    if(nx==0U||ny==0U||nz==0U||!(width>0.0)||!(height>0.0)||!(depth>0.0)) throw std::invalid_argument("invalid box Tet4 mesh dimensions");
    Mesh3D mesh;
    const std::size_t sx=nx+1U,sy=ny+1U;
    const auto id=[sx,sy](std::size_t i,std::size_t j,std::size_t k){return (k*sy+j)*sx+i;};
    mesh.nodes.resize(sx*sy*(nz+1U));
    mesh.boundary_node.resize(mesh.nodes.size(),0U);
    for(std::size_t k=0;k<=nz;++k) for(std::size_t j=0;j<=ny;++j) for(std::size_t i=0;i<=nx;++i) {
        const std::size_t n=id(i,j,k);
        mesh.nodes[n]={width*static_cast<double>(i)/static_cast<double>(nx),height*static_cast<double>(j)/static_cast<double>(ny),depth*static_cast<double>(k)/static_cast<double>(nz)};
        mesh.boundary_node[n]=static_cast<unsigned char>(i==0U||i==nx||j==0U||j==ny||k==0U||k==nz);
    }
    mesh.tetrahedra.reserve(6U*nx*ny*nz);
    for(std::size_t k=0;k<nz;++k) for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const auto n000=id(i,j,k), n100=id(i+1U,j,k), n010=id(i,j+1U,k), n110=id(i+1U,j+1U,k);
        const auto n001=id(i,j,k+1U), n101=id(i+1U,j,k+1U), n011=id(i,j+1U,k+1U), n111=id(i+1U,j+1U,k+1U);
        mesh.tetrahedra.push_back({{n000,n100,n110,n111}});
        mesh.tetrahedra.push_back({{n000,n110,n010,n111}});
        mesh.tetrahedra.push_back({{n000,n010,n011,n111}});
        mesh.tetrahedra.push_back({{n000,n011,n001,n111}});
        mesh.tetrahedra.push_back({{n000,n001,n101,n111}});
        mesh.tetrahedra.push_back({{n000,n101,n100,n111}});
    }
    mesh.validate();
    return mesh;
}

} // namespace cfd::fem
