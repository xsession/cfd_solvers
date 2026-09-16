#include "cfd/fvm/poly_mesh.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfd::fvm {

double magnitude(Vec3 a) noexcept { return std::sqrt(dot(a,a)); }

PolyMesh::PolyMesh(std::vector<Cell> cells,
                   std::vector<Face> faces,
                   std::vector<BoundaryPatch> patches)
    : cells_(std::move(cells)), faces_(std::move(faces)), patches_(std::move(patches)) {
    validate_and_build_adjacency();
}

void PolyMesh::validate_and_build_adjacency() {
    if (cells_.empty()) throw std::invalid_argument("PolyMesh requires at least one cell");
    cell_faces_.assign(cells_.size(), {});
    for (std::size_t c=0;c<cells_.size();++c) {
        if (!(cells_[c].volume > 0.0)) throw std::invalid_argument("PolyMesh cell volume must be positive");
    }
    for (std::size_t f=0;f<faces_.size();++f) {
        const auto& face=faces_[f];
        if (face.owner>=cells_.size()) throw std::invalid_argument("PolyMesh face owner out of range");
        if (!(magnitude(face.area)>0.0)) throw std::invalid_argument("PolyMesh face area must be non-zero");
        cell_faces_[face.owner].push_back(f);
        if (face.boundary()) {
            if (face.patch>=patches_.size()) throw std::invalid_argument("PolyMesh boundary patch out of range");
        } else {
            if (face.neighbour>=cells_.size() || face.neighbour==face.owner) {
                throw std::invalid_argument("PolyMesh face neighbour invalid");
            }
            if (face.patch!=invalid_patch) throw std::invalid_argument("PolyMesh internal face must not have a patch");
            cell_faces_[face.neighbour].push_back(f);
        }
    }
}

std::size_t PolyMesh::boundary_face_count() const noexcept {
    std::size_t n=0;
    for (const auto& f:faces_) if (f.boundary()) ++n;
    return n;
}

std::size_t PolyMesh::patch_index(std::string_view name) const {
    for (std::size_t i=0;i<patches_.size();++i) {
        if (patches_[i].name==name) return i;
    }
    throw std::out_of_range("PolyMesh patch not found: " + std::string(name));
}

PolyMesh make_cartesian_hexa_mesh(std::size_t nx,std::size_t ny,std::size_t nz,
                                  double lx,double ly,double lz) {
    if (nx==0||ny==0||nz==0||!(lx>0.0)||!(ly>0.0)||!(lz>0.0)) {
        throw std::invalid_argument("Cartesian mesh requires positive dimensions and lengths");
    }
    const double dx=lx/static_cast<double>(nx);
    const double dy=ly/static_cast<double>(ny);
    const double dz=lz/static_cast<double>(nz);
    const double volume=dx*dy*dz;
    auto idx=[=](std::size_t i,std::size_t j,std::size_t k){return (k*ny+j)*nx+i;};

    std::vector<Cell> cells;
    cells.reserve(nx*ny*nz);
    for (std::size_t k=0;k<nz;++k)
        for (std::size_t j=0;j<ny;++j)
            for (std::size_t i=0;i<nx;++i)
                cells.push_back({{(static_cast<double>(i)+0.5)*dx,
                                  (static_cast<double>(j)+0.5)*dy,
                                  (static_cast<double>(k)+0.5)*dz},volume});

    std::vector<BoundaryPatch> patches{{"left"},{"right"},{"bottom"},{"top"},{"front"},{"back"}};
    std::vector<Face> faces;
    faces.reserve((nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1));

    // X-normal faces.
    for (std::size_t k=0;k<nz;++k) for (std::size_t j=0;j<ny;++j) {
        faces.push_back({idx(0,j,k),invalid_cell,{0.0,(static_cast<double>(j)+0.5)*dy,(static_cast<double>(k)+0.5)*dz},{-dy*dz,0,0},0});
        for (std::size_t i=1;i<nx;++i) {
            faces.push_back({idx(i-1,j,k),idx(i,j,k),{static_cast<double>(i)*dx,(static_cast<double>(j)+0.5)*dy,(static_cast<double>(k)+0.5)*dz},{dy*dz,0,0},invalid_patch});
        }
        faces.push_back({idx(nx-1,j,k),invalid_cell,{lx,(static_cast<double>(j)+0.5)*dy,(static_cast<double>(k)+0.5)*dz},{dy*dz,0,0},1});
    }
    // Y-normal faces.
    for (std::size_t k=0;k<nz;++k) for (std::size_t i=0;i<nx;++i) {
        faces.push_back({idx(i,0,k),invalid_cell,{(static_cast<double>(i)+0.5)*dx,0.0,(static_cast<double>(k)+0.5)*dz},{0,-dx*dz,0},2});
        for (std::size_t j=1;j<ny;++j) {
            faces.push_back({idx(i,j-1,k),idx(i,j,k),{(static_cast<double>(i)+0.5)*dx,static_cast<double>(j)*dy,(static_cast<double>(k)+0.5)*dz},{0,dx*dz,0},invalid_patch});
        }
        faces.push_back({idx(i,ny-1,k),invalid_cell,{(static_cast<double>(i)+0.5)*dx,ly,(static_cast<double>(k)+0.5)*dz},{0,dx*dz,0},3});
    }
    // Z-normal faces.
    for (std::size_t j=0;j<ny;++j) for (std::size_t i=0;i<nx;++i) {
        faces.push_back({idx(i,j,0),invalid_cell,{(static_cast<double>(i)+0.5)*dx,(static_cast<double>(j)+0.5)*dy,0.0},{0,0,-dx*dy},4});
        for (std::size_t k=1;k<nz;++k) {
            faces.push_back({idx(i,j,k-1),idx(i,j,k),{(static_cast<double>(i)+0.5)*dx,(static_cast<double>(j)+0.5)*dy,static_cast<double>(k)*dz},{0,0,dx*dy},invalid_patch});
        }
        faces.push_back({idx(i,j,nz-1),invalid_cell,{(static_cast<double>(i)+0.5)*dx,(static_cast<double>(j)+0.5)*dy,lz},{0,0,dx*dy},5});
    }

    return PolyMesh(std::move(cells),std::move(faces),std::move(patches));
}

PolyMesh make_sheared_cartesian_hexa_mesh(std::size_t nx,std::size_t ny,std::size_t nz,
                                          double lx,double ly,double lz,double shear_xy) {
    auto base=make_cartesian_hexa_mesh(nx,ny,nz,lx,ly,lz);
    std::vector<Cell> cells=base.cells();
    std::vector<Face> faces=base.faces();
    std::vector<BoundaryPatch> patches=base.patches();
    auto transform_point=[=](Vec3 p) { return Vec3{p.x+shear_xy*p.y,p.y,p.z}; };
    auto transform_area=[=](Vec3 a) { return Vec3{a.x,a.y-shear_xy*a.x,a.z}; };
    for (auto& cell:cells) cell.center=transform_point(cell.center);
    for (auto& face:faces) {
        face.center=transform_point(face.center);
        face.area=transform_area(face.area);
    }
    return PolyMesh(std::move(cells),std::move(faces),std::move(patches));
}

} // namespace cfd::fvm
