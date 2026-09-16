#include "cfd/fvm/operators.hpp"
#include "cfd/core/parallel.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {

void check_scalar_sizes(const PolyMesh& mesh,std::span<const double> cells,std::span<const double> boundary) {
    if (cells.size()!=mesh.cell_count()) throw std::invalid_argument("FVM scalar field size mismatch");
    if (!boundary.empty() && boundary.size()!=mesh.face_count()) throw std::invalid_argument("FVM boundary scalar size mismatch");
}
void check_vector_sizes(const PolyMesh& mesh,std::span<const Vec3> cells,std::span<const Vec3> boundary) {
    if (cells.size()!=mesh.cell_count()) throw std::invalid_argument("FVM vector field size mismatch");
    if (!boundary.empty() && boundary.size()!=mesh.face_count()) throw std::invalid_argument("FVM boundary vector size mismatch");
}

double interpolated_scalar(const PolyMesh& mesh,const Face& f,std::span<const double> values) {
    if (f.boundary()) return values[f.owner];
    const auto& co=mesh.cells()[f.owner].center;
    const auto& cn=mesh.cells()[f.neighbour].center;
    const double dof=magnitude(f.center-co);
    const double dnf=magnitude(cn-f.center);
    const double d=dof+dnf;
    if (!(d>0.0)) throw std::runtime_error("degenerate FVM face interpolation distance");
    return (dnf*values[f.owner]+dof*values[f.neighbour])/d;
}
Vec3 interpolated_vector(const PolyMesh& mesh,const Face& f,std::span<const Vec3> values) {
    if (f.boundary()) return values[f.owner];
    const auto& co=mesh.cells()[f.owner].center;
    const auto& cn=mesh.cells()[f.neighbour].center;
    const double dof=magnitude(f.center-co);
    const double dnf=magnitude(cn-f.center);
    const double d=dof+dnf;
    if (!(d>0.0)) throw std::runtime_error("degenerate FVM face interpolation distance");
    return (values[f.owner]*dnf+values[f.neighbour]*dof)/d;
}

[[nodiscard]] Vec3 outward_area(const Face& face,std::size_t cell) noexcept {
    return face.owner==cell ? face.area : face.area*(-1.0);
}

} // namespace

std::vector<Vec3> gauss_gradient_scalar(const PolyMesh& mesh,
                                         std::span<const double> values,
                                         std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary);
    std::vector<Vec3> result(mesh.cell_count());
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        Vec3 sum{};
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const auto& f=mesh.faces()[fi];
            const double phi_f=f.boundary() && !boundary.empty() ? boundary[fi] : interpolated_scalar(mesh,f,values);
            sum+=outward_area(f,cell)*phi_f;
        }
        result[cell]=sum/mesh.cells()[cell].volume;
    });
    return result;
}

std::vector<double> gauss_divergence_vector(const PolyMesh& mesh,
                                             std::span<const Vec3> values,
                                             std::span<const Vec3> boundary) {
    check_vector_sizes(mesh,values,boundary);
    std::vector<double> result(mesh.cell_count(),0.0);
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double sum=0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const auto& f=mesh.faces()[fi];
            const Vec3 u_f=f.boundary() && !boundary.empty() ? boundary[fi] : interpolated_vector(mesh,f,values);
            sum+=dot(u_f,outward_area(f,cell));
        }
        result[cell]=sum/mesh.cells()[cell].volume;
    });
    return result;
}

std::vector<double> orthogonal_laplacian_scalar(const PolyMesh& mesh,
                                                 std::span<const double> values,
                                                 double diffusivity,
                                                 std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary);
    if (!(diffusivity>0.0)) throw std::invalid_argument("FVM diffusivity must be positive");
    std::vector<double> result(mesh.cell_count(),0.0);
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double sum=0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const auto& f=mesh.faces()[fi];
            const Vec3 sf=outward_area(f,cell);
            const double area=magnitude(sf);
            const Vec3 n=sf/area;
            if (!f.boundary()) {
                const std::size_t other=f.owner==cell ? f.neighbour : f.owner;
                const double distance=std::abs(dot(mesh.cells()[other].center-mesh.cells()[cell].center,n));
                if (!(distance>0.0)) throw std::runtime_error("degenerate FVM internal face distance");
                sum+=diffusivity*area*(values[other]-values[cell])/distance;
            } else if (!boundary.empty()) {
                const double distance=std::abs(dot(f.center-mesh.cells()[cell].center,n));
                if (!(distance>0.0)) throw std::runtime_error("degenerate FVM boundary face distance");
                sum+=diffusivity*area*(boundary[fi]-values[cell])/distance;
            }
        }
        result[cell]=sum/mesh.cells()[cell].volume;
    });
    return result;
}

} // namespace cfd::fvm
