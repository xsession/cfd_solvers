#include "cfd/fvm/operators.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/fvm/pressure_velocity.hpp"

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
    gauss_gradient_scalar_into(mesh,values,result,boundary);
    return result;
}

void gauss_gradient_scalar_into(const PolyMesh& mesh,std::span<const double> values,
                                std::span<Vec3> result,std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary);
    if(result.size()!=mesh.cell_count()) throw std::invalid_argument("gradient output size mismatch");
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        Vec3 sum{};
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const auto& f=mesh.faces()[fi];
            const double phi_f=f.boundary() && !boundary.empty() ? boundary[fi] : interpolated_scalar(mesh,f,values);
            sum+=outward_area(f,cell)*phi_f;
        }
        result[cell]=sum/mesh.cells()[cell].volume;
    });
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

namespace {

Vec3 solve_weighted_gradient_system(double a00,double a01,double a02,
                                    double a11,double a12,double a22,
                                    double b0,double b1,double b2) {
    double a[3][4]={{a00,a01,a02,b0},{a01,a11,a12,b1},{a02,a12,a22,b2}};
    const double trace=std::abs(a00)+std::abs(a11)+std::abs(a22);
    const double reg=std::max(trace*1.0e-14,1.0e-30);
    for(int d=0;d<3;++d)a[d][d]+=reg;
    for(int col=0;col<3;++col){
        int pivot=col;
        for(int r=col+1;r<3;++r)if(std::abs(a[r][col])>std::abs(a[pivot][col]))pivot=r;
        if(std::abs(a[pivot][col])<1.0e-30)continue;
        if(pivot!=col)for(int c=col;c<4;++c)std::swap(a[pivot][c],a[col][c]);
        const double inv=1.0/a[col][col]; for(int c=col;c<4;++c)a[col][c]*=inv;
        for(int r=0;r<3;++r)if(r!=col){const double f=a[r][col];for(int c=col;c<4;++c)a[r][c]-=f*a[col][c];}
    }
    return {a[0][3],a[1][3],a[2][3]};
}

} // namespace

std::vector<Vec3> least_squares_gradient_scalar(const PolyMesh& mesh,
                                                 std::span<const double> values,
                                                 std::span<const double> boundary) {
    check_scalar_sizes(mesh, values, boundary);
    std::vector<Vec3> result(mesh.cell_count());
    least_squares_gradient_scalar_into(mesh,values,result,boundary);
    return result;
}

void least_squares_gradient_scalar_into(const PolyMesh& mesh,std::span<const double> values,
                                        std::span<Vec3> result,std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary);
    if(result.size()!=mesh.cell_count()) throw std::invalid_argument("gradient output size mismatch");
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double a00=0.0,a01=0.0,a02=0.0,a11=0.0,a12=0.0,a22=0.0,b0=0.0,b1=0.0,b2=0.0;
        const auto c0=mesh.cells()[cell].center; const double v0=values[cell];
        for(const std::size_t fi:mesh.cell_faces()[cell]){
            const auto& f=mesh.faces()[fi]; Vec3 d{}; double dv=0.0;
            if(!f.boundary()){
                const std::size_t other=f.owner==cell?f.neighbour:f.owner;
                d=mesh.cells()[other].center-c0; dv=values[other]-v0;
            }else if(!boundary.empty()){
                d=f.center-c0; dv=boundary[fi]-v0;
            }else continue;
            const double d2=dot(d,d); if(!(d2>0.0))continue; const double w=1.0/d2;
            a00+=w*d.x*d.x; a01+=w*d.x*d.y; a02+=w*d.x*d.z;
            a11+=w*d.y*d.y; a12+=w*d.y*d.z; a22+=w*d.z*d.z;
            b0+=w*d.x*dv; b1+=w*d.y*dv; b2+=w*d.z*dv;
        }
        result[cell]=solve_weighted_gradient_system(a00,a01,a02,a11,a12,a22,b0,b1,b2);
    });
}

std::vector<double> corrected_face_normal_gradient_scalar(const PolyMesh& mesh,
                                                            std::span<const double> values,
                                                            std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary);
    const auto grad=least_squares_gradient_scalar(mesh,values,boundary);
    std::vector<double> result(mesh.face_count(),0.0);
    cfd::core::parallel_for(mesh.face_count(),[&](std::size_t fi){
        const auto& f=mesh.faces()[fi]; const auto geom=decompose_face_area(mesh,fi); const double area=magnitude(f.area);
        double delta=0.0; Vec3 gf=grad[f.owner];
        if(!f.boundary()){
            delta=values[f.neighbour]-values[f.owner]; gf=(grad[f.owner]+grad[f.neighbour])*0.5;
        }else if(!boundary.empty()) delta=boundary[fi]-values[f.owner];
        else return;
        const double flux=geom.orthogonal_metric*delta+dot(gf,geom.nonorthogonal_area);
        result[fi]=flux/area;
    });
    return result;
}

std::vector<double> corrected_laplacian_scalar(const PolyMesh& mesh,
                                                std::span<const double> values,
                                                double diffusivity,
                                                std::span<const double> boundary) {
    check_scalar_sizes(mesh,values,boundary); if(!(diffusivity>0.0))throw std::invalid_argument("FVM diffusivity must be positive");
    const auto sn=corrected_face_normal_gradient_scalar(mesh,values,boundary);
    std::vector<double> result(mesh.cell_count(),0.0);
    cfd::core::parallel_for(mesh.cell_count(),[&](std::size_t cell){
        double sum=0.0; for(const std::size_t fi:mesh.cell_faces()[cell]){
            const auto& f=mesh.faces()[fi]; if(f.boundary()&&boundary.empty())continue;
            const double oriented=(f.owner==cell?1.0:-1.0); sum+=diffusivity*oriented*magnitude(f.area)*sn[fi];
        } result[cell]=sum/mesh.cells()[cell].volume;
    });
    return result;
}


} // namespace cfd::fvm
