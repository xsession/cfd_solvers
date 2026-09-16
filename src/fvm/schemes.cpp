#include "cfd/fvm/schemes.hpp"
#include "cfd/fvm/operators.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {

void check_flux(const PolyMesh& mesh, std::span<const double> flux) {
    if (flux.size() != mesh.face_count()) {
        throw std::invalid_argument("FVM face-flux size mismatch");
    }
}

template<class T>
void check_field(const PolyMesh& mesh,
                 std::span<const T> cells,
                 std::span<const T> boundary) {
    if (cells.size() != mesh.cell_count()) {
        throw std::invalid_argument("FVM cell-field size mismatch");
    }
    if (!boundary.empty() && boundary.size() != mesh.face_count()) {
        throw std::invalid_argument("FVM boundary-field size mismatch");
    }
}

template<class T>
T linear_face_value(const PolyMesh& mesh, const Face& face, std::span<const T> values) {
    if (face.boundary()) return values[face.owner];
    const Vec3& co = mesh.cells()[face.owner].center;
    const Vec3& cn = mesh.cells()[face.neighbour].center;
    const double dof = magnitude(face.center - co);
    const double dnf = magnitude(cn - face.center);
    const double total = dof + dnf;
    if (!(total > 0.0)) {
        throw std::runtime_error("degenerate FVM face interpolation distance");
    }
    return values[face.owner] * (dnf / total) + values[face.neighbour] * (dof / total);
}


void bounded_linear_scalar_faces(const PolyMesh& mesh,
                                                               std::span<const double> cells,
                                                               std::span<const double> flux,
                                                               std::span<const double> boundary,
                                                               std::span<double> faces,
                                                               ScalarSchemeWorkspace& workspace) {
    check_field(mesh, cells, boundary);
    check_flux(mesh, flux);
    auto& gradient=workspace.gradient;
    auto& minimum=workspace.minimum;
    auto& maximum=workspace.maximum;
    gauss_gradient_scalar_into(mesh,cells,gradient,boundary);

    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double lo = cells[cell];
        double hi = cells[cell];
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const Face& face = mesh.faces()[fi];
            double neighbour_value = cells[cell];
            if (face.boundary()) {
                if (!boundary.empty()) neighbour_value = boundary[fi];
            } else {
                const std::size_t other = face.owner == cell ? face.neighbour : face.owner;
                neighbour_value = cells[other];
            }
            lo = std::min(lo, neighbour_value);
            hi = std::max(hi, neighbour_value);
        }
        minimum[cell] = lo;
        maximum[cell] = hi;
    });

    auto& limiter=workspace.limiter;
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double alpha = 1.0;
        const double center_value = cells[cell];
        const Vec3 center = mesh.cells()[cell].center;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            const Vec3 delta_x = mesh.faces()[fi].center - center;
            const double delta = dot(gradient[cell], delta_x);
            if (delta > 0.0) {
                const double room = maximum[cell] - center_value;
                alpha = std::min(alpha, room <= 0.0 ? 0.0 : room / delta);
            } else if (delta < 0.0) {
                const double room = minimum[cell] - center_value;
                alpha = std::min(alpha, room >= 0.0 ? 0.0 : room / delta);
            }
        }
        limiter[cell] = std::clamp(alpha, 0.0, 1.0);
    });

    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        const Face& face = mesh.faces()[fi];
        if (face.boundary()) {
            if (flux[fi] < 0.0 && !boundary.empty()) {
                faces[fi] = boundary[fi];
                return;
            }
            const std::size_t upwind = face.owner;
            const double reconstructed = cells[upwind]
                + limiter[upwind] * dot(gradient[upwind], face.center - mesh.cells()[upwind].center);
            faces[fi] = std::clamp(reconstructed, minimum[upwind], maximum[upwind]);
            return;
        }
        const std::size_t upwind = flux[fi] >= 0.0 ? face.owner : face.neighbour;
        const double reconstructed = cells[upwind]
            + limiter[upwind] * dot(gradient[upwind], face.center - mesh.cells()[upwind].center);
        faces[fi] = std::clamp(reconstructed, minimum[upwind], maximum[upwind]);
    });
}

template<class T>
std::vector<T> interpolate_basic(const PolyMesh& mesh,
                                 std::span<const T> cells,
                                 FaceInterpolationScheme scheme,
                                 std::span<const double> flux,
                                 std::span<const T> boundary) {
    check_field(mesh, cells, boundary);
    if (scheme == FaceInterpolationScheme::upwind) check_flux(mesh, flux);

    std::vector<T> faces(mesh.face_count());
    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        const Face& face = mesh.faces()[fi];
        if (face.boundary()) {
            if (!boundary.empty() && (scheme == FaceInterpolationScheme::linear || flux[fi] < 0.0)) {
                faces[fi] = boundary[fi];
            } else {
                faces[fi] = cells[face.owner];
            }
            return;
        }
        if (scheme == FaceInterpolationScheme::upwind) {
            faces[fi] = flux[fi] >= 0.0 ? cells[face.owner] : cells[face.neighbour];
        } else {
            faces[fi] = linear_face_value(mesh, face, cells);
        }
    });
    return faces;
}

[[nodiscard]] double oriented_flux(const Face& face, std::size_t cell, double phi) noexcept {
    return face.owner == cell ? phi : -phi;
}

} // namespace

void ScalarSchemeWorkspace::resize(const PolyMesh& mesh) {
    gradient.resize(mesh.cell_count());minimum.resize(mesh.cell_count());
    maximum.resize(mesh.cell_count());limiter.resize(mesh.cell_count());face_values.resize(mesh.face_count());
}
std::size_t ScalarSchemeWorkspace::allocated_bytes() const noexcept {
    return gradient.capacity()*sizeof(Vec3)+(minimum.capacity()+maximum.capacity()
        +limiter.capacity()+face_values.capacity())*sizeof(double);
}

void interpolate_scalar_to_faces_into(const PolyMesh& mesh,std::span<const double> cells,
    std::span<double> faces,ScalarSchemeWorkspace& workspace,FaceInterpolationScheme scheme,
    std::span<const double> flux,std::span<const double> boundary) {
    check_field(mesh,cells,boundary);
    if(faces.size()!=mesh.face_count()) throw std::invalid_argument("face output size mismatch");
    if(scheme!=FaceInterpolationScheme::linear) check_flux(mesh,flux);
    if(scheme==FaceInterpolationScheme::bounded_linear){
        workspace.resize(mesh);
        bounded_linear_scalar_faces(mesh,cells,flux,boundary,faces,workspace);return;
    }
    const bool tvd=scheme==FaceInterpolationScheme::minmod||scheme==FaceInterpolationScheme::van_leer;
    if(tvd){
        workspace.resize(mesh);
        least_squares_gradient_scalar_into(mesh,cells,workspace.gradient,boundary);
    }
    cfd::core::parallel_for(mesh.face_count(),[&](std::size_t fi){
        const auto& face=mesh.faces()[fi];
        if(face.boundary()){
            faces[fi]=!boundary.empty()&&(scheme==FaceInterpolationScheme::linear||flux[fi]<0.0)
                ?boundary[fi]:cells[face.owner];return;
        }
        if(scheme==FaceInterpolationScheme::linear){faces[fi]=linear_face_value(mesh,face,cells);return;}
        const auto up=flux[fi]>=0.0?face.owner:face.neighbour;
        const auto down=flux[fi]>=0.0?face.neighbour:face.owner;
        faces[fi]=cells[up];
        const double delta=cells[down]-cells[up];
        if(!tvd||std::abs(delta)<1.0e-30)return;
        const Vec3 d=mesh.cells()[down].center-mesh.cells()[up].center;
        const double r=2.0*dot(workspace.gradient[up],d)/delta-1.0;
        const double psi=scheme==FaceInterpolationScheme::minmod?std::clamp(r,0.0,1.0)
            :(r>0.0?2.0/(1.0+1.0/r):0.0);
        const double fraction=std::clamp(dot(face.center-mesh.cells()[up].center,d)/dot(d,d),0.0,1.0);
        faces[fi]=std::clamp(cells[up]+fraction*psi*delta,std::min(cells[up],cells[down]),std::max(cells[up],cells[down]));
    });
}

std::vector<double> interpolate_scalar_to_faces(const PolyMesh& mesh,
                                                 std::span<const double> cells,
                                                 FaceInterpolationScheme scheme,
                                                 std::span<const double> flux,
                                                 std::span<const double> boundary) {
    ScalarSchemeWorkspace workspace;
    std::vector<double> result(mesh.face_count());
    interpolate_scalar_to_faces_into(mesh,cells,result,workspace,scheme,flux,boundary);
    return result;
}

std::vector<Vec3> interpolate_vector_to_faces(const PolyMesh& mesh,
                                               std::span<const Vec3> cells,
                                               FaceInterpolationScheme scheme,
                                               std::span<const double> flux,
                                               std::span<const Vec3> boundary) {
    if (scheme == FaceInterpolationScheme::linear || scheme == FaceInterpolationScheme::upwind) {
        return interpolate_basic(mesh, cells, scheme, flux, boundary);
    }
    check_field(mesh, cells, boundary);
    check_flux(mesh, flux);
    std::vector<double> x(cells.size()), y(cells.size()), z(cells.size());
    std::vector<double> bx, by, bz;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        x[i] = cells[i].x; y[i] = cells[i].y; z[i] = cells[i].z;
    }
    if (!boundary.empty()) {
        bx.resize(boundary.size()); by.resize(boundary.size()); bz.resize(boundary.size());
        for (std::size_t i = 0; i < boundary.size(); ++i) {
            bx[i] = boundary[i].x; by[i] = boundary[i].y; bz[i] = boundary[i].z;
        }
    }
    const auto fx = interpolate_scalar_to_faces(mesh,x,scheme,flux,bx);
    const auto fy = interpolate_scalar_to_faces(mesh,y,scheme,flux,by);
    const auto fz = interpolate_scalar_to_faces(mesh,z,scheme,flux,bz);
    std::vector<Vec3> result(mesh.face_count());
    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t i) {
        result[i] = {fx[i], fy[i], fz[i]};
    });
    return result;
}

std::vector<double> face_flux_from_velocity(const PolyMesh& mesh,
                                             std::span<const Vec3> velocity,
                                             std::span<const Vec3> boundary) {
    const auto face_velocity = interpolate_vector_to_faces(
        mesh, velocity, FaceInterpolationScheme::linear, {}, boundary);
    std::vector<double> flux(mesh.face_count(), 0.0);
    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        flux[fi] = dot(face_velocity[fi], mesh.faces()[fi].area);
    });
    return flux;
}

std::vector<double> convective_divergence_scalar(const PolyMesh& mesh,
                                                  std::span<const double> values,
                                                  std::span<const double> flux,
                                                  FaceInterpolationScheme scheme,
                                                  std::span<const double> boundary) {
    ScalarSchemeWorkspace workspace;
    std::vector<double> result(mesh.cell_count(), 0.0);
    convective_divergence_scalar_into(mesh,values,flux,result,workspace,scheme,boundary);
    return result;
}

void convective_divergence_scalar_into(const PolyMesh& mesh,std::span<const double> values,
    std::span<const double> flux,std::span<double> result,ScalarSchemeWorkspace& workspace,
    FaceInterpolationScheme scheme,std::span<const double> boundary) {
    check_flux(mesh,flux);
    if(result.size()!=mesh.cell_count()) throw std::invalid_argument("divergence output size mismatch");
    workspace.resize(mesh);
    interpolate_scalar_to_faces_into(mesh,values,workspace.face_values,workspace,scheme,flux,boundary);
    const auto& face_values=workspace.face_values;
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double sum = 0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            sum += oriented_flux(mesh.faces()[fi], cell, flux[fi]) * face_values[fi];
        }
        result[cell] = sum / mesh.cells()[cell].volume;
    });
}

std::vector<Vec3> convective_divergence_vector(const PolyMesh& mesh,
                                                std::span<const Vec3> values,
                                                std::span<const double> flux,
                                                FaceInterpolationScheme scheme,
                                                std::span<const Vec3> boundary) {
    check_flux(mesh, flux);
    const auto face_values = interpolate_vector_to_faces(mesh, values, scheme, flux, boundary);
    std::vector<Vec3> result(mesh.cell_count());
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        Vec3 sum{};
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            sum += face_values[fi] * oriented_flux(mesh.faces()[fi], cell, flux[fi]);
        }
        result[cell] = sum / mesh.cells()[cell].volume;
    });
    return result;
}

} // namespace cfd::fvm
