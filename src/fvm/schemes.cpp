#include "cfd/fvm/schemes.hpp"
#include "cfd/core/parallel.hpp"

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

template<class T>
std::vector<T> interpolate_impl(const PolyMesh& mesh,
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

std::vector<double> interpolate_scalar_to_faces(const PolyMesh& mesh,
                                                 std::span<const double> cells,
                                                 FaceInterpolationScheme scheme,
                                                 std::span<const double> flux,
                                                 std::span<const double> boundary) {
    return interpolate_impl(mesh, cells, scheme, flux, boundary);
}

std::vector<Vec3> interpolate_vector_to_faces(const PolyMesh& mesh,
                                               std::span<const Vec3> cells,
                                               FaceInterpolationScheme scheme,
                                               std::span<const double> flux,
                                               std::span<const Vec3> boundary) {
    return interpolate_impl(mesh, cells, scheme, flux, boundary);
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
    check_flux(mesh, flux);
    const auto face_values = interpolate_scalar_to_faces(mesh, values, scheme, flux, boundary);
    std::vector<double> result(mesh.cell_count(), 0.0);
    cfd::core::parallel_for(mesh.cell_count(), [&](std::size_t cell) {
        double sum = 0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            sum += oriented_flux(mesh.faces()[fi], cell, flux[fi]) * face_values[fi];
        }
        result[cell] = sum / mesh.cells()[cell].volume;
    });
    return result;
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
