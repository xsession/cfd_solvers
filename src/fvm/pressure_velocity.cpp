#include "cfd/fvm/pressure_velocity.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/fvm/operators.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace {

void check_patch_sizes(const PolyMesh& mesh,
                       std::span<const VelocityBoundaryCondition> velocity_boundary,
                       std::span<const PressureBoundaryCondition> pressure_boundary) {
    if (velocity_boundary.size() != mesh.patches().size()) {
        throw std::invalid_argument("velocity boundary condition count must match mesh patches");
    }
    if (pressure_boundary.size() != mesh.patches().size()) {
        throw std::invalid_argument("pressure boundary condition count must match mesh patches");
    }
}

[[nodiscard]] double face_interpolate_scalar(const PolyMesh& mesh,
                                              const Face& face,
                                              std::span<const double> values) {
    if (face.boundary()) return values[face.owner];
    const double dof = magnitude(face.center - mesh.cells()[face.owner].center);
    const double dnf = magnitude(mesh.cells()[face.neighbour].center - face.center);
    const double total = dof + dnf;
    if (!(total > 0.0)) throw std::runtime_error("degenerate face interpolation distance");
    return (dnf * values[face.owner] + dof * values[face.neighbour]) / total;
}

[[nodiscard]] Vec3 face_interpolate_vector(const PolyMesh& mesh,
                                            const Face& face,
                                            std::span<const Vec3> values) {
    if (face.boundary()) return values[face.owner];
    const double dof = magnitude(face.center - mesh.cells()[face.owner].center);
    const double dnf = magnitude(mesh.cells()[face.neighbour].center - face.center);
    const double total = dof + dnf;
    if (!(total > 0.0)) throw std::runtime_error("degenerate face interpolation distance");
    return (values[face.owner] * dnf + values[face.neighbour] * dof) / total;
}

[[nodiscard]] Vec3 project_slip(Vec3 owner, Vec3 wall, Vec3 area) {
    const double area_mag = magnitude(area);
    if (!(area_mag > 0.0)) return owner;
    const Vec3 n = area / area_mag;
    const Vec3 relative = owner - wall;
    return wall + relative - n * dot(relative, n);
}

[[nodiscard]] double outward_flux(const Face& face, std::size_t cell, double oriented_flux) noexcept {
    return face.owner == cell ? oriented_flux : -oriented_flux;
}

} // namespace

FaceOrthogonalDecomposition decompose_face_area(const PolyMesh& mesh,
                                                 std::size_t face_index) {
    if (face_index >= mesh.face_count()) throw std::out_of_range("face index out of range");
    const auto& face = mesh.faces()[face_index];
    const Vec3 delta = face.boundary()
        ? face.center - mesh.cells()[face.owner].center
        : mesh.cells()[face.neighbour].center - mesh.cells()[face.owner].center;
    const double d2 = dot(delta, delta);
    if (!(d2 > 0.0)) throw std::runtime_error("degenerate face owner-neighbour displacement");
    const double metric = dot(face.area, delta) / d2;
    if (!(metric > 0.0)) {
        throw std::runtime_error("invalid mesh orientation: face area does not point owner to neighbour/boundary");
    }
    const Vec3 orthogonal = delta * metric;
    return {delta, orthogonal, face.area - orthogonal, metric};
}

std::vector<Vec3> boundary_velocity_values(
    const PolyMesh& mesh,
    std::span<const Vec3> cell_velocity,
    std::span<const VelocityBoundaryCondition> patch_conditions) {
    if (cell_velocity.size() != mesh.cell_count()) throw std::invalid_argument("velocity field size mismatch");
    if (patch_conditions.size() != mesh.patches().size()) throw std::invalid_argument("velocity patch count mismatch");
    std::vector<Vec3> values(mesh.face_count());
    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        const auto& face = mesh.faces()[fi];
        if (!face.boundary()) return;
        const auto& bc = patch_conditions[face.patch];
        switch (bc.type) {
        case VelocityBoundaryType::fixedValue:
            values[fi] = bc.value;
            break;
        case VelocityBoundaryType::zeroGradient:
            values[fi] = cell_velocity[face.owner];
            break;
        case VelocityBoundaryType::slip:
            values[fi] = project_slip(cell_velocity[face.owner], bc.value, face.area);
            break;
        }
    });
    return values;
}

std::vector<double> boundary_pressure_values(
    const PolyMesh& mesh,
    std::span<const double> cell_pressure,
    std::span<const PressureBoundaryCondition> patch_conditions) {
    if (cell_pressure.size() != mesh.cell_count()) throw std::invalid_argument("pressure field size mismatch");
    if (patch_conditions.size() != mesh.patches().size()) throw std::invalid_argument("pressure patch count mismatch");
    std::vector<double> values(mesh.face_count(), 0.0);
    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        const auto& face = mesh.faces()[fi];
        if (!face.boundary()) return;
        const auto& bc = patch_conditions[face.patch];
        values[fi] = bc.type == PressureBoundaryType::fixedValue ? bc.value : cell_pressure[face.owner];
    });
    return values;
}

std::vector<double> rhie_chow_face_flux(
    const PolyMesh& mesh,
    std::span<const Vec3> h_by_a,
    std::span<const double> pressure_mobility,
    std::span<const double> pressure,
    std::span<const VelocityBoundaryCondition> velocity_boundary,
    std::span<const PressureBoundaryCondition> pressure_boundary,
    bool include_nonorthogonal_correction) {
    if (h_by_a.size() != mesh.cell_count() || pressure_mobility.size() != mesh.cell_count() ||
        pressure.size() != mesh.cell_count()) {
        throw std::invalid_argument("Rhie-Chow cell field size mismatch");
    }
    check_patch_sizes(mesh, velocity_boundary, pressure_boundary);
    const auto p_boundary = boundary_pressure_values(mesh, pressure, pressure_boundary);
    const auto h_boundary = boundary_velocity_values(mesh, h_by_a, velocity_boundary);
    const auto grad_p = gauss_gradient_scalar(mesh, pressure, p_boundary);
    std::vector<double> flux(mesh.face_count(), 0.0);

    cfd::core::parallel_for(mesh.face_count(), [&](std::size_t fi) {
        const auto& face = mesh.faces()[fi];
        const auto geometry = decompose_face_area(mesh, fi);
        if (face.boundary() && velocity_boundary[face.patch].type == VelocityBoundaryType::fixedValue) {
            flux[fi] = dot(h_boundary[fi], face.area);
            return;
        }

        const Vec3 hf = face.boundary() ? h_boundary[fi] : face_interpolate_vector(mesh, face, h_by_a);
        const double mobility = face_interpolate_scalar(mesh, face, pressure_mobility);
        double pressure_flux = 0.0;
        if (!face.boundary()) {
            pressure_flux = mobility * geometry.orthogonal_metric *
                (pressure[face.neighbour] - pressure[face.owner]);
            if (include_nonorthogonal_correction) {
                const Vec3 gf = face_interpolate_vector(mesh, face, grad_p);
                pressure_flux += mobility * dot(gf, geometry.nonorthogonal_area);
            }
        } else if (pressure_boundary[face.patch].type == PressureBoundaryType::fixedValue) {
            pressure_flux = mobility * geometry.orthogonal_metric *
                (p_boundary[fi] - pressure[face.owner]);
            if (include_nonorthogonal_correction) {
                pressure_flux += mobility * dot(grad_p[face.owner], geometry.nonorthogonal_area);
            }
        }
        flux[fi] = dot(hf, face.area) - pressure_flux;
    });
    return flux;
}

double flux_divergence_l2(const PolyMesh& mesh, std::span<const double> face_flux) {
    if (face_flux.size() != mesh.face_count()) throw std::invalid_argument("face flux size mismatch");
    const double sum = cfd::core::parallel_sum(mesh.cell_count(), [&](std::size_t cell) {
        double net = 0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            net += outward_flux(mesh.faces()[fi], cell, face_flux[fi]);
        }
        const double div = net / mesh.cells()[cell].volume;
        return div * div;
    });
    return std::sqrt(sum / static_cast<double>(mesh.cell_count()));
}

double flux_divergence_max(const PolyMesh& mesh, std::span<const double> face_flux) {
    if (face_flux.size() != mesh.face_count()) throw std::invalid_argument("face flux size mismatch");
    double result = 0.0;
    for (std::size_t cell = 0; cell < mesh.cell_count(); ++cell) {
        double net = 0.0;
        for (const std::size_t fi : mesh.cell_faces()[cell]) {
            net += outward_flux(mesh.faces()[fi], cell, face_flux[fi]);
        }
        result = std::max(result, std::abs(net / mesh.cells()[cell].volume));
    }
    return result;
}

} // namespace cfd::fvm
