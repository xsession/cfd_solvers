#include "cfd/fvm/resident_sycl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

namespace cfd::fvm {
namespace {

[[nodiscard]] double face_orthogonal_metric(const PolyMesh& mesh, const Face& face) {
    const Vec3 delta = face.boundary()
        ? face.center - mesh.cells()[face.owner].center
        : mesh.cells()[face.neighbour].center - mesh.cells()[face.owner].center;
    const double d2 = dot(delta, delta);
    if (!(d2 > 0.0)) throw std::runtime_error("degenerate FVM face displacement");
    const double metric = dot(face.area, delta) / d2;
    if (!(metric > 0.0)) throw std::runtime_error("invalid FVM face orientation");
    return metric;
}

#if defined(CFD_HAS_SYCL)
[[nodiscard]] double face_normal_metric(const PolyMesh& mesh, const Face& face) {
    const double area = magnitude(face.area);
    if (!(area > 0.0)) throw std::runtime_error("degenerate FVM face area");
    const Vec3 normal = face.area / area;
    const Vec3 delta = face.boundary()
        ? face.center - mesh.cells()[face.owner].center
        : mesh.cells()[face.neighbour].center - mesh.cells()[face.owner].center;
    const double distance = std::abs(dot(delta, normal));
    if (!(distance > 0.0)) throw std::runtime_error("degenerate FVM face-normal distance");
    return area / distance;
}
#endif

[[nodiscard]] double face_mobility(const PolyMesh& mesh,
                                   const Face& face,
                                   std::span<const double> mobility) {
    if (mobility.empty()) return 1.0;
    if (face.boundary()) return mobility[face.owner];
    const double dof = magnitude(face.center - mesh.cells()[face.owner].center);
    const double dnf = magnitude(mesh.cells()[face.neighbour].center - face.center);
    const double d = dof + dnf;
    if (!(d > 0.0)) throw std::runtime_error("degenerate FVM interpolation distance");
    return (dnf * mobility[face.owner] + dof * mobility[face.neighbour]) / d;
}

#if defined(CFD_HAS_SYCL)
[[nodiscard]] cfd::core::CsrMatrix build_cell_neighbour_pattern(const PolyMesh& mesh) {
    cfd::core::CsrBuilder builder(mesh.cell_count(), mesh.cell_count());
    for (std::size_t c = 0U; c < mesh.cell_count(); ++c) builder.add(c, c, 1.0);
    for (const auto& face : mesh.faces()) {
        if (face.boundary()) continue;
        builder.add(face.owner, face.neighbour, -1.0);
        builder.add(face.neighbour, face.owner, -1.0);
    }
    return builder.build();
}

#endif

[[nodiscard]] cfd::core::CsrMatrix build_pressure_matrix_impl(
    const PolyMesh& mesh,
    std::span<const double> mobility,
    bool fixed_value_boundaries,
    std::size_t reference_cell) {
    if (!mobility.empty() && mobility.size() != mesh.cell_count()) {
        throw std::invalid_argument("pressure mobility size mismatch");
    }
    for (double value : mobility) {
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::invalid_argument("pressure mobility must be finite and positive");
        }
    }

    if (reference_cell != invalid_cell && reference_cell >= mesh.cell_count()) {
        throw std::out_of_range("pressure reference cell out of range");
    }
    cfd::core::CsrBuilder builder(mesh.cell_count(), mesh.cell_count());
    for (const auto& face : mesh.faces()) {
        const double coeff = face_mobility(mesh, face, mobility) * face_orthogonal_metric(mesh, face);
        if (!face.boundary()) {
            const bool owner_ref = face.owner == reference_cell;
            const bool neighbour_ref = face.neighbour == reference_cell;
            if (!owner_ref && !neighbour_ref) {
                builder.add(face.owner, face.owner, coeff);
                builder.add(face.owner, face.neighbour, -coeff);
                builder.add(face.neighbour, face.neighbour, coeff);
                builder.add(face.neighbour, face.owner, -coeff);
            } else if (owner_ref && !neighbour_ref) {
                builder.add(face.neighbour, face.neighbour, coeff);
            } else if (!owner_ref && neighbour_ref) {
                builder.add(face.owner, face.owner, coeff);
            }
        } else if (fixed_value_boundaries && face.owner != reference_cell) {
            builder.add(face.owner, face.owner, coeff);
        }
    }
    if (reference_cell != invalid_cell) builder.add(reference_cell, reference_cell, 1.0);
    return builder.build();
}

} // namespace

cfd::core::CsrMatrix build_orthogonal_pressure_matrix(
    const PolyMesh& mesh,
    double mobility,
    bool fixed_value_boundaries,
    std::size_t reference_cell) {
    if (!(mobility > 0.0) || !std::isfinite(mobility)) {
        throw std::invalid_argument("pressure mobility must be finite and positive");
    }
    std::vector<double> values(mesh.cell_count(), mobility);
    return build_pressure_matrix_impl(mesh, values, fixed_value_boundaries, reference_cell);
}

cfd::core::CsrMatrix build_orthogonal_pressure_matrix(
    const PolyMesh& mesh,
    std::span<const double> cell_mobility,
    bool fixed_value_boundaries,
    std::size_t reference_cell) {
    return build_pressure_matrix_impl(mesh, cell_mobility, fixed_value_boundaries, reference_cell);
}

#if defined(CFD_HAS_SYCL)
namespace {

template<class T>
void free_if(T*& ptr, sycl::queue& queue) noexcept {
    if (!ptr) return;
    try { sycl::free(ptr, queue); } catch (...) {}
    ptr = nullptr;
}

void require_pointer(const void* ptr, const char* what) {
    if (!ptr) throw std::invalid_argument(what);
}

} // namespace

ResidentPolyMeshSycl::ResidentPolyMeshSycl(const PolyMesh& mesh, sycl::device device)
    : queue_(device, sycl::property::queue::in_order{}) {
    initialize(mesh);
}

ResidentPolyMeshSycl::ResidentPolyMeshSycl(const PolyMesh& mesh, sycl::queue queue)
    : queue_(std::move(queue)) {
    initialize(mesh);
}

ResidentPolyMeshSycl::~ResidentPolyMeshSycl() noexcept { release(); }

void ResidentPolyMeshSycl::initialize(const PolyMesh& mesh) {
    if (mesh.cell_count() == 0U || mesh.face_count() == 0U) {
        throw std::invalid_argument("resident FVM mesh must be non-empty");
    }
    if (!queue_.get_device().has(sycl::aspect::usm_device_allocations)) {
        throw std::runtime_error("selected SYCL device lacks USM device allocations");
    }

    cell_count_ = mesh.cell_count();
    face_count_ = mesh.face_count();
    std::vector<std::size_t> offsets(cell_count_ + 1U, 0U);
    std::vector<std::size_t> adjacency;
    adjacency.reserve(face_count_ * 2U);
    for (std::size_t cell = 0U; cell < cell_count_; ++cell) {
        adjacency.insert(adjacency.end(), mesh.cell_faces()[cell].begin(), mesh.cell_faces()[cell].end());
        offsets[cell + 1U] = adjacency.size();
    }
    adjacency_count_ = adjacency.size();

    std::vector<std::size_t> owner(face_count_), neighbour(face_count_);
    std::vector<double> cell_geometry(4U * cell_count_, 0.0);
    std::vector<double> face_geometry(9U * face_count_, 0.0);
    for (std::size_t cell = 0U; cell < cell_count_; ++cell) {
        const auto& c = mesh.cells()[cell];
        cell_geometry[cell] = c.center.x;
        cell_geometry[cell_count_ + cell] = c.center.y;
        cell_geometry[2U * cell_count_ + cell] = c.center.z;
        cell_geometry[3U * cell_count_ + cell] = c.volume;
    }
    for (std::size_t fi = 0U; fi < face_count_; ++fi) {
        const auto& f = mesh.faces()[fi];
        owner[fi] = f.owner;
        neighbour[fi] = f.neighbour;
        face_geometry[fi] = f.center.x;
        face_geometry[face_count_ + fi] = f.center.y;
        face_geometry[2U * face_count_ + fi] = f.center.z;
        face_geometry[3U * face_count_ + fi] = f.area.x;
        face_geometry[4U * face_count_ + fi] = f.area.y;
        face_geometry[5U * face_count_ + fi] = f.area.z;
        double owner_weight = 1.0;
        if (!f.boundary()) {
            const double dof = magnitude(f.center - mesh.cells()[f.owner].center);
            const double dnf = magnitude(mesh.cells()[f.neighbour].center - f.center);
            const double d = dof + dnf;
            if (!(d > 0.0)) throw std::runtime_error("degenerate resident FVM face interpolation");
            owner_weight = dnf / d;
        }
        face_geometry[6U * face_count_ + fi] = owner_weight;
        face_geometry[7U * face_count_ + fi] = face_normal_metric(mesh, f);
        face_geometry[8U * face_count_ + fi] = face_orthogonal_metric(mesh, f);
    }

    face_owner_ = sycl::malloc_device<std::size_t>(face_count_, queue_);
    face_neighbour_ = sycl::malloc_device<std::size_t>(face_count_, queue_);
    cell_face_offsets_ = sycl::malloc_device<std::size_t>(cell_count_ + 1U, queue_);
    cell_face_indices_ = sycl::malloc_device<std::size_t>(adjacency_count_, queue_);
    cell_geometry_ = sycl::malloc_device<double>(cell_geometry.size(), queue_);
    face_geometry_ = sycl::malloc_device<double>(face_geometry.size(), queue_);
    if (!face_owner_ || !face_neighbour_ || !cell_face_offsets_ || !cell_face_indices_ ||
        !cell_geometry_ || !face_geometry_) {
        release();
        throw std::bad_alloc{};
    }

    const std::size_t owner_bytes = face_count_ * sizeof(std::size_t);
    const std::size_t offset_bytes = (cell_count_ + 1U) * sizeof(std::size_t);
    const std::size_t adjacency_bytes = adjacency_count_ * sizeof(std::size_t);
    const std::size_t cell_geometry_bytes = cell_geometry.size() * sizeof(double);
    const std::size_t face_geometry_bytes = face_geometry.size() * sizeof(double);
    queue_.memcpy(face_owner_, owner.data(), owner_bytes);
    queue_.memcpy(face_neighbour_, neighbour.data(), owner_bytes);
    queue_.memcpy(cell_face_offsets_, offsets.data(), offset_bytes);
    queue_.memcpy(cell_face_indices_, adjacency.data(), adjacency_bytes);
    queue_.memcpy(cell_geometry_, cell_geometry.data(), cell_geometry_bytes);
    queue_.memcpy(face_geometry_, face_geometry.data(), face_geometry_bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(2U * owner_bytes + offset_bytes + adjacency_bytes +
                                          cell_geometry_bytes + face_geometry_bytes);
    transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::release() noexcept {
    try { queue_.wait_and_throw(); } catch (...) {}
    free_if(face_owner_, queue_);
    free_if(face_neighbour_, queue_);
    free_if(cell_face_offsets_, queue_);
    free_if(cell_face_indices_, queue_);
    free_if(cell_geometry_, queue_);
    free_if(face_geometry_, queue_);
}

std::size_t ResidentPolyMeshSycl::resident_bytes() const noexcept {
    return 2U * face_count_ * sizeof(std::size_t) +
           (cell_count_ + 1U + adjacency_count_) * sizeof(std::size_t) +
           (4U * cell_count_ + 9U * face_count_) * sizeof(double);
}

void ResidentPolyMeshSycl::upload_cell_scalar(std::span<const double> host, double* device) const {
    if (host.size() != cell_count_) throw std::invalid_argument("resident FVM cell scalar size mismatch");
    require_pointer(device, "null resident FVM cell scalar destination");
    const std::size_t bytes = cell_count_ * sizeof(double);
    queue_.memcpy(device, host.data(), bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(bytes); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::download_cell_scalar(const double* device, std::span<double> host) const {
    if (host.size() != cell_count_) throw std::invalid_argument("resident FVM cell scalar size mismatch");
    require_pointer(device, "null resident FVM cell scalar source");
    const std::size_t bytes = cell_count_ * sizeof(double);
    queue_.memcpy(host.data(), device, bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(bytes); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::upload_face_scalar(std::span<const double> host, double* device) const {
    if (host.size() != face_count_) throw std::invalid_argument("resident FVM face scalar size mismatch");
    require_pointer(device, "null resident FVM face scalar destination");
    const std::size_t bytes = face_count_ * sizeof(double);
    queue_.memcpy(device, host.data(), bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(bytes); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::download_face_scalar(const double* device, std::span<double> host) const {
    if (host.size() != face_count_) throw std::invalid_argument("resident FVM face scalar size mismatch");
    require_pointer(device, "null resident FVM face scalar source");
    const std::size_t bytes = face_count_ * sizeof(double);
    queue_.memcpy(host.data(), device, bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(bytes); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::upload_cell_vector(std::span<const Vec3> host, double* device) const {
    if (host.size() != cell_count_) throw std::invalid_argument("resident FVM cell vector size mismatch");
    require_pointer(device, "null resident FVM cell vector destination");
    std::vector<double> soa(3U * cell_count_);
    for (std::size_t i = 0U; i < cell_count_; ++i) {
        soa[i] = host[i].x; soa[cell_count_ + i] = host[i].y; soa[2U * cell_count_ + i] = host[i].z;
    }
    const std::size_t bytes = soa.size() * sizeof(double);
    queue_.memcpy(device, soa.data(), bytes).wait_and_throw();
    transfer_stats_.record_host_to_device(bytes); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::download_cell_vector(const double* device, std::span<Vec3> host) const {
    if (host.size() != cell_count_) throw std::invalid_argument("resident FVM cell vector size mismatch");
    require_pointer(device, "null resident FVM cell vector source");
    std::vector<double> soa(3U * cell_count_);
    const std::size_t bytes = soa.size() * sizeof(double);
    queue_.memcpy(soa.data(), device, bytes).wait_and_throw();
    transfer_stats_.record_device_to_host(bytes); transfer_stats_.record_synchronization();
    for (std::size_t i = 0U; i < cell_count_; ++i) {
        host[i] = {soa[i], soa[cell_count_ + i], soa[2U * cell_count_ + i]};
    }
}

void ResidentPolyMeshSycl::fill(double* device, std::size_t count, double value) const {
    require_pointer(device, "null resident FVM fill destination");
    queue_.parallel_for(sycl::range<1>(count), [=](sycl::id<1> id) { device[id[0]] = value; });
}

void ResidentPolyMeshSycl::set_cell_value(double* device, std::size_t cell, double value) const {
    require_pointer(device, "null resident FVM cell destination");
    if (cell >= cell_count_) throw std::out_of_range("resident FVM cell index out of range");
    queue_.parallel_for(sycl::range<1>(1U), [=](sycl::id<1>) { device[cell] = value; });
}

void ResidentPolyMeshSycl::wait() const {
    queue_.wait_and_throw(); transfer_stats_.record_synchronization();
}

void ResidentPolyMeshSycl::interpolate_scalar_to_faces(const double* cells, double* faces,
                                                        const double* boundary) const {
    require_pointer(cells, "null resident FVM cell scalar"); require_pointer(faces, "null resident FVM face scalar");
    const auto* owner = face_owner_; const auto* neighbour = face_neighbour_; const auto* geom = face_geometry_;
    const std::size_t nf = face_count_; const std::size_t invalid = invalid_cell;
    queue_.parallel_for(sycl::range<1>(nf), [=](sycl::id<1> id) {
        const std::size_t f = id[0], o = owner[f], n = neighbour[f];
        if (n == invalid) { faces[f] = boundary ? boundary[f] : cells[o]; return; }
        const double wo = geom[6U * nf + f]; faces[f] = wo * cells[o] + (1.0 - wo) * cells[n];
    });
}

void ResidentPolyMeshSycl::gauss_gradient_scalar(const double* values, double* gradient,
                                                  const double* boundary) const {
    require_pointer(values, "null resident FVM scalar for gradient"); require_pointer(gradient, "null resident FVM gradient");
    const auto* owner=face_owner_; const auto* neighbour=face_neighbour_; const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_; const auto* cg=cell_geometry_; const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_, nf=face_count_, invalid=invalid_cell;
    queue_.parallel_for(sycl::range<1>(nc), [=](sycl::id<1> id) {
        const std::size_t c=id[0]; double gx=0.0,gy=0.0,gz=0.0;
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f];
            double phi=0.0;
            if(n==invalid) phi=boundary?boundary[f]:values[o];
            else { const double wo=fg[6U*nf+f]; phi=wo*values[o]+(1.0-wo)*values[n]; }
            const double sign=o==c?1.0:-1.0;
            gx+=sign*fg[3U*nf+f]*phi; gy+=sign*fg[4U*nf+f]*phi; gz+=sign*fg[5U*nf+f]*phi;
        }
        const double inv=1.0/cg[3U*nc+c]; gradient[c]=gx*inv; gradient[nc+c]=gy*inv; gradient[2U*nc+c]=gz*inv;
    });
}

void ResidentPolyMeshSycl::gauss_divergence_vector(const double* values, double* divergence,
                                                    const double* boundary) const {
    require_pointer(values, "null resident FVM vector for divergence"); require_pointer(divergence, "null resident FVM divergence");
    const auto* owner=face_owner_; const auto* neighbour=face_neighbour_; const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_; const auto* cg=cell_geometry_; const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    queue_.parallel_for(sycl::range<1>(nc), [=](sycl::id<1> id){
        const std::size_t c=id[0]; double sum=0.0;
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f]; double vx=0.0,vy=0.0,vz=0.0;
            if(n==invalid){
                if(boundary){vx=boundary[f];vy=boundary[nf+f];vz=boundary[2U*nf+f];}
                else {vx=values[o];vy=values[nc+o];vz=values[2U*nc+o];}
            } else {
                const double wo=fg[6U*nf+f]; vx=wo*values[o]+(1.0-wo)*values[n];
                vy=wo*values[nc+o]+(1.0-wo)*values[nc+n]; vz=wo*values[2U*nc+o]+(1.0-wo)*values[2U*nc+n];
            }
            const double sign=o==c?1.0:-1.0;
            sum+=sign*(vx*fg[3U*nf+f]+vy*fg[4U*nf+f]+vz*fg[5U*nf+f]);
        }
        divergence[c]=sum/cg[3U*nc+c];
    });
}

void ResidentPolyMeshSycl::orthogonal_laplacian_scalar(const double* values,double diffusivity,double* laplacian,
                                                        const double* boundary) const {
    require_pointer(values,"null resident FVM scalar for laplacian"); require_pointer(laplacian,"null resident FVM laplacian");
    if(!(diffusivity>0.0)||!std::isfinite(diffusivity))throw std::invalid_argument("resident FVM diffusivity must be positive");
    const auto* owner=face_owner_; const auto* neighbour=face_neighbour_; const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_; const auto* cg=cell_geometry_; const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    queue_.parallel_for(sycl::range<1>(nc), [=](sycl::id<1> id){
        const std::size_t c=id[0]; double sum=0.0;
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f];
            if(n!=invalid){const std::size_t other=o==c?n:o;sum+=diffusivity*fg[7U*nf+f]*(values[other]-values[c]);}
            else if(boundary){sum+=diffusivity*fg[7U*nf+f]*(boundary[f]-values[c]);}
        }
        laplacian[c]=sum/cg[3U*nc+c];
    });
}

void ResidentPolyMeshSycl::predictor_face_flux(const double* velocity,double* flux,const double* boundary) const {
    require_pointer(velocity,"null resident FVM predictor velocity"); require_pointer(flux,"null resident FVM predictor flux");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    queue_.parallel_for(sycl::range<1>(nf),[=](sycl::id<1> id){
        const std::size_t f=id[0],o=owner[f],n=neighbour[f];double vx=0.0,vy=0.0,vz=0.0;
        if(n==invalid){if(boundary){vx=boundary[f];vy=boundary[nf+f];vz=boundary[2U*nf+f];}
                       else {vx=velocity[o];vy=velocity[nc+o];vz=velocity[2U*nc+o];}}
        else {const double wo=fg[6U*nf+f];vx=wo*velocity[o]+(1.0-wo)*velocity[n];
              vy=wo*velocity[nc+o]+(1.0-wo)*velocity[nc+n];vz=wo*velocity[2U*nc+o]+(1.0-wo)*velocity[2U*nc+n];}
        flux[f]=vx*fg[3U*nf+f]+vy*fg[4U*nf+f]+vz*fg[5U*nf+f];
    });
}

void ResidentPolyMeshSycl::pressure_rhs_from_flux(const double* flux,const double* mobility,double* rhs,
                                                   const double* boundary_pressure,
                                                   const std::uint8_t* boundary_pressure_kind,
                                                   const double* nonorthogonal_flux) const {
    require_pointer(flux,"null resident FVM pressure predictor flux");require_pointer(mobility,"null resident FVM pressure mobility");
    require_pointer(rhs,"null resident FVM pressure RHS");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* fg=face_geometry_;const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t c=id[0];double value=0.0;
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f];const double sign=o==c?1.0:-1.0;
            value-=sign*flux[f];
            if(nonorthogonal_flux) value+=sign*nonorthogonal_flux[f];
            if(boundary_pressure&&n==invalid&&(!boundary_pressure_kind||boundary_pressure_kind[f]==fixed)){
                value+=mobility[o]*fg[8U*nf+f]*boundary_pressure[f];
            }
        }
        rhs[c]=value;
    });
}

void ResidentPolyMeshSycl::pressure_correct_face_flux(const double* pressure,const double* mobility,double* flux,
                                                       const double* boundary_pressure,
                                                       const std::uint8_t* boundary_pressure_kind,
                                                       const double* nonorthogonal_flux) const {
    require_pointer(pressure,"null resident FVM pressure");require_pointer(mobility,"null resident FVM pressure mobility");
    require_pointer(flux,"null resident FVM pressure-corrected flux");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* fg=face_geometry_;
    const std::size_t nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nf),[=](sycl::id<1> id){
        const std::size_t f=id[0],o=owner[f],n=neighbour[f];
        if(n!=invalid){const double wo=fg[6U*nf+f];const double mf=wo*mobility[o]+(1.0-wo)*mobility[n];
            flux[f]+=mf*fg[8U*nf+f]*(pressure[o]-pressure[n]);}
        else if(boundary_pressure&&(!boundary_pressure_kind||boundary_pressure_kind[f]==fixed)){
            flux[f]+=mobility[o]*fg[8U*nf+f]*(pressure[o]-boundary_pressure[f]);
        }
        if(nonorthogonal_flux) flux[f]-=nonorthogonal_flux[f];
    });
}

void ResidentPolyMeshSycl::correct_velocity_from_pressure_gradient(const double* h,const double* mobility,
                                                                    const double* gradient,double* velocity) const {
    require_pointer(h,"null resident FVM H/A velocity");require_pointer(mobility,"null resident FVM pressure mobility");
    require_pointer(gradient,"null resident FVM pressure gradient");require_pointer(velocity,"null resident FVM corrected velocity");
    const std::size_t nc=cell_count_;
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const std::size_t c=id[0];const double m=mobility[c];
        velocity[c]=h[c]-m*gradient[c];velocity[nc+c]=h[nc+c]-m*gradient[nc+c];velocity[2U*nc+c]=h[2U*nc+c]-m*gradient[2U*nc+c];});
}

void ResidentPolyMeshSycl::divergence_face_flux(const double* flux,double* divergence) const {
    require_pointer(flux,"null resident FVM face flux");require_pointer(divergence,"null resident FVM flux divergence");
    const auto* owner=face_owner_;const auto* offsets=cell_face_offsets_;const auto* face_indices=cell_face_indices_;const auto* cg=cell_geometry_;
    const std::size_t nc=cell_count_;
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const std::size_t c=id[0];double sum=0.0;
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){const std::size_t f=face_indices[k];sum+=(owner[f]==c?1.0:-1.0)*flux[f];}
        divergence[c]=sum/cg[3U*nc+c];});
}

void ResidentPolyMeshSycl::boundary_velocity_values(const double* cell_velocity,
                                                        const std::uint8_t* boundary_kind,
                                                        const double* boundary_fixed,
                                                        double* boundary_velocity) const {
    require_pointer(cell_velocity,"null resident FVM boundary cell velocity");
    require_pointer(boundary_kind,"null resident FVM boundary kind");
    require_pointer(boundary_fixed,"null resident FVM fixed boundary velocity");
    require_pointer(boundary_velocity,"null resident FVM boundary velocity output");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(VelocityBoundaryType::fixedValue);
    const std::uint8_t slip=static_cast<std::uint8_t>(VelocityBoundaryType::slip);
    queue_.parallel_for(sycl::range<1>(nf),[=](sycl::id<1> id){
        const std::size_t f=id[0];
        if(neighbour[f]!=invalid){boundary_velocity[f]=0.0;boundary_velocity[nf+f]=0.0;boundary_velocity[2U*nf+f]=0.0;return;}
        const std::size_t o=owner[f];double vx=cell_velocity[o],vy=cell_velocity[nc+o],vz=cell_velocity[2U*nc+o];
        if(boundary_kind[f]==fixed){vx=boundary_fixed[f];vy=boundary_fixed[nf+f];vz=boundary_fixed[2U*nf+f];}
        else if(boundary_kind[f]==slip){
            const double wx=boundary_fixed[f],wy=boundary_fixed[nf+f],wz=boundary_fixed[2U*nf+f];
            const double ax=fg[3U*nf+f],ay=fg[4U*nf+f],az=fg[5U*nf+f];
            const double amag=sycl::sqrt(ax*ax+ay*ay+az*az);
            if(amag>0.0){const double nx=ax/amag,ny=ay/amag,nz=az/amag;
                const double rx=vx-wx,ry=vy-wy,rz=vz-wz;const double normal=rx*nx+ry*ny+rz*nz;
                vx=wx+rx-normal*nx;vy=wy+ry-normal*ny;vz=wz+rz-normal*nz;}
        }
        boundary_velocity[f]=vx;boundary_velocity[nf+f]=vy;boundary_velocity[2U*nf+f]=vz;
    });
}

void ResidentPolyMeshSycl::boundary_pressure_values(const double* cell_pressure,
                                                        const std::uint8_t* boundary_kind,
                                                        const double* boundary_fixed,
                                                        double* boundary_pressure) const {
    require_pointer(cell_pressure,"null resident FVM boundary cell pressure");
    require_pointer(boundary_kind,"null resident FVM pressure boundary kind");
    require_pointer(boundary_fixed,"null resident FVM fixed boundary pressure");
    require_pointer(boundary_pressure,"null resident FVM boundary pressure output");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;
    const std::size_t nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nf),[=](sycl::id<1> id){
        const std::size_t f=id[0];
        if(neighbour[f]!=invalid){boundary_pressure[f]=0.0;return;}
        boundary_pressure[f]=boundary_kind[f]==fixed?boundary_fixed[f]:cell_pressure[owner[f]];
    });
}

void ResidentPolyMeshSycl::assemble_momentum_system(const double* time_source,
                                                      const double* pressure_gradient,
                                                      const double* flux,
                                                      const std::uint8_t* boundary_kind,
                                                      const double* boundary_velocity,
                                                      double nu,double dt,bool include_convection,
                                                      double* diagonal,double* mobility,double* rhs) const {
    require_pointer(time_source,"null resident FVM momentum time source");require_pointer(pressure_gradient,"null resident FVM momentum pressure gradient");
    require_pointer(flux,"null resident FVM momentum flux");require_pointer(boundary_kind,"null resident FVM momentum boundary kind");
    require_pointer(boundary_velocity,"null resident FVM momentum boundary velocity");require_pointer(diagonal,"null resident FVM momentum diagonal");
    require_pointer(mobility,"null resident FVM momentum mobility");require_pointer(rhs,"null resident FVM momentum RHS");
    if(!(nu>=0.0)||!std::isfinite(nu)||!(dt>0.0)||!std::isfinite(dt)) throw std::invalid_argument("invalid resident momentum coefficients");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* cg=cell_geometry_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(VelocityBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t c=id[0];double ap=cg[3U*nc+c]/dt;double bx=time_source[c]*ap,by=time_source[nc+c]*ap,bz=time_source[2U*nc+c]*ap;
        bx-=pressure_gradient[c]*cg[3U*nc+c];by-=pressure_gradient[nc+c]*cg[3U*nc+c];bz-=pressure_gradient[2U*nc+c]*cg[3U*nc+c];
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f];const double sign=o==c?1.0:-1.0;
            const double phi=sign*flux[f];const double diff=nu*fg[8U*nf+f];
            if(n!=invalid){const double leaving=include_convection?(phi>0.0?phi:0.0):0.0;ap+=diff+leaving;}
            else if(boundary_kind[f]==fixed){
                const double entering=include_convection?(phi<0.0?-phi:0.0):0.0;const double leaving=include_convection?(phi>0.0?phi:0.0):0.0;
                ap+=diff+leaving;const double coeff=diff+entering;bx+=coeff*boundary_velocity[f];by+=coeff*boundary_velocity[nf+f];bz+=coeff*boundary_velocity[2U*nf+f];
            } else if(include_convection) ap+=phi;
        }
        if(!(ap>0.0)) ap=cg[3U*nc+c]/dt;
        diagonal[c]=ap;mobility[c]=cg[3U*nc+c]/ap;rhs[c]=bx;rhs[nc+c]=by;rhs[2U*nc+c]=bz;
    });
}

void ResidentPolyMeshSycl::momentum_jacobi_sweep(const double* current_velocity,const double* flux,
                                                   const std::uint8_t* boundary_kind,double nu,bool include_convection,
                                                   const double* diagonal,const double* rhs,double* next_velocity) const {
    require_pointer(current_velocity,"null resident FVM momentum iterate");require_pointer(flux,"null resident FVM momentum flux");
    require_pointer(boundary_kind,"null resident FVM momentum boundary kind");require_pointer(diagonal,"null resident FVM momentum diagonal");
    require_pointer(rhs,"null resident FVM momentum RHS");require_pointer(next_velocity,"null resident FVM momentum next iterate");
    (void)boundary_kind;
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* fg=face_geometry_;const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t c=id[0];double sx=rhs[c],sy=rhs[nc+c],sz=rhs[2U*nc+c];
        for(std::size_t k=offsets[c];k<offsets[c+1U];++k){
            const std::size_t f=face_indices[k],o=owner[f],n=neighbour[f];if(n==invalid) continue;
            const std::size_t other=o==c?n:o;const double sign=o==c?1.0:-1.0;const double phi=sign*flux[f];
            const double entering=include_convection?(phi<0.0?-phi:0.0):0.0;const double coupling=nu*fg[8U*nf+f]+entering;
            sx+=coupling*current_velocity[other];sy+=coupling*current_velocity[nc+other];sz+=coupling*current_velocity[2U*nc+other];
        }
        const double inv=1.0/diagonal[c];next_velocity[c]=sx*inv;next_velocity[nc+c]=sy*inv;next_velocity[2U*nc+c]=sz*inv;
    });
}

void ResidentPolyMeshSycl::assemble_momentum_matrix_values(
    const std::size_t* row_offsets,const std::size_t* columns,std::size_t nonzeros,
    const double* flux,const std::uint8_t* boundary_kind,double nu,double dt,
    bool include_convection,double* values) const {
    require_pointer(row_offsets,"null resident momentum row offsets");
    require_pointer(columns,"null resident momentum columns");
    require_pointer(flux,"null resident momentum face flux for matrix assembly");
    require_pointer(boundary_kind,"null resident momentum boundary kind for matrix assembly");
    require_pointer(values,"null resident momentum matrix values");
    if(nonzeros==0U||!(nu>=0.0)||!std::isfinite(nu)||!(dt>0.0)||!std::isfinite(dt))
        throw std::invalid_argument("invalid resident momentum matrix controls");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* cg=cell_geometry_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(VelocityBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t row=id[0];
        for(std::size_t k=row_offsets[row];k<row_offsets[row+1U];++k){
            const std::size_t col=columns[k];double a=0.0;
            if(col==row){
                a=cg[3U*nc+row]/dt;
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];const double sign=o==row?1.0:-1.0;
                    const double phi=sign*flux[f],diff=nu*fg[8U*nf+f];
                    if(n!=invalid){a+=diff+(include_convection&&phi>0.0?phi:0.0);}
                    else if(boundary_kind[f]==fixed){a+=diff+(include_convection&&phi>0.0?phi:0.0);}
                    else if(include_convection){a+=phi;}
                }
                if(!(a>0.0)) a=cg[3U*nc+row]/dt;
            } else {
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n==invalid) continue;
                    const std::size_t other=o==row?n:o;if(other!=col) continue;
                    const double sign=o==row?1.0:-1.0;const double phi=sign*flux[f];
                    const double entering=include_convection&&phi<0.0?-phi:0.0;
                    a=-(nu*fg[8U*nf+f]+entering);break;
                }
            }
            values[k]=a;
        }
    });
}

void ResidentPolyMeshSycl::assemble_scalar_transport_system(
    const std::size_t* row_offsets,const std::size_t* columns,std::size_t nonzeros,
    const double* old_scalar,const double* flux,const std::uint8_t* boundary_kind,
    const double* boundary_fixed,double diffusivity,double dt,const double* source,
    double* values,double* rhs) const {
    require_pointer(row_offsets,"null resident scalar row offsets");require_pointer(columns,"null resident scalar columns");
    require_pointer(old_scalar,"null resident scalar old field");require_pointer(flux,"null resident scalar face flux");
    require_pointer(boundary_kind,"null resident scalar boundary kind");require_pointer(boundary_fixed,"null resident scalar boundary values");
    require_pointer(values,"null resident scalar matrix values");require_pointer(rhs,"null resident scalar RHS");
    if(nonzeros==0U||!(diffusivity>=0.0)||!std::isfinite(diffusivity)||!(dt>0.0)||!std::isfinite(dt))
        throw std::invalid_argument("invalid resident scalar transport controls");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* cg=cell_geometry_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t row=id[0];const double volume=cg[3U*nc+row];
        double b=(volume/dt)*old_scalar[row]+(source?volume*source[row]:0.0);
        for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
            const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n!=invalid) continue;
            const double sign=o==row?1.0:-1.0;const double phi=sign*flux[f];const double diff=diffusivity*fg[8U*nf+f];
            if(boundary_kind[f]==fixed){const double entering=phi<0.0?-phi:0.0;b+=(diff+entering)*boundary_fixed[f];}
        }
        rhs[row]=b;
        for(std::size_t k=row_offsets[row];k<row_offsets[row+1U];++k){
            const std::size_t col=columns[k];double a=0.0;
            if(col==row){
                a=volume/dt;
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];const double sign=o==row?1.0:-1.0;
                    const double phi=sign*flux[f],diff=diffusivity*fg[8U*nf+f];
                    if(n!=invalid){a+=diff+(phi>0.0?phi:0.0);}
                    else if(boundary_kind[f]==fixed){a+=diff+(phi>0.0?phi:0.0);}
                    else {a+=phi>0.0?phi:0.0;}
                }
            } else {
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n==invalid) continue;
                    const std::size_t other=o==row?n:o;if(other!=col) continue;
                    const double sign=o==row?1.0:-1.0;const double phi=sign*flux[f];
                    a=-(diffusivity*fg[8U*nf+f]+(phi<0.0?-phi:0.0));break;
                }
            }
            values[k]=a;
        }
    });
}

void ResidentPolyMeshSycl::assemble_scalar_transport_system_variable(
    const std::size_t* row_offsets,const std::size_t* columns,std::size_t nonzeros,
    const double* old_scalar,const double* flux,const std::uint8_t* boundary_kind,
    const double* boundary_fixed,const double* cell_diffusivity,double dt,const double* source,
    const double* implicit_sink,double* values,double* rhs) const {
    require_pointer(row_offsets,"null resident variable-scalar row offsets");
    require_pointer(columns,"null resident variable-scalar columns");
    require_pointer(old_scalar,"null resident variable-scalar old field");
    require_pointer(flux,"null resident variable-scalar face flux");
    require_pointer(boundary_kind,"null resident variable-scalar boundary kind");
    require_pointer(boundary_fixed,"null resident variable-scalar boundary values");
    require_pointer(cell_diffusivity,"null resident variable-scalar diffusivity");
    require_pointer(values,"null resident variable-scalar matrix values");
    require_pointer(rhs,"null resident variable-scalar RHS");
    if(nonzeros==0U||!(dt>0.0)||!std::isfinite(dt))
        throw std::invalid_argument("invalid resident variable-scalar controls");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;
    const auto* face_indices=cell_face_indices_;const auto* cg=cell_geometry_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t row=id[0];const double volume=cg[3U*nc+row];
        double b=(volume/dt)*old_scalar[row]+(source?volume*source[row]:0.0);
        for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
            const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n!=invalid) continue;
            const double sign=o==row?1.0:-1.0;const double phi=sign*flux[f];
            const double diff=cell_diffusivity[row]*fg[8U*nf+f];
            if(boundary_kind[f]==fixed){const double entering=phi<0.0?-phi:0.0;b+=(diff+entering)*boundary_fixed[f];}
        }
        rhs[row]=b;
        for(std::size_t k=row_offsets[row];k<row_offsets[row+1U];++k){
            const std::size_t col=columns[k];double a=0.0;
            if(col==row){
                a=volume/dt+(implicit_sink?volume*(implicit_sink[row]>0.0?implicit_sink[row]:0.0):0.0);
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];const double sign=o==row?1.0:-1.0;
                    const double phi=sign*flux[f];double gamma=cell_diffusivity[row];
                    if(n!=invalid){const double wo=fg[6U*nf+f];gamma=wo*cell_diffusivity[o]+(1.0-wo)*cell_diffusivity[n];}
                    const double diff=(gamma>0.0?gamma:0.0)*fg[8U*nf+f];
                    if(n!=invalid){a+=diff+(phi>0.0?phi:0.0);}
                    else if(boundary_kind[f]==fixed){a+=diff+(phi>0.0?phi:0.0);}
                    else {a+=phi>0.0?phi:0.0;}
                }
            } else {
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n==invalid) continue;
                    const std::size_t other=o==row?n:o;if(other!=col) continue;
                    const double sign=o==row?1.0:-1.0;const double phi=sign*flux[f];const double wo=fg[6U*nf+f];
                    const double gamma=wo*cell_diffusivity[o]+(1.0-wo)*cell_diffusivity[n];
                    a=-((gamma>0.0?gamma:0.0)*fg[8U*nf+f]+(phi<0.0?-phi:0.0));break;
                }
            }
            values[k]=a;
        }
    });
}

void ResidentPolyMeshSycl::clamp_scalar(double* field,double minimum,double maximum) const {
    require_pointer(field,"null resident scalar clamp field");
    if(!std::isfinite(minimum)||!std::isfinite(maximum)||minimum>maximum)
        throw std::invalid_argument("invalid resident scalar clamp bounds");
    const std::size_t nc=cell_count_;
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t c=id[0];const double v=field[c];field[c]=v<minimum?minimum:(v>maximum?maximum:v);
    });
}

void ResidentPolyMeshSycl::velocity_strain_rate_magnitude(const double* velocity,double* gradient,double* strain) const {
    require_pointer(velocity,"null resident velocity for strain rate");require_pointer(gradient,"null resident velocity-gradient scratch");
    require_pointer(strain,"null resident strain-rate output");const std::size_t nc=cell_count_;
    gauss_gradient_scalar(velocity,gradient,nullptr);
    gauss_gradient_scalar(velocity+nc,gradient+3U*nc,nullptr);
    gauss_gradient_scalar(velocity+2U*nc,gradient+6U*nc,nullptr);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t c=id[0];
        const double dux=gradient[c],duy=gradient[nc+c],duz=gradient[2U*nc+c];
        const double dvx=gradient[3U*nc+c],dvy=gradient[4U*nc+c],dvz=gradient[5U*nc+c];
        const double dwx=gradient[6U*nc+c],dwy=gradient[7U*nc+c],dwz=gradient[8U*nc+c];
        const double sxy=0.5*(duy+dvx),sxz=0.5*(duz+dwx),syz=0.5*(dvz+dwy);
        const double ss=dux*dux+dvy*dvy+dwz*dwz+2.0*(sxy*sxy+sxz*sxz+syz*syz);
        strain[c]=sycl::sqrt(2.0*ss);
    });
}

void ResidentPolyMeshSycl::form_h_by_a(const double* velocity,const double* mobility,const double* pressure_gradient,double* h) const {
    require_pointer(velocity,"null resident FVM momentum velocity");require_pointer(mobility,"null resident FVM momentum mobility");
    require_pointer(pressure_gradient,"null resident FVM momentum pressure gradient");require_pointer(h,"null resident FVM H/A output");
    const std::size_t nc=cell_count_;
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const std::size_t c=id[0],cy=nc+c,cz=2U*nc+c;const double m=mobility[c];
        h[c]=velocity[c]+m*pressure_gradient[c];h[cy]=velocity[cy]+m*pressure_gradient[cy];h[cz]=velocity[cz]+m*pressure_gradient[cz];});
}

void ResidentPolyMeshSycl::assemble_pinned_pressure_values(const std::size_t* row_offsets,const std::size_t* columns,
                                                             std::size_t nonzeros,const double* mobility,std::size_t reference_cell,
                                                             double* values,const std::uint8_t* boundary_pressure_kind,
                                                             bool pin_reference) const {
    require_pointer(row_offsets,"null resident pressure row offsets");require_pointer(columns,"null resident pressure columns");
    require_pointer(mobility,"null resident pressure mobility for assembly");require_pointer(values,"null resident pressure matrix values");
    if(nonzeros==0U||(pin_reference&&reference_cell>=cell_count_)) throw std::invalid_argument("invalid resident pressure matrix assembly controls");
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* offsets=cell_face_offsets_;const auto* face_indices=cell_face_indices_;
    const auto* fg=face_geometry_;const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell,ref=reference_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){
        const std::size_t row=id[0];
        for(std::size_t k=row_offsets[row];k<row_offsets[row+1U];++k){
            const std::size_t col=columns[k];double a=0.0;
            if(pin_reference&&row==ref){a=col==row?1.0:0.0;values[k]=a;continue;}
            if(col==row){
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];
                    if(n==invalid){
                        if(boundary_pressure_kind&&boundary_pressure_kind[f]==fixed) a+=mobility[o]*fg[8U*nf+f];
                        continue;
                    }
                    const double wo=fg[6U*nf+f];const double mf=wo*mobility[o]+(1.0-wo)*mobility[n];a+=mf*fg[8U*nf+f];
                }
            } else {
                for(std::size_t j=offsets[row];j<offsets[row+1U];++j){
                    const std::size_t f=face_indices[j],o=owner[f],n=neighbour[f];if(n==invalid) continue;
                    const std::size_t other=o==row?n:o;if(other==col){
                        if(pin_reference&&col==ref){a=0.0;break;}
                        const double wo=fg[6U*nf+f];const double mf=wo*mobility[o]+(1.0-wo)*mobility[n];a=-mf*fg[8U*nf+f];break;
                    }
                }
            }
            values[k]=a;
        }
    });
}

void ResidentPolyMeshSycl::nonorthogonal_pressure_flux(const double* pressure,const double* mobility,
                                                         double* gradient,double* correction,
                                                         const double* boundary_pressure,
                                                         const std::uint8_t* boundary_pressure_kind) const {
    require_pointer(pressure,"null resident pressure for nonorthogonal flux");require_pointer(mobility,"null resident pressure mobility for nonorthogonal flux");
    require_pointer(gradient,"null resident pressure gradient scratch");require_pointer(correction,"null resident nonorthogonal pressure flux");
    gauss_gradient_scalar(pressure,gradient,boundary_pressure);
    const auto* owner=face_owner_;const auto* neighbour=face_neighbour_;const auto* cg=cell_geometry_;const auto* fg=face_geometry_;
    const std::size_t nc=cell_count_,nf=face_count_,invalid=invalid_cell;
    const std::uint8_t fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    queue_.parallel_for(sycl::range<1>(nf),[=](sycl::id<1> id){
        const std::size_t f=id[0],o=owner[f],n=neighbour[f];
        double mf=mobility[o],gx=gradient[o],gy=gradient[nc+o],gz=gradient[2U*nc+o];
        double dx=fg[f]-cg[o],dy=fg[nf+f]-cg[nc+o],dz=fg[2U*nf+f]-cg[2U*nc+o];
        if(n!=invalid){
            const double wo=fg[6U*nf+f];mf=wo*mobility[o]+(1.0-wo)*mobility[n];
            gx=wo*gradient[o]+(1.0-wo)*gradient[n];gy=wo*gradient[nc+o]+(1.0-wo)*gradient[nc+n];gz=wo*gradient[2U*nc+o]+(1.0-wo)*gradient[2U*nc+n];
            dx=cg[n]-cg[o];dy=cg[nc+n]-cg[nc+o];dz=cg[2U*nc+n]-cg[2U*nc+o];
        } else if(!boundary_pressure||!boundary_pressure_kind||boundary_pressure_kind[f]!=fixed){correction[f]=0.0;return;}
        const double metric=fg[8U*nf+f];const double nx=fg[3U*nf+f]-metric*dx,ny=fg[4U*nf+f]-metric*dy,nz=fg[5U*nf+f]-metric*dz;
        correction[f]=mf*(gx*nx+gy*ny+gz*nz);
    });
}

void ResidentPolyMeshSycl::volume_weighted_vector_norm2(const double* values,double* result) const {
    require_pointer(values,"null resident FVM vector norm field");require_pointer(result,"null resident FVM vector norm reduction");
    *result=0.0;const auto* cg=cell_geometry_;const std::size_t nc=cell_count_;
    auto reduction=sycl::reduction(result,sycl::plus<double>());
    queue_.parallel_for(sycl::range<1>(nc),reduction,[=](sycl::id<1> id,auto& sum){const std::size_t c=id[0];
        const double x=values[c],y=values[nc+c],z=values[2U*nc+c];sum.combine((x*x+y*y+z*z)*cg[3U*nc+c]);}).wait_and_throw();
    transfer_stats_.record_synchronization();
}

ResidentPressureProjectionSycl::ResidentPressureProjectionSycl(
    const PolyMesh& mesh, double mobility, std::size_t reference_cell, sycl::device device)
    : ResidentPressureProjectionSycl(mesh, mobility, reference_cell,
                                     sycl::queue(device, sycl::property::queue::in_order{})) {}

ResidentPressureProjectionSycl::ResidentPressureProjectionSycl(
    const PolyMesh& mesh, double mobility, std::size_t reference_cell, sycl::queue queue)
    : mesh_(mesh, std::move(queue)), reference_cell_(reference_cell) {
    if (reference_cell_ >= mesh.cell_count()) throw std::out_of_range("resident pressure reference cell out of range");
    const auto matrix = build_orthogonal_pressure_matrix(mesh, mobility, false, reference_cell_);
    pressure_solver_ = std::make_unique<cfd::core::SyclCsrLinearAlgebra>(matrix, mesh_.queue());
    allocate(mobility);
}

ResidentPressureProjectionSycl::~ResidentPressureProjectionSycl() noexcept { release(); }

void ResidentPressureProjectionSycl::allocate(double mobility_value) {
    if (!(mobility_value > 0.0) || !std::isfinite(mobility_value))
        throw std::invalid_argument("resident pressure mobility must be finite and positive");
    auto& q=mesh_.queue(); const std::size_t nc=mesh_.cell_count(), nf=mesh_.face_count();
    h_by_a_=sycl::malloc_device<double>(3U*nc,q); pressure_=sycl::malloc_device<double>(nc,q);
    pressure_gradient_=sycl::malloc_device<double>(3U*nc,q); velocity_=sycl::malloc_device<double>(3U*nc,q);
    face_flux_=sycl::malloc_device<double>(nf,q); rhs_=sycl::malloc_device<double>(nc,q);
    mobility_=sycl::malloc_device<double>(nc,q); boundary_velocity_=sycl::malloc_device<double>(3U*nf,q);
    continuity_=sycl::malloc_device<double>(nc,q); reduction_scalar_=sycl::malloc_shared<double>(1U,q);
    if(!h_by_a_||!pressure_||!pressure_gradient_||!velocity_||!face_flux_||!rhs_||!mobility_||
       !boundary_velocity_||!continuity_||!reduction_scalar_){release();throw std::bad_alloc{};}
    mesh_.fill(h_by_a_,3U*nc,0.0);mesh_.fill(pressure_,nc,0.0);mesh_.fill(pressure_gradient_,3U*nc,0.0);
    mesh_.fill(velocity_,3U*nc,0.0);mesh_.fill(face_flux_,nf,0.0);mesh_.fill(rhs_,nc,0.0);
    mesh_.fill(mobility_,nc,mobility_value);mesh_.fill(boundary_velocity_,3U*nf,0.0);mesh_.fill(continuity_,nc,0.0);
    mesh_.wait(); *reduction_scalar_=0.0;
}

void ResidentPressureProjectionSycl::release() noexcept {
    if(!pressure_solver_ && !h_by_a_ && !pressure_ && !pressure_gradient_ && !velocity_ && !face_flux_ &&
       !rhs_ && !mobility_ && !boundary_velocity_ && !continuity_ && !reduction_scalar_) return;
    auto& q=mesh_.queue(); try{q.wait_and_throw();}catch(...){}
    free_if(h_by_a_,q);free_if(pressure_,q);free_if(pressure_gradient_,q);free_if(velocity_,q);
    free_if(face_flux_,q);free_if(rhs_,q);free_if(mobility_,q);free_if(boundary_velocity_,q);
    free_if(continuity_,q);free_if(reduction_scalar_,q);pressure_solver_.reset();
}

void ResidentPressureProjectionSycl::set_predictor(std::span<const Vec3> velocity) { mesh_.upload_cell_vector(velocity,h_by_a_); }
void ResidentPressureProjectionSycl::set_pressure(std::span<const double> pressure) {
    mesh_.upload_cell_scalar(pressure,pressure_);mesh_.set_cell_value(pressure_,reference_cell_,0.0);mesh_.wait();
}

cfd::core::IterativeSolverResult ResidentPressureProjectionSycl::project(
    std::size_t max_iterations,double relative_tolerance) {
    mesh_.predictor_face_flux(h_by_a_,face_flux_,boundary_velocity_);
    mesh_.pressure_rhs_from_flux(face_flux_,mobility_,rhs_);
    mesh_.set_cell_value(rhs_,reference_cell_,0.0);mesh_.set_cell_value(pressure_,reference_cell_,0.0);
    const auto result=pressure_solver_->conjugate_gradient_device(rhs_,pressure_,max_iterations,relative_tolerance);
    mesh_.set_cell_value(pressure_,reference_cell_,0.0);
    mesh_.gauss_gradient_scalar(pressure_,pressure_gradient_);
    mesh_.correct_velocity_from_pressure_gradient(h_by_a_,mobility_,pressure_gradient_,velocity_);
    mesh_.pressure_correct_face_flux(pressure_,mobility_,face_flux_);
    return result;
}

double ResidentPressureProjectionSycl::continuity_l2() {
    mesh_.divergence_face_flux(face_flux_,continuity_);
    *reduction_scalar_=0.0;const std::size_t nc=mesh_.cell_count();const double* div=continuity_;
    auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<double>());
    mesh_.queue().parallel_for(sycl::range<1>(nc),reduction,[=](sycl::id<1> id,auto& sum){const double v=div[id[0]];sum.combine(v*v);}).wait_and_throw();
    return std::sqrt(*reduction_scalar_/static_cast<double>(nc));
}

void ResidentPressureProjectionSycl::download_velocity(std::span<Vec3> velocity) const { mesh_.download_cell_vector(velocity_,velocity); }
void ResidentPressureProjectionSycl::download_pressure(std::span<double> pressure) const { mesh_.download_cell_scalar(pressure_,pressure); }
void ResidentPressureProjectionSycl::reset_transfer_stats() const noexcept {
    mesh_.reset_transfer_stats(); if(pressure_solver_) pressure_solver_->reset_transfer_stats();
}
std::uint64_t ResidentPressureProjectionSycl::hot_loop_host_transfer_bytes() const noexcept {
    const std::uint64_t mesh_bytes=mesh_.transfer_stats().host_transfer_bytes();
    const std::uint64_t solver_bytes=pressure_solver_?pressure_solver_->transfer_stats().host_transfer_bytes():0U;
    return mesh_bytes+solver_bytes;
}
std::size_t ResidentPressureProjectionSycl::resident_bytes() const noexcept {
    const std::size_t nc=mesh_.cell_count(),nf=mesh_.face_count();
    return mesh_.resident_bytes()+(3U*nc+nc+3U*nc+3U*nc+nf+nc+nc+3U*nf+nc)*sizeof(double)+sizeof(double);
}

ResidentIncompressibleSycl::ResidentIncompressibleSycl(
    const PolyMesh& mesh, ResidentIncompressibleConfig config, std::size_t reference_cell, sycl::device device)
    : ResidentIncompressibleSycl(mesh, config, reference_cell,
                                 sycl::queue(device, sycl::property::queue::in_order{})) {}

ResidentIncompressibleSycl::ResidentIncompressibleSycl(
    const PolyMesh& mesh, ResidentIncompressibleConfig config, std::size_t reference_cell, sycl::queue queue)
    : mesh_(mesh, std::move(queue)), config_(config), reference_cell_(reference_cell),
      patch_names_(mesh.patches().size()), patch_faces_(mesh.patches().size()),
      boundary_kind_host_(mesh.face_count(), static_cast<std::uint8_t>(VelocityBoundaryType::zeroGradient)),
      boundary_fixed_host_(3U * mesh.face_count(), 0.0),
      pressure_boundary_kind_host_(mesh.face_count(), static_cast<std::uint8_t>(PressureBoundaryType::zeroGradient)),
      pressure_boundary_fixed_host_(mesh.face_count(), 0.0) {
    if(reference_cell_>=mesh.cell_count()) throw std::out_of_range("resident incompressible pressure reference cell out of range");
    if(!(config_.density>0.0)||!std::isfinite(config_.density)||!(config_.kinematic_viscosity>=0.0)||!std::isfinite(config_.kinematic_viscosity)||
       !(config_.dt>0.0)||!std::isfinite(config_.dt)||config_.momentum_sweeps==0U||config_.momentum_iterations==0U||
       !(config_.momentum_tolerance>0.0)||!std::isfinite(config_.momentum_tolerance)||config_.pressure_iterations==0U||
       !(config_.pressure_tolerance>0.0)||config_.pressure_correctors==0U||config_.outer_correctors==0U||
       !(config_.velocity_relaxation>0.0&&config_.velocity_relaxation<=1.0)||
       !(config_.pressure_relaxation>0.0&&config_.pressure_relaxation<=1.0)) {
        throw std::invalid_argument("invalid resident incompressible controls");
    }
    for(std::size_t p=0U;p<mesh.patches().size();++p) patch_names_[p]=mesh.patches()[p].name;
    for(std::size_t f=0U;f<mesh.face_count();++f) if(mesh.faces()[f].boundary()) patch_faces_[mesh.faces()[f].patch].push_back(f);
    const auto pattern=build_cell_neighbour_pattern(mesh);
    pressure_solver_=std::make_unique<cfd::core::SyclCsrLinearAlgebra>(pattern,mesh_.queue());
    momentum_solver_=std::make_unique<cfd::core::SyclCsrLinearAlgebra>(pattern,mesh_.queue());
    allocate();upload_boundary_state();initialize_uniform();
}

ResidentIncompressibleSycl::~ResidentIncompressibleSycl() noexcept { release(); }

void ResidentIncompressibleSycl::allocate() {
    auto& q=mesh_.queue();const std::size_t nc=mesh_.cell_count(),nf=mesh_.face_count(),nnz=pressure_solver_->nonzeros();
    velocity_=sycl::malloc_device<double>(3U*nc,q);old_velocity_=sycl::malloc_device<double>(3U*nc,q);
    source_velocity_=sycl::malloc_device<double>(3U*nc,q);velocity_before_=sycl::malloc_device<double>(3U*nc,q);
    velocity_work_=sycl::malloc_device<double>(3U*nc,q);h_by_a_=sycl::malloc_device<double>(3U*nc,q);
    pressure_=sycl::malloc_device<double>(nc,q);pressure_before_=sycl::malloc_device<double>(nc,q);
    pressure_gradient_=sycl::malloc_device<double>(3U*nc,q);face_flux_=sycl::malloc_device<double>(nf,q);
    predictor_flux_=sycl::malloc_device<double>(nf,q);nonorthogonal_flux_=sycl::malloc_device<double>(nf,q);
    rhs_=sycl::malloc_device<double>(nc,q);mobility_=sycl::malloc_device<double>(nc,q);
    momentum_diagonal_=sycl::malloc_device<double>(nc,q);momentum_rhs_=sycl::malloc_device<double>(3U*nc,q);
    boundary_velocity_=sycl::malloc_device<double>(3U*nf,q);boundary_kind_=sycl::malloc_device<std::uint8_t>(nf,q);
    boundary_fixed_=sycl::malloc_device<double>(3U*nf,q);
    pressure_boundary_kind_=sycl::malloc_device<std::uint8_t>(nf,q);
    pressure_boundary_fixed_=sycl::malloc_device<double>(nf,q);
    pressure_boundary_values_=sycl::malloc_device<double>(nf,q);
    pressure_matrix_values_=sycl::malloc_device<double>(nnz,q);
    momentum_matrix_values_=sycl::malloc_device<double>(momentum_solver_->nonzeros(),q);
    continuity_=sycl::malloc_device<double>(nc,q);reduction_scalar_=sycl::malloc_shared<double>(1U,q);
    if(!velocity_||!old_velocity_||!source_velocity_||!velocity_before_||!velocity_work_||!h_by_a_||!pressure_||
       !pressure_before_||!pressure_gradient_||!face_flux_||!predictor_flux_||!nonorthogonal_flux_||!rhs_||!mobility_||
       !momentum_diagonal_||!momentum_rhs_||!boundary_velocity_||!boundary_kind_||!boundary_fixed_||
       !pressure_boundary_kind_||!pressure_boundary_fixed_||!pressure_boundary_values_||!pressure_matrix_values_||
       !momentum_matrix_values_||!continuity_||!reduction_scalar_) { release();throw std::bad_alloc{}; }
    mesh_.fill(velocity_,3U*nc,0.0);mesh_.fill(old_velocity_,3U*nc,0.0);mesh_.fill(source_velocity_,3U*nc,0.0);
    mesh_.fill(velocity_before_,3U*nc,0.0);mesh_.fill(velocity_work_,3U*nc,0.0);mesh_.fill(h_by_a_,3U*nc,0.0);
    mesh_.fill(pressure_,nc,0.0);mesh_.fill(pressure_before_,nc,0.0);mesh_.fill(pressure_gradient_,3U*nc,0.0);
    mesh_.fill(face_flux_,nf,0.0);mesh_.fill(predictor_flux_,nf,0.0);mesh_.fill(nonorthogonal_flux_,nf,0.0);
    mesh_.fill(rhs_,nc,0.0);mesh_.fill(mobility_,nc,config_.dt);mesh_.fill(momentum_diagonal_,nc,1.0);
    mesh_.fill(momentum_rhs_,3U*nc,0.0);mesh_.fill(boundary_velocity_,3U*nf,0.0);
    mesh_.fill(pressure_boundary_fixed_,nf,0.0);mesh_.fill(pressure_boundary_values_,nf,0.0);mesh_.fill(continuity_,nc,0.0);
    mesh_.wait();*reduction_scalar_=0.0;
}

void ResidentIncompressibleSycl::release() noexcept {
    auto& q=mesh_.queue();try{q.wait_and_throw();}catch(...){}
    free_if(velocity_,q);free_if(old_velocity_,q);free_if(source_velocity_,q);free_if(velocity_before_,q);free_if(velocity_work_,q);
    free_if(h_by_a_,q);free_if(pressure_,q);free_if(pressure_before_,q);free_if(pressure_gradient_,q);free_if(face_flux_,q);
    free_if(predictor_flux_,q);free_if(nonorthogonal_flux_,q);free_if(rhs_,q);free_if(mobility_,q);free_if(momentum_diagonal_,q);
    free_if(momentum_rhs_,q);free_if(boundary_velocity_,q);free_if(boundary_kind_,q);free_if(boundary_fixed_,q);
    free_if(pressure_boundary_kind_,q);free_if(pressure_boundary_fixed_,q);free_if(pressure_boundary_values_,q);
    free_if(pressure_matrix_values_,q);free_if(momentum_matrix_values_,q);free_if(continuity_,q);free_if(reduction_scalar_,q);
    momentum_solver_.reset();pressure_solver_.reset();
}

void ResidentIncompressibleSycl::upload_boundary_state() {
    auto& q=mesh_.queue();const std::size_t nf=mesh_.face_count();
    const std::size_t kb=nf*sizeof(std::uint8_t),vb=3U*nf*sizeof(double),pb=nf*sizeof(double);
    q.memcpy(boundary_kind_,boundary_kind_host_.data(),kb);q.memcpy(boundary_fixed_,boundary_fixed_host_.data(),vb);
    q.memcpy(pressure_boundary_kind_,pressure_boundary_kind_host_.data(),kb);
    q.memcpy(pressure_boundary_fixed_,pressure_boundary_fixed_host_.data(),pb).wait_and_throw();
    transfer_stats_.record_host_to_device(2U*kb+vb+pb);transfer_stats_.record_synchronization();
}

void ResidentIncompressibleSycl::set_velocity_boundary(std::string_view patch,VelocityBoundaryType type,Vec3 value) {
    std::size_t pi=patch_names_.size();for(std::size_t p=0U;p<patch_names_.size();++p) if(patch_names_[p]==patch){pi=p;break;}
    if(pi==patch_names_.size()) throw std::out_of_range("resident incompressible patch not found");
    const std::size_t nf=mesh_.face_count();const auto kind=static_cast<std::uint8_t>(type);
    for(const std::size_t f:patch_faces_[pi]){boundary_kind_host_[f]=kind;boundary_fixed_host_[f]=value.x;boundary_fixed_host_[nf+f]=value.y;boundary_fixed_host_[2U*nf+f]=value.z;}
    upload_boundary_state();rebuild_face_flux();
}

void ResidentIncompressibleSycl::set_pressure_boundary(std::string_view patch,PressureBoundaryType type,double value) {
    if(!std::isfinite(value)) throw std::invalid_argument("resident incompressible pressure boundary value must be finite");
    std::size_t pi=patch_names_.size();for(std::size_t p=0U;p<patch_names_.size();++p) if(patch_names_[p]==patch){pi=p;break;}
    if(pi==patch_names_.size()) throw std::out_of_range("resident incompressible pressure patch not found");
    const auto kind=static_cast<std::uint8_t>(type);
    for(const std::size_t f:patch_faces_[pi]){pressure_boundary_kind_host_[f]=kind;pressure_boundary_fixed_host_[f]=value;}
    const auto fixed=static_cast<std::uint8_t>(PressureBoundaryType::fixedValue);
    has_fixed_pressure_boundary_=std::any_of(pressure_boundary_kind_host_.begin(),pressure_boundary_kind_host_.end(),
                                             [fixed](std::uint8_t k){return k==fixed;});
    upload_boundary_state();
}

void ResidentIncompressibleSycl::initialize_uniform(Vec3 velocity,double pressure) {
    if(!std::isfinite(pressure)) throw std::invalid_argument("resident incompressible initial pressure must be finite");
    auto& q=mesh_.queue();const std::size_t nc=mesh_.cell_count();double* u=velocity_;double* p=pressure_;
    q.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const std::size_t c=id[0];u[c]=velocity.x;u[nc+c]=velocity.y;u[2U*nc+c]=velocity.z;p[c]=pressure;});
    if(!has_fixed_pressure_boundary_) mesh_.set_cell_value(pressure_,reference_cell_,0.0);
    q.memcpy(old_velocity_,velocity_,3U*nc*sizeof(double));q.memcpy(source_velocity_,velocity_,3U*nc*sizeof(double));
    q.wait_and_throw();transfer_stats_.record_device_to_device(6U*nc*sizeof(double));transfer_stats_.record_synchronization();
    time_=0.0;steps_=0U;rebuild_face_flux();
}

void ResidentIncompressibleSycl::rebuild_face_flux() {
    mesh_.boundary_velocity_values(velocity_,boundary_kind_,boundary_fixed_,boundary_velocity_);
    mesh_.predictor_face_flux(velocity_,face_flux_,boundary_velocity_);
}

void ResidentIncompressibleSycl::momentum_predictor(const double* time_source) {
    mesh_.gauss_gradient_scalar(pressure_,pressure_gradient_);
    mesh_.boundary_velocity_values(velocity_,boundary_kind_,boundary_fixed_,boundary_velocity_);
    mesh_.assemble_momentum_system(time_source,pressure_gradient_,face_flux_,boundary_kind_,boundary_velocity_,
                                   config_.kinematic_viscosity,config_.dt,config_.include_convection,
                                   momentum_diagonal_,mobility_,momentum_rhs_);
    mesh_.assemble_momentum_matrix_values(momentum_solver_->row_offsets_device(),momentum_solver_->column_indices_device(),
                                          momentum_solver_->nonzeros(),face_flux_,boundary_kind_,
                                          config_.kinematic_viscosity,config_.dt,config_.include_convection,
                                          momentum_matrix_values_);
    momentum_solver_->update_values_device(momentum_matrix_values_);

    const std::size_t nc=mesh_.cell_count();
    bool krylov_ok=true;
    for(std::size_t component=0U;component<3U;++component){
        const auto result=momentum_solver_->bicgstab_device(momentum_rhs_+component*nc,velocity_+component*nc,
                                                             config_.momentum_iterations,config_.momentum_tolerance);
        krylov_ok=krylov_ok&&result.converged;
    }
    if(!krylov_ok){
        // Robustness fallback for pathological startup matrices. The primary
        // path is BiCGStab; Jacobi remains useful for recovery/debugging.
        double* current=velocity_;double* next=velocity_work_;
        for(std::size_t sweep=0U;sweep<config_.momentum_sweeps;++sweep){
            mesh_.momentum_jacobi_sweep(current,face_flux_,boundary_kind_,config_.kinematic_viscosity,config_.include_convection,
                                        momentum_diagonal_,momentum_rhs_,next);
            std::swap(current,next);
        }
        if(current!=velocity_){const std::size_t bytes=3U*nc*sizeof(double);mesh_.queue().memcpy(velocity_,current,bytes);transfer_stats_.record_device_to_device(bytes);}
    }
    mesh_.form_h_by_a(velocity_,mobility_,pressure_gradient_,h_by_a_);
}

cfd::core::IterativeSolverResult ResidentIncompressibleSycl::correct_pressure(double pressure_relaxation) {
    auto& q=mesh_.queue();const std::size_t nc=mesh_.cell_count();
    const std::size_t pbytes=nc*sizeof(double);q.memcpy(pressure_before_,pressure_,pbytes);transfer_stats_.record_device_to_device(pbytes);
    mesh_.boundary_velocity_values(h_by_a_,boundary_kind_,boundary_fixed_,boundary_velocity_);
    mesh_.predictor_face_flux(h_by_a_,predictor_flux_,boundary_velocity_);
    mesh_.boundary_pressure_values(pressure_,pressure_boundary_kind_,pressure_boundary_fixed_,pressure_boundary_values_);
    const bool pin_reference=!has_fixed_pressure_boundary_;
    mesh_.assemble_pinned_pressure_values(pressure_solver_->row_offsets_device(),pressure_solver_->column_indices_device(),
                                          pressure_solver_->nonzeros(),mobility_,reference_cell_,pressure_matrix_values_,
                                          pressure_boundary_kind_,pin_reference);
    pressure_solver_->update_values_device(pressure_matrix_values_);
    cfd::core::IterativeSolverResult result{};
    const std::size_t loops=std::max<std::size_t>(1U,config_.nonorthogonal_correctors+1U);
    for(std::size_t nonorth=0U;nonorth<loops;++nonorth){
        mesh_.boundary_pressure_values(pressure_,pressure_boundary_kind_,pressure_boundary_fixed_,pressure_boundary_values_);
        mesh_.nonorthogonal_pressure_flux(pressure_,mobility_,pressure_gradient_,nonorthogonal_flux_,
                                          pressure_boundary_values_,pressure_boundary_kind_);
        mesh_.pressure_rhs_from_flux(predictor_flux_,mobility_,rhs_,pressure_boundary_values_,pressure_boundary_kind_,nonorthogonal_flux_);
        if(pin_reference){mesh_.set_cell_value(rhs_,reference_cell_,0.0);mesh_.set_cell_value(pressure_,reference_cell_,0.0);}
        result=pressure_solver_->conjugate_gradient_device(rhs_,pressure_,config_.pressure_iterations,config_.pressure_tolerance);
        if(!result.converged) break;
    }
    if(result.converged){
        double* p=pressure_;const double* oldp=pressure_before_;const double relax=pressure_relaxation;
        q.parallel_for(sycl::range<1>(nc),[=](sycl::id<1> id){const std::size_t c=id[0];p[c]=oldp[c]+relax*(p[c]-oldp[c]);});
        if(pin_reference) mesh_.set_cell_value(pressure_,reference_cell_,0.0);
        mesh_.boundary_pressure_values(pressure_,pressure_boundary_kind_,pressure_boundary_fixed_,pressure_boundary_values_);
        mesh_.gauss_gradient_scalar(pressure_,pressure_gradient_,pressure_boundary_values_);
        mesh_.correct_velocity_from_pressure_gradient(h_by_a_,mobility_,pressure_gradient_,velocity_work_);
        double* u=velocity_;const double* target=velocity_work_;const double urelax=pressure_relaxation<1.0?config_.velocity_relaxation:1.0;
        q.parallel_for(sycl::range<1>(3U*nc),[=](sycl::id<1> id){const std::size_t i=id[0];u[i]+=urelax*(target[i]-u[i]);});
        mesh_.boundary_velocity_values(h_by_a_,boundary_kind_,boundary_fixed_,boundary_velocity_);
        mesh_.predictor_face_flux(h_by_a_,face_flux_,boundary_velocity_);
        mesh_.nonorthogonal_pressure_flux(pressure_,mobility_,pressure_gradient_,nonorthogonal_flux_,
                                          pressure_boundary_values_,pressure_boundary_kind_);
        mesh_.pressure_correct_face_flux(pressure_,mobility_,face_flux_,pressure_boundary_values_,pressure_boundary_kind_,nonorthogonal_flux_);
    } else {
        q.memcpy(pressure_,pressure_before_,pbytes);transfer_stats_.record_device_to_device(pbytes);
    }
    return result;
}

double ResidentIncompressibleSycl::velocity_rms_change() {
    *reduction_scalar_=0.0;const std::size_t count=3U*mesh_.cell_count();const double* a=velocity_;const double* b=velocity_before_;
    auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<double>());
    mesh_.queue().parallel_for(sycl::range<1>(count),reduction,[=](sycl::id<1> id,auto& sum){const double d=a[id[0]]-b[id[0]];sum.combine(d*d);}).wait_and_throw();
    transfer_stats_.record_synchronization();return std::sqrt(*reduction_scalar_/static_cast<double>(mesh_.cell_count()));
}

ResidentIncompressibleIterationInfo ResidentIncompressibleSycl::coupled_sequence(
    const double* time_source,std::size_t pressure_correctors,double pressure_relaxation) {
    const std::size_t bytes=3U*mesh_.cell_count()*sizeof(double);mesh_.queue().memcpy(velocity_before_,velocity_,bytes);
    transfer_stats_.record_device_to_device(bytes);momentum_predictor(time_source);
    cfd::core::IterativeSolverResult result{};std::size_t solves=0U;
    for(std::size_t corr=0U;corr<pressure_correctors;++corr){result=correct_pressure(pressure_relaxation);++solves;if(!result.converged) break;}
    return {velocity_rms_change(),continuity_l2(),result,solves};
}

ResidentIncompressibleIterationInfo ResidentIncompressibleSycl::iterate_simple() {
    const std::size_t bytes=3U*mesh_.cell_count()*sizeof(double);mesh_.queue().memcpy(source_velocity_,velocity_,bytes);transfer_stats_.record_device_to_device(bytes);
    return coupled_sequence(source_velocity_,1U,config_.pressure_relaxation);
}

ResidentIncompressibleIterationInfo ResidentIncompressibleSycl::step_piso() {
    const std::size_t bytes=3U*mesh_.cell_count()*sizeof(double);mesh_.queue().memcpy(source_velocity_,velocity_,bytes);transfer_stats_.record_device_to_device(bytes);
    auto info=coupled_sequence(source_velocity_,config_.pressure_correctors,1.0);mesh_.queue().memcpy(old_velocity_,source_velocity_,bytes);
    transfer_stats_.record_device_to_device(bytes);time_+=config_.dt;++steps_;return info;
}

ResidentIncompressibleIterationInfo ResidentIncompressibleSycl::step_pimple() {
    const std::size_t bytes=3U*mesh_.cell_count()*sizeof(double);mesh_.queue().memcpy(source_velocity_,velocity_,bytes);transfer_stats_.record_device_to_device(bytes);
    ResidentIncompressibleIterationInfo info{};
    for(std::size_t outer=0U;outer<config_.outer_correctors;++outer){info=coupled_sequence(source_velocity_,config_.pressure_correctors,1.0);if(!info.pressure.converged) break;}
    mesh_.queue().memcpy(old_velocity_,source_velocity_,bytes);transfer_stats_.record_device_to_device(bytes);time_+=config_.dt;++steps_;return info;
}

double ResidentIncompressibleSycl::continuity_l2() {
    mesh_.divergence_face_flux(face_flux_,continuity_);*reduction_scalar_=0.0;const std::size_t nc=mesh_.cell_count();const double* div=continuity_;
    auto reduction=sycl::reduction(reduction_scalar_,sycl::plus<double>());
    mesh_.queue().parallel_for(sycl::range<1>(nc),reduction,[=](sycl::id<1> id,auto& sum){const double v=div[id[0]];sum.combine(v*v);}).wait_and_throw();
    transfer_stats_.record_synchronization();return std::sqrt(*reduction_scalar_/static_cast<double>(nc));
}

double ResidentIncompressibleSycl::kinetic_energy() {
    mesh_.volume_weighted_vector_norm2(velocity_,reduction_scalar_);return 0.5*config_.density*(*reduction_scalar_);
}

void ResidentIncompressibleSycl::download_velocity(std::span<Vec3> velocity) const {
    mesh_.download_cell_vector(velocity_,velocity);
}
void ResidentIncompressibleSycl::download_pressure(std::span<double> pressure) const {
    mesh_.download_cell_scalar(pressure_,pressure);
}
void ResidentIncompressibleSycl::reset_transfer_stats() const noexcept {
    transfer_stats_.reset();mesh_.reset_transfer_stats();if(pressure_solver_)pressure_solver_->reset_transfer_stats();
    if(momentum_solver_)momentum_solver_->reset_transfer_stats();
}
std::uint64_t ResidentIncompressibleSycl::hot_loop_host_transfer_bytes() const noexcept {
    return transfer_stats_.host_transfer_bytes()+mesh_.transfer_stats().host_transfer_bytes()+
           (pressure_solver_?pressure_solver_->transfer_stats().host_transfer_bytes():0U)+
           (momentum_solver_?momentum_solver_->transfer_stats().host_transfer_bytes():0U);
}
std::size_t ResidentIncompressibleSycl::resident_bytes() const noexcept {
    const std::size_t nc=mesh_.cell_count(),nf=mesh_.face_count(),nnz=pressure_solver_?pressure_solver_->nonzeros():0U;
    const std::size_t mnnz=momentum_solver_?momentum_solver_->nonzeros():0U;
    const std::size_t doubles=(3U*nc*8U)+(nc*6U)+(nf*3U)+(3U*nf*2U)+(2U*nf)+nnz+mnnz+nc;
    return mesh_.resident_bytes()+doubles*sizeof(double)+2U*nf*sizeof(std::uint8_t)+sizeof(double);
}

ResidentScalarTransportSycl::ResidentScalarTransportSycl(
    const PolyMesh& mesh,ResidentScalarTransportConfig config,sycl::device device)
    : ResidentScalarTransportSycl(mesh,config,sycl::queue(device,sycl::property::queue::in_order{})) {}

ResidentScalarTransportSycl::ResidentScalarTransportSycl(
    const PolyMesh& mesh,ResidentScalarTransportConfig config,sycl::queue queue)
    : mesh_(mesh,std::move(queue)),config_(config),patch_names_(mesh.patches().size()),patch_faces_(mesh.patches().size()),
      boundary_kind_host_(mesh.face_count(),static_cast<std::uint8_t>(PressureBoundaryType::zeroGradient)),
      boundary_fixed_host_(mesh.face_count(),0.0) {
    if(!(config_.dt>0.0)||!std::isfinite(config_.dt)||!(config_.diffusivity>=0.0)||!std::isfinite(config_.diffusivity)||
       config_.linear_iterations==0U||!(config_.linear_tolerance>0.0)||!std::isfinite(config_.linear_tolerance))
        throw std::invalid_argument("invalid resident scalar transport controls");
    for(std::size_t p=0U;p<mesh.patches().size();++p) patch_names_[p]=mesh.patches()[p].name;
    for(std::size_t f=0U;f<mesh.face_count();++f) if(mesh.faces()[f].boundary()) patch_faces_[mesh.faces()[f].patch].push_back(f);
    const auto pattern=build_cell_neighbour_pattern(mesh);
    solver_=std::make_unique<cfd::core::SyclCsrLinearAlgebra>(pattern,mesh_.queue());
    allocate();upload_boundary_state();initialize_uniform();
}

ResidentScalarTransportSycl::~ResidentScalarTransportSycl() noexcept { release(); }

void ResidentScalarTransportSycl::allocate() {
    auto& q=mesh_.queue();const std::size_t nc=mesh_.cell_count(),nf=mesh_.face_count();
    scalar_=sycl::malloc_device<double>(nc,q);old_scalar_=sycl::malloc_device<double>(nc,q);rhs_=sycl::malloc_device<double>(nc,q);
    matrix_values_=sycl::malloc_device<double>(solver_->nonzeros(),q);boundary_kind_=sycl::malloc_device<std::uint8_t>(nf,q);
    boundary_fixed_=sycl::malloc_device<double>(nf,q);
    if(!scalar_||!old_scalar_||!rhs_||!matrix_values_||!boundary_kind_||!boundary_fixed_){release();throw std::bad_alloc{};}
    mesh_.fill(scalar_,nc,0.0);mesh_.fill(old_scalar_,nc,0.0);mesh_.fill(rhs_,nc,0.0);mesh_.fill(boundary_fixed_,nf,0.0);mesh_.wait();
}

void ResidentScalarTransportSycl::release() noexcept {
    auto& q=mesh_.queue();try{q.wait_and_throw();}catch(...){}
    free_if(scalar_,q);free_if(old_scalar_,q);free_if(rhs_,q);free_if(matrix_values_,q);free_if(boundary_kind_,q);free_if(boundary_fixed_,q);
    solver_.reset();
}

void ResidentScalarTransportSycl::upload_boundary_state() {
    auto& q=mesh_.queue();const std::size_t nf=mesh_.face_count();const std::size_t kb=nf*sizeof(std::uint8_t),vb=nf*sizeof(double);
    q.memcpy(boundary_kind_,boundary_kind_host_.data(),kb);q.memcpy(boundary_fixed_,boundary_fixed_host_.data(),vb).wait_and_throw();
    transfer_stats_.record_host_to_device(kb+vb);transfer_stats_.record_synchronization();
}

void ResidentScalarTransportSycl::set_boundary(std::string_view patch,PressureBoundaryType type,double value) {
    if(!std::isfinite(value)) throw std::invalid_argument("resident scalar boundary value must be finite");
    std::size_t pi=patch_names_.size();for(std::size_t p=0U;p<patch_names_.size();++p) if(patch_names_[p]==patch){pi=p;break;}
    if(pi==patch_names_.size()) throw std::out_of_range("resident scalar patch not found");
    const auto kind=static_cast<std::uint8_t>(type);for(const std::size_t f:patch_faces_[pi]){boundary_kind_host_[f]=kind;boundary_fixed_host_[f]=value;}
    upload_boundary_state();
}

void ResidentScalarTransportSycl::initialize_uniform(double value) {
    if(!std::isfinite(value)) throw std::invalid_argument("resident scalar initial value must be finite");
    mesh_.fill(scalar_,mesh_.cell_count(),value);mesh_.fill(old_scalar_,mesh_.cell_count(),value);mesh_.wait();time_=0.0;steps_=0U;
}

cfd::core::IterativeSolverResult ResidentScalarTransportSycl::step(
    const double* face_flux_device,const double* source_device) {
    if(!face_flux_device) throw std::invalid_argument("resident scalar transport requires face flux");
    const std::size_t bytes=mesh_.cell_count()*sizeof(double);mesh_.queue().memcpy(old_scalar_,scalar_,bytes);
    transfer_stats_.record_device_to_device(bytes);
    mesh_.assemble_scalar_transport_system(solver_->row_offsets_device(),solver_->column_indices_device(),solver_->nonzeros(),
                                           old_scalar_,face_flux_device,boundary_kind_,boundary_fixed_,config_.diffusivity,config_.dt,
                                           source_device,matrix_values_,rhs_);
    solver_->update_values_device(matrix_values_);
    const auto result=solver_->bicgstab_device(rhs_,scalar_,config_.linear_iterations,config_.linear_tolerance);
    if(result.converged){time_+=config_.dt;++steps_;}
    else {mesh_.queue().memcpy(scalar_,old_scalar_,bytes);transfer_stats_.record_device_to_device(bytes);}
    return result;
}

void ResidentScalarTransportSycl::download(std::span<double> scalar) const { mesh_.download_cell_scalar(scalar_,scalar); }
void ResidentScalarTransportSycl::reset_transfer_stats() const noexcept {
    transfer_stats_.reset();mesh_.reset_transfer_stats();if(solver_)solver_->reset_transfer_stats();
}
std::uint64_t ResidentScalarTransportSycl::hot_loop_host_transfer_bytes() const noexcept {
    return transfer_stats_.host_transfer_bytes()+mesh_.transfer_stats().host_transfer_bytes()+
           (solver_?solver_->transfer_stats().host_transfer_bytes():0U);
}
std::size_t ResidentScalarTransportSycl::resident_bytes() const noexcept {
    const std::size_t nc=mesh_.cell_count(),nf=mesh_.face_count(),nnz=solver_?solver_->nonzeros():0U;
    return mesh_.resident_bytes()+(3U*nc+nnz+nf)*sizeof(double)+nf*sizeof(std::uint8_t);
}

#endif

} // namespace cfd::fvm
