#pragma once

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/device_residency.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fvm/pressure_velocity.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>
#include "cfd/core/sycl_sparse.hpp"
#endif

namespace cfd::fvm {

// Symmetric orthogonal pressure operator. Dirichlet boundaries add diagonal
// contributions. With all-Neumann boundaries, set reference_cell to pin a
// zero-pressure reference while preserving an SPD matrix for CG.
[[nodiscard]] cfd::core::CsrMatrix build_orthogonal_pressure_matrix(
    const PolyMesh& mesh,
    double mobility = 1.0,
    bool fixed_value_boundaries = false,
    std::size_t reference_cell = invalid_cell);

[[nodiscard]] cfd::core::CsrMatrix build_orthogonal_pressure_matrix(
    const PolyMesh& mesh,
    std::span<const double> cell_mobility,
    bool fixed_value_boundaries = false,
    std::size_t reference_cell = invalid_cell);

#if defined(CFD_HAS_SYCL)


struct ResidentPolyMeshDeviceView {
    std::size_t cell_count{};
    std::size_t face_count{};
    const std::size_t* face_owner{};
    const std::size_t* face_neighbour{};
    const std::size_t* cell_face_offsets{};
    const std::size_t* cell_face_indices{};
    const double* cell_geometry{};
    const double* face_geometry{};
};

// Device-resident mirror of the geometry/connectivity needed by the core FVM
// operators. Vector fields use SoA layout: [x(0..N), y(0..N), z(0..N)].
class ResidentPolyMeshSycl {
public:
    explicit ResidentPolyMeshSycl(
        const PolyMesh& mesh,
        sycl::device device = sycl::device{sycl::default_selector_v});
    ResidentPolyMeshSycl(const PolyMesh& mesh, sycl::queue queue);
    ~ResidentPolyMeshSycl() noexcept;

    ResidentPolyMeshSycl(const ResidentPolyMeshSycl&) = delete;
    ResidentPolyMeshSycl& operator=(const ResidentPolyMeshSycl&) = delete;
    ResidentPolyMeshSycl(ResidentPolyMeshSycl&&) = delete;
    ResidentPolyMeshSycl& operator=(ResidentPolyMeshSycl&&) = delete;

    [[nodiscard]] std::size_t cell_count() const noexcept { return cell_count_; }
    [[nodiscard]] std::size_t face_count() const noexcept { return face_count_; }
    [[nodiscard]] std::size_t adjacency_count() const noexcept { return adjacency_count_; }
    [[nodiscard]] sycl::queue& queue() noexcept { return queue_; }
    [[nodiscard]] const sycl::queue& queue() const noexcept { return queue_; }
    [[nodiscard]] ResidentPolyMeshDeviceView device_view() const noexcept;

    [[nodiscard]] std::size_t resident_bytes() const noexcept;
    [[nodiscard]] const cfd::core::DeviceTransferStats& transfer_stats() const noexcept {
        return transfer_stats_;
    }
    void reset_transfer_stats() const noexcept { transfer_stats_.reset(); }

    // Explicit setup/output boundary helpers. The hot-loop operators below do
    // not perform host transfers themselves.
    void upload_cell_scalar(std::span<const double> host, double* device) const;
    void download_cell_scalar(const double* device, std::span<double> host) const;
    void upload_face_scalar(std::span<const double> host, double* device) const;
    void download_face_scalar(const double* device, std::span<double> host) const;
    void upload_cell_vector(std::span<const Vec3> host, double* device_soa) const;
    void download_cell_vector(const double* device_soa, std::span<Vec3> host) const;
    void fill(double* device, std::size_t count, double value) const;
    void set_cell_value(double* device, std::size_t cell, double value) const;
    void wait() const;

    void interpolate_scalar_to_faces(const double* cell_values,
                                     double* face_values,
                                     const double* boundary_face_values = nullptr) const;
    void gauss_gradient_scalar(const double* cell_values,
                               double* cell_gradient_soa,
                               const double* boundary_face_values = nullptr) const;
    void gauss_divergence_vector(const double* cell_vector_soa,
                                 double* cell_divergence,
                                 const double* boundary_face_vector_soa = nullptr) const;
    void orthogonal_laplacian_scalar(const double* cell_values,
                                     double diffusivity,
                                     double* cell_laplacian,
                                     const double* boundary_face_values = nullptr) const;

    // Pressure/velocity coupling primitives. predictor_face_flux() forms H/A
    // dot Sf. pressure_rhs_from_flux() forms the unnormalised RHS compatible
    // with build_orthogonal_pressure_matrix(). A non-null boundary pressure
    // array means fixed-value pressure on every boundary face; nullptr means
    // zero-gradient pressure there.
    void predictor_face_flux(const double* h_by_a_soa,
                             double* face_flux,
                             const double* boundary_velocity_soa = nullptr) const;
    void pressure_rhs_from_flux(const double* predictor_face_flux,
                                const double* cell_mobility,
                                double* rhs,
                                const double* boundary_pressure = nullptr,
                                const std::uint8_t* boundary_pressure_kind = nullptr,
                                const double* nonorthogonal_flux = nullptr) const;
    void pressure_correct_face_flux(const double* pressure,
                                    const double* cell_mobility,
                                    double* face_flux_in_out,
                                    const double* boundary_pressure = nullptr,
                                    const std::uint8_t* boundary_pressure_kind = nullptr,
                                    const double* nonorthogonal_flux = nullptr) const;
    void correct_velocity_from_pressure_gradient(const double* h_by_a_soa,
                                                 const double* pressure_mobility,
                                                 const double* pressure_gradient_soa,
                                                 double* velocity_soa) const;
    void divergence_face_flux(const double* face_flux,
                              double* cell_divergence) const;

    // Mixed velocity-boundary handling used by the resident incompressible
    // loop. boundary_kind uses VelocityBoundaryType numeric values per face;
    // fixed/slip values use face-vector SoA layout.
    void boundary_velocity_values(const double* cell_velocity_soa,
                                  const std::uint8_t* boundary_kind,
                                  const double* boundary_fixed_soa,
                                  double* boundary_velocity_soa) const;
    void boundary_pressure_values(const double* cell_pressure,
                                  const std::uint8_t* boundary_kind,
                                  const double* boundary_fixed,
                                  double* boundary_pressure) const;

    // Assemble the first-order Euler/upwind momentum diagonal and RHS entirely
    // on-device. The off-diagonal action is evaluated matrix-free by
    // momentum_jacobi_sweep(), avoiding a host rebuild for nonlinear fluxes.
    void assemble_momentum_system(const double* time_source_soa,
                                  const double* pressure_gradient_soa,
                                  const double* face_flux,
                                  const std::uint8_t* boundary_kind,
                                  const double* boundary_velocity_soa,
                                  double kinematic_viscosity,
                                  double dt,
                                  bool include_convection,
                                  double* diagonal,
                                  double* mobility,
                                  double* rhs_soa,
                                  const double* acceleration_soa = nullptr) const;
    void momentum_jacobi_sweep(const double* current_velocity_soa,
                               const double* face_flux,
                               const std::uint8_t* boundary_kind,
                               double kinematic_viscosity,
                               bool include_convection,
                               const double* diagonal,
                               const double* rhs_soa,
                               double* next_velocity_soa) const;
    void assemble_momentum_matrix_values(const std::size_t* row_offsets,
                                         const std::size_t* columns,
                                         std::size_t nonzeros,
                                         const double* face_flux,
                                         const std::uint8_t* boundary_kind,
                                         double kinematic_viscosity,
                                         double dt,
                                         bool include_convection,
                                         double* values) const;
    void assemble_scalar_transport_system(const std::size_t* row_offsets,
                                          const std::size_t* columns,
                                          std::size_t nonzeros,
                                          const double* old_scalar,
                                          const double* face_flux,
                                          const std::uint8_t* boundary_kind,
                                          const double* boundary_fixed,
                                          double diffusivity,
                                          double dt,
                                          const double* volumetric_source,
                                          double* values,
                                          double* rhs) const;
    // Variable-diffusivity scalar transport with a linearized non-negative
    // sink coefficient. This is the reusable resident equation kernel for
    // turbulence, thermal and reacting-scalar transport.
    void assemble_scalar_transport_system_variable(
        const std::size_t* row_offsets,
        const std::size_t* columns,
        std::size_t nonzeros,
        const double* old_scalar,
        const double* face_flux,
        const std::uint8_t* boundary_kind,
        const double* boundary_fixed,
        const double* cell_diffusivity,
        double dt,
        const double* volumetric_source,
        const double* implicit_sink,
        double* values,
        double* rhs) const;
    void clamp_scalar(double* field, double minimum, double maximum) const;
    void velocity_strain_rate_magnitude(const double* velocity_soa,
                                        double* gradient_scratch_9soa,
                                        double* strain_rate) const;
    void form_h_by_a(const double* velocity_soa,
                     const double* mobility,
                     const double* pressure_gradient_soa,
                     double* h_by_a_soa) const;

    // Variable-mobility pressure support. Numeric CSR values are regenerated
    // on-device against a fixed pinned-Neumann sparsity pattern.
    void assemble_pinned_pressure_values(const std::size_t* row_offsets,
                                         const std::size_t* columns,
                                         std::size_t nonzeros,
                                         const double* cell_mobility,
                                         std::size_t reference_cell,
                                         double* values,
                                         const std::uint8_t* boundary_pressure_kind = nullptr,
                                         bool pin_reference = true) const;
    void nonorthogonal_pressure_flux(const double* pressure,
                                     const double* cell_mobility,
                                     double* pressure_gradient_soa,
                                     double* face_correction,
                                     const double* boundary_pressure = nullptr,
                                     const std::uint8_t* boundary_pressure_kind = nullptr) const;
    void volume_weighted_vector_norm2(const double* vector_soa,
                                      double* shared_result) const;

private:
    mutable sycl::queue queue_;
    std::size_t cell_count_{};
    std::size_t face_count_{};
    std::size_t adjacency_count_{};

    std::size_t* face_owner_{nullptr};
    std::size_t* face_neighbour_{nullptr};
    std::size_t* cell_face_offsets_{nullptr};
    std::size_t* cell_face_indices_{nullptr};
    double* cell_geometry_{nullptr};   // 4*N: cx, cy, cz, volume
    double* face_geometry_{nullptr};   // 9*F: fc xyz, Sf xyz, ownerWeight, normalMetric, orthMetric
    mutable cfd::core::DeviceTransferStats transfer_stats_{};

    void initialize(const PolyMesh& mesh);
    void release() noexcept;
};

// Constant-mobility pressure-projection baseline that keeps all volume fields
// resident. A future resident momentum solver can write H/A directly into
// predictor_device(), call project(), and continue with velocity/face flux
// without a bulk host transfer.
class ResidentPressureProjectionSycl {
public:
    ResidentPressureProjectionSycl(
        const PolyMesh& mesh,
        double mobility = 1.0,
        std::size_t reference_cell = 0U,
        sycl::device device = sycl::device{sycl::default_selector_v});
    ResidentPressureProjectionSycl(const PolyMesh& mesh,
                                   double mobility,
                                   std::size_t reference_cell,
                                   sycl::queue queue);
    ~ResidentPressureProjectionSycl() noexcept;

    ResidentPressureProjectionSycl(const ResidentPressureProjectionSycl&) = delete;
    ResidentPressureProjectionSycl& operator=(const ResidentPressureProjectionSycl&) = delete;
    ResidentPressureProjectionSycl(ResidentPressureProjectionSycl&&) = delete;
    ResidentPressureProjectionSycl& operator=(ResidentPressureProjectionSycl&&) = delete;

    [[nodiscard]] ResidentPolyMeshSycl& mesh() noexcept { return mesh_; }
    [[nodiscard]] const ResidentPolyMeshSycl& mesh() const noexcept { return mesh_; }
    [[nodiscard]] double* predictor_device() noexcept { return h_by_a_; }
    [[nodiscard]] double* pressure_device() noexcept { return pressure_; }
    [[nodiscard]] double* velocity_device() noexcept { return velocity_; }
    [[nodiscard]] double* face_flux_device() noexcept { return face_flux_; }
    [[nodiscard]] double* mobility_device() noexcept { return mobility_; }
    [[nodiscard]] std::size_t reference_cell() const noexcept { return reference_cell_; }

    void set_predictor(std::span<const Vec3> velocity);
    void set_pressure(std::span<const double> pressure);
    [[nodiscard]] cfd::core::IterativeSolverResult project(
        std::size_t max_iterations = 400U,
        double relative_tolerance = 1.0e-10);
    [[nodiscard]] double continuity_l2();
    void download_velocity(std::span<Vec3> velocity) const;
    void download_pressure(std::span<double> pressure) const;

    void reset_transfer_stats() const noexcept;
    [[nodiscard]] std::uint64_t hot_loop_host_transfer_bytes() const noexcept;
    [[nodiscard]] std::size_t resident_bytes() const noexcept;

private:
    ResidentPolyMeshSycl mesh_;
    std::unique_ptr<cfd::core::SyclCsrLinearAlgebra> pressure_solver_;
    std::size_t reference_cell_{};
    double* h_by_a_{nullptr};
    double* pressure_{nullptr};
    double* pressure_gradient_{nullptr};
    double* velocity_{nullptr};
    double* face_flux_{nullptr};
    double* rhs_{nullptr};
    double* mobility_{nullptr};
    double* boundary_velocity_{nullptr};
    double* continuity_{nullptr};
    double* reduction_scalar_{nullptr};

    void allocate(double mobility);
    void release() noexcept;
};

struct ResidentIncompressibleConfig {
    double density{1.0};
    double kinematic_viscosity{1.0e-2};
    double dt{1.0e-2};
    std::size_t momentum_sweeps{40U}; // legacy Jacobi fallback/debug control
    std::size_t momentum_iterations{200U};
    double momentum_tolerance{1.0e-9};
    std::size_t pressure_iterations{600U};
    double pressure_tolerance{1.0e-10};
    std::size_t nonorthogonal_correctors{2U};
    std::size_t pressure_correctors{2U};
    std::size_t outer_correctors{2U};
    double velocity_relaxation{0.7};
    double pressure_relaxation{0.3};
    bool include_convection{true};
};

struct ResidentIncompressibleIterationInfo {
    double velocity_rms_change{};
    double continuity_l2{};
    cfd::core::IterativeSolverResult pressure{};
    std::size_t pressure_solves{};
};

// GPU-resident first-order collocated incompressible baseline. The class keeps
// all hot fields and nonlinear coefficient refreshes on the same in-order SYCL
// queue. Momentum uses matrix-free damped/Jacobi sweeps for the nonsymmetric
// upwind operator; pressure uses the persistent device CSR CG solver.
class ResidentIncompressibleSycl {
public:
    ResidentIncompressibleSycl(
        const PolyMesh& mesh,
        ResidentIncompressibleConfig config = {},
        std::size_t reference_cell = 0U,
        sycl::device device = sycl::device{sycl::default_selector_v});
    ResidentIncompressibleSycl(const PolyMesh& mesh,
                               ResidentIncompressibleConfig config,
                               std::size_t reference_cell,
                               sycl::queue queue);
    ~ResidentIncompressibleSycl() noexcept;

    ResidentIncompressibleSycl(const ResidentIncompressibleSycl&) = delete;
    ResidentIncompressibleSycl& operator=(const ResidentIncompressibleSycl&) = delete;
    ResidentIncompressibleSycl(ResidentIncompressibleSycl&&) = delete;
    ResidentIncompressibleSycl& operator=(ResidentIncompressibleSycl&&) = delete;

    void set_velocity_boundary(std::string_view patch,
                               VelocityBoundaryType type,
                               Vec3 value = {});
    void set_pressure_boundary(std::string_view patch,
                               PressureBoundaryType type,
                               double value = 0.0);
    void initialize_uniform(Vec3 velocity = {}, double pressure = 0.0);

    [[nodiscard]] ResidentIncompressibleIterationInfo iterate_simple();
    [[nodiscard]] ResidentIncompressibleIterationInfo step_piso();
    [[nodiscard]] ResidentIncompressibleIterationInfo step_pimple();

    [[nodiscard]] ResidentPolyMeshSycl& mesh() noexcept { return mesh_; }
    [[nodiscard]] const ResidentPolyMeshSycl& mesh() const noexcept { return mesh_; }
    [[nodiscard]] double* velocity_device() noexcept { return velocity_; }
    [[nodiscard]] double* pressure_device() noexcept { return pressure_; }
    [[nodiscard]] double* face_flux_device() noexcept { return face_flux_; }
    [[nodiscard]] double* mobility_device() noexcept { return mobility_; }
    [[nodiscard]] double* external_acceleration_device() noexcept { return external_acceleration_; }
    void clear_external_acceleration();
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }

    [[nodiscard]] double continuity_l2();
    [[nodiscard]] double kinetic_energy();
    void download_velocity(std::span<Vec3> velocity) const;
    void download_pressure(std::span<double> pressure) const;

    void reset_transfer_stats() const noexcept;
    [[nodiscard]] std::uint64_t hot_loop_host_transfer_bytes() const noexcept;
    [[nodiscard]] std::size_t resident_bytes() const noexcept;

private:
    ResidentPolyMeshSycl mesh_;
    ResidentIncompressibleConfig config_;
    std::unique_ptr<cfd::core::SyclCsrLinearAlgebra> pressure_solver_;
    std::unique_ptr<cfd::core::SyclCsrLinearAlgebra> momentum_solver_;
    std::size_t reference_cell_{};
    std::vector<std::string> patch_names_;
    std::vector<std::vector<std::size_t>> patch_faces_;
    std::vector<std::uint8_t> boundary_kind_host_;
    std::vector<double> boundary_fixed_host_;
    std::vector<std::uint8_t> pressure_boundary_kind_host_;
    std::vector<double> pressure_boundary_fixed_host_;
    bool has_fixed_pressure_boundary_{false};

    double* velocity_{nullptr};
    double* old_velocity_{nullptr};
    double* source_velocity_{nullptr};
    double* velocity_before_{nullptr};
    double* velocity_work_{nullptr};
    double* h_by_a_{nullptr};
    double* pressure_{nullptr};
    double* pressure_before_{nullptr};
    double* pressure_gradient_{nullptr};
    double* face_flux_{nullptr};
    double* predictor_flux_{nullptr};
    double* nonorthogonal_flux_{nullptr};
    double* rhs_{nullptr};
    double* mobility_{nullptr};
    double* momentum_diagonal_{nullptr};
    double* momentum_rhs_{nullptr};
    double* external_acceleration_{nullptr};
    double* boundary_velocity_{nullptr};
    std::uint8_t* boundary_kind_{nullptr};
    double* boundary_fixed_{nullptr};
    std::uint8_t* pressure_boundary_kind_{nullptr};
    double* pressure_boundary_fixed_{nullptr};
    double* pressure_boundary_values_{nullptr};
    double* pressure_matrix_values_{nullptr};
    double* momentum_matrix_values_{nullptr};
    double* continuity_{nullptr};
    double* reduction_scalar_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};
    double time_{};
    std::size_t steps_{};

    void allocate();
    void release() noexcept;
    void upload_boundary_state();
    void rebuild_face_flux();
    void momentum_predictor(const double* time_source);
    [[nodiscard]] cfd::core::IterativeSolverResult correct_pressure(double pressure_relaxation);
    [[nodiscard]] ResidentIncompressibleIterationInfo coupled_sequence(
        const double* time_source,
        std::size_t pressure_correctors,
        double pressure_relaxation);
    [[nodiscard]] double velocity_rms_change();
};

struct ResidentScalarTransportConfig {
    double dt{1.0e-2};
    double diffusivity{1.0e-3};
    std::size_t linear_iterations{300U};
    double linear_tolerance{1.0e-10};
};

// Device-resident first-order implicit scalar advection-diffusion transport.
// A caller may pass ResidentIncompressibleSycl::face_flux_device() directly;
// if both solvers use the same SYCL context, no field staging is required.
class ResidentScalarTransportSycl {
public:
    ResidentScalarTransportSycl(
        const PolyMesh& mesh,
        ResidentScalarTransportConfig config = {},
        sycl::device device = sycl::device{sycl::default_selector_v});
    ResidentScalarTransportSycl(const PolyMesh& mesh,
                                ResidentScalarTransportConfig config,
                                sycl::queue queue);
    ~ResidentScalarTransportSycl() noexcept;

    ResidentScalarTransportSycl(const ResidentScalarTransportSycl&) = delete;
    ResidentScalarTransportSycl& operator=(const ResidentScalarTransportSycl&) = delete;
    ResidentScalarTransportSycl(ResidentScalarTransportSycl&&) = delete;
    ResidentScalarTransportSycl& operator=(ResidentScalarTransportSycl&&) = delete;

    void set_boundary(std::string_view patch, PressureBoundaryType type, double value = 0.0);
    void initialize_uniform(double value = 0.0);
    [[nodiscard]] cfd::core::IterativeSolverResult step(
        const double* face_flux_device,
        const double* volumetric_source_device = nullptr);

    [[nodiscard]] double* scalar_device() noexcept { return scalar_; }
    [[nodiscard]] const double* scalar_device() const noexcept { return scalar_; }
    [[nodiscard]] ResidentPolyMeshSycl& mesh() noexcept { return mesh_; }
    [[nodiscard]] const ResidentPolyMeshSycl& mesh() const noexcept { return mesh_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }
    void download(std::span<double> scalar) const;
    void reset_transfer_stats() const noexcept;
    [[nodiscard]] std::uint64_t hot_loop_host_transfer_bytes() const noexcept;
    [[nodiscard]] std::size_t resident_bytes() const noexcept;

private:
    ResidentPolyMeshSycl mesh_;
    ResidentScalarTransportConfig config_;
    std::unique_ptr<cfd::core::SyclCsrLinearAlgebra> solver_;
    std::vector<std::string> patch_names_;
    std::vector<std::vector<std::size_t>> patch_faces_;
    std::vector<std::uint8_t> boundary_kind_host_;
    std::vector<double> boundary_fixed_host_;
    double* scalar_{nullptr};
    double* old_scalar_{nullptr};
    double* rhs_{nullptr};
    double* matrix_values_{nullptr};
    std::uint8_t* boundary_kind_{nullptr};
    double* boundary_fixed_{nullptr};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};
    double time_{};
    std::size_t steps_{};

    void allocate();
    void release() noexcept;
    void upload_boundary_state();
};

#endif

} // namespace cfd::fvm
