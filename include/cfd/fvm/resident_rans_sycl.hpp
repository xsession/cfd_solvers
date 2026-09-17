#pragma once

#include "cfd/fvm/resident_sycl.hpp"
#include "cfd/solvers/fvm/rans_transport.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

namespace cfd::fvm {

// Shared device-field owner layered on a single ResidentPolyMeshSycl. This is
// intentionally lightweight: coupled equations share one mesh/queue/context
// while the registry owns only named cell fields.
class ResidentFvmFieldRegistrySycl {
public:
    explicit ResidentFvmFieldRegistrySycl(ResidentPolyMeshSycl& mesh) noexcept : mesh_(mesh) {}
    ~ResidentFvmFieldRegistrySycl() noexcept;
    ResidentFvmFieldRegistrySycl(const ResidentFvmFieldRegistrySycl&) = delete;
    ResidentFvmFieldRegistrySycl& operator=(const ResidentFvmFieldRegistrySycl&) = delete;

    double* create_scalar(std::string name, double initial = 0.0);
    double* create_vector(std::string name, Vec3 initial = {});
    double* create_components(std::string name, std::size_t components, double initial = 0.0);
    [[nodiscard]] double* scalar(std::string_view name);
    [[nodiscard]] const double* scalar(std::string_view name) const;
    [[nodiscard]] double* vector(std::string_view name);
    [[nodiscard]] bool contains(std::string_view name) const;
    [[nodiscard]] std::size_t resident_bytes() const noexcept { return resident_bytes_; }
    [[nodiscard]] ResidentPolyMeshSycl& mesh() noexcept { return mesh_; }

private:
    struct Entry { double* data{}; std::size_t components{}; };
    ResidentPolyMeshSycl& mesh_;
    std::unordered_map<std::string, Entry> fields_;
    std::size_t resident_bytes_{};
    double* create(std::string name, std::size_t components);
};

// Generic shared-mesh resident scalar equation. Diffusivity, explicit source
// and linearized sink are device cell fields so source-coupled turbulence,
// thermal and species equations can be chained without host staging.
class ResidentScalarEquationSycl {
public:
    ResidentScalarEquationSycl(const PolyMesh& host_mesh,
                               ResidentFvmFieldRegistrySycl& registry,
                               std::string prefix,
                               double dt,
                               std::size_t linear_iterations = 300U,
                               double linear_tolerance = 1.0e-10);
    ~ResidentScalarEquationSycl() noexcept;
    ResidentScalarEquationSycl(const ResidentScalarEquationSycl&) = delete;
    ResidentScalarEquationSycl& operator=(const ResidentScalarEquationSycl&) = delete;

    void set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void initialize_uniform(double value);
    [[nodiscard]] cfd::core::IterativeSolverResult step(const double* face_flux_device,
                                                         double minimum = -1.0e300,
                                                         double maximum = 1.0e300);

    [[nodiscard]] double* field_device() noexcept { return field_; }
    [[nodiscard]] double* diffusivity_device() noexcept { return diffusivity_; }
    [[nodiscard]] double* source_device() noexcept { return source_; }
    [[nodiscard]] double* sink_device() noexcept { return sink_; }
    [[nodiscard]] const double* field_device() const noexcept { return field_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }
    void download(std::span<double> host) const;
    void reset_transfer_stats() const noexcept;
    [[nodiscard]] std::uint64_t hot_loop_host_transfer_bytes() const noexcept;

private:
    ResidentFvmFieldRegistrySycl& registry_;
    ResidentPolyMeshSycl& mesh_;
    std::string prefix_;
    double dt_{};
    std::size_t linear_iterations_{};
    double linear_tolerance_{};
    std::unique_ptr<cfd::core::SyclCsrLinearAlgebra> solver_;
    std::vector<std::string> patch_names_;
    std::vector<std::vector<std::size_t>> patch_faces_;
    std::vector<std::uint8_t> boundary_kind_host_;
    std::vector<double> boundary_fixed_host_;
    double* field_{};
    double* old_{};
    double* diffusivity_{};
    double* source_{};
    double* sink_{};
    double* rhs_{};
    double* matrix_values_{};
    std::uint8_t* boundary_kind_{};
    double* boundary_fixed_{};
    mutable cfd::core::DeviceTransferStats transfer_stats_{};
    double time_{};
    std::size_t steps_{};
    void upload_boundary_state();
};

class ResidentSpalartAllmarasSycl {
public:
    ResidentSpalartAllmarasSycl(const PolyMesh& host_mesh,
                                ResidentFvmFieldRegistrySycl& registry,
                                SpalartAllmarasConfig config = {},
                                std::string prefix = "sa");
    void initialize_uniform(double nu_tilde);
    void set_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_wall_distance(double distance);
    [[nodiscard]] cfd::core::IterativeSolverResult step(const double* face_flux_device,
                                                         const double* velocity_soa_device);
    [[nodiscard]] double* nu_tilde_device() noexcept { return equation_.field_device(); }
    [[nodiscard]] double* eddy_viscosity_device() noexcept { return nut_; }
    [[nodiscard]] double time() const noexcept { return equation_.time(); }
    void download(std::span<double> host) const { equation_.download(host); }
private:
    ResidentFvmFieldRegistrySycl& registry_;
    SpalartAllmarasConfig config_;
    ResidentScalarEquationSycl equation_;
    double* wall_distance_{}; double* strain_{}; double* velocity_gradient_{}; double* scalar_gradient_{}; double* nut_{};
};

class ResidentKEpsilonSycl {
public:
    ResidentKEpsilonSycl(const PolyMesh& host_mesh,
                         ResidentFvmFieldRegistrySycl& registry,
                         KEpsilonConfig config = {},
                         std::string prefix = "ke");
    void initialize_uniform(double k, double epsilon);
    void set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_epsilon_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    [[nodiscard]] TwoEquationLinearResult step(const double* face_flux_device,
                                                const double* velocity_soa_device);
    [[nodiscard]] double* k_device() noexcept { return k_.field_device(); }
    [[nodiscard]] double* epsilon_device() noexcept { return epsilon_.field_device(); }
    [[nodiscard]] double* eddy_viscosity_device() noexcept { return nut_; }
    [[nodiscard]] double time() const noexcept { return k_.time(); }
private:
    ResidentFvmFieldRegistrySycl& registry_;
    KEpsilonConfig config_;
    ResidentScalarEquationSycl k_, epsilon_;
    double* strain_{}; double* velocity_gradient_{}; double* nut_{};
};

class ResidentKOmegaSSTSycl {
public:
    ResidentKOmegaSSTSycl(const PolyMesh& host_mesh,
                          ResidentFvmFieldRegistrySycl& registry,
                          KOmegaSSTConfig config = {},
                          std::string prefix = "sst");
    void initialize_uniform(double k, double omega);
    void set_k_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_omega_boundary(std::string_view patch, TurbulenceScalarBoundaryType type, double value = 0.0);
    void set_wall_distance(double distance);
    void set_grid_scale(double scale);
    [[nodiscard]] TwoEquationLinearResult step(const double* face_flux_device,
                                                const double* velocity_soa_device);
    [[nodiscard]] double* k_device() noexcept { return k_.field_device(); }
    [[nodiscard]] double* omega_device() noexcept { return omega_.field_device(); }
    [[nodiscard]] double* eddy_viscosity_device() noexcept { return nut_; }
    [[nodiscard]] double* blending_f1_device() noexcept { return f1_; }
    [[nodiscard]] double* blending_f2_device() noexcept { return f2_; }
private:
    ResidentFvmFieldRegistrySycl& registry_;
    KOmegaSSTConfig config_;
    ResidentScalarEquationSycl k_, omega_;
    double* wall_distance_{}; double* grid_scale_{}; double* strain_{}; double* velocity_gradient_{};
    double* grad_k_{}; double* grad_omega_{}; double* nut_{}; double* f1_{}; double* f2_{}; double* hybrid_{};
};

// Device-side one-step reacting source coupling. Outputs are volumetric source
// fields that may be fed directly into shared ResidentScalarEquationSycl
// temperature/species equations.
class ResidentThermoSpeciesSourceSycl {
public:
    ResidentThermoSpeciesSourceSycl(ResidentFvmFieldRegistrySycl& registry,
                                    std::string prefix = "reacting");
    void arrhenius_one_step(const double* temperature,
                            const double* fuel_mass_fraction,
                            double pre_exponential,
                            double activation_temperature,
                            double heat_release);
    [[nodiscard]] double* heat_source_device() noexcept { return heat_source_; }
    [[nodiscard]] double* species_source_device() noexcept { return species_source_; }
    [[nodiscard]] double* reaction_rate_device() noexcept { return rate_; }
    void apply_sources(ResidentScalarEquationSycl& temperature_equation,
                       ResidentScalarEquationSycl& species_equation);
private:
    ResidentFvmFieldRegistrySycl& registry_;
    double* heat_source_{}; double* species_source_{}; double* rate_{};
};

} // namespace cfd::fvm
#endif
