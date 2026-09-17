#pragma once

#include "cfd/core/iterative_solvers.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include "cfd/fvm/temporal.hpp"

#include <cstddef>
#include <functional>
#include <optional>
#include <string_view>
#include <vector>

namespace cfd::fvm {

enum class ScalarBoundaryType { zeroGradient, fixedValue };
struct ScalarBoundaryCondition { ScalarBoundaryType type{ScalarBoundaryType::zeroGradient}; double value{}; };
struct ScalarTransportConfig {
    double dt{1.0e-3}; double diffusivity{1.0e-2}; std::size_t linear_iterations{400};
    std::size_t gmres_restart{30}; double linear_tolerance{1.0e-10};
    TemporalScheme temporal_scheme{TemporalScheme::euler};
    double crank_nicolson_off_centering{1.0};
};

class ScalarTransport {
public:
    explicit ScalarTransport(PolyMesh mesh, ScalarTransportConfig config = {});
    void set_boundary(std::string_view patch, ScalarBoundaryType type, double value = 0.0);
    void initialize(double value);
    void initialize(const std::function<double(Vec3)>& value);
    void set_source(double value);
    void set_source(const std::function<double(Vec3)>& value);
    void set_face_flux(std::vector<double> volumetric_flux);
    cfd::core::IterativeSolverResult step();
    cfd::core::IterativeSolverResult run(std::size_t steps);
    [[nodiscard]] const PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }
    [[nodiscard]] const std::vector<double>& face_flux() const noexcept { return face_flux_; }
    [[nodiscard]] const cfd::core::IterativeSolverResult& linear_result() const noexcept { return linear_result_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double volume_integral() const;
    [[nodiscard]] double minimum() const;
    [[nodiscard]] double maximum() const;
    [[nodiscard]] std::size_t operator_assemblies() const noexcept { return operator_assemblies_; }
private:
    PolyMesh mesh_; ScalarTransportConfig config_; std::vector<ScalarBoundaryCondition> boundary_;
    std::vector<double> values_,old_values_,source_,face_flux_; cfd::core::CsrMatrix spatial_matrix_,matrix_;
    std::optional<cfd::core::Ilu0Preconditioner> preconditioner_;
    std::vector<double> boundary_rhs_,rhs_,candidate_,spatial_action_; bool operator_dirty_{true};
    std::size_t operator_assemblies_{}; std::size_t steps_{}; cfd::core::IterativeSolverResult linear_result_{}; double time_{};
    void assemble_operator();
};
} // namespace cfd::fvm
