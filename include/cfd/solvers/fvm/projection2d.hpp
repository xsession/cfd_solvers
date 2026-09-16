#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/conjugate_gradient.hpp"

#include <cstddef>

namespace cfd::fvm {

struct Projection2DConfig {
    std::size_t nx{64};
    std::size_t ny{64};
    double lx{1.0};
    double ly{1.0};
    double density{1.0};
    double dt{1.0e-3};
    std::size_t pressure_iterations{600};
    double pressure_tolerance{1.0e-10};
};

// Periodic staggered-grid pressure projection. This is the incompressibility
// kernel that later SIMPLE/PISO-style solvers can reuse.
class Projection2D {
public:
    explicit Projection2D(Projection2DConfig config);

    void initialize_divergent(double amplitude = 0.05);
    void initialize_taylor_green(double amplitude = 0.05);
    void project();

    [[nodiscard]] double divergence_l2() const;
    [[nodiscard]] double kinetic_energy() const;

    [[nodiscard]] const cfd::core::AlignedVector<double>& u_faces() const noexcept { return u_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& v_faces() const noexcept { return v_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& pressure() const noexcept { return p_; }
    [[nodiscard]] double dx() const noexcept { return dx_; }
    [[nodiscard]] double dy() const noexcept { return dy_; }
    [[nodiscard]] const Projection2DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const cfd::core::ConjugateGradientResult& pressure_result() const noexcept { return pressure_result_; }

private:
    Projection2DConfig config_;
    double dx_{};
    double dy_{};
    cfd::core::AlignedVector<double> u_;      // (nx+1) * ny vertical faces
    cfd::core::AlignedVector<double> v_;      // nx * (ny+1) horizontal faces
    cfd::core::AlignedVector<double> p_;      // nx * ny cell centers
    cfd::core::AlignedVector<double> rhs_;
    cfd::core::AlignedVector<double> pressure_rhs_;
    cfd::core::ConjugateGradientWorkspace pressure_workspace_;
    cfd::core::ConjugateGradientResult pressure_result_{};

    [[nodiscard]] std::size_t u_index(std::size_t i, std::size_t j) const noexcept;
    [[nodiscard]] std::size_t v_index(std::size_t i, std::size_t j) const noexcept;
    [[nodiscard]] std::size_t p_index(std::size_t i, std::size_t j) const noexcept;
    void enforce_periodic_faces();
};

} // namespace cfd::fvm
