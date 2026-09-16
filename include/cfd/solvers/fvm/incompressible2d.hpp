#pragma once

#include "cfd/core/aligned_allocator.hpp"
#include "cfd/core/conjugate_gradient.hpp"

#include <cstddef>

namespace cfd::fvm {

struct Incompressible2DConfig {
    std::size_t nx{64};
    std::size_t ny{64};
    double lx{1.0};
    double ly{1.0};
    double density{1.0};
    double kinematic_viscosity{1.0e-2};
    double dt{2.5e-4};
    std::size_t pressure_iterations{800};
    double pressure_tolerance{1.0e-10};
};

// Periodic staggered finite-volume incompressible Navier-Stokes baseline.
// Momentum convection is conservative flux-form, diffusion is centered, and
// pressure/velocity coupling uses a projection with matrix-free CG.
class Incompressible2D {
public:
    explicit Incompressible2D(Incompressible2DConfig config);

    void initialize_taylor_green(double amplitude = 0.1);
    void step();
    void run(std::size_t steps);

    [[nodiscard]] double divergence_l2() const;
    [[nodiscard]] double kinetic_energy() const;
    [[nodiscard]] double taylor_green_velocity_error(double initial_amplitude = 0.1) const;
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] std::size_t steps() const noexcept { return steps_; }
    [[nodiscard]] const cfd::core::ConjugateGradientResult& pressure_result() const noexcept { return pressure_result_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& u_faces() const noexcept { return u_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& v_faces() const noexcept { return v_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& pressure() const noexcept { return p_; }
    [[nodiscard]] const Incompressible2DConfig& config() const noexcept { return config_; }

private:
    Incompressible2DConfig config_;
    double dx_{};
    double dy_{};
    double time_{};
    std::size_t steps_{};
    cfd::core::AlignedVector<double> u_;
    cfd::core::AlignedVector<double> v_;
    cfd::core::AlignedVector<double> u_star_;
    cfd::core::AlignedVector<double> v_star_;
    cfd::core::AlignedVector<double> p_;
    cfd::core::AlignedVector<double> pressure_rhs_;
    cfd::core::ConjugateGradientWorkspace pressure_workspace_;
    cfd::core::ConjugateGradientResult pressure_result_{};

    [[nodiscard]] std::size_t u_index(std::size_t i, std::size_t j) const noexcept;
    [[nodiscard]] std::size_t v_index(std::size_t i, std::size_t j) const noexcept;
    [[nodiscard]] std::size_t p_index(std::size_t i, std::size_t j) const noexcept;
    [[nodiscard]] double u_periodic(std::ptrdiff_t i, std::ptrdiff_t j) const noexcept;
    [[nodiscard]] double v_periodic(std::ptrdiff_t i, std::ptrdiff_t j) const noexcept;
    void enforce_periodic(cfd::core::AlignedVector<double>& u, cfd::core::AlignedVector<double>& v);
    void momentum_predictor();
    void project_star_velocity();
};

} // namespace cfd::fvm
