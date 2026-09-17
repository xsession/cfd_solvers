#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace cfd::fvm {

struct IdealGasEquationOfState {
    double gamma{1.4};
    double gas_constant{287.05};
    [[nodiscard]] double pressure(double density,double momentum,double total_energy) const;
    [[nodiscard]] double sound_speed(double density,double pressure_value) const;
    [[nodiscard]] double temperature(double density,double pressure_value) const;
};

struct EulerState1D {
    double density{};
    double momentum{};
    double energy{};
};

struct EulerPrimitive1D {
    double density{};
    double velocity{};
    double pressure{};
};

enum class CompressibleReconstruction {
    muscl_minmod,
    characteristic_weno5
};

enum class CompressibleTimeIntegrator {
    forward_euler,
    ssprk3
};

struct Compressible1DConfig {
    std::size_t cells{200};
    double length{1.0};
    double cfl{0.45};
    IdealGasEquationOfState eos{};
    bool periodic{false};
    CompressibleReconstruction reconstruction{CompressibleReconstruction::muscl_minmod};
    CompressibleTimeIntegrator time_integrator{CompressibleTimeIntegrator::forward_euler};
    double weno_epsilon{1.0e-6};
    double density_floor{1.0e-12};
    double pressure_floor{1.0e-12};
};

class CompressibleEuler1D {
public:
    explicit CompressibleEuler1D(Compressible1DConfig config={});

    void initialize(EulerPrimitive1D left,EulerPrimitive1D right,double interface_fraction=0.5);
    void initialize_uniform(EulerPrimitive1D state);
    void initialize_profile(std::span<const EulerPrimitive1D> states);

    [[nodiscard]] double stable_timestep() const;
    double step(double dt=0.0);
    void run(double duration);

    [[nodiscard]] const std::vector<EulerState1D>& conserved() const noexcept{return state_;}
    [[nodiscard]] EulerPrimitive1D primitive(std::size_t cell) const;
    [[nodiscard]] std::vector<double> pressure_jump_sensor() const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double momentum() const;
    [[nodiscard]] double total_energy() const;

private:
    [[nodiscard]] std::vector<EulerState1D> spatial_residual(std::span<const EulerState1D> state) const;
    [[nodiscard]] bool physical(const EulerState1D& state) const noexcept;
    void require_physical(std::span<const EulerState1D> state) const;

    Compressible1DConfig config_;
    std::vector<EulerState1D> state_,next_;
};

} // namespace cfd::fvm
