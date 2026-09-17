#pragma once
#include <cstddef>
#include <vector>
namespace cfd::fvm {
struct IdealGasEquationOfState {
    double gamma{1.4};
    double gas_constant{287.05};
    [[nodiscard]] double pressure(double density,double momentum,double total_energy) const;
    [[nodiscard]] double sound_speed(double density,double pressure_value) const;
    [[nodiscard]] double temperature(double density,double pressure_value) const;
};
struct EulerState1D {double density{},momentum{},energy{};};
struct EulerPrimitive1D {double density{},velocity{},pressure{};};
struct Compressible1DConfig {std::size_t cells{200};double length{1.0};double cfl{0.45};IdealGasEquationOfState eos{};bool periodic{false};};
class CompressibleEuler1D {
public:
    explicit CompressibleEuler1D(Compressible1DConfig config={});
    void initialize(EulerPrimitive1D left,EulerPrimitive1D right,double interface_fraction=0.5);
    void initialize_uniform(EulerPrimitive1D state);
    [[nodiscard]] double stable_timestep() const;
    double step(double dt=0.0);
    void run(double duration);
    [[nodiscard]] const std::vector<EulerState1D>& conserved() const noexcept{return state_;}
    [[nodiscard]] EulerPrimitive1D primitive(std::size_t cell) const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double total_energy() const;
private:
    Compressible1DConfig config_;
    std::vector<EulerState1D> state_,next_;
};
} // namespace cfd::fvm
