#pragma once

#include "cfd/solvers/fdtd/monitor.hpp"

#include <complex>

namespace cfd::fdtd {

struct WaveAmplitudes1D {
    double forward{};
    double backward{};
};

// TEM-like 1-D port decomposition for the Maxwell1D E_z/H_y convention.
// For +x propagation H_y=-E_z/Z, while the reverse wave has H_y=+E_z/Z.
class WavePort1D {
public:
    WavePort1D(double frequency_hz,double impedance_ohm,int propagation_direction=1);

    void sample(double time,double electric_z,double magnetic_y);
    void clear() noexcept;

    [[nodiscard]] WaveAmplitudes1D instantaneous(double electric_z,double magnetic_y) const noexcept;
    [[nodiscard]] std::complex<double> forward_spectrum() const noexcept { return forward_.integral(); }
    [[nodiscard]] std::complex<double> backward_spectrum() const noexcept { return backward_.integral(); }
    [[nodiscard]] double frequency_hz() const noexcept { return frequency_hz_; }
    [[nodiscard]] double impedance_ohm() const noexcept { return impedance_ohm_; }

private:
    double frequency_hz_{};
    double impedance_ohm_{};
    int direction_{1};
    DftMonitor forward_;
    DftMonitor backward_;
};

struct SParameters1D {
    std::complex<double> s11{};
    std::complex<double> s21{};
};

[[nodiscard]] SParameters1D s_parameters(const WavePort1D& reference,const WavePort1D& transmitted);

} // namespace cfd::fdtd
