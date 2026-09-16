#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace cfd::fdtd {

class TimeProbe {
public:
    void sample(double time,double value);
    void clear() noexcept;
    [[nodiscard]] const std::vector<double>& time() const noexcept{return time_;}
    [[nodiscard]] const std::vector<double>& value() const noexcept{return value_;}
private:
    std::vector<double> time_;
    std::vector<double> value_;
};

class DftMonitor {
public:
    explicit DftMonitor(double frequency_hz);
    void sample(double time,double value);
    void clear() noexcept;
    [[nodiscard]] std::complex<double> integral() const noexcept{return integral_;}
    // Two-sided real-sinusoid amplitude estimate using rectangular integration.
    [[nodiscard]] double amplitude(double duration) const;
    [[nodiscard]] std::size_t samples() const noexcept{return samples_;}
private:
    double frequency_{};
    std::complex<double> integral_{};
    double previous_time_{};
    double previous_value_{};
    bool have_previous_{};
    std::size_t samples_{};
};

} // namespace cfd::fdtd
