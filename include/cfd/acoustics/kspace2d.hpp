#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace cfd::acoustics {

struct KSpace2DConfig {
    std::size_t nx{128};
    std::size_t ny{128};
    double dx{1.0e-4};
    double dy{1.0e-4};
    double dt{0.0};
    double cfl{0.2};
    double reference_sound_speed{0.0};
    std::size_t pml_cells{12};
    double pml_order{3.0};
    double pml_target_reflection{1.0e-8};
};

struct PowerLawAbsorption2D {
    double alpha0_np_per_m{0.0};
    double exponent{1.5};
    double reference_frequency_hz{1.0e6};
};

struct AcousticSensor2D {
    std::size_t cell{};
    std::vector<double> pressure;
};

class KSpaceAcoustic2D {
public:
    explicit KSpaceAcoustic2D(KSpace2DConfig config);

    void set_uniform_medium(double sound_speed_m_s,double density_kg_m3);
    void set_medium(std::span<const double> sound_speed_m_s,std::span<const double> density_kg_m3);
    void set_nonlinearity(std::span<const double> bon_a);
    void set_uniform_nonlinearity(double bon_a);
    void set_power_law_absorption(PowerLawAbsorption2D absorption);
    void clear_absorption() noexcept;

    void initialize_pressure(std::span<const double> pressure_pa);
    void initialize_gaussian_pressure(double center_x_fraction,double center_y_fraction,double sigma_m,double amplitude_pa);
    void add_pressure_source(std::size_t cell,std::function<double(double)> value_pa);
    void add_plane_source_x(std::size_t i,std::function<double(double)> value_pa);
    void add_delayed_source_array(std::span<const std::size_t> cells,std::span<const double> delays_s,std::function<double(double)> waveform_pa);
    std::size_t add_sensor(std::size_t cell);
    std::vector<std::size_t> add_sensor_array(std::span<const std::size_t> cells);
    void clear_sensors();

    // Reconstruct an initial pressure field by enforcing recorded sensor
    // pressure in reverse order as a Dirichlet condition on the sensor mask.
    // Existing sources are ignored during this reconstruction.
    [[nodiscard]] std::vector<double> time_reversal_reconstruct(
        std::span<const std::size_t> sensor_cells,
        std::span<const std::vector<double>> sensor_pressure,
        bool positivity = false);

    void step(std::size_t count=1);
    void run(std::size_t steps){ step(steps); }

    [[nodiscard]] const std::vector<double>& pressure() const noexcept { return pressure_; }
    [[nodiscard]] const std::vector<double>& velocity_x() const noexcept { return ux_; }
    [[nodiscard]] const std::vector<double>& velocity_y() const noexcept { return uy_; }
    [[nodiscard]] const std::vector<double>& density_perturbation() const noexcept { return rho_; }
    [[nodiscard]] const std::vector<AcousticSensor2D>& sensors() const noexcept { return sensors_; }
    [[nodiscard]] double dt() const noexcept { return dt_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double total_energy() const;
    [[nodiscard]] std::size_t cells() const noexcept { return config_.nx*config_.ny; }
    [[nodiscard]] const KSpace2DConfig& config() const noexcept { return config_; }
    [[nodiscard]] const std::vector<double>& sound_speed() const noexcept { return c0_; }
    [[nodiscard]] const std::vector<double>& medium_density() const noexcept { return rho0_; }

private:
    struct Source { std::size_t cell{}; std::function<double(double)> value; };
    KSpace2DConfig config_;
    double dt_{};
    double time_{};
    std::vector<double> c0_,rho0_,bon_a_;
    std::vector<double> ux_,uy_,rho_x_,rho_y_,rho_,pressure_;
    std::vector<double> damp_x_,damp_y_;
    std::vector<double> kx_,ky_,kappa_;
    std::vector<Source> sources_;
    std::vector<AcousticSensor2D> sensors_;
    PowerLawAbsorption2D absorption_{};
    bool absorption_enabled_{false};

    void rebuild_time_step_and_kappa();
    void rebuild_pml();
    void update_pressure_from_density();
    void apply_absorption();
    void spectral_gradient(std::span<const double> input,std::vector<double>& gx,std::vector<double>& gy) const;
    void record_sensors();
};

[[nodiscard]] std::vector<double> fractional_laplacian_periodic_2d(std::span<const double> field,
                                                                    std::size_t nx,std::size_t ny,
                                                                    double dx,double dy,double order);

[[nodiscard]] std::vector<double> power_law_attenuate_periodic_2d(std::span<const double> field,
                                                                   std::size_t nx,std::size_t ny,
                                                                   double dx,double dy,
                                                                   double c_ref,double dt,
                                                                   PowerLawAbsorption2D absorption);

} // namespace cfd::acoustics
