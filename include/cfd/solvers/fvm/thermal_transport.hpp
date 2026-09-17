#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <functional>
#include <span>
#include <string_view>
#include <vector>
namespace cfd::fvm {
struct ThermalBoundary { bool fixed_temperature{};double temperature{}; };
struct ThermalTransportConfig { double dt{1e-4};double density{1.0};double default_heat_capacity{1.0};double default_conductivity{1.0}; };
class ThermalTransport {
public:
    explicit ThermalTransport(PolyMesh mesh,ThermalTransportConfig config={});
    void initialize(double temperature);void initialize(const std::function<double(Vec3)>& temperature);
    void set_boundary(std::string_view patch,ThermalBoundary boundary);
    void set_volumetric_source(std::function<double(Vec3,double)> source){source_=std::move(source);}
    void set_material(std::function<double(Vec3,double)> heat_capacity,std::function<double(Vec3,double)> conductivity){cp_=std::move(heat_capacity);k_=std::move(conductivity);}
    void set_face_mass_flux(std::vector<double> mass_flux);
    void step();void run(std::size_t steps);
    [[nodiscard]] const std::vector<double>& temperature()const noexcept{return temperature_;} [[nodiscard]] double time()const noexcept{return time_;} [[nodiscard]] double total_sensible_enthalpy()const;
private:PolyMesh mesh_;ThermalTransportConfig config_;std::vector<double> temperature_,next_,mass_flux_;std::vector<ThermalBoundary> boundary_;std::function<double(Vec3,double)> source_,cp_,k_;double time_{};
};
[[nodiscard]] std::vector<Vec3> boussinesq_acceleration(const PolyMesh& mesh,std::span<const double> temperature,double reference_temperature,double beta,Vec3 gravity);
} // namespace cfd::fvm
