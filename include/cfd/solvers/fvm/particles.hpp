#pragma once
#include "cfd/fvm/poly_mesh.hpp"
#include <span>
#include <vector>
namespace cfd::fvm {
struct Particle {Vec3 position{},velocity{};double diameter{1e-4};double density{1000.0};};
class ParticleCloud {
public:void add(Particle p);void advance(const PolyMesh& mesh,std::span<const Vec3> fluid_velocity,double dynamic_viscosity,double dt,bool two_way=false);[[nodiscard]] const std::vector<Particle>& particles()const noexcept{return particles_;}[[nodiscard]] const std::vector<Vec3>& fluid_momentum_source()const noexcept{return momentum_source_;}
private:std::vector<Particle>particles_;std::vector<Vec3>momentum_source_;
};

[[nodiscard]] Vec3 saffman_lift_force(Vec3 relative_velocity,Vec3 vorticity,double fluid_density,double kinematic_viscosity,double particle_diameter,double coefficient=1.615);
[[nodiscard]] Vec3 virtual_mass_force(Vec3 fluid_acceleration,Vec3 particle_acceleration,double fluid_density,double particle_diameter,double coefficient=0.5);
void collide_particles_elastic(Particle& a,Particle& b,double restitution=1.0);
[[nodiscard]] Particle coalesce_particles(const Particle& a,const Particle& b);
[[nodiscard]] bool breakup_by_weber(double fluid_density,double relative_speed,double particle_diameter,double surface_tension,double critical_weber=12.0);
[[nodiscard]] double evaporate_d2_law(double diameter,double evaporation_constant,double dt);
[[nodiscard]] std::vector<double> particle_volume_fraction(const PolyMesh& mesh,std::span<const Particle> particles);
}
