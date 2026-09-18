#include "cfd/multibody/cfd_dem_coupling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace cfd::multibody {
namespace {

[[nodiscard]] Vec3 to_mb(cfd::fvm::Vec3 v) noexcept { return {v.x,v.y,v.z}; }

[[nodiscard]] Vec3 cross3(Vec3 a, Vec3 b) noexcept {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}
[[nodiscard]] double mag(Vec3 a) noexcept { return std::sqrt(dot(a,a)); }

[[nodiscard]] std::size_t nearest_cell(const cfd::fvm::PolyMesh& mesh, Vec3 p) {
    if(mesh.cell_count()==0U) throw std::invalid_argument("CFD/DEM coupling requires non-empty mesh");
    std::size_t best=0U; double d2=std::numeric_limits<double>::infinity();
    for(std::size_t c=0;c<mesh.cell_count();++c){
        const auto q=mesh.cells()[c].center;
        const double dx=q.x-p.x,dy=q.y-p.y,dz=q.z-p.z;
        const double r2=dx*dx+dy*dy+dz*dz;
        if(r2<d2){d2=r2;best=c;}
    }
    return best;
}

[[nodiscard]] double sphere_volume(double r) noexcept {
    return (4.0/3.0)*std::numbers::pi*r*r*r;
}

[[nodiscard]] Vec3 schiller_naumann_drag(Vec3 rel,double diameter,double rho,double mu,double void_fraction,double exponent) {
    const double speed=mag(rel);
    if(speed==0.0) return {};
    const double re=rho*speed*diameter/mu;
    const double cd=re<1000.0 ? (24.0/std::max(re,1.0e-30))*(1.0+0.15*std::pow(re,0.687)) : 0.44;
    const double area=0.25*std::numbers::pi*diameter*diameter;
    const double hindrance=std::pow(std::max(void_fraction,1.0e-12),-exponent);
    return rel*(0.5*cd*rho*area*speed*hindrance);
}

[[nodiscard]] Vec3 saffman(Vec3 rel,Vec3 vort,double rho,double nu,double d,double coefficient) {
    const double om=mag(vort);
    if(om==0.0 || mag(rel)==0.0) return {};
    const double scale=coefficient*rho*d*d*std::sqrt(nu*om)/om;
    return cross3(rel,vort)*scale;
}

} // namespace

Vec3 SymmetricStress3::traction(Vec3 n) const noexcept {
    return {xx*n.x+xy*n.y+xz*n.z,
            xy*n.x+yy*n.y+yz*n.z,
            xz*n.x+yz*n.y+zz*n.z};
}

std::vector<RigidSurfaceSample> sphere_surface_quadrature(const RigidBody& body,double radius,std::size_t count) {
    if(!(radius>0.0)||!std::isfinite(radius)) throw std::invalid_argument("sphere quadrature radius must be positive");
    if(count<6U) throw std::invalid_argument("sphere quadrature requires at least six samples");
    std::vector<Vec3> normals;
    normals.reserve(count);
    if(count==6U){
        normals={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    } else {
        // Antipodal Fibonacci pairs preserve zero net area vector exactly.
        if(count%2U!=0U) ++count;
        const std::size_t half=count/2U;
        const double golden=std::numbers::pi*(3.0-std::sqrt(5.0));
        for(std::size_t i=0;i<half;++i){
            const double z=(static_cast<double>(i)+0.5)/static_cast<double>(half);
            const double rr=std::sqrt(std::max(0.0,1.0-z*z));
            const double phi=golden*static_cast<double>(i);
            Vec3 n{rr*std::cos(phi),rr*std::sin(phi),z};
            normals.push_back(n); normals.push_back(n*(-1.0));
        }
    }
    const double area=4.0*std::numbers::pi*radius*radius/static_cast<double>(normals.size());
    std::vector<RigidSurfaceSample> out;out.reserve(normals.size());
    for(const auto n: normals) out.push_back({body.state.position+n*radius,n,area});
    return out;
}

ResolvedSurfaceLoad integrate_resolved_surface_load(
    const cfd::fvm::PolyMesh& mesh,std::span<const double> pressure,
    std::span<const SymmetricStress3> stress,std::span<const RigidSurfaceSample> surface,Vec3 center) {
    if(pressure.size()!=mesh.cell_count()) throw std::invalid_argument("resolved CFD/body pressure size mismatch");
    if(!stress.empty()&&stress.size()!=mesh.cell_count()) throw std::invalid_argument("resolved CFD/body stress size mismatch");
    ResolvedSurfaceLoad result{};
    for(const auto& s:surface){
        if(!(s.area>=0.0)||!std::isfinite(s.area)) throw std::invalid_argument("invalid resolved surface area");
        const auto c=nearest_cell(mesh,s.position);
        Vec3 traction=s.normal*(-pressure[c]);
        if(!stress.empty()) traction+=stress[c].traction(s.normal);
        const Vec3 df=traction*s.area;
        result.force+=df;
        result.torque+=cross3(s.position-center,df);
    }
    return result;
}

void apply_resolved_surface_load(RigidBody& body,const ResolvedSurfaceLoad& load) noexcept {
    body.add_force(load.force); body.add_torque(load.torque);
}

UnresolvedCfdDemResult unresolved_cfd_dem_coupling(
    const cfd::fvm::PolyMesh& mesh,std::span<const RigidBody> bodies,std::span<const SphereShape> spheres,
    std::span<const cfd::fvm::Vec3> fluid_velocity,std::span<const cfd::fvm::Vec3> pressure_gradient,
    std::span<const cfd::fvm::Vec3> vorticity,const UnresolvedCfdDemConfig& cfg) {
    if(fluid_velocity.size()!=mesh.cell_count()) throw std::invalid_argument("CFD/DEM velocity field size mismatch");
    if(!pressure_gradient.empty()&&pressure_gradient.size()!=mesh.cell_count()) throw std::invalid_argument("CFD/DEM pressure-gradient size mismatch");
    if(!vorticity.empty()&&vorticity.size()!=mesh.cell_count()) throw std::invalid_argument("CFD/DEM vorticity size mismatch");
    if(!(cfg.fluid_density>0.0)||!(cfg.dynamic_viscosity>0.0)||!(cfg.minimum_void_fraction>0.0&&cfg.minimum_void_fraction<=1.0)||cfg.void_drag_exponent<0.0)
        throw std::invalid_argument("invalid unresolved CFD/DEM configuration");

    UnresolvedCfdDemResult out;
    out.solid_volume_fraction.assign(mesh.cell_count(),0.0);
    out.void_fraction.assign(mesh.cell_count(),1.0);
    out.particle_force.assign(spheres.size(),{});
    out.fluid_momentum_source.assign(mesh.cell_count(),{});
    std::vector<std::size_t> cell_of_sphere(spheres.size());

    for(std::size_t i=0;i<spheres.size();++i){
        const auto& s=spheres[i];
        if(s.body>=bodies.size()||!(s.radius>0.0)) throw std::invalid_argument("invalid unresolved CFD/DEM sphere");
        const auto c=nearest_cell(mesh,bodies[s.body].state.position); cell_of_sphere[i]=c;
        out.solid_volume_fraction[c]+=sphere_volume(s.radius)/mesh.cells()[c].volume;
    }
    for(std::size_t c=0;c<mesh.cell_count();++c){
        const double alpha=std::clamp(out.solid_volume_fraction[c],0.0,1.0-cfg.minimum_void_fraction);
        out.solid_volume_fraction[c]=alpha; out.void_fraction[c]=std::max(cfg.minimum_void_fraction,1.0-alpha);
    }

    const double nu=cfg.dynamic_viscosity/cfg.fluid_density;
    for(std::size_t i=0;i<spheres.size();++i){
        const auto& s=spheres[i]; const auto& b=bodies[s.body]; const auto c=cell_of_sphere[i];
        const Vec3 uf=to_mb(fluid_velocity[c]); const Vec3 rel=uf-b.state.linear_velocity;
        const double d=2.0*s.radius, vp=sphere_volume(s.radius);
        Vec3 force=schiller_naumann_drag(rel,d,cfg.fluid_density,cfg.dynamic_viscosity,out.void_fraction[c],cfg.void_drag_exponent);
        if(cfg.include_pressure_gradient&&!pressure_gradient.empty()) force-=to_mb(pressure_gradient[c])*vp;
        if(cfg.include_saffman_lift&&!vorticity.empty()) force+=saffman(rel,to_mb(vorticity[c]),cfg.fluid_density,nu,d,cfg.saffman_coefficient);
        out.particle_force[i]=force;
        out.fluid_momentum_source[c]-=force/mesh.cells()[c].volume;
    }
    return out;
}

void apply_unresolved_particle_forces(std::span<RigidBody> bodies,std::span<const SphereShape> spheres,
                                      std::span<const Vec3> forces) {
    if(forces.size()!=spheres.size()) throw std::invalid_argument("unresolved particle-force size mismatch");
    for(std::size_t i=0;i<spheres.size();++i){
        if(spheres[i].body>=bodies.size()) throw std::invalid_argument("unresolved sphere body index out of range");
        bodies[spheres[i].body].add_force(forces[i]);
    }
}

} // namespace cfd::multibody
