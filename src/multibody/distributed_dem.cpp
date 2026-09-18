#include "cfd/multibody/distributed_dem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::multibody {
namespace {

[[nodiscard]] bool finite(Vec3 v) noexcept {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

[[nodiscard]] std::pair<std::uint64_t,std::uint64_t> canonical_pair(std::uint64_t a,std::uint64_t b) {
    if(a==b) throw std::invalid_argument("distributed DEM pair requires distinct global IDs");
    return a<b ? std::pair{a,b} : std::pair{b,a};
}

[[nodiscard]] const DistributedDemParticle* find_particle(
    std::span<const DistributedDemParticle> particles,std::uint64_t id) noexcept {
    for(const auto& particle:particles) if(particle.global_id==id) return &particle;
    return nullptr;
}

[[nodiscard]] DistributedDemContactHistory* find_history(
    std::vector<DistributedDemContactHistory>& history,std::uint64_t a,std::uint64_t b) noexcept {
    for(auto& state:history) if(state.particle_a==a && state.particle_b==b) return &state;
    return nullptr;
}

void validate_particle(const DistributedDemParticle& particle) {
    if(!(particle.mass>0.0) || !std::isfinite(particle.mass) ||
       !(particle.radius>0.0) || !std::isfinite(particle.radius) ||
       !(particle.inertia_diagonal.x>0.0) || !(particle.inertia_diagonal.y>0.0) || !(particle.inertia_diagonal.z>0.0) ||
       !finite(particle.state.position) || !finite(particle.state.linear_velocity) || !finite(particle.state.angular_velocity)) {
        throw std::invalid_argument("invalid distributed DEM particle state");
    }
}

void validate_contact_model(const HertzMindlinContactModel& model) {
    const auto nonnegative=[](double value){return value>=0.0 && std::isfinite(value);};
    if(!nonnegative(model.normal_stiffness) || !nonnegative(model.normal_damping) ||
       !nonnegative(model.tangential_stiffness) || !nonnegative(model.tangential_damping) ||
       !nonnegative(model.friction) || !nonnegative(model.rolling_resistance) ||
       !nonnegative(model.cohesion_force)) {
        throw std::invalid_argument("invalid distributed DEM contact coefficients");
    }
}

void add_contribution(std::vector<DistributedDemForceContribution>& destination,
                      std::uint64_t id,int owner,Vec3 force,Vec3 torque) {
    for(auto& entry:destination) {
        if(entry.global_id==id) {
            if(entry.owner_rank!=owner) throw std::runtime_error("inconsistent distributed DEM force owner");
            entry.force+=force;
            entry.torque+=torque;
            return;
        }
    }
    destination.push_back({id,owner,force,torque});
}

[[nodiscard]] Vec3 point_velocity(const DistributedDemParticle& particle,Vec3 point) noexcept {
    return particle.state.linear_velocity + cross(particle.state.angular_velocity,point-particle.state.position);
}

void evaluate_pair(const DistributedDemParticle& a,
                   const DistributedDemParticle& b,
                   bool b_is_remote,
                   int local_rank,
                   double dt,
                   const HertzMindlinContactModel& model,
                   std::vector<DistributedDemContactHistory>& old_history,
                   std::vector<DistributedDemContactHistory>& next_history,
                   DistributedDemInteractionResult& result) {
    const Vec3 delta=a.state.position-b.state.position;
    const double distance_squared=norm_squared(delta);
    const double radius_sum=a.radius+b.radius;
    if(!(distance_squared<radius_sum*radius_sum)) return;
    if(!(distance_squared>1.0e-30) || !std::isfinite(distance_squared)) {
        throw std::runtime_error("distributed DEM overlapping particle centers are coincident/non-finite");
    }
    const double distance=std::sqrt(distance_squared);
    const Vec3 normal=delta/distance; // force direction on a
    const double penetration=radius_sum-distance;
    const double root=std::sqrt(std::max(0.0,penetration));
    const Vec3 contact_point=a.state.position-normal*(a.radius-0.5*penetration);
    const Vec3 ra=contact_point-a.state.position;
    const Vec3 rb=contact_point-b.state.position;
    const Vec3 relative_velocity=point_velocity(a,contact_point)-point_velocity(b,contact_point);
    const double normal_velocity=dot(relative_velocity,normal);
    const Vec3 tangent_velocity=relative_velocity-normal*normal_velocity;

    const auto [id_a,id_b]=canonical_pair(a.global_id,b.global_id);
    DistributedDemContactHistory state;
    state.particle_a=id_a;
    state.particle_b=id_b;
    if(auto* previous=find_history(old_history,id_a,id_b)) state=*previous;
    state.tangential_displacement+=tangent_velocity*dt;
    state.tangential_displacement-=normal*dot(state.tangential_displacement,normal);
    ++state.age;

    double normal_force=model.normal_stiffness*penetration*root
        -model.normal_damping*root*normal_velocity-model.cohesion_force;
    normal_force=std::max(-model.cohesion_force,normal_force);
    const double kt=model.tangential_stiffness*root;
    const double ct=model.tangential_damping*root;
    Vec3 tangent_force=-(state.tangential_displacement*kt+tangent_velocity*ct);
    const double tangent_magnitude=norm(tangent_force);
    const double friction_limit=model.friction*std::abs(normal_force);
    if(tangent_magnitude>friction_limit && tangent_magnitude>1.0e-30) {
        tangent_force*=friction_limit/tangent_magnitude;
        if(kt>1.0e-30) state.tangential_displacement=-(tangent_force+tangent_velocity*ct)/kt;
    }
    next_history.push_back(state);

    const Vec3 force=normal*normal_force+tangent_force;
    Vec3 torque_a=cross(ra,force);
    Vec3 torque_b=cross(rb,-force);
    const Vec3 relative_omega=a.state.angular_velocity-b.state.angular_velocity;
    const double omega_magnitude=norm(relative_omega);
    if(omega_magnitude>1.0e-30 && model.rolling_resistance>0.0) {
        const double effective_radius=a.radius*b.radius/(a.radius+b.radius);
        const Vec3 rolling=relative_omega*(-model.rolling_resistance*std::abs(normal_force)*effective_radius/omega_magnitude);
        torque_a+=rolling;
        torque_b-=rolling;
    }
    add_contribution(result.local,a.global_id,local_rank,force,torque_a);
    if(b_is_remote) add_contribution(result.remote,b.global_id,b.owner_rank,-force,torque_b);
    else add_contribution(result.local,b.global_id,local_rank,-force,torque_b);
    ++result.contacts;
}

void validate_bond_model(const BondedParticleModel& model) {
    const auto nonnegative=[](double value){return value>=0.0 && std::isfinite(value);};
    if(!nonnegative(model.normal_stiffness) || !nonnegative(model.shear_stiffness) ||
       !nonnegative(model.normal_damping) || !nonnegative(model.shear_damping) ||
       !(model.tensile_failure_force>0.0) || !std::isfinite(model.tensile_failure_force) ||
       !(model.shear_failure_force>0.0) || !std::isfinite(model.shear_failure_force) ||
       !(model.damage_onset_ratio>=0.0 && model.damage_onset_ratio<1.0) || !std::isfinite(model.damage_onset_ratio)) {
        throw std::invalid_argument("invalid distributed DEM bond constitutive parameters");
    }
}

} // namespace

void DemSlabDecomposition::validate() const {
    if(!(x_max>x_min) || !std::isfinite(x_min) || !std::isfinite(x_max)) throw std::invalid_argument("DEM slab bounds must be finite and ordered");
    if(ranks<=0) throw std::invalid_argument("DEM slab decomposition requires at least one rank");
}

double DemSlabDecomposition::slab_width() const {
    validate();
    return (x_max-x_min)/static_cast<double>(ranks);
}

int DemSlabDecomposition::owner_rank(double x) const {
    validate();
    if(!std::isfinite(x)) throw std::invalid_argument("DEM particle coordinate must be finite");
    if(x<=x_min) return 0;
    if(x>=x_max) return ranks-1;
    const double relative=(x-x_min)/(x_max-x_min);
    const auto rank=static_cast<int>(std::floor(relative*static_cast<double>(ranks)));
    return std::clamp(rank,0,ranks-1);
}

double DemSlabDecomposition::slab_min(int rank) const {
    validate();
    if(rank<0 || rank>=ranks) throw std::out_of_range("DEM slab rank out of range");
    return x_min+static_cast<double>(rank)*slab_width();
}

double DemSlabDecomposition::slab_max(int rank) const {
    validate();
    if(rank<0 || rank>=ranks) throw std::out_of_range("DEM slab rank out of range");
    return rank+1==ranks ? x_max : x_min+static_cast<double>(rank+1)*slab_width();
}

DistributedDemParticle make_distributed_dem_particle(std::uint64_t global_id,const RigidBody& body,const SphereShape& sphere,int owner_rank) {
    if(sphere.body==invalid_body) throw std::invalid_argument("distributed DEM sphere must reference a body");
    if(!(sphere.radius>0.0) || !std::isfinite(sphere.radius)) throw std::invalid_argument("distributed DEM sphere radius must be finite and positive");
    body.validate();
    DistributedDemParticle particle;
    particle.global_id=global_id;
    particle.state=body.state;
    particle.mass=body.mass;
    particle.inertia_diagonal=body.inertia_diagonal;
    particle.radius=sphere.radius;
    particle.owner_rank=owner_rank;
    return particle;
}

DistributedDemBond make_distributed_dem_bond(const DistributedDemParticle& a,const DistributedDemParticle& b,BondedParticleModel model,double rest_length) {
    validate_particle(a); validate_particle(b); validate_bond_model(model);
    const auto [id_a,id_b]=canonical_pair(a.global_id,b.global_id);
    const double distance=norm(b.state.position-a.state.position);
    if(rest_length==0.0) rest_length=distance;
    if(!(rest_length>0.0) || !std::isfinite(rest_length)) throw std::invalid_argument("distributed DEM bond rest length must be finite and positive");
    DistributedDemBond bond;
    bond.particle_a=id_a;
    bond.particle_b=id_b;
    bond.rest_length=rest_length;
    bond.model=model;
    return bond;
}

DemExchangePlan plan_dem_slab_exchange(std::span<const DistributedDemParticle> particles,const DemSlabDecomposition& decomposition,int local_rank,double ghost_width) {
    decomposition.validate();
    if(local_rank<0 || local_rank>=decomposition.ranks) throw std::out_of_range("local DEM rank out of range");
    if(!(ghost_width>=0.0) || !std::isfinite(ghost_width) || ghost_width>decomposition.slab_width()) throw std::invalid_argument("DEM ghost width must be finite and no larger than one slab");
    DemExchangePlan plan;
    plan.migrate_indices_by_rank.resize(static_cast<std::size_t>(decomposition.ranks));
    plan.ghost_indices_by_rank.resize(static_cast<std::size_t>(decomposition.ranks));
    const double lower=decomposition.slab_min(local_rank);
    const double upper=decomposition.slab_max(local_rank);
    for(std::size_t i=0;i<particles.size();++i){
        const auto& particle=particles[i];
        validate_particle(particle);
        const int owner=decomposition.owner_rank(particle.state.position.x);
        if(owner!=local_rank){
            plan.migrate_indices_by_rank[static_cast<std::size_t>(owner)].push_back(i);
            continue;
        }
        const double reach=particle.radius+ghost_width;
        if(local_rank>0 && particle.state.position.x-reach<lower){
            plan.ghost_indices_by_rank[static_cast<std::size_t>(local_rank-1)].push_back(i);
        }
        if(local_rank+1<decomposition.ranks && particle.state.position.x+reach>upper){
            plan.ghost_indices_by_rank[static_cast<std::size_t>(local_rank+1)].push_back(i);
        }
    }
    return plan;
}

DistributedDemInteractionResult evaluate_distributed_dem_contacts(
    std::span<const DistributedDemParticle> owned,std::span<const DistributedDemParticle> ghosts,
    std::vector<DistributedDemContactHistory>& history,int local_rank,double dt,const HertzMindlinContactModel& model) {
    if(!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("distributed DEM contact dt must be finite and positive");
    validate_contact_model(model);
    for(const auto& p:owned){validate_particle(p);if(p.owner_rank!=local_rank)throw std::invalid_argument("owned DEM particle has wrong owner rank");}
    for(const auto& p:ghosts) validate_particle(p);

    DistributedDemInteractionResult result;
    std::vector<DistributedDemContactHistory> next_history;
    next_history.reserve(history.size()+owned.size());

    for(std::size_t i=0;i<owned.size();++i){
        for(std::size_t j=i+1;j<owned.size();++j){
            const auto* a=&owned[i]; const auto* b=&owned[j];
            if(a->global_id>b->global_id) std::swap(a,b);
            evaluate_pair(*a,*b,false,local_rank,dt,model,history,next_history,result);
        }
    }
    for(const auto& local:owned){
        for(const auto& ghost:ghosts){
            if(local.global_id>=ghost.global_id) continue; // lower-ID owner computes pair exactly once
            evaluate_pair(local,ghost,true,local_rank,dt,model,history,next_history,result);
        }
    }
    history=std::move(next_history);
    return result;
}

DistributedDemInteractionResult evaluate_distributed_dem_bonds(
    std::span<const DistributedDemParticle> owned,std::span<const DistributedDemParticle> ghosts,
    std::vector<DistributedDemBond>& bonds,int local_rank,double dt) {
    if(!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("distributed DEM bond dt must be finite and positive");
    DistributedDemInteractionResult result;
    for(auto& bond:bonds){
        validate_bond_model(bond.model);
        if(!(bond.particle_a<bond.particle_b)) throw std::invalid_argument("distributed DEM bond IDs must be canonical");
        if(bond.broken){++result.broken_bonds;continue;}
        const auto* a=find_particle(owned,bond.particle_a);
        if(!a) continue; // record belongs on another rank after ownership migration
        if(a->owner_rank!=local_rank) throw std::runtime_error("distributed DEM bond owner mismatch");
        const auto* b=find_particle(owned,bond.particle_b);
        bool remote=false;
        if(!b){b=find_particle(ghosts,bond.particle_b);remote=b!=nullptr;}
        if(!b) continue; // partner not in current halo; caller can enlarge ghost width
        const Vec3 delta=b->state.position-a->state.position;
        const double distance=norm(delta);
        if(!(distance>1.0e-15) || !std::isfinite(distance)) throw std::runtime_error("distributed bonded particle centers are coincident/non-finite");
        const Vec3 n=delta/distance;
        const Vec3 relative_velocity=b->state.linear_velocity-a->state.linear_velocity;
        const double vn=dot(relative_velocity,n);
        const Vec3 vt=relative_velocity-n*vn;
        bond.tangential_displacement+=vt*dt;
        bond.tangential_displacement-=n*dot(bond.tangential_displacement,n);
        const double extension=distance-bond.rest_length;
        const double normal_elastic=bond.model.normal_stiffness*extension;
        const Vec3 shear_elastic=bond.tangential_displacement*bond.model.shear_stiffness;
        const double tensile_ratio=std::max(0.0,normal_elastic)/bond.model.tensile_failure_force;
        const double shear_ratio=norm(shear_elastic)/bond.model.shear_failure_force;
        const double failure_ratio=std::max(tensile_ratio,shear_ratio);
        if(failure_ratio>bond.model.damage_onset_ratio){
            const double target=std::clamp((failure_ratio-bond.model.damage_onset_ratio)/(1.0-bond.model.damage_onset_ratio),0.0,1.0);
            bond.damage=std::max(bond.damage,target);
        }
        if(failure_ratio>=1.0 || bond.damage>=1.0){bond.damage=1.0;bond.broken=true;++result.broken_bonds;continue;}
        ++result.active_bonds;
        const double intact=1.0-bond.damage;
        const double normal_force=intact*normal_elastic+bond.model.normal_damping*vn;
        const Vec3 shear_force=shear_elastic*intact+vt*bond.model.shear_damping;
        const Vec3 force_on_a=n*normal_force+shear_force;
        add_contribution(result.local,a->global_id,local_rank,force_on_a,{});
        if(remote) add_contribution(result.remote,b->global_id,b->owner_rank,-force_on_a,{});
        else add_contribution(result.local,b->global_id,local_rank,-force_on_a,{});
    }
    return result;
}

void accumulate_distributed_dem_forces(std::vector<DistributedDemForceContribution>& destination,
                                        std::span<const DistributedDemForceContribution> source) {
    for(const auto& contribution:source) add_contribution(destination,contribution.global_id,contribution.owner_rank,contribution.force,contribution.torque);
}

void integrate_distributed_dem_particles(std::span<DistributedDemParticle> owned,
                                         std::span<const DistributedDemForceContribution> contributions,
                                         Vec3 gravity,double dt) {
    if(!(dt>0.0) || !std::isfinite(dt) || !finite(gravity)) throw std::invalid_argument("invalid distributed DEM integration parameters");
    for(auto& particle:owned){
        validate_particle(particle);
        Vec3 force=gravity*particle.mass;
        Vec3 torque{};
        for(const auto& contribution:contributions){
            if(contribution.global_id==particle.global_id){force+=contribution.force;torque+=contribution.torque;}
        }
        particle.state.linear_velocity+=force*(dt/particle.mass);
        const Vec3 torque_local=inverse_rotate(particle.state.orientation,torque);
        Vec3 omega_local=inverse_rotate(particle.state.orientation,particle.state.angular_velocity);
        omega_local.x+=dt*torque_local.x/particle.inertia_diagonal.x;
        omega_local.y+=dt*torque_local.y/particle.inertia_diagonal.y;
        omega_local.z+=dt*torque_local.z/particle.inertia_diagonal.z;
        particle.state.angular_velocity=rotate(particle.state.orientation,omega_local);
        particle.state.position+=particle.state.linear_velocity*dt;
        integrate_orientation(particle.state.orientation,particle.state.angular_velocity,dt);
    }
}

} // namespace cfd::multibody
