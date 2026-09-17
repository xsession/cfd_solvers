#include "cfd/particle/transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cfd::particle {
namespace {
constexpr std::size_t npos=static_cast<std::size_t>(-1);

double dot(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
Vec3 add(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
Vec3 scale(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
bool finite(Vec3 v){return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);} 

void validate_world(const TransportWorld& world){
    if(world.particles.empty())throw std::invalid_argument("transport world requires at least one particle definition");
    if(world.materials.empty())throw std::invalid_argument("transport world requires at least one material");
    if(world.regions.empty())throw std::invalid_argument("transport world requires at least one region");
    for(const auto& p:world.particles)if(p.name.empty()||!(p.mass_kg>0.0))throw std::invalid_argument("invalid transport particle definition");
    for(const auto& m:world.materials)if(m.name.empty()||m.density_kg_per_m3<0.0||m.production_cut_energy_ev<0.0)throw std::invalid_argument("invalid transport material");
    for(const auto& r:world.regions){
        if(r.material_index>=world.materials.size())throw std::invalid_argument("transport region references invalid material");
        if(r.sensitive_detector_index!=npos&&r.sensitive_detector_index>=world.sensitive_detectors.size())throw std::invalid_argument("transport region references invalid sensitive detector");
        if(!(r.box.maximum_m.x>r.box.minimum_m.x)||!(r.box.maximum_m.y>r.box.minimum_m.y)||!(r.box.maximum_m.z>r.box.minimum_m.z))throw std::invalid_argument("invalid transport region box");
    }
    for(const auto& detector:world.sensitive_detectors)if(detector.name.empty()||detector.minimum_energy_deposit_ev<0.0)throw std::invalid_argument("invalid transport sensitive detector");
    for(const auto& process:world.processes){
        if(process.name.empty())throw std::invalid_argument("transport process requires a name");
        if(process.particle_index!=npos&&process.particle_index>=world.particles.size())throw std::invalid_argument("transport process references invalid particle");
        if(process.material_index!=npos&&process.material_index>=world.materials.size())throw std::invalid_argument("transport process references invalid material");
        if(process.stopping_power_ev_per_m<0.0||process.physical_interaction_length_m<0.0||process.energy_loss_fraction<0.0||process.secondary_energy_fraction<0.0||process.absorption_threshold_ev<0.0)throw std::invalid_argument("invalid transport process coefficient");
    }
}

void validate_config(const TransportConfig& config){
    if(!(config.max_step_m>0.0)||!(config.boundary_epsilon_m>=0.0)||!(config.light_speed_m_per_s>0.0)||config.energy_cut_ev<0.0)throw std::invalid_argument("invalid transport configuration");
}

} // namespace

Vec3 normalize_direction(Vec3 direction){
    if(!finite(direction))throw std::invalid_argument("invalid direction vector");
    const double n=std::sqrt(dot(direction,direction));
    if(!(n>0.0))throw std::invalid_argument("zero transport direction");
    return {direction.x/n,direction.y/n,direction.z/n};
}

bool contains(const AxisAlignedParticleBox& box,Vec3 p,double tolerance_m){
    return p.x>=box.minimum_m.x-tolerance_m&&p.x<=box.maximum_m.x+tolerance_m&&
           p.y>=box.minimum_m.y-tolerance_m&&p.y<=box.maximum_m.y+tolerance_m&&
           p.z>=box.minimum_m.z-tolerance_m&&p.z<=box.maximum_m.z+tolerance_m;
}

std::size_t locate_region(const TransportWorld& world,Vec3 position,double tolerance_m){
    for(std::size_t i=0;i<world.regions.size();++i)if(contains(world.regions[i].box,position,tolerance_m))return i;
    return npos;
}

double distance_to_box_boundary(const AxisAlignedParticleBox& box,Vec3 p,Vec3 d){
    d=normalize_direction(d);double distance=std::numeric_limits<double>::infinity();
    auto candidate=[&](double dir,double pos,double plane){
        if(dir>0.0){const double t=(plane-pos)/dir;if(t>=0.0)distance=std::min(distance,t);} 
        else if(dir<0.0){const double t=(plane-pos)/dir;if(t>=0.0)distance=std::min(distance,t);} 
    };
    candidate(d.x,p.x,d.x>0.0?box.maximum_m.x:box.minimum_m.x);
    candidate(d.y,p.y,d.y>0.0?box.maximum_m.y:box.minimum_m.y);
    candidate(d.z,p.z,d.z>0.0?box.maximum_m.z:box.minimum_m.z);
    return distance;
}

bool process_applies(const TransportProcess& process,const TransportTrack& track,std::size_t material_index){
    return (process.particle_index==npos||process.particle_index==track.particle_index)&&
           (process.material_index==npos||process.material_index==material_index);
}


TransportWorld apply_physics_list(TransportWorld world,const TransportPhysicsList& physics_list){
    if(physics_list.name.empty())throw std::invalid_argument("transport physics list requires a name");
    if(physics_list.production_cut_energy_ev>=0.0){
        for(auto& material:world.materials)material.production_cut_energy_ev=physics_list.production_cut_energy_ev;
    }
    for(const auto& process:physics_list.processes)world.processes.push_back(process);
    validate_world(world);
    return world;
}

double TransportRandom::uniform_open01(){
    state=state*6364136223846793005ULL+1442695040888963407ULL;
    const std::uint64_t mantissa=(state>>11)|1ULL;
    const double u=static_cast<double>(mantissa)*(1.0/9007199254740992.0);
    return std::min(std::max(u,std::numeric_limits<double>::min()),std::nextafter(1.0,0.0));
}

double sample_exponential_interaction_length_m(double mean_free_path_m,TransportRandom& rng){
    if(!(mean_free_path_m>0.0)||!std::isfinite(mean_free_path_m))throw std::invalid_argument("invalid mean free path");
    const double u=rng.uniform_open01();
    return -mean_free_path_m*std::log(u);
}

namespace {
AxisAlignedParticleBox union_box(AxisAlignedParticleBox a,AxisAlignedParticleBox b){
    return {{std::min(a.minimum_m.x,b.minimum_m.x),std::min(a.minimum_m.y,b.minimum_m.y),std::min(a.minimum_m.z,b.minimum_m.z)},
            {std::max(a.maximum_m.x,b.maximum_m.x),std::max(a.maximum_m.y,b.maximum_m.y),std::max(a.maximum_m.z,b.maximum_m.z)}};
}
Vec3 center(AxisAlignedParticleBox box){return scale(add(box.minimum_m,box.maximum_m),0.5);} 
}

TransportRegionBvh build_region_bvh(const TransportWorld& world,std::size_t leaf_size){
    validate_world(world);if(leaf_size==0U)throw std::invalid_argument("BVH leaf size must be positive");
    TransportRegionBvh bvh;bvh.region_indices.resize(world.regions.size());std::iota(bvh.region_indices.begin(),bvh.region_indices.end(),0U);
    const auto build=[&](auto&& self,std::size_t first,std::size_t count)->std::size_t{
        TransportBvhNode node;node.first=first;node.count=count;AxisAlignedParticleBox bounds=world.regions[bvh.region_indices[first]].box;
        for(std::size_t i=1;i<count;++i)bounds=union_box(bounds,world.regions[bvh.region_indices[first+i]].box);
        node.box=bounds;const std::size_t node_index=bvh.nodes.size();bvh.nodes.push_back(node);
        if(count>leaf_size){
            const Vec3 span={bounds.maximum_m.x-bounds.minimum_m.x,bounds.maximum_m.y-bounds.minimum_m.y,bounds.maximum_m.z-bounds.minimum_m.z};
            const int axis=(span.x>=span.y&&span.x>=span.z)?0:(span.y>=span.z?1:2);
            auto coord=[&](std::size_t region_index){const Vec3 c=center(world.regions[region_index].box);return axis==0?c.x:(axis==1?c.y:c.z);};
            const std::size_t mid=first+count/2U;
            std::nth_element(bvh.region_indices.begin()+static_cast<std::ptrdiff_t>(first),bvh.region_indices.begin()+static_cast<std::ptrdiff_t>(mid),bvh.region_indices.begin()+static_cast<std::ptrdiff_t>(first+count),[&](std::size_t a,std::size_t b){return coord(a)<coord(b);});
            bvh.nodes[node_index].left=self(self,first,mid-first);
            bvh.nodes[node_index].right=self(self,mid,first+count-mid);
            bvh.nodes[node_index].count=0U;
        }
        return node_index;
    };
    if(!bvh.region_indices.empty())(void)build(build,0U,bvh.region_indices.size());
    return bvh;
}

std::size_t locate_region_bvh(const TransportWorld& world,const TransportRegionBvh& bvh,Vec3 position,double tolerance_m){
    if(bvh.nodes.empty())return npos;
    const auto query=[&](auto&& self,std::size_t node_index)->std::size_t{
        const auto& node=bvh.nodes[node_index];if(!contains(node.box,position,tolerance_m))return npos;
        if(node.count>0U){
            for(std::size_t i=0;i<node.count;++i){const std::size_t region=bvh.region_indices[node.first+i];if(contains(world.regions[region].box,position,tolerance_m))return region;}
            return npos;
        }
        const std::size_t left=(node.left==npos)?npos:self(self,node.left);if(left!=npos)return left;
        return node.right==npos?npos:self(self,node.right);
    };
    return query(query,0U);
}

TransportHitCollection collect_transport_hits(const TransportWorld& world,std::span<const TransportStep> steps){
    validate_world(world);TransportHitCollection collection;
    for(const auto& step:steps){
        if(step.region_index>=world.regions.size()||step.material_index>=world.materials.size())continue;
        const std::size_t detector=world.regions[step.region_index].sensitive_detector_index;
        if(detector==npos||detector>=world.sensitive_detectors.size())continue;
        const auto& definition=world.sensitive_detectors[detector];
        if(step.energy_deposit_ev<definition.minimum_energy_deposit_ev&&!(definition.record_zero_deposit_steps&&step.step_length_m>0.0))continue;
        TransportHit hit;hit.detector_index=detector;hit.region_index=step.region_index;hit.material_index=step.material_index;hit.position_m=scale(add(step.pre_position_m,step.post_position_m),0.5);hit.energy_deposit_ev=step.energy_deposit_ev;hit.step_length_m=step.step_length_m;hit.global_time_s=step.global_time_s;hit.limiter=step.limiter;
        if(step.limiting_process_index<world.processes.size())hit.process_name=world.processes[step.limiting_process_index].name;
        collection.total_energy_deposit_ev+=hit.energy_deposit_ev;
        collection.total_weighted_energy_deposit_ev+=hit.energy_deposit_ev*hit.weight;
        collection.hits.push_back(std::move(hit));
    }
    return collection;
}

TransportDoseGrid3D make_transport_dose_grid(std::size_t nx,std::size_t ny,std::size_t nz,AxisAlignedParticleBox box,double density_kg_per_m3){
    if(nx==0U||ny==0U||nz==0U||!(density_kg_per_m3>0.0)||!contains(box,box.minimum_m)||!(box.maximum_m.x>box.minimum_m.x)||!(box.maximum_m.y>box.minimum_m.y)||!(box.maximum_m.z>box.minimum_m.z))throw std::invalid_argument("invalid transport dose grid");
    TransportDoseGrid3D grid;grid.nx=nx;grid.ny=ny;grid.nz=nz;grid.box=box;grid.density_kg_per_m3=density_kg_per_m3;grid.energy_deposit_ev.assign(nx*ny*nz,0.0);grid.dose_gy.assign(nx*ny*nz,0.0);return grid;
}

void score_hits_to_dose_grid(const TransportHitCollection& hits,TransportDoseGrid3D& grid){
    if(grid.nx==0U||grid.ny==0U||grid.nz==0U||grid.energy_deposit_ev.size()!=grid.nx*grid.ny*grid.nz||grid.dose_gy.size()!=grid.energy_deposit_ev.size()||!(grid.density_kg_per_m3>0.0))throw std::invalid_argument("invalid transport dose grid storage");
    const double lx=grid.box.maximum_m.x-grid.box.minimum_m.x,ly=grid.box.maximum_m.y-grid.box.minimum_m.y,lz=grid.box.maximum_m.z-grid.box.minimum_m.z;if(!(lx>0.0)||!(ly>0.0)||!(lz>0.0))throw std::invalid_argument("invalid transport dose grid box");
    for(const auto& hit:hits.hits){
        if(!contains(grid.box,hit.position_m,0.0))continue;
        const double fx=(hit.position_m.x-grid.box.minimum_m.x)/lx,fy=(hit.position_m.y-grid.box.minimum_m.y)/ly,fz=(hit.position_m.z-grid.box.minimum_m.z)/lz;
        const std::size_t ix=std::min(grid.nx-1U,static_cast<std::size_t>(std::floor(fx*static_cast<double>(grid.nx))));
        const std::size_t iy=std::min(grid.ny-1U,static_cast<std::size_t>(std::floor(fy*static_cast<double>(grid.ny))));
        const std::size_t iz=std::min(grid.nz-1U,static_cast<std::size_t>(std::floor(fz*static_cast<double>(grid.nz))));
        grid.energy_deposit_ev[(iz*grid.ny+iy)*grid.nx+ix]+=hit.energy_deposit_ev*hit.weight;
    }
    constexpr double electron_volt_j=1.602176634e-19;
    const double voxel_mass=grid.density_kg_per_m3*(lx/static_cast<double>(grid.nx))*(ly/static_cast<double>(grid.ny))*(lz/static_cast<double>(grid.nz));
    for(std::size_t i=0;i<grid.energy_deposit_ev.size();++i)grid.dose_gy[i]=grid.energy_deposit_ev[i]*electron_volt_j/voxel_mass;
}

std::vector<double> dose_rate_to_sar_w_per_kg(const TransportDoseGrid3D& grid,double exposure_time_s){
    if(!(exposure_time_s>0.0)||grid.dose_gy.empty())throw std::invalid_argument("invalid dose-rate conversion");
    std::vector<double> sar(grid.dose_gy.size());for(std::size_t i=0;i<sar.size();++i){if(grid.dose_gy[i]<0.0||!std::isfinite(grid.dose_gy[i]))throw std::invalid_argument("invalid dose grid value");sar[i]=grid.dose_gy[i]/exposure_time_s;}return sar;
}

std::vector<double> project_dose_grid_to_pennes2d_sar(const TransportDoseGrid3D& grid,double exposure_time_s,std::size_t nx,std::size_t ny){
    if(nx==0U||ny==0U)throw std::invalid_argument("invalid Pennes projection size");
    const auto sar3=dose_rate_to_sar_w_per_kg(grid,exposure_time_s);std::vector<double> out(nx*ny,0.0),count(nx*ny,0.0);
    for(std::size_t iz=0;iz<grid.nz;++iz)for(std::size_t iy=0;iy<grid.ny;++iy)for(std::size_t ix=0;ix<grid.nx;++ix){
        const std::size_t src=(iz*grid.ny+iy)*grid.nx+ix;const std::size_t tx=std::min(nx-1U,ix*nx/grid.nx),ty=std::min(ny-1U,iy*ny/grid.ny),dst=ty*nx+tx;out[dst]+=sar3[src];count[dst]+=1.0;
    }
    for(std::size_t i=0;i<out.size();++i)if(count[i]>0.0)out[i]/=count[i];return out;
}

TransportStep transport_one_step(const TransportWorld& world,TransportTrack& track,const TransportConfig& config){
    validate_world(world);validate_config(config);
    if(track.particle_index>=world.particles.size())throw std::invalid_argument("track references invalid particle");
    if(track.status!=TrackStatus::alive)return {};
    track.direction=normalize_direction(track.direction);
    if(!(track.kinetic_energy_ev>=0.0)||!finite(track.position_m)||!std::isfinite(track.weight))throw std::invalid_argument("invalid transport track");

    const std::size_t region_index=locate_region(world,track.position_m,config.boundary_epsilon_m);
    if(region_index==npos){track.status=TrackStatus::escaped;return {};}
    const auto& region=world.regions[region_index];const std::size_t material_index=region.material_index;
    const auto& material=world.materials[material_index];

    TransportStep step;step.pre_position_m=track.position_m;step.pre_energy_ev=track.kinetic_energy_ev;step.material_index=material_index;step.region_index=region_index;step.global_time_s=track.time_s;
    double selected=config.max_step_m;step.limiter=StepLimiter::user;
    const double boundary=distance_to_box_boundary(region.box,track.position_m,track.direction);
    if(boundary<selected){selected=boundary;step.limiter=StepLimiter::geometry_boundary;}
    for(std::size_t i=0;i<world.processes.size();++i){
        const auto& process=world.processes[i];
        if(!process_applies(process,track,material_index))continue;
        if(process.kind==TransportProcessKind::discrete_interaction&&process.physical_interaction_length_m>0.0&&process.physical_interaction_length_m<selected){selected=process.physical_interaction_length_m;step.limiter=StepLimiter::discrete_process;step.limiting_process_index=i;}
        if(process.kind==TransportProcessKind::absorber&&track.kinetic_energy_ev<=process.absorption_threshold_ev){selected=0.0;step.limiter=StepLimiter::energy_cut;step.limiting_process_index=i;}
    }
    selected=std::max(0.0,selected);
    step.step_length_m=selected;step.post_position_m=add(track.position_m,scale(track.direction,selected));

    double continuous_loss=0.0;
    for(const auto& process:world.processes){
        if(process.kind!=TransportProcessKind::continuous_energy_loss||!process_applies(process,track,material_index))continue;
        continuous_loss+=process.stopping_power_ev_per_m*selected;
    }
    continuous_loss=std::min(track.kinetic_energy_ev,continuous_loss);
    track.kinetic_energy_ev-=continuous_loss;step.energy_deposit_ev+=continuous_loss;

    if(step.limiter==StepLimiter::discrete_process&&step.limiting_process_index<world.processes.size()){
        const auto& process=world.processes[step.limiting_process_index];
        const double loss=std::min(track.kinetic_energy_ev,track.kinetic_energy_ev*std::min(1.0,process.energy_loss_fraction));
        track.kinetic_energy_ev-=loss;step.energy_deposit_ev+=loss;
        const double secondary_energy=track.kinetic_energy_ev*process.secondary_energy_fraction;
        if(secondary_energy>material.production_cut_energy_ev&&secondary_energy>0.0){
            TransportSecondary secondary;secondary.creator_process=process.name;secondary.track=track;secondary.track.position_m=step.post_position_m;secondary.track.kinetic_energy_ev=secondary_energy;secondary.track.track_length_m=0.0;secondary.track.time_s=track.time_s+selected/config.light_speed_m_per_s;secondary.track.weight=track.weight;secondary.track.status=TrackStatus::alive;
            step.secondaries.push_back(secondary);track.kinetic_energy_ev=std::max(0.0,track.kinetic_energy_ev-secondary_energy);
        }
    }

    track.position_m=step.post_position_m;track.track_length_m+=selected;track.time_s+=selected/config.light_speed_m_per_s;
    step.global_time_s=track.time_s;step.post_energy_ev=track.kinetic_energy_ev;
    if(step.limiter==StepLimiter::energy_cut||track.kinetic_energy_ev<=std::max(config.energy_cut_ev,material.production_cut_energy_ev)){track.status=TrackStatus::stopped;}
    if(step.limiter==StepLimiter::geometry_boundary){
        const Vec3 nudged=add(track.position_m,scale(track.direction,config.boundary_epsilon_m*4.0));
        if(locate_region(world,nudged,config.boundary_epsilon_m)==npos)track.status=TrackStatus::escaped;
        else track.position_m=nudged;
    }
    return step;
}

TransportRunResult transport_track(const TransportWorld& world,TransportTrack primary,const TransportConfig& config,std::size_t max_steps){
    validate_world(world);validate_config(config);if(max_steps==0U)throw std::invalid_argument("transport requires positive max steps");
    primary.direction=normalize_direction(primary.direction);
    TransportRunResult result;result.primary=primary;
    for(std::size_t step_index=0;step_index<max_steps&&result.primary.status==TrackStatus::alive;++step_index){
        const auto step=transport_one_step(world,result.primary,config);
        if(step.step_length_m==0.0&&step.secondaries.empty()&&result.primary.status==TrackStatus::alive)break;
        result.scoring.total_energy_deposit_ev+=step.energy_deposit_ev;result.scoring.total_track_length_m+=step.step_length_m;result.scoring.secondaries+=step.secondaries.size();if(step.limiter==StepLimiter::geometry_boundary)++result.scoring.boundary_crossings;++result.scoring.steps;
        for(const auto& secondary:step.secondaries)result.secondaries.push_back(secondary);
        result.steps.push_back(step);
    }
    return result;
}


TransportStep transport_one_step_stochastic(const TransportWorld& world,TransportTrack& track,const TransportConfig& config,TransportRandom& rng){
    TransportWorld sampled=world;
    for(auto& process:sampled.processes){
        if(process.kind==TransportProcessKind::discrete_interaction&&process.physical_interaction_length_m>0.0){
            process.physical_interaction_length_m=sample_exponential_interaction_length_m(process.physical_interaction_length_m,rng);
        }
    }
    return transport_one_step(sampled,track,config);
}

TransportRunResult transport_track_stochastic(const TransportWorld& world,TransportTrack primary,const TransportConfig& config,TransportRandom& rng,std::size_t max_steps){
    validate_world(world);validate_config(config);if(max_steps==0U)throw std::invalid_argument("transport requires positive max steps");
    primary.direction=normalize_direction(primary.direction);TransportRunResult result;result.primary=primary;
    for(std::size_t step_index=0;step_index<max_steps&&result.primary.status==TrackStatus::alive;++step_index){
        const auto step=transport_one_step_stochastic(world,result.primary,config,rng);
        if(step.step_length_m==0.0&&step.secondaries.empty()&&result.primary.status==TrackStatus::alive)break;
        result.scoring.total_energy_deposit_ev+=step.energy_deposit_ev;result.scoring.total_track_length_m+=step.step_length_m;result.scoring.secondaries+=step.secondaries.size();if(step.limiter==StepLimiter::geometry_boundary)++result.scoring.boundary_crossings;++result.scoring.steps;
        for(const auto& secondary:step.secondaries)result.secondaries.push_back(secondary);
        result.steps.push_back(step);
    }
    return result;
}

} // namespace cfd::particle
