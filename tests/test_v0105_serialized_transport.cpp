#include "cfd/particle/electromagnetic.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 
std::size_t flat(std::size_t ix,std::size_t iy,std::size_t iz,std::size_t nx,std::size_t ny){return (iz*ny+iy)*nx+ix;}

double cell_value(std::size_t gx,std::size_t gy,std::size_t gz){return 100.0*static_cast<double>(gx)+10.0*static_cast<double>(gy)+static_cast<double>(gz);} 

cfd::particle::PicParticle3D make_particle(double x,double y,double z,double vx=1.0,double vy=2.0,double vz=3.0){
    cfd::particle::PicParticle3D p;p.position_m={x,y,z};p.velocity_m_per_s={vx,vy,vz};p.charge_c=1.0e-15;p.mass_kg=1.0e-6;p.weight=2.0;return p;
}

struct Fixture {
    cfd::particle::PicDomainGrid3DConfig config;
    std::vector<cfd::particle::PicDomain3D> domains;
    std::vector<std::vector<double>> values;
    std::vector<std::vector<cfd::particle::PicParticle3D>> particles_by_domain;
};

Fixture make_fixture(){
    using namespace cfd::particle;
    Fixture f;f.config.global_nx=6U;f.config.global_ny=4U;f.config.global_nz=3U;f.config.domains_x=3U;f.config.domains_y=2U;f.config.domains_z=1U;f.config.length_x_m=6.0;f.config.length_y_m=4.0;f.config.length_z_m=3.0;f.domains=make_pic_domain_grid_3d(f.config);
    f.values.resize(f.domains.size());
    for(const auto& domain:f.domains){
        f.values[domain.domain_id].resize(domain.cells_x*domain.cells_y*domain.cells_z);
        for(std::size_t iz=0;iz<domain.cells_z;++iz)for(std::size_t iy=0;iy<domain.cells_y;++iy)for(std::size_t ix=0;ix<domain.cells_x;++ix){
            f.values[domain.domain_id][flat(ix,iy,iz,domain.cells_x,domain.cells_y)]=cell_value(domain.first_x+ix,domain.first_y+iy,domain.first_z+iz);
        }
    }
    f.particles_by_domain.resize(f.domains.size());
    f.particles_by_domain[0].push_back(make_particle(0.25,0.25,0.25));
    f.particles_by_domain[0].push_back(make_particle(2.25,0.25,0.25));
    f.particles_by_domain[1].push_back(make_particle(4.25,2.25,1.25,-1.0,0.5,2.5));
    f.particles_by_domain[f.domains.back().domain_id].push_back(make_particle(6.05,0.25,0.25));
    return f;
}

void serialized_roundtrip(){
    using namespace cfd::particle;
    auto fixture=make_fixture();
    const auto migration=pack_particle_migration_messages_3d(fixture.particles_by_domain,fixture.domains,fixture.config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(fixture.values,fixture.domains,fixture.config,1U);
    const auto topology=make_pic_rank_topology_3d(fixture.domains,3U,1U,1U);
    const auto envelopes=make_pic_transport_envelopes_3d(migration.messages,guard_messages,topology);
    const auto serialized=serialize_pic_transport_envelopes_3d(envelopes);
    const auto decoded=deserialize_pic_transport_envelopes_3d(serialized);
    check(decoded.size()==envelopes.size(),"serialized envelope count roundtrip");
    const auto original_diag=summarize_pic_transport_3d(envelopes);
    const auto serialized_diag=summarize_pic_serialized_transport_3d(serialized);
    check(serialized_diag.envelopes==original_diag.envelopes,"serialized diagnostics envelope count");
    check(serialized_diag.particle_payload_count==original_diag.particle_payload_count,"serialized particle payload count");
    check(serialized_diag.scalar_guard_value_count==original_diag.scalar_guard_value_count,"serialized scalar payload count");
    for(std::size_t i=0;i<envelopes.size();++i){
        check(decoded[i].kind==envelopes[i].kind,"serialized kind roundtrip");
        check(decoded[i].source_domain_id==envelopes[i].source_domain_id,"serialized source domain roundtrip");
        check(decoded[i].destination_domain_id==envelopes[i].destination_domain_id,"serialized destination domain roundtrip");
        check(decoded[i].source_rank==envelopes[i].source_rank,"serialized source rank roundtrip");
        check(decoded[i].destination_rank==envelopes[i].destination_rank,"serialized destination rank roundtrip");
        if(decoded[i].kind==PicTransportPayloadKind3D::particle_migration){
            check(decoded[i].particle_migration.particles.size()==envelopes[i].particle_migration.particles.size(),"serialized particle message size");
            for(std::size_t j=0;j<decoded[i].particle_migration.particles.size();++j){
                check(std::abs(decoded[i].particle_migration.particles[j].position_m.x-envelopes[i].particle_migration.particles[j].position_m.x)<1.0e-15,"serialized particle position");
                check(decoded[i].particle_migration.source_indices[j]==envelopes[i].particle_migration.source_indices[j],"serialized particle source index");
            }
        }else{
            check(decoded[i].scalar_guard.values.size()==envelopes[i].scalar_guard.values.size(),"serialized guard message size");
            if(!decoded[i].scalar_guard.values.empty()){
                check(std::abs(decoded[i].scalar_guard.values.front().value-envelopes[i].scalar_guard.values.front().value)<1.0e-15,"serialized guard first value");
            }
        }
    }
}

void serialized_exchange_plan_and_apply(){
    using namespace cfd::particle;
    auto fixture=make_fixture();
    const auto migration=pack_particle_migration_messages_3d(fixture.particles_by_domain,fixture.domains,fixture.config);
    const auto guard_messages=pack_scalar_guard_cell_messages_3d(fixture.values,fixture.domains,fixture.config,1U);
    const auto topology=make_pic_rank_topology_3d(fixture.domains,3U,1U,1U);
    const auto envelopes=make_pic_transport_envelopes_3d(migration.messages,guard_messages,topology);
    const auto plan=plan_serialized_pic_exchange_3d(envelopes,topology.rank_count);
    check(plan.rank_count==3U,"serialized plan rank count");
    std::size_t planned_messages=0U,planned_bytes=0U;
    for(std::size_t rank=0;rank<plan.rank_count;++rank){
        planned_messages+=plan.message_count_by_destination_rank[rank];planned_bytes+=plan.byte_count_by_destination_rank[rank];
        check(plan.envelopes_by_destination_rank[rank].size()==plan.message_count_by_destination_rank[rank],"serialized plan bucket size");
    }
    check(planned_messages==envelopes.size()&&planned_bytes>0U,"serialized plan payload totals");
    std::vector<PicSerializedEnvelope3D> delivered_serialized;
    for(const auto& bucket:plan.envelopes_by_destination_rank)delivered_serialized.insert(delivered_serialized.end(),bucket.begin(),bucket.end());
    const auto delivered=deserialize_pic_transport_envelopes_3d(delivered_serialized);
    std::vector<ParticleMigrationMessage3D> delivered_particles;std::vector<ScalarGuardCellMessage3D> delivered_guards;
    for(const auto& envelope:delivered){
        if(envelope.kind==PicTransportPayloadKind3D::particle_migration)delivered_particles.push_back(envelope.particle_migration);
        else delivered_guards.push_back(envelope.scalar_guard);
    }
    const auto serial_particles=apply_particle_migration_messages_3d(migration.retained_by_domain,migration.messages);
    const auto serialized_particles=apply_particle_migration_messages_3d(migration.retained_by_domain,delivered_particles);
    std::size_t serial_count=0U,serialized_count=0U;
    for(std::size_t domain=0;domain<serial_particles.size();++domain){serial_count+=serial_particles[domain].size();serialized_count+=serialized_particles[domain].size();}
    check(serial_count==serialized_count&&serialized_count==4U,"serialized migration apply parity");
    const auto serial_guarded=apply_scalar_guard_cell_messages_3d(fixture.values,fixture.domains,1U,-1.0,guard_messages);
    const auto serialized_guarded=apply_scalar_guard_cell_messages_3d(fixture.values,fixture.domains,1U,-1.0,delivered_guards);
    check(serial_guarded.size()==serialized_guarded.size(),"serialized guard block count");
    for(std::size_t domain=0;domain<serial_guarded.size();++domain){
        check(serial_guarded[domain].values.size()==serialized_guarded[domain].values.size(),"serialized guard block size");
        for(std::size_t i=0;i<serial_guarded[domain].values.size();++i){
            check(std::abs(serial_guarded[domain].values[i]-serialized_guarded[domain].values[i])<1.0e-15,"serialized guard apply parity");
        }
    }
}

void rejects_malformed_envelopes(){
    using namespace cfd::particle;
    bool caught=false;try{(void)deserialize_pic_transport_envelope_3d(PicSerializedEnvelope3D{});}catch(const std::exception&){caught=true;}check(caught,"empty serialized envelope rejected");
    PicSerializedEnvelope3D bad;bad.bytes={0,1,2,3,4};caught=false;try{(void)deserialize_pic_transport_envelope_3d(bad);}catch(const std::exception&){caught=true;}check(caught,"truncated serialized envelope rejected");
}
}

int main(){
    try{
        serialized_roundtrip();
        serialized_exchange_plan_and_apply();
        rejects_malformed_envelopes();
    }catch(const std::exception& error){
        std::cerr<<"v0.10.5 serialized transport test failed: "<<error.what()<<'\n';
        return 1;
    }
    std::cout<<"v0.10.5 serialized transport tests passed\n";
    return 0;
}
