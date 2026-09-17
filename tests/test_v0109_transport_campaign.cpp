#include "cfd/particle/transport.hpp"
#include "cfd/multiphysics/bioheat.hpp"
#include "cfd/workflow/campaign.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 

cfd::particle::TransportWorld make_transport_world(){
    using namespace cfd::particle;
    TransportWorld world;
    world.particles.push_back({"electron",-1.602176634e-19,9.1093837e-31});
    world.materials.push_back({"argon",1.6,2.0});
    world.sensitive_detectors.push_back({"dose-slab",0.0,false});
    world.regions.push_back({"gas-detector",{{0.0,0.0,0.0},{1.0,1.0,1.0}},0U,0U});
    return world;
}

void stochastic_physics_hits_and_dose(){
    using namespace cfd::particle;
    using namespace cfd::multiphysics;
    auto world=make_transport_world();
    TransportPhysicsList list;list.name="em-standard-lite";list.production_cut_energy_ev=1.0;
    list.processes.push_back({"dEdx",TransportProcessKind::continuous_energy_loss,0U,0U,250.0,0.0,0.0,0.0,0.0});
    list.processes.push_back({"ionization",TransportProcessKind::discrete_interaction,0U,0U,0.0,0.02,0.1,0.15,0.0});
    world=apply_physics_list(std::move(world),list);
    TransportTrack track;track.position_m={0.1,0.5,0.5};track.direction={1.0,0.0,0.0};track.kinetic_energy_ev=1500.0;
    TransportConfig config;config.max_step_m=0.2;config.energy_cut_ev=0.5;
    TransportRandom rng; rng.state=12345U;
    const auto result=transport_track_stochastic(world,track,config,rng,16U);
    check(result.scoring.steps>0U,"stochastic transport records steps");
    check(result.scoring.total_energy_deposit_ev>0.0,"stochastic transport deposits energy");
    check(result.scoring.secondaries>0U,"stochastic discrete process creates secondaries");
    const auto hits=collect_transport_hits(world,result.steps);
    check(!hits.hits.empty(),"sensitive detector receives hits");
    check(hits.total_energy_deposit_ev>0.0,"hit collection accumulates energy");
    auto grid=make_transport_dose_grid(4U,4U,4U,{{0.0,0.0,0.0},{1.0,1.0,1.0}},1000.0);
    score_hits_to_dose_grid(hits,grid);
    const auto sar=dose_rate_to_sar_w_per_kg(grid,0.5);
    check(std::any_of(sar.begin(),sar.end(),[](double v){return v>0.0;}),"dose grid converts to SAR/dose rate");
    const auto pennes_sar=project_dose_grid_to_pennes2d_sar(grid,0.5,6U,6U);
    PennesBioheat2DConfig bio;bio.nx=6U;bio.ny=6U;bio.dt_s=0.01;bio.blood_temperature_k=310.0;PennesBioheat2D solver(bio);solver.initialize(310.0);solver.set_sar(pennes_sar);solver.step();
    check(solver.linear_result().converged,"transport dose projection drives Pennes solve");
}

void bvh_matches_linear_lookup(){
    using namespace cfd::particle;
    TransportWorld world;world.particles.push_back({"gamma",0.0,1.0});world.materials.push_back({"mat",1.0,0.0});
    for(std::size_t i=0;i<8U;++i)world.regions.push_back({"cell",{{static_cast<double>(i),0.0,0.0},{static_cast<double>(i+1U),1.0,1.0}},0U});
    const auto bvh=build_region_bvh(world,2U);
    check(!bvh.nodes.empty(),"BVH has nodes");
    for(std::size_t i=0;i<8U;++i){const cfd::particle::Vec3 p{static_cast<double>(i)+0.25,0.5,0.5};check(locate_region(world,p)==locate_region_bvh(world,bvh,p),"BVH lookup matches linear lookup");}
}

void campaign_doe_adapter_and_gradients(){
    using namespace cfd::workflow;
    const std::vector<CampaignParameter> factorial{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4"}},{"aoa",CampaignParameterKind::discrete,0.0,0.0,{"0","2"}}};
    const auto cases=generate_factorial_campaign(factorial);check(cases.size()==4U,"factorial DOE crosses levels");
    const auto rendered=render_campaign_template("MACH={mach}\n<!-- IF aoa=2 -->HIGH_AOA=1\n<!-- ENDIF -->",cases[1]);
    check(rendered.find("MACH=")!=std::string::npos,"template placeholders render");
    const std::vector<CampaignParameter> lhs{{"re",CampaignParameterKind::continuous,1.0e5,2.0e5,{}},{"model",CampaignParameterKind::discrete,0.0,0.0,{"sa","sst"}}};
    auto latin=generate_latin_hypercube_campaign(lhs,5U,{42U});check(latin.size()==5U&&latin[0].values.contains("re")&&latin[0].values.contains("model"),"LHS DOE carries continuous and discrete values");
    SolverAdapterDescriptor adapter{"su2", "SU2_CFD", "SU2_CFD", {SolverCapability::residuals,SolverCapability::restart,SolverCapability::performance}, {"stop","checkpoint"}};validate_solver_adapter_descriptor(adapter);check(solver_has_capability(adapter,SolverCapability::restart),"adapter capabilities are derived from descriptor");
    auto objective=[](std::span<const double> x){return (x[0]-3.0)*(x[0]-3.0)+(x[1]+1.0)*(x[1]+1.0);};
    std::vector<double> x{2.0,-0.5};const auto grad=finite_difference_gradient(objective,x);check(std::abs(grad[0]+2.0)<1.0e-5&&std::abs(grad[1]-1.0)<1.0e-5,"finite difference gradient approximates analytic gradient");
    const auto next=gradient_descent_update(x,grad,0.25);check(objective(next)<objective(x),"gradient update improves objective");
    std::vector<CampaignRegistryEntry> registry;update_case_status(registry,{cases[0].case_id,CampaignCaseStatus::done,1.5,12U,"ok"});update_case_status(registry,{cases[1].case_id,CampaignCaseStatus::failed,0.0,3U,"diverged"});const auto summary=summarize_campaign_registry(registry);check(summary.done==1U&&summary.failed==1U&&summary.best_case_id==cases[0].case_id,"registry summary tracks outcomes");
}
}

int main(){
    stochastic_physics_hits_and_dose();
    bvh_matches_linear_lookup();
    campaign_doe_adapter_and_gradients();
    std::cout<<"v0.10.9 Geant4/SU2/csauto transport-campaign tests passed\n";
    return 0;
}
