#include "cfd/workflow/campaign.hpp"
#include "cfd/workflow/campaign_control.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);} 
std::filesystem::path unique_temp_root(){
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path()/("cfd_solvers_v0112_"+std::to_string(stamp));
}
}

int main(){
    using namespace cfd::workflow;
    const auto root=unique_temp_root();
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"case_param",CampaignParameterKind::discrete,0.0,0.0,{"A","B","C"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","case={case_param}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    require(persisted.cases_written==3U,"campaign writer should persist three cases");

    SolverAdapterDescriptor adapter{"stub","/bin/echo","stub",{SolverCapability::control,SolverCapability::residuals,SolverCapability::performance},{"stop","extend","checkpoint","flush"}};
    const CampaignControlRequest stop{"case0001",CampaignControlAction::stop,"","operator requested stop"};
    const auto directive=write_campaign_control_directive(root,adapter,stop);
    require(std::filesystem::exists(directive.path),"stop directive should be written");
    const auto discovered=discover_campaign_control_directives(root/"case0001");
    require(discovered.size()>=2U,"stop and latest directives should be discoverable");
    require(discovered.front().case_id=="case0001","discovered directive keeps case id");
    require(campaign_control_action_name(discovered.front().action)=="stop","discovered directive keeps action");

    bool rejected=false;
    try{const auto checkpoint=write_campaign_control_directive(root,adapter,{"case0001",CampaignControlAction::checkpoint,"",""},{".cfd_control",false});(void)checkpoint;}
    catch(const std::exception&){rejected=true;}
    require(!rejected,"different checkpoint directive may coexist with stop directive");

    SolverAdapterDescriptor no_control{"plain","/bin/echo","plain",{},{} };
    bool unsupported=false;
    try{const auto rejected_directive=write_campaign_control_directive(root,no_control,stop);(void)rejected_directive;}catch(const std::invalid_argument&){unsupported=true;}
    require(unsupported,"adapter without control capability should reject directives");

    auto registry=load_campaign_registry(root/"registry.tsv");
    const std::string slurm_text=
        "JOBID CASE STATE ELAPSED\n"
        "101 case0001 RUNNING 12\n"
        "102 case0002 COMPLETED 20\n"
        "103 case0003 CANCELLED 3\n";
    const auto jobs=parse_slurm_queue_table(slurm_text);
    require(jobs.size()==3U,"slurm parser should read three jobs");
    require(jobs[0].state==SchedulerJobState::running,"slurm RUNNING maps to running");
    const auto sched_summary=apply_scheduler_records_to_registry(registry,jobs);
    require(sched_summary.running==1U,"scheduler summary should have one running case");
    require(sched_summary.done==1U,"scheduler summary should have one done case");
    require(sched_summary.stopped==1U,"scheduler summary should have one stopped case");

    {
        std::ofstream out(root/"case0001"/"solver.stdout.log");
        out<<"ITER 7 RESIDUAL 2e-5 FIELD flow\nperf.wall_time 0.4 s\nDONE\n";
    }
    {
        std::ofstream out(root/"case0002"/"history.dat");
        out<<"8 1e-6\n";
    }
    const auto refreshed=refresh_campaign_status_from_outputs(root,registry);
    require(refreshed.done>=2U,"output refresh should mark completed logged cases as done");
    auto case1=registry[0];
    for(const auto& entry:registry)if(entry.case_id=="case0001")case1=entry;
    require(case1.status==CampaignCaseStatus::done,"case0001 should be done after log refresh");
    require(case1.iterations==7U,"case0001 should keep parsed iteration");
    require(case1.objective<1e-4,"case0001 should keep parsed residual objective");

    save_campaign_registry(root/"registry.tsv",registry);
    const auto persisted_summary=refresh_campaign_registry_from_outputs(root);
    require(persisted_summary.done>=2U,"persisted refresh should save summary");

    std::filesystem::remove_all(root);
    return 0;
}
