#include "cfd/workflow/campaign.hpp"
#include "cfd/workflow/deploy.hpp"

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
    return std::filesystem::temp_directory_path()/("cfd_solvers_v0113_"+std::to_string(stamp));
}
std::string read_file(const std::filesystem::path& path){std::ifstream in(path);return std::string((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());} 
}

int main(){
    using namespace cfd::workflow;
    const auto root=unique_temp_root();
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6","0.8","1.0"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    require(persisted.cases_written==5U,"campaign should persist five cases");
    auto registry=load_campaign_registry(root/"registry.tsv");
    registry[3].status=CampaignCaseStatus::done;
    save_campaign_registry(root/"registry.tsv",registry);

    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.10.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:test",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true},
        {"offline","10.10.0.13","cfd","/srv/offline",CampaignServerRuntime::docker,4U,{"gpu"},"cfd-solvers:test",false}
    };
    MultiServerSchedulingOptions opts;opts.pending_only=true;opts.require_online=true;
    const auto plan=plan_multi_server_campaign(root,registry,adapter,servers,opts);
    require(plan.active_servers==2U,"offline server should be filtered");
    require(plan.assignments.size()==4U,"only pending cases should be assigned");
    require(plan.unassigned_cases==0U,"server capacity should fit pending cases");
    require(plan.commands.size()==12U,"each assignment should create mkdir, sync and run commands");
    require(plan.assignments[0].server_name=="cpu-a","tie-breaking should prefer lexical server when equally loaded");
    bool saw_docker=false,saw_ssh=false;
    for(const auto& cmd:plan.commands){
        if(cmd.display_command.find("docker run")!=std::string::npos)saw_docker=true;
        if(cmd.uses_ssh)saw_ssh=true;
    }
    require(saw_docker,"docker server should produce docker run command");
    require(saw_ssh,"remote server should produce ssh/rsync commands");

    const auto written=write_multiserver_plan_files(root,plan);
    require(written.assignments_written==4U,"assignment file should contain four rows");
    require(written.commands_written==12U,"command script should contain commands");
    require(read_file(written.assignments_path).find("case0001")!=std::string::npos,"assignment file should include case id");
    require(read_file(written.commands_path).find("set -eu")!=std::string::npos,"command script should be shell-safe scaffold");

    DockerDeployConfig deploy;deploy.image_name="cfd-solvers:test";deploy.worker_replicas=3U;deploy.worker_cpus=2U;deploy.worker_memory="2g";deploy.enable_gpu_runtime=false;deploy.environment={{"CFD_LOG_LEVEL","info"}};
    const auto files=generate_docker_deploy_files(deploy);
    require(files.size()>=6U,"docker deploy generator should create multiple files");
    const auto deploy_result=write_docker_deploy_system(root/"deploy_stack",deploy,true);
    require(deploy_result.files_written==files.size(),"docker deploy writer should persist every generated file");
    const auto compose=read_file(root/"deploy_stack"/"compose.yaml");
    require(compose.find("campaign-manager")!=std::string::npos,"compose should include manager service");
    require(compose.find("replicas: 3")!=std::string::npos,"compose should include worker replicas");
    require(compose.find("CFD_LOG_LEVEL: info")!=std::string::npos,"compose should include environment values");
    require(std::filesystem::exists(root/"deploy_stack"/"deploy"/"servers.example.tsv"),"server example should be written");

    std::filesystem::remove_all(root);
    return 0;
}
