#include "cfd/workflow/campaign.hpp"
#include "cfd/workflow/deploy.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);} 
std::filesystem::path unique_temp_root(){
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path()/("cfd_solvers_v0114_"+std::to_string(stamp));
}
std::string read_file(const std::filesystem::path& path){std::ifstream in(path);return std::string((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());} 
}

int main(){
    using namespace cfd::workflow;
    const auto root=unique_temp_root();
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4","0.6","0.8"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    require(persisted.cases_written==4U,"campaign should persist four cases");
    auto registry=load_campaign_registry(root/"registry.tsv");

    SolverAdapterDescriptor adapter{"stub","cfd-solve","cfd-solve",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop","checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a","10.10.0.11","cfd","/srv/cfd/campaigns",CampaignServerRuntime::docker,2U,{"gpu","linux"},"cfd-solvers:test",true},
        {"cpu-a","localhost","","/tmp/cfd_remote",CampaignServerRuntime::native,2U,{"cpu","linux"},"",true}
    };
    const auto campaign_plan=plan_multi_server_campaign(root,registry,adapter,servers);
    require(campaign_plan.assignments.size()==4U,"all four cases should be assigned");
    const auto exec_plan=plan_multiserver_execution(campaign_plan,servers);
    require(exec_plan.jobs.size()==4U,"execution plan should create one supervised job per assignment");
    require(exec_plan.health_checks.size()==6U,"three health checks per active server should be generated");
    require(exec_plan.docker_jobs==2U&&exec_plan.native_jobs==2U,"execution plan should preserve docker/native split");
    require(exec_plan.remote_jobs==2U,"remote server jobs should be counted");
    require(exec_plan.jobs.front().launch_command.display_command.find("nohup")!=std::string::npos,"launch command should be supervised");
    require(exec_plan.jobs.front().status_command.display_command.find("printf")!=std::string::npos,"status command should emit a parseable row");
    require(exec_plan.jobs.front().cancel_command.display_command.find("kill")!=std::string::npos,"cancel command should target pid file");

    const auto written=write_multiserver_execution_files(root,exec_plan);
    require(written.health_checks_written==6U,"health script should include six checks");
    require(written.jobs_written==4U,"jobs table should include four jobs");
    require(read_file(written.launch_script_path).find("supervised launch")!=std::string::npos,"launch script should include heading");
    require(read_file(written.status_script_path).find("status case0001")!=std::string::npos,"status script should include case probes");
    require(read_file(written.jobs_path).find("solver.stdout.log")!=std::string::npos,"jobs table should include log paths");

    const auto statuses=parse_remote_job_status_table("case_id\tserver\tstate\texit_code\ncase0001\tcpu-a\trunning\t0\ncase0002\tgpu-a\tdone\t0\ncase0003\tgpu-a\tfailed\t2\ncase0004\tcpu-a\tcancelled\t0\n");
    require(statuses.size()==4U&&statuses[1].state==RemoteJobState::done,"remote status parser should decode table rows");
    const auto summary=apply_remote_job_status_to_registry(registry,statuses);
    require(summary.running==1U&&summary.done==1U&&summary.failed==1U&&summary.stopped==1U,"remote job statuses should update registry states");

    DockerDeployConfig deploy;deploy.image_name="cfd-solvers:test";deploy.worker_replicas=2U;deploy.environment={{"CFD_WORKER_COMMAND","cfd-solve particle-campaign-local-runner"}};
    const auto deploy_result=write_docker_deploy_system(root/"deploy_stack",deploy,true);
    require(deploy_result.files_written>=8U,"deploy writer should include healthcheck and entrypoint scripts");
    const auto compose=read_file(root/"deploy_stack"/"compose.yaml");
    require(compose.find("healthcheck")!=std::string::npos,"compose should contain healthchecks");
    require(std::filesystem::exists(root/"deploy_stack"/"deploy"/"worker-entrypoint.sh"),"worker entrypoint should be written");

    std::filesystem::remove_all(root);
    std::cout<<"v0.11.4 multi-server execution tests passed\n";
    return 0;
}
