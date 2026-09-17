#include "cfd/workflow/campaign.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(_WIN32)
#include <sys/stat.h>
#endif

namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);} 
std::filesystem::path unique_temp_root(){
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path()/("cfd_solvers_v0111_"+std::to_string(stamp));
}
void make_executable(const std::filesystem::path& path){
#if !defined(_WIN32)
    ::chmod(path.c_str(),0755);
#else
    (void)path;
#endif
}
}

int main(){
    using namespace cfd::workflow;
    const auto root=unique_temp_root();
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n"}};
    CampaignWriteOptions write;write.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,write);
    require(persisted.cases_written==2U,"expected two campaign cases");

    const auto script=root/"stub_solver.sh";
    {
        std::ofstream out(script);
        out << "#!/bin/sh\n"
            << "echo 'iter=1 residual=1e-2 equation=flow'\n"
            << "echo 'iter=2 residual=5e-4 equation=flow'\n"
            << "echo 'perf.wall_time 0.02 s'\n"
            << "echo 'DONE'\n"
            << "echo 'ITER 3 RESIDUAL 1e-5 FIELD flow' > residual_history.log\n"
            << "echo 'perf.cells_per_s 1000 cells/s' > timing.log\n"
            << "exit 0\n";
    }
    make_executable(script);

    SolverAdapterDescriptor adapter{"stub",script.string(),"stub",{SolverCapability::residuals,SolverCapability::performance},{"stop"}};
    SolverRunRequest request{adapter,root/"case0001",SolverRuntime::native,1U,1U,{"--dry-run"},{},{},{}};
    const auto doctor=doctor_solver_runtime(request);
    require(doctor.ok,"doctor should accept executable stub and case directory");
    SolverRunRequest missing= request;missing.adapter.native_binary="definitely-not-a-real-cfd-solver";
    require(!doctor_solver_runtime(missing).ok,"doctor should reject missing native executable");

    const auto one=run_local_solver_case(request);
    require(one.status==CampaignCaseStatus::done,"single local run should finish as done");
    require(one.exit_code==0,"single local run should have zero exit code");
    require(one.residuals.size()>=3U,"single local run should parse stdout and file residuals");
    require(one.performance.size()>=2U,"single local run should parse stdout and file performance metrics");
    require(!one.discovered_outputs.residual_files.empty(),"residual file should be discovered");
    require(!one.discovered_outputs.performance_files.empty(),"performance file should be discovered");

    LocalCampaignRunOptions options;options.adapter=adapter;options.runtime=SolverRuntime::native;options.mpi_ranks=1U;options.threads=1U;options.extra_args={"--from-campaign"};
    const auto run=run_local_campaign(root,options);
    require(run.launched==2U,"local campaign runner should launch both pending cases");
    require(run.summary.done==2U,"both cases should be done after local campaign run");
    const auto registry=load_campaign_registry(root/"registry.tsv");
    require(registry.size()==2U,"registry should keep both entries");
    require(registry[0].status==CampaignCaseStatus::done&&registry[1].status==CampaignCaseStatus::done,"registry entries should be done");
    require(registry[0].iterations>=3U,"registry should store parsed iterations");
    require(std::filesystem::exists(root/"case0002"/"solver.stdout.log"),"launcher stdout log should be written");

    std::filesystem::remove_all(root);
    return 0;
}
