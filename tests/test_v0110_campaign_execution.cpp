#include "cfd/workflow/campaign.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

namespace {
void check(bool value,const char* message){if(!value)throw std::runtime_error(message);} 

std::string read_text(const std::filesystem::path& path){
    std::ifstream in(path);if(!in)throw std::runtime_error("failed to read test file");
    return std::string((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
}

void persistent_campaign_folders_and_registry(){
    using namespace cfd::workflow;
    const auto root=std::filesystem::temp_directory_path()/"cfd_v0110_campaign_execution";
    std::filesystem::remove_all(root);
    const std::vector<CampaignParameter> factors{{"mach",CampaignParameterKind::discrete,0.0,0.0,{"0.2","0.4"}},{"aoa",CampaignParameterKind::discrete,0.0,0.0,{"0","2"}}};
    const auto cases=generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg","mach={mach}\n<!-- IF aoa=2 -->\nangle={aoa}\n<!-- ENDIF -->"},{"run.sh","solver --mach {mach}\n"}};
    CampaignWriteOptions options;options.overwrite=true;
    const auto persisted=write_campaign_case_folders(root,cases,templates,options);
    check(persisted.cases_written==4U&&persisted.files_written==12U&&persisted.registry_entries==4U,"campaign folder writer reports case/file/registry counts");
    check(std::filesystem::exists(root/"case0004"/"DATA"/"setup.cfg"),"rendered case setup file exists");
    const auto rendered=read_text(root/"case0004"/"DATA"/"setup.cfg");
    check(rendered.find("mach=0.4")!=std::string::npos&&rendered.find("angle=2")!=std::string::npos,"campaign template renders placeholders and IF blocks");
    auto registry=load_campaign_registry(root/"registry.tsv");
    check(registry.size()==4U&&summarize_campaign_registry(registry).pending==4U,"campaign registry roundtrips pending rows");
    registry[1].status=CampaignCaseStatus::done;registry[1].objective=0.25;registry[1].iterations=18U;registry[1].message="completed\twith note";
    save_campaign_registry(root/"registry.tsv",registry);
    const auto loaded=load_campaign_registry(root/"registry.tsv");
    const auto summary=summarize_campaign_registry(loaded);
    check(summary.done==1U&&summary.best_case_id=="case0002"&&loaded[1].message.find("with note")!=std::string::npos,"campaign registry preserves status, objective and escaped messages");
    std::filesystem::remove_all(root);
}

void solver_command_plans_and_parsers(){
    using namespace cfd::workflow;
    SolverAdapterDescriptor adapter{"stub","stub-solver","stub-solver",{SolverCapability::residuals,SolverCapability::performance,SolverCapability::control},{"stop"}};
    SolverRunRequest native{adapter,"case0001",SolverRuntime::native,4U,2U,{"--steady"},{},{},{}};
    const auto native_plan=build_solver_command_plan(native);
    check(native_plan.argv.size()>=8U&&native_plan.argv[0]=="stub-solver"&&native_plan.display_command.find("--nprocs 4")!=std::string::npos,"native solver command plan is deterministic");
    SolverRunRequest docker{adapter,"/tmp/case0001",SolverRuntime::docker,2U,1U,{"--monitor"},"solver-image:latest",{}, {}};
    const auto docker_plan=build_solver_command_plan(docker);
    check(docker_plan.argv[0]=="docker"&&docker_plan.display_command.find("solver-image:latest")!=std::string::npos,"docker solver command plan wraps adapter binary");
    SolverRunRequest slurm{adapter,"case0001",SolverRuntime::slurm,8U,3U,{"--batch"},{},{},"debug"};
    const auto slurm_plan=build_solver_command_plan(slurm);
    check(slurm_plan.requires_scheduler&&slurm_plan.argv[0]=="sbatch"&&slurm_plan.display_command.find("--partition debug")!=std::string::npos,"slurm solver command plan records scheduler wrapper");
    const std::string log="iter=1 residual=1e-2 equation=flow\niter=2 residual=2.5e-4 equation=flow\nperf.wall_s=12.5 s\nCELLS_PER_S 4000 cell/s\nDONE\n";
    const auto residuals=parse_residual_history(log);
    const auto metrics=parse_performance_metrics(log);
    check(detect_solver_outcome_from_log(log)==CampaignCaseStatus::done,"solver log outcome detection sees completion");
    check(residuals.size()==2U&&residuals.back().iteration==2U&&residuals.back().residual<1.0e-3,"residual parser extracts iteration history");
    check(metrics.size()==2U&&metrics[0].key=="wall_s"&&std::abs(metrics[0].value-12.5)<1.0e-12,"performance parser extracts metrics");
    check(detect_solver_outcome_from_log("FATAL residual diverged")==CampaignCaseStatus::failed,"solver log outcome detection prioritizes failure");
}

void optimization_loop_improves_campaign_objective(){
    using namespace cfd::workflow;
    auto objective=[](std::span<const double> x){return (x[0]-1.5)*(x[0]-1.5)+0.5*(x[1]+0.25)*(x[1]+0.25);};
    CampaignOptimizationConfig config;config.max_iterations=80U;config.step_size=0.35;config.tolerance=1.0e-7;
    const std::vector<double> initial{-2.0,2.0};
    const auto result=run_gradient_descent_campaign(objective,initial,config);
    check(result.objective<objective(initial)*1.0e-3,"campaign gradient loop reduces objective by orders of magnitude");
    check(result.iterations>0U&&!result.parameters.empty()&&!result.last_gradient.empty(),"campaign gradient loop reports diagnostics");
}
}

int main(){
    persistent_campaign_folders_and_registry();
    solver_command_plans_and_parsers();
    optimization_loop_improves_campaign_objective();
    std::cout<<"v0.11.0 campaign execution tests passed\n";
    return 0;
}
