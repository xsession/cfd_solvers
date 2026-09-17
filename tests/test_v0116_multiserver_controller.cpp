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
void require(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }
std::filesystem::path unique_temp_root() {
    const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path() / ("cfd_solvers_v0116_" + std::to_string(stamp));
}
std::string read_file(const std::filesystem::path& path) {
    std::ifstream in(path);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
}
}

int main() {
    using namespace cfd::workflow;
    const auto root = unique_temp_root();
    std::filesystem::remove_all(root);

    const std::vector<CampaignParameter> factors{{"mach", CampaignParameterKind::discrete, 0.0, 0.0, {"0.2", "0.4", "0.6", "0.8"}}};
    const auto cases = generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg", "mach={mach}\n"}};
    CampaignWriteOptions write; write.overwrite = true;
    const auto persisted = write_campaign_case_folders(root, cases, templates, write);
    auto registry = load_campaign_registry(root / "registry.tsv");

    SolverAdapterDescriptor adapter{"stub", "cfd-solve", "cfd-solve", {SolverCapability::residuals, SolverCapability::performance, SolverCapability::control}, {"stop", "checkpoint", "flush"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a", "10.10.0.11", "cfd", "/srv/cfd/campaigns", CampaignServerRuntime::docker, 2U, {"gpu", "linux"}, "cfd-solvers:test", true},
        {"cpu-a", "localhost", "", "/tmp/cfd_remote", CampaignServerRuntime::native, 2U, {"cpu", "linux"}, "", true}
    };
    const auto campaign_plan = plan_multi_server_campaign(root, registry, adapter, servers);
    const auto exec_plan = plan_multiserver_execution(campaign_plan, servers);
    RemoteSupervisionOptions supervision_options;
    supervision_options.log_tail_lines = 24U;
    const auto supervision = plan_multiserver_supervision(exec_plan, servers, supervision_options);

    RemoteControllerConfig controller;
    controller.port = 9090U;
    controller.log_tail_lines = 24U;
    controller.control_actions = {"stop", "checkpoint", "flush"};
    const auto plan = plan_multiserver_controller(supervision, controller);

    require(persisted.cases_written == 4U, "campaign should write four cases");
    require(plan.cases.size() == 4U, "controller should expose all dashboard cases");
    require(plan.routes.size() == 12U, "controller should expose protected production routes");
    require(plan.openapi_json.find("/api/cases/{case_id}/control/{action}") != std::string::npos, "OpenAPI JSON should include control route");
    require(plan.openapi_json.find("9090") != std::string::npos, "OpenAPI JSON should include configured port");
    require(plan.status_json.find("cfd.multiserver.controller.status.v1") != std::string::npos, "status JSON should include schema");
    require(plan.status_json.find("case0001") != std::string::npos, "status JSON should include case IDs");
    require(plan.controller_script.find("CFD_CONTROLLER_TOKEN") != std::string::npos, "controller should support token env var");
    require(plan.controller_script.find("ThreadingHTTPServer") != std::string::npos, "controller should use stdlib HTTP server");
    require(plan.controller_script.find("latest.directive") != std::string::npos, "controller should write control directives");

    const auto written = write_multiserver_controller_files(root, plan);
    require(written.routes_written == 12U, "routes table count mismatch");
    require(written.cases_written == 4U, "case count mismatch");
    require(std::filesystem::exists(written.controller_script_path), "controller script should be written");
    require(std::filesystem::exists(written.start_script_path), "start script should be written");
    require(read_file(written.routes_path).find("requires_token") != std::string::npos, "routes table should have security column");
    require(read_file(written.env_path).find("CFD_CONTROLLER_PORT=9090") != std::string::npos, "env example should include port");
    require(read_file(written.readme_path).find("Multi-server controller") != std::string::npos, "controller readme should be generated");

    RemoteControllerConfig read_only;
    read_only.allow_mutating_actions = false;
    const auto read_only_plan = plan_multiserver_controller(supervision, read_only);
    require(read_only_plan.routes.size() == 11U, "read-only controller should omit mutation route");
    require(read_only_plan.openapi_json.find("/api/cases/{case_id}/control/{action}") == std::string::npos, "read-only OpenAPI should not expose control route");

    std::filesystem::remove_all(root);
    std::cout << "v0.11.6 multi-server controller tests passed\n";
    return 0;
}
