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
    return std::filesystem::temp_directory_path() / ("cfd_solvers_v0115_" + std::to_string(stamp));
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

    SolverAdapterDescriptor adapter{"stub", "cfd-solve", "cfd-solve", {SolverCapability::residuals, SolverCapability::performance, SolverCapability::control}, {"stop", "checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-a", "10.10.0.11", "cfd", "/srv/cfd/campaigns", CampaignServerRuntime::docker, 2U, {"gpu", "linux"}, "cfd-solvers:test", true},
        {"cpu-a", "localhost", "", "/tmp/cfd_remote", CampaignServerRuntime::native, 2U, {"cpu", "linux"}, "", true}
    };
    const auto campaign_plan = plan_multi_server_campaign(root, registry, adapter, servers);
    const auto exec_plan = plan_multiserver_execution(campaign_plan, servers);
    RemoteSupervisionOptions options;
    options.max_attempts = 4U;
    options.retry_delay_seconds = 1U;
    options.log_tail_lines = 40U;
    const auto plan = plan_multiserver_supervision(exec_plan, servers, options);

    require(persisted.cases_written == 4U, "campaign should write four cases");
    require(plan.dashboard_cases.size() == 4U, "dashboard should include one row per job");
    require(plan.retry_launches.size() == 4U, "retry plan should include each job");
    require(plan.log_tails.size() == 4U, "log tail plan should include each job");
    require(plan.access_probes.size() >= 3U, "access probes should include local/ssh/docker/writable checks");
    require(plan.retry_launches.front().command.display_command.find("attempt=1") != std::string::npos, "retry command should include attempt loop");
    require(plan.log_tails.front().stdout_tail_command.display_command.find("tail -n 40") != std::string::npos, "stdout tail command should respect tail length");

    const auto redacted = redact_sensitive_command_display("docker run -e TOKEN=abcd -e PASSWORD=hunter2 image", options.sensitive_markers);
    require(redacted.find("abcd") == std::string::npos && redacted.find("hunter2") == std::string::npos, "redaction should hide secret values");
    require(redacted.find("TOKEN=<redacted>") != std::string::npos, "redaction should keep the marker shape");

    const auto json = build_multiserver_dashboard_json(plan.dashboard_cases);
    require(json.find("cfd.multiserver.dashboard.v1") != std::string::npos, "dashboard JSON should include schema");
    require(json.find("case0001") != std::string::npos, "dashboard JSON should include case IDs");

    const auto written = write_multiserver_supervision_files(root, plan);
    require(written.access_checks_written == plan.access_probes.size(), "access script count mismatch");
    require(written.retry_launches_written == 4U, "retry script count mismatch");
    require(written.log_tail_commands_written == 8U, "tail script should include stdout and stderr per job");
    require(written.dashboard_cases_written == 4U, "dashboard case count mismatch");
    require(read_file(written.access_script_path).find("Docker daemon") != std::string::npos, "access script should include docker probe");
    require(read_file(written.retry_launch_script_path).find("retry launch") != std::string::npos, "retry script should include retry comments");
    require(read_file(written.tail_logs_script_path).find("tail stderr") != std::string::npos, "tail script should include stderr commands");
    require(read_file(written.dashboard_json_path).find("stdout") != std::string::npos, "dashboard JSON file should include log paths");
    require(read_file(written.summary_path).find("case_id\tserver\tstate") != std::string::npos, "summary TSV should include header");

    std::filesystem::remove_all(root);
    std::cout << "v0.11.5 multi-server supervision tests passed\n";
    return 0;
}
