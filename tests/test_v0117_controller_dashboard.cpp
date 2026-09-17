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
    return std::filesystem::temp_directory_path() / ("cfd_solvers_v0117_" + std::to_string(stamp));
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

    const std::vector<CampaignParameter> factors{{"angle", CampaignParameterKind::discrete, 0.0, 0.0, {"0", "4", "8"}}};
    const auto cases = generate_factorial_campaign(factors);
    const std::vector<CampaignTemplateFile> templates{{"DATA/setup.cfg", "angle={angle}\n"}};
    CampaignWriteOptions write; write.overwrite = true;
    const auto persisted = write_campaign_case_folders(root, cases, templates, write);
    auto registry = load_campaign_registry(root / "registry.tsv");

    SolverAdapterDescriptor adapter{"stub", "cfd-solve", "cfd-solve", {SolverCapability::residuals, SolverCapability::performance, SolverCapability::control}, {"stop", "extend", "checkpoint"}};
    const std::vector<CampaignServerDescriptor> servers{
        {"gpu-node", "10.0.0.42", "cfd", "/srv/cfd/campaigns", CampaignServerRuntime::docker, 2U, {"gpu", "linux"}, "cfd-solvers:dash", true},
        {"local-node", "localhost", "", "/tmp/cfd_dashboard", CampaignServerRuntime::native, 1U, {"cpu", "linux"}, "", true}
    };
    const auto campaign_plan = plan_multi_server_campaign(root, registry, adapter, servers);
    const auto exec_plan = plan_multiserver_execution(campaign_plan, servers);
    const auto supervision = plan_multiserver_supervision(exec_plan, servers);

    RemoteControllerConfig controller;
    controller.port = 9191U;
    controller.dashboard_title = "CFD campaign dashboard test";
    controller.dashboard_refresh_seconds = 7U;
    controller.control_actions = {"stop", "extend", "checkpoint"};
    const auto plan = plan_multiserver_controller(supervision, controller);

    require(persisted.cases_written == 3U, "campaign should write three cases");
    require(plan.routes.size() == 12U, "dashboard controller should expose production routes");
    require(plan.index_html.find("CFD campaign dashboard test") != std::string::npos, "dashboard HTML should include title");
    require(plan.index_html.find("data-refresh-seconds=\"7\"") != std::string::npos, "dashboard HTML should include refresh interval");
    require(plan.dashboard_js.find("/api/status") != std::string::npos, "dashboard JS should poll status");
    require(plan.dashboard_js.find("sendControl") != std::string::npos, "dashboard JS should include control action helper");
    require(plan.dashboard_css.find(".state.running") != std::string::npos, "dashboard CSS should style states");
    require(plan.events_ndjson.find("cfd.multiserver.controller.events.v1") != std::string::npos, "events snapshot should include schema");
    require(plan.events_ndjson.find("case0001") != std::string::npos, "events snapshot should include case events");
    require(plan.openapi_json.find("/api/events") != std::string::npos, "OpenAPI should include event snapshot route");
    require(plan.controller_script.find("/api/events") != std::string::npos, "controller script should serve event snapshots");
    require(plan.controller_script.find("dashboard.css") != std::string::npos, "controller script should serve stylesheet asset");

    const auto written = write_multiserver_controller_files(root, plan);
    require(written.routes_written == 12U, "routes count mismatch");
    require(written.static_assets_written == 7U, "generated asset count mismatch");
    require(std::filesystem::exists(written.index_path), "index.html should be written");
    require(std::filesystem::exists(written.dashboard_js_path), "dashboard.js should be written");
    require(std::filesystem::exists(written.dashboard_css_path), "dashboard.css should be written");
    require(std::filesystem::exists(written.events_path), "events.ndjson should be written");
    require(read_file(written.index_path).find("dashboard.js") != std::string::npos, "written HTML should load JS asset");
    require(read_file(written.dashboard_js_path).find("refreshSeconds = 7") != std::string::npos, "written JS should preserve refresh interval");
    require(read_file(written.events_path).find("\"event\":\"case\"") != std::string::npos, "written events should include case records");
    require(read_file(written.readme_path).find("generated live dashboard") != std::string::npos, "controller README should mention dashboard");

    RemoteControllerConfig no_events;
    no_events.expose_event_snapshot = false;
    const auto no_event_plan = plan_multiserver_controller(supervision, no_events);
    require(no_event_plan.routes.size() == 11U, "controller without snapshot route should retain production routes");
    require(no_event_plan.openapi_json.find("\"/api/events\":") == std::string::npos, "snapshot route should be omitted when disabled");

    std::filesystem::remove_all(root);
    std::cout << "v0.11.7 controller dashboard tests passed\n";
    return 0;
}
