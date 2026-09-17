#include "cfd/workflow/deploy.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void require(bool value, const std::string& message) { if (!value) throw std::runtime_error(message); }
std::string read_file(const std::filesystem::path& path) { std::ifstream in(path); return {std::istreambuf_iterator<char>(in), {}}; }
}

int main() {
    using namespace cfd::workflow;
    MultiServerSupervisionPlan supervision;
    supervision.dashboard_cases.push_back({"case0001", "node-a", "running", 0, "case0001/stdout.log", "case0001/stderr.log"});
    RemoteControllerConfig config;
    const auto plan = plan_multiserver_controller(supervision, config);

    require(plan.routes.size() == 12U, "production controller route count mismatch");
    require(plan.openapi_json.find("/api/events/history") != std::string::npos, "history route missing");
    require(plan.openapi_json.find("/api/probes") != std::string::npos, "probe route missing");
    require(plan.controller_script.find("ssl.PROTOCOL_TLS_SERVER") != std::string::npos, "TLS server configuration missing");
    require(plan.controller_script.find("hmac.compare_digest") != std::string::npos, "constant-time token comparison missing");
    require(plan.controller_script.find("required_role='viewer'") != std::string::npos, "role authorization missing");
    require(plan.controller_script.find("events-history.ndjson") != std::string::npos, "durable history missing");
    require(plan.controller_script.find("subprocess.run") != std::string::npos, "supervised live probes missing");
    require(plan.systemd_unit.find("ProtectSystem=strict") != std::string::npos, "systemd hardening missing");
    require(plan.reverse_proxy_config.find("Strict-Transport-Security") != std::string::npos, "reverse proxy security headers missing");
    require(plan.tokens_example_json.find("replace-with-sha256") != std::string::npos, "safe token template missing");

    const auto root = std::filesystem::temp_directory_path() / ("cfd_v0120_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    const auto written = write_multiserver_controller_files(root, plan);
    require(std::filesystem::exists(written.tokens_example_path), "token template not written");
    require(std::filesystem::exists(written.systemd_path), "systemd unit not written");
    require(std::filesystem::exists(written.reverse_proxy_path), "proxy example not written");
    require(read_file(written.env_path).find("change-me") == std::string::npos, "env template must not contain a default credential");
    const std::string compile = "python3 -m py_compile \"" + written.controller_script_path.string() + "\"";
    require(std::system(compile.c_str()) == 0, "generated production controller Python must compile");
    std::filesystem::remove_all(root);

    std::cout << "v0.12.0 production controller tests passed\n";
    return 0;
}
