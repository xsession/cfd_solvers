#include "cfd/workflow/deploy.hpp"

#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    using namespace cfd::workflow;
    MultiServerSupervisionPlan supervision;
    supervision.dashboard_cases.push_back({"case0001", "node-a", "running", 0, "case0001/stdout.log", "case0001/stderr.log"});

    RemoteControllerConfig config;
    config.event_stream_heartbeat_seconds = 3U;
    const auto plan = plan_multiserver_controller(supervision, config);

    require(plan.routes.size() == 12U, "default controller should include SSE and production routes");
    require(plan.openapi_json.find("/api/events/stream") != std::string::npos, "OpenAPI should advertise SSE route");
    require(plan.dashboard_js.find("new EventSource('/api/events/stream'") != std::string::npos, "dashboard should connect with EventSource");
    require(plan.dashboard_js.find("setInterval(refresh") != std::string::npos, "dashboard should retain polling fallback");
    require(plan.controller_script.find("text/event-stream") != std::string::npos, "controller should emit SSE content type");
    require(plan.controller_script.find("EVENT_HEARTBEAT_SECONDS = 3") != std::string::npos, "controller should preserve heartbeat setting");
    require(plan.controller_script.find(": heartbeat") != std::string::npos, "controller should emit heartbeat comments");
    require(plan.controller_script.find("ThreadingHTTPServer((host, port), ProductionHandler)") != std::string::npos, "controller should use production streaming handler");

    const auto script_path = std::filesystem::temp_directory_path() / "cfd_solvers_v0118_controller.py";
    {
        std::ofstream script(script_path, std::ios::trunc);
        require(static_cast<bool>(script), "temporary controller script should open");
        script << plan.controller_script;
    }
    const std::string compile_command = "python3 -m py_compile \"" + script_path.string() + "\"";
    require(std::system(compile_command.c_str()) == 0, "generated Python controller should compile");
    std::filesystem::remove(script_path);
    std::filesystem::remove_all(script_path.parent_path() / "__pycache__");

    RemoteControllerConfig disabled;
    disabled.expose_event_stream = false;
    const auto without_stream = plan_multiserver_controller(supervision, disabled);
    require(without_stream.routes.size() == 11U, "disabled SSE should remove only stream route");
    require(without_stream.openapi_json.find("/api/events/stream") == std::string::npos, "disabled SSE route should be absent from OpenAPI");

    bool rejected = false;
    try {
        RemoteControllerConfig invalid;
        invalid.event_stream_heartbeat_seconds = 0U;
        static_cast<void>(plan_multiserver_controller(supervision, invalid));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "zero heartbeat interval should be rejected");

    std::cout << "v0.11.8 controller SSE tests passed\n";
    return 0;
}
