#pragma once

#include "cfd/workflow/campaign.hpp"

#include <cstddef>
#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::workflow {

enum class CampaignServerRuntime { native, docker };

struct CampaignServerDescriptor {
    std::string name;
    std::string host{"localhost"};
    std::string user;
    std::filesystem::path remote_root;
    CampaignServerRuntime runtime{CampaignServerRuntime::native};
    std::size_t slots{1};
    std::vector<std::string> tags;
    std::string docker_image;
    bool online{true};
};

struct MultiServerSchedulingOptions {
    bool pending_only{true};
    bool require_online{true};
    std::string required_tag;
    std::size_t max_cases_per_server{};
};

struct RemoteCommandPlan {
    std::string server_name;
    std::string description;
    std::vector<std::string> argv;
    std::string display_command;
    bool uses_ssh{};
};

struct ServerCaseAssignment {
    std::string case_id;
    std::string server_name;
    std::filesystem::path local_case_dir;
    std::filesystem::path remote_case_dir;
    SolverCommandPlan solver_plan;
};

struct MultiServerCampaignPlan {
    std::vector<ServerCaseAssignment> assignments;
    std::vector<RemoteCommandPlan> commands;
    std::size_t unassigned_cases{};
    std::size_t active_servers{};
};

struct MultiServerPlanWriteOptions {
    bool overwrite{true};
    std::string assignments_filename{"multiserver_assignments.tsv"};
    std::string commands_filename{"multiserver_commands.sh"};
};

struct MultiServerPlanWriteResult {
    std::filesystem::path assignments_path;
    std::filesystem::path commands_path;
    std::size_t assignments_written{};
    std::size_t commands_written{};
};



enum class RemoteJobState { planned, submitted, running, done, failed, cancelled, lost };

struct RemoteHostHealthCheck {
    std::string server_name;
    std::string check_name;
    RemoteCommandPlan command;
    bool required{true};
};

struct RemoteSupervisedJobPlan {
    std::string case_id;
    std::string server_name;
    std::filesystem::path remote_case_dir;
    std::filesystem::path run_dir;
    std::filesystem::path stdout_path;
    std::filesystem::path stderr_path;
    std::filesystem::path pid_path;
    std::filesystem::path exit_code_path;
    RemoteCommandPlan launch_command;
    RemoteCommandPlan status_command;
    RemoteCommandPlan cancel_command;
    RemoteCommandPlan fetch_logs_command;
    SolverCommandPlan solver_plan;
};

struct MultiServerExecutionPlan {
    std::vector<RemoteHostHealthCheck> health_checks;
    std::vector<RemoteSupervisedJobPlan> jobs;
    std::size_t native_jobs{};
    std::size_t docker_jobs{};
    std::size_t remote_jobs{};
};

struct MultiServerExecutionWriteOptions {
    bool overwrite{true};
    std::string health_filename{"multiserver_health.sh"};
    std::string launch_filename{"multiserver_launch.sh"};
    std::string status_filename{"multiserver_status.sh"};
    std::string fetch_logs_filename{"multiserver_fetch_logs.sh"};
    std::string jobs_filename{"multiserver_jobs.tsv"};
};

struct MultiServerExecutionWriteResult {
    std::filesystem::path health_script_path;
    std::filesystem::path launch_script_path;
    std::filesystem::path status_script_path;
    std::filesystem::path fetch_logs_script_path;
    std::filesystem::path jobs_path;
    std::size_t health_checks_written{};
    std::size_t jobs_written{};
};

struct RemoteJobStatusRecord {
    std::string case_id;
    std::string server_name;
    RemoteJobState state{RemoteJobState::planned};
    int exit_code{};
    std::string raw_state;
};



struct RemoteSupervisionOptions {
    std::size_t max_attempts{3};
    std::size_t retry_delay_seconds{5};
    std::size_t log_tail_lines{120};
    bool include_ssh_probe{true};
    bool include_docker_probe{true};
    bool redact_sensitive_values{true};
    std::vector<std::string> sensitive_markers{"TOKEN=", "PASSWORD=", "SECRET=", "KEY=", "AWS_", "GITHUB_"};
};

struct RemoteAccessProbe {
    std::string server_name;
    std::string probe_name;
    RemoteCommandPlan command;
    bool required{true};
};

struct RemoteRetryCommandPlan {
    std::string case_id;
    std::string server_name;
    std::size_t max_attempts{};
    std::size_t retry_delay_seconds{};
    RemoteCommandPlan command;
};

struct RemoteLogTailPlan {
    std::string case_id;
    std::string server_name;
    RemoteCommandPlan stdout_tail_command;
    RemoteCommandPlan stderr_tail_command;
};

struct RemoteDashboardCase {
    std::string case_id;
    std::string server_name;
    std::string state;
    int exit_code{};
    std::filesystem::path stdout_path;
    std::filesystem::path stderr_path;
};

struct MultiServerSupervisionPlan {
    std::vector<RemoteAccessProbe> access_probes;
    std::vector<RemoteRetryCommandPlan> retry_launches;
    std::vector<RemoteLogTailPlan> log_tails;
    std::vector<RemoteDashboardCase> dashboard_cases;
    std::size_t sensitive_commands_redacted{};
};

struct MultiServerSupervisionWriteOptions {
    bool overwrite{true};
    std::string access_filename{"multiserver_access_checks.sh"};
    std::string retry_launch_filename{"multiserver_retry_launch.sh"};
    std::string tail_logs_filename{"multiserver_tail_logs.sh"};
    std::string dashboard_filename{"multiserver_dashboard.json"};
    std::string summary_filename{"multiserver_supervision.tsv"};
};

struct MultiServerSupervisionWriteResult {
    std::filesystem::path access_script_path;
    std::filesystem::path retry_launch_script_path;
    std::filesystem::path tail_logs_script_path;
    std::filesystem::path dashboard_json_path;
    std::filesystem::path summary_path;
    std::size_t access_checks_written{};
    std::size_t retry_launches_written{};
    std::size_t log_tail_commands_written{};
    std::size_t dashboard_cases_written{};
};


struct RemoteControllerConfig {
    std::string bind_host{"127.0.0.1"};
    std::size_t port{8080};
    std::string token_env_var{"CFD_CONTROLLER_TOKEN"};
    std::string campaign_root_env_var{"CFD_CAMPAIGN_ROOT"};
    std::size_t log_tail_lines{200};
    std::size_t dashboard_refresh_seconds{5};
    std::string dashboard_title{"cfd_solvers campaign dashboard"};
    bool expose_event_snapshot{true};
    bool expose_event_stream{true};
    std::size_t event_stream_heartbeat_seconds{15};
    bool require_authentication{true};
    std::string token_file_env_var{"CFD_CONTROLLER_TOKEN_FILE"};
    bool enable_tls{true};
    std::string tls_cert_env_var{"CFD_CONTROLLER_TLS_CERT"};
    std::string tls_key_env_var{"CFD_CONTROLLER_TLS_KEY"};
    bool enable_event_history{true};
    std::size_t max_event_history_bytes{10485760};
    bool enable_live_probes{true};
    std::size_t probe_timeout_seconds{15};
    bool allow_mutating_actions{true};
    bool require_token_for_mutation{true};
    std::vector<std::string> control_actions{"stop", "extend", "checkpoint", "flush"};
};

struct RemoteControllerRoute {
    std::string method;
    std::string path;
    std::string purpose;
    bool requires_token{};
    bool mutates{};
};

struct MultiServerControllerPlan {
    RemoteControllerConfig config;
    std::vector<RemoteControllerRoute> routes;
    std::vector<RemoteDashboardCase> cases;
    std::string openapi_json;
    std::string status_json;
    std::string controller_script;
    std::string start_script;
    std::string index_html;
    std::string dashboard_js;
    std::string dashboard_css;
    std::string events_ndjson;
    std::string tokens_example_json;
    std::string systemd_unit;
    std::string reverse_proxy_config;
};

struct MultiServerControllerWriteOptions {
    bool overwrite{true};
    std::string controller_dirname{"controller"};
    std::string controller_script_filename{"cfd_controller.py"};
    std::string start_script_filename{"multiserver_controller.sh"};
    std::string openapi_filename{"openapi.json"};
    std::string status_filename{"status.json"};
    std::string routes_filename{"routes.tsv"};
    std::string env_filename{"env.example"};
    std::string readme_filename{"README.md"};
    std::string index_filename{"index.html"};
    std::string dashboard_js_filename{"dashboard.js"};
    std::string dashboard_css_filename{"dashboard.css"};
    std::string events_filename{"events.ndjson"};
    std::string tokens_example_filename{"controller.tokens.example.json"};
    std::string systemd_filename{"cfd-controller.service"};
    std::string reverse_proxy_filename{"Caddyfile.example"};
};

struct MultiServerControllerWriteResult {
    std::filesystem::path controller_script_path;
    std::filesystem::path start_script_path;
    std::filesystem::path openapi_path;
    std::filesystem::path status_path;
    std::filesystem::path routes_path;
    std::filesystem::path env_path;
    std::filesystem::path readme_path;
    std::filesystem::path index_path;
    std::filesystem::path dashboard_js_path;
    std::filesystem::path dashboard_css_path;
    std::filesystem::path events_path;
    std::filesystem::path tokens_example_path;
    std::filesystem::path systemd_path;
    std::filesystem::path reverse_proxy_path;
    std::size_t routes_written{};
    std::size_t cases_written{};
    std::size_t static_assets_written{};
};

struct DockerDeployConfig {
    std::string image_name{"cfd-solvers:latest"};
    std::string service_name{"cfd-campaign-worker"};
    std::string compose_project_name{"cfd_solvers"};
    std::string container_campaign_dir{"/work/campaigns"};
    std::string host_campaign_dir{"./campaigns"};
    std::string worker_command{"cfd-solve particle-campaign-local-runner"};
    std::size_t worker_replicas{1};
    std::size_t worker_cpus{};
    std::string worker_memory;
    bool enable_gpu_runtime{};
    bool include_manager_service{true};
    std::string manager_command{"cfd-solve particle-campaign-control"};
    std::size_t dashboard_port{8080};
    std::vector<std::pair<std::string,std::string>> environment;
};

struct DockerDeployFile {
    std::filesystem::path relative_path;
    std::string contents;
};

struct DockerDeployResult {
    std::filesystem::path root;
    std::vector<DockerDeployFile> files;
    std::size_t files_written{};
};

[[nodiscard]] std::string campaign_server_runtime_name(CampaignServerRuntime runtime);
[[nodiscard]] CampaignServerRuntime campaign_server_runtime_from_name(std::string_view name);
void validate_campaign_server_descriptor(const CampaignServerDescriptor& server);
[[nodiscard]] MultiServerCampaignPlan plan_multi_server_campaign(const std::filesystem::path& local_root,
                                                                 std::span<const CampaignRegistryEntry> registry,
                                                                 const SolverAdapterDescriptor& adapter,
                                                                 std::span<const CampaignServerDescriptor> servers,
                                                                 const MultiServerSchedulingOptions& options={});
[[nodiscard]] MultiServerPlanWriteResult write_multiserver_plan_files(const std::filesystem::path& root,
                                                                      const MultiServerCampaignPlan& plan,
                                                                      const MultiServerPlanWriteOptions& options={});


[[nodiscard]] std::string remote_job_state_name(RemoteJobState state);
[[nodiscard]] RemoteJobState remote_job_state_from_name(std::string_view name);
[[nodiscard]] MultiServerExecutionPlan plan_multiserver_execution(const MultiServerCampaignPlan& campaign_plan,
                                                                  std::span<const CampaignServerDescriptor> servers,
                                                                  std::string run_dirname=".cfd_run");
[[nodiscard]] MultiServerExecutionWriteResult write_multiserver_execution_files(const std::filesystem::path& root,
                                                                                const MultiServerExecutionPlan& plan,
                                                                                const MultiServerExecutionWriteOptions& options={});
[[nodiscard]] std::vector<RemoteJobStatusRecord> parse_remote_job_status_table(std::string_view text);
[[nodiscard]] CampaignSummary apply_remote_job_status_to_registry(std::vector<CampaignRegistryEntry>& registry,
                                                                  std::span<const RemoteJobStatusRecord> records);

[[nodiscard]] std::string redact_sensitive_command_display(std::string_view command, std::span<const std::string> sensitive_markers);
[[nodiscard]] MultiServerSupervisionPlan plan_multiserver_supervision(const MultiServerExecutionPlan& execution_plan,
                                                                      std::span<const CampaignServerDescriptor> servers,
                                                                      const RemoteSupervisionOptions& options={});
[[nodiscard]] std::string build_multiserver_dashboard_json(std::span<const RemoteDashboardCase> cases);
[[nodiscard]] MultiServerSupervisionWriteResult write_multiserver_supervision_files(const std::filesystem::path& root,
                                                                                    const MultiServerSupervisionPlan& plan,
                                                                                    const MultiServerSupervisionWriteOptions& options={});

[[nodiscard]] std::string build_multiserver_controller_openapi_json(std::span<const RemoteControllerRoute> routes,const RemoteControllerConfig& config);
[[nodiscard]] std::string build_multiserver_controller_status_json(std::span<const RemoteDashboardCase> cases);
[[nodiscard]] std::string build_multiserver_controller_events_ndjson(std::span<const RemoteDashboardCase> cases);
[[nodiscard]] std::string build_multiserver_controller_dashboard_html(const RemoteControllerConfig& config);
[[nodiscard]] std::string build_multiserver_controller_dashboard_js(const RemoteControllerConfig& config);
[[nodiscard]] std::string build_multiserver_controller_dashboard_css();
[[nodiscard]] std::string build_multiserver_controller_script(const RemoteControllerConfig& config);
[[nodiscard]] MultiServerControllerPlan plan_multiserver_controller(const MultiServerSupervisionPlan& supervision_plan,const RemoteControllerConfig& config={});
[[nodiscard]] MultiServerControllerWriteResult write_multiserver_controller_files(const std::filesystem::path& root,const MultiServerControllerPlan& plan,const MultiServerControllerWriteOptions& options={});

[[nodiscard]] std::vector<DockerDeployFile> generate_docker_deploy_files(const DockerDeployConfig& config);
[[nodiscard]] DockerDeployResult write_docker_deploy_system(const std::filesystem::path& root,
                                                            const DockerDeployConfig& config={},
                                                            bool overwrite=true);

} // namespace cfd::workflow
