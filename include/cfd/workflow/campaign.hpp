#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <unordered_map>
#include <vector>

namespace cfd::workflow {

enum class CampaignParameterKind { continuous, discrete };

struct CampaignParameter {
    std::string name;
    CampaignParameterKind kind{CampaignParameterKind::continuous};
    double minimum{};
    double maximum{};
    std::vector<std::string> levels;
};

struct CampaignCase {
    std::string case_id;
    std::unordered_map<std::string,std::string> values;
};

struct CampaignRandom {
    std::uint64_t state{0x243f6a8885a308d3ULL};
    [[nodiscard]] double uniform_open01();
};

enum class SolverCapability { residuals, probes, restart, compare, performance, control, gui };

struct SolverAdapterDescriptor {
    std::string name;
    std::string native_binary;
    std::string container_binary;
    std::vector<SolverCapability> capabilities;
    std::vector<std::string> control_actions;
};

enum class CampaignCaseStatus { pending, running, done, failed, stopped };

struct CampaignRegistryEntry {
    std::string case_id;
    CampaignCaseStatus status{CampaignCaseStatus::pending};
    double objective{};
    std::size_t iterations{};
    std::string message;
};

struct CampaignSummary {
    std::size_t pending{},running{},done{},failed{},stopped{};
    double best_objective{};
    std::string best_case_id;
};


struct FiniteDifferenceGradientConfig {
    double absolute_step{1.0e-6};
    double relative_step{1.0e-4};
    bool central{true};
};

struct CampaignTemplateFile {
    std::filesystem::path relative_path;
    std::string contents;
};

struct CampaignWriteOptions {
    bool overwrite{false};
    std::string registry_filename{"registry.tsv"};
    std::string manifest_filename{"template_manifest.tsv"};
    std::string doe_row_filename{"doe_row.csv"};
};

struct CampaignPersistenceResult {
    std::filesystem::path root;
    std::size_t cases_written{};
    std::size_t files_written{};
    std::size_t registry_entries{};
};

enum class SolverRuntime { native, docker, singularity, slurm };

struct SolverRunRequest {
    SolverAdapterDescriptor adapter;
    std::filesystem::path case_directory;
    SolverRuntime runtime{SolverRuntime::native};
    std::size_t mpi_ranks{1};
    std::size_t threads{1};
    std::vector<std::string> extra_args;
    std::string docker_image;
    std::string singularity_image;
    std::string slurm_partition;
};

struct SolverCommandPlan {
    SolverRuntime runtime{SolverRuntime::native};
    std::vector<std::string> argv;
    std::vector<std::pair<std::string,std::string>> environment;
    std::string display_command;
    bool requires_scheduler{};
};

struct ResidualSample {
    std::size_t iteration{};
    std::string equation{"residual"};
    double residual{};
};

struct PerformanceMetric {
    std::string key;
    double value{};
    std::string units;
};


struct RuntimeDoctorCheck {
    std::string name;
    bool ok{};
    std::string message;
};

struct RuntimeDoctorReport {
    std::vector<RuntimeDoctorCheck> checks;
    bool ok{};
};

struct CampaignOutputDiscovery {
    std::vector<std::filesystem::path> log_files;
    std::vector<std::filesystem::path> residual_files;
    std::vector<std::filesystem::path> performance_files;
};

struct LocalSolverRunOptions {
    bool write_launcher_logs{true};
    std::string stdout_filename{"solver.stdout.log"};
    std::string stderr_filename{"solver.stderr.log"};
    std::size_t max_capture_bytes{1U<<20U};
};

struct LocalSolverRunResult {
    std::string case_id;
    CampaignCaseStatus status{CampaignCaseStatus::pending};
    int exit_code{};
    std::string stdout_text;
    std::string stderr_text;
    std::vector<ResidualSample> residuals;
    std::vector<PerformanceMetric> performance;
    CampaignOutputDiscovery discovered_outputs;
};

struct LocalCampaignRunOptions {
    SolverAdapterDescriptor adapter;
    SolverRuntime runtime{SolverRuntime::native};
    std::size_t mpi_ranks{1};
    std::size_t threads{1};
    std::vector<std::string> extra_args;
    bool pending_only{true};
    LocalSolverRunOptions solver_run;
};

struct LocalCampaignRunResult {
    std::vector<LocalSolverRunResult> case_results;
    CampaignSummary summary;
    std::size_t launched{};
};

struct CampaignOptimizationConfig {
    std::size_t max_iterations{32};
    double step_size{0.1};
    double tolerance{1.0e-6};
    FiniteDifferenceGradientConfig finite_difference{};
};

struct CampaignOptimizationResult {
    std::vector<double> parameters;
    double objective{};
    std::size_t iterations{};
    bool converged{};
    std::vector<double> last_gradient;
};

[[nodiscard]] std::vector<CampaignCase> generate_factorial_campaign(std::span<const CampaignParameter> parameters);
[[nodiscard]] std::vector<CampaignCase> generate_latin_hypercube_campaign(std::span<const CampaignParameter> parameters,std::size_t samples,CampaignRandom rng={});
[[nodiscard]] std::string render_campaign_template(std::string_view templ,const CampaignCase& row,bool strict=true);
void validate_solver_adapter_descriptor(const SolverAdapterDescriptor& adapter);
[[nodiscard]] bool solver_has_capability(const SolverAdapterDescriptor& adapter,SolverCapability capability);
void update_case_status(std::vector<CampaignRegistryEntry>& registry,CampaignRegistryEntry entry);
[[nodiscard]] CampaignSummary summarize_campaign_registry(std::span<const CampaignRegistryEntry> registry);

[[nodiscard]] CampaignPersistenceResult write_campaign_case_folders(const std::filesystem::path& root,std::span<const CampaignCase> cases,std::span<const CampaignTemplateFile> templates,const CampaignWriteOptions& options={});
void save_campaign_registry(const std::filesystem::path& path,std::span<const CampaignRegistryEntry> registry);
[[nodiscard]] std::vector<CampaignRegistryEntry> load_campaign_registry(const std::filesystem::path& path);
[[nodiscard]] std::string campaign_status_name(CampaignCaseStatus status);
[[nodiscard]] CampaignCaseStatus campaign_status_from_name(std::string_view name);
[[nodiscard]] SolverCommandPlan build_solver_command_plan(const SolverRunRequest& request);

[[nodiscard]] RuntimeDoctorReport doctor_solver_runtime(const SolverRunRequest& request);
[[nodiscard]] CampaignOutputDiscovery discover_campaign_outputs(const std::filesystem::path& case_directory);
[[nodiscard]] LocalSolverRunResult run_local_solver_case(const SolverRunRequest& request,const LocalSolverRunOptions& options={});
[[nodiscard]] LocalCampaignRunResult run_local_campaign(const std::filesystem::path& root,const LocalCampaignRunOptions& options,const std::string& registry_filename="registry.tsv");
[[nodiscard]] CampaignCaseStatus detect_solver_outcome_from_log(std::string_view log);
[[nodiscard]] std::vector<ResidualSample> parse_residual_history(std::string_view log);
[[nodiscard]] std::vector<PerformanceMetric> parse_performance_metrics(std::string_view log);
[[nodiscard]] CampaignOptimizationResult run_gradient_descent_campaign(const std::function<double(std::span<const double>)>& objective,std::vector<double> initial,const CampaignOptimizationConfig& config={});
[[nodiscard]] std::vector<double> finite_difference_gradient(const std::function<double(std::span<const double>)>& objective,std::span<const double> parameters,const FiniteDifferenceGradientConfig& config={});
[[nodiscard]] std::vector<double> gradient_descent_update(std::span<const double> parameters,std::span<const double> gradient,double step_size);

} // namespace cfd::workflow
