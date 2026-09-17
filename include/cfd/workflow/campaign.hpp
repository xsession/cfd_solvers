#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <string_view>
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

[[nodiscard]] std::vector<CampaignCase> generate_factorial_campaign(std::span<const CampaignParameter> parameters);
[[nodiscard]] std::vector<CampaignCase> generate_latin_hypercube_campaign(std::span<const CampaignParameter> parameters,std::size_t samples,CampaignRandom rng={});
[[nodiscard]] std::string render_campaign_template(std::string_view templ,const CampaignCase& row,bool strict=true);
void validate_solver_adapter_descriptor(const SolverAdapterDescriptor& adapter);
[[nodiscard]] bool solver_has_capability(const SolverAdapterDescriptor& adapter,SolverCapability capability);
void update_case_status(std::vector<CampaignRegistryEntry>& registry,CampaignRegistryEntry entry);
[[nodiscard]] CampaignSummary summarize_campaign_registry(std::span<const CampaignRegistryEntry> registry);
[[nodiscard]] std::vector<double> finite_difference_gradient(const std::function<double(std::span<const double>)>& objective,std::span<const double> parameters,const FiniteDifferenceGradientConfig& config={});
[[nodiscard]] std::vector<double> gradient_descent_update(std::span<const double> parameters,std::span<const double> gradient,double step_size);

} // namespace cfd::workflow
