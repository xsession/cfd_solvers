#pragma once

#include "cfd/workflow/campaign.hpp"

#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::workflow {

enum class CampaignControlAction { stop, extend, checkpoint, flush };

struct CampaignControlRequest {
    std::string case_id;
    CampaignControlAction action{CampaignControlAction::stop};
    std::string value;
    std::string reason;
};

struct CampaignControlOptions {
    std::string control_dirname{".cfd_control"};
    bool overwrite_latest{true};
};

struct CampaignControlDirective {
    std::string case_id;
    CampaignControlAction action{CampaignControlAction::stop};
    std::string value;
    std::string reason;
    std::filesystem::path path;
};

enum class SchedulerJobState { unknown, queued, running, done, failed, cancelled };

struct SchedulerJobRecord {
    std::string job_id;
    std::string case_id;
    SchedulerJobState state{SchedulerJobState::unknown};
    double elapsed_seconds{};
    std::string raw_state;
};

[[nodiscard]] std::string campaign_control_action_name(CampaignControlAction action);
[[nodiscard]] CampaignControlAction campaign_control_action_from_name(std::string_view name);
void validate_campaign_control_request(const SolverAdapterDescriptor& adapter,const CampaignControlRequest& request);
[[nodiscard]] CampaignControlDirective write_campaign_control_directive(const std::filesystem::path& root,const SolverAdapterDescriptor& adapter,const CampaignControlRequest& request,const CampaignControlOptions& options={});
[[nodiscard]] std::vector<CampaignControlDirective> discover_campaign_control_directives(const std::filesystem::path& case_directory,const CampaignControlOptions& options={});

[[nodiscard]] std::string scheduler_job_state_name(SchedulerJobState state);
[[nodiscard]] SchedulerJobState scheduler_job_state_from_name(std::string_view name);
[[nodiscard]] std::vector<SchedulerJobRecord> parse_slurm_queue_table(std::string_view text);
[[nodiscard]] CampaignSummary apply_scheduler_records_to_registry(std::vector<CampaignRegistryEntry>& registry,std::span<const SchedulerJobRecord> records);
[[nodiscard]] CampaignSummary refresh_campaign_status_from_outputs(const std::filesystem::path& root,std::vector<CampaignRegistryEntry>& registry);
[[nodiscard]] CampaignSummary refresh_campaign_registry_from_outputs(const std::filesystem::path& root,const std::string& registry_filename="registry.tsv");

} // namespace cfd::workflow
