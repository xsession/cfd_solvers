#include "cfd/workflow/campaign_control.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace cfd::workflow {
namespace {
std::string trim(std::string value){
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    return value.substr(first,last-first+1U);
}
std::string lower_copy(std::string_view text){
    std::string out(text);
    for(char& ch:out)ch=static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return out;
}
std::string upper_copy(std::string_view text){
    std::string out(text);
    for(char& ch:out)ch=static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    return out;
}
bool adapter_allows_action(const SolverAdapterDescriptor& adapter,std::string_view action){
    if(!solver_has_capability(adapter,SolverCapability::control))return false;
    return std::find(adapter.control_actions.begin(),adapter.control_actions.end(),std::string(action))!=adapter.control_actions.end();
}
std::string escape_value(std::string_view value){
    std::string out;out.reserve(value.size());
    for(char ch:value){
        if(ch=='\\')out+="\\\\";
        else if(ch=='\n')out+="\\n";
        else if(ch=='\r')out+="\\r";
        else out.push_back(ch);
    }
    return out;
}
std::string unescape_value(std::string_view value){
    std::string out;out.reserve(value.size());
    for(std::size_t i=0;i<value.size();++i){
        if(value[i]=='\\'&&i+1U<value.size()){
            const char n=value[++i];
            if(n=='n')out.push_back('\n');
            else if(n=='r')out.push_back('\r');
            else out.push_back(n);
        }else out.push_back(value[i]);
    }
    return out;
}
std::unordered_map<std::string,std::string> parse_key_values(const std::string& text){
    std::unordered_map<std::string,std::string> out;std::stringstream lines(text);std::string line;
    while(std::getline(lines,line)){
        const auto pos=line.find('=');
        if(pos==std::string::npos)continue;
        out[trim(line.substr(0,pos))]=unescape_value(trim(line.substr(pos+1U)));
    }
    return out;
}
std::string read_file(const std::filesystem::path& path){
    std::ifstream in(path,std::ios::binary);if(!in)return {};std::ostringstream out;out<<in.rdbuf();return out.str();
}
std::string control_payload(const CampaignControlRequest& request){
    std::ostringstream out;
    out<<"case_id="<<escape_value(request.case_id)<<'\n';
    out<<"action="<<campaign_control_action_name(request.action)<<'\n';
    out<<"value="<<escape_value(request.value)<<'\n';
    out<<"reason="<<escape_value(request.reason)<<'\n';
    return out.str();
}
bool parse_size_token(std::string_view token,std::size_t& out){
    try{std::size_t pos=0U;auto s=std::string(token);auto v=std::stoull(s,&pos);if(pos!=s.size())return false;out=static_cast<std::size_t>(v);return true;}catch(...){return false;}
}
bool parse_double_token(std::string_view token,double& out){
    try{std::size_t pos=0U;auto s=std::string(token);out=std::stod(s,&pos);return pos==s.size();}catch(...){return false;}
}
std::vector<std::string> tokens(std::string line){
    for(char& ch:line)if(ch==','||ch==';'||ch=='\t')ch=' ';
    std::stringstream ss(line);std::vector<std::string> out;std::string token;while(ss>>token)out.push_back(token);return out;
}
std::size_t max_iteration(std::span<const ResidualSample> residuals){
    std::size_t out=0U;for(const auto& sample:residuals)out=std::max(out,sample.iteration);return out;
}
double objective_from_residuals(std::span<const ResidualSample> residuals,double fallback){
    if(!residuals.empty())return residuals.back().residual;
    return fallback;
}
void update_registry_message(CampaignRegistryEntry& entry,const std::string& prefix,std::size_t residual_count,std::size_t performance_count){
    entry.message=prefix+" residuals="+std::to_string(residual_count)+" perf="+std::to_string(performance_count);
}
}

std::string campaign_control_action_name(CampaignControlAction action){
    switch(action){
        case CampaignControlAction::stop:return "stop";
        case CampaignControlAction::extend:return "extend";
        case CampaignControlAction::checkpoint:return "checkpoint";
        case CampaignControlAction::flush:return "flush";
    }
    return "stop";
}

CampaignControlAction campaign_control_action_from_name(std::string_view name){
    const auto n=lower_copy(name);
    if(n=="stop")return CampaignControlAction::stop;
    if(n=="extend")return CampaignControlAction::extend;
    if(n=="checkpoint")return CampaignControlAction::checkpoint;
    if(n=="flush")return CampaignControlAction::flush;
    throw std::invalid_argument("unknown campaign control action");
}

void validate_campaign_control_request(const SolverAdapterDescriptor& adapter,const CampaignControlRequest& request){
    validate_solver_adapter_descriptor(adapter);
    if(request.case_id.empty())throw std::invalid_argument("campaign control requires case_id");
    const auto action=campaign_control_action_name(request.action);
    if(!adapter_allows_action(adapter,action))throw std::invalid_argument("solver adapter does not support control action: "+action);
}

CampaignControlDirective write_campaign_control_directive(const std::filesystem::path& root,const SolverAdapterDescriptor& adapter,const CampaignControlRequest& request,const CampaignControlOptions& options){
    if(root.empty()||options.control_dirname.empty())throw std::invalid_argument("control directive requires root and control directory");
    validate_campaign_control_request(adapter,request);
    const auto case_dir=root/request.case_id;
    std::error_code ec;
    if(!std::filesystem::exists(case_dir,ec)||!std::filesystem::is_directory(case_dir,ec))throw std::runtime_error("campaign control case directory is missing: "+request.case_id);
    const auto control_dir=case_dir/options.control_dirname;std::filesystem::create_directories(control_dir);
    const auto action=campaign_control_action_name(request.action);
    const auto payload=control_payload(request);
    const auto target=control_dir/(action+".directive");
    if(std::filesystem::exists(target,ec)&&!options.overwrite_latest)throw std::runtime_error("control directive already exists: "+target.string());
    {std::ofstream out(target,std::ios::trunc);if(!out)throw std::runtime_error("failed to write control directive");out<<payload;}
    if(options.overwrite_latest){std::ofstream latest(control_dir/"latest.directive",std::ios::trunc);if(latest)latest<<payload;}
    return {request.case_id,request.action,request.value,request.reason,target};
}

std::vector<CampaignControlDirective> discover_campaign_control_directives(const std::filesystem::path& case_directory,const CampaignControlOptions& options){
    std::vector<CampaignControlDirective> out;
    if(case_directory.empty()||options.control_dirname.empty())return out;
    const auto control_dir=case_directory/options.control_dirname;std::error_code ec;
    if(!std::filesystem::exists(control_dir,ec))return out;
    for(std::filesystem::directory_iterator it(control_dir,ec),end;it!=end&&!ec;it.increment(ec)){
        if(ec||!it->is_regular_file(ec))continue;
        const auto path=it->path();
        if(path.extension()!=".directive")continue;
        auto kv=parse_key_values(read_file(path));
        auto case_it=kv.find("case_id");auto action_it=kv.find("action");
        if(case_it==kv.end()||action_it==kv.end())continue;
        CampaignControlDirective item;item.case_id=case_it->second;item.action=campaign_control_action_from_name(action_it->second);item.path=path;
        if(auto v=kv.find("value");v!=kv.end())item.value=v->second;
        if(auto v=kv.find("reason");v!=kv.end())item.reason=v->second;
        out.push_back(std::move(item));
    }
    std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){return std::tie(a.case_id,a.path)<std::tie(b.case_id,b.path);});
    return out;
}

std::string scheduler_job_state_name(SchedulerJobState state){
    switch(state){
        case SchedulerJobState::unknown:return "unknown";
        case SchedulerJobState::queued:return "queued";
        case SchedulerJobState::running:return "running";
        case SchedulerJobState::done:return "done";
        case SchedulerJobState::failed:return "failed";
        case SchedulerJobState::cancelled:return "cancelled";
    }
    return "unknown";
}

SchedulerJobState scheduler_job_state_from_name(std::string_view name){
    const auto n=upper_copy(name);
    if(n=="PD"||n=="PENDING"||n=="QUEUED"||n=="CONFIGURING")return SchedulerJobState::queued;
    if(n=="R"||n=="RUNNING"||n=="COMPLETING")return SchedulerJobState::running;
    if(n=="CD"||n=="COMPLETED"||n=="DONE"||n=="SUCCESS")return SchedulerJobState::done;
    if(n=="F"||n=="FAILED"||n=="TIMEOUT"||n=="NODE_FAIL"||n=="OUT_OF_MEMORY"||n=="ERROR")return SchedulerJobState::failed;
    if(n=="CA"||n=="CANCELLED"||n=="CANCELED"||n=="STOPPED")return SchedulerJobState::cancelled;
    return SchedulerJobState::unknown;
}

std::vector<SchedulerJobRecord> parse_slurm_queue_table(std::string_view text){
    std::vector<SchedulerJobRecord> out;std::stringstream lines{std::string(text)};std::string line;
    while(std::getline(lines,line)){
        auto row=tokens(line);if(row.empty())continue;
        if(upper_copy(row[0])=="JOBID"||row[0].rfind("#",0U)==0U)continue;
        if(row.size()<3U)continue;
        SchedulerJobRecord record;record.job_id=row[0];record.case_id=row[1];record.raw_state=row[2];record.state=scheduler_job_state_from_name(row[2]);
        if(row.size()>=4U){double elapsed{};if(parse_double_token(row[3],elapsed))record.elapsed_seconds=elapsed;}
        out.push_back(std::move(record));
    }
    return out;
}

CampaignSummary apply_scheduler_records_to_registry(std::vector<CampaignRegistryEntry>& registry,std::span<const SchedulerJobRecord> records){
    for(const auto& record:records){
        if(record.case_id.empty())continue;
        auto it=std::find_if(registry.begin(),registry.end(),[&](const auto& entry){return entry.case_id==record.case_id;});
        if(it==registry.end()){registry.push_back({record.case_id,CampaignCaseStatus::pending,0.0,0U,"scheduler discovered job="+record.job_id});it=std::prev(registry.end());}
        switch(record.state){
            case SchedulerJobState::queued:
            case SchedulerJobState::running:it->status=CampaignCaseStatus::running;break;
            case SchedulerJobState::done:it->status=CampaignCaseStatus::done;break;
            case SchedulerJobState::failed:it->status=CampaignCaseStatus::failed;break;
            case SchedulerJobState::cancelled:it->status=CampaignCaseStatus::stopped;break;
            case SchedulerJobState::unknown:break;
        }
        it->message="scheduler job="+record.job_id+" state="+record.raw_state;
    }
    return summarize_campaign_registry(registry);
}

CampaignSummary refresh_campaign_status_from_outputs(const std::filesystem::path& root,std::vector<CampaignRegistryEntry>& registry){
    if(root.empty())throw std::invalid_argument("campaign output refresh requires root");
    for(auto& entry:registry){
        const auto case_dir=root/entry.case_id;std::error_code ec;
        if(!std::filesystem::exists(case_dir,ec))continue;
        const auto discovery=discover_campaign_outputs(case_dir);
        std::string combined;
        for(const auto& path:discovery.log_files){combined+=read_file(path);combined+='\n';}
        if(combined.empty())continue;
        std::vector<ResidualSample> residuals=parse_residual_history(combined);
        std::vector<PerformanceMetric> perf=parse_performance_metrics(combined);
        for(const auto& path:discovery.residual_files){auto v=parse_residual_history(read_file(path));residuals.insert(residuals.end(),v.begin(),v.end());}
        for(const auto& path:discovery.performance_files){auto v=parse_performance_metrics(read_file(path));perf.insert(perf.end(),v.begin(),v.end());}
        std::sort(residuals.begin(),residuals.end(),[](const auto& a,const auto& b){return std::tie(a.iteration,a.equation,a.residual)<std::tie(b.iteration,b.equation,b.residual);});
        const auto outcome=detect_solver_outcome_from_log(combined);
        if(outcome!=CampaignCaseStatus::running)entry.status=outcome;
        else if(entry.status==CampaignCaseStatus::pending)entry.status=CampaignCaseStatus::running;
        entry.iterations=max_iteration(residuals);
        entry.objective=objective_from_residuals(residuals,entry.objective);
        update_registry_message(entry,"refreshed",residuals.size(),perf.size());
    }
    return summarize_campaign_registry(registry);
}

CampaignSummary refresh_campaign_registry_from_outputs(const std::filesystem::path& root,const std::string& registry_filename){
    if(root.empty()||registry_filename.empty())throw std::invalid_argument("campaign registry refresh requires root and registry filename");
    auto registry=load_campaign_registry(root/registry_filename);
    auto summary=refresh_campaign_status_from_outputs(root,registry);
    save_campaign_registry(root/registry_filename,registry);
    return summary;
}

} // namespace cfd::workflow
