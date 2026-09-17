#include "cfd/workflow/deploy.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace cfd::workflow {
namespace {
std::string lower_copy(std::string_view text){
    std::string out(text);
    for(char& ch:out)ch=static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    return out;
}
std::string shell_quote(std::string_view text){
    if(text.empty())return "''";
    const bool safe=std::all_of(text.begin(),text.end(),[](char ch){
        return std::isalnum(static_cast<unsigned char>(ch))||ch=='_'||ch=='-'||ch=='.'||ch=='/'||ch==':'||ch=='='||ch==',';
    });
    if(safe)return std::string(text);
    std::string out{"'"};
    for(char ch:text){
        if(ch=='\'')out+="'\\''";
        else out.push_back(ch);
    }
    out.push_back('\'');
    return out;
}
std::string join_display(std::span<const std::string> argv){
    std::ostringstream out;
    for(std::size_t i=0;i<argv.size();++i){
        if(i)out<<' ';
        out<<shell_quote(argv[i]);
    }
    return out.str();
}
bool has_tag(const CampaignServerDescriptor& server,const std::string& tag){
    if(tag.empty())return true;
    return std::find(server.tags.begin(),server.tags.end(),tag)!=server.tags.end();
}
std::string ssh_target(const CampaignServerDescriptor& server){
    if(server.user.empty())return server.host;
    return server.user+"@"+server.host;
}
bool is_localhost(std::string_view host){
    return host.empty()||host=="localhost"||host=="127.0.0.1"||host=="::1";
}
RemoteCommandPlan make_command(const CampaignServerDescriptor& server,std::string description,std::vector<std::string> remote_argv){
    RemoteCommandPlan cmd;cmd.server_name=server.name;cmd.description=std::move(description);
    if(is_localhost(server.host)){
        cmd.argv=std::move(remote_argv);
        cmd.uses_ssh=false;
    } else {
        cmd.argv={"ssh",ssh_target(server),join_display(remote_argv)};
        cmd.uses_ssh=true;
    }
    cmd.display_command=join_display(cmd.argv);
    return cmd;
}
RemoteCommandPlan make_rsync_command(const CampaignServerDescriptor& server,const std::filesystem::path& local_case,const std::filesystem::path& remote_case){
    RemoteCommandPlan cmd;cmd.server_name=server.name;cmd.description="sync "+local_case.filename().string()+" to "+server.name;
    const auto source=(local_case.string()+"/");
    if(is_localhost(server.host))cmd.argv={"rsync","-az",source,remote_case.string()+"/"};
    else cmd.argv={"rsync","-az",source,ssh_target(server)+":"+remote_case.string()+"/"};
    cmd.display_command=join_display(cmd.argv);cmd.uses_ssh=!is_localhost(server.host);return cmd;
}
void ensure_relative(const std::filesystem::path& path){
    if(path.empty()||path.is_absolute())throw std::invalid_argument("deploy file path must be relative");
    for(const auto& part:path)if(part.string()=="..")throw std::invalid_argument("deploy file path must not escape root");
}
std::string docker_env_block(const DockerDeployConfig& config,const std::string& indent){
    std::ostringstream out;
    out<<indent<<"CFD_CAMPAIGN_ROOT: "<<config.container_campaign_dir<<"\n";
    out<<indent<<"CFD_DEPLOY_PROJECT: "<<config.compose_project_name<<"\n";
    for(const auto& [key,value]:config.environment)out<<indent<<key<<": "<<value<<"\n";
    return out.str();
}
}

std::string campaign_server_runtime_name(CampaignServerRuntime runtime){
    switch(runtime){
        case CampaignServerRuntime::native:return "native";
        case CampaignServerRuntime::docker:return "docker";
    }
    return "native";
}

CampaignServerRuntime campaign_server_runtime_from_name(std::string_view name){
    const auto n=lower_copy(name);
    if(n=="native")return CampaignServerRuntime::native;
    if(n=="docker")return CampaignServerRuntime::docker;
    throw std::invalid_argument("unknown campaign server runtime");
}

void validate_campaign_server_descriptor(const CampaignServerDescriptor& server){
    if(server.name.empty())throw std::invalid_argument("campaign server requires a name");
    if(server.host.empty())throw std::invalid_argument("campaign server requires a host");
    if(server.remote_root.empty())throw std::invalid_argument("campaign server requires remote root");
    if(server.slots==0U)throw std::invalid_argument("campaign server requires at least one slot");
    if(server.runtime==CampaignServerRuntime::docker&&server.docker_image.empty())throw std::invalid_argument("docker campaign server requires image");
}

MultiServerCampaignPlan plan_multi_server_campaign(const std::filesystem::path& local_root,
                                                    std::span<const CampaignRegistryEntry> registry,
                                                    const SolverAdapterDescriptor& adapter,
                                                    std::span<const CampaignServerDescriptor> servers,
                                                    const MultiServerSchedulingOptions& options){
    if(local_root.empty())throw std::invalid_argument("multi-server campaign planning requires local root");
    validate_solver_adapter_descriptor(adapter);
    if(servers.empty())throw std::invalid_argument("multi-server campaign planning requires servers");
    std::vector<const CampaignServerDescriptor*> active;
    for(const auto& server:servers){
        validate_campaign_server_descriptor(server);
        if(options.require_online&&!server.online)continue;
        if(!has_tag(server,options.required_tag))continue;
        active.push_back(&server);
    }
    if(active.empty())throw std::invalid_argument("no active campaign servers match scheduling options");
    std::vector<std::size_t> assigned(active.size(),0U);
    std::vector<std::size_t> limits(active.size(),0U);
    for(std::size_t i=0;i<active.size();++i){
        limits[i]=active[i]->slots;
        if(options.max_cases_per_server>0U)limits[i]=std::min(limits[i],options.max_cases_per_server);
    }
    MultiServerCampaignPlan plan;plan.active_servers=active.size();
    for(const auto& entry:registry){
        if(options.pending_only&&entry.status!=CampaignCaseStatus::pending)continue;
        std::size_t best=active.size();
        for(std::size_t i=0;i<active.size();++i){
            if(assigned[i]>=limits[i])continue;
            if(best==active.size()||std::tie(assigned[i],active[i]->name)<std::tie(assigned[best],active[best]->name))best=i;
        }
        if(best==active.size()){++plan.unassigned_cases;continue;}
        const auto& server=*active[best];++assigned[best];
        const auto local_case=local_root/entry.case_id;
        const auto remote_case=server.remote_root/entry.case_id;
        SolverRunRequest request;request.adapter=adapter;request.case_directory=remote_case;request.mpi_ranks=1U;request.threads=1U;
        if(server.runtime==CampaignServerRuntime::docker){request.runtime=SolverRuntime::docker;request.docker_image=server.docker_image;}
        else request.runtime=SolverRuntime::native;
        ServerCaseAssignment assignment;assignment.case_id=entry.case_id;assignment.server_name=server.name;assignment.local_case_dir=local_case;assignment.remote_case_dir=remote_case;assignment.solver_plan=build_solver_command_plan(request);
        plan.commands.push_back(make_command(server,"create remote case directory",{"mkdir","-p",remote_case.string()}));
        plan.commands.push_back(make_rsync_command(server,local_case,remote_case));
        plan.commands.push_back(make_command(server,"run "+entry.case_id+" on "+server.name,assignment.solver_plan.argv));
        plan.assignments.push_back(std::move(assignment));
    }
    return plan;
}

MultiServerPlanWriteResult write_multiserver_plan_files(const std::filesystem::path& root,const MultiServerCampaignPlan& plan,const MultiServerPlanWriteOptions& options){
    if(root.empty()||options.assignments_filename.empty()||options.commands_filename.empty())throw std::invalid_argument("multi-server plan writer requires root and filenames");
    std::filesystem::create_directories(root);
    const auto assignments_path=root/options.assignments_filename;
    const auto commands_path=root/options.commands_filename;
    if(!options.overwrite&&(std::filesystem::exists(assignments_path)||std::filesystem::exists(commands_path)))throw std::runtime_error("multi-server plan files already exist");
    {
        std::ofstream out(assignments_path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write multi-server assignments");
        out<<"case_id\tserver\tlocal_case_dir\tremote_case_dir\truntime\tcommand\n";
        for(const auto& a:plan.assignments){
            out<<a.case_id<<'\t'<<a.server_name<<'\t'<<a.local_case_dir.generic_string()<<'\t'<<a.remote_case_dir.generic_string()<<'\t'<<static_cast<int>(a.solver_plan.runtime)<<'\t'<<a.solver_plan.display_command<<'\n';
        }
    }
    {
        std::ofstream out(commands_path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write multi-server commands");
        out<<"#!/usr/bin/env sh\nset -eu\n";
        for(const auto& cmd:plan.commands){
            out<<"\n# "<<cmd.description<<"\n"<<cmd.display_command<<"\n";
        }
    }
    std::filesystem::permissions(commands_path,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
    return {assignments_path,commands_path,plan.assignments.size(),plan.commands.size()};
}



std::string remote_job_state_name(RemoteJobState state){
    switch(state){
        case RemoteJobState::planned:return "planned";
        case RemoteJobState::submitted:return "submitted";
        case RemoteJobState::running:return "running";
        case RemoteJobState::done:return "done";
        case RemoteJobState::failed:return "failed";
        case RemoteJobState::cancelled:return "cancelled";
        case RemoteJobState::lost:return "lost";
    }
    return "planned";
}

RemoteJobState remote_job_state_from_name(std::string_view name){
    const auto n=lower_copy(name);
    if(n=="planned")return RemoteJobState::planned;
    if(n=="submitted"||n=="queued")return RemoteJobState::submitted;
    if(n=="running")return RemoteJobState::running;
    if(n=="done"||n=="complete"||n=="completed"||n=="success")return RemoteJobState::done;
    if(n=="failed"||n=="error")return RemoteJobState::failed;
    if(n=="cancelled"||n=="canceled"||n=="stopped")return RemoteJobState::cancelled;
    if(n=="lost"||n=="missing")return RemoteJobState::lost;
    throw std::invalid_argument("unknown remote job state");
}

namespace {
const CampaignServerDescriptor& find_server_by_name(std::span<const CampaignServerDescriptor> servers,const std::string& name){
    for(const auto& server:servers)if(server.name==name)return server;
    throw std::invalid_argument("assignment references unknown campaign server: "+name);
}
RemoteCommandPlan make_fetch_logs_command(const CampaignServerDescriptor& server,const std::filesystem::path& remote_run,const std::string& case_id){
    RemoteCommandPlan cmd;cmd.server_name=server.name;cmd.description="fetch logs for "+case_id+" from "+server.name;
    const auto destination=(std::filesystem::path("remote_logs")/server.name/case_id).generic_string()+"/";
    const auto source=remote_run.generic_string()+"/";
    if(is_localhost(server.host))cmd.argv={"sh","-lc","mkdir -p "+shell_quote(destination)+" && cp -a "+shell_quote(source)+". "+shell_quote(destination)};
    else cmd.argv={"rsync","-az",ssh_target(server)+":"+source,destination};
    cmd.display_command=join_display(cmd.argv);cmd.uses_ssh=!is_localhost(server.host);return cmd;
}
std::vector<std::string> split_tab_line(std::string_view line){
    std::vector<std::string> out;std::string current;
    for(char ch:line){
        if(ch=='\t'){out.push_back(current);current.clear();}
        else current.push_back(ch);
    }
    out.push_back(current);return out;
}
int parse_int_or_zero(std::string_view text){
    try{return std::stoi(std::string(text));}catch(...){return 0;}
}
void write_command_script(const std::filesystem::path& path,const std::string& title,std::span<const RemoteCommandPlan> commands){
    std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write execution script");
    out<<"#!/usr/bin/env sh\nset -eu\n\n# "<<title<<"\n";
    for(const auto& command:commands)out<<"\n# ["<<command.server_name<<"] "<<command.description<<"\n"<<command.display_command<<"\n";
    std::filesystem::permissions(path,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
}
}

MultiServerExecutionPlan plan_multiserver_execution(const MultiServerCampaignPlan& campaign_plan,
                                                     std::span<const CampaignServerDescriptor> servers,
                                                     std::string run_dirname){
    if(run_dirname.empty())throw std::invalid_argument("multi-server execution requires a run directory name");
    if(campaign_plan.assignments.empty())throw std::invalid_argument("multi-server execution requires assignments");
    if(servers.empty())throw std::invalid_argument("multi-server execution requires server descriptors");
    MultiServerExecutionPlan plan;
    std::vector<std::string> health_added;
    for(const auto& assignment:campaign_plan.assignments){
        const auto& server=find_server_by_name(servers,assignment.server_name);
        validate_campaign_server_descriptor(server);
        if(std::find(health_added.begin(),health_added.end(),server.name)==health_added.end()){
            health_added.push_back(server.name);
            RemoteHostHealthCheck root;root.server_name=server.name;root.check_name="remote-root";
            root.command=make_command(server,"ensure remote campaign root",{"sh","-lc","mkdir -p "+shell_quote(server.remote_root.generic_string())+" && test -d "+shell_quote(server.remote_root.generic_string())});
            plan.health_checks.push_back(std::move(root));
            RemoteHostHealthCheck runtime;runtime.server_name=server.name;runtime.check_name=(server.runtime==CampaignServerRuntime::docker?"docker-runtime":"native-runtime");
            if(server.runtime==CampaignServerRuntime::docker)runtime.command=make_command(server,"verify docker runtime",{"sh","-lc","command -v docker >/dev/null && docker --version >/dev/null"});
            else runtime.command=make_command(server,"verify native solver binary",{"sh","-lc","command -v "+shell_quote(assignment.solver_plan.argv.empty()?std::string{}:assignment.solver_plan.argv.front())+" >/dev/null"});
            plan.health_checks.push_back(std::move(runtime));
            RemoteHostHealthCheck space;space.server_name=server.name;space.check_name="disk-space";space.required=false;
            space.command=make_command(server,"show remote disk space",{"sh","-lc","df -P "+shell_quote(server.remote_root.generic_string())+" | tail -n 1"});
            plan.health_checks.push_back(std::move(space));
        }
        RemoteSupervisedJobPlan job;job.case_id=assignment.case_id;job.server_name=server.name;job.remote_case_dir=assignment.remote_case_dir;job.solver_plan=assignment.solver_plan;
        job.run_dir=assignment.remote_case_dir/run_dirname;
        job.stdout_path=job.run_dir/"solver.stdout.log";job.stderr_path=job.run_dir/"solver.stderr.log";job.pid_path=job.run_dir/"solver.pid";job.exit_code_path=job.run_dir/"exit_code";
        const auto solver_cmd=assignment.solver_plan.display_command;
        std::string launch="mkdir -p "+shell_quote(job.run_dir.generic_string())+
            " && rm -f "+shell_quote(job.exit_code_path.generic_string())+
            " && cd "+shell_quote(job.remote_case_dir.generic_string())+
            " && nohup sh -lc "+shell_quote(solver_cmd+"; echo $? > "+shell_quote(job.exit_code_path.generic_string()))+
            " > "+shell_quote(job.stdout_path.generic_string())+
            " 2> "+shell_quote(job.stderr_path.generic_string())+
            " & echo $! > "+shell_quote(job.pid_path.generic_string());
        job.launch_command=make_command(server,"launch supervised "+assignment.case_id,{"sh","-lc",launch});
        std::string status="state=lost; ec=0; if [ -f "+shell_quote(job.exit_code_path.generic_string())+
            " ]; then ec=$(cat "+shell_quote(job.exit_code_path.generic_string())+
            "); if [ \"$ec\" = 0 ]; then state=done; else state=failed; fi; elif [ -f "+shell_quote(job.pid_path.generic_string())+
            " ] && kill -0 $(cat "+shell_quote(job.pid_path.generic_string())+
            ") 2>/dev/null; then state=running; else state=lost; fi; printf '%s\\t%s\\t%s\\t%s\\n' "+
            shell_quote(job.case_id)+" "+shell_quote(job.server_name)+" \"$state\" \"$ec\"";
        job.status_command=make_command(server,"status "+assignment.case_id,{"sh","-lc",status});
        std::string cancel="if [ -f "+shell_quote(job.pid_path.generic_string())+
            " ]; then kill $(cat "+shell_quote(job.pid_path.generic_string())+
            ") 2>/dev/null || true; fi; echo cancelled > "+shell_quote(job.run_dir.generic_string()+"/cancelled")+" ; true";
        job.cancel_command=make_command(server,"cancel "+assignment.case_id,{"sh","-lc",cancel});
        job.fetch_logs_command=make_fetch_logs_command(server,job.run_dir,job.case_id);
        if(server.runtime==CampaignServerRuntime::docker)++plan.docker_jobs;else ++plan.native_jobs;
        if(!is_localhost(server.host))++plan.remote_jobs;
        plan.jobs.push_back(std::move(job));
    }
    return plan;
}

MultiServerExecutionWriteResult write_multiserver_execution_files(const std::filesystem::path& root,const MultiServerExecutionPlan& plan,const MultiServerExecutionWriteOptions& options){
    if(root.empty())throw std::invalid_argument("execution writer requires a root path");
    if(options.health_filename.empty()||options.launch_filename.empty()||options.status_filename.empty()||options.fetch_logs_filename.empty()||options.jobs_filename.empty())throw std::invalid_argument("execution writer requires filenames");
    std::filesystem::create_directories(root);
    MultiServerExecutionWriteResult result;
    result.health_script_path=root/options.health_filename;result.launch_script_path=root/options.launch_filename;result.status_script_path=root/options.status_filename;result.fetch_logs_script_path=root/options.fetch_logs_filename;result.jobs_path=root/options.jobs_filename;
    if(!options.overwrite){
        for(const auto& path:{result.health_script_path,result.launch_script_path,result.status_script_path,result.fetch_logs_script_path,result.jobs_path})if(std::filesystem::exists(path))throw std::runtime_error("execution file already exists: "+path.string());
    }
    std::vector<RemoteCommandPlan> health,launch,status,fetch;
    for(const auto& check:plan.health_checks)health.push_back(check.command);
    for(const auto& job:plan.jobs){launch.push_back(job.launch_command);status.push_back(job.status_command);fetch.push_back(job.fetch_logs_command);}    
    write_command_script(result.health_script_path,"multi-server health checks",health);
    write_command_script(result.launch_script_path,"multi-server supervised launch",launch);
    write_command_script(result.status_script_path,"multi-server status probes",status);
    write_command_script(result.fetch_logs_script_path,"multi-server log fetch",fetch);
    {
        std::ofstream out(result.jobs_path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write multi-server jobs table");
        out<<"case_id\tserver\tremote_case_dir\trun_dir\tstdout\tstderr\tpid\texit_code\truntime\tlaunch\tstatus\tcancel\tfetch_logs\n";
        for(const auto& job:plan.jobs){
            out<<job.case_id<<'\t'<<job.server_name<<'\t'<<job.remote_case_dir.generic_string()<<'\t'<<job.run_dir.generic_string()<<'\t'
               <<job.stdout_path.generic_string()<<'\t'<<job.stderr_path.generic_string()<<'\t'<<job.pid_path.generic_string()<<'\t'<<job.exit_code_path.generic_string()<<'\t'
               <<static_cast<int>(job.solver_plan.runtime)<<'\t'<<job.launch_command.display_command<<'\t'<<job.status_command.display_command<<'\t'
               <<job.cancel_command.display_command<<'\t'<<job.fetch_logs_command.display_command<<'\n';
        }
    }
    result.health_checks_written=plan.health_checks.size();result.jobs_written=plan.jobs.size();return result;
}

std::vector<RemoteJobStatusRecord> parse_remote_job_status_table(std::string_view text){
    std::vector<RemoteJobStatusRecord> out;std::stringstream input{std::string(text)};std::string line;bool first=true;
    while(std::getline(input,line)){
        if(line.empty())continue;
        auto fields=split_tab_line(line);
        if(first&&fields.size()>=3U&&fields[0]=="case_id"){first=false;continue;}
        first=false;
        if(fields.size()<3U)throw std::runtime_error("malformed remote job status row");
        RemoteJobStatusRecord row;row.case_id=fields[0];row.server_name=fields[1];row.raw_state=fields[2];row.state=remote_job_state_from_name(fields[2]);
        if(fields.size()>=4U)row.exit_code=parse_int_or_zero(fields[3]);
        out.push_back(std::move(row));
    }
    return out;
}

CampaignSummary apply_remote_job_status_to_registry(std::vector<CampaignRegistryEntry>& registry,std::span<const RemoteJobStatusRecord> records){
    for(const auto& record:records){
        CampaignCaseStatus status=CampaignCaseStatus::running;
        switch(record.state){
            case RemoteJobState::planned:
            case RemoteJobState::submitted:
            case RemoteJobState::running:status=CampaignCaseStatus::running;break;
            case RemoteJobState::done:status=CampaignCaseStatus::done;break;
            case RemoteJobState::failed:
            case RemoteJobState::lost:status=CampaignCaseStatus::failed;break;
            case RemoteJobState::cancelled:status=CampaignCaseStatus::stopped;break;
        }
        update_case_status(registry,{record.case_id,status,0.0,0U,"remote "+remote_job_state_name(record.state)+" on "+record.server_name});
    }
    return summarize_campaign_registry(registry);
}


std::string redact_sensitive_command_display(std::string_view command, std::span<const std::string> sensitive_markers){
    std::string out(command);
    for(const auto& marker:sensitive_markers){
        if(marker.empty())continue;
        std::size_t pos=0;
        while((pos=out.find(marker,pos))!=std::string::npos){
            const auto value_start=pos+marker.size();
            auto value_end=value_start;
            while(value_end<out.size()&&!std::isspace(static_cast<unsigned char>(out[value_end]))&&out[value_end]!='\''&&out[value_end]!='"'&&out[value_end]!=';')++value_end;
            out.replace(value_start,value_end-value_start,"<redacted>");
            pos=value_start+10U;
        }
    }
    return out;
}

namespace {
std::string json_escape(std::string_view text){
    std::string out;
    out.reserve(text.size()+8U);
    for(char ch:text){
        switch(ch){
            case '\\':out+="\\\\";break;
            case '"':out+="\\\"";break;
            case '\n':out+="\\n";break;
            case '\r':out+="\\r";break;
            case '\t':out+="\\t";break;
            default:
                if(static_cast<unsigned char>(ch)<0x20U){
                    out+="?";
                } else out.push_back(ch);
        }
    }
    return out;
}
RemoteCommandPlan make_retry_command(const CampaignServerDescriptor& server,
                                     const RemoteSupervisedJobPlan& job,
                                     const RemoteSupervisionOptions& options,
                                     const std::string& launch_display){
    const auto attempts=std::max<std::size_t>(1U,options.max_attempts);
    std::ostringstream script;
    script<<"attempt=1; ec=0; while :; do "
          <<launch_display
          <<" && exit 0; ec=$?; if [ $attempt -ge "<<attempts
          <<" ]; then exit $ec; fi; sleep "<<options.retry_delay_seconds
          <<"; attempt=$((attempt+1)); done";
    return make_command(server,"retry launch "+job.case_id,{"sh","-lc",script.str()});
}
RemoteCommandPlan make_tail_command(const CampaignServerDescriptor& server,
                                    const std::string& description,
                                    const std::filesystem::path& log_path,
                                    std::size_t lines){
    return make_command(server,description,{"sh","-lc","tail -n "+std::to_string(lines)+" "+shell_quote(log_path.generic_string())+" 2>/dev/null || true"});
}
void write_supervision_command_script(const std::filesystem::path& path,const std::string& title,std::span<const RemoteCommandPlan> commands){
    std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write supervision script");
    out<<"#!/usr/bin/env sh\nset -eu\n\n# "<<title<<"\n";
    for(const auto& command:commands)out<<"\n# ["<<command.server_name<<"] "<<command.description<<"\n"<<command.display_command<<"\n";
    std::filesystem::permissions(path,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
}
}

MultiServerSupervisionPlan plan_multiserver_supervision(const MultiServerExecutionPlan& execution_plan,
                                                        std::span<const CampaignServerDescriptor> servers,
                                                        const RemoteSupervisionOptions& options){
    if(execution_plan.jobs.empty())throw std::invalid_argument("multi-server supervision requires jobs");
    if(servers.empty())throw std::invalid_argument("multi-server supervision requires server descriptors");
    if(options.max_attempts==0U)throw std::invalid_argument("retry attempts must be non-zero");
    if(options.log_tail_lines==0U)throw std::invalid_argument("log tail line count must be non-zero");
    MultiServerSupervisionPlan plan;
    std::vector<std::string> probed;
    for(const auto& job:execution_plan.jobs){
        const auto& server=find_server_by_name(servers,job.server_name);
        validate_campaign_server_descriptor(server);
        if(std::find(probed.begin(),probed.end(),server.name)==probed.end()){
            probed.push_back(server.name);
            if(options.include_ssh_probe){
                RemoteAccessProbe probe;probe.server_name=server.name;probe.probe_name=is_localhost(server.host)?"local-access":"ssh-batch-access";probe.required=true;
                if(is_localhost(server.host))probe.command=make_command(server,"verify local execution shell",{"sh","-lc","test -d / && echo local-ok"});
                else probe.command=make_command(server,"verify SSH batch access",{"sh","-lc","echo ssh-ok && test -d "+shell_quote(server.remote_root.generic_string())});
                plan.access_probes.push_back(std::move(probe));
            }
            if(options.include_docker_probe&&server.runtime==CampaignServerRuntime::docker){
                RemoteAccessProbe probe;probe.server_name=server.name;probe.probe_name="docker-daemon";probe.required=true;
                probe.command=make_command(server,"probe Docker daemon",{"sh","-lc","docker info >/dev/null && docker ps --format '{{.ID}}' | head -n 1 >/dev/null"});
                plan.access_probes.push_back(std::move(probe));
            }
            RemoteAccessProbe disk;disk.server_name=server.name;disk.probe_name="campaign-root-writable";disk.required=true;
            disk.command=make_command(server,"verify campaign root is writable",{"sh","-lc","mkdir -p "+shell_quote(server.remote_root.generic_string())+" && test -w "+shell_quote(server.remote_root.generic_string())});
            plan.access_probes.push_back(std::move(disk));
        }
        std::string launch_display=job.launch_command.display_command;
        if(options.redact_sensitive_values){
            const auto redacted=redact_sensitive_command_display(launch_display,options.sensitive_markers);
            if(redacted!=launch_display)++plan.sensitive_commands_redacted;
            launch_display=redacted;
        }
        RemoteRetryCommandPlan retry;retry.case_id=job.case_id;retry.server_name=job.server_name;retry.max_attempts=options.max_attempts;retry.retry_delay_seconds=options.retry_delay_seconds;
        retry.command=make_retry_command(server,job,options,launch_display);
        plan.retry_launches.push_back(std::move(retry));
        RemoteLogTailPlan tail;tail.case_id=job.case_id;tail.server_name=job.server_name;
        tail.stdout_tail_command=make_tail_command(server,"tail stdout for "+job.case_id,job.stdout_path,options.log_tail_lines);
        tail.stderr_tail_command=make_tail_command(server,"tail stderr for "+job.case_id,job.stderr_path,options.log_tail_lines);
        plan.log_tails.push_back(std::move(tail));
        RemoteDashboardCase dash;dash.case_id=job.case_id;dash.server_name=job.server_name;dash.state="planned";dash.exit_code=0;dash.stdout_path=job.stdout_path;dash.stderr_path=job.stderr_path;
        plan.dashboard_cases.push_back(std::move(dash));
    }
    return plan;
}

std::string build_multiserver_dashboard_json(std::span<const RemoteDashboardCase> cases){
    std::ostringstream out;
    out<<"{\n  \"schema\": \"cfd.multiserver.dashboard.v1\",\n  \"cases\": [\n";
    for(std::size_t i=0;i<cases.size();++i){
        const auto& c=cases[i];
        out<<"    {\"case_id\": \""<<json_escape(c.case_id)<<"\", \"server\": \""<<json_escape(c.server_name)
           <<"\", \"state\": \""<<json_escape(c.state)<<"\", \"exit_code\": "<<c.exit_code
           <<", \"stdout\": \""<<json_escape(c.stdout_path.generic_string())<<"\", \"stderr\": \""<<json_escape(c.stderr_path.generic_string())<<"\"}";
        if(i+1U<cases.size())out<<',';
        out<<"\n";
    }
    out<<"  ]\n}\n";
    return out.str();
}

MultiServerSupervisionWriteResult write_multiserver_supervision_files(const std::filesystem::path& root,
                                                                      const MultiServerSupervisionPlan& plan,
                                                                      const MultiServerSupervisionWriteOptions& options){
    if(root.empty())throw std::invalid_argument("supervision writer requires a root path");
    if(options.access_filename.empty()||options.retry_launch_filename.empty()||options.tail_logs_filename.empty()||options.dashboard_filename.empty()||options.summary_filename.empty())throw std::invalid_argument("supervision writer requires filenames");
    std::filesystem::create_directories(root);
    MultiServerSupervisionWriteResult result;
    result.access_script_path=root/options.access_filename;
    result.retry_launch_script_path=root/options.retry_launch_filename;
    result.tail_logs_script_path=root/options.tail_logs_filename;
    result.dashboard_json_path=root/options.dashboard_filename;
    result.summary_path=root/options.summary_filename;
    if(!options.overwrite){
        for(const auto& path:{result.access_script_path,result.retry_launch_script_path,result.tail_logs_script_path,result.dashboard_json_path,result.summary_path})if(std::filesystem::exists(path))throw std::runtime_error("supervision file already exists: "+path.string());
    }
    std::vector<RemoteCommandPlan> access,retry,tail;
    for(const auto& probe:plan.access_probes)access.push_back(probe.command);
    for(const auto& item:plan.retry_launches)retry.push_back(item.command);
    for(const auto& item:plan.log_tails){tail.push_back(item.stdout_tail_command);tail.push_back(item.stderr_tail_command);}    
    write_supervision_command_script(result.access_script_path,"multi-server access and runtime probes",access);
    write_supervision_command_script(result.retry_launch_script_path,"multi-server retry launch commands",retry);
    write_supervision_command_script(result.tail_logs_script_path,"multi-server log tail commands",tail);
    {
        std::ofstream out(result.dashboard_json_path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write dashboard json");
        out<<build_multiserver_dashboard_json(plan.dashboard_cases);
    }
    {
        std::ofstream out(result.summary_path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write supervision summary");
        out<<"case_id\tserver\tstate\tstdout\tstderr\n";
        for(const auto& c:plan.dashboard_cases)out<<c.case_id<<'\t'<<c.server_name<<'\t'<<c.state<<'\t'<<c.stdout_path.generic_string()<<'\t'<<c.stderr_path.generic_string()<<'\n';
    }
    result.access_checks_written=plan.access_probes.size();
    result.retry_launches_written=plan.retry_launches.size();
    result.log_tail_commands_written=plan.log_tails.size()*2U;
    result.dashboard_cases_written=plan.dashboard_cases.size();
    return result;
}


std::string build_multiserver_controller_status_json(std::span<const RemoteDashboardCase> cases){
    std::size_t planned=0U,running=0U,done=0U,failed=0U,cancelled=0U,other=0U;
    for(const auto& c:cases){
        const auto state=lower_copy(c.state);
        if(state=="planned")++planned; else if(state=="running")++running; else if(state=="done")++done; else if(state=="failed")++failed; else if(state=="cancelled"||state=="stopped")++cancelled; else ++other;
    }
    std::ostringstream out;
    out<<"{\n  \"schema\": \"cfd.multiserver.controller.status.v1\",\n";
    out<<"  \"summary\": {\"cases\": "<<cases.size()<<", \"planned\": "<<planned<<", \"running\": "<<running<<", \"done\": "<<done<<", \"failed\": "<<failed<<", \"cancelled\": "<<cancelled<<", \"other\": "<<other<<"},\n";
    out<<"  \"cases\": [\n";
    for(std::size_t i=0;i<cases.size();++i){
        const auto& c=cases[i];
        out<<"    {\"case_id\": \""<<json_escape(c.case_id)<<"\", \"server\": \""<<json_escape(c.server_name)
           <<"\", \"state\": \""<<json_escape(c.state)<<"\", \"exit_code\": "<<c.exit_code
           <<", \"stdout\": \""<<json_escape(c.stdout_path.generic_string())<<"\", \"stderr\": \""<<json_escape(c.stderr_path.generic_string())<<"\"}";
        if(i+1U<cases.size())out<<',';
        out<<"\n";
    }
    out<<"  ]\n}\n";
    return out.str();
}

std::string build_multiserver_controller_events_ndjson(std::span<const RemoteDashboardCase> cases){
    std::ostringstream out;
    out<<"{\"event\":\"summary\",\"schema\":\"cfd.multiserver.controller.events.v1\",\"cases\":"<<cases.size()<<"}\n";
    for(const auto& c:cases){
        out<<"{\"event\":\"case\",\"case_id\":\""<<json_escape(c.case_id)
           <<"\",\"server\":\""<<json_escape(c.server_name)
           <<"\",\"state\":\""<<json_escape(c.state)
           <<"\",\"exit_code\":"<<c.exit_code<<"}\n";
    }
    return out.str();
}

std::string build_multiserver_controller_dashboard_html(const RemoteControllerConfig& config){
    std::ostringstream html;
    html<<"<!doctype html>\n<html lang=\"en\">\n<head>\n";
    html<<"  <meta charset=\"utf-8\">\n  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n";
    html<<"  <title>"<<json_escape(config.dashboard_title)<<"</title>\n";
    html<<"  <link rel=\"stylesheet\" href=\"/dashboard.css\">\n</head>\n<body>\n";
    html<<"  <main data-refresh-seconds=\""<<config.dashboard_refresh_seconds<<"\">\n";
    html<<"    <header><div><p class=\"eyebrow\">cfd_solvers controller</p><h1>"<<json_escape(config.dashboard_title)<<"</h1></div><span id=\"health\">offline</span></header>\n";
    html<<"    <section id=\"summary\" class=\"cards\"></section>\n";
    html<<"    <section><h2>Cases</h2><table><thead><tr><th>Case</th><th>Server</th><th>State</th><th>Exit</th><th>Logs</th><th>Control</th></tr></thead><tbody id=\"cases\"></tbody></table></section>\n";
    html<<"    <section><h2>Log preview</h2><pre id=\"log\">Select stdout or stderr.</pre></section>\n";
    html<<"  </main>\n  <script src=\"/dashboard.js\"></script>\n</body>\n</html>\n";
    return html.str();
}

std::string build_multiserver_controller_dashboard_js(const RemoteControllerConfig& config){
    std::ostringstream js;
    js<<"const refreshSeconds = "<<config.dashboard_refresh_seconds<<";\n";
    js<<"const controlActions = ";
    js<<"[";
    for(std::size_t i=0;i<config.control_actions.size();++i){ if(i)js<<","; js<<"'"<<json_escape(config.control_actions[i])<<"'"; }
    js<<"];\n";
    js<<R"JS(
async function getJson(path) {
  const r = await fetch(path, {cache: 'no-store', headers: authHeaders()});
  if (!r.ok) throw new Error(path + ' -> ' + r.status);
  return await r.json();
}
let controllerToken = sessionStorage.getItem('cfd-controller-token') || '';
function authHeaders() { return controllerToken ? {'Authorization': 'Bearer ' + controllerToken} : {}; }
function ensureToken() {
  if (!controllerToken) controllerToken = window.prompt('Controller access token:', '') || '';
  if (controllerToken) sessionStorage.setItem('cfd-controller-token', controllerToken);
  return controllerToken;
}
function setHealth(text, ok) {
  const el = document.getElementById('health');
  el.textContent = text;
  el.className = ok ? 'ok' : 'bad';
}
let eventSource;
function connectEvents() {
  if (!window.EventSource) return;
  const token = ensureToken();
  eventSource = new EventSource('/api/events/stream' + (token ? '?token=' + encodeURIComponent(token) : ''));
  eventSource.addEventListener('status', event => {
    const status = JSON.parse(event.data);
    setHealth('live', true);
    renderSummary(status.summary || {});
    renderCases(status.cases || []);
  });
  eventSource.onerror = () => setHealth('reconnecting', false);
}
function renderSummary(summary) {
  const root = document.getElementById('summary');
  const keys = ['cases','planned','running','done','failed','cancelled','other'];
  root.innerHTML = keys.map(k => `<article><span>${k}</span><strong>${summary[k] ?? 0}</strong></article>`).join('');
}
async function showLog(caseId, stream) {
  ensureToken();
  const r = await fetch(`/api/cases/${encodeURIComponent(caseId)}/logs/${stream}`, {cache: 'no-store', headers: authHeaders()});
  document.getElementById('log').textContent = await r.text();
}
async function sendControl(caseId, action) {
  ensureToken();
  const r = await fetch(`/api/cases/${encodeURIComponent(caseId)}/control/${action}`, {method: 'POST', headers: authHeaders()});
  const text = await r.text();
  if (!r.ok) window.alert(text);
  await refresh();
}
function renderCases(cases) {
  const body = document.getElementById('cases');
  body.innerHTML = cases.map(c => {
    const logs = `<button onclick="showLog('${c.case_id}','stdout')">stdout</button> <button onclick="showLog('${c.case_id}','stderr')">stderr</button>`;
    const control = controlActions.map(a => `<button onclick="sendControl('${c.case_id}','${a}')">${a}</button>`).join(' ');
    return `<tr><td>${c.case_id}</td><td>${c.server}</td><td><span class="state ${c.state}">${c.state}</span></td><td>${c.exit_code}</td><td>${logs}</td><td>${control}</td></tr>`;
  }).join('');
}
async function refresh() {
  try {
    const health = await getJson('/api/health');
    const status = await getJson('/api/status');
    setHealth(health.ok ? 'online' : 'degraded', !!health.ok);
    renderSummary(status.summary || {});
    renderCases(status.cases || []);
  } catch (err) {
    setHealth('offline: ' + err.message, false);
  }
}
refresh();
setInterval(refresh, refreshSeconds * 1000);
connectEvents();
)JS";
    return js.str();
}

std::string build_multiserver_controller_dashboard_css(){
    return R"CSS(:root { color-scheme: light dark; font-family: Inter, ui-sans-serif, system-ui, sans-serif; }
body { margin: 0; background: #0f172a; color: #e2e8f0; }
main { max-width: 1180px; margin: 0 auto; padding: 28px; }
header { display: flex; justify-content: space-between; align-items: center; gap: 16px; margin-bottom: 24px; }
h1 { margin: 0; font-size: clamp(1.8rem, 4vw, 3rem); }
.eyebrow { margin: 0 0 4px; color: #93c5fd; text-transform: uppercase; letter-spacing: .14em; font-size: .78rem; }
#health { padding: 8px 12px; border-radius: 999px; background: #334155; }
#health.ok { background: #14532d; }
#health.bad { background: #7f1d1d; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 12px; margin-bottom: 24px; }
article { background: #1e293b; border: 1px solid #334155; border-radius: 14px; padding: 14px; }
article span { color: #94a3b8; display: block; }
article strong { display: block; font-size: 2rem; margin-top: 6px; }
section { background: #111827; border: 1px solid #334155; border-radius: 16px; padding: 18px; margin-bottom: 20px; overflow-x: auto; }
table { width: 100%; border-collapse: collapse; }
th, td { padding: 10px; border-bottom: 1px solid #334155; text-align: left; }
button { border: 1px solid #64748b; background: #1e293b; color: #e2e8f0; border-radius: 8px; padding: 6px 10px; margin: 2px; cursor: pointer; }
button:hover { background: #334155; }
.state { padding: 3px 8px; border-radius: 999px; background: #334155; }
.state.running { background: #1d4ed8; }
.state.done { background: #166534; }
.state.failed { background: #991b1b; }
pre { min-height: 160px; white-space: pre-wrap; background: #020617; padding: 14px; border-radius: 12px; }
)CSS";
}

std::string build_multiserver_controller_openapi_json(std::span<const RemoteControllerRoute> routes,const RemoteControllerConfig& config){
    std::ostringstream out;
    out<<"{\n  \"openapi\": \"3.0.3\",\n  \"info\": {\"title\": \"cfd_solvers multi-server controller\", \"version\": \"0.12.0\"},\n";
    out<<"  \"servers\": [{\"url\": \"http://"<<json_escape(config.bind_host)<<":"<<config.port<<"\"}],\n";
    out<<"  \"paths\": {\n";
    for(std::size_t i=0;i<routes.size();++i){
        const auto& r=routes[i];
        out<<"    \""<<json_escape(r.path)<<"\": {\""<<lower_copy(r.method)<<"\": {\"summary\": \""<<json_escape(r.purpose)<<"\", \"x-mutates\": "<<(r.mutates?"true":"false")<<", \"x-requires-token\": "<<(r.requires_token?"true":"false")<<"}}";
        if(i+1U<routes.size())out<<',';
        out<<"\n";
    }
    out<<"  }\n}\n";
    return out.str();
}

std::string build_multiserver_controller_script(const RemoteControllerConfig& config){
    std::ostringstream py;
    py<<"#!/usr/bin/env python3\n";
    py<<"import hashlib, hmac, json, os, pathlib, ssl, subprocess, sys, time, urllib.parse\n";
    py<<"from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer\n\n";
    py<<"HERE = pathlib.Path(__file__).resolve().parent\n";
    py<<"CAMPAIGN_ROOT = pathlib.Path(os.environ.get('"<<config.campaign_root_env_var<<"', HERE.parent)).resolve()\n";
    py<<"TOKEN_ENV = '"<<config.token_env_var<<"'\n";
    py<<"TAIL_LINES = "<<config.log_tail_lines<<"\n";
    py<<"EVENT_HEARTBEAT_SECONDS = "<<config.event_stream_heartbeat_seconds<<"\n";
    py<<"TOKEN_FILE_ENV = '"<<config.token_file_env_var<<"'\n";
    py<<"REQUIRE_AUTH = "<<(config.require_authentication?"True":"False")<<"\n";
    py<<"TLS_ENABLED = "<<(config.enable_tls?"True":"False")<<"\n";
    py<<"TLS_CERT_ENV = '"<<config.tls_cert_env_var<<"'\n";
    py<<"TLS_KEY_ENV = '"<<config.tls_key_env_var<<"'\n";
    py<<"HISTORY_ENABLED = "<<(config.enable_event_history?"True":"False")<<"\n";
    py<<"MAX_HISTORY_BYTES = "<<config.max_event_history_bytes<<"\n";
    py<<"LIVE_PROBES_ENABLED = "<<(config.enable_live_probes?"True":"False")<<"\n";
    py<<"PROBE_TIMEOUT_SECONDS = "<<config.probe_timeout_seconds<<"\n";
    py<<"ALLOW_MUTATION = "<<(config.allow_mutating_actions?"True":"False")<<"\n";
    py<<"REQUIRE_TOKEN = "<<(config.require_token_for_mutation?"True":"False")<<"\n";
    py<<"ACTIONS = set(";
    py<<"[";
    for(std::size_t i=0;i<config.control_actions.size();++i){ if(i)py<<", "; py<<"'"<<config.control_actions[i]<<"'"; }
    py<<"]";
    py<<")\n\n";
    py<<"def read_status():\n    path = HERE / 'status.json'\n    if not path.exists():\n        return {'schema': 'cfd.multiserver.controller.status.v1', 'summary': {'cases': 0}, 'cases': []}\n    return json.loads(path.read_text(encoding='utf-8'))\n\ndef append_history(event, payload):\n    if not HISTORY_ENABLED:\n        return\n    path = HERE / 'events-history.ndjson'\n    record = json.dumps({'time_ns': time.time_ns(), 'event': event, 'payload': payload}, separators=(',', ':')) + '\\n'\n    if path.exists() and path.stat().st_size + len(record.encode('utf-8')) > MAX_HISTORY_BYTES:\n        rotated = HERE / 'events-history.ndjson.1'\n        path.replace(rotated)\n    with path.open('a', encoding='utf-8') as stream:\n        stream.write(record)\n\n";
    py<<"def reply(handler, code, obj, content_type='application/json'):\n    data = json.dumps(obj, indent=2).encode('utf-8') if not isinstance(obj, (bytes, bytearray)) else obj\n    handler.send_response(code)\n    handler.send_header('Content-Type', content_type)\n    handler.send_header('Content-Length', str(len(data)))\n    handler.end_headers()\n    handler.wfile.write(data)\n\n";
    py<<"def authorized(handler, required_role='viewer'):\n    if not REQUIRE_AUTH and not REQUIRE_TOKEN:\n        return True\n    header = handler.headers.get('Authorization', '')\n    query = urllib.parse.parse_qs(urllib.parse.urlparse(handler.path).query)\n    raw = header[7:] if header.startswith('Bearer ') else query.get('token', [''])[0]\n    legacy = os.environ.get(TOKEN_ENV, '')\n    if legacy and hmac.compare_digest(raw, legacy):\n        return True\n    token_path = pathlib.Path(os.environ.get(TOKEN_FILE_ENV, HERE / 'controller.tokens.json'))\n    try:\n        entries = json.loads(token_path.read_text(encoding='utf-8')).get('tokens', [])\n    except (OSError, ValueError):\n        return False\n    digest = hashlib.sha256(raw.encode('utf-8')).hexdigest()\n    levels = {'viewer': 0, 'operator': 1, 'admin': 2}\n    for entry in entries:\n        if hmac.compare_digest(digest, str(entry.get('sha256', ''))) and levels.get(entry.get('role'), -1) >= levels[required_role]:\n            return True\n    return False\n\n";
    py<<"def tail(path_text, lines):\n    path = pathlib.Path(path_text)\n    if not path.is_absolute():\n        path = CAMPAIGN_ROOT / path\n    try:\n        buf = path.read_text(errors='replace').splitlines()\n        return '\\n'.join(buf[-lines:]) + ('\\n' if buf else '')\n    except OSError as exc:\n        return 'log unavailable: %s\\n' % exc\n\n";
    py<<"class Handler(BaseHTTPRequestHandler):\n    def do_GET(self):\n        parsed = urllib.parse.urlparse(self.path)\n        parts = [p for p in parsed.path.split('/') if p]\n        status = read_status()\n        if parsed.path in ('/', '/index.html'):\n            data = (HERE / 'index.html').read_bytes()\n            reply(self, 200, data, 'text/html; charset=utf-8')\n            return\n        if parsed.path == '/dashboard.js':\n            data = (HERE / 'dashboard.js').read_bytes()\n            reply(self, 200, data, 'application/javascript; charset=utf-8')\n            return\n        if parsed.path == '/dashboard.css':\n            data = (HERE / 'dashboard.css').read_bytes()\n            reply(self, 200, data, 'text/css; charset=utf-8')\n            return\n        if parsed.path == '/api/events':\n            path = HERE / 'events.ndjson'\n            data = path.read_bytes() if path.exists() else b''\n            reply(self, 200, data, 'application/x-ndjson; charset=utf-8')\n            return\n        if parsed.path == '/api/health':\n            reply(self, 200, {'ok': True, 'schema': 'cfd.multiserver.controller.health.v1'})\n            return\n        if parsed.path == '/api/status':\n            reply(self, 200, status)\n            return\n        if parsed.path == '/api/cases':\n            reply(self, 200, {'cases': status.get('cases', [])})\n            return\n        if len(parts) == 5 and parts[:2] == ['api', 'cases'] and parts[3] == 'logs':\n            case_id, stream = parts[2], parts[4]\n            for case in status.get('cases', []):\n                if case.get('case_id') == case_id and stream in ('stdout', 'stderr'):\n                    data = tail(case.get(stream, ''), TAIL_LINES).encode('utf-8')\n                    reply(self, 200, data, 'text/plain; charset=utf-8')\n                    return\n            reply(self, 404, {'error': 'case or log stream not found'})\n            return\n        reply(self, 404, {'error': 'not found'})\n\n    def do_POST(self):\n        parsed = urllib.parse.urlparse(self.path)\n        parts = [p for p in parsed.path.split('/') if p]\n        if len(parts) == 5 and parts[:2] == ['api', 'cases'] and parts[3] == 'control':\n            if not ALLOW_MUTATION:\n                reply(self, 403, {'error': 'mutations disabled'})\n                return\n            if not authorized(self):\n                reply(self, 401, {'error': 'missing or invalid controller token'})\n                return\n            case_id, action = parts[2], parts[4]\n            if action not in ACTIONS:\n                reply(self, 400, {'error': 'unsupported action'})\n                return\n            control_dir = CAMPAIGN_ROOT / case_id / '.cfd_control'\n            control_dir.mkdir(parents=True, exist_ok=True)\n            body = 'case_id=%s\\naction=%s\\nreason=controller-api\\n' % (case_id, action)\n            (control_dir / (action + '.directive')).write_text(body, encoding='utf-8')\n            (control_dir / 'latest.directive').write_text(body, encoding='utf-8')\n            reply(self, 202, {'ok': True, 'case_id': case_id, 'action': action})\n            return\n        reply(self, 404, {'error': 'not found'})\n\n    def log_message(self, fmt, *args):\n        sys.stderr.write('controller: ' + fmt % args + '\\n')\n\n";
    py<<"def sse_frame(event, payload, event_id):\n    data = json.dumps(payload, separators=(',', ':'))\n    return ('id: %s\\nevent: %s\\ndata: %s\\n\\n' % (event_id, event, data)).encode('utf-8')\n\n";
    py<<"class StreamingHandler(Handler):\n    protocol_version = 'HTTP/1.1'\n\n    def do_GET(self):\n        if urllib.parse.urlparse(self.path).path != '/api/events/stream':\n            return super().do_GET()\n        self.send_response(200)\n        self.send_header('Content-Type', 'text/event-stream; charset=utf-8')\n        self.send_header('Cache-Control', 'no-cache')\n        self.send_header('Connection', 'keep-alive')\n        self.end_headers()\n        event_id = 0\n        last_payload = None\n        try:\n            while True:\n                payload = read_status()\n                encoded = json.dumps(payload, sort_keys=True)\n                if encoded != last_payload:\n                    event_id += 1\n                    self.wfile.write(sse_frame('status', payload, event_id))\n                    last_payload = encoded\n                else:\n                    self.wfile.write(b': heartbeat\\n\\n')\n                self.wfile.flush()\n                time.sleep(EVENT_HEARTBEAT_SECONDS)\n        except (BrokenPipeError, ConnectionResetError):\n            return\n\n";
    py<<"class ProductionHandler(StreamingHandler):\n    def reject_auth(self):\n        reply(self, 401, {'error': 'authentication required'})\n\n    def do_GET(self):\n        path = urllib.parse.urlparse(self.path).path\n        if path.startswith('/api/') and path != '/api/health' and not authorized(self, 'viewer'):\n            return self.reject_auth()\n        if path == '/api/events/history':\n            history = HERE / 'events-history.ndjson'\n            return reply(self, 200, history.read_bytes() if history.exists() else b'', 'application/x-ndjson; charset=utf-8')\n        if path == '/api/probes':\n            if not LIVE_PROBES_ENABLED:\n                return reply(self, 403, {'error': 'live probes disabled'})\n            script = HERE.parent / 'multiserver_access_checks.sh'\n            try:\n                run = subprocess.run([str(script)], cwd=HERE.parent, capture_output=True, text=True, timeout=PROBE_TIMEOUT_SECONDS, check=False)\n                payload = {'ok': run.returncode == 0, 'exit_code': run.returncode, 'stdout': run.stdout[-8192:], 'stderr': run.stderr[-8192:]}\n            except (OSError, subprocess.TimeoutExpired) as exc:\n                payload = {'ok': False, 'error': str(exc)}\n            append_history('probe', payload)\n            return reply(self, 200 if payload.get('ok') else 503, payload)\n        return super().do_GET()\n\n    def do_POST(self):\n        if not authorized(self, 'operator'):\n            return self.reject_auth()\n        append_history('control-request', {'path': urllib.parse.urlparse(self.path).path})\n        return super().do_POST()\n\n";
    py<<"if __name__ == '__main__':\n    host = os.environ.get('CFD_CONTROLLER_HOST', '"<<config.bind_host<<"')\n    port = int(os.environ.get('CFD_CONTROLLER_PORT', '"<<config.port<<"'))\n    server = ThreadingHTTPServer((host, port), ProductionHandler)\n    if TLS_ENABLED:\n        cert = os.environ.get(TLS_CERT_ENV, '')\n        key = os.environ.get(TLS_KEY_ENV, '')\n        if not cert or not key:\n            raise SystemExit('TLS enabled but certificate/key environment variables are missing')\n        context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)\n        context.minimum_version = ssl.TLSVersion.TLSv1_2\n        context.load_cert_chain(cert, key)\n        server.socket = context.wrap_socket(server.socket, server_side=True)\n    print('serving cfd controller on %s:%d campaign_root=%s tls=%s' % (host, port, CAMPAIGN_ROOT, TLS_ENABLED))\n    server.serve_forever()\n";
    return py.str();
}

MultiServerControllerPlan plan_multiserver_controller(const MultiServerSupervisionPlan& supervision_plan,const RemoteControllerConfig& config){
    if(config.port==0U)throw std::invalid_argument("controller port must be non-zero");
    if(config.token_env_var.empty())throw std::invalid_argument("controller token env var must not be empty");
    if(config.log_tail_lines==0U)throw std::invalid_argument("controller log tail line count must be non-zero");
    if(config.dashboard_refresh_seconds==0U)throw std::invalid_argument("controller dashboard refresh interval must be non-zero");
    if(config.event_stream_heartbeat_seconds==0U)throw std::invalid_argument("controller event heartbeat interval must be non-zero");
    if(config.require_authentication&&config.token_file_env_var.empty())throw std::invalid_argument("controller token file env var must not be empty");
    if(config.enable_tls&&(config.tls_cert_env_var.empty()||config.tls_key_env_var.empty()))throw std::invalid_argument("controller TLS env vars must not be empty");
    if(config.max_event_history_bytes==0U)throw std::invalid_argument("controller event history limit must be non-zero");
    if(config.probe_timeout_seconds==0U)throw std::invalid_argument("controller probe timeout must be non-zero");
    MultiServerControllerPlan plan;plan.config=config;plan.cases=supervision_plan.dashboard_cases;
    plan.routes.push_back({"GET","/","browser dashboard",false,false});
    plan.routes.push_back({"GET","/dashboard.js","dashboard JavaScript asset",false,false});
    plan.routes.push_back({"GET","/dashboard.css","dashboard stylesheet asset",false,false});
    plan.routes.push_back({"GET","/api/health","controller health check",false,false});
    plan.routes.push_back({"GET","/api/status","summary plus case status",false,false});
    plan.routes.push_back({"GET","/api/cases","case table",false,false});
    plan.routes.push_back({"GET","/api/cases/{case_id}/logs/{stream}","tail stdout or stderr metadata path",false,false});
    if(config.expose_event_snapshot){
        plan.routes.push_back({"GET","/api/events","newline-delimited dashboard event snapshot",false,false});
    }
    if(config.expose_event_stream){
        plan.routes.push_back({"GET","/api/events/stream","server-sent status event stream",false,false});
    }
    if(config.enable_event_history){
        plan.routes.push_back({"GET","/api/events/history","durable bounded event history",config.require_authentication,false});
    }
    if(config.enable_live_probes){
        plan.routes.push_back({"GET","/api/probes","run supervised SSH Docker and filesystem probes",config.require_authentication,false});
    }
    if(config.allow_mutating_actions){
        plan.routes.push_back({"POST","/api/cases/{case_id}/control/{action}","write a campaign control directive",config.require_token_for_mutation,true});
    }
    plan.openapi_json=build_multiserver_controller_openapi_json(plan.routes,config);
    plan.status_json=build_multiserver_controller_status_json(plan.cases);
    plan.controller_script=build_multiserver_controller_script(config);
    plan.index_html=build_multiserver_controller_dashboard_html(config);
    plan.dashboard_js=build_multiserver_controller_dashboard_js(config);
    plan.dashboard_css=build_multiserver_controller_dashboard_css();
    plan.events_ndjson=build_multiserver_controller_events_ndjson(plan.cases);
    plan.tokens_example_json="{\n  \"tokens\": [\n    {\"name\": \"dashboard-viewer\", \"role\": \"viewer\", \"sha256\": \"replace-with-sha256-of-random-token\"},\n    {\"name\": \"campaign-operator\", \"role\": \"operator\", \"sha256\": \"replace-with-sha256-of-another-random-token\"}\n  ]\n}\n";
    plan.systemd_unit="[Unit]\nDescription=cfd_solvers campaign controller\nAfter=network-online.target\nWants=network-online.target\n\n[Service]\nType=simple\nUser=cfd-controller\nGroup=cfd-controller\nWorkingDirectory=/var/lib/cfd-controller\nEnvironmentFile=/etc/cfd-controller/controller.env\nExecStart=/usr/bin/python3 /var/lib/cfd-controller/controller/cfd_controller.py\nRestart=on-failure\nRestartSec=3\nNoNewPrivileges=true\nPrivateTmp=true\nProtectSystem=strict\nProtectHome=true\nReadWritePaths=/var/lib/cfd-controller\nCapabilityBoundingSet=\nLockPersonality=true\nMemoryDenyWriteExecute=true\n\n[Install]\nWantedBy=multi-user.target\n";
    plan.reverse_proxy_config=":443 {\n  tls /etc/cfd-controller/tls.crt /etc/cfd-controller/tls.key\n  encode zstd gzip\n  reverse_proxy 127.0.0.1:8080 {\n    flush_interval -1\n  }\n  header {\n    Strict-Transport-Security \"max-age=31536000; includeSubDomains\"\n    X-Content-Type-Options nosniff\n    X-Frame-Options DENY\n    Referrer-Policy no-referrer\n  }\n}\n";
    std::ostringstream sh;
    sh<<"#!/usr/bin/env sh\nset -eu\nSCRIPT_DIR=$(CDPATH= cd -- \"$(dirname -- \"$0\")\" && pwd)\ncd \"$SCRIPT_DIR/controller\"\nexec python3 cfd_controller.py\n";
    plan.start_script=sh.str();
    return plan;
}

MultiServerControllerWriteResult write_multiserver_controller_files(const std::filesystem::path& root,const MultiServerControllerPlan& plan,const MultiServerControllerWriteOptions& options){
    if(root.empty())throw std::invalid_argument("controller writer requires root");
    if(options.controller_dirname.empty()||options.controller_script_filename.empty()||options.start_script_filename.empty()||options.openapi_filename.empty()||options.status_filename.empty()||options.routes_filename.empty()||options.index_filename.empty()||options.dashboard_js_filename.empty()||options.dashboard_css_filename.empty()||options.events_filename.empty()||options.tokens_example_filename.empty()||options.systemd_filename.empty()||options.reverse_proxy_filename.empty())throw std::invalid_argument("controller writer requires filenames");
    std::filesystem::create_directories(root/options.controller_dirname);
    MultiServerControllerWriteResult result;
    const auto controller_dir=root/options.controller_dirname;
    result.controller_script_path=controller_dir/options.controller_script_filename;
    result.start_script_path=root/options.start_script_filename;
    result.openapi_path=controller_dir/options.openapi_filename;
    result.status_path=controller_dir/options.status_filename;
    result.routes_path=controller_dir/options.routes_filename;
    result.env_path=controller_dir/options.env_filename;
    result.readme_path=controller_dir/options.readme_filename;
    result.index_path=controller_dir/options.index_filename;
    result.dashboard_js_path=controller_dir/options.dashboard_js_filename;
    result.dashboard_css_path=controller_dir/options.dashboard_css_filename;
    result.events_path=controller_dir/options.events_filename;
    result.tokens_example_path=controller_dir/options.tokens_example_filename;
    result.systemd_path=controller_dir/options.systemd_filename;
    result.reverse_proxy_path=controller_dir/options.reverse_proxy_filename;
    if(!options.overwrite){
        for(const auto& path:{result.controller_script_path,result.start_script_path,result.openapi_path,result.status_path,result.routes_path,result.env_path,result.readme_path,result.index_path,result.dashboard_js_path,result.dashboard_css_path,result.events_path})if(std::filesystem::exists(path))throw std::runtime_error("controller file already exists: "+path.string());
    }
    auto write_text=[](const std::filesystem::path& path,const std::string& text){std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write controller file: "+path.string());out<<text;};
    write_text(result.controller_script_path,plan.controller_script);
    write_text(result.start_script_path,plan.start_script);
    write_text(result.openapi_path,plan.openapi_json);
    write_text(result.status_path,plan.status_json);
    write_text(result.index_path,plan.index_html);
    write_text(result.dashboard_js_path,plan.dashboard_js);
    write_text(result.dashboard_css_path,plan.dashboard_css);
    write_text(result.events_path,plan.events_ndjson);
    write_text(result.tokens_example_path,plan.tokens_example_json);
    write_text(result.systemd_path,plan.systemd_unit);
    write_text(result.reverse_proxy_path,plan.reverse_proxy_config);
    {
        std::ostringstream out;out<<"method\tpath\trequires_token\tmutates\tpurpose\n";
        for(const auto& r:plan.routes)out<<r.method<<'\t'<<r.path<<'\t'<<(r.requires_token?1:0)<<'\t'<<(r.mutates?1:0)<<'\t'<<r.purpose<<'\n';
        write_text(result.routes_path,out.str());
    }
    {
        std::ostringstream env;env<<"# Store real credentials outside source control.\n"<<plan.config.token_file_env_var<<"=/etc/cfd-controller/controller.tokens.json\nCFD_CONTROLLER_HOST="<<plan.config.bind_host<<"\nCFD_CONTROLLER_PORT="<<plan.config.port<<"\n"<<plan.config.campaign_root_env_var<<"="<<root.generic_string()<<"\n"<<plan.config.tls_cert_env_var<<"=/etc/cfd-controller/tls.crt\n"<<plan.config.tls_key_env_var<<"=/etc/cfd-controller/tls.key\n";
        write_text(result.env_path,env.str());
    }
    {
        std::ostringstream readme;
        readme<<"# Multi-server controller\n\nThis directory contains a dependency-free Python controller for generated multi-server campaign supervision artifacts.\n\n";
        readme<<"Start it from the campaign root with `./"<<options.start_script_filename<<"`. Copy `"<<options.tokens_example_filename<<"` outside the source tree, replace the placeholder hashes with SHA-256 digests of random tokens, and set `"<<plan.config.token_file_env_var<<"`. Configure the TLS certificate/key variables before startup.\n\n";
        readme<<"Routes are documented in `"<<options.routes_filename<<"` and `"<<options.openapi_filename<<"`. Status is served from `"<<options.status_filename<<"`; remote workers still use the generated supervision scripts.\n";
        readme<<"Open `"<<options.index_filename<<"` or `http://127.0.0.1:"<<plan.config.port<<"/` for the generated live dashboard. It prefers the SSE status stream and retains polling as a reconnect fallback. `"<<options.events_filename<<"` stores a newline-delimited event snapshot for downstream dashboards.\n";
        write_text(result.readme_path,readme.str());
    }
    std::filesystem::permissions(result.controller_script_path,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
    std::filesystem::permissions(result.start_script_path,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
    result.routes_written=plan.routes.size();
    result.cases_written=plan.cases.size();
    result.static_assets_written=7U;
    return result;
}

std::vector<DockerDeployFile> generate_docker_deploy_files(const DockerDeployConfig& config){
    if(config.image_name.empty()||config.service_name.empty()||config.container_campaign_dir.empty()||config.worker_command.empty())throw std::invalid_argument("docker deploy config requires image, service, campaign dir and command");
    std::vector<DockerDeployFile> files;
    std::ostringstream dockerfile;
    dockerfile<<"FROM debian:bookworm-slim\n"
              <<"RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates libgomp1 && rm -rf /var/lib/apt/lists/*\n"
              <<"WORKDIR /opt/cfd_solvers\n"
              <<"COPY . /opt/cfd_solvers\n"
              <<"ENV PATH=/opt/cfd_solvers/bin:$PATH\n"
              <<"RUN mkdir -p "<<config.container_campaign_dir<<"\n"
              <<"CMD [\"/bin/sh\",\"-lc\",\""<<config.worker_command<<"\"]\n";
    files.push_back({"Dockerfile",dockerfile.str()});
    std::ostringstream compose;
    compose<<"services:\n";
    if(config.include_manager_service){
        compose<<"  campaign-manager:\n"
               <<"    image: "<<config.image_name<<"\n"
               <<"    build: .\n"
               <<"    command: [\"/bin/sh\", \"-lc\", \""<<config.manager_command<<"\"]\n"
               <<"    environment:\n"<<docker_env_block(config,"      ")
               <<"    volumes:\n"
               <<"      - "<<config.host_campaign_dir<<":"<<config.container_campaign_dir<<"\n"
               <<"    ports:\n"
               <<"      - \""<<config.dashboard_port<<":"<<config.dashboard_port<<"\"\n"
               <<"    healthcheck:\n"
               <<"      test: [\"CMD-SHELL\", \"test -d "<<config.container_campaign_dir<<"\"]\n"
               <<"      interval: 30s\n"
               <<"      timeout: 5s\n"
               <<"      retries: 3\n";
    }
    compose<<"  "<<config.service_name<<":\n"
           <<"    image: "<<config.image_name<<"\n"
           <<"    build: .\n"
           <<"    command: [\"/bin/sh\", \"-lc\", \""<<config.worker_command<<"\"]\n"
           <<"    environment:\n"<<docker_env_block(config,"      ")
           <<"    volumes:\n"
           <<"      - "<<config.host_campaign_dir<<":"<<config.container_campaign_dir<<"\n"
           <<"    healthcheck:\n"
           <<"      test: [\"CMD-SHELL\", \"test -d "<<config.container_campaign_dir<<"\"]\n"
           <<"      interval: 30s\n"
           <<"      timeout: 5s\n"
           <<"      retries: 3\n";
    if(config.worker_replicas>1U||config.worker_cpus>0U||!config.worker_memory.empty()||config.enable_gpu_runtime){
        compose<<"    deploy:\n";
        if(config.worker_replicas>1U)compose<<"      replicas: "<<config.worker_replicas<<"\n";
        if(config.worker_cpus>0U||!config.worker_memory.empty()||config.enable_gpu_runtime){
            compose<<"      resources:\n";
            if(config.worker_cpus>0U||!config.worker_memory.empty()){
                compose<<"        limits:\n";
                if(config.worker_cpus>0U)compose<<"          cpus: '"<<config.worker_cpus<<"'\n";
                if(!config.worker_memory.empty())compose<<"          memory: "<<config.worker_memory<<"\n";
            }
            if(config.enable_gpu_runtime){
                compose<<"        reservations:\n"
                       <<"          devices:\n"
                       <<"            - capabilities: [gpu]\n";
            }
        }
    }
    files.push_back({"compose.yaml",compose.str()});
    files.push_back({".dockerignore","build*/\n.git/\n*.zip\n*.patch\n*.log\n/mnt/\n"});
    files.push_back({"deploy/healthcheck.sh","#!/usr/bin/env sh\nset -eu\ntest -d \"${CFD_CAMPAIGN_ROOT:-/work/campaigns}\"\n"});
    files.push_back({"deploy/worker-entrypoint.sh","#!/usr/bin/env sh\nset -eu\n: \"${CFD_CAMPAIGN_ROOT:=/work/campaigns}\"\nmkdir -p \"$CFD_CAMPAIGN_ROOT\"\nexec /bin/sh -lc \"${CFD_WORKER_COMMAND:-cfd-solve particle-campaign-local-runner}\"\n"});
    std::ostringstream env;
    env<<"CFD_IMAGE="<<config.image_name<<"\n"
       <<"CFD_CAMPAIGN_DIR="<<config.host_campaign_dir<<"\n"
       <<"CFD_WORKER_REPLICAS="<<config.worker_replicas<<"\n";
    files.push_back({"deploy/env.example",env.str()});
    std::ostringstream servers;
    servers<<"name\thost\tuser\tslots\truntime\tremote_root\tdocker_image\ttags\n"
           <<"gpu-a\t10.0.0.11\tcfd\t2\tdocker\t/srv/cfd/campaigns\t"<<config.image_name<<"\tgpu,linux\n"
           <<"cpu-a\t10.0.0.12\tcfd\t4\tnative\t/srv/cfd/campaigns\t\tcpu,linux\n";
    files.push_back({"deploy/servers.example.tsv",servers.str()});
    std::ostringstream script;
    script<<"#!/usr/bin/env sh\nset -eu\n"
          <<"docker compose -p "<<config.compose_project_name<<" build\n"
          <<"docker compose -p "<<config.compose_project_name<<" up -d --scale "<<config.service_name<<"="<<config.worker_replicas<<"\n";
    files.push_back({"deploy/deploy.sh",script.str()});
    std::ostringstream readme;
    readme<<"# Docker deployment\n\n"
          <<"This directory contains a deterministic Docker/Compose deployment scaffold for cfd_solvers campaigns.\n\n"
          <<"1. Copy or generate campaign folders under `"<<config.host_campaign_dir<<"`.\n"
          <<"2. Build the image with `docker compose build`.\n"
          <<"3. Start workers with `docker compose up -d --scale "<<config.service_name<<"="<<config.worker_replicas<<"`.\n"
          <<"4. Use `deploy/servers.example.tsv` as the shape for multi-server placement.\n\n"
          <<"The generated files are an execution scaffold, not a security boundary. Review images, mounts and network exposure before production use.\n";
    files.push_back({"deploy/README.md",readme.str()});
    return files;
}

DockerDeployResult write_docker_deploy_system(const std::filesystem::path& root,const DockerDeployConfig& config,bool overwrite){
    if(root.empty())throw std::invalid_argument("docker deploy writer requires root");
    auto files=generate_docker_deploy_files(config);
    std::filesystem::create_directories(root);
    DockerDeployResult result;result.root=root;result.files=files;
    for(const auto& file:files){
        ensure_relative(file.relative_path);
        const auto target=root/file.relative_path;
        if(std::filesystem::exists(target)&&!overwrite)throw std::runtime_error("docker deploy file already exists: "+target.string());
        std::filesystem::create_directories(target.parent_path());
        std::ofstream out(target,std::ios::trunc);if(!out)throw std::runtime_error("failed to write docker deploy file");
        out<<file.contents;++result.files_written;
        if(target.extension()==".sh")std::filesystem::permissions(target,std::filesystem::perms::owner_exec,std::filesystem::perm_options::add);
    }
    return result;
}

} // namespace cfd::workflow
