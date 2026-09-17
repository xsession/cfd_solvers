#include "cfd/workflow/campaign.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <system_error>
#include <stdexcept>
#include <string_view>
#include <tuple>
#if !defined(_WIN32)
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace cfd::workflow {
namespace {
std::string case_id(std::size_t index){std::ostringstream out;out<<"case"<<std::setw(4)<<std::setfill('0')<<(index+1U);return out.str();}
bool truthy(std::string_view value){return !(value.empty()||value=="false"||value=="False"||value=="FALSE"||value=="0");}
std::string number_to_string(double value){std::ostringstream out;out<<std::setprecision(12)<<value;return out.str();}
std::string trim_ascii(std::string value){
    const auto first=value.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    const auto last=value.find_last_not_of(" \t\r\n");
    return value.substr(first,last-first+1U);
}
void validate_parameter(const CampaignParameter& p){
    if(p.name.empty())throw std::invalid_argument("campaign parameter requires a name");
    if(p.name=="case_id")throw std::invalid_argument("case_id is reserved");
    if(p.kind==CampaignParameterKind::continuous){if(!(p.maximum>=p.minimum)||!std::isfinite(p.minimum)||!std::isfinite(p.maximum))throw std::invalid_argument("invalid continuous campaign parameter");}
    else if(p.levels.empty())throw std::invalid_argument("discrete campaign parameter requires levels");
}
}

double CampaignRandom::uniform_open01(){state=state*2862933555777941757ULL+3037000493ULL;const std::uint64_t mantissa=(state>>11)|1ULL;return std::min(std::nextafter(1.0,0.0),std::max(std::numeric_limits<double>::min(),static_cast<double>(mantissa)*(1.0/9007199254740992.0)));}

std::vector<CampaignCase> generate_factorial_campaign(std::span<const CampaignParameter> parameters){
    if(parameters.empty())throw std::invalid_argument("factorial campaign requires parameters");
    for(const auto& p:parameters){validate_parameter(p);if(p.kind!=CampaignParameterKind::discrete)throw std::invalid_argument("factorial campaign only supports discrete parameters");}
    std::vector<CampaignCase> cases(1U);cases[0].case_id=case_id(0U);
    for(const auto& parameter:parameters){
        std::vector<CampaignCase> next;next.reserve(cases.size()*parameter.levels.size());
        for(const auto& base:cases)for(const auto& level:parameter.levels){auto row=base;row.values[parameter.name]=level;next.push_back(std::move(row));}
        cases=std::move(next);
    }
    for(std::size_t i=0;i<cases.size();++i)cases[i].case_id=case_id(i);
    return cases;
}

std::vector<CampaignCase> generate_latin_hypercube_campaign(std::span<const CampaignParameter> parameters,std::size_t samples,CampaignRandom rng){
    if(parameters.empty()||samples==0U)throw std::invalid_argument("latin hypercube campaign requires parameters and samples");
    for(const auto& p:parameters)validate_parameter(p);
    std::vector<CampaignCase> cases(samples);for(std::size_t i=0;i<samples;++i)cases[i].case_id=case_id(i);
    for(const auto& parameter:parameters){
        if(parameter.kind==CampaignParameterKind::continuous){
            std::vector<double> strata(samples);for(std::size_t i=0;i<samples;++i)strata[i]=(static_cast<double>(i)+rng.uniform_open01())/static_cast<double>(samples);
            for(std::size_t i=samples;i>1U;--i){const std::size_t j=static_cast<std::size_t>(rng.uniform_open01()*static_cast<double>(i));std::swap(strata[i-1U],strata[j]);}
            for(std::size_t i=0;i<samples;++i){const double value=parameter.minimum+(parameter.maximum-parameter.minimum)*strata[i];cases[i].values[parameter.name]=number_to_string(value);} 
        } else {
            for(std::size_t i=0;i<samples;++i)cases[i].values[parameter.name]=parameter.levels[i%parameter.levels.size()];
        }
    }
    return cases;
}

std::string render_campaign_template(std::string_view templ,const CampaignCase& row,bool strict){
    auto replace_placeholders=[&](std::string line){
        std::size_t pos=0U;
        while((pos=line.find('{',pos))!=std::string::npos){const auto end=line.find('}',pos+1U);if(end==std::string::npos)break;const auto key=line.substr(pos+1U,end-pos-1U);auto it=row.values.find(key);if(it==row.values.end()){if(strict)throw std::invalid_argument("missing campaign template placeholder: "+key);pos=end+1U;continue;}line.replace(pos,end-pos+1U,it->second);pos+=it->second.size();}
        return line;
    };
    std::stringstream in{std::string(templ)};std::ostringstream out;std::string line;std::vector<bool> include_stack;
    auto currently_included=[&](){return std::all_of(include_stack.begin(),include_stack.end(),[](bool v){return v;});};
    while(std::getline(in,line)){
        const auto if_pos=line.find("<!-- IF ");
        const auto endif_pos=line.find("<!-- ENDIF -->");
        if(if_pos!=std::string::npos){
            const auto close=line.find("-->",if_pos);if(close==std::string::npos)throw std::invalid_argument("unterminated campaign IF block");std::string expr=trim_ascii(line.substr(if_pos+8U,close-(if_pos+8U)));
            bool include=false;auto neq=expr.find("!=");auto eq=expr.find('=');
            if(neq!=std::string::npos){const auto key=trim_ascii(expr.substr(0,neq));const auto expected=trim_ascii(expr.substr(neq+2U));auto it=row.values.find(key);include=it!=row.values.end()&&it->second!=expected;}
            else if(eq!=std::string::npos){const auto key=trim_ascii(expr.substr(0,eq));const auto expected=trim_ascii(expr.substr(eq+1U));auto it=row.values.find(key);include=it!=row.values.end()&&it->second==expected;}
            else {auto it=row.values.find(expr);include=it!=row.values.end()&&truthy(it->second);}include_stack.push_back(include);continue;
        }
        if(endif_pos!=std::string::npos){if(include_stack.empty())throw std::invalid_argument("campaign ENDIF without IF");include_stack.pop_back();continue;}
        if(currently_included())out<<replace_placeholders(line)<<'\n';
    }
    if(!include_stack.empty())throw std::invalid_argument("unclosed campaign IF block");
    return out.str();
}

void validate_solver_adapter_descriptor(const SolverAdapterDescriptor& adapter){
    if(adapter.name.empty()||adapter.native_binary.empty()||adapter.container_binary.empty())throw std::invalid_argument("solver adapter descriptor requires identity and binaries");
    for(const auto& action:adapter.control_actions)if(action.empty())throw std::invalid_argument("solver control action must not be empty");
}

bool solver_has_capability(const SolverAdapterDescriptor& adapter,SolverCapability capability){return std::find(adapter.capabilities.begin(),adapter.capabilities.end(),capability)!=adapter.capabilities.end();}

void update_case_status(std::vector<CampaignRegistryEntry>& registry,CampaignRegistryEntry entry){
    if(entry.case_id.empty())throw std::invalid_argument("campaign registry entry needs case_id");
    auto it=std::find_if(registry.begin(),registry.end(),[&](const auto& item){return item.case_id==entry.case_id;});
    if(it==registry.end())registry.push_back(std::move(entry));else *it=std::move(entry);
}

CampaignSummary summarize_campaign_registry(std::span<const CampaignRegistryEntry> registry){
    CampaignSummary out;bool have_best=false;
    for(const auto& item:registry){
        switch(item.status){case CampaignCaseStatus::pending:++out.pending;break;case CampaignCaseStatus::running:++out.running;break;case CampaignCaseStatus::done:++out.done;break;case CampaignCaseStatus::failed:++out.failed;break;case CampaignCaseStatus::stopped:++out.stopped;break;}
        if(item.status==CampaignCaseStatus::done&&(!have_best||item.objective<out.best_objective)){have_best=true;out.best_objective=item.objective;out.best_case_id=item.case_id;}
    }
    return out;
}


namespace {
std::string escape_tsv(std::string_view value){
    std::string out;out.reserve(value.size());
    for(char ch:value){
        if(ch=='\\')out+="\\\\";
        else if(ch=='\t')out+="\\t";
        else if(ch=='\n')out+="\\n";
        else out.push_back(ch);
    }
    return out;
}
std::string unescape_tsv(std::string_view value){
    std::string out;out.reserve(value.size());
    for(std::size_t i=0;i<value.size();++i){
        if(value[i]=='\\'&&i+1U<value.size()){
            const char n=value[++i];
            if(n=='t')out.push_back('\t');
            else if(n=='n')out.push_back('\n');
            else out.push_back(n);
        } else out.push_back(value[i]);
    }
    return out;
}
std::vector<std::string_view> split_tabs(std::string_view line){
    std::vector<std::string_view> out;std::size_t begin=0U;
    while(begin<=line.size()){
        const std::size_t pos=line.find('\t',begin);
        if(pos==std::string_view::npos){out.push_back(line.substr(begin));break;}
        out.push_back(line.substr(begin,pos-begin));begin=pos+1U;
    }
    return out;
}
void ensure_safe_relative_path(const std::filesystem::path& relative){
    if(relative.empty()||relative.is_absolute())throw std::invalid_argument("campaign template path must be relative");
    for(const auto& part:relative){
        const auto token=part.string();
        if(token=="..")throw std::invalid_argument("campaign template path must not escape case directory");
    }
}
std::string upper_copy(std::string_view text){
    std::string out(text);for(char& ch:out)ch=static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));return out;
}
std::string join_display(std::span<const std::string> argv){
    std::ostringstream out;
    for(std::size_t i=0;i<argv.size();++i){
        if(i)out<<' ';
        const bool quote=argv[i].find_first_of(" \t\n\"'")!=std::string::npos;
        if(!quote)out<<argv[i];
        else {out<<'\'';for(char ch:argv[i]){if(ch=='\'')out<<"'\\''";else out<<ch;}out<<'\'';}
    }
    return out.str();
}
bool parse_size(std::string_view token,std::size_t& out){
    try{std::size_t pos=0U;const auto s=std::string(token);const unsigned long long v=std::stoull(s,&pos);if(pos!=s.size())return false;out=static_cast<std::size_t>(v);return true;}catch(...){return false;}
}
bool parse_double(std::string_view token,double& out){
    try{std::size_t pos=0U;const auto s=std::string(token);const double v=std::stod(s,&pos);if(pos!=s.size())return false;out=v;return std::isfinite(v);}catch(...){return false;}
}
std::vector<std::string> normalized_tokens(std::string line){
    for(char& ch:line)if(ch=='='||ch==','||ch==';'||ch==':')ch=' ';
    std::stringstream ss(line);std::vector<std::string> tokens;std::string token;while(ss>>token)tokens.push_back(token);return tokens;
}
bool has_path_separator(std::string_view value){return value.find('/')!=std::string_view::npos||value.find('\\')!=std::string_view::npos;}
std::string lower_copy(std::string_view text){
    std::string out(text);for(char& ch:out)ch=static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));return out;
}
bool is_executable_file(const std::filesystem::path& path){
#if defined(_WIN32)
    std::error_code ec;return std::filesystem::exists(path,ec)&&std::filesystem::is_regular_file(path,ec);
#else
    return ::access(path.c_str(),X_OK)==0;
#endif
}
std::filesystem::path find_executable_on_path(const std::string& name){
    if(name.empty())return {};
    if(has_path_separator(name)){return is_executable_file(name)?std::filesystem::path(name):std::filesystem::path{};}
    const char* env=std::getenv("PATH");if(!env)return {};
    std::stringstream paths(env);std::string dir;
    while(std::getline(paths,dir,':')){
        if(dir.empty())dir=".";
        const auto candidate=std::filesystem::path(dir)/name;
        if(is_executable_file(candidate))return candidate;
    }
    return {};
}
std::string read_text_limited(const std::filesystem::path& path,std::size_t max_bytes){
    std::ifstream in(path,std::ios::binary);if(!in)return {};
    std::string out;out.resize(max_bytes);in.read(out.data(),static_cast<std::streamsize>(out.size()));out.resize(static_cast<std::size_t>(in.gcount()));return out;
}
std::string read_text_full(const std::filesystem::path& path){
    std::ifstream in(path,std::ios::binary);if(!in)return {};std::ostringstream out;out<<in.rdbuf();return out.str();
}
void append_unique_path(std::vector<std::filesystem::path>& paths,std::filesystem::path path){
    if(std::find(paths.begin(),paths.end(),path)==paths.end())paths.push_back(std::move(path));
}
double objective_from_residuals(std::span<const ResidualSample> residuals,int exit_code){
    if(!residuals.empty())return residuals.back().residual;
    return exit_code==0?0.0:static_cast<double>(std::max(1,exit_code));
}
std::size_t iterations_from_residuals(std::span<const ResidualSample> residuals){
    std::size_t out=0U;for(const auto& sample:residuals)out=std::max(out,sample.iteration);return out;
}

}

std::string campaign_status_name(CampaignCaseStatus status){
    switch(status){
        case CampaignCaseStatus::pending:return "pending";
        case CampaignCaseStatus::running:return "running";
        case CampaignCaseStatus::done:return "done";
        case CampaignCaseStatus::failed:return "failed";
        case CampaignCaseStatus::stopped:return "stopped";
    }
    return "pending";
}

CampaignCaseStatus campaign_status_from_name(std::string_view name){
    const auto up=upper_copy(name);
    if(up=="PENDING")return CampaignCaseStatus::pending;
    if(up=="RUNNING")return CampaignCaseStatus::running;
    if(up=="DONE"||up=="SUCCESS"||up=="COMPLETED")return CampaignCaseStatus::done;
    if(up=="FAILED"||up=="ERROR")return CampaignCaseStatus::failed;
    if(up=="STOPPED"||up=="CANCELLED")return CampaignCaseStatus::stopped;
    throw std::invalid_argument("unknown campaign case status");
}

void save_campaign_registry(const std::filesystem::path& path,std::span<const CampaignRegistryEntry> registry){
    if(path.empty())throw std::invalid_argument("campaign registry path is empty");
    if(path.has_parent_path())std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("failed to write campaign registry");
    out<<"case_id\tstatus\tobjective\titerations\tmessage\n";
    out<<std::setprecision(17);
    for(const auto& entry:registry){
        out<<escape_tsv(entry.case_id)<<'\t'<<campaign_status_name(entry.status)<<'\t'<<entry.objective<<'\t'<<entry.iterations<<'\t'<<escape_tsv(entry.message)<<'\n';
    }
}

std::vector<CampaignRegistryEntry> load_campaign_registry(const std::filesystem::path& path){
    std::ifstream in(path);if(!in)throw std::runtime_error("failed to read campaign registry");
    std::vector<CampaignRegistryEntry> registry;std::string line;bool header=true;
    while(std::getline(in,line)){
        if(line.empty())continue;
        if(header){header=false;if(line.rfind("case_id\t",0U)==0U)continue;}
        const auto fields=split_tabs(line);if(fields.size()<5U)throw std::runtime_error("malformed campaign registry row");
        CampaignRegistryEntry entry;entry.case_id=unescape_tsv(fields[0]);entry.status=campaign_status_from_name(fields[1]);
        double objective{};std::size_t iterations{};
        if(!parse_double(fields[2],objective)||!parse_size(fields[3],iterations))throw std::runtime_error("malformed numeric campaign registry field");
        entry.objective=objective;entry.iterations=iterations;entry.message=unescape_tsv(fields[4]);registry.push_back(std::move(entry));
    }
    return registry;
}

CampaignPersistenceResult write_campaign_case_folders(const std::filesystem::path& root,std::span<const CampaignCase> cases,std::span<const CampaignTemplateFile> templates,const CampaignWriteOptions& options){
    if(root.empty()||cases.empty()||templates.empty())throw std::invalid_argument("campaign folder writer requires root, cases and templates");
    std::filesystem::create_directories(root);
    std::vector<CampaignRegistryEntry> registry;registry.reserve(cases.size());
    CampaignPersistenceResult result;result.root=root;
    for(const auto& item:templates)ensure_safe_relative_path(item.relative_path);
    for(const auto& row:cases){
        if(row.case_id.empty())throw std::invalid_argument("campaign case requires case_id");
        const auto case_dir=root/row.case_id;
        if(std::filesystem::exists(case_dir)&&!options.overwrite)throw std::runtime_error("campaign case directory already exists: "+case_dir.string());
        std::filesystem::create_directories(case_dir);
        for(const auto& templ:templates){
            const auto target=case_dir/templ.relative_path;std::filesystem::create_directories(target.parent_path());
            std::ofstream out(target,std::ios::trunc);if(!out)throw std::runtime_error("failed to write campaign template file");
            out<<render_campaign_template(templ.contents,row,true);++result.files_written;
        }
        std::vector<std::string> keys;keys.reserve(row.values.size());for(const auto& [key,value]:row.values)keys.push_back(key);std::sort(keys.begin(),keys.end());
        std::ofstream row_out(case_dir/options.doe_row_filename,std::ios::trunc);if(!row_out)throw std::runtime_error("failed to write campaign doe row");
        row_out<<"case_id";for(const auto& key:keys)row_out<<','<<key;row_out<<'\n'<<row.case_id;for(const auto& key:keys)row_out<<','<<row.values.at(key);row_out<<'\n';++result.files_written;
        registry.push_back({row.case_id,CampaignCaseStatus::pending,0.0,0U,"prepared"});++result.cases_written;
    }
    save_campaign_registry(root/options.registry_filename,registry);result.registry_entries=registry.size();
    std::ofstream manifest(root/options.manifest_filename,std::ios::trunc);if(!manifest)throw std::runtime_error("failed to write campaign manifest");
    manifest<<"relative_path\n";for(const auto& templ:templates)manifest<<templ.relative_path.generic_string()<<'\n';
    return result;
}

SolverCommandPlan build_solver_command_plan(const SolverRunRequest& request){
    validate_solver_adapter_descriptor(request.adapter);
    if(request.case_directory.empty()||request.mpi_ranks==0U||request.threads==0U)throw std::invalid_argument("invalid solver run request");
    SolverCommandPlan plan;plan.runtime=request.runtime;plan.environment={{"CFD_SOLVER_ADAPTER",request.adapter.name},{"OMP_NUM_THREADS",std::to_string(request.threads)}};
    const auto case_path=request.case_directory.string();
    auto append_common=[&](std::vector<std::string>& argv,const std::string& case_arg){
        argv.push_back("--case");argv.push_back(case_arg);argv.push_back("--nprocs");argv.push_back(std::to_string(request.mpi_ranks));argv.push_back("--threads");argv.push_back(std::to_string(request.threads));
        argv.insert(argv.end(),request.extra_args.begin(),request.extra_args.end());
    };
    switch(request.runtime){
        case SolverRuntime::native:
            plan.argv.push_back(request.adapter.native_binary);append_common(plan.argv,case_path);break;
        case SolverRuntime::docker:{
            if(request.docker_image.empty())throw std::invalid_argument("docker solver run requires image");
            plan.argv={"docker","run","--rm","-v",case_path+":/case","-w","/case",request.docker_image,request.adapter.container_binary};append_common(plan.argv,"/case");break;
        }
        case SolverRuntime::singularity:{
            if(request.singularity_image.empty())throw std::invalid_argument("singularity solver run requires image");
            plan.argv={"apptainer","exec","--bind",case_path+":/case",request.singularity_image,request.adapter.container_binary};append_common(plan.argv,"/case");break;
        }
        case SolverRuntime::slurm:{
            plan.requires_scheduler=true;std::vector<std::string> inner{request.adapter.native_binary};append_common(inner,case_path);
            plan.argv={"sbatch","--ntasks",std::to_string(request.mpi_ranks),"--cpus-per-task",std::to_string(request.threads)};
            if(!request.slurm_partition.empty()){plan.argv.push_back("--partition");plan.argv.push_back(request.slurm_partition);}plan.argv.push_back("--wrap");plan.argv.push_back(join_display(inner));break;
        }
    }
    plan.display_command=join_display(plan.argv);return plan;
}

RuntimeDoctorReport doctor_solver_runtime(const SolverRunRequest& request){
    RuntimeDoctorReport report;auto add=[&](std::string name,bool ok,std::string message){report.checks.push_back({std::move(name),ok,std::move(message)});};
    try{
        validate_solver_adapter_descriptor(request.adapter);add("adapter",true,"solver adapter descriptor is valid");
    }catch(const std::exception& e){add("adapter",false,e.what());report.ok=false;return report;}
    std::error_code ec;
    const bool case_ok=!request.case_directory.empty()&&std::filesystem::exists(request.case_directory,ec)&&std::filesystem::is_directory(request.case_directory,ec);
    add("case_directory",case_ok,case_ok?"case directory exists":"case directory is missing or not a directory");
    try{
        auto plan=build_solver_command_plan(request);
        const std::string launcher=plan.argv.empty()?std::string{}:plan.argv.front();
        const auto found=find_executable_on_path(launcher);
        add("launcher",!found.empty(),!found.empty()?"found executable: "+found.string():"executable not found: "+launcher);
        if(request.runtime==SolverRuntime::docker)add("docker_image",!request.docker_image.empty(),request.docker_image.empty()?"docker image not configured":"docker image configured");
        if(request.runtime==SolverRuntime::singularity)add("singularity_image",!request.singularity_image.empty(),request.singularity_image.empty()?"singularity image not configured":"singularity image configured");
        if(request.runtime==SolverRuntime::slurm)add("scheduler",plan.requires_scheduler,"slurm command plan requires scheduler submission");
    }catch(const std::exception& e){add("command_plan",false,e.what());}
    report.ok=std::all_of(report.checks.begin(),report.checks.end(),[](const auto& check){return check.ok;});return report;
}

CampaignOutputDiscovery discover_campaign_outputs(const std::filesystem::path& case_directory){
    CampaignOutputDiscovery out;std::error_code ec;
    if(case_directory.empty()||!std::filesystem::exists(case_directory,ec))return out;
    for(std::filesystem::recursive_directory_iterator it(case_directory,ec),end;it!=end&&!ec;it.increment(ec)){
        if(ec||!it->is_regular_file(ec))continue;
        const auto path=it->path();const auto name=lower_copy(path.filename().string());
        const auto ext=lower_copy(path.extension().string());
        const bool text_like=ext==".log"||ext==".out"||ext==".txt"||ext==".dat"||ext==".csv";
        if(!text_like)continue;
        if(ext==".log"||ext==".out"||name.find("stdout")!=std::string::npos||name.find("stderr")!=std::string::npos)append_unique_path(out.log_files,path);
        if(name.find("residual")!=std::string::npos||name.find("history")!=std::string::npos)append_unique_path(out.residual_files,path);
        if(name.find("perf")!=std::string::npos||name.find("timing")!=std::string::npos||name.find("time")!=std::string::npos)append_unique_path(out.performance_files,path);
    }
    std::sort(out.log_files.begin(),out.log_files.end());std::sort(out.residual_files.begin(),out.residual_files.end());std::sort(out.performance_files.begin(),out.performance_files.end());
    return out;
}

LocalSolverRunResult run_local_solver_case(const SolverRunRequest& request,const LocalSolverRunOptions& options){
    auto plan=build_solver_command_plan(request);
    if(plan.requires_scheduler)throw std::invalid_argument("local campaign runner does not submit scheduler jobs");
    auto doctor=doctor_solver_runtime(request);
    if(!doctor.ok){
        std::ostringstream msg;msg<<"runtime doctor failed";for(const auto& check:doctor.checks)if(!check.ok)msg<<"; "<<check.name<<": "<<check.message;throw std::runtime_error(msg.str());
    }
    LocalSolverRunResult result;result.case_id=request.case_directory.filename().string();
    const auto stdout_path=request.case_directory/options.stdout_filename;const auto stderr_path=request.case_directory/options.stderr_filename;
#if defined(_WIN32)
    throw std::runtime_error("local process launcher is not implemented on this platform");
#else
    const int out_flags=O_CREAT|O_TRUNC|O_WRONLY;
    int out_fd=::open(stdout_path.c_str(),out_flags,0644);if(out_fd<0)throw std::runtime_error("failed to open solver stdout log");
    int err_fd=::open(stderr_path.c_str(),out_flags,0644);if(err_fd<0){::close(out_fd);throw std::runtime_error("failed to open solver stderr log");}
    const pid_t pid=::fork();
    if(pid<0){::close(out_fd);::close(err_fd);throw std::runtime_error("failed to fork solver process");}
    if(pid==0){
        ::dup2(out_fd,STDOUT_FILENO);::dup2(err_fd,STDERR_FILENO);::close(out_fd);::close(err_fd);
        ::chdir(request.case_directory.c_str());
        for(const auto& [key,value]:plan.environment)::setenv(key.c_str(),value.c_str(),1);
        std::vector<std::string> argv_storage=plan.argv;std::vector<char*> argv;argv.reserve(argv_storage.size()+1U);
        for(auto& arg:argv_storage)argv.push_back(arg.data());argv.push_back(nullptr);
        ::execvp(argv[0],argv.data());
        std::fprintf(stderr,"execvp failed for %s\n",argv[0]);
        ::_exit(127);
    }
    ::close(out_fd);::close(err_fd);
    int status=0;while(::waitpid(pid,&status,0)<0){if(errno!=EINTR)break;}
    if(WIFEXITED(status))result.exit_code=WEXITSTATUS(status);else if(WIFSIGNALED(status))result.exit_code=128+WTERMSIG(status);else result.exit_code=127;
#endif
    result.stdout_text=read_text_limited(stdout_path,options.max_capture_bytes);
    result.stderr_text=read_text_limited(stderr_path,options.max_capture_bytes);
    result.discovered_outputs=discover_campaign_outputs(request.case_directory);
    const std::string combined=result.stdout_text+"\n"+result.stderr_text;
    auto parsed_res=parse_residual_history(combined);result.residuals.insert(result.residuals.end(),parsed_res.begin(),parsed_res.end());
    auto parsed_perf=parse_performance_metrics(combined);result.performance.insert(result.performance.end(),parsed_perf.begin(),parsed_perf.end());
    for(const auto& file:result.discovered_outputs.residual_files){auto v=parse_residual_history(read_text_full(file));result.residuals.insert(result.residuals.end(),v.begin(),v.end());}
    for(const auto& file:result.discovered_outputs.performance_files){auto v=parse_performance_metrics(read_text_full(file));result.performance.insert(result.performance.end(),v.begin(),v.end());}
    std::sort(result.residuals.begin(),result.residuals.end(),[](const auto& a,const auto& b){return std::tie(a.iteration,a.equation,a.residual)<std::tie(b.iteration,b.equation,b.residual);});
    const auto outcome=detect_solver_outcome_from_log(combined);
    if(result.exit_code==0)result.status=(outcome==CampaignCaseStatus::running?CampaignCaseStatus::done:outcome);
    else result.status=(outcome==CampaignCaseStatus::stopped?CampaignCaseStatus::stopped:CampaignCaseStatus::failed);
    if(!options.write_launcher_logs){std::filesystem::remove(stdout_path);std::filesystem::remove(stderr_path);}
    return result;
}

LocalCampaignRunResult run_local_campaign(const std::filesystem::path& root,const LocalCampaignRunOptions& options,const std::string& registry_filename){
    if(root.empty()||registry_filename.empty())throw std::invalid_argument("local campaign runner requires a root and registry filename");
    validate_solver_adapter_descriptor(options.adapter);
    auto registry=load_campaign_registry(root/registry_filename);
    LocalCampaignRunResult result;
    for(auto& entry:registry){
        if(options.pending_only&&entry.status!=CampaignCaseStatus::pending)continue;
        SolverRunRequest request;request.adapter=options.adapter;request.runtime=options.runtime;request.mpi_ranks=options.mpi_ranks;request.threads=options.threads;request.extra_args=options.extra_args;request.case_directory=root/entry.case_id;
        entry.status=CampaignCaseStatus::running;entry.message="launched";save_campaign_registry(root/registry_filename,registry);
        try{
            auto run=run_local_solver_case(request,options.solver_run);++result.launched;
            entry.status=run.status;entry.objective=objective_from_residuals(run.residuals,run.exit_code);entry.iterations=iterations_from_residuals(run.residuals);entry.message="exit="+std::to_string(run.exit_code)+" residuals="+std::to_string(run.residuals.size())+" perf="+std::to_string(run.performance.size());
            result.case_results.push_back(std::move(run));
        }catch(const std::exception& e){
            entry.status=CampaignCaseStatus::failed;entry.objective=std::numeric_limits<double>::infinity();entry.message=e.what();
        }
        save_campaign_registry(root/registry_filename,registry);
    }
    result.summary=summarize_campaign_registry(registry);return result;
}


CampaignCaseStatus detect_solver_outcome_from_log(std::string_view log){
    const auto text=upper_copy(log);
    if(text.find("FAILED")!=std::string::npos||text.find("FATAL")!=std::string::npos||text.find("ABORT")!=std::string::npos||text.find("ERROR")!=std::string::npos)return CampaignCaseStatus::failed;
    if(text.find("STOPPED")!=std::string::npos||text.find("CANCELLED")!=std::string::npos)return CampaignCaseStatus::stopped;
    if(text.find("DONE")!=std::string::npos||text.find("SUCCESS")!=std::string::npos||text.find("COMPLETED")!=std::string::npos)return CampaignCaseStatus::done;
    return CampaignCaseStatus::running;
}

std::vector<ResidualSample> parse_residual_history(std::string_view log){
    std::vector<ResidualSample> out;std::stringstream lines{std::string(log)};std::string line;
    while(std::getline(lines,line)){
        auto tokens=normalized_tokens(line);if(tokens.empty())continue;ResidualSample sample;bool have_iter=false,have_res=false;
        if(tokens.size()>=2U&&parse_size(tokens[0],sample.iteration)&&parse_double(tokens[1],sample.residual)){sample.equation="residual";out.push_back(sample);continue;}
        for(std::size_t i=0;i+1U<tokens.size();++i){const auto key=upper_copy(tokens[i]);
            if((key=="ITER"||key=="ITERATION")&&parse_size(tokens[i+1U],sample.iteration))have_iter=true;
            else if((key=="RES"||key=="RESIDUAL")&&parse_double(tokens[i+1U],sample.residual))have_res=true;
            else if(key=="EQUATION"||key=="FIELD")sample.equation=tokens[i+1U];
        }
        if(have_iter&&have_res)out.push_back(sample);
    }
    return out;
}

std::vector<PerformanceMetric> parse_performance_metrics(std::string_view log){
    std::vector<PerformanceMetric> out;std::stringstream lines{std::string(log)};std::string line;
    while(std::getline(lines,line)){
        auto tokens=normalized_tokens(line);if(tokens.size()<2U)continue;PerformanceMetric metric;
        std::size_t value_index=1U;std::string key=tokens[0];
        if(key.rfind("perf.",0U)==0U)key=key.substr(5U);
        else if(key.rfind("PERF.",0U)==0U)key=key.substr(5U);
        else if(upper_copy(key)!="WALL_TIME"&&upper_copy(key)!="CPU_TIME"&&upper_copy(key)!="CELLS_PER_S"&&upper_copy(key)!="ITERATIONS_PER_S")continue;
        double value{};if(!parse_double(tokens[value_index],value))continue;
        metric.key=key;metric.value=value;if(tokens.size()>value_index+1U)metric.units=tokens[value_index+1U];out.push_back(metric);
    }
    return out;
}

CampaignOptimizationResult run_gradient_descent_campaign(const std::function<double(std::span<const double>)>& objective,std::vector<double> initial,const CampaignOptimizationConfig& config){
    if(!objective||initial.empty()||config.max_iterations==0U||!(config.step_size>0.0)||!(config.tolerance>0.0))throw std::invalid_argument("invalid campaign optimization request");
    CampaignOptimizationResult result;result.parameters=std::move(initial);result.objective=objective(result.parameters);
    double step=config.step_size;
    for(std::size_t iter=0;iter<config.max_iterations;++iter){
        result.last_gradient=finite_difference_gradient(objective,result.parameters,config.finite_difference);
        double norm2=0.0;for(double g:result.last_gradient)norm2+=g*g;const double norm=std::sqrt(norm2);
        result.iterations=iter+1U;if(norm<config.tolerance){result.converged=true;break;}
        bool accepted=false;
        for(int backtrack=0;backtrack<12;++backtrack){
            auto trial=gradient_descent_update(result.parameters,result.last_gradient,step);const double value=objective(trial);
            if(value<result.objective){result.parameters=std::move(trial);result.objective=value;accepted=true;break;}
            step*=0.5;
        }
        if(!accepted||step<config.tolerance){result.converged=true;break;}
    }
    return result;
}

std::vector<double> finite_difference_gradient(const std::function<double(std::span<const double>)>& objective,std::span<const double> parameters,const FiniteDifferenceGradientConfig& config){
    if(!objective||parameters.empty()||!(config.absolute_step>0.0)||config.relative_step<0.0)throw std::invalid_argument("invalid finite difference gradient request");
    std::vector<double> x(parameters.begin(),parameters.end()),gradient(x.size());
    const double base=objective(x);
    for(std::size_t i=0;i<x.size();++i){const double h=std::max(config.absolute_step,config.relative_step*std::abs(x[i]));auto xp=x;xp[i]+=h;if(config.central){auto xm=x;xm[i]-=h;gradient[i]=(objective(xp)-objective(xm))/(2.0*h);}else gradient[i]=(objective(xp)-base)/h;}
    return gradient;
}

std::vector<double> gradient_descent_update(std::span<const double> parameters,std::span<const double> gradient,double step_size){
    if(parameters.size()!=gradient.size()||parameters.empty()||!(step_size>=0.0))throw std::invalid_argument("invalid gradient descent update");
    std::vector<double> out(parameters.begin(),parameters.end());
    for(std::size_t i=0;i<out.size();++i)out[i]-=step_size*gradient[i];
    return out;
}

} // namespace cfd::workflow
