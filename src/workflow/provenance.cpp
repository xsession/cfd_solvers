#include "cfd/workflow/provenance.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
namespace cfd::workflow {namespace {std::string esc(std::string_view s){std::string o;for(char c:s){if(c=='"'||c=='\\')o+='\\';o+=c;}return o;}}
ResultProvenance make_provenance(std::string solver,std::string version,std::string commit){auto now=std::chrono::system_clock::now();auto tt=std::chrono::system_clock::to_time_t(now);std::tm tm{};
#if defined(_WIN32)
gmtime_s(&tm,&tt);
#else
gmtime_r(&tt,&tm);
#endif
std::ostringstream t;t<<std::put_time(&tm,"%Y-%m-%dT%H:%M:%SZ");return {std::move(solver),std::move(version),std::move(commit),t.str(),{}};}
std::string provenance_json(const ResultProvenance&p){std::ostringstream o;o<<"{\"solver\":\""<<esc(p.solver)<<"\",\"version\":\""<<esc(p.version)<<"\",\"git_commit\":\""<<esc(p.git_commit)<<"\",\"timestamp_utc\":\""<<p.timestamp_utc<<"\",\"metadata\":{";bool first=true;for(auto&kv:p.metadata){if(!first)o<<',';first=false;o<<'"'<<esc(kv.first)<<"\":\""<<esc(kv.second)<<'"';}o<<"}}";return o.str();}
void write_provenance(const std::filesystem::path&path,const ResultProvenance&p){std::ofstream o(path);if(!o)throw std::runtime_error("cannot write provenance");o<<provenance_json(p)<<'\n';}
} // namespace cfd::workflow
