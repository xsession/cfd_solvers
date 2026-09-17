#include "cfd/core/profiler.hpp"
#include <algorithm>
#include <limits>
namespace cfd::core {
Profiler& Profiler::global(){static Profiler p;return p;}
void Profiler::record(std::string_view name,double seconds){if(!enabled_)return;std::scoped_lock lock(mutex_);auto& c=counters_[std::string(name)];if(c.calls==0){c.name=std::string(name);c.minimum_seconds=seconds;c.maximum_seconds=seconds;}else{c.minimum_seconds=std::min(c.minimum_seconds,seconds);c.maximum_seconds=std::max(c.maximum_seconds,seconds);}++c.calls;c.total_seconds+=seconds;}
std::vector<ProfileCounter> Profiler::snapshot() const {std::scoped_lock lock(mutex_);std::vector<ProfileCounter> out;out.reserve(counters_.size());for(const auto& kv:counters_)out.push_back(kv.second);std::sort(out.begin(),out.end(),[](const auto&a,const auto&b){return a.name<b.name;});return out;}
void Profiler::reset(){std::scoped_lock lock(mutex_);counters_.clear();}
ScopedProfile::ScopedProfile(std::string_view name,Profiler& profiler):name_(name),profiler_(&profiler),start_(std::chrono::steady_clock::now()){}
ScopedProfile::~ScopedProfile(){if(profiler_)profiler_->record(name_,std::chrono::duration<double>(std::chrono::steady_clock::now()-start_).count());}
} // namespace cfd::core
