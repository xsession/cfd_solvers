#pragma once
#include <chrono>
#include <cstddef>
#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
namespace cfd::core {
struct ProfileCounter { std::string name; std::size_t calls{}; double total_seconds{}; double minimum_seconds{}; double maximum_seconds{}; };
class Profiler {
public:
    static Profiler& global();
    void set_enabled(bool enabled) noexcept { enabled_=enabled; }
    [[nodiscard]] bool enabled() const noexcept { return enabled_; }
    void record(std::string_view name,double seconds);
    [[nodiscard]] std::vector<ProfileCounter> snapshot() const;
    void reset();
private:
    mutable std::mutex mutex_; std::unordered_map<std::string,ProfileCounter> counters_; bool enabled_{true};
};
class ScopedProfile {
public:
    explicit ScopedProfile(std::string_view name,Profiler& profiler=Profiler::global());
    ~ScopedProfile();
    ScopedProfile(const ScopedProfile&)=delete; ScopedProfile& operator=(const ScopedProfile&)=delete;
private:
    std::string name_; Profiler* profiler_{}; std::chrono::steady_clock::time_point start_{};
};
} // namespace cfd::core
