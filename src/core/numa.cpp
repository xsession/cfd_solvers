#include "cfd/core/numa.hpp"

#include <algorithm>
#include <atomic>
#include <charconv>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#elif defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

namespace cfd::core {
namespace {

[[nodiscard]] unsigned parse_uint(std::string_view token) {
    if (token.empty()) throw std::invalid_argument("empty CPU index");
    unsigned value = 0U;
    const auto* first = token.data();
    const auto* last = token.data() + token.size();
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc{} || ptr != last) throw std::invalid_argument("invalid CPU list token");
    return value;
}

[[nodiscard]] std::vector<unsigned> allowed_cpus() {
#if defined(__linux__)
    long configured = sysconf(_SC_NPROCESSORS_CONF);
    std::size_t capacity = configured > 0 ? static_cast<std::size_t>(configured) : 1024U;
    capacity = std::max<std::size_t>(capacity, 1024U);
    for (int attempt = 0; attempt < 8; ++attempt) {
        cpu_set_t* set = CPU_ALLOC(capacity);
        if (!set) break;
        const std::size_t bytes = CPU_ALLOC_SIZE(capacity);
        CPU_ZERO_S(bytes, set);
        errno = 0;
        const int rc = sched_getaffinity(0, bytes, set);
        if (rc == 0) {
            std::vector<unsigned> cpus;
            for (std::size_t cpu = 0U; cpu < capacity; ++cpu) {
                if (CPU_ISSET_S(cpu, bytes, set)) cpus.push_back(static_cast<unsigned>(cpu));
            }
            CPU_FREE(set);
            if (!cpus.empty()) return cpus;
            break;
        }
        const int error = errno;
        CPU_FREE(set);
        if (error != EINVAL || capacity > (1U << 20U)) break;
        capacity *= 2U;
    }
#endif
    const unsigned count = std::max(1U, std::thread::hardware_concurrency());
    std::vector<unsigned> cpus(count);
    for (unsigned i = 0U; i < count; ++i) cpus[i] = i;
    return cpus;
}

[[nodiscard]] bool contains_cpu(std::span<const unsigned> cpus, unsigned cpu) {
    return std::find(cpus.begin(), cpus.end(), cpu) != cpus.end();
}

} // namespace

std::vector<unsigned> parse_cpu_list(std::string_view text) {
    std::vector<unsigned> cpus;
    std::size_t cursor = 0U;
    while (cursor < text.size()) {
        while (cursor < text.size() && (text[cursor] == ' ' || text[cursor] == '\t' || text[cursor] == '\n')) ++cursor;
        if (cursor >= text.size()) break;
        const std::size_t comma = text.find(',', cursor);
        const std::size_t end = comma == std::string_view::npos ? text.size() : comma;
        auto token = text.substr(cursor, end - cursor);
        while (!token.empty() && (token.back() == ' ' || token.back() == '\t' || token.back() == '\n')) token.remove_suffix(1U);
        const std::size_t dash = token.find('-');
        if (dash == std::string_view::npos) {
            cpus.push_back(parse_uint(token));
        } else {
            const unsigned first = parse_uint(token.substr(0U, dash));
            const unsigned last = parse_uint(token.substr(dash + 1U));
            if (last < first) throw std::invalid_argument("descending CPU range");
            for (unsigned cpu = first;; ++cpu) {
                cpus.push_back(cpu);
                if (cpu == last) break;
                if (cpu == std::numeric_limits<unsigned>::max()) throw std::overflow_error("CPU range overflow");
            }
        }
        if (comma == std::string_view::npos) break;
        cursor = comma + 1U;
    }
    std::sort(cpus.begin(), cpus.end());
    cpus.erase(std::unique(cpus.begin(), cpus.end()), cpus.end());
    return cpus;
}

std::vector<NumaNodeInfo> discover_numa_nodes() {
    const auto allowed = allowed_cpus();
    std::vector<NumaNodeInfo> nodes;
#if defined(__linux__)
    const std::filesystem::path root{"/sys/devices/system/node"};
    std::error_code ec;
    if (std::filesystem::is_directory(root, ec)) {
        for (const auto& entry : std::filesystem::directory_iterator(root, ec)) {
            if (ec) break;
            const auto name = entry.path().filename().string();
            if (name.size() <= 4U || name.rfind("node", 0U) != 0U) continue;
            int node_id = -1;
            const auto suffix = std::string_view{name}.substr(4U);
            const auto [ptr, parse_ec] = std::from_chars(suffix.data(), suffix.data() + suffix.size(), node_id);
            if (parse_ec != std::errc{} || ptr != suffix.data() + suffix.size() || node_id < 0) continue;
            std::ifstream in(entry.path() / "cpulist");
            if (!in) continue;
            std::string line;
            std::getline(in, line);
            auto cpus = parse_cpu_list(line);
            cpus.erase(std::remove_if(cpus.begin(), cpus.end(), [&](unsigned cpu) {
                return !contains_cpu(allowed, cpu);
            }), cpus.end());
            if (!cpus.empty()) nodes.push_back({node_id, std::move(cpus)});
        }
    }
#endif
    if (nodes.empty()) nodes.push_back({0, allowed});
    std::sort(nodes.begin(), nodes.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
    return nodes;
}

NumaPlacementPlan plan_numa_thread_placement(std::size_t thread_count,
                                             NumaPlacementPolicy policy,
                                             std::optional<int> node_id) {
    NumaPlacementPlan plan;
    plan.policy = policy;
    if (policy == NumaPlacementPolicy::inherit) return plan;

    auto nodes = discover_numa_nodes();
    std::size_t available = 0U;
    for (const auto& node : nodes) available += node.cpus.size();
    if (available == 0U) throw std::runtime_error("no CPUs available for NUMA placement");
    if (thread_count == 0U) thread_count = available;

    if (policy == NumaPlacementPolicy::explicit_node) {
        if (!node_id) throw std::invalid_argument("explicit NUMA placement requires a node id");
        const auto it = std::find_if(nodes.begin(), nodes.end(), [&](const auto& node) { return node.id == *node_id; });
        if (it == nodes.end() || it->cpus.empty()) throw std::invalid_argument("requested NUMA node is not available");
        plan.thread_cpus.reserve(thread_count);
        plan.thread_nodes.reserve(thread_count);
        for (std::size_t i = 0U; i < thread_count; ++i) {
            plan.thread_cpus.push_back(it->cpus[i % it->cpus.size()]);
            plan.thread_nodes.push_back(it->id);
        }
        return plan;
    }

    if (policy == NumaPlacementPolicy::compact) {
        std::vector<std::pair<unsigned, int>> flat;
        flat.reserve(available);
        for (const auto& node : nodes) for (unsigned cpu : node.cpus) flat.emplace_back(cpu, node.id);
        plan.thread_cpus.reserve(thread_count);
        plan.thread_nodes.reserve(thread_count);
        for (std::size_t i = 0U; i < thread_count; ++i) {
            const auto& [cpu, node] = flat[i % flat.size()];
            plan.thread_cpus.push_back(cpu);
            plan.thread_nodes.push_back(node);
        }
        return plan;
    }

    if (policy != NumaPlacementPolicy::spread) throw std::invalid_argument("unknown NUMA placement policy");
    plan.thread_cpus.reserve(thread_count);
    plan.thread_nodes.reserve(thread_count);
    std::vector<std::size_t> next(nodes.size(), 0U);
    for (std::size_t i = 0U; i < thread_count; ++i) {
        const std::size_t ni = i % nodes.size();
        const auto& node = nodes[ni];
        const std::size_t ci = next[ni]++ % node.cpus.size();
        plan.thread_cpus.push_back(node.cpus[ci]);
        plan.thread_nodes.push_back(node.id);
    }
    return plan;
}

bool bind_current_thread_to_cpu(unsigned cpu) noexcept {
#if defined(__linux__)
    const std::size_t capacity = std::max<std::size_t>(static_cast<std::size_t>(cpu) + 1U, 1024U);
    cpu_set_t* set = CPU_ALLOC(capacity);
    if (!set) return false;
    const std::size_t bytes = CPU_ALLOC_SIZE(capacity);
    CPU_ZERO_S(bytes, set);
    CPU_SET_S(static_cast<std::size_t>(cpu), bytes, set);
    const bool ok = sched_setaffinity(0, bytes, set) == 0;
    CPU_FREE(set);
    return ok;
#elif defined(_WIN32)
    constexpr unsigned bits = static_cast<unsigned>(sizeof(DWORD_PTR) * 8U);
    if (cpu >= bits) return false;
    const DWORD_PTR mask = static_cast<DWORD_PTR>(1) << cpu;
    return SetThreadAffinityMask(GetCurrentThread(), mask) != 0;
#else
    (void)cpu;
    return false;
#endif
}

bool apply_numa_thread_placement(const NumaPlacementPlan& plan) noexcept {
    if (!plan.binds_threads()) return true;
#if defined(CFD_HAS_OPENMP)
    const auto limited = std::min<std::size_t>(plan.thread_count(), static_cast<std::size_t>(std::numeric_limits<int>::max()));
    const int requested = static_cast<int>(limited);
    std::atomic<bool> ok{true};
    std::atomic<int> actual{0};
    #pragma omp parallel num_threads(requested)
    {
        actual.store(omp_get_num_threads(), std::memory_order_relaxed);
        const auto tid = static_cast<std::size_t>(omp_get_thread_num());
        if (tid >= plan.thread_cpus.size() || !bind_current_thread_to_cpu(plan.thread_cpus[tid])) {
            ok.store(false, std::memory_order_relaxed);
        }
    }
    return ok.load(std::memory_order_relaxed) && actual.load(std::memory_order_relaxed) == requested;
#else
    return plan.thread_count() == 1U && bind_current_thread_to_cpu(plan.thread_cpus.front());
#endif
}

} // namespace cfd::core
