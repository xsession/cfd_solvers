#pragma once

#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

#ifndef CFD_SOLVERS_VERSION
#define CFD_SOLVERS_VERSION "unknown"
#endif

namespace cfd::examples {

struct Options {
    bool quick{};
    std::size_t scale{1U};
};

inline Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        if (arg == "--quick") o.quick = true;
        else if (arg == "--scale" && i + 1 < argc) {
            o.scale = static_cast<std::size_t>(std::strtoull(argv[++i], nullptr, 10));
            if (o.scale == 0U) o.scale = 1U;
        }
    }
    return o;
}

class Timer {
public:
    Timer() : start_(clock::now()) {}
    [[nodiscard]] double milliseconds() const {
        return std::chrono::duration<double, std::milli>(clock::now() - start_).count();
    }
private:
    using clock = std::chrono::steady_clock;
    clock::time_point start_;
};

struct BenchmarkRecord {
    std::string feature_set;
    std::string case_name;
    std::string backend{"cpu"};
    std::string unit{"work_items"};
    std::size_t problem_size{};
    std::size_t steps{};
    double setup_ms{};
    double simulation_ms{};
    double work_units{};
    double checksum{};
};

inline void emit(const BenchmarkRecord& r) {
    const double seconds = r.simulation_ms * 1.0e-3;
    const double throughput = seconds > 0.0 ? r.work_units / seconds : 0.0;
    std::cout << std::setprecision(17)
              << "CFD_BENCH {"
              << "\"schema\":1,"
              << "\"version\":\"" << CFD_SOLVERS_VERSION << "\","
              << "\"feature_set\":\"" << r.feature_set << "\","
              << "\"case\":\"" << r.case_name << "\","
              << "\"backend\":\"" << r.backend << "\","
              << "\"problem_size\":" << r.problem_size << ','
              << "\"steps\":" << r.steps << ','
              << "\"setup_ms\":" << r.setup_ms << ','
              << "\"simulation_ms\":" << r.simulation_ms << ','
              << "\"total_ms\":" << (r.setup_ms + r.simulation_ms) << ','
              << "\"work_units\":" << r.work_units << ','
              << "\"work_unit_name\":\"" << r.unit << "\","
              << "\"throughput\":" << throughput << ','
              << "\"checksum\":" << r.checksum
              << "}\n";
}

} // namespace cfd::examples
