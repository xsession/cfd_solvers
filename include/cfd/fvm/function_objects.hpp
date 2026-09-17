#pragma once
#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
namespace cfd::fvm {
struct SolverReport { std::size_t step{}; double time{}; std::unordered_map<std::string,double> residuals; };
using FunctionObject=std::function<void(const SolverReport&)>;
class FunctionObjectManager {
public:
    void add(FunctionObject object);
    void execute(const SolverReport& report) const;
    [[nodiscard]] std::size_t size() const noexcept { return objects_.size(); }
private: std::vector<FunctionObject> objects_;
};
class ResidualHistory {
public:
    void observe(const SolverReport& report);
    [[nodiscard]] const std::vector<SolverReport>& records() const noexcept { return records_; }
    [[nodiscard]] double latest(std::string_view field) const;
private:std::vector<SolverReport> records_;
};
} // namespace cfd::fvm
