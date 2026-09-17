#pragma once
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <vector>
namespace cfd::workflow {
struct SweepResult {double parameter{},objective{};};
[[nodiscard]] std::vector<SweepResult> parameter_sweep(std::span<const double> values,const std::function<double(double)>& objective);
struct InverseResult {std::vector<double> parameters;double objective{};std::size_t iterations{};bool converged{};};
[[nodiscard]] InverseResult inverse_coordinate_search(const std::function<double(std::span<const double>)>& objective,std::vector<double> initial,std::vector<double> step,std::size_t max_iterations=200,double tolerance=1e-8);
class RestartableWorkflow {
public:
    using Action=std::function<void()>;
    std::size_t add(std::string name,Action action,std::vector<std::size_t> dependencies={});
    void run(const std::filesystem::path& checkpoint);
    [[nodiscard]] const std::vector<std::string>& completed()const noexcept{return completed_names_;}
private:
    struct Node{std::string name;Action action;std::vector<std::size_t> deps;};std::vector<Node>nodes_;std::vector<std::string>completed_names_;
};
}
