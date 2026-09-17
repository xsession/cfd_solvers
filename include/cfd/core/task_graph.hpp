#pragma once
#include <cstddef>
#include <functional>
#include <string>
#include <vector>
namespace cfd::core {
class TaskGraph {
public:
    using TaskId=std::size_t;TaskId add(std::string name,std::function<void()> task,std::vector<TaskId> dependencies={});
    void run(bool asynchronous=true);[[nodiscard]] std::size_t size()const noexcept{return nodes_.size();}
private:struct Node{std::string name;std::function<void()> fn;std::vector<TaskId> deps;};std::vector<Node>nodes_;[[nodiscard]] std::vector<TaskId> topological_order()const;
};
} // namespace cfd::core
