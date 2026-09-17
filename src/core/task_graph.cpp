#include "cfd/core/task_graph.hpp"
#include <future>
#include <stdexcept>
namespace cfd::core {
TaskGraph::TaskId TaskGraph::add(std::string name,std::function<void()> task,std::vector<TaskId> deps){if(!task)throw std::invalid_argument("task graph requires callable");for(auto d:deps)if(d>=nodes_.size())throw std::out_of_range("task dependency does not exist");nodes_.push_back({std::move(name),std::move(task),std::move(deps)});return nodes_.size()-1;}
std::vector<TaskGraph::TaskId> TaskGraph::topological_order()const{std::vector<unsigned>state(nodes_.size());std::vector<TaskId>out;std::function<void(TaskId)>visit=[&](TaskId n){if(state[n]==1)throw std::runtime_error("task graph cycle");if(state[n]==2)return;state[n]=1;for(auto d:nodes_[n].deps)visit(d);state[n]=2;out.push_back(n);};for(TaskId i=0;i<nodes_.size();++i)visit(i);return out;}
void TaskGraph::run(bool async){const auto order=topological_order();if(!async){for(auto id:order)nodes_[id].fn();return;}std::vector<std::shared_future<void>>f(nodes_.size());for(auto id:order){auto deps=nodes_[id].deps;auto fn=nodes_[id].fn;f[id]=std::async(std::launch::async,[deps=std::move(deps),fn=std::move(fn),&f]{for(auto d:deps)f[d].get();fn();}).share();}for(auto id:order)f[id].get();}
} // namespace cfd::core
