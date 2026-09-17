#pragma once
#include "cfd/solvers/optics/sequential.hpp"
#include <string>
namespace cfd::optics {
[[nodiscard]] std::string sequential_system_to_json(const SequentialOpticalSystem& system,double object_space_index=1.0);
[[nodiscard]] SequentialOpticalSystem sequential_system_from_json(const std::string& json);
} // namespace cfd::optics
