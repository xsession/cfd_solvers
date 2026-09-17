#pragma once
#include <string>
#include <vector>
namespace cfd::workflow {struct ValidationCase {std::string name,family,command,criterion;};[[nodiscard]] const std::vector<ValidationCase>& validation_catalog();}
