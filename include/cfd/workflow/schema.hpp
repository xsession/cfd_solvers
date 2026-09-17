#pragma once
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
namespace cfd::workflow {
struct UnitValue { double value{}; std::string unit; };
struct ParameterRule { std::string key;std::string unit;std::optional<double> minimum,maximum;bool required{true}; };
class CaseConfig {
public:
    void set(std::string key,UnitValue value){values_[std::move(key)]=std::move(value);} [[nodiscard]] const UnitValue& at(std::string_view key)const;
    void validate(const std::vector<ParameterRule>& rules)const; [[nodiscard]] std::string to_json()const;
private:std::unordered_map<std::string,UnitValue> values_;
};
struct Material {std::string name;std::unordered_map<std::string,UnitValue> property;};
class MaterialDatabase {
public:void add(Material material);[[nodiscard]] const Material& at(std::string_view name)const;[[nodiscard]] std::string to_json()const;
private:std::unordered_map<std::string,Material> materials_;
};

[[nodiscard]] CaseConfig case_config_from_json(std::string_view json);
[[nodiscard]] CaseConfig case_config_from_yaml(std::string_view yaml);
} // namespace cfd::workflow
