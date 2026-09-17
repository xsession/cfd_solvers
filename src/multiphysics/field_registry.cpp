#include "cfd/multiphysics/field_registry.hpp"

#include <algorithm>
#include <stdexcept>

namespace cfd::multiphysics {

void FieldRegistry::add(FieldMetadata metadata, double initial_value) {
    if (metadata.name.empty() || metadata.components == 0U || metadata.entities == 0U) {
        throw std::invalid_argument("invalid field metadata");
    }
    if (!(metadata.units.scale_to_si > 0.0)) throw std::invalid_argument("field unit scale must be positive");
    const std::string key = metadata.name;
    if (fields_.contains(key)) throw std::invalid_argument("field already exists: " + key);
    Record record;
    record.metadata = std::move(metadata);
    record.data.assign(record.metadata.components * record.metadata.entities, initial_value);
    fields_.emplace(key, std::move(record));
}

bool FieldRegistry::contains(std::string_view name) const noexcept {
    return fields_.find(std::string(name)) != fields_.end();
}

const FieldMetadata& FieldRegistry::metadata(std::string_view name) const {
    const auto it = fields_.find(std::string(name));
    if (it == fields_.end()) throw std::out_of_range("unknown field: " + std::string(name));
    return it->second.metadata;
}

std::span<double> FieldRegistry::data(std::string_view name) {
    auto it = fields_.find(std::string(name));
    if (it == fields_.end()) throw std::out_of_range("unknown field: " + std::string(name));
    return it->second.data;
}

std::span<const double> FieldRegistry::data(std::string_view name) const {
    const auto it = fields_.find(std::string(name));
    if (it == fields_.end()) throw std::out_of_range("unknown field: " + std::string(name));
    return it->second.data;
}

std::vector<std::string> FieldRegistry::names() const {
    std::vector<std::string> out;
    out.reserve(fields_.size());
    for (const auto& item : fields_) out.push_back(item.first);
    std::sort(out.begin(), out.end());
    return out;
}

UnitSignature dimensionless_units() { return {}; }
UnitSignature kelvin_units() { return {{{0,0,0,0,1,0,0}}, 1.0, "K"}; }
UnitSignature metre_units() { return {{{0,1,0,0,0,0,0}}, 1.0, "m"}; }
UnitSignature watt_per_cubic_metre_units() {
    // W/m^3 = kg m^-1 s^-3
    return {{{1,-1,-3,0,0,0,0}}, 1.0, "W/m^3"};
}

} // namespace cfd::multiphysics
