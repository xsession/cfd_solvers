#pragma once

#include <array>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace cfd::multiphysics {

enum class FieldLocation {
    cell,
    face,
    node,
    edge,
    gridPoint,
    surface,
    global
};

enum class FieldTopology {
    structured,
    unstructured,
    optical,
    meshIndependent
};

// SI base-dimension exponents in the order kg, m, s, A, K, mol, cd.
struct UnitSignature {
    std::array<int, 7> exponent{};
    double scale_to_si{1.0};
    std::string symbol{"1"};

    [[nodiscard]] bool dimensionally_compatible(const UnitSignature& other) const noexcept {
        return exponent == other.exponent;
    }
};

struct FieldMetadata {
    std::string name;
    std::size_t components{1U};
    std::size_t entities{};
    FieldLocation location{FieldLocation::cell};
    FieldTopology topology{FieldTopology::unstructured};
    UnitSignature units{};
    std::string producer;
};

class FieldRegistry {
public:
    void add(FieldMetadata metadata, double initial_value = 0.0);
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] const FieldMetadata& metadata(std::string_view name) const;
    [[nodiscard]] std::span<double> data(std::string_view name);
    [[nodiscard]] std::span<const double> data(std::string_view name) const;
    [[nodiscard]] std::vector<std::string> names() const;

private:
    struct Record {
        FieldMetadata metadata;
        std::vector<double> data;
    };
    std::unordered_map<std::string, Record> fields_;
};

[[nodiscard]] UnitSignature dimensionless_units();
[[nodiscard]] UnitSignature kelvin_units();
[[nodiscard]] UnitSignature metre_units();
[[nodiscard]] UnitSignature watt_per_cubic_metre_units();

} // namespace cfd::multiphysics
