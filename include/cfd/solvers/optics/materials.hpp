#pragma once

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::optics {

struct SellmeierMaterial {
    std::string name;
    std::array<double,3> b{};
    std::array<double,3> c_um2{};
    [[nodiscard]] double refractive_index_nm(double wavelength_nm) const;
};

class MaterialCatalog {
public:
    MaterialCatalog();
    void add(SellmeierMaterial material);
    [[nodiscard]] const SellmeierMaterial& at(std::string_view name) const;
    [[nodiscard]] const std::vector<SellmeierMaterial>& materials() const noexcept{return materials_;}
private:
    std::vector<SellmeierMaterial> materials_;
};

struct WavelengthSample {
    double wavelength_nm{550.0};
    double weight{1.0};
};

struct OpticalField {
    double x{};
    double y{};
    bool angular{true};
};

} // namespace cfd::optics
