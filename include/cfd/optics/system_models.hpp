#pragma once
#include "cfd/solvers/optics/ray.hpp"
#include <cstddef>
#include <stdexcept>
#include <vector>
namespace cfd::optics {
enum class FieldType { angle, object_height };
struct FieldPoint { double x{}; double y{}; double weight{1.0}; };
class FieldModel {
public:
    explicit FieldModel(FieldType type=FieldType::angle):type_(type){}
    void add(FieldPoint field){if(!(field.weight>=0.0))throw std::invalid_argument("invalid optical field weight");fields_.push_back(field);}
    [[nodiscard]] FieldType type() const noexcept{return type_;} [[nodiscard]] const std::vector<FieldPoint>& fields() const noexcept{return fields_;}
private:FieldType type_;std::vector<FieldPoint> fields_;
};
struct Wavelength { double nm{}; double weight{1.0}; bool primary{}; };
class WavelengthModel {
public:
    void add(Wavelength w){if(!(w.nm>0.0)||!(w.weight>=0.0))throw std::invalid_argument("invalid optical wavelength");if(w.primary)for(auto&x:wavelengths_)x.primary=false;wavelengths_.push_back(w);}
    [[nodiscard]] double primary_nm() const {for(const auto&w:wavelengths_)if(w.primary)return w.nm;if(wavelengths_.empty())throw std::runtime_error("wavelength model is empty");return wavelengths_.front().nm;}
    [[nodiscard]] const std::vector<Wavelength>& wavelengths() const noexcept{return wavelengths_;}
private:std::vector<Wavelength> wavelengths_;
};
enum class ApertureKind { circular, rectangular, elliptical };
struct Aperture { ApertureKind kind{ApertureKind::circular};double x_radius{1.0};double y_radius{1.0};bool stop{};
    [[nodiscard]] bool contains(double x,double y) const {if(!(x_radius>0.0)||!(y_radius>0.0))return false;switch(kind){case ApertureKind::circular:return x*x+y*y<=x_radius*x_radius;case ApertureKind::rectangular:return std::abs(x)<=x_radius&&std::abs(y)<=y_radius;case ApertureKind::elliptical:return x*x/(x_radius*x_radius)+y*y/(y_radius*y_radius)<=1.0;}return false;}
};
} // namespace cfd::optics
