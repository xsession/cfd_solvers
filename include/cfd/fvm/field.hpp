#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
namespace cfd::fvm {
struct Dimensions { std::array<int,7> exponent{}; constexpr bool operator==(const Dimensions&) const noexcept=default; };
inline constexpr Dimensions dimless{};
inline constexpr Dimensions dim_length{{1,0,0,0,0,0,0}};
inline constexpr Dimensions dim_mass{{0,1,0,0,0,0,0}};
inline constexpr Dimensions dim_time{{0,0,1,0,0,0,0}};
inline constexpr Dimensions dim_temperature{{0,0,0,1,0,0,0}};
[[nodiscard]] constexpr Dimensions operator+(Dimensions a,Dimensions b) noexcept {for(std::size_t i=0;i<7;++i)a.exponent[i]+=b.exponent[i];return a;}
[[nodiscard]] constexpr Dimensions operator-(Dimensions a,Dimensions b) noexcept {for(std::size_t i=0;i<7;++i)a.exponent[i]-=b.exponent[i];return a;}
enum class FieldLocation { cell, face, point };
template<class T> class DimensionedField {
public:
    DimensionedField()=default;DimensionedField(std::string name,Dimensions dimensions,FieldLocation location,std::size_t size,T initial={})
      :name_(std::move(name)),dimensions_(dimensions),location_(location),values_(size,initial){}
    [[nodiscard]] const std::string& name()const noexcept{return name_;} [[nodiscard]] Dimensions dimensions()const noexcept{return dimensions_;} [[nodiscard]] FieldLocation location()const noexcept{return location_;}
    [[nodiscard]] std::vector<T>& values()noexcept{return values_;} [[nodiscard]] const std::vector<T>& values()const noexcept{return values_;}
    [[nodiscard]] std::size_t size()const noexcept{return values_.size();}
private:std::string name_;Dimensions dimensions_{};FieldLocation location_{FieldLocation::cell};std::vector<T> values_;
};
template<class T> class FieldHistory {
public:
    explicit FieldHistory(std::size_t levels=2):levels_(levels){if(levels==0)throw std::invalid_argument("field history requires at least one level");}
    void push(std::vector<T> state){if(!states_.empty()&&state.size()!=states_.front().size())throw std::invalid_argument("field history size mismatch");states_.insert(states_.begin(),std::move(state));if(states_.size()>levels_)states_.pop_back();}
    [[nodiscard]] const std::vector<T>& current()const{if(states_.empty())throw std::runtime_error("empty field history");return states_[0];}
    [[nodiscard]] const std::vector<T>& old(std::size_t level=1)const{if(level>=states_.size())throw std::out_of_range("field history level unavailable");return states_[level];}
    [[nodiscard]] std::size_t stored_levels()const noexcept{return states_.size();}
private:std::size_t levels_;std::vector<std::vector<T>> states_;
};
struct Bdf2Coefficients { double current{1.5},old{-2.0},older{0.5}; };
struct CrankNicolsonWeights { double implicit{0.5},explicit_old{0.5}; };
} // namespace cfd::fvm
