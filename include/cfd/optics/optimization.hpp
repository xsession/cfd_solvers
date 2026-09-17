#pragma once
#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>
namespace cfd::optics {
using MeritFunction=std::function<double(std::span<const double>)>;
struct OptimizationResult { std::vector<double> variables;double merit{};std::size_t iterations{};bool converged{}; };
[[nodiscard]] OptimizationResult coordinate_descent(MeritFunction merit,std::vector<double> initial,std::vector<double> step,std::size_t max_iterations=200,double tolerance=1e-10);
[[nodiscard]] OptimizationResult differential_evolution(MeritFunction merit,std::span<const double> lower,std::span<const double> upper,std::size_t population=32,std::size_t generations=200,double mutation=0.8,double crossover=0.9,std::uint64_t seed=1);
[[nodiscard]] std::vector<double> finite_difference_sensitivity(const MeritFunction& merit,std::span<const double> point,double relative_step=1e-6);
struct MonteCarloStats { double mean{},standard_deviation{},minimum{},maximum{};std::size_t samples{}; };
[[nodiscard]] MonteCarloStats monte_carlo_tolerance(const MeritFunction& merit,std::span<const double> nominal,std::span<const double> sigma,std::size_t samples,std::uint64_t seed=1);

template<class T> struct Dual {
    T value{};
    T derivative{};
    constexpr Dual()=default;
    constexpr Dual(T v,T d=0):value(v),derivative(d){}
    friend constexpr Dual operator+(Dual a,Dual b){return {a.value+b.value,a.derivative+b.derivative};}
    friend constexpr Dual operator-(Dual a,Dual b){return {a.value-b.value,a.derivative-b.derivative};}
    friend constexpr Dual operator*(Dual a,Dual b){return {a.value*b.value,a.derivative*b.value+a.value*b.derivative};}
    friend constexpr Dual operator/(Dual a,Dual b){return {a.value/b.value,(a.derivative*b.value-a.value*b.derivative)/(b.value*b.value)};}
};
template<class T> constexpr Dual<T> square(Dual<T> x){return x*x;}
template<class Function> [[nodiscard]] double autodiff_derivative(Function&& function,double point){Dual<double>x{point,1.0};return function(x).derivative;}
} // namespace cfd::optics
