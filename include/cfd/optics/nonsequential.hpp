#pragma once
#include "cfd/solvers/optics/ray.hpp"
#include <cstddef>
#include <vector>
#include <span>
namespace cfd::optics {
struct NsPlane {Vec3 point{};Vec3 normal{0,0,1};double refractive_index_before{1.0},refractive_index_after{1.5};double reflectance{};double transmittance{1.0};double absorptance{};double scatter_fraction{};};
struct NsRay {Ray ray{};double power{1.0};std::size_t bounces{};};
class NonSequentialScene {
public:void add_plane(NsPlane p);[[nodiscard]] std::vector<NsRay> trace(NsRay ray,std::size_t max_bounces=4,double min_power=1e-6)const;
private:std::vector<NsPlane>planes_;
};

struct StrayLightSummary {double total_power{},reflected_power{},transmitted_power{},scattered_power{},absorbed_power{};};
[[nodiscard]] StrayLightSummary summarize_stray_light(std::span<const NsRay> rays,double launched_power=1.0);
struct GhostPath {std::size_t bounces{};double power{};Ray ray{};};
[[nodiscard]] std::vector<GhostPath> rank_ghost_paths(std::span<const NsRay> rays,std::size_t minimum_bounces=2);

} // namespace cfd::optics
