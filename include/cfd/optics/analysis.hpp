#pragma once
#include "cfd/solvers/optics/sequential.hpp"
#include <vector>
namespace cfd::optics {
struct SpotPoint { double x{},y{}; };
struct SpotDiagram { std::vector<SpotPoint> points;double centroid_x{},centroid_y{},rms_radius{}; };
[[nodiscard]] SpotDiagram spot_diagram_at_plane(const std::vector<RayTraceResult>& rays,double image_z);
struct RayFanPoint { double pupil_coordinate{};double transverse_aberration{}; };
[[nodiscard]] std::vector<RayFanPoint> tangential_ray_fan(const std::vector<RayTraceResult>& rays,double image_z,double reference_x=0.0);
[[nodiscard]] double distortion_percent(double actual_image_height,double ideal_image_height);
} // namespace cfd::optics
