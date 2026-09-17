#include "cfd/optics/analysis.hpp"
#include <cmath>
#include <stdexcept>
namespace cfd::optics {
SpotDiagram spot_diagram_at_plane(const std::vector<RayTraceResult>& rays,double z){SpotDiagram d;for(const auto&r:rays){if(!r.valid||std::abs(r.ray.direction.z)<1e-14)continue;const double t=(z-r.ray.origin.z)/r.ray.direction.z;if(t<0)continue;d.points.push_back({r.ray.origin.x+t*r.ray.direction.x,r.ray.origin.y+t*r.ray.direction.y});}if(d.points.empty())throw std::runtime_error("no rays reach spot plane");for(auto p:d.points){d.centroid_x+=p.x;d.centroid_y+=p.y;}d.centroid_x/=d.points.size();d.centroid_y/=d.points.size();for(auto p:d.points){const double x=p.x-d.centroid_x,y=p.y-d.centroid_y;d.rms_radius+=x*x+y*y;}d.rms_radius=std::sqrt(d.rms_radius/d.points.size());return d;}
std::vector<RayFanPoint> tangential_ray_fan(const std::vector<RayTraceResult>& rays,double z,double ref){std::vector<RayFanPoint> out;out.reserve(rays.size());double maxp=0;for(const auto&r:rays)if(r.valid)maxp=std::max(maxp,std::abs(r.ray.origin.y));if(maxp==0)maxp=1;for(const auto&r:rays){if(!r.valid||std::abs(r.ray.direction.z)<1e-14)continue;const double t=(z-r.ray.origin.z)/r.ray.direction.z;if(t<0)continue;out.push_back({r.ray.origin.y/maxp,r.ray.origin.y+t*r.ray.direction.y-ref});}return out;}
double distortion_percent(double actual,double ideal){if(ideal==0.0)throw std::invalid_argument("ideal image height must be nonzero");return 100.0*(actual-ideal)/ideal;}
} // namespace cfd::optics
