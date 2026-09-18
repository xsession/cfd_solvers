#include "cfd/acoustics/postprocess.hpp"

#include <cmath>
#include <stdexcept>
#include <limits>

namespace cfd::acoustics {

double plane_wave_intensity_w_m2(double p,double rho,double c){
    if(!(p>=0.0)||!(rho>0.0)||!(c>0.0)||!std::isfinite(p)||!std::isfinite(rho)||!std::isfinite(c)) throw std::invalid_argument("invalid acoustic intensity inputs");
    return p*p/(rho*c);
}

double radiation_pressure_pa(double intensity,double c,double reflection){
    if(!(intensity>=0.0)||!(c>0.0)||!(reflection>=0.0&&reflection<=1.0)||!std::isfinite(intensity)||!std::isfinite(c)||!std::isfinite(reflection)) throw std::invalid_argument("invalid radiation pressure inputs");
    return (1.0+reflection)*intensity/c;
}

double acoustic_heating_w_m3(double intensity,double alpha){
    if(!(intensity>=0.0)||!(alpha>=0.0)||!std::isfinite(intensity)||!std::isfinite(alpha)) throw std::invalid_argument("invalid acoustic heating inputs");
    return 2.0*alpha*intensity;
}

std::vector<double> delay_and_sum_beamform(std::span<const BeamformSample> sensors,double fx,double fy,double dt,double c,std::size_t nout){
    if(sensors.empty()||!(dt>0.0)||!(c>0.0)||nout==0U) throw std::invalid_argument("invalid beamforming controls");
    std::vector<double> out(nout,0.0);
    double dmin=std::numeric_limits<double>::infinity();
    std::vector<double> distance(sensors.size());
    for(std::size_t s=0;s<sensors.size();++s){distance[s]=std::hypot(sensors[s].x_m-fx,sensors[s].y_m-fy);dmin=std::min(dmin,distance[s]);}
    for(std::size_t s=0;s<sensors.size();++s){
        if(sensors[s].pressure.empty()) continue;
        const double delay=(distance[s]-dmin)/c/dt;
        for(std::size_t n=0;n<nout;++n){
            const double src=static_cast<double>(n)+delay;
            const auto i0=static_cast<std::size_t>(std::floor(src));
            if(i0+1U>=sensors[s].pressure.size()) continue;
            const double frac=src-static_cast<double>(i0);
            out[n]+=(1.0-frac)*sensors[s].pressure[i0]+frac*sensors[s].pressure[i0+1U];
        }
    }
    const double scale=1.0/static_cast<double>(sensors.size());
    for(double& value:out)value*=scale;
    return out;
}

} // namespace cfd::acoustics
