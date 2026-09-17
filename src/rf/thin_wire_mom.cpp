#include "cfd/rf/thin_wire_mom.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>
#include <stdexcept>

namespace cfd::rf {
namespace {
using Complex = std::complex<double>;
constexpr double eps0 = 8.8541878128e-12;
constexpr double mu0 = 1.25663706212e-6;

struct Quadrature {
    std::vector<double> x;
    std::vector<double> w;
};

Quadrature gauss_legendre(std::size_t n) {
    if (n < 2U || n > 32U) throw std::invalid_argument("thin-wire MoM quadrature order must be 2..32");
    Quadrature q; q.x.resize(n); q.w.resize(n);
    const std::size_t half=(n+1U)/2U;
    for(std::size_t i=0;i<half;++i){
        double z=std::cos(std::numbers::pi*(static_cast<double>(i)+0.75)/(static_cast<double>(n)+0.5));
        double previous=0.0,derivative=0.0;
        for(std::size_t iteration=0;iteration<50U;++iteration){
            double p0=1.0,p1=z;
            if(n==1U){derivative=1.0;}
            else{
                for(std::size_t k=2;k<=n;++k){
                    const double kd=static_cast<double>(k);
                    const double p=((2.0*kd-1.0)*z*p1-(kd-1.0)*p0)/kd;
                    p0=p1;p1=p;
                }
                derivative=static_cast<double>(n)*(z*p1-p0)/(z*z-1.0);
                previous=z; z-=p1/derivative;
                if(std::abs(z-previous)<1e-15)break;
                continue;
            }
        }
        double p0=1.0,p1=z;
        for(std::size_t k=2;k<=n;++k){
            const double kd=static_cast<double>(k);
            const double p=((2.0*kd-1.0)*z*p1-(kd-1.0)*p0)/kd;
            p0=p1;p1=p;
        }
        derivative=static_cast<double>(n)*(z*p1-p0)/(z*z-1.0);
        const double weight=2.0/((1.0-z*z)*derivative*derivative);
        q.x[i]=-z;q.x[n-1U-i]=z;q.w[i]=weight;q.w[n-1U-i]=weight;
    }
    return q;
}

template<class T>
std::vector<T> solve_dense(std::vector<T> a,std::vector<T> b,std::size_t n){
    for(std::size_t col=0;col<n;++col){
        std::size_t pivot=col;double best=std::abs(a[col*n+col]);
        for(std::size_t row=col+1U;row<n;++row){const double v=std::abs(a[row*n+col]);if(v>best){best=v;pivot=row;}}
        if(best<1e-24)throw std::runtime_error("singular thin-wire MoM matrix");
        if(pivot!=col){for(std::size_t j=col;j<n;++j)std::swap(a[pivot*n+j],a[col*n+j]);std::swap(b[pivot],b[col]);}
        const T d=a[col*n+col];
        for(std::size_t j=col;j<n;++j)a[col*n+j]/=d;
        b[col]/=d;
        for(std::size_t row=0;row<n;++row){
            if(row==col)continue;
            const T f=a[row*n+col];
            if(std::abs(f)==0.0)continue;
            for(std::size_t j=col;j<n;++j)a[row*n+j]-=f*a[col*n+j];
            b[row]-=f*b[col];
        }
    }
    return b;
}

WirePoint3 add(WirePoint3 a,WirePoint3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
WirePoint3 sub(WirePoint3 a,WirePoint3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
WirePoint3 mul(WirePoint3 a,double s){return {a.x*s,a.y*s,a.z*s};}
double dot(WirePoint3 a,WirePoint3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double norm(WirePoint3 a){return std::sqrt(dot(a,a));}
bool finite(WirePoint3 a){return std::isfinite(a.x)&&std::isfinite(a.y)&&std::isfinite(a.z);}

double segment_distance(WirePoint3 p1,WirePoint3 q1,WirePoint3 p2,WirePoint3 q2){
    const WirePoint3 d1=sub(q1,p1),d2=sub(q2,p2),r=sub(p1,p2);
    const double a=dot(d1,d1),e=dot(d2,d2),f=dot(d2,r);
    double s=0.0,t=0.0;
    constexpr double tiny=1.0e-30;
    if(a<=tiny&&e<=tiny)return norm(r);
    if(a<=tiny){t=std::clamp(f/e,0.0,1.0);}
    else{
        const double c=dot(d1,r);
        if(e<=tiny)s=std::clamp(-c/a,0.0,1.0);
        else{
            const double b=dot(d1,d2),denom=a*e-b*b;
            if(std::abs(denom)>tiny)s=std::clamp((b*f-c*e)/denom,0.0,1.0);
            t=(b*s+f)/e;
            if(t<0.0){t=0.0;s=std::clamp(-c/a,0.0,1.0);}
            else if(t>1.0){t=1.0;s=std::clamp((b-c)/a,0.0,1.0);}
        }
    }
    return norm(sub(add(p1,mul(d1,s)),add(p2,mul(d2,t))));
}

struct JunctionEndpoint {
    std::size_t wire{};
    bool start{};
};

std::vector<std::vector<JunctionEndpoint>> endpoint_junctions(const OrientedWireMomConfig& config) {
    if(config.junction_tolerance_m==0.0)return {};
    if(!(config.junction_tolerance_m>0.0)||!std::isfinite(config.junction_tolerance_m))
        throw std::invalid_argument("junction tolerance must be finite and nonnegative");
    struct Endpoint { WirePoint3 point; std::size_t wire; bool start; };
    std::vector<Endpoint> endpoints;endpoints.reserve(2U*config.wires.size());
    double minimum_length=std::numeric_limits<double>::infinity();
    for(std::size_t wire=0;wire<config.wires.size();++wire){
        endpoints.push_back({config.wires[wire].start_m,wire,true});
        endpoints.push_back({config.wires[wire].end_m,wire,false});
        minimum_length=std::min(minimum_length,norm(sub(config.wires[wire].end_m,config.wires[wire].start_m)));
    }
    if(config.junction_tolerance_m>=0.25*minimum_length)
        throw std::invalid_argument("junction tolerance is too large relative to wire length");
    std::vector<std::size_t> parent(endpoints.size());std::iota(parent.begin(),parent.end(),0U);
    const auto find_root=[&](std::size_t x){while(parent[x]!=x){parent[x]=parent[parent[x]];x=parent[x];}return x;};
    const auto unite=[&](std::size_t a,std::size_t b){
        a=find_root(a);b=find_root(b);if(a!=b)parent[b]=a;
    };
    for(std::size_t a=0;a<endpoints.size();++a)for(std::size_t b=a+1U;b<endpoints.size();++b){
        if(endpoints[a].wire==endpoints[b].wire)continue;
        if(norm(sub(endpoints[a].point,endpoints[b].point))<=config.junction_tolerance_m)unite(a,b);
    }
    std::vector<std::vector<JunctionEndpoint>> groups;
    std::vector<std::size_t> roots;
    for(std::size_t i=0;i<endpoints.size();++i){
        const auto root=find_root(i);auto found=std::find(roots.begin(),roots.end(),root);
        std::size_t group{};
        if(found==roots.end()){roots.push_back(root);groups.emplace_back();group=groups.size()-1U;}
        else group=static_cast<std::size_t>(found-roots.begin());
        groups[group].push_back({endpoints[i].wire,endpoints[i].start});
    }
    groups.erase(std::remove_if(groups.begin(),groups.end(),[](const auto& group){return group.size()<2U;}),groups.end());
    return groups;
}

bool endpoint_pair_is_junction(const OrientedThinWire& a,const OrientedThinWire& b,double tolerance) {
    if(!(tolerance>0.0))return false;
    return norm(sub(a.start_m,b.start_m))<=tolerance||norm(sub(a.start_m,b.end_m))<=tolerance
        ||norm(sub(a.end_m,b.start_m))<=tolerance||norm(sub(a.end_m,b.end_m))<=tolerance;
}

bool endpoint_junction_overlaps(const OrientedThinWire& a,const OrientedThinWire& b,double tolerance) {
    struct Candidate { WirePoint3 ap,bp; bool as,bs; };
    const Candidate candidates[]={{a.start_m,b.start_m,true,true},{a.start_m,b.end_m,true,false},
                                  {a.end_m,b.start_m,false,true},{a.end_m,b.end_m,false,false}};
    for(const auto& candidate:candidates){
        if(norm(sub(candidate.ap,candidate.bp))>tolerance)continue;
        const WirePoint3 ad=candidate.as?sub(a.end_m,a.start_m):sub(a.start_m,a.end_m);
        const WirePoint3 bd=candidate.bs?sub(b.end_m,b.start_m):sub(b.start_m,b.end_m);
        const double cosine=dot(ad,bd)/(norm(ad)*norm(bd));
        if(cosine>1.0-1.0e-10)return true;
    }
    return false;
}

Complex kernel_integral_oriented(const OrientedWireSegment& observation,
                                 const OrientedWireSegment& source,
                                 double source_radius,bool same_wire,
                                 double wave_number,double omega,const Quadrature& q){
    const Complex j{0.0,1.0};
    const double tangent_dot=dot(observation.tangent,source.tangent);
    Complex sum{};
    for(std::size_t k=0;k<q.x.size();++k){
        const WirePoint3 source_point=add(source.center_m,mul(source.tangent,0.5*source.length_m*q.x[k]));
        const WirePoint3 displacement=sub(observation.center_m,source_point);
        const double smooth2=same_wire?source_radius*source_radius:0.0;
        const double r=std::sqrt(dot(displacement,displacement)+smooth2);
        if(!(r>0.0)||!std::isfinite(r))throw std::runtime_error("singular oriented thin-wire interaction");
        const double inv_r=1.0/r,inv_r2=inv_r*inv_r,inv_r3=inv_r2*inv_r;
        const double obs_radial=dot(observation.tangent,displacement)*inv_r;
        const double src_radial=dot(source.tangent,displacement)*inv_r;
        const Complex phase=std::exp(-j*wave_number*r)/(4.0*std::numbers::pi);
        const Complex radial=(-wave_number*wave_number*inv_r+3.0*j*wave_number*inv_r2+3.0*inv_r3)
            *obs_radial*src_radial;
        const Complex isotropic=(-j*wave_number*inv_r2-inv_r3)*tangent_dot;
        const Complex vector_term=wave_number*wave_number*tangent_dot*inv_r;
        sum+=q.w[k]*phase*(vector_term+radial+isotropic);
    }
    return j/(omega*eps0)*(0.5*source.length_m)*sum;
}

OrientedWireSegment pec_image_segment(const OrientedWireSegment& source,double plane_z){
    OrientedWireSegment image=source;
    image.center_m.z=2.0*plane_z-source.center_m.z;
    // Electric-current image over a PEC plane: tangential current reverses,
    // normal current keeps its sign. Reversing the mirrored geometric tangent
    // spans the same source segment while providing exactly that image current.
    image.tangent={-source.tangent.x,-source.tangent.y,source.tangent.z};
    return image;
}

struct LineZParameters { Complex z11{},z12{},z21{},z22{}; };
LineZParameters line_z_parameters(const ThinWireTransmissionLineLoad& load){
    if(!(load.length_m>0.0)||!std::isfinite(load.length_m)
       ||!std::isfinite(load.characteristic_impedance_ohm.real())||!std::isfinite(load.characteristic_impedance_ohm.imag())
       ||!(load.characteristic_impedance_ohm.real()>0.0)
       ||!std::isfinite(load.propagation_constant_per_m.real())||!std::isfinite(load.propagation_constant_per_m.imag())
       ||(load.first_current_sign!=1&&load.first_current_sign!=-1)
       ||(load.second_current_sign!=1&&load.second_current_sign!=-1))
        throw std::invalid_argument("invalid thin-wire transmission-line load");
    const Complex gl=load.propagation_constant_per_m*load.length_m;
    const Complex sh=std::sinh(gl),ch=std::cosh(gl);
    if(std::abs(sh)<1.0e-14)throw std::invalid_argument("thin-wire transmission-line section is singular at this electrical length");
    const Complex diagonal=load.characteristic_impedance_ohm*ch/sh;
    const Complex transfer=load.characteristic_impedance_ohm/sh;
    return {diagonal,transfer,transfer,diagonal};
}

std::size_t feed_segment(std::size_t segments,const ThinWireDeltaGapFeed& feed){
    if(feed.segment==std::numeric_limits<std::size_t>::max()){
        if((segments%2U)==0U)throw std::invalid_argument("center delta-gap feed requires odd wire segmentation");
        return segments/2U;
    }
    if(feed.segment>=segments)throw std::out_of_range("thin-wire feed segment out of range");
    return feed.segment;
}
}

OrientedWireMomResult solve_oriented_thin_wires(const OrientedWireMomConfig& config){
    if(!(config.frequency_hz>0.0)||!std::isfinite(config.frequency_hz)||config.wires.empty()||config.feeds.empty())
        throw std::invalid_argument("invalid oriented thin-wire MoM configuration");
    std::vector<std::size_t> offsets(config.wires.size()+1U,0U);
    std::vector<double> lengths(config.wires.size());
    for(std::size_t w=0;w<config.wires.size();++w){
        const auto& wire=config.wires[w];const double length=norm(sub(wire.end_m,wire.start_m));lengths[w]=length;
        if(!finite(wire.start_m)||!finite(wire.end_m)||!(length>0.0)||!(wire.radius_m>0.0)
           ||!std::isfinite(wire.radius_m)||wire.radius_m>=0.1*length||wire.segments<5U
           ||wire.conductivity_s_per_m<0.0||!std::isfinite(wire.conductivity_s_per_m)
           ||!(wire.relative_permeability>0.0)||!std::isfinite(wire.relative_permeability)
           ||wire.dielectric_loss_ohm_per_m<0.0||!std::isfinite(wire.dielectric_loss_ohm_per_m))
            throw std::invalid_argument("invalid oriented thin wire");
        if(config.ground_model==ThinWireGroundModel::pec_plane){
            const double minimum_z=std::min(wire.start_m.z,wire.end_m.z);
            if(!(minimum_z-config.ground_plane_z_m>wire.radius_m))
                throw std::invalid_argument("PEC-ground thin wires must remain above the ground plane by more than their radius");
        }
        const double ds=length/static_cast<double>(wire.segments);
        if(wire.radius_m>=0.5*ds)throw std::invalid_argument("oriented thin-wire MoM requires radius < half segment length");
        if(offsets[w]>std::numeric_limits<std::size_t>::max()-wire.segments)throw std::overflow_error("thin-wire segment count overflow");
        offsets[w+1U]=offsets[w]+wire.segments;
    }
    for(std::size_t a=0;a<config.wires.size();++a)for(std::size_t b=a+1U;b<config.wires.size();++b){
        const double separation=segment_distance(config.wires[a].start_m,config.wires[a].end_m,
                                                  config.wires[b].start_m,config.wires[b].end_m);
        if(!(separation>config.wires[a].radius_m+config.wires[b].radius_m)){
            const bool endpoint_junction=endpoint_pair_is_junction(config.wires[a],config.wires[b],config.junction_tolerance_m);
            if(!endpoint_junction||endpoint_junction_overlaps(config.wires[a],config.wires[b],config.junction_tolerance_m))
                throw std::invalid_argument("oriented thin-wire axes overlap/cross outside an enabled endpoint junction");
        }
    }
    const auto junctions=endpoint_junctions(config);
    const std::size_t n=offsets.back();
    OrientedWireMomResult result;result.segments.reserve(n);
    for(std::size_t w=0;w<config.wires.size();++w){
        const auto& wire=config.wires[w];const WirePoint3 delta=sub(wire.end_m,wire.start_m);
        const double length=lengths[w],ds=length/static_cast<double>(wire.segments);const WirePoint3 tangent=mul(delta,1.0/length);
        for(std::size_t s=0;s<wire.segments;++s){
            const double fraction=(static_cast<double>(s)+0.5)/static_cast<double>(wire.segments);
            result.segments.push_back({w,s,add(wire.start_m,mul(delta,fraction)),tangent,ds});
        }
    }
    const double c=1.0/std::sqrt(mu0*eps0),omega=2.0*std::numbers::pi*config.frequency_hz,wave_number=omega/c;
    const auto quadrature=gauss_legendre(config.quadrature_order);
    std::vector<Complex> matrix(n*n),rhs(n,Complex{});
    std::vector<double> conductor_resistance_per_m(config.wires.size(),0.0);
    for(std::size_t wire=0;wire<config.wires.size();++wire){
        const auto& definition=config.wires[wire];
        if(definition.conductivity_s_per_m>0.0){
            const double surface_resistance=std::sqrt(std::numbers::pi*config.frequency_hz*mu0*definition.relative_permeability
                                                      /definition.conductivity_s_per_m);
            conductor_resistance_per_m[wire]=surface_resistance/(2.0*std::numbers::pi*definition.radius_m);
        }
    }
    for(std::size_t row=0;row<n;++row){
        const auto& observation=result.segments[row];
        for(std::size_t col=0;col<n;++col){
            const auto& source=result.segments[col];
            matrix[row*n+col]=kernel_integral_oriented(observation,source,config.wires[source.wire].radius_m,
                                                       observation.wire==source.wire,wave_number,omega,quadrature);
            if(config.ground_model==ThinWireGroundModel::pec_plane){
                const auto image=pec_image_segment(source,config.ground_plane_z_m);
                matrix[row*n+col]+=kernel_integral_oriented(observation,image,config.wires[source.wire].radius_m,
                                                            false,wave_number,omega,quadrature);
            }
        }
        const auto& wire=config.wires[observation.wire];
        matrix[row*n+row]+=Complex{conductor_resistance_per_m[observation.wire]+wire.dielectric_loss_ohm_per_m,0.0};
    }
    for(const auto& load:config.loads){
        if(load.wire>=config.wires.size()||load.segment>=config.wires[load.wire].segments
           ||!std::isfinite(load.impedance_ohm.real())||!std::isfinite(load.impedance_ohm.imag()))
            throw std::invalid_argument("invalid thin-wire series load");
        const std::size_t index=offsets[load.wire]+load.segment;
        matrix[index*n+index]+=load.impedance_ohm/result.segments[index].length_m;
    }
    for(const auto& load:config.transmission_line_loads){
        if(load.first_wire>=config.wires.size()||load.second_wire>=config.wires.size()
           ||load.first_segment>=config.wires[load.first_wire].segments
           ||load.second_segment>=config.wires[load.second_wire].segments)
            throw std::invalid_argument("thin-wire transmission-line load segment out of range");
        const auto z=line_z_parameters(load);
        const std::size_t first=offsets[load.first_wire]+load.first_segment;
        const std::size_t second=offsets[load.second_wire]+load.second_segment;
        if(first==second)throw std::invalid_argument("thin-wire transmission-line ports must use distinct segments");
        const double s1=static_cast<double>(load.first_current_sign),s2=static_cast<double>(load.second_current_sign);
        matrix[first*n+first]+=s1*s1*z.z11/result.segments[first].length_m;
        matrix[first*n+second]+=s1*s2*z.z12/result.segments[first].length_m;
        matrix[second*n+first]+=s2*s1*z.z21/result.segments[second].length_m;
        matrix[second*n+second]+=s2*s2*z.z22/result.segments[second].length_m;
    }
    std::vector<std::size_t> feed_indices;feed_indices.reserve(config.feeds.size());
    for(const auto& feed:config.feeds){
        if(feed.wire>=config.wires.size()||std::abs(feed.voltage_v)==0.0
           ||!std::isfinite(feed.voltage_v.real())||!std::isfinite(feed.voltage_v.imag()))
            throw std::invalid_argument("invalid thin-wire delta-gap feed");
        const std::size_t local=feed_segment(config.wires[feed.wire].segments,feed);
        const std::size_t index=offsets[feed.wire]+local;feed_indices.push_back(index);
        rhs[index]-=feed.voltage_v/result.segments[index].length_m;
    }
    if(junctions.empty()) result.current_a=solve_dense(std::move(matrix),std::move(rhs),n);
    else {
        const std::size_t dimension=n+junctions.size();
        std::vector<Complex> constrained(dimension*dimension,Complex{});
        std::vector<Complex> constrained_rhs(dimension,Complex{});
        double constraint_scale=0.0;
        for(std::size_t row=0;row<n;++row){
            constrained_rhs[row]=rhs[row];
            constraint_scale=std::max(constraint_scale,std::abs(matrix[row*n+row]));
            for(std::size_t column=0;column<n;++column)constrained[row*dimension+column]=matrix[row*n+column];
        }
        constraint_scale=std::max(constraint_scale,1.0);
        for(std::size_t junc=0;junc<junctions.size();++junc){
            const std::size_t row=n+junc;
            for(const auto& endpoint:junctions[junc]){
                const std::size_t segment=endpoint.start?offsets[endpoint.wire]:offsets[endpoint.wire+1U]-1U;
                const double sign=endpoint.start?1.0:-1.0;
                const Complex coefficient{constraint_scale*sign,0.0};
                constrained[row*dimension+segment]=coefficient;
                constrained[segment*dimension+row]=coefficient;
            }
        }
        auto solution=solve_dense(std::move(constrained),std::move(constrained_rhs),dimension);
        solution.resize(n);result.current_a=std::move(solution);
    }
    const Complex first_impedance=config.feeds.front().voltage_v/result.current_a[feed_indices.front()];
    if(first_impedance.real()<0.0)for(auto& current:result.current_a)current=-current;
    result.feed_impedance_ohm.reserve(config.feeds.size());result.feed_current_a.reserve(config.feeds.size());
    for(std::size_t i=0;i<config.feeds.size();++i){
        const Complex current=result.current_a[feed_indices[i]];
        if(std::abs(current)<1e-30)throw std::runtime_error("thin-wire MoM zero feed current");
        result.feed_impedance_ohm.push_back(config.feeds[i].voltage_v/current);
        result.feed_current_a.push_back(std::abs(current));
        result.accepted_power_w+=0.5*std::real(config.feeds[i].voltage_v*std::conj(current));
    }
    for(std::size_t segment=0;segment<result.segments.size();++segment){
        const auto& geometry=result.segments[segment];const auto& wire=config.wires[geometry.wire];
        const double current2=std::norm(result.current_a[segment]);
        result.conductor_loss_w+=0.5*conductor_resistance_per_m[geometry.wire]*geometry.length_m*current2;
        result.dielectric_loss_w+=0.5*wire.dielectric_loss_ohm_per_m*geometry.length_m*current2;
    }
    for(const auto& load:config.loads){
        if(load.impedance_ohm.real()>0.0){
            const std::size_t segment=offsets[load.wire]+load.segment;
            result.lumped_load_loss_w+=0.5*load.impedance_ohm.real()*std::norm(result.current_a[segment]);
        }
    }
    for(const auto& load:config.transmission_line_loads){
        const auto z=line_z_parameters(load);
        const std::size_t first=offsets[load.first_wire]+load.first_segment;
        const std::size_t second=offsets[load.second_wire]+load.second_segment;
        const Complex i1=static_cast<double>(load.first_current_sign)*result.current_a[first];
        const Complex i2=static_cast<double>(load.second_current_sign)*result.current_a[second];
        const Complex v1=z.z11*i1+z.z12*i2,v2=z.z21*i1+z.z22*i2;
        const double power=0.5*std::real(v1*std::conj(i1)+v2*std::conj(i2));
        result.transmission_line_loss_w+=std::max(0.0,power);
    }
    const double dissipated=result.conductor_loss_w+result.dielectric_loss_w+result.lumped_load_loss_w+result.transmission_line_loss_w;
    result.radiated_power_w=std::max(0.0,result.accepted_power_w-dissipated);
    result.radiation_efficiency=result.accepted_power_w>0.0
        ?std::clamp(result.radiated_power_w/result.accepted_power_w,0.0,1.0):0.0;
    return result;
}


ComplexMatrix oriented_wire_port_impedance_matrix(const OrientedWireMomConfig& config){
    if(config.feeds.empty())throw std::invalid_argument("wire multiport extraction requires at least one feed");
    const std::size_t ports=config.feeds.size();
    std::vector<std::size_t> offsets(config.wires.size()+1U,0U);
    for(std::size_t w=0;w<config.wires.size();++w)offsets[w+1U]=offsets[w]+config.wires[w].segments;
    std::vector<std::size_t> indices(ports);
    for(std::size_t p=0;p<ports;++p){
        const auto& feed=config.feeds[p];
        if(feed.wire>=config.wires.size())throw std::invalid_argument("wire multiport feed wire out of range");
        indices[p]=offsets[feed.wire]+feed_segment(config.wires[feed.wire].segments,feed);
    }
    ComplexMatrix admittance(ports);
    for(std::size_t excited=0;excited<ports;++excited){
        OrientedWireMomConfig single=config;
        single.feeds.clear();
        auto feed=config.feeds[excited];
        feed.voltage_v={1.0,0.0};
        single.feeds.push_back(feed);
        const auto solved=solve_oriented_thin_wires(single);
        for(std::size_t observed=0;observed<ports;++observed)
            admittance(observed,excited)=solved.current_a[indices[observed]];
    }
    return inverse(admittance);
}

ComplexMatrix oriented_wire_port_s_parameters(const OrientedWireMomConfig& config,
                                               std::span<const double> reference_impedance){
    return z_to_s(oriented_wire_port_impedance_matrix(config),reference_impedance);
}

ParallelWireMomResult solve_parallel_thin_wires(const ParallelWireMomConfig& config){
    OrientedWireMomConfig oriented;oriented.frequency_hz=config.frequency_hz;oriented.feeds=config.feeds;
    oriented.loads=config.loads;oriented.quadrature_order=config.quadrature_order;oriented.wires.reserve(config.wires.size());
    for(const auto& wire:config.wires){
        oriented.wires.push_back({{wire.x_m,wire.y_m,wire.center_z_m-0.5*wire.length_m},
                                  {wire.x_m,wire.y_m,wire.center_z_m+0.5*wire.length_m},
                                  wire.radius_m,wire.segments});
    }
    const auto solved=solve_oriented_thin_wires(oriented);
    ParallelWireMomResult result;result.current_a=solved.current_a;result.feed_impedance_ohm=solved.feed_impedance_ohm;
    result.feed_current_a=solved.feed_current_a;result.segments.reserve(solved.segments.size());
    for(const auto& segment:solved.segments)
        result.segments.push_back({segment.wire,segment.local_segment,segment.center_m.x,segment.center_m.y,
                                   segment.center_m.z,segment.length_m});
    return result;
}

ThinWireMomResult solve_center_fed_thin_wire(const ThinWireMomConfig& config){
    if(!(config.length_m>0.0)||!(config.radius_m>0.0)||!(config.frequency_hz>0.0)
       ||!std::isfinite(config.length_m)||!std::isfinite(config.radius_m)||!std::isfinite(config.frequency_hz)
       ||config.radius_m>=0.1*config.length_m||config.segments<5U||(config.segments%2U)==0U
       ||!(config.feed_voltage_v!=0.0)||!std::isfinite(config.feed_voltage_v))
        throw std::invalid_argument("invalid thin-wire MoM configuration");
    ParallelWireMomConfig multi;
    multi.frequency_hz=config.frequency_hz;multi.quadrature_order=config.quadrature_order;
    multi.wires.push_back({0.0,0.0,0.0,config.length_m,config.radius_m,config.segments});
    multi.feeds.push_back({0U,std::numeric_limits<std::size_t>::max(),Complex{config.feed_voltage_v,0.0}});
    const auto solved=solve_parallel_thin_wires(multi);
    ThinWireMomResult result;result.current_a=solved.current_a;result.segment_center_z_m.reserve(solved.segments.size());
    for(const auto& segment:solved.segments)result.segment_center_z_m.push_back(segment.center_z_m);
    result.feed_impedance_ohm=solved.feed_impedance_ohm.front();
    result.feed_current_a=solved.feed_current_a.front();
    return result;
}

} // namespace cfd::rf
