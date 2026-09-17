#include "cfd/solvers/fem/magnetostatics2d.hpp"

#include "cfd/fem/reference_element.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfd::fem {

Magnetostatics2D::Magnetostatics2D(Mesh2D mesh,ScalarDiffusion2DConfig config):solver_(std::move(mesh),config){}
void Magnetostatics2D::set_boundary(int patch,ScalarBoundaryCondition2D condition){solver_.set_boundary(patch,std::move(condition));}
void Magnetostatics2D::solve(const std::function<double(Node2)>& reluctivity,const std::function<double(Node2)>& current_density_z){
    if(!reluctivity||!current_density_z)throw std::invalid_argument("magnetostatic callback missing");
    solver_.solve(reluctivity,[](Node2){return 0.0;},current_density_z);
}
std::vector<MagneticFluxDensity2> Magnetostatics2D::element_flux_density() const{
    std::vector<MagneticFluxDensity2> out;out.reserve(mesh().triangles.size());
    for(const auto&tri:mesh().triangles){
        std::array<Point3,3> node{{{mesh().nodes[tri.node[0]].x,mesh().nodes[tri.node[0]].y,0.0},
                                  {mesh().nodes[tri.node[1]].x,mesh().nodes[tri.node[1]].y,0.0},
                                  {mesh().nodes[tri.node[2]].x,mesh().nodes[tri.node[2]].y,0.0}}};
        const auto iso=evaluate_isoparametric(ElementType::tri3,node,{1.0/3.0,1.0/3.0,0.0});
        double gx=0.0,gy=0.0;
        for(std::size_t i=0;i<3U;++i){gx+=vector_potential()[tri.node[i]]*iso.gradient_physical[i].x;gy+=vector_potential()[tri.node[i]]*iso.gradient_physical[i].y;}
        out.push_back({gy,-gx});
    }
    return out;
}


Magnetostatics2D::NonlinearResult Magnetostatics2D::solve_nonlinear(
    const std::function<double(double)>& law,const std::function<double(Node2)>& current,
    std::size_t max_iterations,double tolerance,double relaxation){
    if(!law||!current||max_iterations==0U||!(tolerance>0.0)||!(relaxation>0.0&&relaxation<=1.0))throw std::invalid_argument("invalid nonlinear magnetostatic controls");
    std::vector<double> nu(mesh().triangles.size(),law(0.0));for(double value:nu)if(!(value>0.0)||!std::isfinite(value))throw std::invalid_argument("invalid B-H reluctivity law");
    auto element_for_point=[&](Node2 p){
        std::size_t nearest=0U;double best=std::numeric_limits<double>::infinity();
        for(std::size_t e=0;e<mesh().triangles.size();++e){const auto& t=mesh().triangles[e];const auto&a=mesh().nodes[t.node[0]],&b=mesh().nodes[t.node[1]],&c=mesh().nodes[t.node[2]];const double det=(b.y-c.y)*(a.x-c.x)+(c.x-b.x)*(a.y-c.y);const double l0=((b.y-c.y)*(p.x-c.x)+(c.x-b.x)*(p.y-c.y))/det;const double l1=((c.y-a.y)*(p.x-c.x)+(a.x-c.x)*(p.y-c.y))/det;const double l2=1.0-l0-l1;if(l0>=-1.0e-10&&l1>=-1.0e-10&&l2>=-1.0e-10)return e;const double cx=(a.x+b.x+c.x)/3.0,cy=(a.y+b.y+c.y)/3.0,d2=(p.x-cx)*(p.x-cx)+(p.y-cy)*(p.y-cy);if(d2<best){best=d2;nearest=e;}}
        return nearest;
    };
    NonlinearResult result;
    for(std::size_t iteration=0;iteration<max_iterations;++iteration){solver_.solve([&](Node2 p){return nu[element_for_point(p)];},[](Node2){return 0.0;},current);const auto b=element_flux_density();double max_change=0.0;std::vector<double> updated(nu.size());for(std::size_t e=0;e<nu.size();++e){const double magnitude=std::hypot(b[e].x,b[e].y),target=law(magnitude);if(!(target>0.0)||!std::isfinite(target))throw std::runtime_error("nonlinear B-H law returned invalid reluctivity");updated[e]=(1.0-relaxation)*nu[e]+relaxation*target;max_change=std::max(max_change,std::abs(updated[e]-nu[e])/std::max(std::abs(updated[e]),1.0e-30));}nu=std::move(updated);result={iteration+1U,max_change,max_change<=tolerance};if(result.converged)return result;}
    return result;
}

PiecewiseLinearBHCurve::PiecewiseLinearBHCurve(std::vector<double> b,std::vector<double> h):b_(std::move(b)),h_(std::move(h)){
    if(b_.size()<2U||b_.size()!=h_.size()||b_.front()!=0.0||h_.front()!=0.0)throw std::invalid_argument("B-H curve must start at origin and contain at least two samples");for(std::size_t i=1;i<b_.size();++i)if(!(b_[i]>b_[i-1])||!(h_[i]>h_[i-1]))throw std::invalid_argument("B-H samples must be strictly increasing");
}
double PiecewiseLinearBHCurve::field_strength(double b) const {if(!(b>=0.0)||!std::isfinite(b))throw std::invalid_argument("invalid B-H query");if(b>=b_.back()){const std::size_t n=b_.size();const double slope=(h_[n-1]-h_[n-2])/(b_[n-1]-b_[n-2]);return h_.back()+slope*(b-b_.back());}const auto it=std::upper_bound(b_.begin(),b_.end(),b);const std::size_t hi=static_cast<std::size_t>(it-b_.begin()),lo=hi-1U;const double t=(b-b_[lo])/(b_[hi]-b_[lo]);return h_[lo]+t*(h_[hi]-h_[lo]);}
double PiecewiseLinearBHCurve::reluctivity(double b) const {if(b<=1.0e-15)return h_[1]/b_[1];return field_strength(b)/b;}

double RectangularStrandedCoil2D::area_m2() const {if(!(xmax>xmin&&ymax>ymin)||!std::isfinite(turns)||!std::isfinite(current_a))throw std::invalid_argument("invalid rectangular stranded coil");return (xmax-xmin)*(ymax-ymin);}
bool RectangularStrandedCoil2D::contains(Node2 p) const noexcept{return p.x>=xmin&&p.x<=xmax&&p.y>=ymin&&p.y<=ymax;}
double RectangularStrandedCoil2D::current_density_a_per_m2(Node2 p) const {const double area=area_m2();return contains(p)?turns*current_a/area:0.0;}

double coil_flux_linkage_wb_turn(const Mesh2D& mesh,std::span<const double> a,const RectangularStrandedCoil2D& coil){mesh.validate();if(a.size()!=mesh.node_count())throw std::invalid_argument("coil flux-linkage field size mismatch");const double area=coil.area_m2();double integral=0.0;for(const auto& tri:mesh.triangles){const Node2 c{(mesh.nodes[tri.node[0]].x+mesh.nodes[tri.node[1]].x+mesh.nodes[tri.node[2]].x)/3.0,(mesh.nodes[tri.node[0]].y+mesh.nodes[tri.node[1]].y+mesh.nodes[tri.node[2]].y)/3.0};if(!coil.contains(c))continue;const double twice_area=std::abs((mesh.nodes[tri.node[1]].x-mesh.nodes[tri.node[0]].x)*(mesh.nodes[tri.node[2]].y-mesh.nodes[tri.node[0]].y)-(mesh.nodes[tri.node[1]].y-mesh.nodes[tri.node[0]].y)*(mesh.nodes[tri.node[2]].x-mesh.nodes[tri.node[0]].x));integral+=(a[tri.node[0]]+a[tri.node[1]]+a[tri.node[2]])/3.0*0.5*twice_area;}return coil.turns*integral/area;}
double coil_inductance_h(const Mesh2D& mesh,std::span<const double> a,const RectangularStrandedCoil2D& coil){if(std::abs(coil.current_a)<=1.0e-30)throw std::invalid_argument("coil inductance requires non-zero current");return coil_flux_linkage_wb_turn(mesh,a,coil)/coil.current_a;}

MaxwellForceTorque2D maxwell_stress_boundary_force_2d(const Mesh2D& mesh,std::span<const MagneticFluxDensity2> b,int patch,double mu,Node2 origin){mesh.validate();if(b.size()!=mesh.element_count()||!(mu>0.0))throw std::invalid_argument("invalid Maxwell-stress inputs");MaxwellForceTorque2D result;for(const auto& edge:mesh.boundary_edges){if(edge.patch!=patch)continue;std::size_t owner=mesh.element_count();for(std::size_t e=0;e<mesh.triangles.size();++e){bool first=false,second=false;for(auto node:mesh.triangles[e].node){first=first||node==edge.node[0];second=second||node==edge.node[1];}if(first&&second){owner=e;break;}}if(owner==mesh.element_count())throw std::runtime_error("boundary edge has no adjacent triangle");const auto&p0=mesh.nodes[edge.node[0]],&p1=mesh.nodes[edge.node[1]];const double dx=p1.x-p0.x,dy=p1.y-p0.y,length=std::hypot(dx,dy);double nx=dy/length,ny=-dx/length;const auto&t=mesh.triangles[owner];const Node2 centroid{(mesh.nodes[t.node[0]].x+mesh.nodes[t.node[1]].x+mesh.nodes[t.node[2]].x)/3.0,(mesh.nodes[t.node[0]].y+mesh.nodes[t.node[1]].y+mesh.nodes[t.node[2]].y)/3.0};const Node2 midpoint{0.5*(p0.x+p1.x),0.5*(p0.y+p1.y)};if(nx*(centroid.x-midpoint.x)+ny*(centroid.y-midpoint.y)>0.0){nx=-nx;ny=-ny;}const double bx=b[owner].x,by=b[owner].y,bn=bx*nx+by*ny,b2=bx*bx+by*by;const double tx=(bx*bn-0.5*b2*nx)/mu,ty=(by*bn-0.5*b2*ny)/mu,fx=tx*length,fy=ty*length;result.force_x_n_per_m+=fx;result.force_y_n_per_m+=fy;result.torque_z_n+=(midpoint.x-origin.x)*fy-(midpoint.y-origin.y)*fx;}return result;}

} // namespace cfd::fem
