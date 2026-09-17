#include "cfd/em/adaptivity3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace cfd::em {
namespace {
using Face=std::array<std::size_t,3>;
using Edge=std::array<std::size_t,2>;
constexpr std::array<std::array<std::size_t,3>,4> local_faces{{{{1U,2U,3U}},{{0U,2U,3U}},{{0U,1U,3U}},{{0U,1U,2U}}}};
constexpr std::array<std::array<std::size_t,2>,6> local_edges{{{{0U,1U}},{{0U,2U}},{{0U,3U}},{{1U,2U}},{{1U,3U}},{{2U,3U}}}};
Face canonical_face(Face f){std::sort(f.begin(),f.end());return f;}
Edge canonical_edge(std::size_t a,std::size_t b){return a<b?Edge{{a,b}}:Edge{{b,a}};}
std::array<double,3> sub(cfd::fem::Point3 a,cfd::fem::Point3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
std::array<double,3> cross(std::array<double,3> a,std::array<double,3> b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(std::array<double,3> a,std::array<double,3> b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(std::array<double,3> a){return std::sqrt(dot(a,a));}
double distance2(cfd::fem::Point3 a,cfd::fem::Point3 b){const auto d=sub(a,b);return dot(d,d);}
double tangential_jump2(ComplexVec3 a,ComplexVec3 b,std::array<double,3> n){ComplexVec3 d{{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};const Complex dn=d[0]*n[0]+d[1]*n[1]+d[2]*n[2];for(std::size_t k=0;k<3U;++k)d[k]-=dn*n[k];return std::norm(d[0])+std::norm(d[1])+std::norm(d[2]);}

std::map<Face,std::vector<std::size_t>> face_owners(const cfd::fem::Mesh3D& mesh){std::map<Face,std::vector<std::size_t>> owners;for(std::size_t e=0;e<mesh.tetrahedra.size();++e){const auto& t=mesh.tetrahedra[e];for(const auto& lf:local_faces){Face face{{t.node[lf[0]],t.node[lf[1]],t.node[lf[2]]}};owners[canonical_face(face)].push_back(e);}}return owners;}
std::set<Edge> boundary_edges(const cfd::fem::Mesh3D& mesh){std::set<Edge> edges;for(const auto& [face,owners]:face_owners(mesh))if(owners.size()==1U)for(std::size_t i=0;i<3U;++i)for(std::size_t j=i+1U;j<3U;++j)edges.insert(canonical_edge(face[i],face[j]));return edges;}
}

std::vector<double> maxwell_face_jump_indicators_3d(const cfd::fem::Mesh3D& mesh,const EdgeMaxwell3DResult& field,double impedance){
    mesh.validate();if(field.electric_centroid_v_per_m.size()!=mesh.element_count()||field.magnetic_centroid_a_per_m.size()!=mesh.element_count()||!(impedance>0.0))throw std::invalid_argument("invalid Maxwell jump-estimator inputs");std::vector<double> indicator(mesh.element_count(),0.0);
    for(const auto& [face,owners]:face_owners(mesh)){if(owners.size()!=2U)continue;const auto& p0=mesh.nodes[face[0]],&p1=mesh.nodes[face[1]],&p2=mesh.nodes[face[2]];auto normal=cross(sub(p1,p0),sub(p2,p0));const double twice_area=norm(normal);if(!(twice_area>0.0))continue;for(double& value:normal)value/=twice_area;const double e_jump=tangential_jump2(field.electric_centroid_v_per_m[owners[0]],field.electric_centroid_v_per_m[owners[1]],normal);const double h_jump=tangential_jump2(field.magnetic_centroid_a_per_m[owners[0]],field.magnetic_centroid_a_per_m[owners[1]],normal);const double eta=0.5*twice_area*(e_jump+impedance*impedance*h_jump);indicator[owners[0]]+=0.5*eta;indicator[owners[1]]+=0.5*eta;}
    for(double& value:indicator)value=std::sqrt(std::max(0.0,value));return indicator;
}

std::vector<unsigned char> mark_maxwell_dorfler(std::span<const double> indicators,double fraction){
    if(indicators.empty()||!(fraction>0.0&&fraction<=1.0))throw std::invalid_argument("invalid Maxwell Dorfler controls");double total=0.0;for(double value:indicators){if(!(value>=0.0)||!std::isfinite(value))throw std::invalid_argument("invalid Maxwell error indicator");total+=value*value;}std::vector<unsigned char> marked(indicators.size(),0U);if(!(total>0.0))return marked;std::vector<std::size_t> order(indicators.size());std::iota(order.begin(),order.end(),0U);std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return indicators[a]>indicators[b];});double selected=0.0;for(auto i:order){marked[i]=1U;selected+=indicators[i]*indicators[i];if(selected>=fraction*total)break;}return marked;
}

cfd::fem::Mesh3D refine_tet4_marked_longest_edges(const cfd::fem::Mesh3D& mesh,std::span<const unsigned char> marked){
    mesh.validate();if(marked.size()!=mesh.element_count())throw std::invalid_argument("Tet4 refinement marking size mismatch");std::set<Edge> selected;for(std::size_t e=0;e<mesh.element_count();++e)if(marked[e]){const auto& tet=mesh.tetrahedra[e];double best=-1.0;Edge longest{};for(const auto& le:local_edges){const Edge edge=canonical_edge(tet.node[le[0]],tet.node[le[1]]);const double length2=distance2(mesh.nodes[edge[0]],mesh.nodes[edge[1]]);if(length2>best){best=length2;longest=edge;}}selected.insert(longest);}if(selected.empty())return mesh;
    cfd::fem::Mesh3D out=mesh;std::map<Edge,std::size_t> midpoint;
    for(const auto& edge:selected){std::vector<std::size_t> incident;for(std::size_t e=0;e<out.tetrahedra.size();++e){bool a=false,b=false;for(auto node:out.tetrahedra[e].node){a=a||node==edge[0];b=b||node==edge[1];}if(a&&b)incident.push_back(e);}if(incident.empty())continue;const auto boundary=boundary_edges(out);std::size_t mid{};auto it=midpoint.find(edge);if(it==midpoint.end()){const auto& a=out.nodes[edge[0]],&b=out.nodes[edge[1]];mid=out.nodes.size();out.nodes.push_back({0.5*(a.x+b.x),0.5*(a.y+b.y),0.5*(a.z+b.z)});out.boundary_node.push_back(static_cast<unsigned char>(boundary.contains(edge)));midpoint.emplace(edge,mid);}else mid=it->second;
        std::vector<cfd::fem::Tet4> rebuilt;rebuilt.reserve(out.tetrahedra.size()+incident.size());std::vector<unsigned char> split(out.tetrahedra.size(),0U);for(auto e:incident)split[e]=1U;for(std::size_t e=0;e<out.tetrahedra.size();++e){if(!split[e]){rebuilt.push_back(out.tetrahedra[e]);continue;}const auto& t=out.tetrahedra[e];std::array<std::size_t,2> other{};std::size_t count=0U;for(auto node:t.node)if(node!=edge[0]&&node!=edge[1])other[count++]=node;if(count!=2U)throw std::runtime_error("invalid Tet4 edge-star refinement");rebuilt.push_back({{edge[0],mid,other[0],other[1]}});rebuilt.push_back({{mid,edge[1],other[0],other[1]}});}out.tetrahedra=std::move(rebuilt);
    }
    out.validate();return out;
}

} // namespace cfd::em
