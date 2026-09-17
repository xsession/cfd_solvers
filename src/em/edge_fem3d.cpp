#include "cfd/em/edge_fem3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <utility>

namespace cfd::em {
namespace {
constexpr double epsilon0=8.8541878128e-12;
constexpr double mu0=1.25663706212e-6;
constexpr std::array<std::array<std::size_t,2>,6> local_edges{{
    {{0U,1U}},{{0U,2U}},{{0U,3U}},{{1U,2U}},{{1U,3U}},{{2U,3U}}
}};
constexpr std::array<std::array<std::size_t,3>,4> local_faces{{
    {{1U,2U,3U}},{{0U,2U,3U}},{{0U,1U,3U}},{{0U,1U,2U}}
}};

RealVec3 sub(cfd::fem::Point3 a,cfd::fem::Point3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
RealVec3 cross(RealVec3 a,RealVec3 b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double dot(RealVec3 a,RealVec3 b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm(RealVec3 a){return std::sqrt(dot(a,a));}
Complex cdot(ComplexVec3 a,RealVec3 b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double norm2(ComplexVec3 a){return std::norm(a[0])+std::norm(a[1])+std::norm(a[2]);}

std::pair<std::array<std::size_t,2>,int> canonical_edge(std::size_t a,std::size_t b){
    if(a<b)return {{{a,b}},1};
    return {{{b,a}},-1};
}
std::array<std::size_t,3> canonical_face(std::array<std::size_t,3> face){std::sort(face.begin(),face.end());return face;}

cfd::fem::Point3 barycentric_point(const cfd::fem::Mesh3D& mesh,const cfd::fem::Tet4& tet,
                                    const std::array<double,4>& lambda){
    cfd::fem::Point3 p{};
    for(std::size_t i=0;i<4U;++i){const auto& q=mesh.nodes[tet.node[i]];p.x+=lambda[i]*q.x;p.y+=lambda[i]*q.y;p.z+=lambda[i]*q.z;}
    return p;
}

struct Topology3D {
    std::vector<std::array<std::size_t,2>> edges;
    std::map<std::array<std::size_t,2>,std::size_t> edge_index;
    std::vector<std::array<std::size_t,6>> tet_edges;
    std::vector<std::array<int,6>> tet_signs;
    std::vector<unsigned char> boundary_edge;
    std::map<std::array<std::size_t,3>,std::pair<std::size_t,std::size_t>> boundary_faces; // face -> element, count
};

Topology3D build_topology(const cfd::fem::Mesh3D& mesh){
    Topology3D topology;topology.tet_edges.resize(mesh.tetrahedra.size());topology.tet_signs.resize(mesh.tetrahedra.size());
    std::map<std::array<std::size_t,3>,std::pair<std::size_t,std::size_t>> faces;
    for(std::size_t e=0;e<mesh.tetrahedra.size();++e){
        const auto& tet=mesh.tetrahedra[e];
        for(std::size_t le=0;le<6U;++le){const auto [edge,sign]=canonical_edge(tet.node[local_edges[le][0]],tet.node[local_edges[le][1]]);auto [it,inserted]=topology.edge_index.emplace(edge,topology.edges.size());if(inserted)topology.edges.push_back(edge);topology.tet_edges[e][le]=it->second;topology.tet_signs[e][le]=sign;}
        for(std::size_t lf=0;lf<4U;++lf){std::array<std::size_t,3> face{{tet.node[local_faces[lf][0]],tet.node[local_faces[lf][1]],tet.node[local_faces[lf][2]]}};face=canonical_face(face);auto [it,inserted]=faces.emplace(face,std::pair<std::size_t,std::size_t>{e,1U});if(!inserted)++it->second.second;}
    }
    topology.boundary_edge.assign(topology.edges.size(),0U);
    for(const auto& [face,owner_count]:faces){if(owner_count.second!=1U)continue;topology.boundary_faces.emplace(face,owner_count);for(std::size_t i=0;i<3U;++i)for(std::size_t j=i+1U;j<3U;++j)topology.boundary_edge[topology.edge_index.at(canonical_edge(face[i],face[j]).first)]=1U;}
    return topology;
}

ComplexVec3 scale_add(ComplexVec3 out,Complex value,RealVec3 basis){for(std::size_t d=0;d<3U;++d)out[d]+=value*basis[d];return out;}
}

NedelecTet4Geometry nedelec_tet4_geometry(const cfd::fem::Mesh3D& mesh,const cfd::fem::Tet4& tet){
    const auto& p0=mesh.nodes.at(tet.node[0]);const auto& p1=mesh.nodes.at(tet.node[1]);const auto& p2=mesh.nodes.at(tet.node[2]);const auto& p3=mesh.nodes.at(tet.node[3]);
    const RealVec3 v1=sub(p1,p0),v2=sub(p2,p0),v3=sub(p3,p0);
    const double det=dot(v1,cross(v2,v3));if(!(std::abs(det)>1.0e-30))throw std::invalid_argument("degenerate Nedelec tetrahedron");
    // Rows of M^{-1}, M=[v1 v2 v3], are the gradients of lambda_1..lambda_3.
    const RealVec3 g1=cross(v2,v3);const RealVec3 g2=cross(v3,v1);const RealVec3 g3=cross(v1,v2);
    NedelecTet4Geometry geometry;geometry.volume_m3=std::abs(det)/6.0;
    geometry.gradient_lambda[1]={g1[0]/det,g1[1]/det,g1[2]/det};
    geometry.gradient_lambda[2]={g2[0]/det,g2[1]/det,g2[2]/det};
    geometry.gradient_lambda[3]={g3[0]/det,g3[1]/det,g3[2]/det};
    for(std::size_t d=0;d<3U;++d)geometry.gradient_lambda[0][d]=-(geometry.gradient_lambda[1][d]+geometry.gradient_lambda[2][d]+geometry.gradient_lambda[3][d]);
    return geometry;
}

RealVec3 nedelec_tet4_basis(const NedelecTet4Geometry& geometry,std::size_t local_edge,std::array<double,4> lambda){
    if(local_edge>=6U)throw std::out_of_range("Nedelec Tet4 local edge index");
    const auto i=local_edges[local_edge][0],j=local_edges[local_edge][1];
    RealVec3 value{};
    for(std::size_t d=0;d<3U;++d)value[d]=lambda[i]*geometry.gradient_lambda[j][d]-lambda[j]*geometry.gradient_lambda[i][d];
    return value;
}

RealVec3 nedelec_tet4_curl(const NedelecTet4Geometry& geometry,std::size_t local_edge){
    if(local_edge>=6U)throw std::out_of_range("Nedelec Tet4 local edge index");
    const auto i=local_edges[local_edge][0],j=local_edges[local_edge][1];
    auto value=cross(geometry.gradient_lambda[i],geometry.gradient_lambda[j]);
    for(double& component:value)component*=2.0;
    return value;
}

EdgeMaxwell3DResult solve_driven_edge_maxwell_3d(const cfd::fem::Mesh3D& mesh,const EdgeMaxwell3DConfig& config,
                                                  const CurrentDensity3D& source,std::span<const EdgeCurrentPort3D> ports){
    mesh.validate();if(!(config.frequency_hz>0.0)||!(config.relative_permittivity>0.0)||!(config.relative_permeability>0.0)||!(config.conductivity_s_per_m>=0.0))throw std::invalid_argument("invalid 3-D edge-Maxwell configuration");
    const auto topology=build_topology(mesh);std::vector<std::size_t> free_index(topology.edges.size(),static_cast<std::size_t>(-1));std::size_t free_count=0U;for(std::size_t e=0;e<topology.edges.size();++e)if(!topology.boundary_edge[e])free_index[e]=free_count++;
    if(free_count==0U)throw std::runtime_error("3-D edge-Maxwell mesh has no free interior edges");
    cfd::core::ComplexCsrBuilder matrix(free_count);std::vector<Complex> rhs(free_count,Complex{});
    const double omega=2.0*std::numbers::pi*config.frequency_hz,mu=mu0*config.relative_permeability,epsilon=epsilon0*config.relative_permittivity;
    constexpr double a=0.5854101966249685,b=0.1381966011250105;
    constexpr std::array<std::array<double,4>,4> quadrature{{{{a,b,b,b}},{{b,a,b,b}},{{b,b,a,b}},{{b,b,b,a}}}};
    for(std::size_t element=0;element<mesh.tetrahedra.size();++element){
        const auto& tet=mesh.tetrahedra[element];const auto geometry=nedelec_tet4_geometry(mesh,tet);std::array<std::array<double,6>,6> mass{};std::array<Complex,6> load{};
        for(const auto& lambda:quadrature){std::array<RealVec3,6> basis{};for(std::size_t le=0;le<6U;++le)basis[le]=nedelec_tet4_basis(geometry,le,lambda);const double weight=geometry.volume_m3/4.0;ComplexVec3 current{};if(source)current=source(barycentric_point(mesh,tet,lambda));for(std::size_t i=0;i<6U;++i){load[i]+=weight*cdot(current,basis[i]);for(std::size_t j=0;j<6U;++j)mass[i][j]+=weight*dot(basis[i],basis[j]);}}
        std::array<RealVec3,6> curl{};for(std::size_t le=0;le<6U;++le)curl[le]=nedelec_tet4_curl(geometry,le);
        for(std::size_t i=0;i<6U;++i){const auto gi=topology.tet_edges[element][i];if(topology.boundary_edge[gi])continue;const std::size_t ri=free_index[gi];const double si=static_cast<double>(topology.tet_signs[element][i]);rhs[ri]+=si*Complex{0.0,-omega}*load[i];for(std::size_t j=0;j<6U;++j){const auto gj=topology.tet_edges[element][j];if(topology.boundary_edge[gj])continue;const double sj=static_cast<double>(topology.tet_signs[element][j]);const Complex local=geometry.volume_m3*dot(curl[i],curl[j])/mu+Complex{-omega*omega*epsilon,omega*config.conductivity_s_per_m}*mass[i][j];matrix.add(ri,free_index[gj],si*sj*local);}}
    }
    for(const auto& port:ports){if(port.node_a==port.node_b||port.node_a>=mesh.node_count()||port.node_b>=mesh.node_count())throw std::invalid_argument("invalid 3-D edge current port");const auto [edge,sign]=canonical_edge(port.node_a,port.node_b);const auto it=topology.edge_index.find(edge);if(it==topology.edge_index.end())throw std::invalid_argument("3-D edge current port is not a mesh edge");if(topology.boundary_edge[it->second])throw std::invalid_argument("3-D edge current port cannot excite a PEC boundary edge");rhs[free_index[it->second]]+=static_cast<double>(sign)*Complex{0.0,-omega}*port.current_a;}
    const auto solved=cfd::core::solve_complex_sparse(matrix,rhs,config.linear);if(!solved.linear_result.converged)throw std::runtime_error("3-D edge-Maxwell GMRES failed to converge");
    EdgeMaxwell3DResult result;result.edges=topology.edges;result.boundary_edge=topology.boundary_edge;result.edge_voltage_v.assign(topology.edges.size(),Complex{});result.linear_result=solved.linear_result;for(std::size_t e=0;e<topology.edges.size();++e)if(!topology.boundary_edge[e])result.edge_voltage_v[e]=solved.solution[free_index[e]];
    result.electric_centroid_v_per_m.resize(mesh.tetrahedra.size());result.magnetic_centroid_a_per_m.resize(mesh.tetrahedra.size());constexpr std::array<double,4> centroid{{0.25,0.25,0.25,0.25}};
    for(std::size_t element=0;element<mesh.tetrahedra.size();++element){const auto geometry=nedelec_tet4_geometry(mesh,mesh.tetrahedra[element]);ComplexVec3 electric{},curl_e{};for(std::size_t le=0;le<6U;++le){const Complex dof=static_cast<double>(topology.tet_signs[element][le])*result.edge_voltage_v[topology.tet_edges[element][le]];electric=scale_add(electric,dof,nedelec_tet4_basis(geometry,le,centroid));curl_e=scale_add(curl_e,dof,nedelec_tet4_curl(geometry,le));}result.electric_centroid_v_per_m[element]=electric;const Complex factor{0.0,1.0/(omega*mu)};for(std::size_t d=0;d<3U;++d)result.magnetic_centroid_a_per_m[element][d]=factor*curl_e[d];}
    return result;
}

ResonatorQuality3D resonator_quality_3d(const cfd::fem::Mesh3D& mesh,const EdgeMaxwell3DConfig& config,
                                        const EdgeMaxwell3DResult& field,double wall_sigma,double wall_mu_r){
    mesh.validate();if(field.electric_centroid_v_per_m.size()!=mesh.element_count()||field.magnetic_centroid_a_per_m.size()!=mesh.element_count())throw std::invalid_argument("3-D resonator field/mesh size mismatch");if(!(config.frequency_hz>0.0)||!(config.relative_permittivity>0.0)||!(config.relative_permeability>0.0)||!(config.conductivity_s_per_m>=0.0)||!(wall_sigma>=0.0)||!(wall_mu_r>0.0))throw std::invalid_argument("invalid resonator quality controls");
    const double epsilon=epsilon0*config.relative_permittivity,mu=mu0*config.relative_permeability,omega=2.0*std::numbers::pi*config.frequency_hz;ResonatorQuality3D out;
    for(std::size_t e=0;e<mesh.element_count();++e){const double volume=nedelec_tet4_geometry(mesh,mesh.tetrahedra[e]).volume_m3;const double e2=norm2(field.electric_centroid_v_per_m[e]),h2=norm2(field.magnetic_centroid_a_per_m[e]);out.electric_energy_j+=0.25*epsilon*e2*volume;out.magnetic_energy_j+=0.25*mu*h2*volume;out.dielectric_loss_w+=0.5*config.conductivity_s_per_m*e2*volume;}
    if(wall_sigma>0.0){const auto topology=build_topology(mesh);const double rs=std::sqrt(std::numbers::pi*config.frequency_hz*mu0*wall_mu_r/wall_sigma);for(const auto& [face,owner_count]:topology.boundary_faces){const auto& p0=mesh.nodes[face[0]],&p1=mesh.nodes[face[1]],&p2=mesh.nodes[face[2]];const RealVec3 normal_raw=cross(sub(p1,p0),sub(p2,p0));const double twice_area=norm(normal_raw);if(!(twice_area>0.0))continue;RealVec3 normal=normal_raw;for(double& v:normal)v/=twice_area;const auto& h=field.magnetic_centroid_a_per_m[owner_count.first];const Complex hn=h[0]*normal[0]+h[1]*normal[1]+h[2]*normal[2];ComplexVec3 ht=h;for(std::size_t d=0;d<3U;++d)ht[d]-=hn*normal[d];out.conductor_loss_w+=0.5*rs*norm2(ht)*(0.5*twice_area);}}
    const double loss=out.dielectric_loss_w+out.conductor_loss_w,stored=out.electric_energy_j+out.magnetic_energy_j;out.quality_factor=loss>0.0?omega*stored/loss:std::numeric_limits<double>::infinity();return out;
}

} // namespace cfd::em
