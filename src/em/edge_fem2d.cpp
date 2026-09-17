#include "cfd/em/edge_fem2d.hpp"

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
constexpr std::array<std::array<std::size_t,2>,3> local_edges{{{{0U,1U}},{{1U,2U}},{{2U,0U}}}};

std::vector<Complex> solve_dense(std::vector<Complex> matrix,std::vector<Complex> rhs,std::size_t n){
    if(matrix.size()!=n*n||rhs.size()!=n)throw std::invalid_argument("edge-Maxwell dense solve size mismatch");
    for(std::size_t k=0;k<n;++k){std::size_t pivot=k;double best=std::abs(matrix[k*n+k]);for(std::size_t row=k+1U;row<n;++row){const double value=std::abs(matrix[row*n+k]);if(value>best){best=value;pivot=row;}}if(!(best>100.0*std::numeric_limits<double>::epsilon()))throw std::runtime_error("singular edge-Maxwell system");if(pivot!=k){for(std::size_t col=k;col<n;++col)std::swap(matrix[k*n+col],matrix[pivot*n+col]);std::swap(rhs[k],rhs[pivot]);}const Complex diagonal=matrix[k*n+k];for(std::size_t row=k+1U;row<n;++row){const Complex factor=matrix[row*n+k]/diagonal;if(std::abs(factor)==0.0)continue;matrix[row*n+k]={};for(std::size_t col=k+1U;col<n;++col)matrix[row*n+col]-=factor*matrix[k*n+col];rhs[row]-=factor*rhs[k];}}
    std::vector<Complex> x(n);for(std::size_t back=0;back<n;++back){const std::size_t row=n-1U-back;Complex value=rhs[row];for(std::size_t col=row+1U;col<n;++col)value-=matrix[row*n+col]*x[col];x[row]=value/matrix[row*n+row];}return x;
}

std::pair<std::array<std::size_t,2>,int> canonical_edge(std::size_t from,std::size_t to){
    if(from<to)return {{{from,to}},1};
    return {{{to,from}},-1};
}

cfd::fem::Node2 barycentric_point(const cfd::fem::Mesh2D& mesh,const cfd::fem::Tri3& triangle,std::array<double,3> lambda){
    cfd::fem::Node2 point{};for(std::size_t i=0;i<3U;++i){point.x+=lambda[i]*mesh.nodes[triangle.node[i]].x;point.y+=lambda[i]*mesh.nodes[triangle.node[i]].y;}return point;
}
}

NedelecTri3Geometry nedelec_tri3_geometry(const cfd::fem::Mesh2D& mesh,const cfd::fem::Tri3& triangle){
    const auto& a=mesh.nodes.at(triangle.node[0]);const auto& b=mesh.nodes.at(triangle.node[1]);const auto& c=mesh.nodes.at(triangle.node[2]);
    const double twice_signed_area=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);if(!(std::abs(twice_signed_area)>1.0e-30))throw std::invalid_argument("degenerate Nedelec triangle");
    NedelecTri3Geometry geometry;geometry.area=0.5*std::abs(twice_signed_area);
    geometry.gradient_lambda[0]={(b.y-c.y)/twice_signed_area,(c.x-b.x)/twice_signed_area};
    geometry.gradient_lambda[1]={(c.y-a.y)/twice_signed_area,(a.x-c.x)/twice_signed_area};
    geometry.gradient_lambda[2]={(a.y-b.y)/twice_signed_area,(b.x-a.x)/twice_signed_area};
    return geometry;
}

std::array<double,2> nedelec_tri3_basis(const NedelecTri3Geometry& geometry,std::size_t local_edge,std::array<double,3> lambda){
    if(local_edge>=3U)throw std::out_of_range("Nedelec local edge index");
    const auto i=local_edges[local_edge][0],j=local_edges[local_edge][1];
    return {lambda[i]*geometry.gradient_lambda[j][0]-lambda[j]*geometry.gradient_lambda[i][0],
            lambda[i]*geometry.gradient_lambda[j][1]-lambda[j]*geometry.gradient_lambda[i][1]};
}

double nedelec_tri3_curl(const NedelecTri3Geometry& geometry,std::size_t local_edge){
    if(local_edge>=3U)throw std::out_of_range("Nedelec local edge index");
    const auto i=local_edges[local_edge][0],j=local_edges[local_edge][1];
    const auto& gi=geometry.gradient_lambda[i];const auto& gj=geometry.gradient_lambda[j];
    return 2.0*(gi[0]*gj[1]-gi[1]*gj[0]);
}

EdgeMaxwell2DResult solve_driven_edge_maxwell_2d(const cfd::fem::Mesh2D& mesh,const EdgeMaxwell2DConfig& config,
                                                  const CurrentDensity2D& source){
    mesh.validate();if(!(config.frequency_hz>0.0)||!(config.relative_permittivity>0.0)||!(config.relative_permeability>0.0)||!(config.conductivity_s_per_m>=0.0)||!source)throw std::invalid_argument("invalid edge-Maxwell configuration/source");
    std::map<std::array<std::size_t,2>,std::size_t> edge_index;std::vector<std::array<std::size_t,2>> edges;
    std::vector<std::array<std::size_t,3>> triangle_edges(mesh.triangles.size());std::vector<std::array<int,3>> triangle_signs(mesh.triangles.size());
    for(std::size_t element=0;element<mesh.triangles.size();++element){const auto& tri=mesh.triangles[element];for(std::size_t local=0;local<3U;++local){const auto [edge,sign]=canonical_edge(tri.node[local_edges[local][0]],tri.node[local_edges[local][1]]);auto [it,inserted]=edge_index.emplace(edge,edges.size());if(inserted)edges.push_back(edge);triangle_edges[element][local]=it->second;triangle_signs[element][local]=sign;}}
    const std::size_t n=edges.size();std::vector<Complex> matrix(n*n,Complex{}),rhs(n,Complex{});const double omega=2.0*std::numbers::pi*config.frequency_hz,mu=mu0*config.relative_permeability,epsilon=epsilon0*config.relative_permittivity;
    constexpr std::array<std::array<double,3>,3> quadrature{{{{2.0/3.0,1.0/6.0,1.0/6.0}},{{1.0/6.0,2.0/3.0,1.0/6.0}},{{1.0/6.0,1.0/6.0,2.0/3.0}}}};
    for(std::size_t element=0;element<mesh.triangles.size();++element){const auto& tri=mesh.triangles[element];const auto geometry=nedelec_tri3_geometry(mesh,tri);std::array<std::array<double,3>,3> mass{};std::array<Complex,3> load{};
        for(const auto& lambda:quadrature){std::array<std::array<double,2>,3> basis{};for(std::size_t a=0;a<3U;++a)basis[a]=nedelec_tri3_basis(geometry,a,lambda);const auto point=barycentric_point(mesh,tri,lambda);const auto current=source(point);const double weight=geometry.area/3.0;for(std::size_t a=0;a<3U;++a){load[a]+=weight*(basis[a][0]*current[0]+basis[a][1]*current[1]);for(std::size_t b=0;b<3U;++b)mass[a][b]+=weight*(basis[a][0]*basis[b][0]+basis[a][1]*basis[b][1]);}}
        std::array<double,3> curl{};for(std::size_t a=0;a<3U;++a)curl[a]=nedelec_tri3_curl(geometry,a);
        for(std::size_t a=0;a<3U;++a){const std::size_t ga=triangle_edges[element][a];const double sa=static_cast<double>(triangle_signs[element][a]);rhs[ga]+=sa*Complex{0.0,-omega}*load[a];for(std::size_t b=0;b<3U;++b){const std::size_t gb=triangle_edges[element][b];const double sb=static_cast<double>(triangle_signs[element][b]);const Complex local=geometry.area*curl[a]*curl[b]/mu+Complex{-omega*omega*epsilon,omega*config.conductivity_s_per_m}*mass[a][b];matrix[ga*n+gb]+=sa*sb*local;}}
    }
    std::set<std::array<std::size_t,2>> boundary;for(const auto& item:mesh.boundary_edges)boundary.insert(canonical_edge(item.node[0],item.node[1]).first);
    for(const auto& edge:boundary){const std::size_t fixed=edge_index.at(edge);for(std::size_t column=0;column<n;++column){matrix[fixed*n+column]={};matrix[column*n+fixed]={};}matrix[fixed*n+fixed]=1.0;rhs[fixed]={};}
    EdgeMaxwell2DResult result;result.edges=edges;result.edge_voltage_v=solve_dense(std::move(matrix),std::move(rhs),n);result.electric_centroid_v_per_m.resize(mesh.triangles.size());constexpr std::array<double,3> centroid{{1.0/3.0,1.0/3.0,1.0/3.0}};
    for(std::size_t element=0;element<mesh.triangles.size();++element){const auto geometry=nedelec_tri3_geometry(mesh,mesh.triangles[element]);ComplexVec2 field{};for(std::size_t local=0;local<3U;++local){const auto basis=nedelec_tri3_basis(geometry,local,centroid);const Complex dof=static_cast<double>(triangle_signs[element][local])*result.edge_voltage_v[triangle_edges[element][local]];field[0]+=dof*basis[0];field[1]+=dof*basis[1];}result.electric_centroid_v_per_m[element]=field;}
    return result;
}

} // namespace cfd::em
