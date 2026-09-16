#include "cfd/fem/adaptivity.hpp"

#include "cfd/fem/reference_element.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace cfd::fem {
namespace {

using Edge = std::pair<std::size_t, std::size_t>;

Edge edge_key(std::size_t a, std::size_t b) { return a < b ? Edge{a,b} : Edge{b,a}; }

double signed_area(Node2 a, Node2 b, Node2 c) {
    return 0.5 * ((b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x));
}

double edge_length(Node2 a, Node2 b) { return std::hypot(b.x-a.x,b.y-a.y); }

std::array<double,2> triangle_gradient(const Mesh2D& mesh, const Tri3& tri, std::span<const double> u) {
    const auto&a=mesh.nodes[tri.node[0]],&b=mesh.nodes[tri.node[1]],&c=mesh.nodes[tri.node[2]];
    const double det=(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    if (std::abs(det) <= 1.0e-30) throw std::runtime_error("degenerate triangle in error estimator");
    const double inv=1.0/det;
    const std::array<std::array<double,2>,3> grad{{
        {{(b.y-c.y)*inv,(c.x-b.x)*inv}},
        {{(c.y-a.y)*inv,(a.x-c.x)*inv}},
        {{(a.y-b.y)*inv,(b.x-a.x)*inv}}}};
    std::array<double,2> g{};
    for(std::size_t i=0;i<3U;++i){g[0]+=u[tri.node[i]]*grad[i][0];g[1]+=u[tri.node[i]]*grad[i][1];}
    return g;
}

void add_oriented_triangle(Mesh2D& out, std::size_t a, std::size_t b, std::size_t c) {
    if (std::abs(signed_area(out.nodes[a],out.nodes[b],out.nodes[c])) <= 1.0e-30) {
        throw std::runtime_error("adaptive refinement produced degenerate triangle");
    }
    if (signed_area(out.nodes[a],out.nodes[b],out.nodes[c]) < 0.0) std::swap(b,c);
    out.triangles.push_back({{{a,b,c}}});
}

} // namespace

PoissonErrorEstimate2D estimate_poisson_error_tri3(const Mesh2D& mesh,
                                                    std::span<const double> nodal_solution,
                                                    const std::function<double(Node2)>& source,
                                                    double diffusivity) {
    mesh.validate();
    if (nodal_solution.size()!=mesh.node_count() || !source || !(diffusivity>0.0)) {
        throw std::invalid_argument("invalid Poisson estimator input");
    }
    const std::size_t ne=mesh.triangles.size();
    PoissonErrorEstimate2D result;
    result.element_indicator.assign(ne,0.0);
    std::vector<double> eta2(ne,0.0);
    std::vector<std::array<double,2>> grad(ne);
    struct Adj {std::size_t triangle{}; std::size_t a{}; std::size_t b{};};
    std::map<Edge,std::vector<Adj>> adjacency;

    for(std::size_t e=0;e<ne;++e){
        const auto&tri=mesh.triangles[e];
        grad[e]=triangle_gradient(mesh,tri,nodal_solution);
        const Node2 a=mesh.nodes[tri.node[0]],b=mesh.nodes[tri.node[1]],c=mesh.nodes[tri.node[2]];
        const double area=std::abs(signed_area(a,b,c));
        const double h=std::max({edge_length(a,b),edge_length(b,c),edge_length(c,a)});
        const Node2 center{(a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0};
        const double f=source(center);
        eta2[e]+=h*h*area*f*f;
        adjacency[edge_key(tri.node[0],tri.node[1])].push_back({e,tri.node[0],tri.node[1]});
        adjacency[edge_key(tri.node[1],tri.node[2])].push_back({e,tri.node[1],tri.node[2]});
        adjacency[edge_key(tri.node[2],tri.node[0])].push_back({e,tri.node[2],tri.node[0]});
    }

    for(const auto& item:adjacency){
        const auto& adj=item.second;
        if(adj.size()!=2U) continue;
        const auto p0=mesh.nodes[item.first.first],p1=mesh.nodes[item.first.second];
        const double dx=p1.x-p0.x,dy=p1.y-p0.y;
        const double length=std::hypot(dx,dy);
        if(!(length>0.0)) continue;
        const double nx=dy/length,ny=-dx/length;
        const auto&t0=grad[adj[0].triangle];
        const auto&t1=grad[adj[1].triangle];
        const double jump=diffusivity*((t0[0]-t1[0])*nx+(t0[1]-t1[1])*ny);
        const double contribution=0.5*length*jump*jump;
        eta2[adj[0].triangle]+=contribution;
        eta2[adj[1].triangle]+=contribution;
    }

    double total=0.0;
    for(std::size_t e=0;e<ne;++e){
        result.element_indicator[e]=std::sqrt(std::max(0.0,eta2[e]));
        total+=eta2[e];
    }
    result.global_indicator=std::sqrt(total);
    return result;
}

std::vector<unsigned char> mark_dorfler(std::span<const double> indicator, double target_fraction) {
    if(indicator.empty() || !(target_fraction>0.0) || target_fraction>1.0) {
        throw std::invalid_argument("invalid Dorfler marking input");
    }
    std::vector<std::size_t> order(indicator.size());
    std::iota(order.begin(),order.end(),0U);
    std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return indicator[a]>indicator[b];});
    double total=0.0; for(double v:indicator){if(!(v>=0.0)||!std::isfinite(v)) throw std::invalid_argument("invalid error indicator"); total+=v*v;}
    std::vector<unsigned char> marked(indicator.size(),0U);
    if(total==0.0){marked[order.front()]=1U; return marked;}
    const double target=target_fraction*total;
    double accumulated=0.0;
    for(const auto e:order){marked[e]=1U; accumulated+=indicator[e]*indicator[e]; if(accumulated>=target) break;}
    return marked;
}

Mesh2D refine_tri3_longest_edges(const Mesh2D& mesh, std::span<const unsigned char> marked_element) {
    mesh.validate();
    if(marked_element.size()!=mesh.element_count()) throw std::invalid_argument("adaptive mark size mismatch");
    std::set<Edge> split_edges;
    for(std::size_t e=0;e<mesh.triangles.size();++e){
        if(!marked_element[e]) continue;
        const auto&tri=mesh.triangles[e];
        const std::array<Edge,3> edge{{edge_key(tri.node[0],tri.node[1]),edge_key(tri.node[1],tri.node[2]),edge_key(tri.node[2],tri.node[0])}};
        std::array<double,3> len{{edge_length(mesh.nodes[edge[0].first],mesh.nodes[edge[0].second]),
                                  edge_length(mesh.nodes[edge[1].first],mesh.nodes[edge[1].second]),
                                  edge_length(mesh.nodes[edge[2].first],mesh.nodes[edge[2].second])}};
        const auto it=std::max_element(len.begin(),len.end());
        split_edges.insert(edge[static_cast<std::size_t>(std::distance(len.begin(),it))]);
    }
    if(split_edges.empty()) return mesh;

    Mesh2D out;
    out.nodes=mesh.nodes;
    std::map<Edge,std::size_t> midpoint;
    for(const auto&e:split_edges){
        const auto a=mesh.nodes[e.first],b=mesh.nodes[e.second];
        midpoint[e]=out.nodes.size();
        out.nodes.push_back({0.5*(a.x+b.x),0.5*(a.y+b.y)});
    }

    for(const auto&tri:mesh.triangles){
        const std::size_t a=tri.node[0],b=tri.node[1],c=tri.node[2];
        const Edge e0=edge_key(a,b),e1=edge_key(b,c),e2=edge_key(c,a);
        const bool s0=midpoint.contains(e0),s1=midpoint.contains(e1),s2=midpoint.contains(e2);
        const unsigned count=static_cast<unsigned>(s0)+static_cast<unsigned>(s1)+static_cast<unsigned>(s2);
        if(count==0U){add_oriented_triangle(out,a,b,c); continue;}
        const std::size_t m0=s0?midpoint[e0]:0U,m1=s1?midpoint[e1]:0U,m2=s2?midpoint[e2]:0U;
        if(count==1U){
            if(s0){add_oriented_triangle(out,a,m0,c);add_oriented_triangle(out,m0,b,c);}
            else if(s1){add_oriented_triangle(out,b,m1,a);add_oriented_triangle(out,m1,c,a);}
            else {add_oriented_triangle(out,c,m2,b);add_oriented_triangle(out,m2,a,b);}
        } else if(count==2U){
            if(s0&&s1){add_oriented_triangle(out,b,m1,m0);add_oriented_triangle(out,m0,m1,c);add_oriented_triangle(out,a,m0,c);}
            else if(s1&&s2){add_oriented_triangle(out,c,m2,m1);add_oriented_triangle(out,m1,m2,a);add_oriented_triangle(out,b,m1,a);}
            else {add_oriented_triangle(out,a,m0,m2);add_oriented_triangle(out,m2,m0,b);add_oriented_triangle(out,c,m2,b);}
        } else {
            add_oriented_triangle(out,a,m0,m2);add_oriented_triangle(out,m0,b,m1);add_oriented_triangle(out,m2,m1,c);add_oriented_triangle(out,m0,m1,m2);
        }
    }

    for(const auto&edge:mesh.boundary_edges){
        const Edge key=edge_key(edge.node[0],edge.node[1]);
        const auto it=midpoint.find(key);
        if(it==midpoint.end()) out.boundary_edges.push_back(edge);
        else {
            out.boundary_edges.push_back({{{edge.node[0],it->second}},edge.patch});
            out.boundary_edges.push_back({{{it->second,edge.node[1]}},edge.patch});
        }
    }
    out.boundary_node.assign(out.nodes.size(),0U);
    for(const auto&edge:out.boundary_edges){out.boundary_node[edge.node[0]]=1U;out.boundary_node[edge.node[1]]=1U;}
    out.validate();
    return out;
}

} // namespace cfd::fem
