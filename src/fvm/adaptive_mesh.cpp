#include "cfd/fvm/adaptive_mesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace cfd::fvm {
namespace {

constexpr double eps = 1.0e-12;

bool finite(Vec3 p) {
    return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
}

double cell_volume(const AdaptiveHexCell& c) {
    return (c.upper.x-c.lower.x)*(c.upper.y-c.lower.y)*(c.upper.z-c.lower.z);
}

Vec3 cell_center(const AdaptiveHexCell& c) {
    return {(c.lower.x+c.upper.x)*0.5,(c.lower.y+c.upper.y)*0.5,(c.lower.z+c.upper.z)*0.5};
}

bool near(double a,double b) {
    return std::abs(a-b) <= eps*std::max({1.0,std::abs(a),std::abs(b)});
}

double overlap_1d(double a0,double a1,double b0,double b1) {
    return std::max(0.0,std::min(a1,b1)-std::max(a0,b0));
}

double overlap_volume(const AdaptiveHexCell& a,const AdaptiveHexCell& b) {
    return overlap_1d(a.lower.x,a.upper.x,b.lower.x,b.upper.x)
          *overlap_1d(a.lower.y,a.upper.y,b.lower.y,b.upper.y)
          *overlap_1d(a.lower.z,a.upper.z,b.lower.z,b.upper.z);
}

struct FaceContact {
    int axis{-1};
    double plane{};
    double lo0{};
    double hi0{};
    double lo1{};
    double hi1{};
    int direction{}; // +1 when b lies in positive axis direction from a
};

FaceContact face_contact(const AdaptiveHexCell& a,const AdaptiveHexCell& b) {
    FaceContact c;
    if (near(a.upper.x,b.lower.x) || near(b.upper.x,a.lower.x)) {
        const double oy=overlap_1d(a.lower.y,a.upper.y,b.lower.y,b.upper.y);
        const double oz=overlap_1d(a.lower.z,a.upper.z,b.lower.z,b.upper.z);
        if (oy>eps && oz>eps) {
            c.axis=0;c.direction=near(a.upper.x,b.lower.x)?1:-1;
            c.plane=c.direction>0?a.upper.x:a.lower.x;
            c.lo0=std::max(a.lower.y,b.lower.y);c.hi0=std::min(a.upper.y,b.upper.y);
            c.lo1=std::max(a.lower.z,b.lower.z);c.hi1=std::min(a.upper.z,b.upper.z);
            return c;
        }
    }
    if (near(a.upper.y,b.lower.y) || near(b.upper.y,a.lower.y)) {
        const double ox=overlap_1d(a.lower.x,a.upper.x,b.lower.x,b.upper.x);
        const double oz=overlap_1d(a.lower.z,a.upper.z,b.lower.z,b.upper.z);
        if (ox>eps && oz>eps) {
            c.axis=1;c.direction=near(a.upper.y,b.lower.y)?1:-1;
            c.plane=c.direction>0?a.upper.y:a.lower.y;
            c.lo0=std::max(a.lower.x,b.lower.x);c.hi0=std::min(a.upper.x,b.upper.x);
            c.lo1=std::max(a.lower.z,b.lower.z);c.hi1=std::min(a.upper.z,b.upper.z);
            return c;
        }
    }
    if (near(a.upper.z,b.lower.z) || near(b.upper.z,a.lower.z)) {
        const double ox=overlap_1d(a.lower.x,a.upper.x,b.lower.x,b.upper.x);
        const double oy=overlap_1d(a.lower.y,a.upper.y,b.lower.y,b.upper.y);
        if (ox>eps && oy>eps) {
            c.axis=2;c.direction=near(a.upper.z,b.lower.z)?1:-1;
            c.plane=c.direction>0?a.upper.z:a.lower.z;
            c.lo0=std::max(a.lower.x,b.lower.x);c.hi0=std::min(a.upper.x,b.upper.x);
            c.lo1=std::max(a.lower.y,b.lower.y);c.hi1=std::min(a.upper.y,b.upper.y);
            return c;
        }
    }
    return c;
}

std::array<AdaptiveHexCell,8> split_cell(const AdaptiveHexCell& c) {
    const Vec3 m=cell_center(c);
    std::array<AdaptiveHexCell,8> out{};
    for (std::size_t child=0;child<8;++child) {
        const bool hx=(child&1U)!=0U,hy=(child&2U)!=0U,hz=(child&4U)!=0U;
        auto& q=out[child];
        q.lower={hx?m.x:c.lower.x,hy?m.y:c.lower.y,hz?m.z:c.lower.z};
        q.upper={hx?c.upper.x:m.x,hy?c.upper.y:m.y,hz?c.upper.z:m.z};
        q.level=c.level+1U;q.root=c.root;q.lineage=(c.lineage<<3U)|child;
    }
    return out;
}

AdaptiveHexCell merge_siblings(std::span<const AdaptiveHexCell> siblings) {
    if (siblings.size()!=8U) throw std::invalid_argument("hex coarsening requires eight siblings");
    AdaptiveHexCell parent=siblings.front();
    parent.lower=siblings.front().lower;parent.upper=siblings.front().upper;
    for (const auto& c:siblings) {
        parent.lower.x=std::min(parent.lower.x,c.lower.x);parent.lower.y=std::min(parent.lower.y,c.lower.y);parent.lower.z=std::min(parent.lower.z,c.lower.z);
        parent.upper.x=std::max(parent.upper.x,c.upper.x);parent.upper.y=std::max(parent.upper.y,c.upper.y);parent.upper.z=std::max(parent.upper.z,c.upper.z);
    }
    parent.level=siblings.front().level-1U;
    parent.root=siblings.front().root;
    parent.lineage=siblings.front().lineage>>3U;
    return parent;
}

AdaptiveHexMesh refine_once(const AdaptiveHexMesh& mesh,std::span<const unsigned char> marked,std::size_t max_level) {
    if (marked.size()!=mesh.cell_count()) throw std::invalid_argument("adaptive refinement mark size mismatch");
    std::vector<AdaptiveHexCell> out;
    out.reserve(mesh.cell_count()*2U);
    for (std::size_t i=0;i<mesh.cell_count();++i) {
        const auto& c=mesh.cells()[i];
        if (marked[i] && c.level<max_level) {
            const auto children=split_cell(c);
            out.insert(out.end(),children.begin(),children.end());
        } else out.push_back(c);
    }
    return AdaptiveHexMesh(std::move(out),mesh.domain_lower(),mesh.domain_upper());
}

AdaptiveHexMesh enforce_balance(AdaptiveHexMesh mesh,std::size_t max_level) {
    for (std::size_t pass=0;pass<64U;++pass) {
        std::vector<unsigned char> mark(mesh.cell_count(),0U);
        bool any=false;
        for (std::size_t i=0;i<mesh.cell_count();++i) for (std::size_t j=i+1;j<mesh.cell_count();++j) {
            if (face_contact(mesh.cells()[i],mesh.cells()[j]).axis<0) continue;
            const auto li=mesh.cells()[i].level,lj=mesh.cells()[j].level;
            if (li+1U<lj) { if (li>=max_level) throw std::runtime_error("2:1 balance requires refinement above maximum level"); mark[i]=1U;any=true; }
            else if (lj+1U<li) { if (lj>=max_level) throw std::runtime_error("2:1 balance requires refinement above maximum level"); mark[j]=1U;any=true; }
        }
        if (!any) return mesh;
        mesh=refine_once(mesh,mark,max_level);
    }
    throw std::runtime_error("2:1 balancing did not converge");
}

AdaptiveHexMesh coarsen_once(const AdaptiveHexMesh& mesh,std::span<const unsigned char> marked,std::size_t min_level) {
    if (marked.size()!=mesh.cell_count()) throw std::invalid_argument("adaptive coarsening mark size mismatch");
    using Key=std::tuple<std::size_t,std::size_t,std::uint64_t>;
    std::map<Key,std::vector<std::size_t>> groups;
    for (std::size_t i=0;i<mesh.cell_count();++i) {
        const auto& c=mesh.cells()[i];
        if (c.level>min_level) groups[{c.root,c.level,c.lineage>>3U}].push_back(i);
    }
    std::vector<unsigned char> consumed(mesh.cell_count(),0U);
    std::vector<AdaptiveHexCell> parents;
    for (const auto& [key,indices]:groups) {
        if (indices.size()!=8U) continue;
        bool all=true;std::array<unsigned char,8> child_seen{};
        for (const auto idx:indices) {
            all=all&&marked[idx];
            const auto child=static_cast<std::size_t>(mesh.cells()[idx].lineage&7U);
            child_seen[child]=1U;
        }
        for (auto seen:child_seen) all=all&&seen;
        if (!all) continue;
        std::array<AdaptiveHexCell,8> siblings{};
        for (std::size_t k=0;k<indices.size();++k) { siblings[k]=mesh.cells()[indices[k]];consumed[indices[k]]=1U; }
        parents.push_back(merge_siblings(siblings));
    }
    std::vector<AdaptiveHexCell> out;
    out.reserve(mesh.cell_count());
    for (std::size_t i=0;i<mesh.cell_count();++i) if (!consumed[i]) out.push_back(mesh.cells()[i]);
    out.insert(out.end(),parents.begin(),parents.end());
    return AdaptiveHexMesh(std::move(out),mesh.domain_lower(),mesh.domain_upper());
}

std::vector<unsigned char> low_error_marks(const AdaptiveHexMesh& mesh,std::span<const double> indicator,const AdaptiveAmrConfig& cfg) {
    std::vector<unsigned char> mark(mesh.cell_count(),0U);
    if (indicator.empty()) return mark;
    const double peak=*std::max_element(indicator.begin(),indicator.end());
    const double threshold=cfg.coarsen_relative_threshold*peak;
    for (std::size_t i=0;i<mesh.cell_count();++i) {
        if (mesh.cells()[i].level>cfg.minimum_level && indicator[i]<=threshold) mark[i]=1U;
    }
    return mark;
}

} // namespace

AdaptiveHexMesh::AdaptiveHexMesh(std::vector<AdaptiveHexCell> cells,Vec3 domain_lower,Vec3 domain_upper)
    :cells_(std::move(cells)),domain_lower_(domain_lower),domain_upper_(domain_upper){validate();}

AdaptiveHexMesh AdaptiveHexMesh::cartesian(std::size_t nx,std::size_t ny,std::size_t nz,double lx,double ly,double lz) {
    if(nx==0U||ny==0U||nz==0U||!(lx>0.0)||!(ly>0.0)||!(lz>0.0))throw std::invalid_argument("adaptive Cartesian mesh requires positive dimensions and lengths");
    const double dx=lx/static_cast<double>(nx),dy=ly/static_cast<double>(ny),dz=lz/static_cast<double>(nz);
    std::vector<AdaptiveHexCell> cells;cells.reserve(nx*ny*nz);std::size_t root=0U;
    for(std::size_t k=0;k<nz;++k)for(std::size_t j=0;j<ny;++j)for(std::size_t i=0;i<nx;++i,++root){
        const double x0=static_cast<double>(i)*dx;
        const double y0=static_cast<double>(j)*dy;
        const double z0=static_cast<double>(k)*dz;
        cells.push_back({{x0,y0,z0},{x0+dx,y0+dy,z0+dz},0U,root,0U});
    }
    return AdaptiveHexMesh(std::move(cells),{0,0,0},{lx,ly,lz});
}

void AdaptiveHexMesh::validate() const {
    if(cells_.empty())throw std::invalid_argument("AdaptiveHexMesh requires at least one cell");
    if(!finite(domain_lower_)||!finite(domain_upper_)||!(domain_upper_.x>domain_lower_.x)||!(domain_upper_.y>domain_lower_.y)||!(domain_upper_.z>domain_lower_.z))throw std::invalid_argument("invalid adaptive mesh domain");
    double volume=0.0;
    for(const auto& c:cells_){
        if(!finite(c.lower)||!finite(c.upper)||!(c.upper.x>c.lower.x)||!(c.upper.y>c.lower.y)||!(c.upper.z>c.lower.z))throw std::invalid_argument("invalid adaptive cell bounds");
        if(c.lower.x<domain_lower_.x-eps||c.lower.y<domain_lower_.y-eps||c.lower.z<domain_lower_.z-eps||c.upper.x>domain_upper_.x+eps||c.upper.y>domain_upper_.y+eps||c.upper.z>domain_upper_.z+eps)throw std::invalid_argument("adaptive cell lies outside domain");
        if(c.level>20U)throw std::invalid_argument("adaptive refinement level exceeds lineage capacity");
        volume+=cell_volume(c);
    }
    const double domain=(domain_upper_.x-domain_lower_.x)*(domain_upper_.y-domain_lower_.y)*(domain_upper_.z-domain_lower_.z);
    if(std::abs(volume-domain)>1.0e-9*std::max(1.0,domain))throw std::invalid_argument("adaptive cells do not conserve domain volume");
    for(std::size_t i=0;i<cells_.size();++i)for(std::size_t j=i+1;j<cells_.size();++j){
        if(overlap_volume(cells_[i],cells_[j])>1.0e-13*domain)throw std::invalid_argument("adaptive cells overlap");
    }
}

double AdaptiveHexMesh::total_volume() const noexcept {double v=0.0;for(const auto& c:cells_)v+=cell_volume(c);return v;}

PolyMesh AdaptiveHexMesh::poly_mesh() const {
    std::vector<Cell> cells;cells.reserve(cells_.size());
    for(const auto& c:cells_)cells.push_back({cell_center(c),cell_volume(c)});
    std::vector<BoundaryPatch> patches{{"left"},{"right"},{"bottom"},{"top"},{"front"},{"back"}};
    std::vector<Face> faces;
    for(std::size_t i=0;i<cells_.size();++i)for(std::size_t j=i+1;j<cells_.size();++j){
        const auto contact=face_contact(cells_[i],cells_[j]);if(contact.axis<0)continue;
        Vec3 center{},area{};const double a=(contact.hi0-contact.lo0)*(contact.hi1-contact.lo1);
        if(contact.axis==0){center={contact.plane,(contact.lo0+contact.hi0)*0.5,(contact.lo1+contact.hi1)*0.5};area={contact.direction*a,0,0};}
        else if(contact.axis==1){center={(contact.lo0+contact.hi0)*0.5,contact.plane,(contact.lo1+contact.hi1)*0.5};area={0,contact.direction*a,0};}
        else {center={(contact.lo0+contact.hi0)*0.5,(contact.lo1+contact.hi1)*0.5,contact.plane};area={0,0,contact.direction*a};}
        faces.push_back({i,j,center,area,invalid_patch});
    }
    for(std::size_t i=0;i<cells_.size();++i){const auto& c=cells_[i];const double dx=c.upper.x-c.lower.x,dy=c.upper.y-c.lower.y,dz=c.upper.z-c.lower.z;
        if(near(c.lower.x,domain_lower_.x))faces.push_back({i,invalid_cell,{c.lower.x,(c.lower.y+c.upper.y)*.5,(c.lower.z+c.upper.z)*.5},{-dy*dz,0,0},0});
        if(near(c.upper.x,domain_upper_.x))faces.push_back({i,invalid_cell,{c.upper.x,(c.lower.y+c.upper.y)*.5,(c.lower.z+c.upper.z)*.5},{dy*dz,0,0},1});
        if(near(c.lower.y,domain_lower_.y))faces.push_back({i,invalid_cell,{(c.lower.x+c.upper.x)*.5,c.lower.y,(c.lower.z+c.upper.z)*.5},{0,-dx*dz,0},2});
        if(near(c.upper.y,domain_upper_.y))faces.push_back({i,invalid_cell,{(c.lower.x+c.upper.x)*.5,c.upper.y,(c.lower.z+c.upper.z)*.5},{0,dx*dz,0},3});
        if(near(c.lower.z,domain_lower_.z))faces.push_back({i,invalid_cell,{(c.lower.x+c.upper.x)*.5,(c.lower.y+c.upper.y)*.5,c.lower.z},{0,0,-dx*dy},4});
        if(near(c.upper.z,domain_upper_.z))faces.push_back({i,invalid_cell,{(c.lower.x+c.upper.x)*.5,(c.lower.y+c.upper.y)*.5,c.upper.z},{0,0,dx*dy},5});
    }
    return PolyMesh(std::move(cells),std::move(faces),std::move(patches));
}

TopologyChangeMap build_topology_change_map(const AdaptiveHexMesh& old_mesh,const AdaptiveHexMesh& new_mesh){
    TopologyChangeMap map;map.old_cell_count=old_mesh.cell_count();map.new_cell_count=new_mesh.cell_count();map.new_from_old.resize(new_mesh.cell_count());
    for(std::size_t j=0;j<new_mesh.cell_count();++j){double covered=0.0;for(std::size_t i=0;i<old_mesh.cell_count();++i){const double v=overlap_volume(old_mesh.cells()[i],new_mesh.cells()[j]);if(v>eps){map.new_from_old[j].push_back({i,v});covered+=v;}}
        const double target=cell_volume(new_mesh.cells()[j]);if(std::abs(covered-target)>1.0e-9*std::max(1.0,target))throw std::runtime_error("topology change map does not fully cover target cell");}
    return map;
}

std::vector<double> conservative_adaptive_remap(const AdaptiveHexMesh& old_mesh,
                                                      std::span<const double> old_field,
                                                      const AdaptiveHexMesh& new_mesh) {
    if (old_field.size()!=old_mesh.cell_count()) {
        throw std::invalid_argument("adaptive remap field size mismatch");
    }
    for (double value:old_field) {
        if (!std::isfinite(value)) throw std::invalid_argument("adaptive remap requires finite field values");
    }
    const auto map=build_topology_change_map(old_mesh,new_mesh);
    std::vector<double> out(new_mesh.cell_count(),0.0);
    for (std::size_t j=0;j<out.size();++j) {
        const double volume=cell_volume(new_mesh.cells()[j]);
        for (const auto& weight:map.new_from_old[j]) {
            out[j]+=old_field[weight.source_cell]*weight.volume/volume;
        }
    }
    return out;
}

MarkedHexRefinement::MarkedHexRefinement(std::vector<unsigned char> marked,
                                         std::size_t maximum_level,
                                         bool enforce_two_to_one)
    :marked_(std::move(marked)),maximum_level_(maximum_level),enforce_two_to_one_(enforce_two_to_one) {
    if (maximum_level_>20U) throw std::invalid_argument("maximum adaptive level exceeds lineage capacity");
}

AdaptiveTopologyResult MarkedHexRefinement::apply(const AdaptiveHexMesh& mesh) const {
    auto out=refine_once(mesh,marked_,maximum_level_);
    if (enforce_two_to_one_) out=enforce_balance(std::move(out),maximum_level_);
    return {out,build_topology_change_map(mesh,out)};
}

MarkedHexCoarsening::MarkedHexCoarsening(std::vector<unsigned char> marked,
                                         std::size_t minimum_level,
                                         bool enforce_two_to_one)
    :marked_(std::move(marked)),minimum_level_(minimum_level),enforce_two_to_one_(enforce_two_to_one) {}

AdaptiveTopologyResult MarkedHexCoarsening::apply(const AdaptiveHexMesh& mesh) const {
    auto out=coarsen_once(mesh,marked_,minimum_level_);
    if (enforce_two_to_one_) {
        std::size_t max_level=0U;
        for (const auto& cell:mesh.cells()) max_level=std::max(max_level,cell.level);
        out=enforce_balance(std::move(out),max_level);
    }
    return {out,build_topology_change_map(mesh,out)};
}

bool satisfies_two_to_one(const AdaptiveHexMesh& mesh) {
    for (std::size_t i=0;i<mesh.cell_count();++i) {
        for (std::size_t j=i+1;j<mesh.cell_count();++j) {
            if (face_contact(mesh.cells()[i],mesh.cells()[j]).axis<0) continue;
            const auto a=mesh.cells()[i].level;
            const auto b=mesh.cells()[j].level;
            if (a>b+1U || b>a+1U) return false;
        }
    }
    return true;
}

std::vector<double> scalar_jump_error_indicator(const AdaptiveHexMesh& mesh,
                                                std::span<const double> field) {
    if (field.size()!=mesh.cell_count()) throw std::invalid_argument("AMR indicator field size mismatch");
    for (double value:field) {
        if (!std::isfinite(value)) throw std::invalid_argument("AMR indicator requires finite field values");
    }
    std::vector<double> eta(mesh.cell_count(),0.0);
    const auto poly=mesh.poly_mesh();
    for (const auto& face:poly.faces()) {
        if (face.boundary()) continue;
        const auto owner=face.owner;
        const auto neighbour=face.neighbour;
        const double distance=magnitude(poly.cells()[neighbour].center-poly.cells()[owner].center);
        if (!(distance>0.0)) continue;
        const double jump=field[neighbour]-field[owner];
        const double contribution=magnitude(face.area)*jump*jump/distance;
        eta[owner]+=contribution;
        eta[neighbour]+=contribution;
    }
    for (double& value:eta) value=std::sqrt(std::max(0.0,value));
    return eta;
}

std::vector<unsigned char> mark_dorfler(std::span<const double> indicators,double theta) {
    if (!(theta>0.0&&theta<=1.0)||!std::isfinite(theta)) {
        throw std::invalid_argument("Dorfler theta must lie in (0,1]");
    }
    std::vector<unsigned char> mark(indicators.size(),0U);
    std::vector<std::size_t> order(indicators.size());
    std::iota(order.begin(),order.end(),0U);
    double total=0.0;
    for (double indicator:indicators) {
        if (!(indicator>=0.0)||!std::isfinite(indicator)) {
            throw std::invalid_argument("Dorfler indicators must be finite non-negative");
        }
        total+=indicator*indicator;
    }
    if (total<=0.0) return mark;
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b){return indicators[a]>indicators[b];});
    double selected=0.0;
    for (const auto i:order) {
        mark[i]=1U;
        selected+=indicators[i]*indicators[i];
        if (selected>=theta*total) break;
    }
    return mark;
}

AdaptiveAmrResult adapt_scalar_field(const AdaptiveHexMesh& mesh,
                                     std::span<const double> field,
                                     const AdaptiveAmrConfig& cfg) {
    if (field.size()!=mesh.cell_count()) throw std::invalid_argument("AMR field size mismatch");
    if (cfg.minimum_level>cfg.maximum_level||cfg.maximum_level>20U) {
        throw std::invalid_argument("invalid AMR level bounds");
    }
    if (!(cfg.refine_theta>0.0&&cfg.refine_theta<=1.0)
        || !(cfg.coarsen_relative_threshold>=0.0&&cfg.coarsen_relative_threshold<=1.0)) {
        throw std::invalid_argument("invalid AMR marking controls");
    }

    auto current=mesh;
    std::vector<double> values(field.begin(),field.end());
    auto indicator=scalar_jump_error_indicator(current,values);
    const auto before_coarsen=current.cell_count();

    auto coarsen_mark=low_error_marks(current,indicator,cfg);
    MarkedHexCoarsening coarsen(std::move(coarsen_mark),cfg.minimum_level,cfg.enforce_two_to_one);
    auto coarsened=coarsen.apply(current);
    if (coarsened.mesh.cell_count()!=current.cell_count()) {
        values=conservative_adaptive_remap(current,values,coarsened.mesh);
    }
    current=std::move(coarsened.mesh);
    const auto after_coarsen=current.cell_count();

    indicator=scalar_jump_error_indicator(current,values);
    auto refine_mark=mark_dorfler(indicator,cfg.refine_theta);
    for (std::size_t i=0;i<refine_mark.size();++i) {
        if (current.cells()[i].level>=cfg.maximum_level) refine_mark[i]=0U;
    }
    const auto before_refine=current.cell_count();
    MarkedHexRefinement refine(std::move(refine_mark),cfg.maximum_level,cfg.enforce_two_to_one);
    auto refined=refine.apply(current);
    if (refined.mesh.cell_count()!=current.cell_count()) {
        values=conservative_adaptive_remap(current,values,refined.mesh);
    }
    current=std::move(refined.mesh);
    indicator=scalar_jump_error_indicator(current,values);

    AdaptiveAmrResult result{std::move(current),std::move(values),std::move(indicator),0U,0U};
    result.coarsened_cell_count=before_coarsen>after_coarsen?before_coarsen-after_coarsen:0U;
    result.refined_cell_count=result.mesh.cell_count()>before_refine?result.mesh.cell_count()-before_refine:0U;
    return result;
}

} // namespace cfd::fvm
