#include "cfd/fvm/adaptive_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}

double integral(const cfd::fvm::AdaptiveHexMesh& mesh,const std::vector<double>& field){double sum=0.0;for(std::size_t i=0;i<mesh.cell_count();++i){const auto& c=mesh.cells()[i];sum+=field[i]*(c.upper.x-c.lower.x)*(c.upper.y-c.lower.y)*(c.upper.z-c.lower.z);}return sum;}

void test_local_refinement_and_nonconformal_faces(){
    const auto base=cfd::fvm::AdaptiveHexMesh::cartesian(2,1,1);
    std::vector<unsigned char> mark(base.cell_count(),0U);mark[0]=1U;
    cfd::fvm::MarkedHexRefinement op(mark,3U,true);
    const cfd::fvm::MeshTopologyOperation& iface=op;
    const auto result=iface.apply(base);
    require(result.mesh.cell_count()==9U,"one hex split should replace one cell by eight children");
    require(std::abs(result.mesh.total_volume()-1.0)<1e-12,"refinement must conserve volume");
    const auto poly=result.mesh.poly_mesh();
    std::size_t coarse_fine_faces=0U;
    for(const auto& face:poly.faces())if(!face.boundary()){
        const auto lo=std::min(result.mesh.cells()[face.owner].level,result.mesh.cells()[face.neighbour].level);
        const auto hi=std::max(result.mesh.cells()[face.owner].level,result.mesh.cells()[face.neighbour].level);
        if(lo==0U&&hi==1U)++coarse_fine_faces;
    }
    require(coarse_fine_faces==4U,"coarse/fine hex interface should be represented by four conservative subfaces");
    require(result.mapping.new_cell_count==9U&&result.mapping.old_cell_count==2U,"topology map should describe old/new cell counts");
}

void test_refine_coarsen_roundtrip_and_conservative_remap(){
    const auto base=cfd::fvm::AdaptiveHexMesh::cartesian(2,1,1);
    std::vector<double> original{2.0,5.0};
    std::vector<unsigned char> mark(base.cell_count(),1U);
    auto refined=cfd::fvm::MarkedHexRefinement(mark,2U,true).apply(base).mesh;
    auto fine_field=cfd::fvm::conservative_adaptive_remap(base,original,refined);
    require(std::abs(integral(base,original)-integral(refined,fine_field))<1e-12,"refinement remap must conserve scalar integral exactly");
    std::vector<unsigned char> coarsen_mark(refined.cell_count(),1U);
    auto coarsened=cfd::fvm::MarkedHexCoarsening(coarsen_mark,0U,true).apply(refined).mesh;
    auto coarse_field=cfd::fvm::conservative_adaptive_remap(refined,fine_field,coarsened);
    require(coarsened.cell_count()==base.cell_count(),"coarsening all complete sibling groups should restore root mesh");
    require(std::abs(coarse_field[0]-2.0)<1e-12&&std::abs(coarse_field[1]-5.0)<1e-12,"refine/coarsen remap should preserve piecewise-constant root fields");
}

void test_two_to_one_balancing(){
    auto mesh=cfd::fvm::AdaptiveHexMesh::cartesian(1,1,1);
    mesh=cfd::fvm::MarkedHexRefinement(std::vector<unsigned char>(1U,1U),4U,true).apply(mesh).mesh;
    std::vector<unsigned char> one(mesh.cell_count(),0U);one[0]=1U;mesh=cfd::fvm::MarkedHexRefinement(one,4U,true).apply(mesh).mesh;
    one.assign(mesh.cell_count(),0U);
    auto finest=std::max_element(mesh.cells().begin(),mesh.cells().end(),[](const auto& a,const auto& b){return a.level<b.level;});
    one[static_cast<std::size_t>(std::distance(mesh.cells().begin(),finest))]=1U;
    mesh=cfd::fvm::MarkedHexRefinement(one,4U,true).apply(mesh).mesh;
    require(cfd::fvm::satisfies_two_to_one(mesh),"balanced refinement must keep face-neighbour levels within one");
}

void test_error_indicator_and_dorfler_amr(){
    auto mesh=cfd::fvm::AdaptiveHexMesh::cartesian(8,1,1);
    std::vector<double> field(mesh.cell_count());
    for(std::size_t i=0;i<mesh.cell_count();++i)field[i]=mesh.cells()[i].upper.x<=0.5?0.0:1.0;
    const auto indicator=cfd::fvm::scalar_jump_error_indicator(mesh,field);
    require(*std::max_element(indicator.begin(),indicator.end())>0.0,"jump estimator should detect a discontinuity");
    const auto marked=cfd::fvm::mark_dorfler(indicator,0.5);
    require(std::count(marked.begin(),marked.end(),1U)>0,"Dorfler marking must select nonzero-error cells");
    cfd::fvm::AdaptiveAmrConfig cfg;cfg.maximum_level=2U;cfg.refine_theta=0.6;cfg.coarsen_relative_threshold=0.0;
    const auto adapted=cfd::fvm::adapt_scalar_field(mesh,field,cfg);
    require(adapted.mesh.cell_count()>mesh.cell_count(),"error-driven AMR should refine the jump region");
    bool refined_near_interface=false;
    for(const auto& c:adapted.mesh.cells())if(c.level>0U&&c.lower.x<0.5+1e-12&&c.upper.x>0.5-1e-12)refined_near_interface=true;
    require(refined_near_interface,"AMR should place refined cells at the scalar discontinuity");
    require(std::abs(integral(mesh,field)-integral(adapted.mesh,adapted.field))<1e-12,"AMR field transfer must conserve the scalar integral");
}

void test_low_error_coarsening(){
    auto mesh=cfd::fvm::AdaptiveHexMesh::cartesian(1,1,1);
    mesh=cfd::fvm::MarkedHexRefinement(std::vector<unsigned char>(1U,1U),2U,true).apply(mesh).mesh;
    std::vector<double> field(mesh.cell_count(),3.0);
    cfd::fvm::AdaptiveAmrConfig cfg;cfg.maximum_level=2U;cfg.minimum_level=0U;cfg.refine_theta=0.5;cfg.coarsen_relative_threshold=1.0;
    const auto adapted=cfd::fvm::adapt_scalar_field(mesh,field,cfg);
    require(adapted.mesh.cell_count()==1U,"zero-error sibling group should coarsen back to its parent");
    require(std::abs(adapted.field.front()-3.0)<1e-12,"coarsening must preserve uniform field");
}
}

int main(){try{test_local_refinement_and_nonconformal_faces();test_refine_coarsen_roundtrip_and_conservative_remap();test_two_to_one_balancing();test_error_indicator_and_dorfler_amr();test_low_error_coarsening();std::cout<<"v0.15.3 FVM AMR tests passed\n";return 0;}catch(const std::exception& e){std::cerr<<"v0.15.3 FVM AMR test failure: "<<e.what()<<'\n';return 1;}}
