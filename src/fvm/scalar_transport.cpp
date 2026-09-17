#include "cfd/solvers/fvm/scalar_transport.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/fvm/pressure_velocity.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::fvm {
namespace { double outward_flux(const Face& face,std::size_t cell,double phi) noexcept { return face.owner==cell?phi:-phi; } }
ScalarTransport::ScalarTransport(PolyMesh mesh,ScalarTransportConfig config)
    :mesh_(std::move(mesh)),config_(config),boundary_(mesh_.patches().size()),values_(mesh_.cell_count(),0.0),
     old_values_(mesh_.cell_count(),0.0),source_(mesh_.cell_count(),0.0),face_flux_(mesh_.face_count(),0.0),
     boundary_rhs_(mesh_.cell_count(),0.0),rhs_(mesh_.cell_count(),0.0),candidate_(mesh_.cell_count(),0.0),spatial_action_(mesh_.cell_count(),0.0) {
    if(!(config_.dt>0.0)||!(config_.diffusivity>=0.0)||config_.linear_iterations==0U||config_.gmres_restart==0U||!(config_.linear_tolerance>0.0))
        throw std::invalid_argument("invalid scalar transport controls");
    if(config_.temporal_scheme==TemporalScheme::crank_nicolson)
        (void)crank_nicolson_implicit_weight(config_.crank_nicolson_off_centering);
}
void ScalarTransport::set_boundary(std::string_view patch,ScalarBoundaryType type,double value){boundary_[mesh_.patch_index(patch)]={type,value};operator_dirty_=true;}
void ScalarTransport::initialize(double value){std::fill(values_.begin(),values_.end(),value);old_values_=values_;time_=0.0;steps_=0U;}
void ScalarTransport::initialize(const std::function<double(Vec3)>& value){for(std::size_t c=0;c<mesh_.cell_count();++c)values_[c]=value(mesh_.cells()[c].center);old_values_=values_;time_=0.0;steps_=0U;}
void ScalarTransport::set_source(double value){std::fill(source_.begin(),source_.end(),value);}
void ScalarTransport::set_source(const std::function<double(Vec3)>& value){for(std::size_t c=0;c<mesh_.cell_count();++c)source_[c]=value(mesh_.cells()[c].center);}
void ScalarTransport::set_face_flux(std::vector<double> flux){if(flux.size()!=mesh_.face_count())throw std::invalid_argument("scalar transport face-flux size mismatch");face_flux_=std::move(flux);operator_dirty_=true;}
void ScalarTransport::assemble_operator(){
    const std::size_t n=mesh_.cell_count(); cfd::core::CsrBuilder builder(n,n);
    for(std::size_t cell=0;cell<n;++cell){double ap=0.0,b=0.0;
        for(const std::size_t fi:mesh_.cell_faces()[cell]){const auto& face=mesh_.faces()[fi];const double phi=outward_flux(face,cell,face_flux_[fi]);
            if(!face.boundary()){const std::size_t other=face.owner==cell?face.neighbour:face.owner;const double diff=config_.diffusivity*decompose_face_area(mesh_,fi).orthogonal_metric;
                ap+=diff+std::max(phi,0.0);const double off=-diff+std::min(phi,0.0);if(off!=0.0)builder.add(cell,other,off);
            }else{const auto& bc=boundary_[face.patch];if(bc.type==ScalarBoundaryType::fixedValue){const double distance=magnitude(face.center-mesh_.cells()[cell].center);if(!(distance>0.0))throw std::runtime_error("degenerate scalar boundary distance");const double diff=config_.diffusivity*magnitude(face.area)/distance;ap+=diff+std::max(phi,0.0);b+=diff*bc.value-std::min(phi,0.0)*bc.value;}else ap+=phi;}}
        if(!std::isfinite(ap)) throw std::runtime_error("non-finite scalar transport spatial diagonal");
        if(ap!=0.0)builder.add(cell,cell,ap);
        boundary_rhs_[cell]=b;
    }
    spatial_matrix_=builder.build();operator_dirty_=false;++operator_assemblies_;
}
cfd::core::IterativeSolverResult ScalarTransport::step(){
    if(operator_dirty_)assemble_operator();
    const std::size_t n=mesh_.cell_count();
    double alpha0=1.0,theta=1.0;
    const bool use_bdf2=config_.temporal_scheme==TemporalScheme::backward_bdf2&&steps_>0U;
    if(use_bdf2)alpha0=1.5;
    if(config_.temporal_scheme==TemporalScheme::crank_nicolson)theta=crank_nicolson_implicit_weight(config_.crank_nicolson_off_centering);

    cfd::core::CsrBuilder builder(n,n);
    const auto& ro=spatial_matrix_.row_offsets();const auto& ci=spatial_matrix_.column_indices();const auto& av=spatial_matrix_.values();
    for(std::size_t row=0;row<n;++row){
        for(std::size_t k=ro[row];k<ro[row+1U];++k)builder.add(row,ci[k],theta*av[k]);
        builder.add(row,row,alpha0*mesh_.cells()[row].volume/config_.dt);
    }
    matrix_=builder.build();preconditioner_.emplace(matrix_);
    if(theta<1.0)spatial_matrix_.multiply(values_,spatial_action_);else std::fill(spatial_action_.begin(),spatial_action_.end(),0.0);
    candidate_=values_;
    for(std::size_t c=0;c<n;++c){
        const double vdt=mesh_.cells()[c].volume/config_.dt;
        const double temporal=use_bdf2?(2.0*values_[c]-0.5*old_values_[c])*vdt:values_[c]*vdt;
        rhs_[c]=temporal+source_[c]*mesh_.cells()[c].volume+boundary_rhs_[c]-(1.0-theta)*spatial_action_[c];
    }
    linear_result_=cfd::core::restarted_gmres(rhs_,candidate_,[&](std::span<const double> in,std::span<double> out){matrix_.multiply(in,out);},[&](std::span<const double> in,std::span<double> out){(*preconditioner_)(in,out);},config_.linear_iterations,config_.gmres_restart,config_.linear_tolerance);
    if(!linear_result_.converged) throw std::runtime_error("scalar transport linear solve did not converge");
    old_values_=values_;values_.swap(candidate_);++steps_;time_+=config_.dt;
    return linear_result_;
}
cfd::core::IterativeSolverResult ScalarTransport::run(std::size_t steps){if(steps==0U)return linear_result_;for(std::size_t i=0;i<steps;++i)step();return linear_result_;}
double ScalarTransport::volume_integral() const{return cfd::core::parallel_sum(mesh_.cell_count(),[&](std::size_t c){return values_[c]*mesh_.cells()[c].volume;});}
double ScalarTransport::minimum() const{return *std::min_element(values_.begin(),values_.end());}
double ScalarTransport::maximum() const{return *std::max_element(values_.begin(),values_.end());}
} // namespace cfd::fvm
