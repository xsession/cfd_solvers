#include "cfd/multibody/advanced_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cfd::multibody {
namespace {

[[nodiscard]] std::vector<double> solve_dense(std::vector<double> a,std::vector<double> b,std::size_t n){
    if(a.size()!=n*n||b.size()!=n) throw std::invalid_argument("dense solve size mismatch");
    for(std::size_t col=0;col<n;++col){
        std::size_t pivot=col;
        for(std::size_t row=col+1;row<n;++row) if(std::abs(a[row*n+col])>std::abs(a[pivot*n+col])) pivot=row;
        if(!(std::abs(a[pivot*n+col])>1.0e-14)) throw std::runtime_error("singular generalized-coordinate matrix");
        if(pivot!=col){
            for(std::size_t j=col;j<n;++j) std::swap(a[pivot*n+j],a[col*n+j]);
            std::swap(b[pivot],b[col]);
        }
        const double diagonal=a[col*n+col];
        for(std::size_t j=col;j<n;++j) a[col*n+j]/=diagonal;
        b[col]/=diagonal;
        for(std::size_t row=0;row<n;++row){
            if(row==col) continue;
            const double f=a[row*n+col];
            if(f==0.0) continue;
            for(std::size_t j=col;j<n;++j) a[row*n+j]-=f*a[col*n+j];
            b[row]-=f*b[col];
        }
    }
    return b;
}

[[nodiscard]] double inf_norm(std::span<const double> x) noexcept {
    double value=0.0;
    for(double v:x) value=std::max(value,std::abs(v));
    return value;
}

[[nodiscard]] Vec3 world_inertia_mul(const RigidBody& body,Vec3 world_vector){
    const Vec3 local=inverse_rotate(body.state.orientation,world_vector);
    return rotate(body.state.orientation,component_mul(body.inertia_diagonal,local));
}

} // namespace

ImplicitSecondOrderConfig ImplicitSecondOrderConfig::newmark(double dt,double beta,double gamma){
    ImplicitSecondOrderConfig c;c.dt=dt;c.beta=beta;c.gamma=gamma;c.alpha_m=0.0;c.alpha_f=0.0;return c;
}

ImplicitSecondOrderConfig ImplicitSecondOrderConfig::generalized_alpha(double dt,double rho){
    if(!(rho>=0.0&&rho<=1.0)||!std::isfinite(rho)) throw std::invalid_argument("generalized-alpha spectral radius must be in [0,1]");
    ImplicitSecondOrderConfig c;c.dt=dt;
    c.alpha_m=(2.0*rho-1.0)/(rho+1.0);
    c.alpha_f=rho/(rho+1.0);
    c.gamma=0.5-c.alpha_m+c.alpha_f;
    const double term=1.0-c.alpha_m+c.alpha_f;
    c.beta=0.25*term*term;
    return c;
}

ImplicitSecondOrderConfig ImplicitSecondOrderConfig::hht(double dt,double alpha){
    if(!(alpha>=-1.0/3.0&&alpha<=0.0)||!std::isfinite(alpha)) throw std::invalid_argument("HHT alpha must be in [-1/3,0]");
    ImplicitSecondOrderConfig c;c.dt=dt;c.alpha_m=0.0;c.alpha_f=-alpha;
    c.gamma=0.5-alpha;
    c.beta=0.25*(1.0-alpha)*(1.0-alpha);
    return c;
}

ImplicitGeneralizedIntegrator::ImplicitGeneralizedIntegrator(std::vector<double> mass_matrix,std::size_t dofs,
                                                             GeneralizedForceFunction force,ImplicitSecondOrderConfig config)
    :mass_(std::move(mass_matrix)),n_(dofs),force_(std::move(force)),config_(config),q_(dofs),v_(dofs),a_(dofs){
    if(n_==0||mass_.size()!=n_*n_) throw std::invalid_argument("implicit multibody mass matrix size mismatch");
    if(!force_) throw std::invalid_argument("implicit multibody force callback required");
    if(!(config_.dt>0.0&&config_.beta>0.0&&config_.gamma>0.0&&config_.newton_iterations>0&&config_.nonlinear_tolerance>0.0&&config_.finite_difference_epsilon>0.0))
        throw std::invalid_argument("invalid implicit multibody controls");
    if(!(config_.alpha_m<1.0&&config_.alpha_f<1.0)) throw std::invalid_argument("generalized-alpha interpolation must remain implicit");
}

void ImplicitGeneralizedIntegrator::initialize(std::vector<double> q,std::vector<double> v,double time){
    if(q.size()!=n_) throw std::invalid_argument("implicit multibody q size mismatch");
    if(v.empty()) v.assign(n_,0.0);
    if(v.size()!=n_) throw std::invalid_argument("implicit multibody v size mismatch");
    if(!std::isfinite(time)) throw std::invalid_argument("implicit multibody time must be finite");
    q_=std::move(q);v_=std::move(v);time_=time;
    std::vector<double> f(n_);force_(q_,v_,time_,f);
    a_=solve_dense(mass_,std::move(f),n_);
    last_newton_iterations_=0;last_residual_norm_=0.0;
}

std::vector<double> ImplicitGeneralizedIntegrator::residual(std::span<const double> a1,
                                                             std::span<const double> q0,
                                                             std::span<const double> v0,
                                                             std::span<const double> a0) const {
    const double dt=config_.dt,b=config_.beta,g=config_.gamma;
    std::vector<double> q1(n_),v1(n_),qaf(n_),vaf(n_),aam(n_),f(n_),r(n_);
    for(std::size_t i=0;i<n_;++i){
        const double qp=q0[i]+dt*v0[i]+dt*dt*(0.5-b)*a0[i];
        const double vp=v0[i]+dt*(1.0-g)*a0[i];
        q1[i]=qp+b*dt*dt*a1[i];
        v1[i]=vp+g*dt*a1[i];
        qaf[i]=(1.0-config_.alpha_f)*q1[i]+config_.alpha_f*q0[i];
        vaf[i]=(1.0-config_.alpha_f)*v1[i]+config_.alpha_f*v0[i];
        aam[i]=(1.0-config_.alpha_m)*a1[i]+config_.alpha_m*a0[i];
    }
    const double teval=(1.0-config_.alpha_f)*(time_+dt)+config_.alpha_f*time_;
    force_(qaf,vaf,teval,f);
    for(std::size_t i=0;i<n_;++i){
        double ma=0.0;
        for(std::size_t j=0;j<n_;++j) ma+=mass_[i*n_+j]*aam[j];
        r[i]=ma-f[i];
    }
    return r;
}

void ImplicitGeneralizedIntegrator::step(){
    const auto q0=q_,v0=v_,a0=a_;
    std::vector<double> a1=a0;
    last_newton_iterations_=0;
    for(std::size_t iteration=0;iteration<config_.newton_iterations;++iteration){
        auto r=residual(a1,q0,v0,a0);
        last_residual_norm_=inf_norm(r);
        if(last_residual_norm_<=config_.nonlinear_tolerance){ last_newton_iterations_=iteration;break; }
        std::vector<double> jac(n_*n_);
        for(std::size_t col=0;col<n_;++col){
            auto perturbed=a1;
            const double h=config_.finite_difference_epsilon*std::max(1.0,std::abs(a1[col]));
            perturbed[col]+=h;
            const auto rp=residual(perturbed,q0,v0,a0);
            for(std::size_t row=0;row<n_;++row) jac[row*n_+col]=(rp[row]-r[row])/h;
        }
        for(double& value:r) value=-value;
        const auto delta=solve_dense(std::move(jac),std::move(r),n_);
        for(std::size_t i=0;i<n_;++i) a1[i]+=delta[i];
        last_newton_iterations_=iteration+1;
        if(inf_norm(delta)<=config_.nonlinear_tolerance*std::max(1.0,inf_norm(a1))){
            const auto rr=residual(a1,q0,v0,a0);last_residual_norm_=inf_norm(rr);break;
        }
    }
    if(last_residual_norm_>100.0*config_.nonlinear_tolerance) throw std::runtime_error("implicit multibody Newton solve did not converge");
    const double dt=config_.dt,b=config_.beta,g=config_.gamma;
    for(std::size_t i=0;i<n_;++i){
        const double qp=q0[i]+dt*v0[i]+dt*dt*(0.5-b)*a0[i];
        const double vp=v0[i]+dt*(1.0-g)*a0[i];
        q_[i]=qp+b*dt*dt*a1[i];
        v_[i]=vp+g*dt*a1[i];
    }
    a_=std::move(a1);time_+=dt;
}
void ImplicitGeneralizedIntegrator::run(std::size_t steps){for(std::size_t i=0;i<steps;++i)step();}

ArticulatedSystem::ArticulatedSystem(Vec3 gravity):gravity_(gravity){}

std::size_t ArticulatedSystem::add_link(ArticulatedLink link){
    if(link.parent!=invalid_body&&link.parent>=links_.size()) throw std::invalid_argument("articulated parent must precede child");
    if(norm(link.axis_parent)<1.0e-14) throw std::invalid_argument("articulated joint axis must be nonzero");
    link.axis_parent=normalized(link.axis_parent);link.body.validate();
    const auto id=links_.size();links_.push_back(std::move(link));
    joint_origin_world_.resize(links_.size());joint_axis_world_.resize(links_.size());forward_kinematics();return id;
}
void ArticulatedSystem::set_effort(std::size_t link,double effort){links_.at(link).applied_effort=effort;}

bool ArticulatedSystem::is_ancestor(std::size_t joint,std::size_t body) const {
    std::size_t current=body;
    while(current!=invalid_body){if(current==joint)return true;current=links_[current].parent;}
    return false;
}

void ArticulatedSystem::forward_kinematics(){
    for(std::size_t i=0;i<links_.size();++i){
        auto& link=links_[i];
        Vec3 parent_position{};Quaternion parent_orientation{};
        if(link.parent!=invalid_body){parent_position=links_[link.parent].body.state.position;parent_orientation=links_[link.parent].body.state.orientation;}
        const Vec3 axis_world=rotate(parent_orientation,link.axis_parent);
        Vec3 joint_origin=parent_position+rotate(parent_orientation,link.parent_joint_offset);
        Quaternion child_orientation{};
        if(link.joint_type==ReducedJointType::revolute){
            child_orientation=normalized(parent_orientation*quaternion_from_axis_angle(link.axis_parent,link.q)*link.reference_orientation);
        }else{
            joint_origin+=axis_world*link.q;
            child_orientation=normalized(parent_orientation*link.reference_orientation);
        }
        joint_origin_world_[i]=joint_origin;joint_axis_world_[i]=normalized(axis_world);
        link.body.state.orientation=child_orientation;
        link.body.state.position=joint_origin+rotate(child_orientation,link.joint_to_com_local);
    }
    // Velocities are reconstructed exactly from the generalized-coordinate Jacobians.
    for(std::size_t body=0;body<links_.size();++body){
        Vec3 v{},w{};
        for(std::size_t joint=0;joint<links_.size();++joint){
            if(!is_ancestor(joint,body)) continue;
            v+=linear_jacobian(joint,body)*links_[joint].qd;
            w+=angular_jacobian(joint,body)*links_[joint].qd;
        }
        links_[body].body.state.linear_velocity=v;links_[body].body.state.angular_velocity=w;
    }
}

Vec3 ArticulatedSystem::linear_jacobian(std::size_t joint,std::size_t body) const {
    if(!is_ancestor(joint,body)) return {};
    if(links_[joint].joint_type==ReducedJointType::prismatic) return joint_axis_world_[joint];
    return cross(joint_axis_world_[joint],links_[body].body.state.position-joint_origin_world_[joint]);
}
Vec3 ArticulatedSystem::angular_jacobian(std::size_t joint,std::size_t body) const {
    if(!is_ancestor(joint,body)||links_[joint].joint_type==ReducedJointType::prismatic) return {};
    return joint_axis_world_[joint];
}

std::vector<double> ArticulatedSystem::mass_matrix() const {
    const std::size_t n=links_.size();std::vector<double> matrix(n*n);
    for(std::size_t body=0;body<n;++body){
        const auto& rb=links_[body].body;
        for(std::size_t i=0;i<n;++i){
            if(!is_ancestor(i,body)) continue;
            const Vec3 jvi=linear_jacobian(i,body),jwi=angular_jacobian(i,body);
            for(std::size_t j=0;j<n;++j){
                if(!is_ancestor(j,body)) continue;
                const Vec3 jvj=linear_jacobian(j,body),jwj=angular_jacobian(j,body);
                matrix[i*n+j]+=rb.mass*dot(jvi,jvj)+dot(jwi,world_inertia_mul(rb,jwj));
            }
        }
    }
    return matrix;
}

std::vector<double> ArticulatedSystem::generalized_forces() const {
    const std::size_t n=links_.size();std::vector<double> force(n);
    for(std::size_t i=0;i<n;++i){
        force[i]=links_[i].applied_effort-links_[i].damping*links_[i].qd;
        for(std::size_t body=0;body<n;++body){
            if(!is_ancestor(i,body)) continue;
            force[i]+=links_[body].body.mass*dot(linear_jacobian(i,body),gravity_);
            force[i]+=dot(angular_jacobian(i,body),links_[body].body.torque);
            force[i]+=dot(linear_jacobian(i,body),links_[body].body.force);
        }
    }
    return force;
}

void ArticulatedSystem::step(double dt){
    if(!(dt>0.0)||!std::isfinite(dt)) throw std::invalid_argument("articulated dt must be finite and positive");
    forward_kinematics();const auto m=mass_matrix();auto rhs=generalized_forces();const auto qdd=solve_dense(m,std::move(rhs),links_.size());
    for(std::size_t i=0;i<links_.size();++i){links_[i].qdd=qdd[i];links_[i].qd+=dt*qdd[i];links_[i].q+=dt*links_[i].qd;links_[i].body.clear_accumulators();}
    forward_kinematics();
}
void ArticulatedSystem::run(std::size_t steps,double dt){for(std::size_t i=0;i<steps;++i)step(dt);}
double ArticulatedSystem::kinetic_energy() const noexcept {double e=0.0;for(const auto& l:links_){const auto& b=l.body;e+=0.5*b.mass*norm_squared(b.state.linear_velocity)+0.5*dot(b.state.angular_velocity,world_inertia_mul(b,b.state.angular_velocity));}return e;}

FlexibleBodyInterface::FlexibleBodyInterface(std::vector<FlexibleInterfaceNode> nodes,std::vector<FlexibleMode> modes)
    :nodes_(std::move(nodes)),modes_(std::move(modes)),q_(modes_.size()),qd_(modes_.size()){
    if(nodes_.empty()) throw std::invalid_argument("flexible interface requires nodes");
    for(const auto& mode:modes_){
        if(mode.shape_local.size()!=nodes_.size()) throw std::invalid_argument("flexible mode shape size mismatch");
        if(!(mode.modal_mass>0.0)||!std::isfinite(mode.modal_mass)) throw std::invalid_argument("flexible modal mass must be finite and positive");
    }
}
void FlexibleBodyInterface::set_modal_state(std::vector<double> q,std::vector<double> qd){
    if(q.size()!=modes_.size()) throw std::invalid_argument("flexible modal q size mismatch");
    if(qd.empty()) qd.assign(modes_.size(),0.0);
    if(qd.size()!=modes_.size()) throw std::invalid_argument("flexible modal qd size mismatch");
    q_=std::move(q);qd_=std::move(qd);
}
Vec3 FlexibleBodyInterface::local_deformation(std::size_t node) const noexcept {Vec3 d{};for(std::size_t m=0;m<modes_.size();++m)d+=modes_[m].shape_local[node]*q_[m];return d;}
Vec3 FlexibleBodyInterface::local_deformation_velocity(std::size_t node) const noexcept {Vec3 d{};for(std::size_t m=0;m<modes_.size();++m)d+=modes_[m].shape_local[node]*qd_[m];return d;}
std::vector<Vec3> FlexibleBodyInterface::node_positions(const RigidBody& frame) const {std::vector<Vec3> out(nodes_.size());for(std::size_t i=0;i<nodes_.size();++i)out[i]=frame.state.position+rotate(frame.state.orientation,nodes_[i].reference_local+local_deformation(i));return out;}
std::vector<Vec3> FlexibleBodyInterface::node_velocities(const RigidBody& frame) const {std::vector<Vec3> out(nodes_.size());for(std::size_t i=0;i<nodes_.size();++i){const Vec3 r=rotate(frame.state.orientation,nodes_[i].reference_local+local_deformation(i));out[i]=frame.state.linear_velocity+cross(frame.state.angular_velocity,r)+rotate(frame.state.orientation,local_deformation_velocity(i));}return out;}
FlexibleInterfaceLoad FlexibleBodyInterface::project_nodal_forces(const RigidBody& frame,std::span<const Vec3> nodal_forces) const {
    if(nodal_forces.size()!=nodes_.size()) throw std::invalid_argument("flexible nodal force size mismatch");
    FlexibleInterfaceLoad load;load.modal_force.assign(modes_.size(),0.0);const auto positions=node_positions(frame);
    for(std::size_t i=0;i<nodes_.size();++i){
        const Vec3 f=nodal_forces[i];load.body_force+=f;load.body_torque+=cross(positions[i]-frame.state.position,f);
        for(std::size_t m=0;m<modes_.size();++m) load.modal_force[m]+=dot(rotate(frame.state.orientation,modes_[m].shape_local[i]),f);
    }
    return load;
}
FlexibleInterfaceLoad FlexibleBodyInterface::apply_nodal_forces(RigidBody& frame,std::span<const Vec3> nodal_forces) const {auto load=project_nodal_forces(frame,nodal_forces);frame.add_force(load.body_force);frame.add_torque(load.body_torque);return load;}

} // namespace cfd::multibody
