#include "cfd/multibody/system.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace cfd::multibody {

namespace {
[[nodiscard]] double clamp_value(double x,double lo,double hi) noexcept { return std::max(lo,std::min(hi,x)); }
[[nodiscard]] Vec3 body_point(const RigidBody& b,Vec3 local) { return b.state.position+rotate(b.state.orientation,local); }
[[nodiscard]] Vec3 local_point(const RigidBody& b,Vec3 world) { return inverse_rotate(b.state.orientation,world-b.state.position); }
[[nodiscard]] Vec3 world_axis(const RigidBody& b,Vec3 local) { return normalized(rotate(b.state.orientation,local)); }
[[nodiscard]] Quaternion relative_orientation(const RigidBody& a,const RigidBody& b) {
    return normalized(conjugate(normalized(a.state.orientation))*normalized(b.state.orientation));
}
[[nodiscard]] Vec3 orientation_error(Quaternion current,Quaternion reference) {
    Quaternion e=normalized(current*conjugate(normalized(reference)));
    if(e.w<0.0){ e.w=-e.w;e.x=-e.x;e.y=-e.y;e.z=-e.z; }
    return {2.0*e.x,2.0*e.y,2.0*e.z};
}
[[nodiscard]] Vec3 velocity_at_contact(const std::vector<RigidBody>& bodies,std::size_t id,Vec3 p) {
    return id==invalid_body?Vec3{}:bodies.at(id).point_velocity(p);
}
[[nodiscard]] double contact_effective_mass(const std::vector<RigidBody>& bodies,const ContactPoint& c,Vec3 dir) {
    double k=0.0;
    const auto accumulate=[&](std::size_t id,double sign){
        if(id==invalid_body) return;
        const auto& b=bodies.at(id);
        if(b.fixed) return;
        const Vec3 r=c.point-b.state.position;
        const Vec3 ang=cross(r,sign*dir);
        k+=b.inverse_mass()*dot(dir,dir)+dot(ang,b.world_inverse_inertia_mul(ang));
    };
    accumulate(c.body_a,1.0);
    accumulate(c.body_b,-1.0);
    return k;
}
void apply_contact_impulse(std::vector<RigidBody>& bodies,const ContactPoint& c,Vec3 impulse_on_a) {
    if(c.body_a!=invalid_body) bodies.at(c.body_a).apply_impulse( impulse_on_a,c.point);
    if(c.body_b!=invalid_body) bodies.at(c.body_b).apply_impulse(-impulse_on_a,c.point);
}
[[nodiscard]] Vec3 angular_velocity_relative(const std::vector<RigidBody>& bodies,const ContactPoint& c) {
    const Vec3 wa=c.body_a==invalid_body?Vec3{}:bodies.at(c.body_a).state.angular_velocity;
    const Vec3 wb=c.body_b==invalid_body?Vec3{}:bodies.at(c.body_b).state.angular_velocity;
    return wa-wb;
}
[[nodiscard]] double sphere_radius_for_body(const std::vector<SphereShape>& spheres,std::size_t body) {
    for(const auto& s:spheres) if(s.body==body) return s.radius;
    return 0.0;
}
[[nodiscard]] std::array<Vec3,3> basis_axes() noexcept { return {Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}}; }
void add_point_rows(std::vector<ConstraintRow>& rows,const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,
                    Vec3 pa,Vec3 pb,double dt,double beta,const std::array<Vec3,3>& axes=basis_axes()) {
    const Vec3 ra=pa-bodies.at(a).state.position;
    const Vec3 rb=pb-bodies.at(b).state.position;
    const Vec3 err=pb-pa;
    for(Vec3 n:axes){
        n=normalized(n);
        ConstraintRow r;
        r.body_a=a;r.body_b=b;
        r.linear_a=-n;r.angular_a=-cross(ra,n);
        r.linear_b=n;r.angular_b=cross(rb,n);
        r.bias=(beta/dt)*dot(err,n);
        rows.push_back(r);
    }
}
void add_angular_lock_rows(std::vector<ConstraintRow>& rows,const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,
                           Quaternion reference,double dt,double beta,const std::vector<Vec3>& world_axes) {
    const Vec3 local_err=orientation_error(relative_orientation(bodies.at(a),bodies.at(b)),reference);
    const Vec3 world_err=rotate(bodies.at(a).state.orientation,local_err);
    for(Vec3 n:world_axes){
        n=normalized(n);
        ConstraintRow r;
        r.body_a=a;r.body_b=b;
        r.angular_a=-n;r.angular_b=n;
        r.bias=(beta/dt)*dot(world_err,n);
        rows.push_back(r);
    }
}
[[nodiscard]] double constraint_velocity(const std::vector<RigidBody>& bodies,const ConstraintRow& r) {
    double v=0.0;
    if(r.body_a!=invalid_body){ const auto& b=bodies.at(r.body_a); v+=dot(r.linear_a,b.state.linear_velocity)+dot(r.angular_a,b.state.angular_velocity); }
    if(r.body_b!=invalid_body){ const auto& b=bodies.at(r.body_b); v+=dot(r.linear_b,b.state.linear_velocity)+dot(r.angular_b,b.state.angular_velocity); }
    return v;
}
[[nodiscard]] double constraint_diag(const std::vector<RigidBody>& bodies,const ConstraintRow& r) {
    double d=0.0;
    if(r.body_a!=invalid_body){ const auto& b=bodies.at(r.body_a); d+=b.inverse_mass()*dot(r.linear_a,r.linear_a)+dot(r.angular_a,b.world_inverse_inertia_mul(r.angular_a)); }
    if(r.body_b!=invalid_body){ const auto& b=bodies.at(r.body_b); d+=b.inverse_mass()*dot(r.linear_b,r.linear_b)+dot(r.angular_b,b.world_inverse_inertia_mul(r.angular_b)); }
    return d;
}
void apply_constraint_delta(std::vector<RigidBody>& bodies,const ConstraintRow& r,double dl) {
    if(r.body_a!=invalid_body){
        auto& b=bodies.at(r.body_a);
        if(!b.fixed){ b.state.linear_velocity+=r.linear_a*(dl*b.inverse_mass()); b.state.angular_velocity+=b.world_inverse_inertia_mul(r.angular_a*dl); }
    }
    if(r.body_b!=invalid_body){
        auto& b=bodies.at(r.body_b);
        if(!b.fixed){ b.state.linear_velocity+=r.linear_b*(dl*b.inverse_mass()); b.state.angular_velocity+=b.world_inverse_inertia_mul(r.angular_b*dl); }
    }
}


struct SupportPoint {
    Vec3 point{};
    Vec3 point_a{};
    Vec3 point_b{};
};

[[nodiscard]] Vec3 convex_support_world(const RigidBody& body,const ConvexHullShape& shape,Vec3 direction) {
    if(shape.vertices.empty()) throw std::invalid_argument("convex hull must contain vertices");
    const Vec3 local_direction=inverse_rotate(body.state.orientation,direction);
    const Vec3* best=&shape.vertices.front();
    double best_dot=dot(*best,local_direction);
    for(const auto& vertex:shape.vertices){
        const double value=dot(vertex,local_direction);
        if(value>best_dot){ best=&vertex;best_dot=value; }
    }
    return body.state.position+rotate(body.state.orientation,*best);
}

[[nodiscard]] SupportPoint minkowski_support(const std::vector<RigidBody>& bodies,const ConvexHullShape& a,const ConvexHullShape& b,Vec3 direction) {
    if(norm_squared(direction)<1.0e-28) direction={1.0,0.0,0.0};
    const Vec3 pa=convex_support_world(bodies.at(a.body),a,direction);
    const Vec3 pb=convex_support_world(bodies.at(b.body),b,-direction);
    return {pa-pb,pa,pb};
}

[[nodiscard]] Vec3 triple_cross(Vec3 a,Vec3 b,Vec3 c) noexcept { return cross(cross(a,b),c); }
[[nodiscard]] bool same_direction(Vec3 a,Vec3 b) noexcept { return dot(a,b)>0.0; }

[[nodiscard]] bool update_gjk_simplex(std::vector<SupportPoint>& simplex,Vec3& direction) {
    const SupportPoint a=simplex.back();
    const Vec3 ao=-a.point;
    if(simplex.size()==2U){
        const SupportPoint b=simplex[0];
        const Vec3 ab=b.point-a.point;
        if(same_direction(ab,ao)){
            direction=triple_cross(ab,ao,ab);
            if(norm_squared(direction)<1.0e-28) direction=any_perpendicular(ab);
        }else{
            simplex={a};direction=ao;
        }
        return false;
    }
    if(simplex.size()==3U){
        const SupportPoint b=simplex[1],c=simplex[0];
        const Vec3 ab=b.point-a.point,ac=c.point-a.point;
        const Vec3 abc=cross(ab,ac);
        const Vec3 ac_out=cross(abc,ac);
        if(same_direction(ac_out,ao)){
            if(same_direction(ac,ao)){
                simplex={c,a};direction=triple_cross(ac,ao,ac);
            }else if(same_direction(ab,ao)){
                simplex={b,a};direction=triple_cross(ab,ao,ab);
            }else{
                simplex={a};direction=ao;
            }
            return false;
        }
        const Vec3 ab_out=cross(ab,abc);
        if(same_direction(ab_out,ao)){
            if(same_direction(ab,ao)){
                simplex={b,a};direction=triple_cross(ab,ao,ab);
            }else{
                simplex={a};direction=ao;
            }
            return false;
        }
        if(same_direction(abc,ao)) direction=abc;
        else { simplex={b,c,a};direction=-abc; }
        return false;
    }
    if(simplex.size()==4U){
        const SupportPoint b=simplex[2],c=simplex[1],d=simplex[0];
        const auto outside_face=[&](SupportPoint p1,SupportPoint p2,SupportPoint opposite,std::vector<SupportPoint> kept)->bool{
            Vec3 n=cross(p1.point-a.point,p2.point-a.point);
            if(dot(n,opposite.point-a.point)>0.0) n=-n;
            if(dot(n,ao)>0.0){ simplex=std::move(kept);direction=n;return true; }
            return false;
        };
        if(outside_face(b,c,d,{c,b,a})) return false;
        if(outside_face(c,d,b,{d,c,a})) return false;
        if(outside_face(d,b,c,{b,d,a})) return false;
        return true;
    }
    direction=ao;
    return false;
}

struct EpaFace {
    std::size_t a{},b{},c{};
    Vec3 normal{};
    double distance{};
};

[[nodiscard]] EpaFace make_epa_face(const std::vector<SupportPoint>& vertices,std::size_t ia,std::size_t ib,std::size_t ic) {
    EpaFace f{ia,ib,ic,{},{}};
    Vec3 n=cross(vertices[ib].point-vertices[ia].point,vertices[ic].point-vertices[ia].point);
    const double length=norm(n);
    if(length<1.0e-14) return {ia,ib,ic,{1.0,0.0,0.0},std::numeric_limits<double>::infinity()};
    n=n/length;
    double d=dot(n,vertices[ia].point);
    if(d<0.0){ std::swap(f.b,f.c);n=-n;d=-d; }
    f.normal=n;f.distance=d;return f;
}

[[nodiscard]] std::array<double,3> barycentric_triangle(Vec3 p,Vec3 a,Vec3 b,Vec3 c) noexcept {
    const Vec3 v0=b-a,v1=c-a,v2=p-a;
    const double d00=dot(v0,v0),d01=dot(v0,v1),d11=dot(v1,v1),d20=dot(v2,v0),d21=dot(v2,v1);
    const double denom=d00*d11-d01*d01;
    if(std::abs(denom)<1.0e-24) return {1.0,0.0,0.0};
    const double v=(d11*d20-d01*d21)/denom;
    const double w=(d00*d21-d01*d20)/denom;
    return {1.0-v-w,v,w};
}

[[nodiscard]] bool same_contact_key(const ContactPoint& a,const ContactPoint& b) noexcept {
    return a.body_a==b.body_a && a.body_b==b.body_b && a.feature_id==b.feature_id;
}

[[nodiscard]] Vec3 triangle_vertex_world(const std::vector<RigidBody>& bodies,const TriangleMeshShape& mesh,Vec3 local) {
    if(mesh.body==invalid_body) return local;
    const auto& body=bodies.at(mesh.body);
    return body.state.position+rotate(body.state.orientation,local);
}

[[nodiscard]] Vec3 closest_point_on_triangle(Vec3 p,Vec3 a,Vec3 b,Vec3 c) noexcept {
    const Vec3 ab=b-a,ac=c-a,ap=p-a;
    const double d1=dot(ab,ap),d2=dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return a;

    const Vec3 bp=p-b;
    const double d3=dot(ab,bp),d4=dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return b;

    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0){
        const double v=d1/(d1-d3);
        return a+ab*v;
    }

    const Vec3 cp=p-c;
    const double d5=dot(ab,cp),d6=dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return c;

    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0){
        const double w=d2/(d2-d6);
        return a+ac*w;
    }

    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0){
        const Vec3 bc=c-b;
        const double w=(d4-d3)/((d4-d3)+(d5-d6));
        return b+bc*w;
    }

    const double denom=1.0/(va+vb+vc);
    const double v=vb*denom,w=vc*denom;
    return a+ab*v+ac*w;
}
}

Vec3 RigidBody::world_inverse_inertia_mul(Vec3 world_vector) const {
    if(fixed) return {};
    const Vec3 local=inverse_rotate(state.orientation,world_vector);
    const Vec3 local_result{
        local.x/inertia_diagonal.x,
        local.y/inertia_diagonal.y,
        local.z/inertia_diagonal.z
    };
    return rotate(state.orientation,local_result);
}
Vec3 RigidBody::point_velocity(Vec3 world_point) const noexcept {
    return state.linear_velocity+cross(state.angular_velocity,world_point-state.position);
}
void RigidBody::add_force_at_point(Vec3 f,Vec3 world_point) noexcept { force+=f;torque+=cross(world_point-state.position,f); }
void RigidBody::apply_impulse(Vec3 impulse,Vec3 world_point) noexcept {
    if(fixed) return;
    state.linear_velocity+=impulse*inverse_mass();
    state.angular_velocity+=world_inverse_inertia_mul(cross(world_point-state.position,impulse));
}
void RigidBody::validate() const {
    if(!(mass>0.0) || !std::isfinite(mass)) throw std::invalid_argument("rigid body mass must be finite and positive");
    if(!(inertia_diagonal.x>0.0 && inertia_diagonal.y>0.0 && inertia_diagonal.z>0.0)) throw std::invalid_argument("rigid body inertia must be positive");
    if(!std::isfinite(norm_squared(state.position)) || !std::isfinite(norm_squared(state.linear_velocity)) || !std::isfinite(norm_squared(state.angular_velocity))) throw std::invalid_argument("rigid body state must be finite");
    (void)normalized(state.orientation);
}
void integrate_rigid_body(RigidBody& body,double dt,RigidBodyIntegrator integrator) {
    if(body.fixed){ body.clear_accumulators(); return; }
    body.validate();
    if(!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("rigid body dt must be finite and positive");
    const Vec3 a=body.force*body.inverse_mass();
    const Vec3 alpha=body.world_inverse_inertia_mul(body.torque-cross(body.state.angular_velocity,rotate(body.state.orientation,component_mul(body.inertia_diagonal,inverse_rotate(body.state.orientation,body.state.angular_velocity)))));
    if(integrator==RigidBodyIntegrator::velocity_verlet){
        body.state.position+=body.state.linear_velocity*dt+a*(0.5*dt*dt);
        body.state.linear_velocity+=a*dt;
        body.state.angular_velocity+=alpha*dt;
        integrate_orientation(body.state.orientation,body.state.angular_velocity-alpha*(0.5*dt),dt);
    } else {
        body.state.linear_velocity+=a*dt;
        body.state.angular_velocity+=alpha*dt;
        body.state.position+=body.state.linear_velocity*dt;
        integrate_orientation(body.state.orientation,body.state.angular_velocity,dt);
    }
    body.clear_accumulators();
}

DistanceConstraint make_distance_constraint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 wa,Vec3 wb,double distance) {
    const double d=distance>=0.0?distance:norm(wb-wa);
    return {a,b,local_point(bodies.at(a),wa),local_point(bodies.at(b),wb),d};
}
SphericalJoint make_spherical_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 w) {
    return {a,b,local_point(bodies.at(a),w),local_point(bodies.at(b),w)};
}
RevoluteJoint make_revolute_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 w,Vec3 axis) {
    axis=normalized(axis);
    return {a,b,local_point(bodies.at(a),w),local_point(bodies.at(b),w),inverse_rotate(bodies.at(a).state.orientation,axis),inverse_rotate(bodies.at(b).state.orientation,axis)};
}
PrismaticJoint make_prismatic_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 w,Vec3 axis) {
    axis=normalized(axis);
    return {a,b,local_point(bodies.at(a),w),local_point(bodies.at(b),w),inverse_rotate(bodies.at(a).state.orientation,axis),relative_orientation(bodies.at(a),bodies.at(b))};
}
FixedJoint make_fixed_joint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 w) {
    return {a,b,local_point(bodies.at(a),w),local_point(bodies.at(b),w),relative_orientation(bodies.at(a),bodies.at(b))};
}
GearConstraint make_gear_constraint(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 aa,Vec3 ab,double ratio) {
    if(!std::isfinite(ratio) || ratio==0.0) throw std::invalid_argument("gear ratio must be finite and non-zero");
    aa=normalized(aa);ab=normalized(ab);
    return {a,b,inverse_rotate(bodies.at(a).state.orientation,aa),inverse_rotate(bodies.at(b).state.orientation,ab),ratio};
}
LinearMotor make_linear_motor(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 axis,double target_speed,double max_force) {
    if(!std::isfinite(target_speed) || !(max_force>=0.0) || std::isnan(max_force)) throw std::invalid_argument("invalid linear motor parameters");
    axis=normalized(axis);
    return {a,b,inverse_rotate(bodies.at(a).state.orientation,axis),target_speed,max_force};
}
AngularMotor make_angular_motor(const std::vector<RigidBody>& bodies,std::size_t a,std::size_t b,Vec3 axis,double target_speed,double max_torque) {
    if(!std::isfinite(target_speed) || !(max_torque>=0.0) || std::isnan(max_torque)) throw std::invalid_argument("invalid angular motor parameters");
    axis=normalized(axis);
    return {a,b,inverse_rotate(bodies.at(a).state.orientation,axis),target_speed,max_torque};
}

std::vector<ConstraintRow> build_constraint_rows(const std::vector<RigidBody>& bodies,const Joint& joint,double dt,double beta) {
    if(!(dt>0.0) || !(beta>=0.0)) throw std::invalid_argument("invalid constraint timestep/Baumgarte factor");
    std::vector<ConstraintRow> rows;
    std::visit([&](const auto& j){
        using T=std::decay_t<decltype(j)>;
        const auto& a=bodies.at(j.body_a); const auto& b=bodies.at(j.body_b);
        if constexpr(std::is_same_v<T,DistanceConstraint>){
            const Vec3 pa=body_point(a,j.local_anchor_a),pb=body_point(b,j.local_anchor_b);
            const Vec3 d=pb-pa; const double l=norm(d); const Vec3 n=l>1e-12?d/l:Vec3{1,0,0};
            ConstraintRow r; r.body_a=j.body_a;r.body_b=j.body_b;
            const Vec3 ra=pa-a.state.position,rb=pb-b.state.position;
            r.linear_a=-n;r.angular_a=-cross(ra,n);r.linear_b=n;r.angular_b=cross(rb,n);
            r.bias=(beta/dt)*(l-j.distance); rows.push_back(r);
        } else if constexpr(std::is_same_v<T,SphericalJoint>){
            add_point_rows(rows,bodies,j.body_a,j.body_b,body_point(a,j.local_anchor_a),body_point(b,j.local_anchor_b),dt,beta);
        } else if constexpr(std::is_same_v<T,RevoluteJoint>){
            add_point_rows(rows,bodies,j.body_a,j.body_b,body_point(a,j.local_anchor_a),body_point(b,j.local_anchor_b),dt,beta);
            const Vec3 aa=world_axis(a,j.local_axis_a),ab=world_axis(b,j.local_axis_b);
            Vec3 t1=any_perpendicular(aa),t2=normalized(cross(aa,t1));
            const Vec3 err=cross(aa,ab);
            for(Vec3 t:{t1,t2}){ ConstraintRow r;r.body_a=j.body_a;r.body_b=j.body_b;r.angular_a=-t;r.angular_b=t;r.bias=(beta/dt)*dot(err,t);rows.push_back(r); }
        } else if constexpr(std::is_same_v<T,PrismaticJoint>){
            const Vec3 pa=body_point(a,j.local_anchor_a),pb=body_point(b,j.local_anchor_b),axis=world_axis(a,j.local_axis_a);
            const Vec3 t1=any_perpendicular(axis),t2=normalized(cross(axis,t1));
            const Vec3 ra=pa-a.state.position,rb=pb-b.state.position,err=pb-pa;
            for(Vec3 t:{t1,t2}){ ConstraintRow r;r.body_a=j.body_a;r.body_b=j.body_b;r.linear_a=-t;r.angular_a=-cross(ra,t);r.linear_b=t;r.angular_b=cross(rb,t);r.bias=(beta/dt)*dot(err,t);rows.push_back(r); }
            add_angular_lock_rows(rows,bodies,j.body_a,j.body_b,j.reference_relative_orientation,dt,beta,{Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}});
        } else if constexpr(std::is_same_v<T,FixedJoint>){
            add_point_rows(rows,bodies,j.body_a,j.body_b,body_point(a,j.local_anchor_a),body_point(b,j.local_anchor_b),dt,beta);
            add_angular_lock_rows(rows,bodies,j.body_a,j.body_b,j.reference_relative_orientation,dt,beta,{Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}});
        } else if constexpr(std::is_same_v<T,GearConstraint>){
            ConstraintRow r;r.body_a=j.body_a;r.body_b=j.body_b;
            r.angular_a=world_axis(a,j.local_axis_a);r.angular_b=world_axis(b,j.local_axis_b)*j.ratio;rows.push_back(r);
        } else if constexpr(std::is_same_v<T,LinearMotor>){
            const Vec3 axis=world_axis(a,j.local_axis_a);
            ConstraintRow r;r.body_a=j.body_a;r.body_b=j.body_b;
            r.linear_a=-axis;r.linear_b=axis;r.bias=-j.target_speed;
            const double limit=std::isfinite(j.max_force)?j.max_force*dt:std::numeric_limits<double>::infinity();
            r.lower=-limit;r.upper=limit;rows.push_back(r);
        } else if constexpr(std::is_same_v<T,AngularMotor>){
            const Vec3 axis=world_axis(a,j.local_axis_a);
            ConstraintRow r;r.body_a=j.body_a;r.body_b=j.body_b;
            r.angular_a=-axis;r.angular_b=axis;r.bias=-j.target_speed;
            const double limit=std::isfinite(j.max_torque)?j.max_torque*dt:std::numeric_limits<double>::infinity();
            r.lower=-limit;r.upper=limit;rows.push_back(r);
        }
    },joint);
    return rows;
}
void solve_constraint_rows(std::vector<RigidBody>& bodies,std::vector<ConstraintRow>& rows,std::size_t iterations) {
    for(std::size_t it=0;it<iterations;++it){
        for(auto& r:rows){
            const double d=constraint_diag(bodies,r); if(d<=1e-18) continue;
            const double old=r.impulse;
            const double dl=-(constraint_velocity(bodies,r)+r.bias)/d;
            r.impulse=clamp_value(old+dl,r.lower,r.upper);
            apply_constraint_delta(bodies,r,r.impulse-old);
        }
    }
}

Aabb sphere_aabb(const RigidBody& body,const SphereShape& sphere) {
    if(!(sphere.radius>0.0)) throw std::invalid_argument("sphere radius must be positive");
    const Vec3 r{sphere.radius,sphere.radius,sphere.radius}; return {body.state.position-r,body.state.position+r};
}
Aabb convex_aabb(const RigidBody& body,const ConvexHullShape& shape) {
    if(shape.vertices.empty()) throw std::invalid_argument("convex hull must contain vertices");
    Vec3 first=body.state.position+rotate(body.state.orientation,shape.vertices.front());
    Aabb box{first,first};
    for(const auto& local:shape.vertices){
        const Vec3 p=body.state.position+rotate(body.state.orientation,local);
        box.minimum={std::min(box.minimum.x,p.x),std::min(box.minimum.y,p.y),std::min(box.minimum.z,p.z)};
        box.maximum={std::max(box.maximum.x,p.x),std::max(box.maximum.y,p.y),std::max(box.maximum.z,p.z)};
    }
    return box;
}

GjkEpaResult convex_gjk_epa_contact(const std::vector<RigidBody>& bodies,const ConvexHullShape& a,const ConvexHullShape& b,std::size_t max_gjk_iterations,std::size_t max_epa_iterations,double tolerance) {
    if(a.body>=bodies.size() || b.body>=bodies.size() || a.body==b.body) throw std::invalid_argument("invalid convex hull body attachment");
    if(a.vertices.size()<4U || b.vertices.size()<4U) throw std::invalid_argument("GJK/EPA convex hull requires at least four vertices per shape");
    if(max_gjk_iterations==0U || max_epa_iterations==0U || !(tolerance>0.0)) throw std::invalid_argument("invalid GJK/EPA iteration/tolerance settings");

    Vec3 direction=bodies.at(a.body).state.position-bodies.at(b.body).state.position;
    if(norm_squared(direction)<1.0e-24) direction={1.0,0.0,0.0};
    std::vector<SupportPoint> simplex;
    simplex.reserve(4);
    simplex.push_back(minkowski_support(bodies,a,b,direction));
    direction=-simplex.back().point;

    std::size_t gjk_iterations=0U;
    bool intersect=false;
    for(;gjk_iterations<max_gjk_iterations;++gjk_iterations){
        if(norm_squared(direction)<tolerance*tolerance){ intersect=true;break; }
        const SupportPoint pnt=minkowski_support(bodies,a,b,direction);
        if(dot(pnt.point,direction)<0.0) return {false,{},gjk_iterations+1U,0U};
        simplex.push_back(pnt);
        if(update_gjk_simplex(simplex,direction)){ intersect=true;++gjk_iterations;break; }
    }
    if(!intersect || simplex.size()<4U) return {false,{},gjk_iterations,0U};

    std::vector<SupportPoint> vertices=simplex;
    std::vector<EpaFace> faces;
    faces.reserve(64);
    faces.push_back(make_epa_face(vertices,0,1,2));
    faces.push_back(make_epa_face(vertices,0,3,1));
    faces.push_back(make_epa_face(vertices,0,2,3));
    faces.push_back(make_epa_face(vertices,1,3,2));

    const auto finish_contact=[&](const EpaFace& face,std::size_t epa_iterations)->GjkEpaResult{
        const Vec3 projected=face.normal*face.distance;
        auto w=barycentric_triangle(projected,vertices[face.a].point,vertices[face.b].point,vertices[face.c].point);
        for(double& wi:w) wi=clamp_value(wi,0.0,1.0);
        const double ws=w[0]+w[1]+w[2];
        if(ws>1.0e-15){ w[0]/=ws;w[1]/=ws;w[2]/=ws; }
        const Vec3 pa=vertices[face.a].point_a*w[0]+vertices[face.b].point_a*w[1]+vertices[face.c].point_a*w[2];
        const Vec3 pb=vertices[face.a].point_b*w[0]+vertices[face.b].point_b*w[1]+vertices[face.c].point_b*w[2];
        ContactPoint contact;
        contact.body_a=a.body;contact.body_b=b.body;
        contact.point=(pa+pb)*0.5;
        contact.normal=-face.normal;
        contact.penetration=std::max(0.0,face.distance);
        return {true,contact,gjk_iterations,epa_iterations};
    };

    for(std::size_t iter=0;iter<max_epa_iterations;++iter){
        auto closest=std::min_element(faces.begin(),faces.end(),[](const EpaFace& x,const EpaFace& y){return x.distance<y.distance;});
        if(closest==faces.end() || !std::isfinite(closest->distance)) break;
        const EpaFace best=*closest;
        const SupportPoint pnt=minkowski_support(bodies,a,b,best.normal);
        const double support_distance=dot(pnt.point,best.normal);
        if(support_distance-best.distance<=tolerance) return finish_contact(best,iter+1U);

        const std::size_t new_index=vertices.size();
        vertices.push_back(pnt);
        std::vector<std::pair<std::size_t,std::size_t>> horizon;
        const auto add_edge=[&](std::size_t x,std::size_t y){
            const auto rev=std::find(horizon.begin(),horizon.end(),std::pair<std::size_t,std::size_t>{y,x});
            if(rev!=horizon.end()) horizon.erase(rev); else horizon.emplace_back(x,y);
        };
        std::vector<EpaFace> kept;
        kept.reserve(faces.size()+8U);
        for(const auto& face:faces){
            if(dot(face.normal,pnt.point-vertices[face.a].point)>tolerance){
                add_edge(face.a,face.b);add_edge(face.b,face.c);add_edge(face.c,face.a);
            }else kept.push_back(face);
        }
        if(horizon.empty()) return finish_contact(best,iter+1U);
        for(const auto& edge:horizon){
            EpaFace f=make_epa_face(vertices,edge.first,edge.second,new_index);
            if(std::isfinite(f.distance)) kept.push_back(f);
        }
        faces=std::move(kept);
    }
    if(!faces.empty()){
        const auto closest=std::min_element(faces.begin(),faces.end(),[](const EpaFace& x,const EpaFace& y){return x.distance<y.distance;});
        if(closest!=faces.end() && std::isfinite(closest->distance)) return finish_contact(*closest,max_epa_iterations);
    }
    return {false,{},gjk_iterations,max_epa_iterations};
}

std::vector<std::pair<std::size_t,std::size_t>> broad_phase_sweep_and_prune(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres) {
    struct Entry{std::size_t sphere{};Aabb box{};};
    std::vector<Entry> e; e.reserve(spheres.size());
    for(std::size_t i=0;i<spheres.size();++i) e.push_back({i,sphere_aabb(bodies.at(spheres[i].body),spheres[i])});
    std::sort(e.begin(),e.end(),[](const Entry& a,const Entry& b){return a.box.minimum.x<b.box.minimum.x;});
    std::vector<std::pair<std::size_t,std::size_t>> out;
    for(std::size_t i=0;i<e.size();++i) for(std::size_t j=i+1;j<e.size() && e[j].box.minimum.x<=e[i].box.maximum.x;++j){
        const auto& a=e[i].box;const auto& b=e[j].box;
        if(a.maximum.y<b.minimum.y||b.maximum.y<a.minimum.y||a.maximum.z<b.minimum.z||b.maximum.z<a.minimum.z) continue;
        if(spheres[e[i].sphere].body==spheres[e[j].sphere].body) continue;
        out.emplace_back(e[i].sphere,e[j].sphere);
    }
    return out;
}

std::vector<std::pair<std::size_t,std::size_t>> broad_phase_bvh(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres) {
    struct Node { Aabb box{}; int left{-1}; int right{-1}; std::size_t begin{}; std::size_t end{}; };
    if(spheres.size()<2U) return {};
    std::vector<std::size_t> order(spheres.size());
    std::iota(order.begin(),order.end(),0U);
    std::vector<Aabb> boxes(spheres.size());
    for(std::size_t i=0;i<spheres.size();++i) boxes[i]=sphere_aabb(bodies.at(spheres[i].body),spheres[i]);
    auto merged=[](Aabb a,Aabb b){return Aabb{{std::min(a.minimum.x,b.minimum.x),std::min(a.minimum.y,b.minimum.y),std::min(a.minimum.z,b.minimum.z)},
                                                {std::max(a.maximum.x,b.maximum.x),std::max(a.maximum.y,b.maximum.y),std::max(a.maximum.z,b.maximum.z)}};};
    auto overlaps=[](Aabb a,Aabb b){return !(a.maximum.x<b.minimum.x||b.maximum.x<a.minimum.x||a.maximum.y<b.minimum.y||b.maximum.y<a.minimum.y||a.maximum.z<b.minimum.z||b.maximum.z<a.minimum.z);};
    std::vector<Node> nodes;nodes.reserve(2*spheres.size());
    const auto build=[&](auto&& self,std::size_t begin,std::size_t end)->int{
        Aabb box=boxes[order[begin]];for(std::size_t i=begin+1;i<end;++i)box=merged(box,boxes[order[i]]);
        const int id=static_cast<int>(nodes.size());nodes.push_back({box,-1,-1,begin,end});
        if(end-begin<=2U)return id;
        const Vec3 ext=box.maximum-box.minimum;const int axis=(ext.y>ext.x&&ext.y>=ext.z)?1:((ext.z>ext.x&&ext.z>ext.y)?2:0);
        const std::size_t mid=begin+(end-begin)/2U;
        auto center=[&](std::size_t si){const auto& b=boxes[si];return axis==0?0.5*(b.minimum.x+b.maximum.x):(axis==1?0.5*(b.minimum.y+b.maximum.y):0.5*(b.minimum.z+b.maximum.z));};
        std::nth_element(order.begin()+static_cast<std::ptrdiff_t>(begin),order.begin()+static_cast<std::ptrdiff_t>(mid),order.begin()+static_cast<std::ptrdiff_t>(end),[&](std::size_t a,std::size_t b){return center(a)<center(b);});
        nodes[id].left=self(self,begin,mid);nodes[id].right=self(self,mid,end);return id;
    };
    const int root=build(build,0U,order.size());
    std::vector<std::pair<std::size_t,std::size_t>> out;
    const auto emit=[&](auto&& self,int na,int nb){
        if(na<0||nb<0||!overlaps(nodes[static_cast<std::size_t>(na)].box,nodes[static_cast<std::size_t>(nb)].box))return;
        const Node& a=nodes[static_cast<std::size_t>(na)];const Node& b=nodes[static_cast<std::size_t>(nb)];
        const bool la=a.left<0,lb=b.left<0;
        if(la&&lb){for(std::size_t ia=a.begin;ia<a.end;++ia)for(std::size_t ib=b.begin;ib<b.end;++ib){const std::size_t sa=order[ia],sb=order[ib];if(sa>=sb)continue;if(spheres[sa].body==spheres[sb].body)continue;if(overlaps(boxes[sa],boxes[sb]))out.emplace_back(sa,sb);}return;}
        if(na==nb){self(self,a.left,a.left);self(self,a.left,a.right);self(self,a.right,a.right);return;}
        if(lb||(!la&&(a.end-a.begin)>=(b.end-b.begin))){self(self,a.left,nb);self(self,a.right,nb);}else{self(self,na,b.left);self(self,na,b.right);}
    };
    emit(emit,root,root);
    std::sort(out.begin(),out.end());out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}

std::vector<std::pair<std::size_t,std::size_t>> broad_phase_cell_linked(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,double cell_size) {
    if(!(cell_size>0.0) || !std::isfinite(cell_size)) throw std::invalid_argument("cell-linked broad-phase cell size must be finite and positive");
    if(spheres.empty()) return {};
    Vec3 lo{std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
    Vec3 hi{-lo.x,-lo.y,-lo.z};
    for(const auto& sphere:spheres){
        if(sphere.body>=bodies.size() || !(sphere.radius>0.0) || 2.0*sphere.radius>cell_size) throw std::invalid_argument("cell-linked broad phase requires valid spheres with diameter <= cell size");
        const auto& p=bodies[sphere.body].state.position;
        lo.x=std::min(lo.x,p.x);lo.y=std::min(lo.y,p.y);lo.z=std::min(lo.z,p.z);
        hi.x=std::max(hi.x,p.x);hi.y=std::max(hi.y,p.y);hi.z=std::max(hi.z,p.z);
    }
    lo-=Vec3{cell_size,cell_size,cell_size};hi+=Vec3{cell_size,cell_size,cell_size};
    const auto extent=[&](double a,double b){return std::max<std::size_t>(1U,static_cast<std::size_t>(std::ceil((b-a)/cell_size)));};
    const std::size_t nx=extent(lo.x,hi.x),ny=extent(lo.y,hi.y),nz=extent(lo.z,hi.z),cells=nx*ny*nz;
    std::vector<std::vector<std::size_t>> buckets(cells);
    const auto cell_of=[&](Vec3 p){
        const auto axis=[&](double v,double origin,std::size_t n){const auto raw=static_cast<long long>((v-origin)/cell_size);return static_cast<std::size_t>(std::clamp<long long>(raw,0,static_cast<long long>(n)-1));};
        const std::size_t x=axis(p.x,lo.x,nx),y=axis(p.y,lo.y,ny),z=axis(p.z,lo.z,nz);return (z*ny+y)*nx+x;
    };
    for(std::size_t i=0;i<spheres.size();++i) buckets[cell_of(bodies[spheres[i].body].state.position)].push_back(i);
    std::vector<std::pair<std::size_t,std::size_t>> pairs;
    for(std::size_t z=0;z<nz;++z)for(std::size_t y=0;y<ny;++y)for(std::size_t x=0;x<nx;++x){
        const auto& source=buckets[(z*ny+y)*nx+x];
        for(const auto i:source){
            for(int dz=-1;dz<=1;++dz){const long long zz=static_cast<long long>(z)+dz;if(zz<0||zz>=static_cast<long long>(nz))continue;
                for(int dy=-1;dy<=1;++dy){const long long yy=static_cast<long long>(y)+dy;if(yy<0||yy>=static_cast<long long>(ny))continue;
                    for(int dx=-1;dx<=1;++dx){const long long xx=static_cast<long long>(x)+dx;if(xx<0||xx>=static_cast<long long>(nx))continue;
                        for(const auto j:buckets[(static_cast<std::size_t>(zz)*ny+static_cast<std::size_t>(yy))*nx+static_cast<std::size_t>(xx)]) if(j>i){
                            const auto ai=sphere_aabb(bodies[spheres[i].body],spheres[i]);const auto aj=sphere_aabb(bodies[spheres[j].body],spheres[j]);
                            const bool overlap=ai.minimum.x<=aj.maximum.x && ai.maximum.x>=aj.minimum.x && ai.minimum.y<=aj.maximum.y && ai.maximum.y>=aj.minimum.y && ai.minimum.z<=aj.maximum.z && ai.maximum.z>=aj.minimum.z;
                            if(overlap) pairs.emplace_back(i,j);
                        }
                    }
                }
            }
        }
    }
    std::sort(pairs.begin(),pairs.end());pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());return pairs;
}

std::vector<ContactPoint> sphere_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const std::vector<std::pair<std::size_t,std::size_t>>& pairs) {
    std::vector<ContactPoint> out;
    for(const auto& [ia,ib]:pairs){
        const auto& sa=spheres.at(ia);const auto& sb=spheres.at(ib);const auto& a=bodies.at(sa.body);const auto& b=bodies.at(sb.body);
        const Vec3 d=a.state.position-b.state.position;const double l=norm(d);const double pen=sa.radius+sb.radius-l;
        if(pen<=0.0) continue;
        const Vec3 n=l>1e-12?d/l:Vec3{1,0,0};
        const Vec3 p=a.state.position-n*(sa.radius-0.5*pen); out.push_back({sa.body,sb.body,p,n,pen,0U});
    }
    return out;
}
std::vector<ContactPoint> sphere_plane_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const std::vector<PlaneShape>& planes) {
    std::vector<ContactPoint> out;
    for(const auto& s:spheres){
        const auto& b=bodies.at(s.body);
        for(std::size_t pi=0;pi<planes.size();++pi){
            const auto& pr=planes[pi];const Vec3 n=normalized(pr.normal);const double gap=dot(n,b.state.position)-pr.offset-s.radius;
            if(gap<0.0) out.push_back({s.body,invalid_body,b.state.position-n*s.radius,n,-gap,static_cast<std::uint64_t>(pi)+1U});
        }
    }
    return out;
}

std::vector<ContactPoint> sphere_triangle_mesh_contacts(const std::vector<RigidBody>& bodies,const std::vector<SphereShape>& spheres,const TriangleMeshShape& mesh) {
    if(mesh.body!=invalid_body && mesh.body>=bodies.size()) throw std::invalid_argument("triangle mesh references invalid body");
    if(mesh.vertices.empty() || mesh.triangles.empty()) return {};
    std::vector<ContactPoint> out;
    for(std::size_t si=0;si<spheres.size();++si){
        const auto& sphere=spheres[si];
        if(sphere.body>=bodies.size()) throw std::invalid_argument("sphere references invalid body");
        if(mesh.body==sphere.body) continue;
        const Vec3 center=bodies.at(sphere.body).state.position;
        for(std::size_t ti=0;ti<mesh.triangles.size();++ti){
            const auto tri=mesh.triangles[ti];
            if(tri[0]>=mesh.vertices.size() || tri[1]>=mesh.vertices.size() || tri[2]>=mesh.vertices.size()) throw std::invalid_argument("triangle mesh index out of range");
            const Vec3 a=triangle_vertex_world(bodies,mesh,mesh.vertices[tri[0]]);
            const Vec3 b=triangle_vertex_world(bodies,mesh,mesh.vertices[tri[1]]);
            const Vec3 c=triangle_vertex_world(bodies,mesh,mesh.vertices[tri[2]]);
            const Vec3 cp=closest_point_on_triangle(center,a,b,c);
            const Vec3 delta=center-cp;
            const double distance=norm(delta);
            if(distance>=sphere.radius) continue;
            Vec3 normal{};
            if(distance>1.0e-12) normal=delta/distance;
            else {
                normal=normalized(cross(b-a,c-a));
                const Vec3 centroid=(a+b+c)/3.0;
                if(dot(normal,center-centroid)<0.0) normal=-normal;
            }
            const double penetration=sphere.radius-distance;
            const std::uint64_t feature=(static_cast<std::uint64_t>(mesh.feature_namespace)<<32U) | static_cast<std::uint64_t>(ti+1U);
            out.push_back({sphere.body,mesh.body,cp,normal,penetration,feature});
        }
    }
    return out;
}

void ContactManifoldCache::update(const std::vector<ContactPoint>& contacts) {
    std::vector<PersistentContactState> next;
    next.reserve(contacts.size());
    for(const auto& contact:contacts){
        auto it=std::find_if(states_.begin(),states_.end(),[&](const PersistentContactState& state){return same_contact_key(state.contact,contact);});
        PersistentContactState state;
        if(it!=states_.end()) state=*it;
        state.contact=contact;
        state.tangential_displacement-=contact.normal*dot(state.tangential_displacement,contact.normal);
        state.age=(it==states_.end())?1U:(it->age+1U);
        next.push_back(state);
    }
    states_=std::move(next);
}
void apply_penalty_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const PenaltyContactModel& m) {
    for(const auto& c:contacts){
        const Vec3 rel=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);const double vn=dot(rel,c.normal);
        const double fn=std::max(0.0,m.normal_stiffness*c.penetration-m.normal_damping*vn); const Vec3 vt=rel-c.normal*vn; Vec3 ft{};
        const double vtmag=norm(vt);if(vtmag>1e-14) ft=vt*(-std::min(m.friction*fn,m.tangential_damping*vtmag)/vtmag);
        const Vec3 f=c.normal*fn+ft; if(c.body_a!=invalid_body)bodies.at(c.body_a).add_force_at_point(f,c.point);if(c.body_b!=invalid_body)bodies.at(c.body_b).add_force_at_point(-f,c.point);
    }
}
void apply_hertz_mindlin_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const std::vector<SphereShape>& spheres,const HertzMindlinContactModel& m) {
    for(const auto& c:contacts){
        const Vec3 rel=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);const double vn=dot(rel,c.normal);const double root=std::sqrt(std::max(0.0,c.penetration));
        double fn=m.normal_stiffness*c.penetration*root-m.normal_damping*root*vn-m.cohesion_force; fn=std::max(-m.cohesion_force,fn);
        const Vec3 vt=rel-c.normal*vn;Vec3 ft{};const double vtmag=norm(vt);if(vtmag>1e-14)ft=vt*(-std::min(m.friction*std::abs(fn),m.tangential_damping*root*vtmag)/vtmag);
        const Vec3 f=c.normal*fn+ft;if(c.body_a!=invalid_body)bodies.at(c.body_a).add_force_at_point(f,c.point);if(c.body_b!=invalid_body)bodies.at(c.body_b).add_force_at_point(-f,c.point);
        const Vec3 wr=angular_velocity_relative(bodies,c);const double wmag=norm(wr);if(wmag>1e-14 && m.rolling_resistance>0.0){const double ra=sphere_radius_for_body(spheres,c.body_a);const double rb=sphere_radius_for_body(spheres,c.body_b);const double reff=(rb>0.0)?(ra*rb/(ra+rb)):ra;const Vec3 tr=wr*(-m.rolling_resistance*std::abs(fn)*reff/wmag);if(c.body_a!=invalid_body)bodies.at(c.body_a).add_torque(tr);if(c.body_b!=invalid_body)bodies.at(c.body_b).add_torque(-tr);}
    }
}

void apply_history_dependent_mindlin_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,const std::vector<SphereShape>& spheres,ContactManifoldCache& cache,double dt,const HertzMindlinContactModel& m) {
    if(!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("Mindlin contact dt must be finite and positive");
    if(!(m.tangential_stiffness>=0.0) || !std::isfinite(m.tangential_stiffness)) throw std::invalid_argument("Mindlin tangential stiffness must be finite and non-negative");
    cache.update(contacts);
    for(auto& state:cache.states()){
        const auto& c=state.contact;
        const Vec3 rel=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);
        const double vn=dot(rel,c.normal);
        const double root=std::sqrt(std::max(0.0,c.penetration));
        double fn=m.normal_stiffness*c.penetration*root-m.normal_damping*root*vn-m.cohesion_force;
        fn=std::max(-m.cohesion_force,fn);
        const Vec3 vt=rel-c.normal*vn;
        state.tangential_displacement+=vt*dt;
        state.tangential_displacement-=c.normal*dot(state.tangential_displacement,c.normal);
        const double kt=m.tangential_stiffness*root;
        const double ct=m.tangential_damping*root;
        Vec3 ft=state.tangential_displacement*(-kt)+vt*(-ct);
        const double max_ft=m.friction*std::abs(fn);
        const double ft_mag=norm(ft);
        if(ft_mag>max_ft && ft_mag>1.0e-14){
            ft*=max_ft/ft_mag;
            if(kt>1.0e-14) state.tangential_displacement=-(ft+vt*ct)/kt;
        }
        const Vec3 f=c.normal*fn+ft;
        if(c.body_a!=invalid_body)bodies.at(c.body_a).add_force_at_point(f,c.point);
        if(c.body_b!=invalid_body)bodies.at(c.body_b).add_force_at_point(-f,c.point);
        const Vec3 wr=angular_velocity_relative(bodies,c);const double wmag=norm(wr);
        if(wmag>1e-14 && m.rolling_resistance>0.0){
            const double ra=sphere_radius_for_body(spheres,c.body_a);const double rb=sphere_radius_for_body(spheres,c.body_b);
            const double reff=(rb>0.0)?(ra*rb/(ra+rb)):ra;
            const Vec3 tr=wr*(-m.rolling_resistance*std::abs(fn)*reff/wmag);
            if(c.body_a!=invalid_body)bodies.at(c.body_a).add_torque(tr);
            if(c.body_b!=invalid_body)bodies.at(c.body_b).add_torque(-tr);
        }
    }
}
ParticleBond make_particle_bond(const std::vector<RigidBody>& bodies,std::size_t body_a,std::size_t body_b,BondedParticleModel model,double rest_length) {
    if(body_a>=bodies.size() || body_b>=bodies.size() || body_a==body_b) throw std::invalid_argument("invalid bonded-particle body pair");
    const auto valid_nonnegative=[](double x){ return x>=0.0 && std::isfinite(x); };
    if(!valid_nonnegative(model.normal_stiffness) || !valid_nonnegative(model.shear_stiffness) ||
       !valid_nonnegative(model.normal_damping) || !valid_nonnegative(model.shear_damping) ||
       !(model.tensile_failure_force>0.0) || !std::isfinite(model.tensile_failure_force) ||
       !(model.shear_failure_force>0.0) || !std::isfinite(model.shear_failure_force) ||
       !(model.damage_onset_ratio>=0.0 && model.damage_onset_ratio<1.0) || !std::isfinite(model.damage_onset_ratio)) {
        throw std::invalid_argument("invalid bonded-particle constitutive parameters");
    }
    const double distance=norm(bodies[body_b].state.position-bodies[body_a].state.position);
    if(rest_length==0.0) rest_length=distance;
    if(!(rest_length>0.0) || !std::isfinite(rest_length)) throw std::invalid_argument("bond rest length must be finite and positive");
    ParticleBond bond;
    bond.body_a=body_a;
    bond.body_b=body_b;
    bond.rest_length=rest_length;
    bond.model=model;
    return bond;
}

BondUpdateStats apply_particle_bonds(std::vector<RigidBody>& bodies,std::vector<ParticleBond>& bonds,double dt) {
    if(!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("bond update dt must be finite and positive");
    BondUpdateStats stats;
    for(auto& bond:bonds){
        if(bond.body_a>=bodies.size() || bond.body_b>=bodies.size() || bond.body_a==bond.body_b) throw std::out_of_range("bond references invalid body");
        if(bond.broken){ ++stats.broken_total; continue; }
        ++stats.active;
        auto& a=bodies[bond.body_a];
        auto& b=bodies[bond.body_b];
        const Vec3 delta=b.state.position-a.state.position;
        const double distance=norm(delta);
        if(!(distance>1.0e-15) || !std::isfinite(distance)) throw std::runtime_error("bonded particle centers are coincident/non-finite");
        const Vec3 n=delta/distance;
        const Vec3 relative_velocity=b.state.linear_velocity-a.state.linear_velocity;
        const double vn=dot(relative_velocity,n);
        const Vec3 vt=relative_velocity-n*vn;
        bond.tangential_displacement+=vt*dt;
        bond.tangential_displacement-=n*dot(bond.tangential_displacement,n);

        const double extension=distance-bond.rest_length;
        const double normal_elastic=bond.model.normal_stiffness*extension;
        const Vec3 shear_elastic=bond.tangential_displacement*bond.model.shear_stiffness;
        const double tensile_ratio=std::max(0.0,normal_elastic)/bond.model.tensile_failure_force;
        const double shear_ratio=norm(shear_elastic)/bond.model.shear_failure_force;
        const double failure_ratio=std::max(tensile_ratio,shear_ratio);
        if(failure_ratio>bond.model.damage_onset_ratio){
            const double target=clamp_value((failure_ratio-bond.model.damage_onset_ratio)/(1.0-bond.model.damage_onset_ratio),0.0,1.0);
            bond.damage=std::max(bond.damage,target);
        }
        if(failure_ratio>=1.0 || bond.damage>=1.0){
            bond.damage=1.0;
            bond.broken=true;
            --stats.active;
            ++stats.broken_total;
            ++stats.broken_this_step;
            continue;
        }
        if(bond.damage>0.0) ++stats.damaged;
        const double intact=1.0-bond.damage;
        const double normal_force=intact*normal_elastic + bond.model.normal_damping*vn;
        const Vec3 shear_force=shear_elastic*intact + vt*bond.model.shear_damping;
        const Vec3 force_on_a=n*normal_force+shear_force;
        if(!a.fixed) a.add_force(force_on_a);
        if(!b.fixed) b.add_force(-force_on_a);
    }
    return stats;
}

void resolve_impulse_contacts(std::vector<RigidBody>& bodies,const std::vector<ContactPoint>& contacts,double dt,const ImpulseContactModel& m) {
    if(!(dt>0.0)) throw std::invalid_argument("contact dt must be positive");
    std::vector<double> accumulated(contacts.size(),0.0);
    std::vector<double> target_velocity(contacts.size(),0.0);
    for(std::size_t ci=0;ci<contacts.size();++ci){
        const auto& c=contacts[ci];
        const Vec3 rel=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);
        const double vn=dot(rel,c.normal);
        target_velocity[ci]=std::max((m.baumgarte/dt)*c.penetration,-m.restitution*std::min(vn,0.0));
    }
    for(std::size_t it=0;it<m.iterations;++it){
        for(std::size_t ci=0;ci<contacts.size();++ci){
            const auto& c=contacts[ci];
            const Vec3 rel=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);
            const double vn=dot(rel,c.normal);
            const double k=contact_effective_mass(bodies,c,c.normal);
            if(k<=1e-18)continue;
            const double old=accumulated[ci];
            const double now=std::max(0.0,old+(target_velocity[ci]-vn)/k);
            const double djn=now-old;accumulated[ci]=now;apply_contact_impulse(bodies,c,c.normal*djn);
            const Vec3 rel2=velocity_at_contact(bodies,c.body_a,c.point)-velocity_at_contact(bodies,c.body_b,c.point);
            const Vec3 vt=rel2-c.normal*dot(rel2,c.normal);const double vm=norm(vt);
            if(vm>1e-14){const Vec3 t=vt/vm;const double kt=contact_effective_mass(bodies,c,t);if(kt>1e-18){const double jt=clamp_value(-vm/kt,-m.friction*now,m.friction*now);apply_contact_impulse(bodies,c,t*jt);}}
        }
    }
}


RigidBodySystem::RigidBodySystem(MultibodyConfig c):config_(c){ if(!(c.dt>0.0))throw std::invalid_argument("multibody dt must be positive"); }
std::size_t RigidBodySystem::add_body(RigidBody b){b.validate();b.state.orientation=normalized(b.state.orientation);bodies_.push_back(b);return bodies_.size()-1;}
void RigidBodySystem::add_sphere(std::size_t body,double radius){if(body>=bodies_.size()||!(radius>0.0))throw std::invalid_argument("invalid sphere attachment");spheres_.push_back({body,radius});}
void RigidBodySystem::add_plane(PlaneShape p){p.normal=normalized(p.normal);planes_.push_back(p);}
void RigidBodySystem::add_triangle_mesh(TriangleMeshShape mesh){
    if(mesh.body!=invalid_body && mesh.body>=bodies_.size()) throw std::invalid_argument("invalid triangle mesh body attachment");
    if(mesh.feature_namespace==0U) mesh.feature_namespace=static_cast<std::uint32_t>(triangle_meshes_.size()+1U);
    triangle_meshes_.push_back(std::move(mesh));
}
void RigidBodySystem::add_joint(Joint j){joints_.push_back(std::move(j));}
void RigidBodySystem::step(){
    for(auto& b:bodies_) if(!b.fixed) b.add_force(config_.gravity*b.mass);
    auto pairs=broad_phase_sweep_and_prune(bodies_,spheres_);auto contacts=sphere_contacts(bodies_,spheres_,pairs);auto plane=sphere_plane_contacts(bodies_,spheres_,planes_);contacts.insert(contacts.end(),plane.begin(),plane.end());for(const auto& mesh:triangle_meshes_){auto mesh_contacts=sphere_triangle_mesh_contacts(bodies_,spheres_,mesh);contacts.insert(contacts.end(),mesh_contacts.begin(),mesh_contacts.end());}
    if(config_.contact_method==ContactMethod::smooth_penalty) apply_penalty_contacts(bodies_,contacts,config_.penalty);
    for(auto& b:bodies_) integrate_rigid_body(b,config_.dt,config_.integrator);
    std::vector<ConstraintRow> rows;for(const auto& j:joints_){auto r=build_constraint_rows(bodies_,j,config_.dt,config_.constraint_baumgarte);rows.insert(rows.end(),r.begin(),r.end());}solve_constraint_rows(bodies_,rows,config_.constraint_iterations);
    if(config_.contact_method==ContactMethod::nonsmooth_impulse) resolve_impulse_contacts(bodies_,contacts,config_.dt,config_.impulse);
    time_+=config_.dt;
}
void RigidBodySystem::run(std::size_t n){for(std::size_t i=0;i<n;++i)step();}
double RigidBodySystem::total_kinetic_energy() const noexcept {double e=0.0;for(const auto& b:bodies_){if(b.fixed)continue;e+=0.5*b.mass*norm_squared(b.state.linear_velocity);const Vec3 wl=inverse_rotate(b.state.orientation,b.state.angular_velocity);e+=0.5*(b.inertia_diagonal.x*wl.x*wl.x+b.inertia_diagonal.y*wl.y*wl.y+b.inertia_diagonal.z*wl.z*wl.z);}return e;}
Vec3 RigidBodySystem::linear_momentum() const noexcept {Vec3 p{};for(const auto& b:bodies_)if(!b.fixed)p+=b.state.linear_velocity*b.mass;return p;}

ExplicitDemSystem::ExplicitDemSystem(DemConfig c):config_(c){if(!(c.dt>0.0))throw std::invalid_argument("DEM dt must be positive");}
std::size_t ExplicitDemSystem::add_particle(Vec3 pos,double r,double density,Vec3 vel){if(!(r>0.0&&density>0.0))throw std::invalid_argument("DEM radius/density must be positive");constexpr double pi=3.14159265358979323846;const double m=density*(4.0/3.0)*pi*r*r*r;const double I=0.4*m*r*r;RigidBody b;b.state.position=pos;b.state.linear_velocity=vel;b.mass=m;b.inertia_diagonal={I,I,I};const auto id=bodies_.size();bodies_.push_back(b);spheres_.push_back({id,r});return id;}
void ExplicitDemSystem::add_plane(PlaneShape p){p.normal=normalized(p.normal);planes_.push_back(p);}
void ExplicitDemSystem::add_triangle_mesh(TriangleMeshShape mesh){
    if(mesh.body!=invalid_body && mesh.body>=bodies_.size()) throw std::invalid_argument("invalid DEM triangle mesh body attachment");
    if(mesh.feature_namespace==0U) mesh.feature_namespace=static_cast<std::uint32_t>(triangle_meshes_.size()+1U);
    triangle_meshes_.push_back(std::move(mesh));
}
std::size_t ExplicitDemSystem::add_bond(std::size_t body_a,std::size_t body_b,BondedParticleModel model,double rest_length){
    bonds_.push_back(make_particle_bond(bodies_,body_a,body_b,model,rest_length));
    return bonds_.size()-1U;
}
void ExplicitDemSystem::step(){for(auto& b:bodies_)b.add_force(config_.gravity*b.mass);last_bond_stats_=apply_particle_bonds(bodies_,bonds_,config_.dt);auto pairs=broad_phase_sweep_and_prune(bodies_,spheres_);auto contacts=sphere_contacts(bodies_,spheres_,pairs);auto plane=sphere_plane_contacts(bodies_,spheres_,planes_);contacts.insert(contacts.end(),plane.begin(),plane.end());for(const auto& mesh:triangle_meshes_){auto mesh_contacts=sphere_triangle_mesh_contacts(bodies_,spheres_,mesh);contacts.insert(contacts.end(),mesh_contacts.begin(),mesh_contacts.end());}if(config_.history_dependent_tangential)apply_history_dependent_mindlin_contacts(bodies_,contacts,spheres_,contact_history_,config_.dt,config_.contact);else{contact_history_.clear();apply_hertz_mindlin_contacts(bodies_,contacts,spheres_,config_.contact);}for(auto& b:bodies_)integrate_rigid_body(b,config_.dt,RigidBodyIntegrator::semi_implicit_euler);time_+=config_.dt;}
void ExplicitDemSystem::run(std::size_t n){for(std::size_t i=0;i<n;++i)step();}
Vec3 stokes_drag_force(const RigidBody& body,double radius,Vec3 fluid_velocity,double mu){if(!(radius>0.0&&mu>=0.0))throw std::invalid_argument("invalid Stokes drag parameters");constexpr double pi=3.14159265358979323846;return (fluid_velocity-body.state.linear_velocity)*(6.0*pi*mu*radius);}

} // namespace cfd::multibody
