#include "cfd/solvers/fvm/particles.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace cfd::fvm {namespace {
std::size_t nearest(const PolyMesh&m,Vec3 p){std::size_t best=0;double d2=1e300;for(std::size_t c=0;c<m.cell_count();++c){auto d=m.cells()[c].center-p;double v=dot(d,d);if(v<d2){d2=v;best=c;}}return best;}
Vec3 cross(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double volume(const Particle&p){return std::numbers::pi*p.diameter*p.diameter*p.diameter/6.0;}
}
void ParticleCloud::add(Particle p){if(!(p.diameter>0)||!(p.density>0))throw std::invalid_argument("invalid particle");particles_.push_back(p);}
void ParticleCloud::advance(const PolyMesh&m,std::span<const Vec3>u,double mu,double dt,bool two){if(u.size()!=m.cell_count()||!(mu>0)||!(dt>0))throw std::invalid_argument("invalid particle advance");momentum_source_.assign(m.cell_count(),{});for(auto&p:particles_){auto c=nearest(m,p.position);Vec3 old=p.velocity;double tau=p.density*p.diameter*p.diameter/(18.0*mu),decay=std::exp(-dt/tau);p.velocity=u[c]+(p.velocity-u[c])*decay;p.position+=p.velocity*dt;if(two){double mass=p.density*volume(p);momentum_source_[c]-=(p.velocity-old)*(mass/(dt*m.cells()[c].volume));}}}
Vec3 saffman_lift_force(Vec3 rel,Vec3 vort,double rho,double nu,double d,double c){if(!(rho>0)||!(nu>0)||!(d>0)||!(c>=0))throw std::invalid_argument("invalid Saffman lift state");const double om=magnitude(vort);if(om==0)return {};const double scale=c*rho*d*d*std::sqrt(nu*om);return cross(rel,vort)*(scale/om);}
Vec3 virtual_mass_force(Vec3 af,Vec3 ap,double rho,double d,double c){if(!(rho>0)||!(d>0)||!(c>=0))throw std::invalid_argument("invalid virtual-mass state");const double displaced=rho*std::numbers::pi*d*d*d/6.0;return (af-ap)*(c*displaced);}
void collide_particles_elastic(Particle&a,Particle&b,double e){if(!(e>=0&&e<=1))throw std::invalid_argument("invalid restitution");Vec3 n=b.position-a.position;double dist=magnitude(n);if(dist==0)return;n=n/dist;double rel=dot(b.velocity-a.velocity,n);if(rel>=0)return;double ma=a.density*volume(a),mb=b.density*volume(b);double j=-(1+e)*rel/(1/ma+1/mb);a.velocity-=n*(j/ma);b.velocity+=n*(j/mb);}
Particle coalesce_particles(const Particle&a,const Particle&b){double va=volume(a),vb=volume(b),ma=a.density*va,mb=b.density*vb,m=ma+mb,v=va+vb;Particle p;p.diameter=std::cbrt(6*v/std::numbers::pi);p.density=m/v;p.position=(a.position*ma+b.position*mb)/m;p.velocity=(a.velocity*ma+b.velocity*mb)/m;return p;}
bool breakup_by_weber(double rho,double u,double d,double sigma,double critical){if(!(rho>0)||!(d>0)||!(sigma>0)||!(critical>0))throw std::invalid_argument("invalid breakup state");return rho*u*u*d/sigma>critical;}
double evaporate_d2_law(double d,double k,double dt){if(!(d>=0)||!(k>=0)||!(dt>=0))throw std::invalid_argument("invalid evaporation state");return std::sqrt(std::max(0.0,d*d-k*dt));}
std::vector<double> particle_volume_fraction(const PolyMesh&m,std::span<const Particle>p){std::vector<double>a(m.cell_count());for(const auto&x:p){auto c=nearest(m,x.position);a[c]+=volume(x)/m.cells()[c].volume;}for(double&v:a)v=std::clamp(v,0.0,1.0);return a;}
}
