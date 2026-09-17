#include "cfd/particle/electromagnetic.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);} 

double max_abs_on_boundary(const std::vector<double>& values,std::size_t nx,std::size_t ny){
    double peak=0.0;
    for(std::size_t iy=0;iy<ny;++iy){peak=std::max(peak,std::abs(values[iy*nx]));peak=std::max(peak,std::abs(values[iy*nx+nx-1U]));}
    for(std::size_t ix=0;ix<nx;++ix){peak=std::max(peak,std::abs(values[ix]));peak=std::max(peak,std::abs(values[(ny-1U)*nx+ix]));}
    return peak;
}

double max_abs(const std::vector<double>& values){double out=0.0;for(double value:values)out=std::max(out,std::abs(value));return out;}

void staggered_vacuum_tm_mode_is_bounded(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig config;config.nx=32U;config.ny=20U;config.length_x_m=1.0;config.length_y_m=0.75;
    const double dx=config.length_x_m/static_cast<double>(config.nx);const double dy=config.length_y_m/static_cast<double>(config.ny);
    config.dt_s=0.045/(299792458.0*std::sqrt(1.0/(dx*dx)+1.0/(dy*dy)));
    StaggeredElectromagneticPic2D solver(config);solver.initialize_tm_z_mode(0.8,1U,1U);
    const double initial=solver.diagnostics().field_energy_j;require(initial>0.0,"staggered Yee PIC initial field energy");
    solver.step(16U);const auto diag=solver.diagnostics();
    require(diag.particle_count==0U,"staggered Yee vacuum remains particle-free");
    require(std::isfinite(diag.field_energy_j)&&diag.field_energy_j>0.60*initial&&diag.field_energy_j<1.45*initial,"staggered Yee vacuum TM mode bounded energy");
}

void electric_wall_boundary_zeroes_tangential_fields(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig config;config.nx=12U;config.ny=10U;config.length_x_m=1.0;config.length_y_m=0.8;config.field_boundary.mode=GridBoundaryMode2D::electric_wall;config.periodic_particles=false;config.dt_s=1.0e-12;
    StaggeredElectromagneticPic2D solver(config);const std::size_t n=config.nx*config.ny;std::vector<double> one(n,1.0),two(n,2.0),three(n,3.0),bx(n,0.1),by(n,0.2),bz(n,0.3);
    solver.set_fields(one,two,three,bx,by,bz);solver.step(1U);
    double peak_ez=max_abs_on_boundary(solver.electric_z(),config.nx,config.ny);
    double peak_ex_y=0.0,peak_ey_x=0.0;
    for(std::size_t ix=0;ix<config.nx;++ix){peak_ex_y=std::max(peak_ex_y,std::abs(solver.electric_x()[ix]));peak_ex_y=std::max(peak_ex_y,std::abs(solver.electric_x()[(config.ny-1U)*config.nx+ix]));}
    for(std::size_t iy=0;iy<config.ny;++iy){peak_ey_x=std::max(peak_ey_x,std::abs(solver.electric_y()[iy*config.nx]));peak_ey_x=std::max(peak_ey_x,std::abs(solver.electric_y()[iy*config.nx+config.nx-1U]));}
    require(peak_ez==0.0&&peak_ex_y==0.0&&peak_ey_x==0.0,"staggered electric wall zeroes tangential E");
}

void absorbing_sponge_reduces_edge_energy(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig config;config.nx=18U;config.ny=14U;config.length_x_m=1.0;config.length_y_m=0.7;config.field_boundary.mode=GridBoundaryMode2D::absorbing_sponge;config.field_boundary.sponge_cells=4U;config.field_boundary.sponge_strength=2.5;config.periodic_particles=false;config.dt_s=1.0e-13;
    StaggeredElectromagneticPic2D solver(config);const std::size_t n=config.nx*config.ny;std::vector<double> zero(n,0.0),ez(n,0.0);
    for(std::size_t iy=0;iy<config.ny;++iy)for(std::size_t ix=0;ix<config.nx;++ix)if(ix<3U||iy<3U||ix+3U>=config.nx||iy+3U>=config.ny)ez[iy*config.nx+ix]=1.0;
    solver.set_fields(zero,zero,ez,zero,zero,zero);const double initial=solver.diagnostics().field_energy_j;solver.step(1U);const double after=solver.diagnostics().field_energy_j;
    require(after<initial,"staggered absorbing sponge damps edge-localized fields");
}

void particle_walls_absorb_and_reflect(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig absorb;absorb.nx=10U;absorb.ny=8U;absorb.length_x_m=1.0;absorb.length_y_m=0.8;absorb.dt_s=2.0e-11;absorb.field_boundary.mode=GridBoundaryMode2D::electric_wall;absorb.periodic_particles=false;absorb.particle_boundary=ParticleWallMode::absorb;absorb.enable_secondary_yield_report=true;
    StaggeredElectromagneticPic2D absorbing(absorb);PicParticle2D p;p.mass_kg=1.0;p.charge_c=1.0e-12;p.weight=2.0;p.position_m={0.999,0.4};p.velocity_m_per_s={1.0e8,0.0,0.0};absorbing.set_particles({p});absorbing.step(1U);auto diag=absorbing.diagnostics();
    require(diag.particle_count==0U&&diag.absorbed_particles==1U,"staggered absorbing particle wall removes escaping particle");

    StaggeredElectromagneticPic2DConfig reflect=absorb;reflect.particle_boundary=ParticleWallMode::specular_reflect;reflect.enable_secondary_yield_report=false;StaggeredElectromagneticPic2D reflecting(reflect);reflecting.set_particles({p});reflecting.step(1U);diag=reflecting.diagnostics();
    require(diag.particle_count==1U&&reflecting.particles().front().velocity_m_per_s.x<0.0,"staggered reflecting wall mirrors normal particle velocity");
}

void transverse_current_and_collision_couple(){
    using namespace cfd::particle;
    StaggeredElectromagneticPic2DConfig config;config.nx=18U;config.ny=14U;config.length_x_m=1.0;config.length_y_m=0.7;config.dt_s=1.0e-12;
    StaggeredElectromagneticPic2D solver(config);PicParticle2D p;p.mass_kg=9.1093837139e-31;p.charge_c=-1.602176634e-19;p.position_m={0.35,0.41};p.velocity_m_per_s={2.0e5,-1.0e5,3.0e5};solver.set_particles({p});NeutralCollisionModel gas;gas.neutral_density_per_m3=1.0e31;gas.elastic_cross_section_m2=1.0e-18;gas.random_seed=42U;solver.step(2U,&gas);const auto diag=solver.diagnostics();
    require(max_abs(solver.current_z())>0.0,"staggered Yee PIC deposits transverse current");
    require(max_abs(solver.electric_z())>0.0,"staggered Yee PIC transverse current drives Ez");
    require(diag.collision_statistics.elastic_events>0U,"staggered Yee PIC invokes MCC coupling");
}
}

int main(){
    try{
        staggered_vacuum_tm_mode_is_bounded();
        electric_wall_boundary_zeroes_tangential_fields();
        absorbing_sponge_reduces_edge_energy();
        particle_walls_absorb_and_reflect();
        transverse_current_and_collision_couple();
        std::cout<<"v0.9.6 staggered 2-D Yee EM-PIC tests passed\n";
        return 0;
    }catch(const std::exception& error){
        std::cerr<<"v0.9.6 staggered PIC2D test failure: "<<error.what()<<'\n';
        return 1;
    }
}
