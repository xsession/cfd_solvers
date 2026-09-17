#include "cfd/solvers/lbm/visualization.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace cfd::lbm {namespace {std::size_t id(std::size_t x,std::size_t y,std::size_t nx){return y*nx+x;}Vector2 sample(std::size_t nx,std::size_t ny,std::span<const Vector2>u,Vector2 p){if(p.x<0||p.y<0||p.x>static_cast<double>(nx-1)||p.y>static_cast<double>(ny-1))return {};auto x=static_cast<std::size_t>(std::clamp(p.x,0.0,static_cast<double>(nx-1))),y=static_cast<std::size_t>(std::clamp(p.y,0.0,static_cast<double>(ny-1)));return u[id(x,y,nx)];}}
void write_vtk_structured_2d(const std::filesystem::path&p,std::size_t nx,std::size_t ny,std::span<const double>s,std::span<const Vector2>u,const char*name){if(nx*ny!=s.size()||(!u.empty()&&u.size()!=s.size()))throw std::invalid_argument("LBM VTK field size mismatch");std::ofstream o(p);if(!o)throw std::runtime_error("cannot create LBM VTK");o<<"# vtk DataFile Version 3.0\ncfd_solvers LBM\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<nx<<' '<<ny<<" 1\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA "<<s.size()<<"\nSCALARS "<<name<<" double 1\nLOOKUP_TABLE default\n"<<std::setprecision(17);for(double v:s)o<<v<<'\n';if(!u.empty()){o<<"VECTORS velocity double\n";for(auto v:u)o<<v.x<<' '<<v.y<<" 0\n";}}
std::vector<double> vorticity_2d(std::size_t nx,std::size_t ny,std::span<const Vector2>u,double dx,double dy){if(nx<3||ny<3||u.size()!=nx*ny||!(dx>0)||!(dy>0))throw std::invalid_argument("invalid LBM vorticity grid");std::vector<double>w(nx*ny);for(std::size_t y=1;y+1<ny;++y)for(std::size_t x=1;x+1<nx;++x){double dvdx=(u[id(x+1,y,nx)].y-u[id(x-1,y,nx)].y)/(2*dx),dudy=(u[id(x,y+1,nx)].x-u[id(x,y-1,nx)].x)/(2*dy);w[id(x,y,nx)]=dvdx-dudy;}return w;}
std::vector<double> q_criterion_3d(std::size_t nx,std::size_t ny,std::size_t nz,std::span<const Vector3>u,double dx,double dy,double dz){
    if(nx<3||ny<3||nz<3||u.size()!=nx*ny*nz||!(dx>0)||!(dy>0)||!(dz>0))throw std::invalid_argument("invalid LBM Q-criterion grid");
    const auto at=[&](std::size_t x,std::size_t y,std::size_t z)->const Vector3&{return u[(z*ny+y)*nx+x];};
    std::vector<double>q(nx*ny*nz,0.0);
    for(std::size_t z=1;z+1<nz;++z)for(std::size_t y=1;y+1<ny;++y)for(std::size_t x=1;x+1<nx;++x){
        const Vector3 xp=at(x+1,y,z),xm=at(x-1,y,z),yp=at(x,y+1,z),ym=at(x,y-1,z),zp=at(x,y,z+1),zm=at(x,y,z-1);
        const double g[3][3]={{(xp.x-xm.x)/(2*dx),(yp.x-ym.x)/(2*dy),(zp.x-zm.x)/(2*dz)},
                              {(xp.y-xm.y)/(2*dx),(yp.y-ym.y)/(2*dy),(zp.y-zm.y)/(2*dz)},
                              {(xp.z-xm.z)/(2*dx),(yp.z-ym.z)/(2*dy),(zp.z-zm.z)/(2*dz)}};
        double s2=0.0,o2=0.0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){const double sij=0.5*(g[i][j]+g[j][i]),oij=0.5*(g[i][j]-g[j][i]);s2+=sij*sij;o2+=oij*oij;}
        q[(z*ny+y)*nx+x]=0.5*(o2-s2);
    }
    return q;
}
std::vector<double> horizontal_slice(std::size_t nx,std::size_t ny,std::span<const double>f,std::size_t y){if(f.size()!=nx*ny||y>=ny)throw std::invalid_argument("invalid LBM slice");return {f.begin()+static_cast<std::ptrdiff_t>(y*nx),f.begin()+static_cast<std::ptrdiff_t>((y+1)*nx)};}
std::vector<Vector2> streamline_2d(std::size_t nx,std::size_t ny,std::span<const Vector2>u,Vector2 p,double h,std::size_t max){if(u.size()!=nx*ny||!(h>0)||max==0)throw std::invalid_argument("invalid streamline controls");std::vector<Vector2>out;for(std::size_t i=0;i<max;++i){if(p.x<0||p.y<0||p.x>=static_cast<double>(nx)||p.y>=static_cast<double>(ny))break;out.push_back(p);auto v=sample(nx,ny,u,p);double m=std::hypot(v.x,v.y);if(m<1e-14)break;p.x+=h*v.x/m;p.y+=h*v.y/m;}return out;}
void write_html_scalar_viewer(const std::filesystem::path&p,std::size_t nx,std::size_t ny,std::span<const double>f){
    if(f.size()!=nx*ny||f.empty())throw std::invalid_argument("invalid HTML viewer field");
    std::ofstream o(p);if(!o)throw std::runtime_error("cannot create LBM HTML viewer");
    o<<"<!doctype html><meta charset=\"utf-8\"><title>cfd_solvers LBM viewer</title>\n"
       "<style>html,body{margin:0;height:100%;background:#111;color:#ddd;font:14px sans-serif;overflow:hidden}#hud{position:fixed;left:12px;top:10px;background:#000b;padding:8px;border-radius:6px}canvas{width:100%;height:100%;image-rendering:pixelated;cursor:grab}</style>\n"
       "<canvas id=\"view\"></canvas><div id=\"hud\">wheel: zoom · drag: pan<br><span id=\"value\"></span></div><script>\n";
    o<<"const nx="<<nx<<",ny="<<ny<<",field=["<<std::setprecision(17);for(std::size_t i=0;i<f.size();++i){if(i)o<<',';o<<f[i];}o<<"];\n";
    o<<R"JS(const c=document.getElementById('view'),g=c.getContext('2d'),v=document.getElementById('value');let zoom=1,ox=0,oy=0,drag=false,lx=0,ly=0;
const lo=Math.min(...field),hi=Math.max(...field),span=hi-lo||1,off=document.createElement('canvas');off.width=nx;off.height=ny;const og=off.getContext('2d'),im=og.createImageData(nx,ny);
for(let i=0;i<field.length;i++){const t=Math.max(0,Math.min(1,(field[i]-lo)/span)),r=Math.round(255*t),b=Math.round(255*(1-t)),q=4*i;im.data[q]=r;im.data[q+1]=Math.round(255*(1-Math.abs(2*t-1)));im.data[q+2]=b;im.data[q+3]=255;}og.putImageData(im,0,0);
function resize(){c.width=innerWidth*devicePixelRatio;c.height=innerHeight*devicePixelRatio;draw()}function draw(){g.setTransform(1,0,0,1,0,0);g.clearRect(0,0,c.width,c.height);const base=Math.min(c.width/nx,c.height/ny)*zoom,w=nx*base,h=ny*base,x=(c.width-w)/2+ox,y=(c.height-h)/2+oy;g.imageSmoothingEnabled=false;g.drawImage(off,x,y,w,h);c._map={x,y,w,h}}
c.addEventListener('wheel',e=>{e.preventDefault();zoom*=Math.exp(-e.deltaY*.001);zoom=Math.max(.1,Math.min(100,zoom));draw()},{passive:false});c.addEventListener('pointerdown',e=>{drag=true;lx=e.clientX*devicePixelRatio;ly=e.clientY*devicePixelRatio;c.setPointerCapture(e.pointerId)});c.addEventListener('pointermove',e=>{const x=e.clientX*devicePixelRatio,y=e.clientY*devicePixelRatio;if(drag){ox+=x-lx;oy+=y-ly;lx=x;ly=y;draw()}const m=c._map,ix=Math.floor((x-m.x)/m.w*nx),iy=Math.floor((y-m.y)/m.h*ny);v.textContent=(ix>=0&&iy>=0&&ix<nx&&iy<ny)?`(${ix}, ${iy}) = ${field[iy*nx+ix]}`:''});c.addEventListener('pointerup',()=>drag=false);addEventListener('resize',resize);resize();
</script>)JS";
}
void write_pgm_scalar(const std::filesystem::path&p,std::size_t nx,std::size_t ny,std::span<const double>f){if(f.size()!=nx*ny||f.empty())throw std::invalid_argument("invalid PGM field");auto [mn,mx]=std::minmax_element(f.begin(),f.end());double span=*mx-*mn;std::ofstream o(p);if(!o)throw std::runtime_error("cannot create PGM");o<<"P2\n"<<nx<<' '<<ny<<"\n255\n";for(double v:f){int q=span>0?static_cast<int>(std::lround(255.0*(v-*mn)/span)):0;o<<std::clamp(q,0,255)<<' ';}}
} // namespace cfd::lbm
