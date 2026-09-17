#include "cfd/python/c_api.h"
#include "cfd/io/mesh_io.hpp"
#include <algorithm>
#include <cstring>
#include <exception>
#include <vector>

namespace {void copy_error(char* out,size_t cap,const char* text){if(!out||cap==0U)return;const auto n=std::min(cap-1U,std::strlen(text));std::memcpy(out,text,n);out[n]='\0';}}
extern "C" const char* cfd_solvers_version(void){return "0.16.3";}
extern "C" int cfd_write_xdmf_tri3(const char* path,const double* xy,size_t nodes,const size_t* tri,size_t triangles,const double* scalar,const char* name,char* error,size_t cap){
    try{
        if(!path||!xy||!tri||nodes==0U||triangles==0U)throw std::invalid_argument("path, coordinates and triangles are required");
        cfd::fem::Mesh2D mesh;mesh.nodes.reserve(nodes);mesh.boundary_node.assign(nodes,0U);for(size_t i=0;i<nodes;++i)mesh.nodes.push_back({xy[2U*i],xy[2U*i+1U]});mesh.triangles.reserve(triangles);for(size_t i=0;i<triangles;++i)mesh.triangles.push_back({{tri[3U*i],tri[3U*i+1U],tri[3U*i+2U]}});
        std::vector<double> values;if(scalar)values.assign(scalar,scalar+nodes);cfd::io::write_xdmf_inline(path,mesh,values,name?name:"field");return 0;
    }catch(const std::exception& exc){copy_error(error,cap,exc.what());return 1;}catch(...){copy_error(error,cap,"unknown C++ exception");return 2;}
}
