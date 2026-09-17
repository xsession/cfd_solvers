#include "cfd/io/mesh_io.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

#ifdef CFD_HAS_HDF5
#include <hdf5.h>
#endif

namespace cfd::io {
namespace {
void validate(const cfd::fem::Mesh2D& mesh,const std::vector<double>& scalar){
    mesh.validate();
    if(!scalar.empty()&&scalar.size()!=mesh.node_count())throw std::invalid_argument("XDMF scalar size mismatch");
}
void header(std::ostream& out){out<<"<?xml version=\"1.0\" ?>\n<Xdmf Version=\"3.0\">\n <Domain>\n  <Grid Name=\"mesh\" GridType=\"Uniform\">\n";}
void footer(std::ostream& out){out<<"  </Grid>\n </Domain>\n</Xdmf>\n";}
}

void write_xdmf_inline(const std::filesystem::path& path,const cfd::fem::Mesh2D& mesh,const std::vector<double>& scalar,std::string name){
    validate(mesh,scalar);std::ofstream out(path,std::ios::trunc);if(!out)throw std::runtime_error("cannot create XDMF");out<<std::setprecision(17);header(out);
    out<<"   <Topology TopologyType=\"Triangle\" NumberOfElements=\""<<mesh.element_count()<<"\">\n    <DataItem Dimensions=\""<<mesh.element_count()<<" 3\" NumberType=\"UInt\" Precision=\"8\" Format=\"XML\">\n";
    for(const auto& tri:mesh.triangles)out<<tri.node[0]<<' '<<tri.node[1]<<' '<<tri.node[2]<<'\n';
    out<<"    </DataItem>\n   </Topology>\n   <Geometry GeometryType=\"XY\">\n    <DataItem Dimensions=\""<<mesh.node_count()<<" 2\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n";
    for(const auto& node:mesh.nodes)out<<node.x<<' '<<node.y<<'\n';
    out<<"    </DataItem>\n   </Geometry>\n";
    if(!scalar.empty()){out<<"   <Attribute Name=\""<<name<<"\" AttributeType=\"Scalar\" Center=\"Node\">\n    <DataItem Dimensions=\""<<scalar.size()<<"\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\">\n";for(double value:scalar)out<<value<<' ';out<<"\n    </DataItem>\n   </Attribute>\n";}
    footer(out);
}

bool hdf5_output_available() noexcept{
#ifdef CFD_HAS_HDF5
    return true;
#else
    return false;
#endif
}

void write_hdf5_xdmf(const std::filesystem::path& xdmf_path,const std::filesystem::path& hdf5_path,const cfd::fem::Mesh2D& mesh,const std::vector<double>& scalar,std::string name){
    validate(mesh,scalar);
#ifndef CFD_HAS_HDF5
    static_cast<void>(xdmf_path);static_cast<void>(hdf5_path);static_cast<void>(name);
    throw std::runtime_error("HDF5 output requires a build with CFD_ENABLE_HDF5 and an HDF5 library");
#else
    const hid_t file=H5Fcreate(hdf5_path.string().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);if(file<0)throw std::runtime_error("cannot create HDF5 file");
    auto close_file=[&](){H5Fclose(file);};
    std::vector<double> points;points.reserve(mesh.node_count()*2U);for(const auto& n:mesh.nodes){points.push_back(n.x);points.push_back(n.y);}
    std::vector<unsigned long long> cells;cells.reserve(mesh.element_count()*3U);for(const auto& t:mesh.triangles)for(auto i:t.node)cells.push_back(static_cast<unsigned long long>(i));
    const hid_t mesh_group=H5Gcreate2(file,"/Mesh",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(mesh_group<0){close_file();throw std::runtime_error("cannot create HDF5 mesh group");}H5Gclose(mesh_group);
    const hid_t fields_group=H5Gcreate2(file,"/Fields",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(fields_group<0){close_file();throw std::runtime_error("cannot create HDF5 fields group");}H5Gclose(fields_group);
    auto dataset=[&](const std::string& path,hid_t type,const void* data,const std::vector<hsize_t>& dims){hid_t space=H5Screate_simple(static_cast<int>(dims.size()),dims.data(),nullptr);hid_t set=H5Dcreate2(file,path.c_str(),type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);if(set<0||H5Dwrite(set,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data)<0){if(set>=0)H5Dclose(set);H5Sclose(space);throw std::runtime_error("cannot write HDF5 dataset");}H5Dclose(set);H5Sclose(space);};
    try{dataset("/Mesh/Points",H5T_NATIVE_DOUBLE,points.data(),{mesh.node_count(),2});dataset("/Mesh/Cells",H5T_NATIVE_ULLONG,cells.data(),{mesh.element_count(),3});if(!scalar.empty())dataset("/Fields/"+name,H5T_NATIVE_DOUBLE,scalar.data(),{scalar.size()});}catch(...){close_file();throw;}close_file();
    std::ofstream out(xdmf_path,std::ios::trunc);if(!out)throw std::runtime_error("cannot create XDMF sidecar");const auto ref=hdf5_path.filename().generic_string();header(out);
    out<<"   <Topology TopologyType=\"Triangle\" NumberOfElements=\""<<mesh.element_count()<<"\"><DataItem Dimensions=\""<<mesh.element_count()<<" 3\" NumberType=\"UInt\" Precision=\"8\" Format=\"HDF\">"<<ref<<":/Mesh/Cells</DataItem></Topology>\n";
    out<<"   <Geometry GeometryType=\"XY\"><DataItem Dimensions=\""<<mesh.node_count()<<" 2\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">"<<ref<<":/Mesh/Points</DataItem></Geometry>\n";
    if(!scalar.empty())out<<"   <Attribute Name=\""<<name<<"\" AttributeType=\"Scalar\" Center=\"Node\"><DataItem Dimensions=\""<<scalar.size()<<"\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">"<<ref<<":/Fields/"<<name<<"</DataItem></Attribute>\n";
    footer(out);
#endif
}
} // namespace cfd::io
