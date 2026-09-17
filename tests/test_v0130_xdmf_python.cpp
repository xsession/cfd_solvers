#include "cfd/fem/mesh2d.hpp"
#include "cfd/io/mesh_io.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}std::string read(const std::filesystem::path& p){std::ifstream in(p);return {std::istreambuf_iterator<char>(in),{}};}}
int main(){
    const auto root=std::filesystem::temp_directory_path()/("cfd_v0130_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));std::filesystem::create_directories(root);
    const auto mesh=cfd::fem::make_rectangle_tri_mesh(2,2);std::vector<double> scalar(mesh.node_count());for(std::size_t i=0;i<scalar.size();++i)scalar[i]=static_cast<double>(i);
    const auto path=root/"field.xdmf";cfd::io::write_xdmf_inline(path,mesh,scalar,"temperature");const auto xml=read(path);
    require(xml.find("TopologyType=\"Triangle\"")!=std::string::npos,"triangle topology missing");require(xml.find("GeometryType=\"XY\"")!=std::string::npos,"XY geometry missing");require(xml.find("Name=\"temperature\"")!=std::string::npos,"nodal attribute missing");
    bool rejected=false;try{cfd::io::write_xdmf_inline(root/"bad.xdmf",mesh,{1.0},"bad");}catch(const std::invalid_argument&){rejected=true;}require(rejected,"scalar mismatch should be rejected");
    if(!cfd::io::hdf5_output_available()){bool unavailable=false;try{cfd::io::write_hdf5_xdmf(root/"field-hdf.xdmf",root/"field.h5",mesh,scalar,"temperature");}catch(const std::runtime_error&){unavailable=true;}require(unavailable,"non-HDF5 build should fail explicitly");}
    std::filesystem::remove_all(root);std::cout<<"v0.13.0 XDMF/HDF5 boundary tests passed\n";return 0;
}
