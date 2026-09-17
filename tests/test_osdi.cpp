#include "cfd/circuit/osdi_loader.hpp"
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv){
    try{
        if(argc!=2)throw std::runtime_error("expected fixture path");
        cfd::circuit::OsdiLibrary library(argv[1]);
        if(library.version_major()!=0U||library.version_minor()!=4U||library.descriptor_count()!=2U||!library.raw_descriptors())
            throw std::runtime_error("OSDI metadata mismatch");
        cfd::circuit::OsdiLibrary moved(std::move(library));
        if(moved.descriptor_count()!=2U||!moved.raw_descriptors())throw std::runtime_error("OSDI move semantics");
        std::cout<<"OSDI loader seam passed\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"OSDI test failure: "<<e.what()<<'\n';return 1;}
}
