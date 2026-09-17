#include "cfd/circuit/osdi_loader.hpp"

#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace cfd::circuit {
namespace {

#if defined(_WIN32)
void* open_library(const std::string& path) {
    return reinterpret_cast<void*>(LoadLibraryA(path.c_str()));
}
void close_library(void* handle) noexcept {
    if (handle) FreeLibrary(reinterpret_cast<HMODULE>(handle));
}
void* load_symbol(void* handle,const char* name) {
    return reinterpret_cast<void*>(GetProcAddress(reinterpret_cast<HMODULE>(handle),name));
}
std::string loader_error() { return "Windows dynamic-loader error"; }
#else
void* open_library(const std::string& path) {
    dlerror();
    return dlopen(path.c_str(),RTLD_NOW|RTLD_LOCAL);
}
void close_library(void* handle) noexcept {
    if (handle) dlclose(handle);
}
void* load_symbol(void* handle,const char* name) {
    dlerror();
    return dlsym(handle,name);
}
std::string loader_error() {
    const char* error=dlerror();
    return error?std::string(error):std::string("dynamic-loader error");
}
#endif

template<class T>
const T* required_variable(void* handle,const char* name) {
    void* symbol=load_symbol(handle,name);
    if(!symbol) throw std::runtime_error(std::string("OSDI module missing symbol ")+name+": "+loader_error());
    return static_cast<const T*>(symbol);
}

} // namespace

OsdiLibrary::OsdiLibrary(const std::string& path):path_(path) {
    if(path.empty()) throw std::invalid_argument("OSDI module path must not be empty");
    handle_=open_library(path);
    if(!handle_) throw std::runtime_error("cannot load OSDI module "+path+": "+loader_error());
    try {
        version_major_=*required_variable<std::uint32_t>(handle_,"OSDI_VERSION_MAJOR");
        version_minor_=*required_variable<std::uint32_t>(handle_,"OSDI_VERSION_MINOR");
        descriptor_count_=*required_variable<std::uint32_t>(handle_,"OSDI_NUM_DESCRIPTORS");
        descriptors_=load_symbol(handle_,"OSDI_DESCRIPTORS");
        if(descriptor_count_>0U&&!descriptors_)
            throw std::runtime_error("OSDI module declares descriptors but exports no OSDI_DESCRIPTORS symbol");
    } catch (...) {
        close();
        throw;
    }
}

OsdiLibrary::~OsdiLibrary(){close();}

OsdiLibrary::OsdiLibrary(OsdiLibrary&& other) noexcept
    :handle_(std::exchange(other.handle_,nullptr)),
     version_major_(other.version_major_),version_minor_(other.version_minor_),
     descriptor_count_(other.descriptor_count_),descriptors_(other.descriptors_),
     path_(std::move(other.path_)) {
    other.descriptors_=nullptr;other.descriptor_count_=0U;
}

OsdiLibrary& OsdiLibrary::operator=(OsdiLibrary&& other) noexcept {
    if(this==&other)return *this;
    close();
    handle_=std::exchange(other.handle_,nullptr);
    version_major_=other.version_major_;version_minor_=other.version_minor_;
    descriptor_count_=other.descriptor_count_;descriptors_=other.descriptors_;
    path_=std::move(other.path_);
    other.descriptors_=nullptr;other.descriptor_count_=0U;
    return *this;
}

void OsdiLibrary::close() noexcept {
    if(handle_){close_library(handle_);handle_=nullptr;}
    descriptors_=nullptr;
}

} // namespace cfd::circuit
