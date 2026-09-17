#include "cfd/core/backend.hpp"
#include <stdexcept>
namespace cfd::core {
const char* backend_name(RuntimeBackend b) noexcept {switch(b){case RuntimeBackend::serial:return "serial";case RuntimeBackend::openmp:return "openmp";case RuntimeBackend::mpi:return "mpi";case RuntimeBackend::sycl:return "sycl";}return "unknown";}
bool backend_available(RuntimeBackend b) noexcept {switch(b){case RuntimeBackend::serial:return true;case RuntimeBackend::openmp:
#if defined(CFD_HAS_OPENMP)
return true;
#else
return false;
#endif
case RuntimeBackend::mpi:
#if defined(CFD_HAS_MPI)
return true;
#else
return false;
#endif
case RuntimeBackend::sycl:
#if defined(CFD_HAS_SYCL)
return true;
#else
return false;
#endif
}return false;}
std::vector<BackendInfo> compiled_backends(){std::vector<BackendInfo> out;for(auto b:{RuntimeBackend::serial,RuntimeBackend::openmp,RuntimeBackend::mpi,RuntimeBackend::sycl})out.push_back({b,backend_name(b),backend_available(b)});return out;}
RuntimeBackend parse_backend(std::string_view s){if(s=="serial"||s=="cpu-serial")return RuntimeBackend::serial;if(s=="openmp"||s=="omp"||s=="cpu")return RuntimeBackend::openmp;if(s=="mpi")return RuntimeBackend::mpi;if(s=="sycl"||s=="gpu")return RuntimeBackend::sycl;throw std::invalid_argument("unknown runtime backend");}
RuntimeBackend select_backend(std::string_view request){if(request=="auto"){for(auto b:{RuntimeBackend::sycl,RuntimeBackend::mpi,RuntimeBackend::openmp,RuntimeBackend::serial})if(backend_available(b))return b;}
    const auto b=parse_backend(request);if(!backend_available(b))throw std::runtime_error("requested runtime backend was not compiled in");return b;}
} // namespace cfd::core
