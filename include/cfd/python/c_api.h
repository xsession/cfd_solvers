#pragma once
#include <stddef.h>
#ifdef _WIN32
#define CFD_PY_API __declspec(dllexport)
#else
#define CFD_PY_API __attribute__((visibility("default")))
#endif
#ifdef __cplusplus
extern "C" {
#endif
CFD_PY_API const char* cfd_solvers_version(void);
CFD_PY_API int cfd_write_xdmf_tri3(const char* path,const double* xy,size_t node_count,const size_t* tri,size_t triangle_count,const double* scalar,const char* scalar_name,char* error,size_t error_capacity);
#ifdef __cplusplus
}
#endif
