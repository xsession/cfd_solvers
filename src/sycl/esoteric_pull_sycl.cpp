#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"

#if defined(CFD_HAS_SYCL)
namespace cfd::lbm {

template class EsotericPullSyclSolver<D2Q9InPlaceDescriptor>;
template class EsotericPullSyclSolver<D3Q19Descriptor>;
template class EsotericPullSyclSolver<D3Q27Descriptor>;

} // namespace cfd::lbm
#endif
