#include "cfd/solvers/lbm/esoteric_pull.hpp"

namespace cfd::lbm {

template class EsotericPullSolver<D2Q9InPlaceDescriptor>;
template class EsotericPullSolver<D3Q19Descriptor>;
template class EsotericPullSolver<D3Q27Descriptor>;

} // namespace cfd::lbm
