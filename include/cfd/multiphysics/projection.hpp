#pragma once
#include "cfd/fem/mesh2d.hpp"
#include "cfd/fvm/poly_mesh.hpp"
#include <cstddef>
#include <span>
#include <vector>
namespace cfd::multiphysics {
// Piecewise-linear Tri3 interpolation from FEM nodal values to FVM cell centers.
[[nodiscard]] std::vector<double> fem_nodes_to_fvm_cells(const cfd::fem::Mesh2D& source,
    std::span<const double> nodal_values,const cfd::fvm::PolyMesh& target);
// Inverse-distance reconstruction from FVM cell-centered values to FEM nodes.
// Exact for constants; `neighbours` controls the nearest-cell stencil size.
[[nodiscard]] std::vector<double> fvm_cells_to_fem_nodes(const cfd::fvm::PolyMesh& source,
    std::span<const double> cell_values,const cfd::fem::Mesh2D& target,std::size_t neighbours=8);
struct StructuredScalarGrid3D {
    std::size_t nx{},ny{},nz{}; double lx{1.0},ly{1.0},lz{1.0}; std::vector<double> values;
    void validate() const;
    [[nodiscard]] double sample(double x,double y,double z) const;
};
[[nodiscard]] std::vector<double> structured_to_fvm_cells(const StructuredScalarGrid3D& source,const cfd::fvm::PolyMesh& target);
[[nodiscard]] std::vector<double> structured_to_fem_nodes(const StructuredScalarGrid3D& source,const cfd::fem::Mesh2D& target,double z=0.0);
} // namespace cfd::multiphysics
