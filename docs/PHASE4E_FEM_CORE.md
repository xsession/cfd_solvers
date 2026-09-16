# Phase 4E - General FEM Core

This checkpoint expands the Elmer-class branch from isolated proof solvers into a reusable element/assembly layer.

Implemented and validated:

- reference topology for Line2, Tri3, Quad4, Tet4, Hex8, Prism6 and Pyramid5;
- partition-of-unity shape functions and reference gradients;
- first/second-order Gaussian quadrature rules appropriate to each reference cell;
- isoparametric point mapping, Jacobians and physical shape gradients;
- Tri3 global CSR Laplace assembly plus a matrix-free element action with CSR parity test;
- patch-based scalar Dirichlet, Neumann and Robin conditions;
- variable-coefficient/reaction scalar diffusion on Tri3 meshes;
- 3-D Tet4 Poisson with manufactured sine solution and mesh-refinement convergence;
- plane-stress/plane-strain and axisymmetric small-strain elasticity patch tests;
- transient heat with a consistent P1 mass matrix;
- electrostatic potential/electric-field and DC conduction/current-density wrappers using the same scalar assembly.

Current limitations include no high-order elements, no mixed saddle-point element spaces, no nonlinear material laws, no adaptivity and no distributed FEM assembly yet.
