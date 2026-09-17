# v0.15.4 upstream documentation review - characteristic WENO compressible transport

This checkpoint is a clean-room implementation guided by equations and numerical-method descriptions in public documentation and literature. No upstream source code was copied.

## References reviewed

1. Jiang, G.-S. and Shu, C.-W., *Efficient Implementation of Weighted ENO Schemes*, Journal of Computational Physics 126 (1996), DOI 10.1006/jcph.1996.0130. NASA NTRS mirror: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960007052.pdf
2. Clawpack 5.10.x, one-dimensional Euler / Shu-Osher example: https://www.clawpack.org/gallery/pyclaw/gallery/shocksine.html
3. Clawpack SharpClaw solver documentation: https://www.clawpack.org/v5.10.x/pyclaw/solvers.html
4. SU2 governing-equations documentation for the compressible Euler system: https://su2code.github.io/docs_v7/Theory/

## Clean-room lessons applied

- The Euler working state is conservative: density, momentum and total energy. Pressure is recovered through a perfect-gas closure.
- High-order reconstruction must not be applied blindly component-by-component across strong coupled waves. The v0.15.4 path freezes a Roe-averaged Euler eigensystem at each face, projects the five-cell stencils into characteristic coordinates, performs scalar Jiang-Shu WENO5 reconstruction there, then transforms the two interface states back to conserved variables.
- The existing Rusanov numerical flux remains the robust Riemann-flux baseline. v0.15.4 changes reconstruction and time accuracy without pretending to introduce an exact or contact-resolving Riemann solver.
- Fifth-order spatial reconstruction requires a compatible method-of-lines time integrator for smooth transient tests. The new optional SSPRK3 path is used by the WENO validation while the legacy forward-Euler path remains available for backward compatibility.
- Near a shock, reconstructed states can be nonphysical even if neighboring cell averages are physical. The implementation therefore validates density and pressure and locally falls back to the adjacent cell average before evaluating a flux.
- A normalized pressure-jump sensor is exposed independently from reconstruction so the already implemented AMR layer can consume shock information later without coupling mesh policy into the Euler solver.

## Validation mapping

`tests/test_v0154_compressible_weno.cpp` validates three different failure modes:

- periodic smooth entropy-wave transport, including measured refinement order and conservation of mass, momentum and total energy;
- exact preservation of a uniform Euler state;
- Sod shock-tube positivity plus detection by the pressure-jump sensor.

The smooth-wave profile is initialized as an exact cell average of a sinusoidal entropy wave rather than point samples. That distinction is necessary when using a finite-volume reconstruction convergence test.

## Deliberate scope limits

This closes the tracker item for a characteristic high-order/WENO compressible scheme, but it is still a 1-D ideal-gas Euler baseline. It does not claim multidimensional characteristic reconstruction, positivity-preserving Zhang-Shu scaling, HLLC/Roe fluxes, viscous compressible Navier-Stokes, reacting shocks, or unstructured WENO reconstruction. Those remain independent extensions.
