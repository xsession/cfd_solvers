# Phase 14 acoustics clean-room research

The Phase-14 implementation uses public k-Wave documentation and published method descriptions as capability and validation references only. No k-Wave source code is transplanted.

The public k-Wave documentation describes a first-order time-domain acoustic system with Fourier-collocation spatial derivatives, a k-space-corrected finite-difference time update, heterogeneous sound speed and density, nonlinear propagation, frequency power-law absorption, and a split-field PML. It also documents photoacoustic initial-value propagation, point/transducer sources, sensor recording, and related ultrasound workflows.

Primary references:

- https://www.k-wave.org/documentation.php
- https://k-wave.org/documentation/kspaceFirstOrder2D.php
- https://www.k-wave.org/documentation/example_na_modelling_absorption.php
- Treeby & Cox, *k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields*, J. Biomed. Opt. 15(2), 2010.
- Treeby et al., *Modeling nonlinear ultrasound propagation in heterogeneous media with power law absorption using a k-space pseudospectral method*, JASA 131(6), 2012.
- Treeby & Cox, *Modeling power law absorption and dispersion for acoustic propagation using the fractional Laplacian*, JASA 127(5), 2010.

## v0.17.0 clean-room boundaries

Implemented independently in v0.17.0:

- radix-2 1-D/2-D complex FFT;
- Fourier spectral derivatives and fractional Laplacian;
- 2-D first-order pressure/particle-velocity acoustic solver;
- heterogeneous density/sound-speed fields;
- k-space sinc correction;
- split density components plus directional PML damping;
- nonlinear B/A pressure-density term;
- frequency power-law spectral attenuation baseline;
- initial-pressure/photoacoustic propagation;
- point/plane/delayed-array pressure sources;
- pressure sensors and delay-and-sum beamforming;
- intensity, radiation-pressure and heating helpers.

The current power-law attenuation path is intentionally narrower than the complete causal k-Wave absorption/dispersion model. It damps each spatial Fourier mode according to a user-defined frequency power law; it does not yet reproduce the full Kramers-Kronig-consistent fractional constitutive relation in heterogeneous media. Likewise, the PML is a readable split-density/directional-damping baseline, not a claim of coefficient-level parity with k-Wave.

## Next validation/implementation steps

1. time-reversal reconstruction;
2. causal fractional absorption and associated dispersion;
3. 3-D and axisymmetric k-space solvers;
4. ultrasound transducer aperture/directivity models;
5. acoustic heating -> existing Pennes solver coupling;
6. structural and piezoelectric acoustic coupling;
7. inverse/optimization workflows;
8. resident SYCL FFT/k-space execution and distributed FFT support.

## v0.17.1 completion notes

The time-reversal workflow follows the documented k-Wave formulation conceptually: recorded pressure is enforced on the sensor mask in reverse time order as a Dirichlet boundary condition. The implementation is independent C++ and does not copy k-Wave source.

The multiphysics closure deliberately uses narrow reusable interfaces rather than merging solver internals: pressure traction for structures, volumetric heating for Pennes, radiation-force acceleration for FVM, a small coupled piezoelectric FEM baseline, and trace-misfit optimization through the common workflow layer.
