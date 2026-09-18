# cfd_solvers v0.17.0 - k-space acoustics and photoacoustics

v0.17.0 starts Phase 14 with an independent Fourier-pseudospectral acoustics backend inspired by the public numerical-method descriptions around k-Wave.

## Added

- reusable radix-2 complex FFT and 2-D FFT;
- periodic Fourier spectral derivatives;
- periodic fractional-Laplacian operator;
- first-order 2-D k-space acoustic pressure/velocity solver;
- heterogeneous sound speed and density;
- k-space temporal sinc correction;
- split-density directional PML baseline;
- nonlinear B/A equation-of-state term;
- spectral frequency power-law attenuation baseline;
- photoacoustic initial-pressure initialization;
- point, plane and delayed-array pressure sources;
- sensor arrays;
- delay-and-sum beamforming;
- acoustic intensity, radiation pressure and heating helpers.

## Validation

The focused regression checks spectral derivative accuracy, fractional-Laplacian eigenmodes, stronger attenuation of higher spatial frequencies, bounded lossless propagation energy error, heterogeneous nonlinear photoacoustic propagation, source/sensor arrays, beamforming and post-processing formulas.

The complete default project matrix passes **141/141 CTest targets** and Python ABI **0.17.0**. Phase 14 starts at **15/21 = 71.4%**; the expanded overall tracker is **664/777 = 85.5%**. The full result is recorded in the packaged `ctest_v0170_full.log`.

## Scope limits

- current solver is 2-D Cartesian and power-of-two in each FFT dimension;
- current power-law attenuation is a spectral attenuation baseline rather than the complete causal fractional absorption/dispersion constitutive model;
- time reversal and direct Pennes/structural/piezoelectric coupling remain open;
- FFT workspaces are CPU-resident and not yet SYCL/distributed.
