# Release 0.13.0 - Phase 9 interoperability completion

## Summary

This release completes the two remaining Phase 9 capabilities: HDF5/XDMF result output and Python bindings. Both are designed around a narrow, stable boundary rather than exposing internal C++ object layouts.

## HDF5/XDMF output

- `write_xdmf_inline(...)` writes Tri3 topology, XY geometry and an optional nodal scalar directly in portable XDMF XML.
- `write_hdf5_xdmf(...)` writes `/Mesh/Points`, `/Mesh/Cells` and `/Fields/<name>` datasets plus an XDMF sidecar when HDF5 is available.
- `hdf5_output_available()` lets callers inspect the compiled capability.
- Builds without HDF5 retain inline XDMF and fail explicitly if native HDF5 output is requested.
- `CFD_ENABLE_HDF5` controls discovery and integration.

## Python bindings

- `cfd_solvers_capi` provides a stable shared C ABI.
- `python/cfd_solvers` loads that ABI with Python's standard `ctypes` module.
- The first bound operation writes Tri3 meshes and nodal fields through the same validated C++ XDMF path.
- `CFD_SOLVERS_LIBRARY` can select the shared library at runtime.
- No pybind11 or CPython development headers are required.

## Validation

- `cfd-v0130-xdmf-python-tests` validates topology, geometry, scalar output, invalid input handling and the non-HDF5 fail-closed boundary.
- `cfd-v0130-python-binding-smoke` loads the produced shared library and writes a real XDMF file through Python.
- Native HDF5 dataset execution requires an HDF5-equipped CI/build environment; this runtime did not provide HDF5.

## Status

- Overall tracker: 580/719 = 80.7%.
- Phase 9: 50/50 = 100.0%.
