# Release 0.10.1 - PIC particle sorting and guard-halo decomposition foundation

## Added

- Stable 3-D PIC particle sorting by cell with cell offsets, cell counts, original-index recovery and occupied-cell diagnostics.
- 3-D guard-halo particle classification for x/y/z minus/plus face exchange; edge and corner particles can appear in multiple lists.
- Optional sorted particle storage in `StaggeredElectromagneticPic3D` with configurable interval and diagnostics.
- `particle-sort-halo3d` CLI smoke case.
- `cfd-v0101-particle-sorting-tests` focused regression target.

## Validation

- Full debug/CI-speed CPU CTest matrix: 74/74 passed.
- Focused ASan + UBSan + leak detection: passed for the v0.10.1 particle sorting and guard-halo slice.

## Tracker

- Overall: 521/663 capabilities, 78.6%.
- Phase 12 CST-class expansion: 55/75 capabilities, 73.3%.

## Limitations

This release is not yet MPI particle exchange, guard-cell field exchange or GPU particle tiling. It is the validated serial data-layout foundation for those next roadmap steps.
