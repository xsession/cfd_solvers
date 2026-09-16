# Contributing

## Rules

- Keep solver/kernel code C++20.
- A normal CPU build must not require a GPU SDK.
- Add a manufactured/analytic regression test for every numerical operator where practical.
- Do not introduce virtual dispatch, allocation or string processing inside hot loops without benchmark evidence.
- Do not paste or mechanically translate source from GPL or FluidX3D-restricted references. Use mathematical literature/public APIs and independently written tests.
- Record benchmark hardware, compiler, backend and precision when claiming performance improvements.
- Keep reference and optimized implementations comparable so performance work cannot silently change physics.

## Local verification

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```
