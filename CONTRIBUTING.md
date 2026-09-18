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

## Formatting hook

Install the repository Git hooks once per clone:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/install_git_hooks.ps1
```

or on POSIX shells:

```bash
sh scripts/install_git_hooks.sh
```

The pre-commit hook runs `clang-format` on staged C/C++ files and re-stages the formatted files. If a staged file also has unstaged edits, the hook stops so local work is not mixed into the commit.
