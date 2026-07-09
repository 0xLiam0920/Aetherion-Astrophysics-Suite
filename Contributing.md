# Contributing to the Aetherion Astrophysics Suite

Thanks for your interest in contributing! This document explains how to build
the project, run the tests, and submit changes.

## Code of Conduct

This project and everyone participating in it is governed by the
[Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to
uphold it. Please report unacceptable behaviour to the maintainer.

## Ways to contribute

- **Report bugs** by opening an issue. Include your OS, compiler version, the
  exact build command you ran, and the full error output or a screenshot.
- **Suggest features or physics improvements** by opening an issue describing
  the use case.
- **Submit code** via a pull request (see below).

## Building from source

The suite is C++17 and uses CMake (>= 3.10). It depends on Qt6, SFML 3.x, GLM,
OpenGL, and (on non-Apple platforms) GLEW. Dear ImGui and ImGui-SFML are
vendored under `external/`. See [DEPENDENCIES.md](DEPENDENCIES.md) for the full
list and per-platform install instructions.

Convenience scripts wrap the configure/build steps:

```bash
./rebuild_macos.sh      # macOS (Homebrew toolchain)
./rebuild_linux.sh      # Linux
./rebuild_windows.ps1   # Windows (MSVC + vcpkg)
```

Or configure manually:

```bash
mkdir build && cd build
cmake ..
cmake --build . -j
```

## Running the tests

The physics regression suite is wired into CTest. After building:

```bash
cd build
ctest --output-on-failure
```

This runs the `physics-regression` target, which validates the geodesic
integrators against known analytic results (e.g. the ~1.75″ solar-limb
deflection and long-run energy/angular-momentum conservation). Please make sure
the suite passes before opening a pull request, and add a new test case when you
change physics code or fix a physics bug.

## Pull request checklist

1. Fork the repository and create a topic branch off `main`.
2. Keep changes focused; unrelated fixes belong in separate PRs.
3. Match the surrounding code style (the physics core favours small, header-only
   structs in `src/2D/2D-physics/` with the scientific API kept
   rendering-independent).
4. Build cleanly and ensure `ctest` passes.
5. Describe **what** changed and **why** in the PR description. For physics
   changes, cite the reference or derivation you used.

## Project layout

- `src/2D/2D-physics/` — rendering-independent scientific core (metric,
  geodesics, integrators, units). Reused by both the 2D and 3D front-ends.
- `src/2D/2D-simulation/` — simulation orchestration and CSV/FITS export.
- `src/3D/` — the 3D viewer (carries a self-contained copy of the Schwarzschild
  helpers in `physics_overlay.hpp`).
- `src/MAIN-MENU/`, `src/QT-LAUNCHER/` — UI host and launcher.
- `tests/` — physics regression tests.

## Licensing

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE) that covers the project.
