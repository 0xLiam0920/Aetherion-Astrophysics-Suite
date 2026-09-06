# Contributing to Aetherion 

Welcome, and thank you for your interest in contributing. Below will explain how to do it properly, and what can/can’t be done.

## Code of Conduct

This project and everything about it is tied to the **CODE OF CONDUCT** (`CODE_OF_CONDUCT.md`). By contributing, you are expected to go by it. Please report anything that violates it in either a private message or a github issue.

## Ways to contribute

- **Report bugs** by opening an issue. Include your OS, compiler version, the
  exact build command you ran, and the full error output or a screenshot. The more details the better; I am one person after all.
- **Suggest features or physics improvements** by opening an issue describing
  the use case and justification. 
A note on this one: Eventually, I will add a plugin system that allows users to add extra parameters and other black holes via JSON format. If you suggest a black hole be added, I may use it as a demo for such a system later.
- When **submitting code**, do it via a pull request (see below).

## Building from source

The suite uses C++ version 20 and requires a CMAKE version greater than or equal to v3.10. It depends on Qt6, SFML 3x, GLM,
OpenGL, and (*NON APPLE ONLY*) GLEW. Dear ImGui and ImGui-SFML are
Already provided under `external/`. See the DEPENDENCIES file (`DEPENDENCIES.md`) for the full list and platform install tutorials.

For user convenience, I have written scripts to help with configuration/building-related steps:

```bash
./rebuild_macos.sh      # macOS (Homebrew toolchain)
./rebuild_linux.sh      # Linux
./rebuild_windows.ps1   # Windows (MSVC + vcpkg)
```

Alternatively, you may configure manually:

```bash
mkdir build && cd build
cmake ..
cmake --build . -j
```

## Running the tests

The physics regression suite is wired into CTest and is automatically run after every commit, PR, etc. If you wish to run it manually, then do the following after building:

```bash
cd build
ctest --output-on-failure
```

This runs the `physics-regression` C++ file, with the sole purpose of validating geodesic logic against known sources (e.g. the ~1.75″ solar-limb
deflection and angular-momentum conservation laws). Kindly ensure
it passes fine before opening a pull request. Add a new test case if you aim to change physics or fix a flaw in existing physics. Failure to do so will result in a closure to your PR and a comment regarding it.

## Pull request checklist

1. Fork this repository and make a new branch.
2. Keep changes to one task; other aspects can go in a different PR for convenience.
3. Try to keep the style and logic consistent; the code itself favors smaller, header-based structure in paths like `src/2D/2D-physics/` with all science and physics logic kept contained internally.
4. Build, give the relevant files (i.e, pkgs) as needed, and ensure all tests pass properly. Ignore JS warnings about deprecation.
5. In your actual Pull Request, give a changelog (added, changed, fixed, removed, etc) within the description or comments.
** FOR CHANGES IN PHYSICS AND MATHEMATICAL CONSTRUCTS: PLEASE CITE YOUR WORK AND/OR PROVIDE MATHEMATICAL BREAKDOWNS OF YOUR LOGIC **

## Essential project layout

* explained below are some important things to quickly know (summarized regarding the layout of this project: *

- `src/2D/2D-physics/`: Backend of scientific logic and structure: contains metric,
  geodesics, integrators, and units/their relevant conversions. Reused by both subprograms’ frontends
- `src/2D/2D-simulation/`:  2D subprogram and succeedinglogic
- `src/3D/`:  the 3D subprogram (has a copy of the Schwarzschild
  Logic found within the file `physics_overlay.hpp`).
- `src/MAIN-MENU/`, `src/QT-LAUNCHER/`:  UI program and main launcher.
- `tests/`: Testing files for physics, accretion, regression, disk color/temp, etc.

## Licensing

By contributing, you agree that your contributions will be licensed under the
[MIT License] (LICENSE) that covers the project. See `LICENSE.md` for more details.

