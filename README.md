# Aetherion: Astrophysics Simulation Suite

[![Version](https://img.shields.io/badge/version-0.2.1--beta2-blue)]
[![Status](https://img.shields.io/badge/status-beta-orange)]
[![License](https://img.shields.io/badge/license-MIT-green)]

Aetherion is an open source astrophysics simulation and analysis suite for experimenting with, visualising, and analysing BLACK HOLES.

DISCLAIMER!
                                                       
As of 17 April 2026, Aetherion is in active beta. The core simulation, visualization,
analysis, and export workflows are implemented and continue to undergo testing and
refinement. Features and interfaces may change between releases, hence the project
is not ready to be considered production-grade.

If you encounter a crash, incorrect result, or other reproducible issue, please open
a GitHub issue with the relevant platform, build configuration, inputs / scenario,
and steps required to reproduce the problem. This information is particularly useful
while the project’s regression coverage and platform support continue to expand.

**Physics Model Limitations:**
- **Temperature/Sound Speed Model:** The gas temperature calculations assume ideal gas behavior, which is optimistic for hot accreting plasma near black holes where radiation pressure and magnetic effects dominate. This is fine for Bondi radius order-of-magnitude estimates but should not be relied upon for precision plasma dynamics.
---

## About

Aetherion sits between a raw physics simulator and a data-analysis tool. You can run a scenario, watch what happens, and then export the data to look at it properly.

The program focuses on:
- Deterministic simulation runs (same inputs, same outputs)
- A modular codebase so individual pieces can be swapped out
- Clean data export for external or in-app analysis
- Built-in visualisation

The goal is reproducible results without having to handle inputs and test outputs manually. 

---

## Core Features

- **Simulation Engine**
  - Orbital mechanics
  - Photon trajectory modeling
  - Relativistic effects (e.g., deflection, precession)

- **Data Analysis Suite**
  - CSV and FITS export
  - Interactive graphing (orbit paths, energy conservation, etc.)
  - Multi-dataset comparison

- **Visualization**
  - Real-time graphing visualization
  - 2D and 3D rendering of orbits and photon paths
  - Multi-tab analysis interface
  - Clean UI for rapid inspection

- **Object Library**
  - Reusable custom bodies (i.e, pulsar)
  - Configurable parameters (mass, gas temp, etc, still WIP)

- **Architecture**
  - Modularity for easy extension (NOT YET IMPLEMENTED)
  - Plugin system for custom modularity (NOT YET IMPLEMENTED)
  - Configurable for keybinds, scenarios, etc

---

## Installation

### macOS (Homebrew)

```bash
brew install cmake sfml glm
git clone https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite.git
cd Aetherion-Astrophysics-Suite
mkdir build && cd build
cmake ..
make -j$(nproc)
./blackhole-sim
```

### Linux (Ubuntu / Debian)

```bash
sudo apt install cmake g++ \
    libsfml-dev \
    libglm-dev \
    libglew-dev \
    qt6-base-dev qt6-charts-dev \
    qt6-qpa-plugins           # XCB platform plugin, required at runtime
git clone https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite.git
cd Aetherion-Astrophysics-Suite
mkdir build && cd build
cmake ..
make -j$(nproc)
./blackhole-sim
```

> **Note:** SFML 3.x is required. As of May 2026, Ubuntu/Debian's `libsfml-dev`
> still ships **SFML 2.6** on most stable releases (24.04 LTS included), this
> will not work; the build will fail at `cmake` with a clear error. Either build
> SFML 3.0 from source ([SFML 3.0.0 release](https://github.com/SFML/SFML/releases/tag/3.0.0))
> or use the Flatpak build below. Ubuntu 25.04+ / Debian trixie+ ship SFML 3 and
> work out of the box.

### Linux (Fedora / RHEL)

```bash
sudo dnf install cmake gcc-c++ \
    SFML-devel \
    glm-devel \
    glew-devel \
    qt6-qtbase-devel qt6-qtcharts-devel
```

### Linux (Flatpak, recommended for end users)

```bash
flatpak install flathub org.kde.Platform//6.8 org.kde.Sdk//6.8
flatpak-builder --user --install --force-clean build-flatpak \
    flatpak/io.github.0xLiam0920.AetherionSuite.json
flatpak run io.github.0xLiam0920.AetherionSuite
```

### Windows (MSVC + vcpkg)

Prerequisites (all one-time):

1. **Git**, **CMake ≥ 3.21**, and **Visual Studio 2022** with the
   *Desktop development with C++* workload. From an elevated PowerShell:

   ```powershell
   winget install --id Git.Git -e
   winget install --id Kitware.CMake -e
   # Full IDE:
   winget install --id Microsoft.VisualStudio.2022.Community -e `
     --override "--quiet --wait --add Microsoft.VisualStudio.Workload.NativeDesktop --includeRecommended"
   # …or Build Tools only (no IDE, ~6 GB smaller):
   winget install --id Microsoft.VisualStudio.2022.BuildTools -e `
     --override "--quiet --wait --add Microsoft.VisualStudio.Workload.VCTools --includeRecommended"
   ```

2. **vcpkg**, bootstrapped and exported on the environment:

   ```powershell
   git clone https://github.com/microsoft/vcpkg C:\vcpkg
   C:\vcpkg\bootstrap-vcpkg.bat
   [Environment]::SetEnvironmentVariable('VCPKG_ROOT', 'C:\vcpkg', 'User')
   $env:VCPKG_ROOT = 'C:\vcpkg'
   ```

3. **Dependencies** (x64-windows triplet):

   ```powershell
   C:\vcpkg\vcpkg.exe install sfml glm glew qtbase qtcharts --triplet x64-windows
   ```

Then clone and build:

```powershell
git clone https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite.git
cd Aetherion-Astrophysics-Suite
.\rebuild_windows.ps1 -Config Release
.\build\Release\blackhole-sim.exe
```

The build script auto-detects `VCPKG_ROOT`, invokes the Visual Studio 2022
generator, and deploys the required Qt/SFML runtime DLLs plus the
`platforms\qwindows.dll` plugin next to each executable. Output binaries:

- `build\Release\blackhole-sim.exe` — main Qt launcher
- `build\Release\blackhole-2D.exe` — standalone 2D sim
- `build\Release\blackhole-3D.exe` — standalone 3D sim

> **Note:** `rebuild_windows.ps1` must run from a shell where the MSVC toolchain
> is on `PATH`. Either launch *Developer PowerShell for VS 2022* from the Start
> menu, or dot-source `Launch-VsDevShell.ps1 -Arch amd64` from a regular
> PowerShell first.

### Physics regression tests

```bash
cd build
ctest --output-on-failure -R physics-regression
```

---

## Known Issues

- SFML 2.x is not supported, 3.x API is required.
- On Wayland, the app forces X11/XWayland via `QT_QPA_PLATFORM=xcb` (required by SFML). Set `QT_QPA_PLATFORM=xcb` manually if needed.
- None of the packages are certified/signed for their relevant platforms; your AV may flag them as potentially dangerous. Just whitelist the app (or binaries if building from source) as needed
- Various scenarios and objects are still in progress; 2d subprogram logic is largely complete, but the 3d ones are about 40% through as of `v0.2.1`
- Windows: MSVC + vcpkg via `rebuild_windows.ps1` is the supported path (see *Installation → Windows*). Packaging an `.exe` installer is handled by `make_exe.sh`.
- If you build from source, the background file (`background.png`) may not be copied to the build directory. If you see a white background, either the copy is corrupted or missing, so just copy it over. Will be fixed as live backgrounds based on star data are added ("BETA 3").
---

## License

See [LICENSE](LICENSE) for full terms.

