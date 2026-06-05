# Aetherion Astrophysics Suite — `v0.2.0-beta2`

> **The 2D Sub-App Maturity Release**
>
> *Release date:* June 2026
> *Codename:* **Equatorial**
> *Prior release:* [`BETA` (v0.1.0)](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/releases/tag/BETA), April 21 2026
> *Window covered:* 70 days · ~30 k lines added · ~18 k removed · 110+ files touched

[![Version](https://img.shields.io/badge/version-0.2.0--beta2-blue)]()
[![Status](https://img.shields.io/badge/status-beta-orange)]()
[![Platforms](https://img.shields.io/badge/platforms-macOS%20%7C%20Linux%20%7C%20Windows-lightgrey)]()
[![License](https://img.shields.io/badge/license-MIT-green)]()

---

## TL;DR

Beta 2 brings the 2D sub-application to ~95 % completion. The headline work is a
unified Qt + SFML launcher, full preset library of real-world black holes, custom
black-hole and custom-body authoring with on-disk persistence, BH–BH / BH–NS /
BH–Pulsar / BH–Star / BH–WD merger events with kind-specific dynamics, a Kerr
equatorial geodesic overlay, a thermal disk photon emitter, a sticky high-resolution
lensing mode, dark/light HUD theming, an animated parallax starfield, an
auto-updater, and full Windows + Linux/Flatpak parity alongside macOS.

The 3D sub-app received infrastructure work (modular `bh3d_core`, `physics_overlay`,
`imgui_qt_adapter`, photoreal disk shader) but its feature-completion target is
**beta 3**.

---

## Table of Contents

1. [Highlights](#highlights)
2. [What's new — 2D](#whats-new--2d)
3. [What's new — 3D](#whats-new--3d)
4. [What's new — Launcher & Suite](#whats-new--launcher--suite)
5. [Platform & packaging](#platform--packaging)
6. [Physics & numerics](#physics--numerics)
7. [Bug fixes](#bug-fixes)
8. [Breaking & behavioural changes](#breaking--behavioural-changes)
9. [Known limitations](#known-limitations)
10. [Deferred to beta 3 / 4](#deferred-to-beta-3--4)
11. [Full commit changelog](#full-commit-changelog)
12. [Upgrade notes](#upgrade-notes)
13. [Acknowledgements](#acknowledgements)

---

## Highlights

| Area | Headline |
|------|----------|
| **2D sub-app** | Now feature-complete to the beta-2 95 % bar. |
| **Black hole library** | Every mainstream stellar / intermediate / supermassive candidate (Sgr A\*, M87\*, Gaia BH1/BH2/BH3, Cyg X-1, GW190521, …) ships as a preset with a "Learn more" reference URL. |
| **Custom BH authoring** | Full menu-driven authoring of custom black holes and custom orbital bodies; presets persist across runs in `saves/`. |
| **Mergers** | Five secondary kinds (BH / NS / Pulsar / Star / WD) with kind-specific post-merger dynamics (GW kick + ring, X-ray burst, tidal-disruption stream). Three time-scale modes (cinematic / default / realtime). Auto-follow camera. |
| **Relativistic visuals** | Kerr equatorial geodesic overlay (J), thermal disk photon emitter (O), sticky high-resolution lensing near `b_crit` (K). |
| **HUD** | Dark / light theming, animated parallax starfield, View-reset key, transient-event cleanup on preset switch. |
| **Suite** | Unified Qt launcher, hosted 2D/3D widgets, multi-tab analysis, auto-updater. |
| **Platforms** | First-class **Windows** support (vcpkg + VS 2022, MSI/EXE installer), continued macOS `.app`/`.pkg`, refreshed Flatpak manifest at `io.github.0xLiam0920.AetherionSuite`. |

---

## What's new — 2D

### Black hole presets

- A full table of physically-grounded presets in `src/2D/2D-utils/presets_2d.hpp`,
  including stellar-mass, intermediate-mass, and supermassive candidates plus
  recent LIGO/Virgo events. Each entry carries:
  - mass in solar units (geometrized internally),
  - influence/tidal/disk zones in geometric `M`,
  - a `learnMoreUrl` opened by the `L` key.
- A galaxy-system spawn path for galactic-centre presets (Sgr A\*, M87\*, etc.)
  populates representative orbiting bodies on activation.

### Custom black hole authoring

- New custom-BH dialog (`src/QT-LAUNCHER/custom_bh_dialog.{h,cpp}`, ~580 lines)
  with full parameter sweep — mass, spin, disk geometry, lensing strength,
  preview thumbnail rendered live via `bh_preview_widget` (~2 k lines).
- Authored presets are written to
  `<userdata>/Aetherion/saves/custom_presets.tsv`
  and reloaded at startup. Tab-separated for forward-compatibility and to
  avoid pulling in a JSON dependency.

### Custom bodies (key `N`)

- In-sim menu to spawn arbitrary orbital bodies (seven curated types, e.g.
  *G-type star*, *neutron star*, *pulsar*, *brown dwarf*, *gas giant*). Three
  fields: type (Left/Right), semi-major axis in `M`, eccentricity.
- `Simulation::addCustomBody(type, smM, ecc)` tags spawned bodies as galaxy
  bodies so existing render and event-handling code paths work unchanged.

### Merger events

- `MergerSecondaryKind`: `BlackHole`, `NeutronStar`, `Pulsar`, `Star`,
  `WhiteDwarf`.
- `TidalEventKind::XrayBurst` for NS/Pulsar inspirals.
- `reinitBodiesPostMerger()` branches per kind:
  - **BH** → GW kick + photon ring,
  - **NS / Pulsar** → fast X-ray burst, no kick,
  - **Star / WD** → tidal-disruption fallback stream centred at origin.
- `MERGER_SECONDARY_PRESETS[]` ships canonical secondaries
  (NS, *PSR J0108-1431*, Sun, *Sirius B*).
- Three time-scale dials via `M` (cinematic 0.35×, default 1.0×, realtime
  ~3.0× — approximate, full Peters `M_tot` scaling is beta-3 work).
- Auto-follow camera: `pixelsPerM` smoothly lerps so the inspiral binary
  stays framed at ~40 % of the viewport. Capped so the horizon can't
  blow up past ~250 px as `r → 0`. Wall-clock flash and remnant shot
  keep the pre-zoom view so the cinematic beat reads cleanly.

### Relativistic visualisation

- **Kerr equatorial overlay** (`J`): 24 null geodesics swept across
  `b ∈ [−8M, 8M]`, integrated in the equatorial plane of a Kerr metric with
  `a = 0.7 M` (configurable via `kerrSpin`). Drawn in amber over the
  Schwarzschild fan for direct visual comparison.
- **Disk photon emitter** (`O`): Monte-Carlo thermal disk emission from
  `r ∈ [r_ISCO, 12 M]`, with rest-frame `λ(r) ∝ r^{3/4}`. Rays use
  `Photon::computeEmissionPath` and feed the lensing analytics histogram.
- **High-resolution lensing** (`K`): adds 20 extra rays per side, log-spaced
  near the critical impact parameter `b_crit = 3√3·M`, with finer angular
  step `dφ = 0.001`. Now a **sticky** mode — preserved across preset
  switches, zoom changes, and test-scenario swaps. Reset via `R`.

### HUD, controls, theming

- **Light/dark mode** (`B`): re-tints clear colour, starfield, info text,
  data panel, controls panel, merger menu, and custom-body menu. Persisted
  in `QSettings("Aetherion","AetherionSuite")` and applied to every open
  2D tab in the Qt launcher.
- **Animated starfield**: three parallax layers (120 / 90 / 60 stars,
  drift 1.5 / 0.6 / 0.25 px/s, twinkle 0.6 / 0.4 / 0.25 Hz). Independent
  per-star phase. Runs while paused, theme-aware.
- **View reset** (`V`): snaps `pixelsPerM` back to the preset's horizon
  target (or `defaultPixelsPerM` if no preset is active) without
  disturbing physics state.
- **Transient cleanup**: `Simulation::clearTransientEvents()` is invoked
  on preset switch, scenario reset, and `R`, so stale "MERGER REMNANT"
  or "X-RAY BURST" overlays can no longer leak after a mid-flash
  switch.
- **Preset combo in the launcher** ("Open Preset Tab") now actually
  applies the chosen preset to the spawned tab instead of opening
  blank.

### Keybindings

The 2D keymap is now fully exposed in `key_config_2d.cfg` with new entries
for every beta-2 feature:

```
toggle_high_res_lensing   = K
toggle_kerr_overlay       = J
toggle_disk_emitter       = O
toggle_light_mode         = B
reset_view                = V
cycle_merger_speed        = M
learn_more                = L
open_custom_body_menu     = N
```

Re-binding via the in-app keybind editor writes back to the same file.

---

## What's new — 3D

The 3D sub-app received an infrastructure pass to set up for its beta-3
feature-completion push. Highlights:

- **Modularisation** — the monolithic 3D entry point was split into
  `bh3d_core.hpp` (~1150 lines), `physics_overlay.hpp` (~580 lines),
  `hud_imgui.hpp` / `hud_panel.hpp`, and a new `imgui_qt_adapter.hpp`
  so the 3D view can host an embedded ImGui HUD inside the Qt widget.
- **Photoreal disk shader** — `BlackHole3D_PhotorealDisk.frag`
  (~2.7 k lines) now ships with the build; selectable as a render mode.
- **Resource manager** — texture / shader lifetime separated from the
  GL context so multiple 3D tabs no longer fight over the same handles.
- **Camera controller** — rewritten with stable spherical interpolation
  and a configurable focus point for galactic-centre views.
- **Transparent-orbiter regression fixed** — orbiting bodies in front
  of the black hole no longer render translucent.

The 3D simulation's broader feature completeness remains
**deferred to beta 3** by design.

---

## What's new — Launcher & Suite

- **Unified Qt launcher** hosting both sub-apps:
  - `Simulation2DWidget` (embedded SFML via `QSFMLCanvas`),
  - `Simulation3DWidget` (embedded OpenGL surface),
  - tabbed multi-instance editing — open multiple 2D simulations side by side
    with independent state.
- **Multi-tab analysis** with per-tab CSV / FITS / binary export.
- **Page registry** in `src/MAIN-MENU/` for menu pages
  (Overview, Simulation, Object Library, Data Analysis, Settings, Export,
  Breakdown). Pages are now discoverable from a single registry instead
  of being open-coded into `launcher.hpp`.
- **Black-hole preview widget** (`bh_preview_widget`) — live thumbnail of
  any preset or custom BH, used in both the object library and the
  custom-BH dialog.
- **Keybind UI** (`keybindbuttonwidget`) rewritten with theming and
  conflict-detection.
- **Auto-updater** (`updater.{h,cpp}`) checks the GitHub Releases API on
  startup (opt-in), surfaces a notification, and can download and stage
  the update. Skip / Install Later / Install Now flow.
- **Workspace persistence** — `exportFormat`, `lightMode`, and per-tab
  simulation state round-trip through `QSettings` and the workspace JSON.

---

## Platform & packaging

| Platform | Status | Notes |
|----------|--------|-------|
| **macOS (Apple Silicon + Intel)** | ✅ Primary | `make_pkg.sh` produces `.pkg`; `bundle_app.sh` produces an ad-hoc-signed `.app`. Homebrew dependencies. |
| **Linux (Flatpak)** | ✅ | Manifest renamed to `io.github.0xLiam0920.AetherionSuite.{json,desktop,metainfo.xml}`. Checksums helper at `flatpak/get-checksums.sh`. `rebuild_linux.sh` for ad-hoc builds. |
| **Windows (x64)** | ✅ **NEW** | vcpkg toolchain + Visual Studio 2022. `rebuild_windows.ps1` / `rebuild_windows.sh` (Git Bash). `make_exe.sh` builds an installer. SmartScreen is unsigned — see README. |
| **Linux (native build)** | ✅ | `rebuild_linux.sh`; system Qt6 / SFML 3 / OpenGL. |

Other platform work:

- `src/platform.hpp` — `platformUserDataDir()` + `platformOpenUrl()`
  with proper branches per OS. Raw `getenv("HOME")` use has been removed
  from the codebase; all save / config paths route through the helper.
- `Info.plist.in` updated for both the main app and the Qt launcher
  (camera-free, network entitlement for the updater only).
- `CMakeLists.txt` cleaned up: single source of truth for all four
  executable targets (`blackhole-sim`, `blackhole-2D`, `blackhole-3D`,
  `physics-regression-tests`) plus the `aetherion_imgui` static lib.
- `imgui-sfml` is now a tracked submodule (`.gitmodules`) pinned to a
  known-good upstream master.

---

## Physics & numerics

- **`Photon::computePath`** (`src/2D/2D-simulation/photon.hpp`)
  - Now scales its integration cutoff `rMax = max(rMax, 200·M)` so it
    behaves identically across all preset mass scales. Previously the
    default `1e5` was below typical periapsis at SMBH scale, causing
    the ray loop to terminate on step 0.
- **`Photon::computeEmissionPath`**
  - Same M-scaled clamp on `rMax`; the disk emitter now renders
    full ray paths at every preset mass.
- **`Kerr::integrateKerrEquatorial`** (`src/2D/2D-physics/kerr.hpp`)
  - `rMax = max(rMax, 200·M)` and `dlam = max(dlam, 0.25·M)`. The Kerr
    overlay now produces correctly-sampled trajectories at preset
    scales (previously the default `rMax = 1e3` was below `r0 = 60M`
    at SMBH scales).
- **`Simulation::Params`** gains a sticky `highResLensing` flag so
  `rebuildPhotons()` honours the user's K toggle even when called
  without the explicit `highRes` argument (the case for every
  preset / scenario / zoom rebuild path).
- **Schwarzschild metric** (`schwarzschild.hpp`) — refactored:
  `criticalImpact()`, `findPeriapsis(b)`, `f(r)`, `isco()`, `horizon()`
  are all single-source-of-truth helpers; periapsis bisection runs
  200 iterations (machine precision) for numerical headroom.
- **Geodesic integrator** (`geodesic.hpp`, `integrator.hpp`) —
  documented Binet-equation conventions, time-budget guard
  (40–50 ms per ray) so a near-critical photon can no longer stall
  the frame.
- **Pulsar orbital module** (`pulsar_orbital.hpp`, ~720 lines)
  rewritten with corrected secular precession and a stable
  bright-core + dual-beam visualiser shared between the live pulsar
  body and the merger secondary view.
- **Lensing analytics** (`research_data.hpp`) — emitter rays are
  tagged `fromEmitter = true` and excluded from the global lensing
  histogram so the spectrum, deflection curve, and critical-impact
  estimate stay clean.

### Regression suite

`tests/physics_regression_tests.cpp` runs as part of the build. Latest
green output:

```
[PASS] solar deflection: 0.000486012 deg (target 0.000486111 deg)
[PASS] conservation:     max dE/E = 2.27e-16,  max dL/L = 1.13e-08
Physics regression tests PASSED
```

The solar light-deflection test reproduces the canonical
`1.75″ ≈ 4.86 × 10⁻⁴°` value to four significant figures.

---

## Bug fixes

This window resolves dozens of issues touched by ~30 k lines of diff.
Selected highlights:

- **HUD-disable crash in 3D** (V0.1.1) — fixed.
- **Mergers** (V0.1.4 / V0.1.5) — secondary-body gravity, post-merger
  reset, flash-phase wall-clock, kick direction and trail rendering.
- **Gaia BH bodies** (V0.1.4) — mass/spectral parameters corrected for
  Gaia BH1, BH2, BH3.
- **Auto-updater** (V0.1.4) — version compare, download integrity.
- **Primary-BH cosmetic wobble** — now gated to `q ≥ 0.05` so equal-mass
  pairs no longer jitter the primary visually in both render sites.
- **Transparent orbiters in 3D** — orbital bodies in front of the
  horizon no longer render translucent.
- **Preset-mode visual regression** — `K` (high-res lensing), `J`
  (Kerr overlay), and `O` (disk emitter) all now render correctly
  when a preset is active and survive preset switches.
- **Preset combo opening blank tab** — fixed; now propagates the
  chosen preset to the spawned 2D tab.
- **Light-mode persistence** — `lightMode` now round-trips through
  workspace state.
- **Export format wiring** — the `exportFormat` choice is honoured by
  newly-opened tabs (was being dropped from `QSettings`).
- **Key conflict** — `K` was double-bound to high-res lensing and
  another action; resolved.
- **Stale merger overlay** — switching presets mid-flash no longer
  leaks "MERGER REMNANT" / "X-RAY BURST" labels.
- **`-Wswitch` warning** at `BlackHole2D.cpp:86` — cleared. Build is
  warning-clean.
- **Cross-cloud user-data paths** — all `getenv("HOME")` removed; uses
  `platformUserDataDir()`.
- **`imgui-sfml` tracking** (V0.1.5) — switched from embedded copy
  to a real submodule pinned to upstream master, fixing intermittent
  rebuild churn.

---

## Breaking & behavioural changes

- **Default keybinds added/changed.** Pre-beta-2 configs may not have
  entries for `J / K / O / B / V / M / L / N`. The loader treats
  missing entries as defaults (see [Keybindings](#keybindings)), so
  no migration step is required, but rebinding `K` is recommended if
  your prior config bound it elsewhere — it now toggles high-res
  lensing.
- **`Photon::computePath` / `computeEmissionPath` `rMax`.** The
  effective cutoff is now `max(caller_rMax, 200·M)`. Caller-passed
  values smaller than `200·M` are ignored — this is intentional and
  required for the integrator to function at SMBH mass scales.
- **`Simulation::Params::highResLensing`** field added. ABI-incompatible
  for any external linker; this is a beta release, no stable ABI is
  promised.
- **Flatpak app-id renamed** to `io.github.0xLiam0920.AetherionSuite`
  (was `io.github.aetherion` in early Flatpak drafts). Reinstall, do
  not upgrade in place.
- **`fuckingrebuild.sh`** was renamed to `rebuild.sh` (then split into
  `rebuild_linux.sh` / `rebuild_macos.sh` / `rebuild_windows.{sh,ps1}`).

---

## Known limitations

- **Windows code-signing** — the executable is unsigned. SmartScreen
  will warn on first run. See README §Install.
- **macOS notarization** — ad-hoc signed only. Gatekeeper will warn.
- **Merger realtime mode** is an approximate `3×` factor; full
  Peters `M_tot` scaling is **deferred to beta 3**.
- **Primary BH stays at origin** in physics (`TODO(physics-honesty)`,
  `simulation.hpp:427`). Cosmetic wobble only. True barycentric
  integration is deferred to beta 3.
- **Hot-plasma temperature model** assumes ideal-gas behaviour
  (`schwarzschild.hpp:151`). Acceptable for Bondi-radius
  order-of-magnitude estimates; not a precision plasma-dynamics model.
- **Custom-body type list** is curated to seven types; arbitrary
  spectral classes are a beta-3 polish item.
- **`research_data.hpp` full FITS plan** (dwarf-galaxy / loss-cone /
  full Bondi-Hoyle export) is staged but not shipped. See
  `docs/2d-research-roadmap.md` for the full plan.

---

## Deferred to beta 3 / 4

| Item | Target |
|------|--------|
| 3D sub-app feature completeness | beta 3 |
| Full Peters `M_tot` merger scaling | beta 3 |
| Barycentric primary-BH integration | beta 3 |
| `research_data.hpp` FITS dwarf-galaxy plan | beta 3 |
| "Mathematical logic" page | beta 4 |
| Beginner tutorial page | beta 4 |
| Sources / attribution page | beta 4 |
| Plugin / addon system (custom BHs, custom objects via SDK) | beta 4 / V1 |
| Final cross-suite polish & V1 readiness | V1 |

---

## Full commit changelog

Commits since [`BETA` (v0.1.0)](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/releases/tag/BETA),
ordered oldest → newest:

| Date | SHA | Summary |
|------|-----|---------|
| 2026-04-21 | `7814dd1` | README hotfix for v0.1.1 |
| 2026-04-21 | `9f2a975` | README update |
| 2026-04-22 | `6f64b18` | README update |
| 2026-04-24 | `f5bdd96` | README disclaimer refresh |
| 2026-04-24 | `d741d30` | **V0.1.2** — restored old docs, alpha prototype of custom BH system, bugfixes |
| 2026-04-24 | `89f797a` | Merge `origin/main` → V0.1.2 |
| 2026-04-25 | `0a1373b` | **V0.1.3** — added mainstream black-hole library to 2D/3D, configurability pass, bugfixes |
| 2026-05-04 | `831df25` | **V0.1.4** — auto-updater, Gaia BH fixes, physics rewrites, content patches |
| 2026-05-04 | `684d017` | README update |
| 2026-05-25 | `cd7ae57` | WIP — LQP shading, Linux compatibility |
| 2026-05-30 | `f1ce2e6` | **V0.1.4 follow-up** — merger fixes, mathematical-logic improvements |
| 2026-05-31 | `aaa03c9` | **V0.1.5** — merger physics, secondary-body gravity, `imgui-sfml` submodule tracking |
| 2026-05-31 | `282a96e` | Re-add `imgui-sfml` as submodule, bump to upstream master |
| 2026-06-01 | `7840c8e` | Pre-beta-2 patches, Windows compatibility pipeline starts |
| 2026-06-01 | `b1acc6b` | Convert `imgui-sfml` from embedded repo to regular folder |
| 2026-06-02 | `3deb94d` | Windows compatibility — best-effort baseline (Phase 1 complete) |
| **HEAD** | `unstaged` | Beta-2 final polish — preset visuals (K/J/O), HUD theming, M-scaled integrators, sticky high-res mode, sim-2D widget parity |

Tagged releases on the path to beta 2:
[`V0.1.3`](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/releases/tag/V0.1.3)
· [`V0.1.4`](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/releases/tag/V0.1.4)
· [`V0.1.5`](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/releases/tag/V0.1.5).

---

## Upgrade notes

### From v0.1.x

```bash
git pull
git submodule update --init --recursive    # imgui-sfml is now a submodule
```

Then rebuild for your platform:

```bash
./rebuild_macos.sh       # macOS
./rebuild_linux.sh       # Linux
./rebuild_windows.ps1    # Windows (PowerShell)
```

Or via CMake directly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

User data and saves are preserved automatically — paths now route through
`platformUserDataDir()` on every OS, but the underlying directories
(`~/Library/Application Support/Aetherion/`,
`$XDG_DATA_HOME/Aetherion/`, `%APPDATA%\Aetherion\`) are unchanged.

### From BETA (v0.1.0)

Same as above. Note the binary names — the unified launcher
(`blackhole-sim`) is the user-facing target; `blackhole-2D` and
`blackhole-3D` remain as standalone debugging entry points.

### Custom-preset migration

Pre-beta-2 custom BHs (if any) were not persisted across runs. They
now live in `<userdata>/Aetherion/saves/custom_presets.tsv`. No
migration is needed; re-author from the custom-BH dialog and they
will reload on next launch.

---

## Acknowledgements

- **SFML 3.x** for the rendering / windowing layer.
- **Qt 6** (Core, Gui, Widgets, Charts, Network) for the launcher and
  analysis surface.
- **Dear ImGui** + **imgui-sfml** for the in-sim debug HUD.
- The **LIGO / Virgo / KAGRA Collaboration** for the merger-event
  parameters used as preset reference values.
- The **Event Horizon Telescope Collaboration** for the M87\* and
  Sgr A\* parameter references.
- The **Gaia mission** for the BH1 / BH2 / BH3 detections.
- Everyone who filed an issue or sent a crash report against beta 1 —
  most of the bug-fix list above traces back to one of you.

---

## What's next

Beta 3 will focus on bringing the 3D sub-app up to the same 95 %
maturity bar that 2D now meets:

- Low-quality preset rendering rewrite
- Full custom-object pipeline in 3D
- Photoreal disk shader exposure in the UI
- Camera + HUD polish
- 3D-side analytics export

After that, beta 4 / V1 closes the long-tail polish, ships the
plugin/addon SDK, and lands the in-app explainer pages
(mathematical logic, tutorial, sources/attribution).

If you find a way to crash it spectacularly, [file an issue](https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite/issues).
The collection is growing.

— *Aetherion maintainers, June 2026*
