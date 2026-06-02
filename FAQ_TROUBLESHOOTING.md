# FAQ & Troubleshooting

Common questions, common mistakes, and the occasional rant. If your problem
isn't in here, congratulations, you've found a new and exciting one. Please
file an issue with actual error output, not a screenshot of a screenshot.

---

## Build Issues

### CMake can't find SFML

**Error:**
```
Could not find a package configuration file provided by "SFML"
```

**Solution:**
- Make sure SFML **3.x** is installed (not 2.x).
- macOS: `brew install sfml`
- Linux: `sudo apt install libsfml-dev`, check that your distro provides SFML 3.x. If not, build SFML from source.
- Windows: Use vcpkg or build SFML 3.x from source. Install via vcpkg:
  ```bash
  vcpkg install sfml glm
  ```
- FreeBSD: yeah no I ain't adding compat for that bud, best you switch to something mainstream

- If SFML is installed in a non-standard location, pass its path to CMake:
  ```bash
  cmake -DCMAKE_PREFIX_PATH=/path/to/sfml ..
  ```

### CMake can't find GLM

**Error:**
```
Could not find a package configuration file provided by "GLM"
```

**Solution:**
- macOS: `brew install glm`
- Linux: `sudo apt install libglm-dev`
- The CMake config tries both `glm` (lowercase) and `GLM` (uppercase) package names. If neither works, set the include path manually:
  ```bash
  cmake -DGLM_INCLUDE_DIRS=/path/to/glm/include ..
  ```

### SFML 2.x vs 3.x API errors

**Error (examples):**
```
error: no member named 'antiAliasingLevel' in 'sf::ContextSettings'
error: no matching constructor for 'sf::VideoMode'
```

**Cause:** The project requires SFML 3.x. SFML 2.x has a different API.

**Solution:**
- Upgrade to SFML 3.x.
- macOS: `brew upgrade sfml` (Homebrew provides 3.x).
- Linux: If your package manager only has 2.x, build SFML 3.x from source (https://github.com/SFML/SFML).
- Windows: Use vcpkg or build SFML 3.x from source. Install via vcpkg:
  ```bash
  vcpkg install sfml glm
  ```

### OpenGL deprecation warnings on macOS

**Warning:**
```
'glGenFramebuffers' is deprecated: first deprecated in macOS 10.14
```

**Cause:** Apple deprecated OpenGL in favour of Metal. The project defines `GL_SILENCE_DEPRECATION` to suppress these, but some may still appear depending on compiler flags.

**Solution:** These warnings are cosmetic and can be safely ignored. The simulator works fine with macOS's OpenGL implementation.

### `build-essential` / no compiler found (Linux)

**Error:**
```
No CMAKE_CXX_COMPILER could be found.
```

**Cause:** You are trying to compile C++ without a C++ compiler, don't be a doofus.

**Solution:**
```bash
sudo apt install build-essential
```

---

## Runtime Issues

### Black screen on launch (`blackhole-3D`)

**Possible causes:**
1. **OpenGL 3.3 not supported**, The simulator requires OpenGL 3.3 Core. Check your GPU/driver support:
   ```bash
   glxinfo | grep "OpenGL version"   # Linux
   ```
   On macOS, OpenGL 3.3 is supported on all Macs from 2012 onwards.
   On windows, check your GPU specs online, or use a tool like GPU-Z. Corrupted drivers can also cause this, so use DDU or something along that line.

2. **Shader files missing**, The build copies `.frag` files next to the executable. If you moved the binary, make sure `BlackHole3D.frag` and `BlackHole3D_PhotorealDisk.frag` are in the same directory as the executable, or in the project's `src/3D/` directory.

3. **Missing textures**, Background or disk textures are optional but the shader may render a black disk without them. See the [Textures](#textures-not-loading) section.

### Shaders fail to compile at runtime

**Error (in terminal output):**
```
Shader compilation failed: ...
```

**Solution:**
- Ensure your GPU supports GLSL `#version 330 core`.
- If you edited a `.frag` file, check for syntax errors. The simulator loads shaders from disk at startup, so edits take effect on next launch without recompiling.
- On Linux with Intel integrated graphics, install the latest Mesa drivers:
  ```bash
  sudo apt install mesa-utils libgl1-mesa-dri
  ```

### Textures not loading

**Symptoms:** Black or missing accretion disk, no background stars.

**Solution:**
- **Disk texture:** Run the generator script (requires Python 3 + Pillow):
  ```bash
  pip install Pillow
  python makedisktexture.py
  ```
  This creates `src/disk_texture.png`. Rebuild to copy it to the build directory.
- **Background texture:** Place a `background.png` in the project root or `build/` directory. The build system copies it automatically.

### Low FPS / poor performance

**Tips:**
- The 3D simulator is GPU-bound (full-screen ray marching every frame). A discrete GPU is recommended for good performance.
- Reduce the window size, don't use fullscreen to avoid fps stuttering.
- Toggle off expensive features:
  - **J**, disable relativistic jets
  - **G**, disable Broad Line Region
  - **V**, disable Doppler beaming
- Check that your system isn't falling back to software rendering (integrated GPU), for windows sometimes task manager priority going to very high/realtime can fix this.

### Mouse cursor stuck / input not working

- Press **Esc** to release the mouse cursor from the window.
- Right-click and hold to enter freelook mode; release to stop.
- Press **F** to toggle between freelook and orbit camera modes.

---

## macOS-Specific Issues

### `.app` bundle won't open ("damaged" or "unidentified developer")

**Cause:** The bundle isn't code-signed, because Apple charges $100/year for the privilege of telling macOS that I am, in fact, a real person and not a Russian botnet.

**Solution:**
```bash
# Remove quarantine attribute
xattr -cr build/Aetherion.app
```

Or right-click the app → Open → confirm in the dialog.

### `bundle_app.sh` warnings about missing dylibs

**Warning:**
```
WARNING: /opt/homebrew/opt/glm/lib/libglm.dylib not found, skipping.
```

**Cause:** GLM is header-only on some Homebrew versions and doesn't produce a `.dylib`.

**Solution:** This warning is harmless for GLM since it's header-only at runtime. For other missing dylibs (SFML, FreeType, libpng), reinstall the relevant package:
```bash
brew reinstall sfml freetype libpng
```

### Homebrew on Apple Silicon vs Intel

The build system checks both `/opt/homebrew` (Apple Silicon) and `/usr/local` (Intel) prefix paths. If CMake still can't find packages:
```bash
cmake -DCMAKE_PREFIX_PATH="$(brew --prefix)" ..
```

---

## Linux-Specific Issues

### App crashes immediately ("Could not find the Qt platform plugin xcb")

**Error:**
```
qt.qpa.plugin: Could not find the Qt platform plugin "xcb" in ""
This application failed to start because no Qt platform plugin could be initialized.
```

**Cause:** The Qt6 XCB platform plugin is missing. The app forces `QT_QPA_PLATFORM=xcb` on Linux so SFML can embed into a Qt widget via X11/XWayland. If the `xcb` plugin isn't installed, Qt can't open any window.

**Solution:**
```bash
# Ubuntu / Debian
sudo apt install qt6-qpa-plugins

# Fedora
sudo dnf install qt6-qtbase-gui

# Arch Linux
sudo pacman -S qt6-base
```
If you're using the bundled tarball (`bundle_app.sh --platform linux`), the plugin is already included under `plugins/platforms/libqxcb.so`.

### No OpenGL context / `GLXBadFBConfig`

**Cause:** X11/Wayland driver issue or missing OpenGL support.

**Solution:**
```bash
sudo apt install libgl1-mesa-dri libgl1-mesa-glx mesa-utils
```
Verify with `glxinfo | grep "OpenGL version"`, you need 3.3+.

### SFML linking errors with system packages

**Error:**
```
undefined reference to `sf::RenderWindow::...'
```

**Cause:** Mismatched SFML version (system has 2.x, project needs 3.x).

**Solution:** Build SFML 3.x from source:
```bash
git clone https://github.com/SFML/SFML.git
cd SFML && mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local ..
cmake --build . && sudo cmake --install .
```
Then rebuild the project with:
```bash
cmake -DCMAKE_PREFIX_PATH=/usr/local ..
```

---

## Windows-Specific Issues

So you want to build Aetherion on Windows. First, deep breath. Second, follow this in order, because skipping a step here is the source of roughly 90% of the "doesn't work" issues people DM me about.

### What you actually need before running anything

The `rebuild_windows.ps1` / `rebuild_windows.sh` scripts assume you have already done all of the following. They will not install these for you. They are scripts, not your IT department.

1. **Visual Studio 2022** with the **"Desktop development with C++"** workload. Not VS Code. Not "I have Build Tools installed". The full IDE with that exact workload checked. The Community edition is free.
2. **vcpkg** cloned and bootstrapped somewhere sane (e.g. `C:\vcpkg`):
   ```powershell
   git clone https://github.com/microsoft/vcpkg.git C:\vcpkg
   C:\vcpkg\bootstrap-vcpkg.bat
   ```
3. **The actual libraries**, installed via vcpkg with the **x64-windows** triplet (the default `x86` triplet will silently install 32-bit binaries and then nothing will link, and you will be sad):
   ```powershell
   C:\vcpkg\vcpkg.exe install qt6-base qt6-charts sfml glm glew --triplet x64-windows
   ```
   Go make a coffee. Go make lunch. This step takes a while because Qt is, structurally, a small operating system.
4. **`VCPKG_ROOT` set in your environment** (or pass `--vcpkg-root` every time, which gets old fast):
   ```powershell
   [Environment]::SetEnvironmentVariable('VCPKG_ROOT', 'C:\vcpkg', 'User')
   ```
   Open a new shell after setting this. PowerShell will not magically pick it up in the current session.
5. A **"Developer PowerShell for VS 2022"** prompt (Start menu → search "Developer PowerShell"). This is the one that has `cmake.exe`, `cl.exe`, and friends on PATH. A regular PowerShell window does not, and you will get baffling "CMAKE_CXX_COMPILER not found" errors.

### Now actually build it

From the repo root in that Developer PowerShell window:
```powershell
.\rebuild_windows.ps1
```
Or from Git Bash / MSYS2 spawned from that same Developer prompt:
```bash
./rebuild_windows.sh
```
The output executables land in `build\Release\` as `blackhole-sim.exe`, `blackhole-2D.exe`, and `blackhole-3D.exe`.

### "CMake Error: Could not find package configuration file for Qt6" (or SFML, or GLM…)

**Cause:** Either `VCPKG_ROOT` isn't set, or you ran the script without passing `--vcpkg-root`, or you installed the wrong triplet. Almost always the triplet, in my experience.

**Solution:** Check that all three of these are true:
- `vcpkg list` shows the package with `:x64-windows` after its name (not `:x86-windows`).
- `echo $env:VCPKG_ROOT` in PowerShell prints a real path with a `scripts\buildsystems\vcpkg.cmake` file inside it.
- You are in a Developer PowerShell for VS 2022 window.

If all three check out and it still fails, blow away the `build\` folder and configure from scratch. CMake caches are like cats; once they decide they don't like you, you replace them.

### "The application failed to start because Qt6Core.dll was not found"

**Cause:** Qt's runtime DLLs aren't next to the `.exe`. The build links against vcpkg's Qt, but Windows won't go hunting through vcpkg's install tree at runtime.

**Solution:** Deploy them with `windeployqt`:
```powershell
windeployqt build\Release\blackhole-sim.exe
windeployqt build\Release\blackhole-2D.exe
windeployqt build\Release\blackhole-3D.exe
```
`windeployqt` lives inside the vcpkg Qt install, usually under `vcpkg\installed\x64-windows\tools\Qt6\bin\windeployqt.exe`. Add that folder to PATH, or call it by full path. If you ran `bundle_app.sh --platform windows` or `make_exe.sh`, this is already done for you inside the bundled folder.

### `make_exe.sh` says "makensis was not found"

**Cause:** NSIS (the installer compiler) isn't installed. The script is honest about this and falls back to the portable `.zip` so you still get *something* useful.

**Solution:** Install NSIS, then re-run with `--no-bundle` so you don't have to wait for the bundle step again:
```powershell
choco install nsis            # if you have Chocolatey
winget install NSIS.NSIS      # or with winget
```
Then:
```bash
./make_exe.sh --no-bundle
```

### SmartScreen says "Windows protected your PC" when running the installer

**Cause:** The installer isn't code-signed. See the README rant. TL;DR: code-signing certificates cost real money for a hobby/research project, so the installer ships unsigned.

**Solution:** Click **"More info"** → **"Run anyway"**. If you don't trust the binary, that's fair, build it from source — that's literally why this is open source.

### `bash: ./rebuild_windows.sh: /usr/bin/env: bad interpreter`

**Cause:** Git for Windows checked out the script with CRLF line endings, and the `#!/usr/bin/env bash` shebang now ends in `\r`, which is not a real path.

**Solution:**
```bash
git config --global core.autocrlf input
git rm --cached rebuild_windows.sh make_exe.sh bundle_app.sh
git checkout -- rebuild_windows.sh make_exe.sh bundle_app.sh
```
Or just run the `.ps1` instead and avoid the whole class of problem.

### "Antivirus quarantined Aetherion.exe"

**Cause:** You compiled a brand-new unsigned executable that does GPU things and reads files from disk. To a heuristic AV, that is indistinguishable from approximately every piece of malware ever written.

**Solution:** Add an exclusion for the build folder, or submit a false-positive report to your AV vendor. I cannot whitelist binaries with Windows Defender from a README, sadly.

---

---------- RANT STARTS HERE ----------

Note: No, I am not adding support for other versions of SFML or any graphical library. 
MacOS only supports OpenGL 3.3, so anything 2x for SFML is not going to happen for cross-compatibility sake.
If you want to use an older or different version/alternative to any library here, fork this project and change the code to support it yourself. I will not be doing it for you.
If you want to add support directly to this project, make a Pull Request. Any help here is appreciated, but my priority is in maintaining the codebase and adding research 
grade features. Quality over quantity.


---------- RANT OVER ----------

---

## General FAQ

### What are the three executables?

| Executable | Description |
|---|---|
| `blackhole-sim` | Main launcher menu (ImGui interface to choose which simulator to run) |
| `blackhole-2D` | 2D Schwarzschild null-geodesic simulator | if behaves as expected should be able to run either standalone (terminal) or in the aetherion app
| `blackhole-3D` | 3D Kerr-metric ray-marched black hole with accretion disk, jets, bloom |, same case as above.

### Can I edit shaders without recompiling?

The `.frag` shader files are loaded from disk at runtime. Edit `BlackHole3D.frag` or `BlackHole3D_PhotorealDisk.frag` in the build directory and relaunch the simulator. Do so at your own risk; GLSL is a fickle beast and "why is the screen entirely magenta" is almost always your fault, not mine.

### What units does the simulation use?

Simulation units are in standarized Schwarzschild radii (Rs = 1). For reference, TON 618's Schwarzschild radius is approximately 1,300 AU (1.95×10¹⁴ m).

### How do I switch between black hole presets?

The 3D simulator includes presets for TON 618, Sgr A*, M87*, and others. Check the HUD overlay (**H** to toggle) and key bindings for preset switching. Presets are defined in [src/3D/presets.hpp](src/3D/presets.hpp).
  Note: You can change the keybinds to whatever you want at any time. Just edit the keybinds in the aetherion app, though try not to create conflicts.

### Does this work on Windows?

Yes. Windows is now a supported platform alongside macOS and Linux. The cross-platform C++17 codebase builds with MSVC (recommended) or MinGW. You will need:
- Visual Studio 2022 with the "Desktop development with C++" workload (or MinGW with C++17 support)
- [vcpkg](https://vcpkg.io/) with: `vcpkg install qt6-base qt6-charts sfml glm glew --triplet x64-windows`
- CMake ≥ 3.10
- A GPU with OpenGL 3.3 support

From a "Developer PowerShell for VS 2022" prompt run `rebuild_windows.ps1`, or from Git Bash / MSYS2 run `./rebuild_windows.sh`. To produce a redistributable installer, run `./make_exe.sh` (uses NSIS if available, otherwise falls back to a portable `.zip`).

### How do I regenerate the accretion disk texture?

do this:
```bash
pip install Pillow
python makedisktexture.py
```

This writes the `src/disk_texture.png`. Rebuild to have it copied to `build/`.

### Can you add `insert blackhole here` 

**Answer:** Maybe. If it's a significant enough discovery and is highly requested, I can add a preset for it.
Otherwise, try making a custom preset based on properties yourself. In future, I may make a plugin system
for adding more properties, and maybe a hub of sorts to share custom black hole files if this somehow gets enough traction.


### FreeBSD?

**Answer:** I have zero idea why I got asked by several people for this, but absolutely not.
I have zero interest in supporting FreeBSD. If you want to run this on FreeBSD, 
you can try building it yourself and fixing any issues that come up, but I will not be supporting you on this.
This is an open source research tool, not a game.