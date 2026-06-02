#!/usr/bin/env bash
# =============================================================================
# rebuild_windows.sh - clean build of Aetherion on Windows (MSVC + vcpkg).
#
# Bash counterpart to rebuild_windows.ps1, for the noble souls who installed
# Git Bash specifically so they would never have to type `$env:VAR = "value"`
# ever again. Functionally identical to the .ps1: points CMake at vcpkg's
# toolchain file and lets the Visual Studio 2022 generator do its thing.
#
# Prerequisites:
#   1. Visual Studio 2022 with the "Desktop development with C++" workload.
#      Launch this script from a "Developer Command Prompt for VS 2022" shell
#      (e.g. Git Bash spawned from one) so cmake.exe / cl.exe are on PATH.
#   2. vcpkg installed somewhere (e.g. C:\vcpkg). Set VCPKG_ROOT in the env
#      or pass --vcpkg-root.
#   3. The required vcpkg packages installed (x64-windows triplet):
#         vcpkg install qt6-base qt6-charts sfml glm glew --triplet x64-windows
#
# Usage (from the repo root, in Git Bash / MSYS2):
#   ./rebuild_windows.sh
#   ./rebuild_windows.sh --vcpkg-root /c/dev/vcpkg --config Debug
# =============================================================================
set -euo pipefail

VCPKG_ROOT_ARG="${VCPKG_ROOT:-}"
CONFIG="Release"
TRIPLET="x64-windows"
BUILD_DIR="build"
GENERATOR="Visual Studio 17 2022"
ARCH="x64"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcpkg-root)    VCPKG_ROOT_ARG="${2:?'--vcpkg-root requires a path'}"; shift 2 ;;
        --vcpkg-root=*)  VCPKG_ROOT_ARG="${1#--vcpkg-root=}"; shift ;;
        --config)        CONFIG="${2:?'--config requires Release|Debug|RelWithDebInfo'}"; shift 2 ;;
        --config=*)      CONFIG="${1#--config=}"; shift ;;
        --triplet)       TRIPLET="${2:?'--triplet requires a value'}"; shift 2 ;;
        --triplet=*)     TRIPLET="${1#--triplet=}"; shift ;;
        --build-dir)     BUILD_DIR="${2:?'--build-dir requires a path'}"; shift 2 ;;
        --build-dir=*)   BUILD_DIR="${1#--build-dir=}"; shift ;;
        --generator)     GENERATOR="${2:?'--generator requires a value'}"; shift 2 ;;
        --generator=*)   GENERATOR="${1#--generator=}"; shift ;;
        -h|--help)
            sed -n '2,21p' "$0" | sed 's/^# \?//'
            exit 0 ;;
        *) printf 'ERROR: Unknown argument: %s\nRun with --help for usage.\n' "$1" >&2; exit 1 ;;
    esac
done

case "$CONFIG" in
    Release|Debug|RelWithDebInfo) ;;
    *) printf 'ERROR: --config must be Release, Debug, or RelWithDebInfo (got "%s").\n' "$CONFIG" >&2; exit 1 ;;
esac

# ── Sanity: we should be on Windows (Git Bash / MSYS2 / Cygwin) ─────────────
# If you're running this on Linux and got here by accident, you probably
# wanted rebuild_linux.sh. We just warn instead of exit because if you really
# insist on doing something cursed, that is your prerogative.
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) ;;
    *) printf 'WARNING: rebuild_windows.sh is intended for Git Bash / MSYS2 on Windows.\n' >&2 ;;
esac

if [[ -z "${VCPKG_ROOT_ARG}" ]]; then
    printf 'ERROR: VCPKG_ROOT is not set. Pass --vcpkg-root <path> or export VCPKG_ROOT.\n' >&2
    exit 1
fi

# Normalize MSYS-style paths (/c/vcpkg) to Windows-style for CMake.
# CMake on Windows does not speak forward-slash-with-no-drive-letter. This is
# a hill it will die on, so we translate.
to_windows_path() {
    local p="$1"
    if command -v cygpath &>/dev/null; then
        cygpath -w "$p"
    else
        printf '%s' "$p"
    fi
}

VCPKG_WIN="$(to_windows_path "${VCPKG_ROOT_ARG}")"
TOOLCHAIN="${VCPKG_ROOT_ARG}/scripts/buildsystems/vcpkg.cmake"
if [[ ! -f "${TOOLCHAIN}" ]]; then
    printf 'ERROR: vcpkg toolchain file not found at %s\n' "${TOOLCHAIN}" >&2
    exit 1
fi
TOOLCHAIN_WIN="$(to_windows_path "${TOOLCHAIN}")"

mkdir -p "${BUILD_DIR}"

# Ensure submodules (external/imgui-sfml etc.) are populated, mirrors rebuild_linux.sh.
if [[ -d .git ]] && command -v git &>/dev/null; then
    git submodule update --init --recursive
fi

echo "Configuring (vcpkg=${VCPKG_WIN}, triplet=${TRIPLET}, config=${CONFIG}) …"
cmake -S . -B "${BUILD_DIR}" \
    -G "${GENERATOR}" -A "${ARCH}" \
    "-DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_WIN}" \
    "-DVCPKG_TARGET_TRIPLET=${TRIPLET}"

echo "Building (${CONFIG}) …"
cmake --build "${BUILD_DIR}" --config "${CONFIG}" --parallel

cat <<EOF

Build complete. Binaries are in:
  ${BUILD_DIR}/${CONFIG}/blackhole-sim.exe
  ${BUILD_DIR}/${CONFIG}/blackhole-2D.exe
  ${BUILD_DIR}/${CONFIG}/blackhole-3D.exe

On first run you may need to deploy Qt runtime DLLs next to the exe:
  windeployqt ${BUILD_DIR}/${CONFIG}/blackhole-sim.exe

Next steps:
  ./bundle_app.sh --platform windows     # portable folder + zip
  ./make_exe.sh                          # NSIS .exe installer (if makensis is on PATH)
EOF
