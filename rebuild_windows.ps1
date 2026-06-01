# ─────────────────────────────────────────────────────────────────────────────
# rebuild_windows.ps1 — clean build of Aetherion on Windows (MSVC + vcpkg).
#
# Prerequisites:
#   1. Visual Studio 2022 with the "Desktop development with C++" workload.
#      Launch this script from a "Developer PowerShell for VS 2022" prompt so
#      cl.exe / cmake.exe are on PATH.
#   2. vcpkg installed somewhere (e.g. C:\vcpkg). Set the VCPKG_ROOT env var
#      or pass -VcpkgRoot to this script.
#   3. The following vcpkg packages installed (x64-windows triplet):
#         vcpkg install qt6-base qt6-charts sfml glm glew
#
# Usage (from the repo root):
#   pwsh -File .\rebuild_windows.ps1
#   pwsh -File .\rebuild_windows.ps1 -VcpkgRoot C:\dev\vcpkg -Config Debug
# ─────────────────────────────────────────────────────────────────────────────
[CmdletBinding()]
param(
    [string]$VcpkgRoot = $env:VCPKG_ROOT,
    [ValidateSet('Release', 'Debug', 'RelWithDebInfo')]
    [string]$Config   = 'Release',
    [string]$Triplet  = 'x64-windows',
    [string]$BuildDir = 'build'
)

$ErrorActionPreference = 'Stop'

if (-not $VcpkgRoot) {
    throw "VCPKG_ROOT is not set. Pass -VcpkgRoot <path> or set the env var."
}
$toolchain = Join-Path $VcpkgRoot 'scripts\buildsystems\vcpkg.cmake'
if (-not (Test-Path $toolchain)) {
    throw "vcpkg toolchain file not found at $toolchain"
}

if (-not (Test-Path $BuildDir)) {
    New-Item -ItemType Directory -Path $BuildDir | Out-Null
}

Write-Host "Configuring (vcpkg=$VcpkgRoot, triplet=$Triplet, config=$Config)..."
cmake -S . -B $BuildDir `
    -G "Visual Studio 17 2022" -A x64 `
    "-DCMAKE_TOOLCHAIN_FILE=$toolchain" `
    "-DVCPKG_TARGET_TRIPLET=$Triplet"

Write-Host "Building..."
cmake --build $BuildDir --config $Config --parallel

Write-Host ""
Write-Host "Build complete. Binaries are in:"
Write-Host "  $BuildDir\$Config\blackhole-sim.exe"
Write-Host "  $BuildDir\$Config\blackhole-2D.exe"
Write-Host "  $BuildDir\$Config\blackhole-3D.exe"
Write-Host ""
Write-Host "On first run you may need to deploy Qt runtime DLLs next to the exe:"
Write-Host "  windeployqt $BuildDir\$Config\blackhole-sim.exe"
