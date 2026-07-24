cmake_minimum_required(VERSION 3.10)
project(blackhole-sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# ── Windows / MSVC global compile flags ───────────────────────
# _USE_MATH_DEFINES: expose M_PI and friends from <cmath> (MSVC omits them).
# NOMINMAX:          prevent <windows.h> from defining min/max macros that
#                    break std::min / std::max and ADL-resolved calls.
# _CRT_SECURE_NO_WARNINGS: silence MSVC's deprecation warnings about
#                          standard C functions (getenv, etc.) used cross-platform.
# WIN32_LEAN_AND_MEAN:     trim down <windows.h> include cost.
if(WIN32)
    add_compile_definitions(
        _USE_MATH_DEFINES
        NOMINMAX
        WIN32_LEAN_AND_MEAN
    )
    if(MSVC)
        add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
    endif()
endif()

# ── Platform-specific package search hints ────────────────────
if(APPLE)
    # Help CMake find Homebrew-provided packages when configured via IDE tools.
    if(EXISTS "/opt/homebrew")
        list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew")
    endif()
    if(EXISTS "/usr/local")
        list(APPEND CMAKE_PREFIX_PATH "/usr/local")
    endif()
elseif(UNIX)
    # Common Linux distro install prefixes for Qt6, SFML, GLEW, GLM.
    # Adjust or remove entries that don't exist on your distro.
    foreach(_prefix "/usr" "/usr/local" "/usr/lib/x86_64-linux-gnu"
                    "/usr/lib/aarch64-linux-gnu")
        if(EXISTS "${_prefix}")
            list(APPEND CMAKE_PREFIX_PATH "${_prefix}")
        endif()
    endforeach()
elseif(WIN32)
    # ── Windows: MSYS2 / mingw-w64 ─────────────────────────────
    # When building with MSYS2 + mingw-w64, CMake is typically invoked from
    # the MSYS2 shell and already has the right paths. These are fallbacks.
    foreach(_prefix
            "C:/msys64/mingw64"
            "C:/msys64/ucrt64"
            "C:/Qt/6.9.0/mingw_64"
            "C:/Qt/6.9.0/msvc2022_64")
        if(EXISTS "${_prefix}")
            list(APPEND CMAKE_PREFIX_PATH "${_prefix}")
        endif()
    endforeach()
    # ── Windows: vcpkg ──────────────────────────────────────────
    # If using vcpkg, pass -DCMAKE_TOOLCHAIN_FILE=<vcpkg>/scripts/buildsystems/vcpkg.cmake
    # at configure time. No extra hints needed here in that case.
endif()


# Find Qt6
find_package(Qt6 COMPONENTS Core Gui Widgets Charts Network REQUIRED)

# SFML 3.x is required. The API and component naming changed between 2.x and 3.x:
#   SFML 2.x components: graphics window system     (lowercase)
#   SFML 3.x components: Graphics Window System     (capitalized)
# Some distros (e.g. Ubuntu 24.04 libsfml-dev) still ship SFML 2.x, which would
# fail the capitalized component lookup with a cryptic error. Probe explicitly
# so we can emit an actionable message.
find_package(SFML 3 COMPONENTS Graphics Window System QUIET)
if(NOT SFML_FOUND)
    # Was a 2.x install detected? If so, tell the user instead of letting the
    # find_package error spam scroll past.
    find_package(SFML 2 COMPONENTS graphics window system QUIET)
    if(SFML_FOUND)
        message(FATAL_ERROR
            "Found SFML ${SFML_VERSION}, but Aetherion requires SFML 3.x.\n"
            "  On Ubuntu/Debian < 25.04, libsfml-dev ships SFML 2.x, build SFML 3 "
            "from source (https://github.com/SFML/SFML/releases/tag/3.0.0) or use "
            "the Flatpak build (see README.md).")
    endif()
    # Neither 3.x nor 2.x found, fall through to the canonical REQUIRED call
    # so CMake produces its normal "package not found" diagnostic.
    find_package(SFML 3 COMPONENTS Graphics Window System REQUIRED)
endif()

# Find GLM (package name differs between distros)
find_package(glm CONFIG QUIET)
if(NOT glm_FOUND)
    find_package(GLM REQUIRED)
endif()

# Find OpenGL
find_package(OpenGL REQUIRED)

# ── GLEW (Linux / Windows, not needed on macOS) ─────────────
# macOS uses the system OpenGL framework directly (Core Profile via
# OpenGL/gl3.h). All other platforms need a runtime loader.
if(NOT APPLE)
    find_package(GLEW REQUIRED)
endif()

# ── Unified Qt6-based Application with Embedded Simulations ──
# All windows and simulations are tied to one executable instance
# WIN32 tells MSVC/MinGW to use the Windows GUI subsystem (no console window).
# NOTE: WIN32 must be a literal CMake keyword (parsed before generator
# expressions are evaluated). It is silently ignored on non-Windows platforms.
add_executable(blackhole-sim WIN32
    # Qt Launcher & Main Window
    src/QT-LAUNCHER/main.cpp
    src/QT-LAUNCHER/mainwindow.h
    src/QT-LAUNCHER/mainwindow.cpp
    
    # SFML Canvas for embedding SFML in Qt
    src/QT-LAUNCHER/sfml_canvas.h
    src/QT-LAUNCHER/sfml_canvas.cpp
    
    # Simulation Widgets
    src/QT-LAUNCHER/simulation_2d_widget.h
    src/QT-LAUNCHER/simulation_2d_widget.cpp
    src/QT-LAUNCHER/simulation_3d_widget.h
    src/QT-LAUNCHER/simulation_3d_widget.cpp

    # Keybind remapping widget (requires AUTOMOC for Q_OBJECT)
    src/QT-LAUNCHER/keybindbuttonwidget.h
    src/QT-LAUNCHER/keybindbuttonwidget.cpp

    # Custom black hole builder dialog
    src/QT-LAUNCHER/custom_bh_dialog.h
    src/QT-LAUNCHER/custom_bh_dialog.cpp

    # Auto-updater (checks GitHub Releases)
    src/QT-LAUNCHER/updater.h
    src/QT-LAUNCHER/updater.cpp

    # Object Library preview widget
    src/QT-LAUNCHER/bh_preview_widget.h
    src/QT-LAUNCHER/bh_preview_widget.cpp

    # 2D/3D integration files should be linked as libraries or refactored modules,
    # not standalone entrypoints with their own main().
    src/3D/shader_utils.cpp
    src/3D/texture_utils.cpp
)

target_include_directories(blackhole-sim PRIVATE 
    ${CMAKE_SOURCE_DIR}/src/QT-LAUNCHER
    ${CMAKE_SOURCE_DIR}/src/2D
    ${CMAKE_SOURCE_DIR}/src/3D
    ${CMAKE_SOURCE_DIR}/src          # exposes platform.hpp as "platform.hpp"
    ${CMAKE_SOURCE_DIR}/external/json
    ${CMAKE_SOURCE_DIR}
)

target_link_libraries(blackhole-sim PRIVATE
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::Network
    SFML::Graphics SFML::Window SFML::System
    OpenGL::GL
    aetherion_imgui
)

# ── Application version (must match the GitHub release tag, without 'v') ──
set(AETHERION_VERSION "0.1.3")
target_compile_definitions(blackhole-sim PRIVATE
    AETHERION_VERSION="${AETHERION_VERSION}"
    AETHERION_QT_HOST
)

# GLEW for non-Apple platforms, because MacOS mainly uses their own proprietary framework for OpenGL.
if(NOT APPLE)
    target_link_libraries(blackhole-sim PRIVATE GLEW::GLEW)
endif()

if(glm_FOUND)
    target_link_libraries(blackhole-sim PRIVATE glm::glm)
else()
    target_include_directories(blackhole-sim PRIVATE ${GLM_INCLUDE_DIRS})
endif()

# Link OpenGL on macOS (also silences Apple deprecation warnings via platform.hpp define)
if(APPLE)
    find_library(OpenGL_LIBRARY OpenGL)
    target_link_libraries(blackhole-sim PRIVATE ${OpenGL_LIBRARY})
endif()

# GCC < 9 ships <filesystem> in a separate library; warn-proof the build.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
    target_link_libraries(blackhole-sim PRIVATE stdc++fs)
endif()

# ── Windows: deploy Qt platform plugin ────────────────────────────────────
# Qt apps fail at startup with "could not find or load the Qt platform plugin
# 'windows'" unless qwindows.dll sits in a `platforms/` subfolder next to the
# exe. vcpkg copies the core Qt6*.dll files automatically but does NOT copy
# plugins. Locate qwindows.dll from the imported Qt6::QWindowsIntegrationPlugin
# target (preferred) or fall back to probing the Qt6_DIR install tree.
if(WIN32)
    set(_qwindows_dll "")
    if(TARGET Qt6::QWindowsIntegrationPlugin)
        set(_qwindows_dll "$<TARGET_FILE:Qt6::QWindowsIntegrationPlugin>")
    elseif(DEFINED Qt6_DIR)
        # Qt6_DIR is typically <prefix>/share/Qt6 or <prefix>/lib/cmake/Qt6.
        # Walk up to the install prefix and look in the standard plugin path.
        get_filename_component(_qt_prefix "${Qt6_DIR}/../../.." ABSOLUTE)
        foreach(_candidate
                "${_qt_prefix}/Qt6/plugins/platforms/qwindows.dll"
                "${_qt_prefix}/plugins/platforms/qwindows.dll"
                "${_qt_prefix}/share/qt6/plugins/platforms/qwindows.dll")
            if(EXISTS "${_candidate}")
                set(_qwindows_dll "${_candidate}")
                break()
            endif()
        endforeach()
    endif()
    if(_qwindows_dll)
        add_custom_command(TARGET blackhole-sim POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E make_directory
                    $<TARGET_FILE_DIR:blackhole-sim>/platforms
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${_qwindows_dll}"
                    $<TARGET_FILE_DIR:blackhole-sim>/platforms/
            COMMENT "Deploying Qt platform plugin qwindows.dll")
    else()
        message(WARNING
            "Could not locate qwindows.dll. The app will fail at startup with "
            "'could not find or load the Qt platform plugin \"windows\"'. "
            "Manually copy <qt-prefix>/plugins/platforms/qwindows.dll into "
            "$<TARGET_FILE_DIR:blackhole-sim>/platforms/.")
    endif()
endif()

# Copy shader files for 3D simulation to build output
add_custom_command(TARGET blackhole-sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/src/3D/BlackHole3D.frag
            $<TARGET_FILE_DIR:blackhole-sim>/BlackHole3D.frag)

add_custom_command(TARGET blackhole-sim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_SOURCE_DIR}/src/3D/BlackHole3D_PhotorealDisk.frag
            $<TARGET_FILE_DIR:blackhole-sim>/BlackHole3D_PhotorealDisk.frag)

# ── Standalone simulation executables (restored) ───────────────────────────
add_executable(blackhole-2D src/2D/BlackHole2D.cpp)
target_include_directories(blackhole-2D PRIVATE
    ${CMAKE_SOURCE_DIR}/src/2D
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/external/json   # future 2D catalog adapter
)
target_link_libraries(blackhole-2D PRIVATE SFML::Graphics SFML::Window SFML::System)

# ── Dear ImGui static library ─────────────────────────────────────────────
# Contains imgui core + imgui-SFML (standalone path) + imgui_impl_opengl3
# (Qt-embedded path). Both backends coexist: the SFML backend is used by
# blackhole-3D; the raw GL3 backend is used by blackhole-sim (AETHERION_QT_HOST).
add_library(aetherion_imgui STATIC
    external/imgui/imgui.cpp
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui-sfml/imgui-SFML.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)
# These third-party sources don't use Qt; skip the meta-object compiler.
set_target_properties(aetherion_imgui PROPERTIES
    AUTOMOC OFF AUTORCC OFF AUTOUIC OFF
)
target_include_directories(aetherion_imgui PUBLIC
    ${CMAKE_SOURCE_DIR}/external/imgui
    ${CMAKE_SOURCE_DIR}/external/imgui-sfml
    ${CMAKE_SOURCE_DIR}/external/imgui/backends
)
target_compile_definitions(aetherion_imgui PUBLIC
    IMGUI_USER_CONFIG="imconfig-SFML.h"
)
target_link_libraries(aetherion_imgui PUBLIC
    SFML::Graphics SFML::Window SFML::System OpenGL::GL
)
if(NOT APPLE)
    target_link_libraries(aetherion_imgui PUBLIC GLEW::GLEW)
endif()

add_executable(blackhole-3D
    src/3D/BlackHole3D.cpp
    src/3D/shader_utils.cpp
    src/3D/texture_utils.cpp
)
target_include_directories(blackhole-3D PRIVATE
    ${CMAKE_SOURCE_DIR}/src/3D
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/external/json
    ${CMAKE_SOURCE_DIR}
)
target_link_libraries(blackhole-3D PRIVATE
    SFML::Graphics SFML::Window SFML::System OpenGL::GL
    aetherion_imgui
)
if(NOT APPLE)
    target_link_libraries(blackhole-3D PRIVATE GLEW::GLEW)
endif()
if(glm_FOUND)
    target_link_libraries(blackhole-3D PRIVATE glm::glm)
else()
    target_include_directories(blackhole-3D PRIVATE ${GLM_INCLUDE_DIRS})
endif()

# ── Shader copy: file-driven so any .frag edit triggers a re-copy ────────────
# Using add_custom_command with OUTPUT means cmake tracks the .frag as a real
# dependency, editing the source file causes a rebuild of the copy target,
# regardless of whether any C++ source changed.
foreach(_frag BlackHole3D.frag BlackHole3D_PhotorealDisk.frag)
    add_custom_command(
        OUTPUT  ${CMAKE_BINARY_DIR}/${_frag}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/src/3D/${_frag}
                ${CMAKE_BINARY_DIR}/${_frag}
        DEPENDS ${CMAKE_SOURCE_DIR}/src/3D/${_frag}
        COMMENT "Syncing shader ${_frag} to build directory"
    )
endforeach()
add_custom_target(sync_shaders_3D ALL
    DEPENDS
        ${CMAKE_BINARY_DIR}/BlackHole3D.frag
        ${CMAKE_BINARY_DIR}/BlackHole3D_PhotorealDisk.frag
)
add_dependencies(blackhole-3D sync_shaders_3D)

# If a disk texture exists, copy it next to the executable for convenience.
if(EXISTS ${CMAKE_SOURCE_DIR}/src/disk_texture.png)
    add_custom_command(TARGET blackhole-sim POST_BUILD # we live in a cruel world, a world where I have to write fifty lines of CMake to copy a single file
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/src/disk_texture.png
                $<TARGET_FILE_DIR:blackhole-sim>/disk_texture.png)
    add_custom_command(TARGET blackhole-3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/src/disk_texture.png
                $<TARGET_FILE_DIR:blackhole-3D>/disk_texture.png)
endif()

# Copy the black hole JSON catalog next to each executable so ResourceManager's
# findDir("catalog") resolves it in dev/run-from-build-dir workflows.
foreach(_bh_target blackhole-sim blackhole-3D)
    add_custom_command(TARGET ${_bh_target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/src/common/resources/catalog
                $<TARGET_FILE_DIR:${_bh_target}>/catalog)
endforeach()

# macOS: also land the catalog inside the .app bundle Resources (mirrors the
# icon copy). The bundle's Resources dir is on the ResourceManager search path in this case
if(APPLE)
    add_custom_command(TARGET blackhole-sim POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_SOURCE_DIR}/src/common/resources/catalog
                ${CMAKE_BINARY_DIR}/Aetherion.app/Contents/Resources/catalog)
endif()

# Copy app icon next to the executable (sidebar logo display at runtime).
if(EXISTS ${CMAKE_SOURCE_DIR}/icon.png)
    add_custom_command(TARGET blackhole-sim POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/icon.png
                $<TARGET_FILE_DIR:blackhole-sim>/icon.png)
    if(APPLE)
        add_custom_command(TARGET blackhole-sim POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${CMAKE_SOURCE_DIR}/icon.png
                    ${CMAKE_BINARY_DIR}/Aetherion.app/Contents/Resources/icon.png)
    endif()
endif()
# Pre-baked transparent version of the icon (preferred by the sidebar logo loader).
if(EXISTS ${CMAKE_SOURCE_DIR}/icon_transparent.png)
    add_custom_command(TARGET blackhole-sim POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/icon_transparent.png
                $<TARGET_FILE_DIR:blackhole-sim>/icon_transparent.png)
    if(APPLE)
        add_custom_command(TARGET blackhole-sim POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${CMAKE_SOURCE_DIR}/icon_transparent.png
                    ${CMAKE_BINARY_DIR}/Aetherion.app/Contents/Resources/icon_transparent.png)
    endif()
endif()

# Background texture (optional) is often kept next to the executable in build/.
if(EXISTS ${CMAKE_SOURCE_DIR}/build/background.png)
    add_custom_command(TARGET blackhole-sim POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/build/background.png
                $<TARGET_FILE_DIR:blackhole-sim>/background.png)
    add_custom_command(TARGET blackhole-3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/build/background.png
                $<TARGET_FILE_DIR:blackhole-3D>/background.png)
elseif(EXISTS ${CMAKE_SOURCE_DIR}/background.png)
    add_custom_command(TARGET blackhole-sim POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/background.png
                $<TARGET_FILE_DIR:blackhole-sim>/background.png)
    add_custom_command(TARGET blackhole-3D POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_SOURCE_DIR}/background.png
                $<TARGET_FILE_DIR:blackhole-3D>/background.png)
endif()

# ── Install targets ───────────────────────────────────────────
if(APPLE)
    install(TARGETS blackhole-sim blackhole-2D blackhole-3D RUNTIME DESTINATION bin BUNDLE DESTINATION .)
else()
    install(TARGETS blackhole-sim blackhole-2D blackhole-3D RUNTIME DESTINATION bin)
endif()

# Runtime data files (shaders, textures) loaded via ResourceManager::find().
# Installed to share/aetherionsuite/ so XDG_DATA_DIRS search (Flatpak/system) finds them.
install(FILES
    src/3D/BlackHole3D.frag
    src/3D/BlackHole3D_PhotorealDisk.frag
    DESTINATION share/aetherionsuite)

if(EXISTS ${CMAKE_SOURCE_DIR}/src/disk_texture.png)
    install(FILES src/disk_texture.png DESTINATION share/aetherionsuite)
endif()

# Black hole JSON catalog is at share/aetherionsuite/catalog (covered by the XDG
# search paths ResourceManager already probes on Linux/Flatpak).
install(DIRECTORY src/common/resources/catalog DESTINATION share/aetherionsuite)
if(EXISTS ${CMAKE_SOURCE_DIR}/build/background.png)
    install(FILES build/background.png DESTINATION share/aetherionsuite)
elseif(EXISTS ${CMAKE_SOURCE_DIR}/background.png)
    install(FILES background.png DESTINATION share/aetherionsuite)
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/icon.png)
    install(FILES icon.png DESTINATION share/aetherionsuite)
    install(FILES icon.png
        DESTINATION share/icons/hicolor/256x256/apps
        RENAME io.github.0xLiam0920.AetherionSuite.png)
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/icon_transparent.png)
    install(FILES icon_transparent.png DESTINATION share/aetherionsuite)
endif()

# Desktop entry and AppStream metadata (used by Flatpak and system package managers).
if(EXISTS ${CMAKE_SOURCE_DIR}/flatpak/io.github.0xLiam0920.AetherionSuite.desktop)
    install(FILES flatpak/io.github.0xLiam0920.AetherionSuite.desktop
        DESTINATION share/applications)
endif()
if(EXISTS ${CMAKE_SOURCE_DIR}/flatpak/io.github.0xLiam0920.AetherionSuite.metainfo.xml)
    install(FILES flatpak/io.github.0xLiam0920.AetherionSuite.metainfo.xml
        DESTINATION share/metainfo)
endif()

# ── macOS .app bundle target ──────────────────────────────────
if(APPLE)
    add_custom_target(bundle
        COMMAND bash ${CMAKE_SOURCE_DIR}/bundle_app.sh
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        DEPENDS blackhole-sim blackhole-2D blackhole-3D
        COMMENT "Assembling Aetherion.app bundle"
    )
endif()

# ── Physics regression tests (beta validation) ──────────────────────────────
enable_testing()

add_executable(physics-regression-tests
    tests/physics_regression_tests.cpp
)
target_include_directories(physics-regression-tests PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/2D
    ${CMAKE_SOURCE_DIR}/external/json
)

# Enable the JSON-catalog round-trip test when glm is available (it needs
# config.hpp). The test header-guards itself on <glm/glm.hpp> + <nlohmann/json.hpp>
# so the lightweight CI compile (no glm) simply skips it.
if(glm_FOUND)
    target_link_libraries(physics-regression-tests PRIVATE glm::glm)
elseif(DEFINED GLM_INCLUDE_DIRS)
    target_include_directories(physics-regression-tests PRIVATE ${GLM_INCLUDE_DIRS})
endif()
target_compile_definitions(physics-regression-tests PRIVATE
    AETHERION_CATALOG_DIR="${CMAKE_SOURCE_DIR}/src/common/resources/catalog")

add_test(
    NAME physics-regression
    COMMAND physics-regression-tests
)
