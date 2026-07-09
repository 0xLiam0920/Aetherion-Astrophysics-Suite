// ============================================================
// bh3d_core.hpp
// ============================================================
// The shared core of the 3D black-hole simulation.
//
// This header consolidates everything that used to be duplicated between
//   src/3D/BlackHole3D.cpp                    (standalone window)
//   src/QT-LAUNCHER/simulation_3d_widget.cpp  (Qt-embedded widget)
//
// Both call-sites now keep only their platform glue (window creation,
// event source, viewport size lookup) and delegate the actual simulation,
// rendering, HUD, and input handling here.
//
// The module is header-only with `inline` storage to avoid touching CMake.
// I might wanna split this into multiple files later, it's a bit monolithic.
// ============================================================
#pragma once

#include "platform.hpp"
#include <SFML/Window/Keyboard.hpp>
#include <SFML/Window/Mouse.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Graphics/Font.hpp>
#include <glm/glm.hpp>

#include "gl_types.hpp"
#include <cstdint>
#include <cmath>
#include "config.hpp"
#include "resource_manager.hpp"
#include "bh3d_shadersrcs.hpp"
#include "bh3d_shaderutilities.hpp"
#include "bh3d_textureutils.hpp"
#include "bh3d_bloompipeline.hpp"
#include "bh3d_camera.hpp"
#include "bh3d_input.hpp"
#include "bh3d_simulationstates.hpp"
#include "bh3d_presets.hpp"
#include "orbital_body.hpp"
#include "gl_font.hpp"
#include "hud_panel.hpp"
#include "physics_overlay.hpp"

#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace bh3d { // TODO: this header has grown large; consider splitting it into
                // separate files as the core functionality keeps expanding.
/*
Planned continuation of split:
- bh3d_core.hpp: Core simulation state, physics, and rendering logic (this file)
-bh3d_mergers.hpp: Merger logic, states, and sub rendering.

Also need to add more of the logic to other files


*/

// ────────────────────────────────────────────────────────────
//  Small string helpers
// ────────────────────────────────────────────────────────────
inline std::string trimCopy(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

inline std::string upperCopy(std::string s) {
    for (char& c : s) c = char(std::toupper((unsigned char)c));
    return s;
}

inline const char* keyToString(sf::Keyboard::Key k) {
    switch (k) {
        case sf::Keyboard::Key::A: return "A"; case sf::Keyboard::Key::B: return "B";
        case sf::Keyboard::Key::C: return "C"; case sf::Keyboard::Key::D: return "D";
        case sf::Keyboard::Key::E: return "E"; case sf::Keyboard::Key::F: return "F";
        case sf::Keyboard::Key::G: return "G"; case sf::Keyboard::Key::H: return "H";
        case sf::Keyboard::Key::I: return "I"; case sf::Keyboard::Key::J: return "J";
        case sf::Keyboard::Key::K: return "K"; case sf::Keyboard::Key::L: return "L";
        case sf::Keyboard::Key::M: return "M"; case sf::Keyboard::Key::N: return "N";
        case sf::Keyboard::Key::O: return "O"; case sf::Keyboard::Key::P: return "P";
        case sf::Keyboard::Key::Q: return "Q"; case sf::Keyboard::Key::R: return "R";
        case sf::Keyboard::Key::S: return "S"; case sf::Keyboard::Key::T: return "T";
        case sf::Keyboard::Key::U: return "U"; case sf::Keyboard::Key::V: return "V";
        case sf::Keyboard::Key::W: return "W"; case sf::Keyboard::Key::X: return "X";
        case sf::Keyboard::Key::Y: return "Y"; case sf::Keyboard::Key::Z: return "Z";
        case sf::Keyboard::Key::Num0: return "NUM0"; case sf::Keyboard::Key::Num1: return "NUM1";
        case sf::Keyboard::Key::Num2: return "NUM2"; case sf::Keyboard::Key::Num3: return "NUM3";
        case sf::Keyboard::Key::Num4: return "NUM4"; case sf::Keyboard::Key::Num5: return "NUM5";
        case sf::Keyboard::Key::Num6: return "NUM6"; case sf::Keyboard::Key::Num7: return "NUM7";
        case sf::Keyboard::Key::Num8: return "NUM8"; case sf::Keyboard::Key::Num9: return "NUM9";
        case sf::Keyboard::Key::Space:    return "SPACE";
        case sf::Keyboard::Key::LControl: return "LCTRL";
        case sf::Keyboard::Key::LShift:   return "LSHIFT";
        case sf::Keyboard::Key::Escape:   return "ESCAPE";
        case sf::Keyboard::Key::Equal:    return "=";
        case sf::Keyboard::Key::Hyphen:   return "-";
        default: return "UNKNOWN";
    }
}

inline bool keyFromString(const std::string& raw, sf::Keyboard::Key& out) {
    const std::string s = upperCopy(trimCopy(raw));
    if (s.empty()) return false;
    if (s=="A"){out=sf::Keyboard::Key::A;return true;} if (s=="B"){out=sf::Keyboard::Key::B;return true;}
    if (s=="C"){out=sf::Keyboard::Key::C;return true;} if (s=="D"){out=sf::Keyboard::Key::D;return true;}
    if (s=="E"){out=sf::Keyboard::Key::E;return true;} if (s=="F"){out=sf::Keyboard::Key::F;return true;}
    if (s=="G"){out=sf::Keyboard::Key::G;return true;} if (s=="H"){out=sf::Keyboard::Key::H;return true;}
    if (s=="I"){out=sf::Keyboard::Key::I;return true;} if (s=="J"){out=sf::Keyboard::Key::J;return true;}
    if (s=="K"){out=sf::Keyboard::Key::K;return true;} if (s=="L"){out=sf::Keyboard::Key::L;return true;}
    if (s=="M"){out=sf::Keyboard::Key::M;return true;} if (s=="N"){out=sf::Keyboard::Key::N;return true;}
    if (s=="O"){out=sf::Keyboard::Key::O;return true;} if (s=="P"){out=sf::Keyboard::Key::P;return true;}
    if (s=="Q"){out=sf::Keyboard::Key::Q;return true;} if (s=="R"){out=sf::Keyboard::Key::R;return true;}
    if (s=="S"){out=sf::Keyboard::Key::S;return true;} if (s=="T"){out=sf::Keyboard::Key::T;return true;}
    if (s=="U"){out=sf::Keyboard::Key::U;return true;} if (s=="V"){out=sf::Keyboard::Key::V;return true;}
    if (s=="W"){out=sf::Keyboard::Key::W;return true;} if (s=="X"){out=sf::Keyboard::Key::X;return true;}
    if (s=="Y"){out=sf::Keyboard::Key::Y;return true;} if (s=="Z"){out=sf::Keyboard::Key::Z;return true;}
    if (s=="NUM0"||s=="0"){out=sf::Keyboard::Key::Num0;return true;}
    if (s=="NUM1"||s=="1"){out=sf::Keyboard::Key::Num1;return true;}
    if (s=="NUM2"||s=="2"){out=sf::Keyboard::Key::Num2;return true;}
    if (s=="NUM3"||s=="3"){out=sf::Keyboard::Key::Num3;return true;}
    if (s=="NUM4"||s=="4"){out=sf::Keyboard::Key::Num4;return true;}
    if (s=="NUM5"||s=="5"){out=sf::Keyboard::Key::Num5;return true;}
    if (s=="NUM6"||s=="6"){out=sf::Keyboard::Key::Num6;return true;}
    if (s=="NUM7"||s=="7"){out=sf::Keyboard::Key::Num7;return true;}
    if (s=="NUM8"||s=="8"){out=sf::Keyboard::Key::Num8;return true;}
    if (s=="NUM9"||s=="9"){out=sf::Keyboard::Key::Num9;return true;}
    if (s=="SPACE"){out=sf::Keyboard::Key::Space;return true;}
    if (s=="LCTRL"||s=="LCONTROL"||s=="CTRL"){out=sf::Keyboard::Key::LControl;return true;}
    if (s=="LSHIFT"||s=="SHIFT"){out=sf::Keyboard::Key::LShift;return true;}
    if (s=="ESC"||s=="ESCAPE"){out=sf::Keyboard::Key::Escape;return true;}
    if (s=="="||s=="EQUAL"||s=="PLUS"||s=="+"){out=sf::Keyboard::Key::Equal;return true;}
    if (s=="-"||s=="HYPHEN"||s=="MINUS"){out=sf::Keyboard::Key::Hyphen;return true;}
    return false;
}

// ────────────────────────────────────────────────────────────
//  Action keybinds (editable via on-disk config)
// ────────────────────────────────────────────────────────────
struct ActionKeybinds {
    sf::Keyboard::Key toggleFreelook   = sf::Keyboard::Key::F;
    sf::Keyboard::Key toggleJets       = sf::Keyboard::Key::J;
    sf::Keyboard::Key toggleBLR        = sf::Keyboard::Key::G;
    sf::Keyboard::Key toggleOrbBody    = sf::Keyboard::Key::O;
    sf::Keyboard::Key toggleDoppler    = sf::Keyboard::Key::V;
    sf::Keyboard::Key toggleBlueshift  = sf::Keyboard::Key::U;
    sf::Keyboard::Key toggleHostGalaxy = sf::Keyboard::Key::Num1;
    sf::Keyboard::Key toggleLAB        = sf::Keyboard::Key::Num2;
    sf::Keyboard::Key toggleCGM        = sf::Keyboard::Key::Num3;
    sf::Keyboard::Key cycleAnimSpeed   = sf::Keyboard::Key::Y;
    sf::Keyboard::Key toggleHUD        = sf::Keyboard::Key::H;
    sf::Keyboard::Key toggleDebugHUD   = sf::Keyboard::Key::B;
    sf::Keyboard::Key resetTilt        = sf::Keyboard::Key::R;
    sf::Keyboard::Key nextProfile      = sf::Keyboard::Key::N;
    sf::Keyboard::Key toggleRender     = sf::Keyboard::Key::P;
    sf::Keyboard::Key releaseMouse     = sf::Keyboard::Key::Escape;
    sf::Keyboard::Key toggleOverlays   = sf::Keyboard::Key::M;
    sf::Keyboard::Key speedUp          = sf::Keyboard::Key::Equal;   // '+' / '='
    sf::Keyboard::Key speedDown        = sf::Keyboard::Key::Hyphen;  // '-'
    sf::Keyboard::Key toggleRK4Orbits  = sf::Keyboard::Key::K;
    sf::Keyboard::Key toggleRK4Photons = sf::Keyboard::Key::L;
    sf::Keyboard::Key toggleSpacetime  = sf::Keyboard::Key::T;
    sf::Keyboard::Key openMergerMenu   = sf::Keyboard::Key::C;
};

inline std::filesystem::path keybindConfigPath() {
    return std::filesystem::path(platformUserDataDir()) / "blackhole3d_keybinds.cfg";
}

inline void enforceKeybindConflicts(ActionKeybinds& keys) {
    const ActionKeybinds defaults{};
    auto isMovementKey = [](sf::Keyboard::Key k) {
        return k == sf::Keyboard::Key::W  || k == sf::Keyboard::Key::A  ||
               k == sf::Keyboard::Key::S  || k == sf::Keyboard::Key::D  ||
               k == sf::Keyboard::Key::Q  || k == sf::Keyboard::Key::E  ||
               k == sf::Keyboard::Key::Space    ||
               k == sf::Keyboard::Key::LControl ||
               k == sf::Keyboard::Key::LShift;
    };
    struct Entry { const char* name; sf::Keyboard::Key* key; sf::Keyboard::Key def; };
    std::array<Entry, 23> entries = {{
        {"toggle_freelook",   &keys.toggleFreelook,   defaults.toggleFreelook},
        {"toggle_jets",       &keys.toggleJets,       defaults.toggleJets},
        {"toggle_blr",        &keys.toggleBLR,        defaults.toggleBLR},
        {"toggle_orb_body",   &keys.toggleOrbBody,    defaults.toggleOrbBody},
        {"toggle_doppler",    &keys.toggleDoppler,    defaults.toggleDoppler},
        {"toggle_blueshift",  &keys.toggleBlueshift,  defaults.toggleBlueshift},
        {"toggle_host_galaxy",&keys.toggleHostGalaxy, defaults.toggleHostGalaxy},
        {"toggle_lab",        &keys.toggleLAB,        defaults.toggleLAB},
        {"toggle_cgm",        &keys.toggleCGM,        defaults.toggleCGM},
        {"cycle_anim_speed",  &keys.cycleAnimSpeed,   defaults.cycleAnimSpeed},
        {"toggle_hud",        &keys.toggleHUD,        defaults.toggleHUD},
        {"toggle_debug_hud",  &keys.toggleDebugHUD,   defaults.toggleDebugHUD},
        {"reset_tilt",        &keys.resetTilt,        defaults.resetTilt},
        {"next_profile",      &keys.nextProfile,      defaults.nextProfile},
        {"toggle_render",     &keys.toggleRender,     defaults.toggleRender},
        {"release_mouse",     &keys.releaseMouse,     defaults.releaseMouse},
        {"toggle_overlays",   &keys.toggleOverlays,   defaults.toggleOverlays},
        {"speed_up",          &keys.speedUp,          defaults.speedUp},
        {"speed_down",        &keys.speedDown,        defaults.speedDown},
        {"toggle_rk4_orbits", &keys.toggleRK4Orbits,  defaults.toggleRK4Orbits},
        {"toggle_rk4_photons",&keys.toggleRK4Photons, defaults.toggleRK4Photons},
        {"toggle_spacetime",  &keys.toggleSpacetime,  defaults.toggleSpacetime},
        {"open_merger_menu",  &keys.openMergerMenu,   defaults.openMergerMenu}
    }};
    for (auto& e : entries) {
        if (isMovementKey(*e.key)) *e.key = e.def;
    }
    std::unordered_map<int, const char*> seen;
    for (auto& e : entries) {
        const int code = static_cast<int>(*e.key);
        auto it = seen.find(code);
        if (it != seen.end()) *e.key = e.def;
        seen[static_cast<int>(*e.key)] = e.name;
    }
}

inline void writeDefaultKeybindFile(const std::filesystem::path& path, const ActionKeybinds& k) {
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) return;
    out << "# Aetherion 3D keybinds\n# Format: action=KEY\n\n";
    out << "toggle_freelook="   << keyToString(k.toggleFreelook)   << "\n";
    out << "toggle_jets="       << keyToString(k.toggleJets)       << "\n";
    out << "toggle_blr="        << keyToString(k.toggleBLR)        << "\n";
    out << "toggle_orb_body="   << keyToString(k.toggleOrbBody)    << "\n";
    out << "toggle_doppler="    << keyToString(k.toggleDoppler)    << "\n";
    out << "toggle_blueshift="  << keyToString(k.toggleBlueshift)  << "\n";
    out << "toggle_host_galaxy="<< keyToString(k.toggleHostGalaxy) << "\n";
    out << "toggle_lab="        << keyToString(k.toggleLAB)        << "\n";
    out << "toggle_cgm="        << keyToString(k.toggleCGM)        << "\n";
    out << "cycle_anim_speed="  << keyToString(k.cycleAnimSpeed)   << "\n";
    out << "toggle_hud="        << keyToString(k.toggleHUD)        << "\n";
    out << "toggle_debug_hud="  << keyToString(k.toggleDebugHUD)   << "\n";
    out << "reset_tilt="        << keyToString(k.resetTilt)        << "\n";
    out << "next_profile="      << keyToString(k.nextProfile)      << "\n";
    out << "toggle_render="     << keyToString(k.toggleRender)     << "\n";
    out << "release_mouse="     << keyToString(k.releaseMouse)     << "\n";
    out << "toggle_overlays="   << keyToString(k.toggleOverlays)   << "\n";
    out << "speed_up="          << keyToString(k.speedUp)          << "\n";
    out << "speed_down="        << keyToString(k.speedDown)        << "\n";
    out << "toggle_rk4_orbits=" << keyToString(k.toggleRK4Orbits)  << "\n";
    out << "toggle_rk4_photons="<< keyToString(k.toggleRK4Photons) << "\n";
    out << "toggle_spacetime="  << keyToString(k.toggleSpacetime)  << "\n";
    out << "open_merger_menu="  << keyToString(k.openMergerMenu)   << "\n";
}

inline ActionKeybinds loadActionKeybinds() {
    ActionKeybinds keys;
    const auto cfgPath = keybindConfigPath();
    if (!std::filesystem::exists(cfgPath)) {
        writeDefaultKeybindFile(cfgPath, keys);
        return keys;
    }
    std::ifstream in(cfgPath);
    if (!in) return keys;
    std::unordered_map<std::string, sf::Keyboard::Key*> targets = {
        {"toggle_freelook",   &keys.toggleFreelook},
        {"toggle_jets",       &keys.toggleJets},
        {"toggle_blr",        &keys.toggleBLR},
        {"toggle_orb_body",   &keys.toggleOrbBody},
        {"toggle_doppler",    &keys.toggleDoppler},
        {"toggle_blueshift",  &keys.toggleBlueshift},
        {"toggle_host_galaxy",&keys.toggleHostGalaxy},
        {"toggle_lab",        &keys.toggleLAB},
        {"toggle_cgm",        &keys.toggleCGM},
        {"cycle_anim_speed",  &keys.cycleAnimSpeed},
        {"toggle_hud",        &keys.toggleHUD},
        {"toggle_debug_hud",  &keys.toggleDebugHUD},
        {"reset_tilt",        &keys.resetTilt},
        {"next_profile",      &keys.nextProfile},
        {"toggle_render",     &keys.toggleRender},
        {"release_mouse",     &keys.releaseMouse},
        {"toggle_overlays",   &keys.toggleOverlays},
        {"speed_up",          &keys.speedUp},
        {"speed_down",        &keys.speedDown},
        {"toggle_rk4_orbits", &keys.toggleRK4Orbits},
        {"toggle_rk4_photons",&keys.toggleRK4Photons},
        {"toggle_spacetime",  &keys.toggleSpacetime},
        {"open_merger_menu",  &keys.openMergerMenu}
    };
    std::string line;
    while (std::getline(in, line)) {
        const auto trimmed = trimCopy(line);
        if (trimmed.empty() || trimmed[0] == '#') continue;
        const auto eqPos = trimmed.find('=');
        if (eqPos == std::string::npos) continue;
        const std::string key   = trimCopy(trimmed.substr(0, eqPos));
        const std::string value = trimCopy(trimmed.substr(eqPos + 1));
        auto it = targets.find(key);
        if (it == targets.end()) continue;
        sf::Keyboard::Key parsed{};
        if (keyFromString(value, parsed)) *it->second = parsed;
    }
    enforceKeybindConflicts(keys);
    return keys;
}

// ────────────────────────────────────────────────────────────
//  Scene shader uniform locations
// ────────────────────────────────────────────────────────────
struct SceneUniforms {
    GLint resolution, cameraPos, cameraDir, cameraUp, fov;
    GLint blackHoleRadius, blackHolePos, backgroundTex, diskTex;
    GLint diskInnerRadius, diskOuterRadius, diskHalfThickness;
    GLint spinParameter;
    GLint jetRadius, jetLength, jetColor;
    GLint showJets, showBLR, showDoppler, showBlueshift;
    GLint blrInnerRadius, blrOuterRadius, blrThickness, blrStrength;
    GLint uTime;
    GLint showOrbBody;
    GLint orbBodyPos[10];
    GLint orbBodyRadius[10];
    GLint orbBodyColor[10];
    GLint orbBodyType[10];
    GLint numOrbBodies;
    GLint showHostGalaxy, hostGalaxyRadius, hostGalaxyColor;
    GLint showLAB, labPos, labRadius, labColor;
    GLint showCGM, cgmRadius, cgmColor;
    GLint maxStepsOverride;
    GLint diskPeakTemp, diskDisplayTempInner, diskDisplayTempOuter;
    GLint diskSatBoostInner, diskSatBoostOuter;

    // Inspiralling merger BH identities
    GLint secBHActive, secBHDiskNormal, secBHSpin;
    GLint secBHDiskInner, secBHDiskOuter, secBHDiskStrength;
    GLint secBHColorInner, secBHColorOuter;
    GLint secBHShowJets, secBHJetColor, secBHJetRadius, secBHJetLength;

    static SceneUniforms lookup(GLProgram& prog) {
        SceneUniforms u;
        u.resolution           = prog.uniform("resolution");
        u.cameraPos            = prog.uniform("cameraPos");
        u.cameraDir            = prog.uniform("cameraDir");
        u.cameraUp             = prog.uniform("cameraUp");
        u.fov                  = prog.uniform("fov");
        u.blackHoleRadius      = prog.uniform("blackHoleRadius");
        u.blackHolePos         = prog.uniform("blackHolePos");
        u.backgroundTex        = prog.uniform("backgroundTex");
        u.diskTex              = prog.uniform("diskTex");
        u.diskInnerRadius      = prog.uniform("diskInnerRadius");
        u.diskOuterRadius      = prog.uniform("diskOuterRadius");
        u.diskHalfThickness    = prog.uniform("diskHalfThickness");
        u.spinParameter        = prog.uniform("spinParameter");
        u.jetRadius            = prog.uniform("jetRadius");
        u.jetLength            = prog.uniform("jetLength");
        u.jetColor             = prog.uniform("jetColor");
        u.showJets             = prog.uniform("showJets");
        u.showBLR              = prog.uniform("showBLR");
        u.showDoppler          = prog.uniform("showDoppler");
        u.showBlueshift        = prog.uniform("showBlueshift");
        u.blrInnerRadius       = prog.uniform("blrInnerRadius");
        u.blrOuterRadius       = prog.uniform("blrOuterRadius");
        u.blrThickness         = prog.uniform("blrThickness");
        u.blrStrength          = prog.uniform("blrStrength");
        u.uTime                = prog.uniform("uTime");
        u.showOrbBody          = prog.uniform("showOrbBody");
        for (int i = 0; i < 10; ++i) {
            char buf[64];
            snprintf(buf, sizeof(buf), "orbBodyPos[%d]",    i); u.orbBodyPos[i]    = prog.uniform(buf);
            snprintf(buf, sizeof(buf), "orbBodyRadius[%d]", i); u.orbBodyRadius[i] = prog.uniform(buf);
            snprintf(buf, sizeof(buf), "orbBodyColor[%d]",  i); u.orbBodyColor[i]  = prog.uniform(buf);
            snprintf(buf, sizeof(buf), "orbBodyType[%d]",   i); u.orbBodyType[i]   = prog.uniform(buf);
        }
        u.numOrbBodies         = prog.uniform("numOrbBodies");
        u.showHostGalaxy       = prog.uniform("showHostGalaxy");
        u.hostGalaxyRadius     = prog.uniform("hostGalaxyRadius");
        u.hostGalaxyColor      = prog.uniform("hostGalaxyColor");
        u.showLAB              = prog.uniform("showLAB");
        u.labPos               = prog.uniform("labPos");
        u.labRadius            = prog.uniform("labRadius");
        u.labColor             = prog.uniform("labColor");
        u.showCGM              = prog.uniform("showCGM");
        u.cgmRadius            = prog.uniform("cgmRadius");
        u.cgmColor             = prog.uniform("cgmColor");
        u.maxStepsOverride     = prog.uniform("maxStepsOverride");
        u.diskPeakTemp         = prog.uniform("diskPeakTemp");
        u.diskDisplayTempInner = prog.uniform("diskDisplayTempInner");
        u.diskDisplayTempOuter = prog.uniform("diskDisplayTempOuter");
        u.diskSatBoostInner    = prog.uniform("diskSatBoostInner");
        u.diskSatBoostOuter    = prog.uniform("diskSatBoostOuter");
        u.secBHActive          = prog.uniform("secBHActive");
        u.secBHDiskNormal      = prog.uniform("secBHDiskNormal");
        u.secBHSpin            = prog.uniform("secBHSpin");
        u.secBHDiskInner       = prog.uniform("secBHDiskInner");
        u.secBHDiskOuter       = prog.uniform("secBHDiskOuter");
        u.secBHDiskStrength    = prog.uniform("secBHDiskStrength");
        u.secBHColorInner      = prog.uniform("secBHColorInner");
        u.secBHColorOuter      = prog.uniform("secBHColorOuter");
        u.secBHShowJets        = prog.uniform("secBHShowJets");
        u.secBHJetColor        = prog.uniform("secBHJetColor");
        u.secBHJetRadius       = prog.uniform("secBHJetRadius");
        u.secBHJetLength       = prog.uniform("secBHJetLength");
        return u;
    }
};

inline void setSceneUniforms(GLProgram& prog, const SceneUniforms& u,
                             const PhysicsSnapshot& snap,
                             GLuint bgTexId, GLuint diskTexId)
{
    prog.use();
    glUniform2f(u.resolution,     float(snap.windowW), float(snap.windowH));
    glUniform3f(u.cameraPos,      snap.cameraPos.x, snap.cameraPos.y, snap.cameraPos.z);
    glUniform3f(u.cameraDir,      snap.cameraDir.x, snap.cameraDir.y, snap.cameraDir.z);
    glUniform3f(u.cameraUp,       snap.cameraUp.x,  snap.cameraUp.y,  snap.cameraUp.z);
    glUniform1f(u.fov,            snap.fov);
    glUniform1f(u.blackHoleRadius, snap.bhRadius);
    glUniform3f(u.blackHolePos,   snap.bhPosition.x, snap.bhPosition.y, snap.bhPosition.z);

    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, bgTexId);   glUniform1i(u.backgroundTex, 0);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, diskTexId); glUniform1i(u.diskTex, 1);

    glUniform1f(u.diskInnerRadius,       snap.diskInnerRadius);
    glUniform1f(u.diskOuterRadius,       snap.diskOuterRadius);
    glUniform1f(u.diskHalfThickness,     snap.diskHalfThickness);
    glUniform1f(u.diskPeakTemp,          snap.diskPeakTemp);
    glUniform1f(u.diskDisplayTempInner,  snap.diskDisplayTempInner);
    glUniform1f(u.diskDisplayTempOuter,  snap.diskDisplayTempOuter);
    glUniform1f(u.diskSatBoostInner,     snap.diskSatBoostInner);
    glUniform1f(u.diskSatBoostOuter,     snap.diskSatBoostOuter);
    glUniform1i(u.secBHActive,           snap.secBHActive ? 1 : 0);
    glUniform3f(u.secBHDiskNormal,       snap.secBHDiskNormal.x, snap.secBHDiskNormal.y, snap.secBHDiskNormal.z);
    glUniform1f(u.secBHSpin,             snap.secBHSpin);
    glUniform1f(u.secBHDiskInner,        snap.secBHDiskInner);
    glUniform1f(u.secBHDiskOuter,        snap.secBHDiskOuter);
    glUniform1f(u.secBHDiskStrength,     snap.secBHDiskStrength);
    glUniform3f(u.secBHColorInner,       snap.secBHColorInner.x, snap.secBHColorInner.y, snap.secBHColorInner.z);
    glUniform3f(u.secBHColorOuter,       snap.secBHColorOuter.x, snap.secBHColorOuter.y, snap.secBHColorOuter.z);
    glUniform1i(u.secBHShowJets,         snap.secBHShowJets ? 1 : 0);
    glUniform3f(u.secBHJetColor,         snap.secBHJetColor.x, snap.secBHJetColor.y, snap.secBHJetColor.z);
    glUniform1f(u.secBHJetRadius,        snap.secBHJetRadius);
    glUniform1f(u.secBHJetLength,        snap.secBHJetLength);
    glUniform1f(u.spinParameter,         snap.bhSpin);
    glUniform1f(u.jetRadius,             snap.jetRadius);
    glUniform1f(u.jetLength,             snap.jetLength);
    glUniform3f(u.jetColor,              snap.jetColor.x, snap.jetColor.y, snap.jetColor.z);
    glUniform1i(u.showJets,              snap.jetsEnabled    ? 1 : 0);
    glUniform1i(u.showBLR,               snap.blrEnabled     ? 1 : 0);
    glUniform1i(u.showDoppler,           snap.dopplerEnabled  ? 1 : 0);
    glUniform1i(u.showBlueshift,         snap.blueshiftEnabled ? 1 : 0);
    glUniform1f(u.blrInnerRadius,        snap.blrInnerRadius);
    glUniform1f(u.blrOuterRadius,        snap.blrOuterRadius);
    glUniform1f(u.blrThickness,          snap.blrThickness);
    glUniform1f(u.blrStrength,           snap.blrStrength);
    glUniform1f(u.uTime,                 snap.animTime);

    glUniform1i(u.showOrbBody,           snap.orbBodyEnabled ? 1 : 0);
    glUniform1i(u.numOrbBodies,          (GLint)snap.orbBodyPositions.size());
    for (size_t i = 0; i < snap.orbBodyPositions.size() && i < 10; ++i) {
        glUniform3f(u.orbBodyPos[i],    snap.orbBodyPositions[i].x, snap.orbBodyPositions[i].y, snap.orbBodyPositions[i].z);
        glUniform1f(u.orbBodyRadius[i], snap.orbBodyRadii[i]);
        glUniform3f(u.orbBodyColor[i],  snap.orbBodyColors[i].x,    snap.orbBodyColors[i].y,    snap.orbBodyColors[i].z);
        glUniform1i(u.orbBodyType[i],   i < snap.orbBodyTypes.size() ? snap.orbBodyTypes[i] : 0);
    }
    glUniform1i(u.showHostGalaxy,   snap.hostGalaxyEnabled ? 1 : 0);
    glUniform1f(u.hostGalaxyRadius, 10.0f);
    glUniform3f(u.hostGalaxyColor,  0.8f, 0.6f, 0.4f);
    glUniform1i(u.showLAB,          snap.labEnabled ? 1 : 0);
    glUniform3f(u.labPos,           50.0f, 0.0f, 0.0f);
    glUniform1f(u.labRadius,        100.0f);
    glUniform3f(u.labColor,         0.4f, 0.7f, 1.0f);
    glUniform1i(u.showCGM,          snap.cgmEnabled ? 1 : 0);
    glUniform1f(u.cgmRadius,        200.0f);
    glUniform3f(u.cgmColor,         0.3f, 0.5f, 0.8f);
    glUniform1i(u.maxStepsOverride, snap.maxSteps);
}

// ────────────────────────────────────────────────────────────
//  Build the OrbitalBody set for a given profile
// ────────────────────────────────────────────────────────────
inline void rebuildOrbBodies(std::vector<OrbitalBody>& out,
                             const BlackHoleProfile&  prof,
                             const cfg::SimConfig&    cfgRef)
{
    out.clear();
    if (!prof.galaxyBodies.empty()) {
        out.reserve(prof.galaxyBodies.size());
        for (size_t i = 0; i < prof.galaxyBodies.size(); ++i) {
            const auto& gb = prof.galaxyBodies[i];
            cfg::OrbitalConfig oc = profiles::makeOrbitalBody(gb);
            // Preserve per-profile time-scale tuning from legacy orbital config
            oc.timeScale     = cfgRef.orbital.timeScale;
            oc.fastTimeScale = cfgRef.orbital.fastTimeScale;
            out.emplace_back(oc);
            // Stagger initial mean anomaly so high-eccentricity orbits don't
            // all start clustered at periapsis.
            constexpr float goldenAngle = 2.39996323f;
            out.back().setInitialPhase(static_cast<float>(i) * goldenAngle);
        }
    } else {
        out.emplace_back(cfgRef.orbital);
    }
}

// ────────────────────────────────────────────────────────────
//  All per-instance simulation state
// ────────────────────────────────────────────────────────────
struct State {
    // Profiles
    std::array<BlackHoleProfile, profiles::NUM_PROFILES> profilesArr;
    int            profileIdx = 0;
    cfg::SimConfig config;

    // Orbital bodies
    std::vector<OrbitalBody> orbBodies;

    // GL programs and cached uniform locations
    GLProgram     photorealProgram;
    GLProgram     simpleProgram;
    bool          havePhotoreal = false;
    bool          haveSimple    = false;
    SceneUniforms photorealUfs{};
    SceneUniforms simpleUfs{};
    GLProgram*    activeProgram = nullptr;
    SceneUniforms activeUfs{};
    bool          cinematicMode = false;

    // Geometry
    GLVertexArray quad;

    // Textures
    GLTexture2D bgTex;
    GLTexture2D photorealDiskTex;
    GLTexture2D simpleDiskTex;

    // Bloom
    BloomPipeline bloom;
    int           lastW = 0, lastH = 0;

    // HUD
    GLBitmapFont       glFont;
    bool               fontLoaded = false;
    hud::OverlaysState overlays;
    hud::HudFrame      hudFrame;
    hud::PresetMenuState presetMenu;
    // When true, the host (e.g. standalone SFML window with ImGui::SFML
    // initialized) draws the preset menu via ImGui instead of the bespoke
    // GLBitmapFont overlay. The Qt-embedded build leaves this false.
    bool               useImGuiHud = false;

    // Camera/input
    CameraController camera;
    KeyState         keys;

    // Toggles
    bool jetsEnabled       = false;
    bool blrEnabled        = false;
    bool orbBodyEnabled    = false;
    bool dopplerEnabled    = false;
    bool hostGalaxyEnabled = false;
    bool labEnabled        = false;
    bool cgmEnabled        = false;
    bool blueshiftEnabled  = true;
    bool showHUD           = true;
    bool showDebugHUD      = false;

    // Animation
    static constexpr float animSpeedPresets[3] = { 1.0f, 3.0f, 8.0f };
    static constexpr int   animSpeedCount      = 3;
    int   animSpeedIdx = 1;
    float animSpeed    = 3.0f;
    float animTime     = 0.0f;
    float totalTime    = 0.0f;

    // FPS
    sf::Clock frameClock;
    sf::Clock fpsClock;
    int       fpsFrameCount = 0;
    float     fpsDisplay    = 0.0f;

    // Keybinds
    ActionKeybinds actionKeys;

    // Persistent snapshot
    PhysicsSnapshot snap{};

    // RK4 orbit + null-geodesic photon-ray overlay (world-space lines).
    // Cached vertex buffer; rebuilt on profile switch / toggle change.
    PhysicsOverlay physOverlay;

    // Tidal disruption event (3D)
    struct TidalEvent3D {
        bool   active        = false;
        float  flashTimer    = 0.0f;
        glm::vec3 eventPos   = glm::vec3(0.0f);
        std::vector<glm::vec3> debrisPos;
        std::vector<glm::vec3> debrisVel;
        std::vector<float>     debrisLife;
        std::vector<float>     debrisMaxLife;
        std::vector<bool>      isFallback;
        static constexpr float FLASH_DURATION  = 0.45f;
        static constexpr float DEBRIS_LIFETIME = 6.0f;
        static constexpr float STREAM_LIFETIME = 9.0f;
    } tde3D;
    std::vector<bool> orbBodyDisrupted; // sized to match orbBodies

    // ── Visual black hole merger (mirrors the 2D simulator) ──────────────
    // A secondary object spirals into the primary on a tilted orbit, both
    // wobble about the barycenter, then coalesce in a flash that leaves a
    // grown remnant and a gravitational-wave shockwave ring. During the
    // inspiral the secondary also gravitationally perturbs the orbiting bodies
    // and tidally disrupts those that stray too close (see applyMergerInfluence).
    // The coalescence itself is still a *visual* effect (no Peters realism).
    struct MergerState3D {
        bool  active = false;   // inspiral or aftermath running
        bool  merged = false;   // reached coalescence (flash/aftermath phase)
        MergerSecondaryKind3D kind = MergerSecondaryKind3D::BlackHole;

        // Cinematic timing dial, mirrors the 2D MergerState::TimeScale.
        enum class TimeScale { Cinematic = 0, Default = 1, Realtime = 2 } timeScale = TimeScale::Default;
        float timeScaleFactor() const {
            switch (timeScale) {
                case TimeScale::Cinematic: return 0.35f;
                case TimeScale::Realtime:  return 3.0f;
                default:                   return 1.0f;
            }
        }
        const char* timeScaleName() const {
            switch (timeScale) {
                case TimeScale::Cinematic: return "Cinematic";
                case TimeScale::Realtime:  return "Realtime";
                default:                   return "Default";
            }
        }

        // Inspiral geometry (all separations in Rs world units).
        float  r    = 16.0f;   // current separation
        float  phi  = 0.0f;    // orbital angle
        float  incl = 0.35f;   // orbital-plane tilt [rad]
        float  bhR0 = 1.0f;    // primary Schwarzschild radius captured at start
        double massSolar = 10.0; // secondary mass [Msun]
        double q    = 0.5;     // m2 / (m1 + m2)
        float  secRadius = 1.0f; // secondary visual radius [Rs]

        // ── Secondary visual identity ────────────────────────────────────
        // Grabbed when the merger starts so the incoming secondary looks the
        // same as it would as a standalone primary (disk / spin / jets), as long
        // as the preset points at one. The disk/jet sizes are ratios of secRadius
        // applied in the shader, so they scale with the secondary no matter how
        // big it ends up on screen.
        struct SecVisual {
            bool  hasDisk = false;   // draw an accretion disk around the secondary
            bool  hasJets = false;   // draw relativistic jets
            float spin    = 0.0f;    // a* → photon-ring tightness
            float diskInnerRatio = 2.2f;   // × secRadius
            float diskOuterRatio = 12.0f;  // × secRadius
            float diskStrength   = 0.0f;    // overall disk brightness (0 = none)
            glm::vec3 colorInner = glm::vec3(1.0f, 0.85f, 0.55f); // inner disk RGB
            glm::vec3 colorOuter = glm::vec3(0.9f, 0.45f, 0.20f); // outer disk RGB
            glm::vec3 jetColor   = glm::vec3(0.4f, 0.7f, 1.0f);
            float jetRadiusRatio = 0.15f;  // × secRadius
            float jetLengthRatio = 10.0f;  // × secRadius
            glm::vec3 diskNormal = glm::vec3(0.0f, 1.0f, 0.0f);   // disk-plane normal
        } secVis;

        // Coalescence + remnant growth.
        float flashTimer    = 0.0f;
        float growthT       = 0.0f;  // 0..1 growth animation
        float preRadius     = 1.0f,  targetRadius  = 1.0f;
        float preDiskInner  = 2.0f,  targetDiskInner  = 2.0f;
        float preDiskOuter  = 20.0f, targetDiskOuter  = 20.0f;

        // Gravitational-wave shockwave ring.
        bool      ringActive = false;
        float     ringT      = 0.0f;  // 0..1
        glm::vec3 eventPos   = glm::vec3(0.0f);
        float     remnantTimer = 0.0f;

        // Death-spiral trail (secondary world positions, oldest → newest).
        std::deque<glm::vec3> trail;

        static constexpr float START_SEP_RS    = 16.0f;
        static constexpr float MERGE_SEP_RS    = 2.6f;
        static constexpr float INSPIRAL_K      = 1600.0f;
        static constexpr float OMEGA_K         = 14.0f;
        static constexpr float FLASH_DURATION  = 0.6f;
        static constexpr float GROWTH_DURATION = 0.6f;
        static constexpr float RING_DURATION   = 1.6f;
        static constexpr float REMNANT_LABEL   = 3.2f;
        static constexpr size_t MAX_TRAIL      = 160;

        // Unit direction from barycenter toward the secondary in the tilted
        // orbital plane (orbit in XZ, tilted about X by `incl`).
        glm::vec3 unitDir() const {
            float cp = std::cos(phi), sp = std::sin(phi);
            float ci = std::cos(incl), si = std::sin(incl);
            return glm::vec3(cp, sp * si, sp * ci);
        }
        // Secondary sits at +(1−q)·sep, primary recoils to −q·sep.
        glm::vec3 secondaryWorld(const glm::vec3& bary) const {
            return bary + unitDir() * (r * float(1.0 - q));
        }
        glm::vec3 primaryOffset() const {
            return -unitDir() * (r * float(q));
        }
    } merger;

    // Merger selector menu (ImGui), mirrors presetMenu.
    struct MergerMenuState {
        bool  open     = false;
        int   selected = 0;
        bool  scrollToSelected = false; // request keyboard-driven auto-scroll
        float t        = 0.0f;   // fade/slide 0..1
        static constexpr float transitionSec = 0.18f;
    } mergerMenu;

    // Mouse-look state (for standalone window cursor grab; widget tracks its own)
    bool  looking = false;
    float lastMouseX = 0.0f;
    float lastMouseY = 0.0f;

    explicit State(const cfg::CameraConfig& camCfg) : camera(camCfg) {
        snap.orbBodyPositions.reserve(10);
        snap.orbBodyRadii.reserve(10);
        snap.orbBodyColors.reserve(10);
        snap.orbBodyTypes.reserve(10);
        snap.orbBodyLabels.reserve(10);
    }
};

// ────────────────────────────────────────────────────────────
//  Initialisation helpers
// ────────────────────────────────────────────────────────────

// Load profiles, set initial profile, build orbital bodies, set toggle defaults.
inline void initProfiles(State& s, int initialIdx = 0) {
    auto bh = profiles::allProfiles();
    s.profilesArr = std::move(bh);
    s.profileIdx  = std::max(0, std::min<int>(initialIdx, profiles::NUM_PROFILES - 1));
    s.config      = s.profilesArr[s.profileIdx].config;
    rebuildOrbBodies(s.orbBodies, s.profilesArr[s.profileIdx], s.config);
    s.orbBodyDisrupted.assign(s.orbBodies.size(), false);
    s.tde3D = State::TidalEvent3D{};
    const auto& prof = s.profilesArr[s.profileIdx];
    s.jetsEnabled       = prof.defaultJets;
    s.blrEnabled        = prof.defaultBLR;
    s.orbBodyEnabled    = prof.defaultOrbBody;
    s.dopplerEnabled    = prof.defaultDoppler;
    s.hostGalaxyEnabled = prof.defaultHostGalaxy;
    s.labEnabled        = prof.defaultLAB;
    s.cgmEnabled        = prof.defaultCGM;
}

// Try the standard list of monospace font paths; populate glFont.
// Caller is responsible for activating the GL context first.
inline void initFont(State& s, sf::Font& tmpFont) {
    auto tryFont = [&](const std::string& path) -> bool {
        std::ifstream probe(path);
        return probe.good() && tmpFont.openFromFile(path);
    };
    bool loaded =
        tryFont("/usr/share/fonts/dejavu/DejaVuSansMono.ttf")                       ||
        tryFont("/usr/share/fonts/liberation-fonts/LiberationMono-Regular.ttf")     ||
        tryFont("/run/host/fonts/truetype/dejavu/DejaVuSansMono.ttf")               ||
        tryFont("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf")              ||
        tryFont("/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf")  ||
        tryFont("/System/Library/Fonts/Menlo.ttc")                                  ||
        tryFont("/System/Library/Fonts/SFNSMono.ttf")                               ||
        tryFont("C:/Windows/Fonts/consola.ttf")                                     ||
        tryFont("C:/Windows/Fonts/cour.ttf")                                        ||
        tryFont("C:/Windows/Fonts/lucon.ttf");
    if (loaded) s.fontLoaded = s.glFont.initFromFont(tmpFont, 14);
}

// Compile both the photoreal and the simple shader from disk.
// Returns false if neither compiled (caller should bail out).
inline bool initShaders(State& s, ResourceManager& res) {
    GLShader sharedVS(GL_VERTEX_SHADER, vertexShaderSrc);
    if (!sharedVS) { std::cerr << "[bh3d] vertex shader failed\n"; return false; }
    {
        std::string fragPath = res.find("BlackHole3D_PhotorealDisk.frag");
        if (!fragPath.empty()) {
            std::string src = readFile(fragPath.c_str());
            if (!src.empty()) {
                GLShader fs(GL_FRAGMENT_SHADER, src.c_str());
                if (fs && s.photorealProgram.link(sharedVS.id(), fs.id())) {
                    s.havePhotoreal = true;
                    std::cerr << "[bh3d] photoreal shader compiled: " << fragPath << "\n";
                }
            }
        }
    }
    {
        std::string fragPath = res.find("BlackHole3D.frag");
        if (!fragPath.empty()) {
            std::string src = readFile(fragPath.c_str());
            if (!src.empty()) {
                GLShader fs(GL_FRAGMENT_SHADER, src.c_str());
                if (fs && s.simpleProgram.link(sharedVS.id(), fs.id())) {
                    s.haveSimple = true;
                    std::cerr << "[bh3d] simple shader compiled: " << fragPath << "\n";
                }
            }
        }
    }
    if (!s.havePhotoreal && !s.haveSimple) return false;
    if (s.havePhotoreal) s.photorealUfs = SceneUniforms::lookup(s.photorealProgram);
    if (s.haveSimple)    s.simpleUfs    = SceneUniforms::lookup(s.simpleProgram);
    s.cinematicMode  = s.havePhotoreal;
    s.activeProgram  = s.cinematicMode ? &s.photorealProgram : &s.simpleProgram;
    s.activeUfs      = s.cinematicMode ? s.photorealUfs      : s.simpleUfs;
    return true;
}

inline void initTextures(State& s, ResourceManager& res) {
    {
        std::string bgPath = res.find("background.png");
        if (!bgPath.empty()) s.bgTex = loadTexture2D(bgPath.c_str());
        if (!s.bgTex) s.bgTex = createFallbackWhiteTexture();
    }
    s.photorealDiskTex = createFallbackWhiteTexture();
    {
        std::string diskPath = res.find("disk_texture.png");
        if (!diskPath.empty()) s.simpleDiskTex = loadTexture2D(diskPath.c_str());
        if (!s.simpleDiskTex) s.simpleDiskTex = createDiskTextureProcedural(512);
    }
}

// ────────────────────────────────────────────────────────────
//  Visual black hole merger, Basically ported from 2D but with the 3d tilt and other parameters
// ────────────────────────────────────────────────────────────

// Map a secondary kind to a shader body archetype. Black holes use the
// dedicated BODY_BLACKHOLE (=7); the rest reuse existing GalaxyBody3DType
// ints so the shader renders them through the normal orbiting-body path.
inline int mergerBodyShaderType(MergerSecondaryKind3D k) {
    switch (k) {
        case MergerSecondaryKind3D::BlackHole:   return 7; // BODY_BLACKHOLE
        case MergerSecondaryKind3D::NeutronStar: return static_cast<int>(GalaxyBody3DType::NeutronStar);
        case MergerSecondaryKind3D::Pulsar:      return static_cast<int>(GalaxyBody3DType::NeutronStar);
        case MergerSecondaryKind3D::Star:        return static_cast<int>(GalaxyBody3DType::Star);
        case MergerSecondaryKind3D::WhiteDwarf:  return static_cast<int>(GalaxyBody3DType::WhiteDwarf);
    }
    return 7;
}

inline const char* mergerKindLabel(MergerSecondaryKind3D k) {
    switch (k) {
        case MergerSecondaryKind3D::BlackHole:   return "Black Hole";
        case MergerSecondaryKind3D::NeutronStar: return "Neutron Star";
        case MergerSecondaryKind3D::Pulsar:      return "Pulsar";
        case MergerSecondaryKind3D::Star:        return "Star";
        case MergerSecondaryKind3D::WhiteDwarf:  return "White Dwarf";
    }
    return "Secondary";
}

inline glm::vec3 mergerKindColor(MergerSecondaryKind3D k) {
    switch (k) {
        case MergerSecondaryKind3D::BlackHole:   return glm::vec3(0.02f, 0.02f, 0.03f);
        case MergerSecondaryKind3D::NeutronStar: return glm::vec3(1.00f, 0.35f, 0.90f);
        case MergerSecondaryKind3D::Pulsar:      return glm::vec3(0.70f, 0.85f, 1.00f);
        case MergerSecondaryKind3D::Star:        return glm::vec3(1.00f, 0.93f, 0.78f);
        case MergerSecondaryKind3D::WhiteDwarf:  return glm::vec3(1.00f, 1.00f, 0.95f);
    }
    return glm::vec3(1.0f);
}

// Approximate blackbody colour (Tanner Helland algorithm), mirroring the
// shader's kelvinToRGB so a merger secondary's disk colours can be precomputed
// on the CPU and shipped as uniforms (the fast shader has no kelvinToRGB).
inline glm::vec3 bh3dKelvinToRGB(float kelvin) {
    float t = std::clamp(kelvin, 1000.0f, 40000.0f) / 100.0f;
    float r, g, b;
    if (t <= 66.0f) {
        r = 255.0f;
        g = std::clamp(99.4708025861f * std::log(std::max(t, 1.0f)) - 161.1195681661f, 0.0f, 255.0f);
    } else {
        r = std::clamp(329.698727446f * std::pow(t - 60.0f, -0.1332047592f), 0.0f, 255.0f);
        g = std::clamp(288.1221695283f * std::pow(t - 60.0f, -0.0755148492f), 0.0f, 255.0f);
    }
    if (t >= 66.0f)      b = 255.0f;
    else if (t <= 19.0f) b = 0.0f;
    else                 b = std::clamp(138.5177312231f * std::log(t - 10.0f) - 305.0447927307f, 0.0f, 255.0f);
    return glm::vec3(r, g, b) / 255.0f;
}

// Begin a merger: place the secondary at the outer separation and seed the
// inspiral. The cinematic time-scale choice persists across re-triggers. Do note that the secondary mass is NOT clamped, 
//meaning that something like a 100 Msun secondary can be merged into a 10 Msun primary, which is not physically realistic but is visually interesting.
inline void startMerger3D(State& s, MergerSecondaryKind3D kind, double massSolar,
                          const char* secProfileName = nullptr) { 
    const auto prevScale = s.merger.timeScale;
    // Undo any remnant growth left by a previous, already-coalesced merger so
    // repeated mergers grow from the profile baseline instead of compounding
    // (1.8× per merger). A profile switch already resets merger + config, so
    // this only matters for back-to-back mergers on the same profile.
    if (s.merger.merged && s.merger.growthT > 0.0f) {
        s.config.blackHole.radius = s.merger.preRadius;
        s.config.disk.innerRadius = s.merger.preDiskInner;
        s.config.disk.outerRadius = s.merger.preDiskOuter;
    }
    s.merger = State::MergerState3D{};
    s.merger.timeScale = prevScale;

    auto& m = s.merger;
    m.active    = true;
    m.merged    = false;
    m.kind      = kind;
    m.massSolar = massSolar;
    m.bhR0      = std::max(0.05f, s.config.blackHole.radius);
    m.r         = State::MergerState3D::START_SEP_RS * m.bhR0;
    m.phi       = 0.0f;
    m.incl      = 0.32f;   // gentle tilt so the inspiral reads as 3D

    const double m1 = std::max(1e-6, s.profilesArr[s.profileIdx].massSolar);
    const double m2 = std::max(1e-9, massSolar);
    m.q = m2 / (m1 + m2);

    switch (kind) {
        case MergerSecondaryKind3D::BlackHole: {
            // Visual radius grows with mass but is compressed so an SMBH
            // secondary stays on-screen rather than swallowing the frame.
            double rel = std::cbrt(m2 / m1);            // equal-mass → 1
            m.secRadius = (float)std::clamp(rel, 0.25, 2.6) * m.bhR0;
            break;
        }
        case MergerSecondaryKind3D::NeutronStar:
        case MergerSecondaryKind3D::Pulsar:     m.secRadius = 0.42f * m.bhR0; break;
        case MergerSecondaryKind3D::WhiteDwarf: m.secRadius = 0.52f * m.bhR0; break;
        case MergerSecondaryKind3D::Star:       m.secRadius = 0.85f * m.bhR0; break;
    }

    // ── Secondary visual identity ────────────────────────────────────────
    // Hand a black-hole secondary the same disk / spin / jet look it would have
    // as a standalone primary (when the preset points at one). The other kinds
    // (NS/pulsar/star/WD) just keep their own surface shading and don't get a
    // disk.
    m.secVis = State::MergerState3D::SecVisual{};
    if (kind == MergerSecondaryKind3D::BlackHole) {
        auto& v = m.secVis;
        // Tilt the secondary's disk so it reads as a 3D object rather than a
        // face-on ring: normal derived from the inspiral inclination.
        {
            float ci = std::cos(m.incl), si = std::sin(m.incl);
            v.diskNormal = glm::normalize(glm::vec3(0.18f, ci, si * 0.85f));
        }
        const ::BlackHoleProfile* prof =
            secProfileName ? profiles::findProfileByName(secProfileName) : nullptr;
        if (prof) {
            const auto& c = prof->config;
            float pR0 = std::max(0.05f, c.blackHole.radius);
            v.hasDisk        = true;
            v.spin           = c.blackHole.spinParameter;
            v.diskInnerRatio = std::clamp(c.disk.innerRadius / pR0, 1.4f, 4.0f);
            v.diskOuterRatio = std::clamp(c.disk.outerRadius / pR0, 6.0f, 26.0f);
            v.colorInner     = bh3dKelvinToRGB(c.disk.displayTempInner);
            v.colorOuter     = bh3dKelvinToRGB(c.disk.displayTempOuter);
            // Saturation boost baked into the precomputed colours so the shader
            // stays simple.
            auto boost = [](glm::vec3 col, float sat) {
                float lum = glm::dot(col, glm::vec3(0.2126f, 0.7152f, 0.0722f));
                return glm::max(glm::mix(glm::vec3(lum), col, sat), glm::vec3(0.0f));
            };
            v.colorInner = boost(v.colorInner, c.disk.saturationBoostInner);
            v.colorOuter = boost(v.colorOuter, c.disk.saturationBoostOuter);
            v.diskStrength   = 1.0f;
            v.hasJets        = prof->defaultJets;
            v.jetColor       = c.jet.color;
            v.jetRadiusRatio = std::clamp(c.jet.radius / pR0, 0.05f, 0.6f);
            v.jetLengthRatio = std::clamp(c.jet.length / pR0, 4.0f, 22.0f);
        } else {
            // Generic archetype with no dedicated profile: a tasteful warm
            // disk derived from mass (heavier → cooler/redder, quasar-like
            // secondaries stay hot). No jets.
            float logM = (float)std::log10(std::max(1e-6, massSolar));
            float tInner = (float)std::clamp(9000.0 - logM * 600.0, 3200.0, 9000.0);
            float tOuter = tInner * 0.42f;
            v.hasDisk        = true;
            v.spin           = 0.5f;
            v.diskInnerRatio = 2.2f;
            v.diskOuterRatio = 12.0f;
            v.colorInner     = bh3dKelvinToRGB(tInner);
            v.colorOuter     = bh3dKelvinToRGB(tOuter);
            v.diskStrength   = 0.85f;
            v.hasJets        = false;
        }
    }

    s.mergerMenu.open = false;
    // Reset any prior orbit perturbation / disruption so each merger starts clean.
    for (auto& b : s.orbBodies) b.resetPerturbation();
    s.orbBodyDisrupted.assign(s.orbBodies.size(), false);
    std::cerr << "[bh3d] merger start: " << mergerKindLabel(kind)
              << " (" << massSolar << " Msun, q=" << m.q << ")\n";
}

// On coalescence, spawn matter ejecta for non-BH secondaries into the shared
// tidal-event particle system (already advanced + drawn elsewhere). Black-hole
// mergers stay clean: only the gravitational-wave ring, no matter debris.
inline void spawnMergerAftermath(State& s) {
    auto& m = s.merger;
    if (m.kind == MergerSecondaryKind3D::BlackHole) return;

    s.tde3D = State::TidalEvent3D{};
    s.tde3D.active     = true;
    s.tde3D.eventPos   = m.eventPos;
    s.tde3D.flashTimer = State::TidalEvent3D::FLASH_DURATION;

    auto rnd = [](float a, float b) {
        return a + (b - a) * (float)std::rand() / (float)RAND_MAX;
    };

    const bool  compact = (m.kind == MergerSecondaryKind3D::NeutronStar
                        || m.kind == MergerSecondaryKind3D::Pulsar);
    const int   nDebris = compact ? 40 : 64;
    const float speed   = compact ? 6.0f : 3.5f;   // NS ejecta is faster/hotter
    const float life    = compact ? State::TidalEvent3D::DEBRIS_LIFETIME * 0.7f
                                  : State::TidalEvent3D::STREAM_LIFETIME;

    for (int i = 0; i < nDebris; ++i) {
        float a  = rnd(0.0f, 6.2831853f);
        float el = rnd(-0.5f, 0.5f);
        glm::vec3 dir = glm::vec3(std::cos(a) * std::cos(el),
                                  std::sin(el),
                                  std::sin(a) * std::cos(el));
        float v = speed * rnd(0.5f, 1.0f);
        s.tde3D.debrisPos.push_back(m.eventPos + dir * (m.bhR0 * rnd(1.0f, 2.5f)));
        s.tde3D.debrisVel.push_back(dir * v);
        float ml = life * rnd(0.6f, 1.0f);
        s.tde3D.debrisLife.push_back(ml);
        s.tde3D.debrisMaxLife.push_back(ml);
        // NS/pulsar → hot blue burst (all fallback); star/WD → mixed stream.
        s.tde3D.isFallback.push_back(compact ? true : (i % 3 == 0));
    }
}

// Spawn a detailed tidal-disruption stream when the in-spiralling merger
// secondary shreds an orbiting body. Unlike spawnMergerAftermath this APPENDS
// to the shared tidal particle system (several bodies can be torn apart during
// one merger) and lays the debris out as a stretched, two-tailed stream along
// the body's orbit with an infalling component toward the secondary, the
// classic spaghettification look, rendered in detail by the HUD overlay.
inline void spawnBodyDisruption(State& s, const glm::vec3& pos,
                                const glm::vec3& tangent,
                                const glm::vec3& towardSec,
                                int count, bool hot) {
    s.tde3D.active     = true;
    s.tde3D.eventPos   = pos;
    s.tde3D.flashTimer = std::max(s.tde3D.flashTimer,
                                  State::TidalEvent3D::FLASH_DURATION);

    auto rnd = [](float a, float b) {
        return a + (b - a) * (float)std::rand() / (float)RAND_MAX;
    };
    glm::vec3 t  = (glm::length(tangent)   > 1e-4f) ? glm::normalize(tangent)   : glm::vec3(1, 0, 0);
    glm::vec3 ts = (glm::length(towardSec) > 1e-4f) ? glm::normalize(towardSec) : glm::vec3(0, 1, 0);
    glm::vec3 up = glm::cross(t, ts);
    up = (glm::length(up) > 1e-4f) ? glm::normalize(up) : glm::vec3(0, 1, 0);

    for (int i = 0; i < count; ++i) {
        // Two-tailed: roughly half the matter leads, half trails along the orbit.
        float along  = rnd(-1.0f, 1.0f);
        float infall = rnd(0.0f, 1.0f);
        glm::vec3 jitter = up * rnd(-0.35f, 0.35f)
                         + ts * rnd(-0.15f, 0.15f)
                         + t  * rnd(-0.15f, 0.15f);
        // Position: elongated along the orbit, drawn slightly toward the secondary.
        glm::vec3 p = pos + t * (along * rnd(0.6f, 2.2f)) + ts * (infall * 0.6f) + jitter * 0.6f;
        // Velocity: streams away along ±tangent, with an infall pull and spread.
        glm::vec3 v = t * (along * rnd(0.6f, 1.6f))
                    + ts * (infall * rnd(0.4f, 1.2f))
                    + jitter * rnd(0.3f, 0.9f);
        s.tde3D.debrisPos.push_back(p);
        s.tde3D.debrisVel.push_back(v);
        float ml = State::TidalEvent3D::STREAM_LIFETIME * rnd(0.55f, 1.0f);
        s.tde3D.debrisLife.push_back(ml);
        s.tde3D.debrisMaxLife.push_back(ml);
        // Hot (NS/WD) → mostly blue-white fallback; star/gas → mixed warm stream.
        s.tde3D.isFallback.push_back(hot ? (i % 4 != 0) : (i % 3 == 0));
    }
}

// Advance the merger one frame. Drives the inspiral, triggers coalescence,
// then animates the flash, remnant growth, and gravitational-wave ring.
inline void updateMerger3D(State& s, float dt) {
    auto& m = s.merger;
    if (!m.active) return;

    // Merger runs on its own cinematic clock, independent of animSpeed, so the
    // choreography is actually predictable regardless of the global time dial. 
    // Clamp the frame delta so a stutter can't make the explicit integrator
    // overshoot the merge, and so the choreography stays bounded on slow frames.
    const float mdt = std::min(dt * m.timeScaleFactor(), 0.10f);
    const glm::vec3 bary = s.config.blackHole.position;

    if (!m.merged) {
        // ── Inspiral ────────────────────────────────────────────────
        // Coalesce when the horizons meet, not at a fixed separation, so a
        // large / heavy secondary doesn't visibly interpenetrate the primary.
        const float mergeSep = std::max(State::MergerState3D::MERGE_SEP_RS * m.bhR0,
                                        (m.bhR0 + m.secRadius) * 1.05f);
        // Symmetric-mass-ratio factor (η/0.25), floored so even extreme
        // ratios still coalesce in a watchable time (visual, not EMRI).
        const float eta4     = 4.0f * (float)(m.q * (1.0 - m.q));
        const float etaScale = std::clamp(eta4, 0.2f, 1.0f);

        // Integrate in units of the primary radius (ρ = r / bhR0) so the inspiral
        // rate is independent of the primary's size, and advance in bounded
        // sub-steps so the trajectory is frame-rate independent.
        const float invR0    = 1.0f / std::max(m.bhR0, 1e-3f);
        const float rhoMerge = mergeSep * invR0;
        float       rho      = m.r * invR0;
        const int   nSub     = std::clamp((int)std::ceil(mdt / 0.01f), 1, 32);
        const float sub      = mdt / (float)nSub;
        bool        reached  = false;
        for (int k = 0; k < nSub; ++k) {
            const float r3 = rho * rho * rho;
            rho   -= (State::MergerState3D::INSPIRAL_K * etaScale * sub) / std::max(r3, 1e-3f);
            m.phi += (State::MergerState3D::OMEGA_K * sub) / std::pow(std::max(rho, 0.2f), 1.5f);
            if (rho <= rhoMerge) { rho = rhoMerge; reached = true; break; }
        }
        m.r = rho * m.bhR0;

        m.trail.push_back(m.secondaryWorld(bary));
        while (m.trail.size() > State::MergerState3D::MAX_TRAIL) m.trail.pop_front();

        if (reached) {
            // ── Coalescence trigger ─────────────────────────────────
            m.merged       = true;
            m.flashTimer   = State::MergerState3D::FLASH_DURATION;
            m.ringActive   = true;
            m.ringT        = 0.0f;
            m.eventPos     = bary;
            m.remnantTimer = State::MergerState3D::REMNANT_LABEL;

            // Bounded remnant growth so the framing doesn't blow up.
            float gFactor      = std::clamp(1.0f + (float)m.q * 1.6f, 1.0f, 1.8f);
            m.preRadius        = s.config.blackHole.radius;
            m.targetRadius     = m.preRadius * gFactor;
            m.preDiskInner     = s.config.disk.innerRadius;
            m.targetDiskInner  = m.preDiskInner * gFactor;
            m.preDiskOuter     = s.config.disk.outerRadius;
            m.targetDiskOuter  = m.preDiskOuter * gFactor;
            m.growthT          = 0.0f;

            spawnMergerAftermath(s);
        }
    } else {
        // ── Aftermath: flash decay, remnant growth, ring expansion ───
        if (m.flashTimer > 0.0f) m.flashTimer -= mdt;

        if (m.growthT < 1.0f) {
            m.growthT = std::min(1.0f, m.growthT + mdt / State::MergerState3D::GROWTH_DURATION);
            float e = m.growthT * m.growthT * (3.0f - 2.0f * m.growthT); // smoothstep
            s.config.blackHole.radius = m.preRadius    + (m.targetRadius    - m.preRadius)    * e;
            s.config.disk.innerRadius = m.preDiskInner + (m.targetDiskInner - m.preDiskInner) * e;
            s.config.disk.outerRadius = m.preDiskOuter + (m.targetDiskOuter - m.preDiskOuter) * e;
        }

        if (m.ringActive) {
            m.ringT += mdt / State::MergerState3D::RING_DURATION;
            if (m.ringT >= 1.0f) { m.ringT = 1.0f; m.ringActive = false; }
        }
        if (m.remnantTimer > 0.0f) m.remnantTimer -= mdt;

        // Dissolve the death-spiral trail after coalescence, dt-scaled so the
        // fade rate is independent of frame rate and pacing.
        if (!m.trail.empty()) {
            const float frac  = std::clamp(mdt / 0.5f, 0.0f, 1.0f);
            const size_t drop = (size_t)std::ceil((float)m.trail.size() * frac);
            for (size_t k = 0; k < drop && !m.trail.empty(); ++k) m.trail.pop_front();
        }

        if (m.flashTimer <= 0.0f && !m.ringActive && m.remnantTimer <= 0.0f) {
            m.active = false;
        }
    }
}

// While the secondary is spiralling in, it gravitationally tugs every orbiting
// body, bending their Keplerian paths, and tidally shreds (or swallows whole)
// those that stray within reach. A more massive secondary (e.g. a TON 618
// SMBH) perturbs and disrupts bodies from much farther out. Disruption reuses
// the shared tidal particle system for the debris streams. Inspiral phase only.
inline void applyMergerInfluence(State& s, float dt) {
    auto& m = s.merger;
    if (!m.active || m.merged) return;
    if (s.orbBodies.empty()) return;

    if (s.orbBodyDisrupted.size() != s.orbBodies.size())
        s.orbBodyDisrupted.assign(s.orbBodies.size(), false);

    const glm::vec3 bary   = s.config.blackHole.position;
    const glm::vec3 secPos = m.secondaryWorld(bary);
    const float q    = (float)std::clamp(m.q, 0.0, 0.999);
    const float mu2  = std::clamp(q / std::max(1e-3f, 1.0f - q), 0.02f, 8.0f);
    const float capR = std::max(m.secRadius * 0.9f, m.bhR0 * 0.4f);
    const float maxOffset = 8.0f * m.bhR0;
    // Perturbation shares the merger's cinematic clock (not the global animation
    // speed) so a body's total orbit-bending is the same regardless of pacing,
    // and can't blow up quadratically when animSpeed is high.
    constexpr float PERTURB_K = 4.8f;
    const float infl = dt * std::max(0.05f, m.timeScaleFactor());
    // A heavier secondary tidally reaches farther.
    const float massReach = std::clamp(0.7f + 0.5f * mu2, 0.7f, 3.0f);

    for (size_t i = 0; i < s.orbBodies.size(); ++i) {
        auto& body = s.orbBodies[i];
        if (body.disrupted()) continue;

        const glm::vec3 p = body.position();
        glm::vec3 d = secPos - p;
        float dist = glm::length(d);
        if (dist < 1e-4f) continue;
        const glm::vec3 dir = d / dist;
        const int type = body.bodyType();

        // Swallowed whole: the body falls inside the secondary's capture radius.
        if (dist < capR) {
            body.setDisrupted();
            s.orbBodyDisrupted[i] = true;
            glm::vec3 tangent = glm::cross(glm::vec3(0, 1, 0), p - bary);
            spawnBodyDisruption(s, p, tangent, dir, 28, true);
            continue;
        }

        // Tidal reach depends on body type (compact remnants hold together
        // longer) and the secondary's mass + size. Extended systems
        // (clusters, dwarf galaxies) are only perturbed, never shredded here.
        const bool shreddable =
            type == static_cast<int>(GalaxyBody3DType::Star)          ||
            type == static_cast<int>(GalaxyBody3DType::CompanionStar) ||
            type == static_cast<int>(GalaxyBody3DType::GasCloud)      ||
            type == static_cast<int>(GalaxyBody3DType::NeutronStar)   ||
            type == static_cast<int>(GalaxyBody3DType::WhiteDwarf);
        if (shreddable) {
            float tFactor;
            if      (type == static_cast<int>(GalaxyBody3DType::NeutronStar)) tFactor = 1.8f;
            else if (type == static_cast<int>(GalaxyBody3DType::WhiteDwarf))  tFactor = 2.0f;
            else if (type == static_cast<int>(GalaxyBody3DType::GasCloud))    tFactor = 3.4f;
            else                                                             tFactor = 3.0f; // Star/Companion
            float tidalR = std::max(m.secRadius * tFactor * massReach, m.bhR0 * 1.3f);
            if (dist < tidalR) {
                body.setDisrupted();
                s.orbBodyDisrupted[i] = true;
                glm::vec3 tangent = glm::cross(glm::vec3(0, 1, 0), p - bary);
                const bool hot = (type == static_cast<int>(GalaxyBody3DType::NeutronStar)
                               || type == static_cast<int>(GalaxyBody3DType::WhiteDwarf));
                const int count = (type == static_cast<int>(GalaxyBody3DType::GasCloud)) ? 72 : 56;
                spawnBodyDisruption(s, p, tangent, dir, count, hot);
                continue;
            }
        }

        // Otherwise: a gravitational tug that visibly bends the orbit inward.
        const float invd2 = 1.0f / (dist * dist);
        const glm::vec3 accel = dir * (mu2 * PERTURB_K * invd2);
        body.applyExternalAccel(accel, infl, maxOffset);
    }
}

// ────────────────────────────────────────────────────────────
//  Per-frame physics AND snapshot construction
// ────────────────────────────────────────────────────────────
inline void tickPhysics(State& s, float dt) {
    s.totalTime += dt;
    s.animTime  += dt * s.animSpeed;
    s.camera.update(dt, s.keys);
    // Orbiting bodies share the global animSpeed multiplier so +/- (and Y)
    // also speed/slow Keplerian orbits, matching the 2D simulator's behavior.
    const float orbDt = dt * s.animSpeed;
    for (auto& body : s.orbBodies) body.update(orbDt);

    // Visual black hole merger (independent cinematic clock).
    updateMerger3D(s, dt);
    // Gravitational perturbation + tidal disruption of orbiting bodies by the
    // in-spiralling secondary (orbits bend, nearby stars/gas get shredded).
    applyMergerInfluence(s, dt);

    // ---- Tidal disruption detection ----
    // NOTE: All 3D orbBodies are on permanent Keplerian orbits (no energy loss /
    // radiation-reaction inspiral). Triggering a TDE when a high-eccentricity body
    // passes periapsis would fire every single orbit, physically nonsensical for
    // stable S-stars, pulsars, etc. TDE in 3D is reserved for bodies that are
    // genuinely accreting (future feature). Detection loop intentionally disabled.
    {

        // Advance existing TDE particles. While a merger is running they share
        // its cinematic clock so ejecta stay in sync with the flash / ring.
        if (s.tde3D.active) {
            const float tdt = s.merger.active ? dt * s.merger.timeScaleFactor() : dt;
            if (s.tde3D.flashTimer > 0.0f)
                s.tde3D.flashTimer -= tdt;

            for (size_t i = 0; i < s.tde3D.debrisPos.size(); ) {
                s.tde3D.debrisLife[i] -= tdt;
                if (s.tde3D.debrisLife[i] <= 0.0f) {
                    // Swap-erase
                    size_t last = s.tde3D.debrisPos.size() - 1;
                    s.tde3D.debrisPos[i]     = s.tde3D.debrisPos[last];
                    s.tde3D.debrisVel[i]     = s.tde3D.debrisVel[last];
                    s.tde3D.debrisLife[i]    = s.tde3D.debrisLife[last];
                    s.tde3D.debrisMaxLife[i] = s.tde3D.debrisMaxLife[last];
                    s.tde3D.isFallback[i]    = s.tde3D.isFallback[last];
                    s.tde3D.debrisPos.pop_back();
                    s.tde3D.debrisVel.pop_back();
                    s.tde3D.debrisLife.pop_back();
                    s.tde3D.debrisMaxLife.pop_back();
                    s.tde3D.isFallback.pop_back();
                } else {
                    s.tde3D.debrisPos[i] += s.tde3D.debrisVel[i] * tdt;
                    ++i;
                }
            }

            if (s.tde3D.flashTimer <= 0.0f && s.tde3D.debrisPos.empty())
                s.tde3D.active = false;
        }
    }

    // Preset-menu transition (open=fade in, closed=fade out).
    const float speed = 1.0f / hud::PresetMenuState::transitionSec;
    if (s.presetMenu.open) s.presetMenu.t = std::min(1.0f, s.presetMenu.t + dt * speed);
    else                   s.presetMenu.t = std::max(0.0f, s.presetMenu.t - dt * speed);

    // Merger-menu transition (mirrors the preset menu).
    const float mSpeed = 1.0f / State::MergerMenuState::transitionSec;
    if (s.mergerMenu.open) s.mergerMenu.t = std::min(1.0f, s.mergerMenu.t + dt * mSpeed);
    else                   s.mergerMenu.t = std::max(0.0f, s.mergerMenu.t - dt * mSpeed);
}

inline void tickFPS(State& s) { // exponential moving average so the readout doesn't jump around every frame
    float dt = s.fpsClock.restart().asSeconds();
    if (dt == 0.f) return; // Prevent div by zero

    float instantFPS = 1.0f / dt;
    float smoothing = 0.05f; // Lower = smoother/slower, Higher = more responsive

    s.fpsDisplay = (instantFPS * smoothing) + (s.fpsDisplay * (1.0f - smoothing));
}

inline void buildSnapshot(State& s, int w, int h, float dt) {
    auto& snap = s.snap;
    snap.cameraPos    = s.camera.position();
    snap.cameraDir    = s.camera.direction();
    snap.cameraUp     = s.camera.up();
    snap.fov          = s.camera.fov();
    snap.roll         = s.camera.roll();
    snap.freelook     = s.camera.isFreelook();
    snap.totalTime    = s.totalTime;
    snap.animTime     = s.animTime;
    snap.animSpeed    = s.animSpeed;
    snap.dt           = dt;
    snap.jetsEnabled       = s.jetsEnabled;
    snap.blrEnabled        = s.blrEnabled;
    snap.dopplerEnabled    = s.dopplerEnabled;
    snap.blueshiftEnabled  = s.blueshiftEnabled;
    snap.cinematicMode     = s.cinematicMode;
    snap.bhRadius          = s.config.blackHole.radius;
    snap.bhPosition        = s.config.blackHole.position;
    snap.bhSpin            = s.config.blackHole.spinParameter;
    snap.secBHActive       = false;  // set true only when a BH secondary inspirals
    snap.diskInnerRadius   = s.config.disk.innerRadius;
    snap.diskOuterRadius   = s.config.disk.outerRadius;
    snap.diskHalfThickness = s.config.disk.halfThickness;
    snap.diskPeakTemp          = s.config.disk.peakTemp;
    snap.diskDisplayTempInner  = s.config.disk.displayTempInner;
    snap.diskDisplayTempOuter  = s.config.disk.displayTempOuter;
    snap.diskSatBoostInner     = s.config.disk.saturationBoostInner;
    snap.diskSatBoostOuter     = s.config.disk.saturationBoostOuter;
    snap.jetRadius         = s.config.jet.radius;
    snap.jetLength         = s.config.jet.length;
    snap.jetColor          = s.config.jet.color;
    snap.blrInnerRadius    = s.config.blr.innerRadius;
    snap.blrOuterRadius    = s.config.blr.outerRadius;
    snap.blrThickness      = s.config.blr.thickness;
    snap.blrStrength       = s.config.blr.strength;

    snap.orbBodyPositions.clear();
    snap.orbBodyRadii.clear();
    snap.orbBodyColors.clear();
    snap.orbBodyTypes.clear();
    snap.orbBodyLabels.clear();
    {
        const auto& curProf = s.profilesArr[s.profileIdx];

        // ── Barycentric binary mode (Gaia BH1/2/3) ──────────────────────────
        // For these systems the companion mass is large enough that the BH
        // noticeably orbits the centre of mass, exactly the astrometric signal
        // Gaia measured.  We keep the Keplerian integrator running in
        // BH-centred coords but, at snapshot time, shift everything so the
        // *barycenter* sits at world origin:
        //   r_BH_bary      = -r_comp × q           (q = M_comp / (M_BH + M_comp))
        //   r_comp_bary    =  r_comp × (1 − q)     (= r_comp × M_BH/(M_BH+M_comp))
        snap.barycentricMode = curProf.isBinaryWithBarycenter
                               && !s.orbBodies.empty()
                               && s.orbBodyEnabled;
        snap.barycenterPos   = glm::vec3(0.0f);  // always world origin

        glm::vec3 bhOffset(0.0f); // FIXME: if we ever add non-binary profiles with barycenterMode=true, we'll need to compute an actual barycenter offset for those too
        float     baryScale = 1.0f;
        if (snap.barycentricMode) {
            double M_bh   = curProf.massSolar;
            double M_comp = curProf.companionMassSolar;
            double q      = M_comp / (M_bh + M_comp);
            baryScale     = (float)(1.0 - q);   // companion scale

            // Companion world position from the integrator (BH-centred)
            glm::vec3 compPos = s.orbBodies[0].position();
            // BH offset: moves *opposite* the companion by the mass fraction
            bhOffset = -compPos * (float)q;
        }

        // Apply BH offset
        snap.bhPosition = s.config.blackHole.position + bhOffset;

        for (size_t i = 0; i < s.orbBodies.size(); ++i) {
            const auto& body = s.orbBodies[i];
            // Bodies tidally disrupted / swallowed by a merger secondary have
            // become a debris stream, don't draw them as intact bodies.
            if (body.disrupted()) continue;
            // Scale companion position to barycenter-centred world coords
            float scale = (snap.barycentricMode && i == 0) ? baryScale : 1.0f;
            snap.orbBodyPositions.push_back(body.position() * scale);
            snap.orbBodyRadii.push_back(body.bodyRadius());
            snap.orbBodyColors.push_back(body.bodyColor());
            snap.orbBodyTypes.push_back(body.bodyType());
            std::string lbl;
            if (i < curProf.galaxyBodies.size() && !curProf.galaxyBodies[i].label.empty()) {
                lbl = curProf.galaxyBodies[i].label;
            } else {
                lbl = hud::defaultLabelForType(body.bodyType());
            }
            snap.orbBodyLabels.push_back(std::move(lbl));
        }

        // ── Visual merger: inject the inspiralling secondary ────────────────
        // During the inspiral the secondary rides into the primary as an extra
        // orbiting body, and the primary recoils opposite it about the
        // barycenter. This reuses the existing orbBody shader path entirely.
        if (s.merger.active && !s.merger.merged) {
            const glm::vec3 bary   = s.config.blackHole.position;
            const glm::vec3 secPos = s.merger.secondaryWorld(bary);
            // If the profile has no orbiting bodies, show ONLY the secondary so
            // an incidental legacy body doesn't appear mid-merger.
            if (!s.orbBodyEnabled) {
                snap.orbBodyPositions.clear();
                snap.orbBodyRadii.clear();
                snap.orbBodyColors.clear();
                snap.orbBodyTypes.clear();
                snap.orbBodyLabels.clear();
                snap.barycentricMode = false;
            }
            // Always show the inspiralling secondary: if the profile already
            // fills every shader body slot, evict the last one so the merger's
            // secondary is never silently dropped.
            if (snap.orbBodyPositions.size() >= 10) {
                snap.orbBodyPositions.pop_back();
                snap.orbBodyRadii.pop_back();
                snap.orbBodyColors.pop_back();
                snap.orbBodyTypes.pop_back();
                snap.orbBodyLabels.pop_back();
            }
            snap.orbBodyPositions.push_back(secPos);
            snap.orbBodyRadii.push_back(s.merger.secRadius);
            snap.orbBodyColors.push_back(mergerKindColor(s.merger.kind));
            snap.orbBodyTypes.push_back(mergerBodyShaderType(s.merger.kind));
            snap.orbBodyLabels.push_back(mergerKindLabel(s.merger.kind));
            snap.bhPosition += s.merger.primaryOffset();

            // Secondary black-hole visual identity (accretion disk / spin /
            // jets) so the inspiralling BH looks like its profile counterpart.
            const auto& v = s.merger.secVis;
            snap.secBHActive       = (s.merger.kind == MergerSecondaryKind3D::BlackHole) && v.hasDisk;
            snap.secBHDiskNormal   = v.diskNormal;
            snap.secBHSpin         = v.spin;
            snap.secBHDiskInner    = v.diskInnerRatio;
            snap.secBHDiskOuter    = v.diskOuterRatio;
            snap.secBHDiskStrength = v.diskStrength;
            snap.secBHColorInner   = v.colorInner;
            snap.secBHColorOuter   = v.colorOuter;
            snap.secBHShowJets     = v.hasJets;
            snap.secBHJetColor     = v.jetColor;
            snap.secBHJetRadius    = v.jetRadiusRatio;
            snap.secBHJetLength    = v.jetLengthRatio;
        } else {
            snap.secBHActive = false;
        }
    }
    // Force orbiting bodies on while the secondary is inspiralling, even if the
    // active profile normally hides them.
    snap.orbBodyEnabled    = s.orbBodyEnabled || (s.merger.active && !s.merger.merged);
    snap.hostGalaxyEnabled = s.hostGalaxyEnabled;
    snap.labEnabled        = s.labEnabled;
    snap.cgmEnabled        = s.cgmEnabled;
    snap.maxSteps          = s.cinematicMode ? 300 : 200;
    snap.profileName       = s.profilesArr[s.profileIdx].name;
    snap.windowW           = w;
    snap.windowH           = h;
    snap.fps               = s.fpsDisplay;

    // ---- Tidal disruption event data for HUD overlay ----
    snap.tdeActive     = s.tde3D.active;
    snap.tdeFlashAlpha = (s.tde3D.flashTimer > 0.0f)
        ? (s.tde3D.flashTimer / State::TidalEvent3D::FLASH_DURATION)
        : 0.0f;
    snap.tdeEventPos   = s.tde3D.eventPos;
    snap.tdeDebrisPos         = s.tde3D.debrisPos;
    snap.tdeDebrisIsFallback  = s.tde3D.isFallback;
    snap.tdeDebrisLifeF.clear();
    snap.tdeDebrisLifeF.reserve(s.tde3D.debrisLife.size());
    for (size_t i = 0; i < s.tde3D.debrisLife.size(); ++i) {
        float maxL = s.tde3D.debrisMaxLife[i];
        snap.tdeDebrisLifeF.push_back(maxL > 0.0f ? s.tde3D.debrisLife[i] / maxL : 0.0f);
    }

    // ---- Visual merger data for HUD overlay ----
    snap.mergerActive     = s.merger.active;
    snap.mergerInspiral   = s.merger.active && !s.merger.merged;
    snap.mergerFlashAlpha = (s.merger.flashTimer > 0.0f)
        ? std::min(1.0f, s.merger.flashTimer / State::MergerState3D::FLASH_DURATION)
        : 0.0f;
    snap.mergerEventPos       = s.merger.eventPos;
    snap.mergerRingActive     = s.merger.ringActive;
    snap.mergerRingT          = s.merger.ringT;
    snap.mergerSecondaryKind  = static_cast<int>(s.merger.kind);
    snap.mergerSecondaryLabel = mergerKindLabel(s.merger.kind);
    snap.mergerSecondaryPos   = s.merger.secondaryWorld(s.config.blackHole.position);
    snap.mergerSepRs          = (s.merger.bhR0 > 0.0f) ? (s.merger.r / s.merger.bhR0) : s.merger.r;
    snap.mergerRemnantAlpha   = (s.merger.remnantTimer > 0.0f)
        ? std::min(1.0f, s.merger.remnantTimer / State::MergerState3D::REMNANT_LABEL)
        : 0.0f;
    snap.mergerTrail.assign(s.merger.trail.begin(), s.merger.trail.end());
}

// ────────────────────────────────────────────────────────────
//  Render: scene + (optional) bloom + HUD overlay
// ────────────────────────────────────────────────────────────

// Render the scene. Caller must have an active GL context and a default
// framebuffer bound (this function will switch FBOs internally for bloom).
inline void renderScene(State& s, int w, int h) {
    if (w != s.lastW || h != s.lastH) {
        glViewport(0, 0, w, h);
        s.bloom.resize(w, h);
        s.lastW = w; s.lastH = h;
    }
    if (s.cinematicMode) { 
        s.bloom.sceneFBO.bind();
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        setSceneUniforms(*s.activeProgram, s.activeUfs, s.snap,
                         s.bgTex.id(), s.photorealDiskTex.id());
        s.quad.drawQuad();
        s.bloom.execute(s.quad, cfg::cinematicBloom(), true, s.totalTime);
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, w, h);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);
        setSceneUniforms(*s.activeProgram, s.activeUfs, s.snap,
                         s.bgTex.id(), s.simpleDiskTex.id());
        s.quad.drawQuad();
    }
    // bloom.execute() saves+restores prevFBO which ends up as sceneFBO; force
    // back to default before any HUD or display(), required on macOS Metal-GL.
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, w, h);

    // RK4 orbit + photon-ray overlays (world-space lines, drawn over the
    // ray-marched scene). Lazy-init the GL resources on first draw so
    // tickPhysics doesn't have to know about GL context state.
    if (s.physOverlay.orbitsEnabled || s.physOverlay.photonsEnabled ||
        s.physOverlay.spacetimeEnabled) {
        s.physOverlay.init();
        s.physOverlay.notifyScale(s.snap.bhRadius, s.snap.diskOuterRadius);
        if (s.physOverlay.dirty()) {
            s.physOverlay.rebuild(s.snap, s.orbBodies);
        }
        s.physOverlay.draw(s.snap);
    }
}

// Render the panel-based HUD + label-view overlay + overlays toggle panel.
// Updates s.hudFrame so subsequent left-clicks can be hit-tested.
inline void renderHUD(State& s, int w, int h) {
    const bool drawAnyHUD = ((s.showHUD || s.showDebugHUD) || s.overlays.labelView) && s.fontLoaded;
    if (!drawAnyHUD) return;

    for (int i = 0; i < 8; ++i) {
        glActiveTexture(GL_TEXTURE0 + i);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const float rollDeg = s.snap.roll * (180.f / 3.14159265f);
    s.hudFrame.clear();

    // 1) In-scene labels first so HUD panels render on top
    if (s.overlays.labelView && s.snap.orbBodyEnabled) {
        for (size_t i = 0; i < s.snap.orbBodyPositions.size(); ++i) {
            float sx = 0.f, sy = 0.f;
            if (hud::projectToScreen(s.snap.orbBodyPositions[i],
                                     s.snap.cameraPos, s.snap.cameraDir, s.snap.cameraUp,
                                     s.snap.fov, s.snap.windowW, s.snap.windowH, sx, sy)) {
                if (sx < -200.f || sx > s.snap.windowW + 200.f) continue;
                if (sy < -200.f || sy > s.snap.windowH + 200.f) continue;
                const std::string& nm = (i < s.snap.orbBodyLabels.size())
                    ? s.snap.orbBodyLabels[i] : std::string("Body");
                hud::drawBodyLabel(s.glFont, nm, sx, sy, s.snap.windowW, s.snap.windowH);
            }
        }

        // Barycenter crosshair for Gaia binary systems
        if (s.snap.barycentricMode) {
            float bsx = 0.f, bsy = 0.f;
            if (hud::projectToScreen(s.snap.barycenterPos,
                                     s.snap.cameraPos, s.snap.cameraDir, s.snap.cameraUp,
                                     s.snap.fov, s.snap.windowW, s.snap.windowH, bsx, bsy)) {
                hud::drawBarycenterMarker(s.glFont, bsx, bsy,
                                         s.snap.windowW, s.snap.windowH);
            }
        }
    }

    if (s.showHUD) {
        if (!s.useImGuiHud) {
            hud::drawStatusPanel(s.glFont, s.snap, rollDeg, s.snap.windowW, s.snap.windowH);
            hud::drawControlsHint(s.glFont, s.snap.windowW, s.snap.windowH);
            if (s.overlays.panelOpen) {
                hud::drawOverlaysPanel(s.glFont, s.hudFrame, s.overlays,
                                       s.snap.windowW, s.snap.windowH);
            }
        }
    }
    if (s.showDebugHUD && !s.useImGuiHud) {
        hud::DebugPanelInfo dbg{};
        dbg.cinematic    = s.snap.cinematicMode;
        dbg.freelook     = s.snap.freelook;
        dbg.jets         = s.snap.jetsEnabled;
        dbg.blr          = s.snap.blrEnabled;
        dbg.doppler      = s.snap.dopplerEnabled;
        dbg.havePhotoreal= s.havePhotoreal;
        dbg.haveSimple   = s.haveSimple;
        dbg.animSpeed    = s.snap.animSpeed;
        dbg.fps          = s.snap.fps;
        dbg.totalTime    = s.snap.totalTime;
        dbg.animTime     = s.snap.animTime;
        dbg.maxSteps     = s.snap.maxSteps;
        dbg.camPos       = s.snap.cameraPos;
        dbg.camDir       = s.snap.cameraDir;
        dbg.rollDeg      = rollDeg;
        hud::drawDebugPanel(s.glFont, dbg, s.snap.windowW, s.snap.windowH);
    }

    // Preset menu draws on top of all other HUD elements (modal). It also
    // renders during fade-out (presetMenu.t > 0 even after open=false).
    // Skipped entirely when an ImGui-based HUD is active in the host.
    if (!s.useImGuiHud && s.presetMenu.t > 0.001f) {
        // Build parallel name/description arrays from the loaded profiles.
        std::array<std::string, profiles::NUM_PROFILES> names;
        std::array<std::string, profiles::NUM_PROFILES> descs;
        for (int i = 0; i < profiles::NUM_PROFILES; ++i) {
            names[i] = s.profilesArr[i].name;
            descs[i] = s.profilesArr[i].description;
        }
        hud::drawPresetMenu(s.glFont, s.hudFrame, s.presetMenu,
                            names.data(), descs.data(),
                            profiles::NUM_PROFILES, s.profileIdx,
                            s.snap.windowW, s.snap.windowH);
    }

    glDisable(GL_BLEND);
    (void)w; (void)h;
}

// ────────────────────────────────────────────────────────────
//  Profile switching
// ────────────────────────────────────────────────────────────
inline void switchToProfile(State& s, int newIdx) {
    s.profileIdx = ((newIdx % profiles::NUM_PROFILES) + profiles::NUM_PROFILES) % profiles::NUM_PROFILES;
    const auto& prof = s.profilesArr[s.profileIdx];
    s.config            = prof.config;
    s.jetsEnabled       = prof.defaultJets;
    s.blrEnabled        = prof.defaultBLR;
    s.orbBodyEnabled    = prof.defaultOrbBody;
    s.dopplerEnabled    = prof.defaultDoppler;
    s.hostGalaxyEnabled = prof.defaultHostGalaxy;
    s.labEnabled        = prof.defaultLAB;
    s.cgmEnabled        = prof.defaultCGM;
    s.camera            = CameraController(prof.config.camera);
    rebuildOrbBodies(s.orbBodies, prof, s.config);
    s.orbBodyDisrupted.assign(s.orbBodies.size(), false);
    s.tde3D = State::TidalEvent3D{};
    // Cancel any in-flight merger (and its grown remnant) on profile switch.
    s.merger      = State::MergerState3D{};
    s.mergerMenu  = State::MergerMenuState{};
    s.physOverlay.markDirty();
    std::cerr << "[bh3d] profile: " << prof.name
              << ", " << prof.description << "\n";
}

// ────────────────────────────────────────────────────────────
//  Input dispatch
// ────────────────────────────────────────────────────────────

// Dispatch one action keypress. Movement keys (WASD/QE/Space/Shift/Ctrl)
// should already have been forwarded to s.keys.onKeyPressed by the caller.
inline void onActionKey(State& s, sf::Keyboard::Key code) {
    const auto& kb = s.actionKeys;

    // ── Preset menu intercepts navigation when open ──
    if (s.presetMenu.open) {
        if (s.useImGuiHud) {
            // ImGui owns Up/Down/Enter/Esc, we only handle the toggle key
            // here so the user can still close the menu with N.
            if (code == kb.nextProfile) {
                s.presetMenu.open = false;
                return;
            }
            // All other keys fall through (camera/HUD toggles still work).
        } else {
            if (code == sf::Keyboard::Key::Up) {
                s.presetMenu.selected = std::max(0, s.presetMenu.selected - 1);
                return;
            }
            if (code == sf::Keyboard::Key::Down) {
                s.presetMenu.selected = std::min(profiles::NUM_PROFILES - 1,
                                                 s.presetMenu.selected + 1);
                return;
            }
            if (code == sf::Keyboard::Key::Enter) {
                switchToProfile(s, s.presetMenu.selected);
                s.presetMenu.open = false;
                return;
            }
            if (code == sf::Keyboard::Key::Escape || code == kb.nextProfile) {
                s.presetMenu.open = false;
                return;
            }
            // Other keys fall through (so e.g. M still toggles overlays etc.)
        }
    }

    if      (code == kb.toggleFreelook)    s.camera.toggleMode();
    else if (code == kb.toggleJets)        s.jetsEnabled       = !s.jetsEnabled;
    else if (code == kb.toggleBLR)         s.blrEnabled        = !s.blrEnabled;
    else if (code == kb.toggleOrbBody)     s.orbBodyEnabled    = !s.orbBodyEnabled;
    else if (code == kb.toggleDoppler)     s.dopplerEnabled    = !s.dopplerEnabled;
    else if (code == kb.toggleBlueshift)   s.blueshiftEnabled  = !s.blueshiftEnabled;
    else if (code == kb.toggleHostGalaxy)  s.hostGalaxyEnabled = !s.hostGalaxyEnabled;
    else if (code == kb.toggleLAB)         s.labEnabled        = !s.labEnabled;
    else if (code == kb.toggleCGM)         s.cgmEnabled        = !s.cgmEnabled;
    else if (code == kb.cycleAnimSpeed) {
        s.animSpeedIdx = (s.animSpeedIdx + 1) % State::animSpeedCount;
        s.animSpeed    = State::animSpeedPresets[s.animSpeedIdx];
    }
    else if (code == kb.speedUp) {
        // Mirror 2D: continuous doubling with caps. Decouple from preset cycle.
        s.animSpeed = std::min(64.0f, s.animSpeed * 2.0f);
    }
    else if (code == kb.speedDown) {
        s.animSpeed = std::max(0.125f, s.animSpeed * 0.5f);
    }
    else if (code == kb.toggleHUD)         s.showHUD      = !s.showHUD;
    else if (code == kb.toggleDebugHUD)    s.showDebugHUD = !s.showDebugHUD;
    else if (code == kb.toggleOverlays)    s.overlays.panelOpen = !s.overlays.panelOpen;
    else if (code == kb.toggleRK4Orbits) {
        s.physOverlay.orbitsEnabled = !s.physOverlay.orbitsEnabled;
        s.physOverlay.markDirty();
    }
    else if (code == kb.toggleRK4Photons) {
        s.physOverlay.photonsEnabled = !s.physOverlay.photonsEnabled;
        s.physOverlay.markDirty();
    }
    else if (code == kb.toggleSpacetime) {
        s.physOverlay.spacetimeEnabled = !s.physOverlay.spacetimeEnabled;
        s.physOverlay.markDirty();
    }
    else if (code == kb.resetTilt)         s.camera.resetRoll();
    else if (code == kb.nextProfile) {
        // Open the preset menu (or close it if already open). The menu-open
        // branch above also handles toggling-off via this same key.
        s.presetMenu.open     = true;
        s.presetMenu.selected = s.profileIdx;
    }
    else if (code == kb.openMergerMenu) {
        // Toggle the merger selector. Opening it closes the preset menu so the
        // two modal panels never overlap.
        s.mergerMenu.open = !s.mergerMenu.open;
        if (s.mergerMenu.open) s.presetMenu.open = false;
    }
    else if (code == kb.toggleRender) {
        if (s.havePhotoreal && s.haveSimple) {
            s.cinematicMode  = !s.cinematicMode;
            s.activeProgram  = s.cinematicMode ? &s.photorealProgram : &s.simpleProgram;
            s.activeUfs      = s.cinematicMode ? s.photorealUfs      : s.simpleUfs;
        }
    }
    else if (code == kb.releaseMouse)      s.looking = false;
}

// Hit-test a left click against the open overlays panel. Returns true if
// the click hit a button (and the matching overlay was toggled).
inline bool onLeftClickOverlays(State& s, float x, float y) {
    // Preset menu (modal) takes priority while it is open.
    if (s.presetMenu.open) {
        const int hit = s.hudFrame.hitTestPresetMenu(x, y);
        if (hit == -2) {
            s.presetMenu.scroll = std::max(0, s.presetMenu.scroll - 1);
            return true;
        }
        if (hit == -3) {
            s.presetMenu.scroll = s.presetMenu.scroll + 1; // clamped on next draw
            return true;
        }
        if (hit >= 0) {
            switchToProfile(s, hit);
            s.presetMenu.open = false;
            return true;
        }
        // Click outside any item closes the menu.
        s.presetMenu.open = false;
        return true;
    }

    if (!(s.overlays.panelOpen && s.showHUD)) return false;
    const hud::HudButtonId id = s.hudFrame.hitTest(x, y);
    switch (id) {
        case hud::BTN_LABEL_VIEW:    s.overlays.labelView       = !s.overlays.labelView;       return true;
        case hud::BTN_SPACETIME_VIZ: s.overlays.spacetimeViz    = !s.overlays.spacetimeViz;    return true;
        case hud::BTN_GRAV_REDSHIFT: s.overlays.gravRedshift    = !s.overlays.gravRedshift;    return true;
        case hud::BTN_KERR_SHADOW:   s.overlays.kerrShadowGuide = !s.overlays.kerrShadowGuide; return true;
        default: break;
    }
    return false;
}

} // namespace bh3d