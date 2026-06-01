// ============================================================
// bh3d_core.hpp
// ============================================================
// Single source of truth for the 3D black-hole simulation.
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
#include "shader_sources.hpp"
#include "shader_utils.hpp"
#include "texture_utils.hpp"
#include "bloom_pipeline.hpp"
#include "camera_controller.hpp"
#include "input.hpp"
#include "simulation_state.hpp"
#include "presets.hpp"
#include "orbital_body.hpp"
#include "gl_font.hpp"
#include "hud_panel.hpp"
#include "physics_overlay.hpp"

#include <array>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace bh3d {

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
    std::array<Entry, 22> entries = {{
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
        {"toggle_spacetime",  &keys.toggleSpacetime,  defaults.toggleSpacetime}
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
        {"toggle_spacetime",  &keys.toggleSpacetime}
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
//  Per-frame physics + snapshot construction
// ────────────────────────────────────────────────────────────
inline void tickPhysics(State& s, float dt) {
    s.totalTime += dt;
    s.animTime  += dt * s.animSpeed;
    s.camera.update(dt, s.keys);
    // Orbiting bodies share the global animSpeed multiplier so +/- (and Y)
    // also speed/slow Keplerian orbits, matching the 2D simulator's behavior.
    const float orbDt = dt * s.animSpeed;
    for (auto& body : s.orbBodies) body.update(orbDt);

    // ---- Tidal disruption detection ----
    // NOTE: All 3D orbBodies are on permanent Keplerian orbits (no energy loss /
    // radiation-reaction inspiral). Triggering a TDE when a high-eccentricity body
    // passes periapsis would fire every single orbit, physically nonsensical for
    // stable S-stars, pulsars, etc. TDE in 3D is reserved for bodies that are
    // genuinely accreting (future feature). Detection loop intentionally disabled.
    {

        // Advance existing TDE particles
        if (s.tde3D.active) {
            if (s.tde3D.flashTimer > 0.0f)
                s.tde3D.flashTimer -= dt;

            for (size_t i = 0; i < s.tde3D.debrisPos.size(); ) {
                s.tde3D.debrisLife[i] -= dt;
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
                    s.tde3D.debrisPos[i] += s.tde3D.debrisVel[i] * dt;
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
}

inline void tickFPS(State& s) { // much better than executing the previous value by firing squad once every frame
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
    }
    snap.orbBodyEnabled    = s.orbBodyEnabled;
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

} // namespace stuff