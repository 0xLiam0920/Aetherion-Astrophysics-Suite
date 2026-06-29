#pragma once
// ============================================================
// hud_panel.hpp: Panel-based HUD, overlays menu, body labels
//
// Reuses GLBitmapFont for both text and solid-rect drawing.
// All public functions are stateless; the caller passes in
// state structs and collects button rects for hit-testing.
// ============================================================

#include "gl_font.hpp"
#include "bh3d_presets.hpp"
#include "bh3d_simulationstates.hpp"

#include <glm/glm.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>

namespace hud {

/*--------- Overlay menu state (persistent across frames) ---------*/
struct OverlaysState {
    bool panelOpen        = false; // toggled by 'M' key
    bool labelView        = false; // FUNCTIONAL, draws name tags over orbiting bodies
    bool spacetimeViz     = false; // placeholder
    bool gravRedshift     = false; // placeholder
    bool ISCOSpinParam    = false; // placeholder until I can actually connect it to something
    bool FluxLuminosity   = false; // placeholder
    bool kerrShadowGuide  = false; // placeholder
};

/*--------- Per-frame button hit-test record ---------*/
enum HudButtonId {
    BTN_NONE = 0,
    BTN_LABEL_VIEW,
    BTN_SPACETIME_VIZ,
    BTN_GRAV_REDSHIFT,
    BTN_KERR_SHADOW,
};

struct HudButtonRect {
    float x, y, w, h;
    HudButtonId id;
};

/*--------- Preset menu hit area (id is the profile index) ---------*/
struct PresetItemRect {
    float x, y, w, h;
    int   profileIdx;
};

/*--------- Preset menu state (open/closed + selection + animation) ---------*/
struct PresetMenuState {
    bool  open       = false;
    int   selected   = 0;
    int   scroll     = 0;
    float t          = 0.0f;       // 0..1 transition progress
    static constexpr int   visibleCount = 8;
    static constexpr float transitionSec = 0.18f;
};

struct HudFrame {
    std::vector<HudButtonRect>  buttons;
    std::vector<PresetItemRect> presetItems;
    PresetItemRect              presetScrollUp{0,0,0,0,-1};
    PresetItemRect              presetScrollDown{0,0,0,0,-1};
    bool                        presetScrollUpVisible   = false;
    bool                        presetScrollDownVisible = false;
    void clear() {
        buttons.clear();
        presetItems.clear();
        presetScrollUpVisible   = false;
        presetScrollDownVisible = false;
    }
    HudButtonId hitTest(float mx, float my) const {
        for (const auto& b : buttons) {
            if (mx >= b.x && mx <= b.x + b.w && my >= b.y && my <= b.y + b.h) return b.id;
        }
        return BTN_NONE;
    }
    // Returns profile index, -2 for scroll-up, -3 for scroll-down, -1 for miss.
    int hitTestPresetMenu(float mx, float my) const {
        if (presetScrollUpVisible &&
            mx >= presetScrollUp.x && mx <= presetScrollUp.x + presetScrollUp.w &&
            my >= presetScrollUp.y && my <= presetScrollUp.y + presetScrollUp.h)
            return -2;
        if (presetScrollDownVisible &&
            mx >= presetScrollDown.x && mx <= presetScrollDown.x + presetScrollDown.w &&
            my >= presetScrollDown.y && my <= presetScrollDown.y + presetScrollDown.h)
            return -3;
        for (const auto& it : presetItems) {
            if (mx >= it.x && mx <= it.x + it.w && my >= it.y && my <= it.y + it.h)
                return it.profileIdx;
        }
        return -1;
    }
};

/*--------- Color palette (RGBA, 0..1) ---------*/
struct Color { float r, g, b, a; };
namespace col {
    constexpr Color panelBg     = { 0.04f, 0.06f, 0.10f, 0.82f };
    constexpr Color panelBorder = { 0.35f, 0.55f, 0.85f, 0.95f };
    constexpr Color titleBar    = { 0.10f, 0.16f, 0.28f, 0.95f };
    constexpr Color divider     = { 0.20f, 0.28f, 0.42f, 0.85f };
    constexpr Color textPrimary = { 0.96f, 0.97f, 0.99f, 1.00f };
    constexpr Color textDim     = { 0.70f, 0.74f, 0.82f, 1.00f };
    constexpr Color textAccent  = { 0.55f, 0.85f, 1.00f, 1.00f };
    constexpr Color valueOn     = { 0.55f, 0.95f, 0.55f, 1.00f };
    constexpr Color valueOff    = { 0.85f, 0.55f, 0.55f, 1.00f };
    constexpr Color btnBg       = { 0.08f, 0.12f, 0.20f, 0.95f };
    constexpr Color btnBgActive = { 0.18f, 0.45f, 0.30f, 0.95f };
    constexpr Color btnBorder   = { 0.40f, 0.55f, 0.80f, 0.95f };
    constexpr Color btnPlaceholder = { 0.10f, 0.10f, 0.12f, 0.85f };
    constexpr Color labelBg     = { 1.00f, 1.00f, 1.00f, 0.92f };
    constexpr Color labelText   = { 0.05f, 0.05f, 0.05f, 1.00f };
    constexpr Color labelLine   = { 1.00f, 1.00f, 1.00f, 0.85f };
    constexpr Color labelDot    = { 1.00f, 1.00f, 1.00f, 1.00f };
}

/*--------- Primitive: 1px-style border via 4 thin rects ---------*/
inline void drawBorder(const GLBitmapFont& font, float x, float y, float w, float h,
                       Color c, int sw, int sh, float t = 1.0f) {
    font.drawRect(x,         y,         w, t,    c.r, c.g, c.b, c.a, sw, sh);
    font.drawRect(x,         y + h - t, w, t,    c.r, c.g, c.b, c.a, sw, sh);
    font.drawRect(x,         y,         t, h,    c.r, c.g, c.b, c.a, sw, sh);
    font.drawRect(x + w - t, y,         t, h,    c.r, c.g, c.b, c.a, sw, sh);
}

inline void drawFilled(const GLBitmapFont& font, float x, float y, float w, float h,
                       Color c, int sw, int sh) {
    font.drawRect(x, y, w, h, c.r, c.g, c.b, c.a, sw, sh);
}

inline void drawString(const GLBitmapFont& font, const std::string& s, float x, float y,
                       Color c, int sw, int sh) {
    font.drawText(s, x, y, c.r, c.g, c.b, c.a, sw, sh);
}

/*--------- Panel: titled rectangle with bg + border + title bar ---------*/
inline void drawPanel(const GLBitmapFont& font, float x, float y, float w, float h,
                      const std::string& title, int sw, int sh) {
    drawFilled(font, x, y, w, h, col::panelBg, sw, sh);
    const float titleH = 22.0f;
    drawFilled(font, x, y, w, titleH, col::titleBar, sw, sh);
    drawString(font, title, x + 8.0f, y + 4.0f, col::textAccent, sw, sh);
    drawBorder(font, x, y, w, h, col::panelBorder, sw, sh, 1.0f);
    // Divider under title
    drawFilled(font, x, y + titleH, w, 1.0f, col::divider, sw, sh);
}

/*--------- Key/value row inside a panel ---------*/
inline void drawKVRow(const GLBitmapFont& font, float x, float y, float w,
                      const std::string& key, const std::string& value,
                      Color valueColor, int sw, int sh) {
    drawString(font, key, x, y, col::textDim, sw, sh);
    float vw = font.textWidth(value);
    drawString(font, value, x + w - vw - 4.0f, y, valueColor, sw, sh);
}

/*--------- Toggle button (rectangle + label + state pill) ---------*/
inline void drawToggleButton(const GLBitmapFont& font, HudFrame& frame,
                             HudButtonId id, float x, float y, float w, float h,
                             const std::string& label, bool active, bool placeholder,
                             int sw, int sh) {
    Color bg = placeholder ? col::btnPlaceholder
                           : (active ? col::btnBgActive : col::btnBg);
    drawFilled(font, x, y, w, h, bg, sw, sh);
    drawBorder(font, x, y, w, h, col::btnBorder, sw, sh, 1.0f);

    // Status pill on the left
    const float pillSize = 10.0f;
    const float pillX = x + 6.0f;
    const float pillY = y + (h - pillSize) * 0.5f;
    Color pillCol = placeholder ? Color{0.45f, 0.45f, 0.50f, 1.0f}
                                : (active ? col::valueOn : col::valueOff);
    drawFilled(font, pillX, pillY, pillSize, pillSize, pillCol, sw, sh);

    // Label text
    drawString(font, label, x + 6.0f + pillSize + 6.0f, y + 4.0f,
               placeholder ? col::textDim : col::textPrimary, sw, sh);

    // "(soon)" suffix for placeholder buttons
    if (placeholder) {
        const std::string suf = "(soon)";
        float sw2 = font.textWidth(suf);
        drawString(font, suf, x + w - sw2 - 6.0f, y + 4.0f, col::textDim, sw, sh);
    }

    frame.buttons.push_back({x, y, w, h, id});
}

/*--------- Default label fallback per body type ---------*/
inline std::string defaultLabelForType(int bodyTypeEnum) {
    switch (static_cast<GalaxyBody3DType>(bodyTypeEnum)) {
        case GalaxyBody3DType::Star:           return "Star";
        case GalaxyBody3DType::GasCloud:       return "Gas Cloud";
        case GalaxyBody3DType::StellarCluster: return "Stellar Cluster";
        case GalaxyBody3DType::DwarfGalaxy:    return "Dwarf Galaxy";
        case GalaxyBody3DType::NeutronStar:    return "Neutron Star";
        case GalaxyBody3DType::WhiteDwarf:     return "White Dwarf";
        case GalaxyBody3DType::CompanionStar:  return "Companion Star";
    }
    return "Body";
}

/*--------- 3D world position -> screen pixel coords ---------*/
// Returns true if the point is in front of the camera AND inside the viewport.
inline bool projectToScreen(const glm::vec3& worldPos,
                            const glm::vec3& camPos,
                            const glm::vec3& camDir,
                            const glm::vec3& camUp,
                            float fovDeg, int sw, int sh,
                            float& outX, float& outY) {
    glm::vec3 forward = glm::normalize(camDir);
    glm::vec3 rightV  = glm::normalize(glm::cross(forward, glm::normalize(camUp)));
    glm::vec3 upV     = glm::normalize(glm::cross(rightV, forward));
    glm::vec3 rel     = worldPos - camPos;
    float zCam = glm::dot(rel, forward);
    if (zCam < 0.05f) return false; // behind camera or too close
    float xCam = glm::dot(rel, rightV);
    float yCam = glm::dot(rel, upV);

    const float aspect = (sh > 0) ? float(sw) / float(sh) : 1.0f;
    const float halfFov = 0.5f * fovDeg * 3.14159265f / 180.0f;
    const float tanHF = std::tan(halfFov);
    if (tanHF < 1e-6f) return false;

    float ndcX = (xCam / zCam) / (tanHF * aspect);
    float ndcY = (yCam / zCam) / tanHF;

    outX = (ndcX * 0.5f + 0.5f) * float(sw);
    outY = (1.0f - (ndcY * 0.5f + 0.5f)) * float(sh);
    return true;
}

/*--------- Barycenter crosshair (Gaia binary systems) ---------*/
// Draws a small "+" crosshair with a "Center of Mass" label at the given
// screen position.  Called only when snap.barycentricMode is true.
inline void drawBarycenterMarker(const GLBitmapFont& font,
                                 float sx, float sy,
                                 int sw, int sh) {
    const float armLen  = 7.0f;
    const float armThk  = 1.5f;
    const Color comc    = {0.78f, 1.00f, 0.78f, 0.80f}; // pale green
    const Color textCol = {0.78f, 1.00f, 0.78f, 0.90f};

    // Horizontal arm
    drawFilled(font, sx - armLen, sy - armThk * 0.5f, armLen * 2.0f, armThk, comc, sw, sh);
    // Vertical arm
    drawFilled(font, sx - armThk * 0.5f, sy - armLen, armThk, armLen * 2.0f, comc, sw, sh);

    // "Center of Mass" label to the right
    const std::string lbl = "Center of Mass";
    const float lx = sx + armLen + 4.0f;
    const float ly = sy - float(font.lineHeight()) * 0.5f;
    drawString(font, lbl, lx, ly, textCol, sw, sh);
}

/*--------- Westbrook-style label tag with leader line ---------*/
inline void drawBodyLabel(const GLBitmapFont& font, const std::string& text,
                          float bodyScreenX, float bodyScreenY,
                          int sw, int sh) {
    // Anchor dot on the body
    const float dotR = 3.0f;
    drawFilled(font, bodyScreenX - dotR, bodyScreenY - dotR, dotR * 2.0f, dotR * 2.0f,
               col::labelDot, sw, sh);

    // Label box hovers above the body. Leader line drops vertically.
    const float padX = 6.0f, padY = 3.0f;
    const float labelW = font.textWidth(text) + padX * 2.0f;
    const float labelH = float(font.lineHeight()) + padY * 2.0f - 2.0f;
    const float leaderLen = 38.0f;

    float boxX = bodyScreenX - labelW * 0.5f;
    float boxY = bodyScreenY - leaderLen - labelH;

    // Clamp box inside viewport
    if (boxX < 4.0f) boxX = 4.0f;
    if (boxX + labelW > sw - 4.0f) boxX = sw - 4.0f - labelW;
    if (boxY < 4.0f) {
        // If above goes off-screen, place below the body instead
        boxY = bodyScreenY + leaderLen;
    }

    // Leader line: 1px vertical from body up to the bottom-center of the box
    float lineX = bodyScreenX;
    float lineY0 = std::min(boxY + labelH, bodyScreenY);
    float lineY1 = std::max(boxY, bodyScreenY);
    drawFilled(font, lineX - 0.5f, lineY0, 1.0f, lineY1 - lineY0,
               col::labelLine, sw, sh);

    // Label background + thin border
    drawFilled(font, boxX, boxY, labelW, labelH, col::labelBg, sw, sh);
    drawBorder(font, boxX, boxY, labelW, labelH,
               Color{0.85f, 0.85f, 0.85f, 0.95f}, sw, sh, 1.0f);
    drawString(font, text, boxX + padX, boxY + padY - 1.0f, col::labelText, sw, sh);
}

/*--------- Status HUD: profile + camera + toggles ---------*/
inline void drawStatusPanel(const GLBitmapFont& font,
                            const PhysicsSnapshot& snap,
                            float rollDeg, int sw, int sh) {
    const float x = 12.0f, y = 12.0f;
    const float w = 280.0f;
    const float titleH = 22.0f;
    const float rowH   = 16.0f;
    const float pad    = 8.0f;

    // Section sizes
    int profileLines = 1; // profile row
    int cameraLines  = 2; // mode + tilt
    int toggleLines  = 9; // jets, blr, orbBody, doppler, blueshift, host, lab, cgm, anim
    int perfLines    = 2; // fps + render
    float panelH = titleH + pad
                 + profileLines * rowH + pad + 1
                 + cameraLines  * rowH + pad + 1
                 + toggleLines  * rowH + pad + 1
                 + perfLines    * rowH + pad;

    drawPanel(font, x, y, w, panelH, "Status", sw, sh);

    float cy = y + titleH + pad - 2.0f;

    auto section = [&](const std::string& heading) {
        drawString(font, heading, x + pad, cy, col::textAccent, sw, sh);
        cy += rowH;
    };
    auto kv = [&](const std::string& k, const std::string& v, Color vc) {
        drawKVRow(font, x + pad, cy, w - pad * 2.0f, k, v, vc, sw, sh);
        cy += rowH;
    };
    auto onoff = [&](const std::string& k, bool on) {
        kv(k, on ? "ON" : "OFF", on ? col::valueOn : col::valueOff);
    };
    auto divider = [&]() {
        drawFilled(font, x + 6.0f, cy + 1.0f, w - 12.0f, 1.0f, col::divider, sw, sh);
        cy += pad;
    };

    // Profile
    kv("Profile", snap.profileName, col::textPrimary);
    divider();

    // Camera
    onoff("Freelook", snap.freelook);
    {
        char buf[32]; std::snprintf(buf, sizeof(buf), "%d deg", int(rollDeg + 0.5f));
        kv("Tilt", buf, col::textPrimary);
    }
    divider();

    // Toggles
    onoff("Jets",        snap.jetsEnabled);
    onoff("BLR",         snap.blrEnabled);
    onoff("Orb Bodies",  snap.orbBodyEnabled);
    onoff("Doppler",     snap.dopplerEnabled);
    onoff("Blueshift",   snap.blueshiftEnabled);
    onoff("Host Galaxy", snap.hostGalaxyEnabled);
    onoff("LAB",         snap.labEnabled);
    onoff("CGM",         snap.cgmEnabled);
    {
        char buf[32];
        if (snap.animSpeed >= 1.0f) std::snprintf(buf, sizeof(buf), "%dx", int(snap.animSpeed + 0.5f));
        else                       std::snprintf(buf, sizeof(buf), "%.2fx", snap.animSpeed);
        kv("Anim Speed", buf, col::textPrimary);
    }
    divider();

    // Perf
    kv("Render", snap.cinematicMode ? "CINEMATIC" : "FAST",
       snap.cinematicMode ? col::valueOn : col::textPrimary);
    {
        char buf[32]; std::snprintf(buf, sizeof(buf), "%d", int(snap.fps + 0.5f));
        kv("FPS", buf, col::textPrimary);
    }
}

/*--------- Controls hint strip (bottom of screen) ---------*/
inline void drawControlsHint(const GLBitmapFont& font, int sw, int sh) {
    const std::string txt =
        "[WASD] Move  [Q/E] Tilt  [Space/LCTRL] Up/Down  [LSHIFT] Fast  "
        "[RMB] Look  [N] Next Profile  [H] Hide HUD  [B] Debug  [M] Overlays";
    float tw = font.textWidth(txt);
    float bh = float(font.lineHeight()) + 8.0f;
    float bx = (sw - tw) * 0.5f - 10.0f;
    float by = sh - bh - 10.0f;
    drawFilled(font, bx, by, tw + 20.0f, bh, col::panelBg, sw, sh);
    drawBorder(font, bx, by, tw + 20.0f, bh, col::panelBorder, sw, sh, 1.0f);
    drawString(font, txt, bx + 10.0f, by + 4.0f, col::textPrimary, sw, sh);
}

/*--------- Debug panel (right side) ---------*/
struct DebugPanelInfo {
    bool cinematic, freelook, jets, blr, doppler;
    bool havePhotoreal, haveSimple;
    float animSpeed, fps, totalTime, animTime, rollDeg;
    int maxSteps;
    glm::vec3 camPos, camDir;
};

inline void drawDebugPanel(const GLBitmapFont& font, const DebugPanelInfo& d, int sw, int sh) {
    const float w = 280.0f;
    const float titleH = 22.0f;
    const float rowH = 16.0f;
    const float pad = 8.0f;

    // Layout: 5 sections (Shader/Camera/Render/Animation/Window) with rows
    int rows = 1 /*shader*/ + 2 /*compiled*/
             + 1 /*mode*/   + 3 /*pos/dir/tilt*/
             + 1 /*render*/ + 2 /*steps/cinematic*/
             + 1 /*anim*/   + 3 /*speed/animT/totalT*/
             + 1 /*window*/ + 1 /*size*/;
    int dividers = 4;
    float panelH = titleH + pad + rows * rowH + dividers * pad + dividers * 1.0f + pad;

    const float x = sw - w - 12.0f;
    const float y = 12.0f;
    drawPanel(font, x, y, w, panelH, "Debug", sw, sh);

    float cy = y + titleH + pad - 2.0f;
    auto kv = [&](const std::string& k, const std::string& v, Color vc) {
        drawKVRow(font, x + pad, cy, w - pad * 2.0f, k, v, vc, sw, sh);
        cy += rowH;
    };
    auto sect = [&](const std::string& s) {
        drawString(font, s, x + pad, cy, col::textAccent, sw, sh);
        cy += rowH;
    };
    auto divider = [&]() {
        drawFilled(font, x + 6.0f, cy + 1.0f, w - 12.0f, 1.0f, col::divider, sw, sh);
        cy += pad;
    };
    char buf[64];

    sect("Shader");
    kv("Photoreal", d.havePhotoreal ? "compiled" : "N/A",
       d.havePhotoreal ? col::valueOn : col::valueOff);
    kv("Simple", d.haveSimple ? "compiled" : "N/A",
       d.haveSimple ? col::valueOn : col::valueOff);
    divider();

    sect("Camera");
    std::snprintf(buf, sizeof(buf), "(%.1f, %.1f, %.1f)", d.camPos.x, d.camPos.y, d.camPos.z);
    kv("Pos", buf, col::textPrimary);
    std::snprintf(buf, sizeof(buf), "(%.2f, %.2f, %.2f)", d.camDir.x, d.camDir.y, d.camDir.z);
    kv("Dir", buf, col::textPrimary);
    std::snprintf(buf, sizeof(buf), "%.1f deg", d.rollDeg);
    kv("Tilt", buf, col::textPrimary);
    divider();

    sect("Render");
    std::snprintf(buf, sizeof(buf), "%d", d.maxSteps);
    kv("Max Steps", buf, col::textPrimary);
    kv("Cinematic", d.cinematic ? "ON" : "OFF",
       d.cinematic ? col::valueOn : col::valueOff);
    divider();

    sect("Animation");
    if (d.animSpeed >= 1.0f) std::snprintf(buf, sizeof(buf), "%.0fx", d.animSpeed);
    else                     std::snprintf(buf, sizeof(buf), "%.2fx", d.animSpeed);
    kv("Speed", buf, col::textPrimary);
    std::snprintf(buf, sizeof(buf), "%.1fs", d.animTime);
    kv("Anim Time", buf, col::textPrimary);
    std::snprintf(buf, sizeof(buf), "%.1fs", d.totalTime);
    kv("Total Time", buf, col::textPrimary);
    divider();

    sect("Window");
    std::snprintf(buf, sizeof(buf), "%dx%d  @ %.0f FPS", sw, sh, d.fps);
    kv("Size", buf, col::textPrimary);
}

/*--------- Overlays panel: placeholder feature toggles ---------*/
inline void drawOverlaysPanel(const GLBitmapFont& font, HudFrame& frame,
                              const OverlaysState& ov, int sw, int sh) {
    const float w = 260.0f;
    const float titleH = 22.0f;
    const float pad    = 8.0f;
    const float btnH   = 26.0f;
    const float btnGap = 6.0f;

    const int nButtons = 4;
    float panelH = titleH + pad + nButtons * btnH + (nButtons - 1) * btnGap + pad
                 + 14.0f /*hint row*/ + pad;

    const float x = 12.0f;
    const float y = sh - panelH - 12.0f - 30.0f; // leave room for the bottom controls strip
    drawPanel(font, x, y, w, panelH, "Overlays", sw, sh);

    float cy = y + titleH + pad;
    const float bx = x + pad;
    const float bw = w - pad * 2.0f;

    drawToggleButton(font, frame, BTN_LABEL_VIEW,    bx, cy, bw, btnH,
                     "Label View",            ov.labelView,       false, sw, sh);
    cy += btnH + btnGap;
    drawToggleButton(font, frame, BTN_SPACETIME_VIZ, bx, cy, bw, btnH,
                     "Spacetime Visualization", ov.spacetimeViz,  true,  sw, sh);
    cy += btnH + btnGap;
    drawToggleButton(font, frame, BTN_GRAV_REDSHIFT, bx, cy, bw, btnH,
                     "Gravitational Redshift",  ov.gravRedshift,  true,  sw, sh);
    cy += btnH + btnGap;
    drawToggleButton(font, frame, BTN_KERR_SHADOW,   bx, cy, bw, btnH,
                     "Kerr Shadow Guide",       ov.kerrShadowGuide, true, sw, sh);
    cy += btnH + btnGap;

    drawString(font, "Click a button to toggle.", bx, cy + 2.0f, col::textDim, sw, sh);
}

/*--------- Preset menu: modal sleek list with smooth transition ---------
 * Renders a centered "Available Presets" panel with scroll arrows, item
 * highlight, and a fading dim backdrop. The caller is responsible for
 * advancing `state.t` (0 closed → 1 open) per frame. When `state.t <= 0`
 * the function early-outs without drawing or registering hit areas.
 *
 * `profileNames` and `profileDescriptions` must have size == numProfiles.
 * `currentIdx` is the profile that is presently active (drawn with an
 * accent dot). `state.selected` is the highlighted item the user is
 * navigating with arrow keys.
 *
 * After the call, `frame.presetItems`, `presetScrollUp`/`Down` are
 * populated for hit-testing.
 */
inline float easeOutCubic(float x) {
    const float inv = 1.0f - x;
    return 1.0f - inv * inv * inv;
}

inline void drawPresetMenu(const GLBitmapFont& font, HudFrame& frame,
                           PresetMenuState& state,
                           const std::string* profileNames,
                           const std::string* profileDescriptions,
                           int  numProfiles,
                           int  currentIdx,
                           int  sw, int sh) {
    if (state.t <= 0.001f) return;

    const float eased = easeOutCubic(std::min(std::max(state.t, 0.0f), 1.0f));

    // Side-panel layout (anchored to right edge, slides in from off-screen).
    const float panelW   = 300.0f;
    const float headerH  = 38.0f;
    const float itemH    = 28.0f;
    const float descH    = 56.0f;   // selected-item description block
    const float footerH  = 26.0f;
    const float padX     = 12.0f;
    const float margin   = 16.0f;

    const int   visible  = std::min(numProfiles, PresetMenuState::visibleCount);
    const float panelH   = headerH + visible * itemH + descH + footerH;

    // Slide in from the right: at t=0 panel sits fully off-screen, at t=1 it
    // rests at `margin` from the right edge.
    const float restX    = float(sw) - panelW - margin;
    const float panelX   = std::round(restX + (1.0f - eased) * (panelW + margin + 8.0f));
    const float panelY   = std::round((float(sh) - panelH) * 0.5f);

    // Drop shadow
    {
        const float a = 0.45f * eased;
        font.drawRect(panelX + 4.0f, panelY + 6.0f, panelW, panelH,
                      0.0f, 0.0f, 0.0f, a, sw, sh);
    }

    // Panel body + header bar
    Color bg     = col::panelBg;     bg.a     *= eased;
    Color border = col::panelBorder; border.a *= eased;
    Color title  = col::titleBar;    title.a  *= eased;
    drawFilled(font, panelX, panelY, panelW, panelH, bg, sw, sh);
    drawFilled(font, panelX, panelY, panelW, headerH, title, sw, sh);
    drawBorder(font, panelX, panelY, panelW, panelH, border, sw, sh, 1.0f);
    drawFilled(font, panelX, panelY + headerH, panelW, 1.0f,
               { col::divider.r, col::divider.g, col::divider.b, col::divider.a * eased },
               sw, sh);

    // Header text, left-aligned in the side panel
    {
        Color h = col::textAccent; h.a *= eased;
        drawString(font, "Available Presets",
                   panelX + padX,
                   panelY + (headerH - float(font.lineHeight())) * 0.5f - 1.0f,
                   h, sw, sh);
    }

    // Clamp scroll into range
    const int maxScroll = std::max(0, numProfiles - visible);
    if (state.scroll < 0) state.scroll = 0;
    if (state.scroll > maxScroll) state.scroll = maxScroll;
    if (state.selected < state.scroll) state.scroll = state.selected;
    if (state.selected >= state.scroll + visible)
        state.scroll = state.selected - visible + 1;

    // Item rows (single-line: index + name)
    const float itemX = panelX + padX * 0.5f;
    const float itemW = panelW - padX;
    for (int row = 0; row < visible; ++row) {
        const int idx = state.scroll + row;
        if (idx >= numProfiles) break;
        const float iy = panelY + headerH + row * itemH + 2.0f;
        const float ih = itemH - 4.0f;

        const bool isSel = (idx == state.selected);
        const bool isCur = (idx == currentIdx);

        if (isSel) {
            Color sel = col::btnBgActive; sel.a *= eased;
            drawFilled(font, itemX, iy, itemW, ih, sel, sw, sh);
            Color selB = col::btnBorder; selB.a *= eased;
            drawBorder(font, itemX, iy, itemW, ih, selB, sw, sh, 1.0f);
        }

        // Index label on the left ("01 ", dim)
        char idxBuf[8];
        std::snprintf(idxBuf, sizeof(idxBuf), "%02d", idx + 1);
        {
            Color ic = col::textDim; ic.a *= eased;
            drawString(font, idxBuf,
                       itemX + 8.0f,
                       iy + (ih - float(font.lineHeight())) * 0.5f,
                       ic, sw, sh);
        }

        // Name
        {
            Color tc = col::textPrimary; tc.a *= eased;
            drawString(font, profileNames[idx],
                       itemX + 36.0f,
                       iy + (ih - float(font.lineHeight())) * 0.5f,
                       tc, sw, sh);
        }

        // Active-profile dot on the right
        if (isCur) {
            const float dotR = 3.5f;
            Color dot = col::valueOn; dot.a *= eased;
            drawFilled(font, itemX + itemW - 14.0f, iy + ih * 0.5f - dotR,
                       dotR * 2.0f, dotR * 2.0f, dot, sw, sh);
        }

        frame.presetItems.push_back({itemX, iy, itemW, ih, idx});
    }

    // Description block for selected item (word-wrapped naive: just one line
    //, descriptions in profiles are already concise.)
    {
        const float dy = panelY + headerH + visible * itemH;
        drawFilled(font, panelX, dy, panelW, 1.0f,
                   { col::divider.r, col::divider.g, col::divider.b, col::divider.a * eased },
                   sw, sh);
        if (state.selected >= 0 && state.selected < numProfiles) {
            const std::string& nm = profileNames[state.selected];
            const std::string& ds = profileDescriptions[state.selected];
            Color nc = col::textAccent; nc.a *= eased;
            Color dc = col::textDim;    dc.a *= eased;
            drawString(font, nm,
                       panelX + padX, dy + 8.0f, nc, sw, sh);
            // Description (truncate visually if too wide; HUD has no wrap).
            const float maxW   = panelW - padX * 2.0f;
            std::string trimmed = ds;
            while (!trimmed.empty() && font.textWidth(trimmed) > maxW) {
                trimmed.pop_back();
            }
            if (trimmed.size() < ds.size() && trimmed.size() > 1) {
                trimmed.replace(trimmed.size() - 1, 1, ".");
                if (trimmed.size() > 2) trimmed[trimmed.size() - 2] = '.';
                if (trimmed.size() > 3) trimmed[trimmed.size() - 3] = '.';
            }
            drawString(font, trimmed,
                       panelX + padX,
                       dy + 8.0f + float(font.lineHeight()) + 4.0f,
                       dc, sw, sh);
        }
    }

    // Scroll arrows (right-edge mini buttons)
    const float arrowW = 28.0f;
    const float arrowH = 20.0f;
    const float arrowX = panelX + panelW - arrowW - 6.0f;
    if (state.scroll > 0) {
        const float ay = panelY + headerH + 4.0f;
        Color ab = col::btnBg;     ab.a *= eased * 0.7f;
        Color ac = col::textAccent; ac.a *= eased;
        drawFilled(font, arrowX, ay, arrowW, arrowH, ab, sw, sh);
        drawString(font, " ^",
                   arrowX + (arrowW - font.textWidth(" ^")) * 0.5f,
                   ay + (arrowH - float(font.lineHeight())) * 0.5f,
                   ac, sw, sh);
        frame.presetScrollUp = {arrowX, ay, arrowW, arrowH, -2};
        frame.presetScrollUpVisible = true;
    }
    if (state.scroll < maxScroll) {
        const float ay = panelY + headerH + visible * itemH - arrowH - 4.0f;
        Color ab = col::btnBg;     ab.a *= eased * 0.7f;
        Color ac = col::textAccent; ac.a *= eased;
        drawFilled(font, arrowX, ay, arrowW, arrowH, ab, sw, sh);
        drawString(font, " v",
                   arrowX + (arrowW - font.textWidth(" v")) * 0.5f,
                   ay + (arrowH - float(font.lineHeight())) * 0.5f,
                   ac, sw, sh);
        frame.presetScrollDown = {arrowX, ay, arrowW, arrowH, -3};
        frame.presetScrollDownVisible = true;
    }

    // Footer hint (compact)
    {
        Color hc = col::textDim; hc.a *= eased;
        const std::string hint = "Up/Down  Enter  Esc";
        const float tw = font.textWidth(hint);
        drawString(font, hint,
                   panelX + (panelW - tw) * 0.5f,
                   panelY + panelH - footerH + (footerH - float(font.lineHeight())) * 0.5f - 1.0f,
                   hc, sw, sh);
    }
}

} // namespace hud
