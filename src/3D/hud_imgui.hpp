// ============================================================
// hud_imgui.hpp — ImGui-based HUD overlays for the standalone
// 3D black hole viewer.
//
//  Phase 1: preset selector menu.
//  Phase 2: status, debug, controls hint, overlays toggles.
//
// In-scene body labels still use the GLBitmapFont path because
// they require world→screen projection from the camera state.
// ============================================================
#pragma once

#include <cstdio>
#include "imgui.h"
#include "bh3d_core.hpp"

namespace hud_im {

// ────────────────────────────────────────────────────────────
// Apply the Aetherion HUD palette/typography to ImGui.
// Call once after ImGui::SFML::Init().
// ────────────────────────────────────────────────────────────
inline void applyAetherionStyle() {
    ImGuiStyle& s = ImGui::GetStyle();

    // Geometry
    s.WindowRounding   = 6.0f;
    s.ChildRounding    = 6.0f;
    s.FrameRounding    = 4.0f;
    s.GrabRounding     = 4.0f;
    s.PopupRounding    = 6.0f;
    s.ScrollbarRounding= 6.0f;
    s.WindowBorderSize = 1.0f;
    s.FrameBorderSize  = 0.0f;
    s.WindowPadding    = ImVec2(14, 12);
    s.FramePadding     = ImVec2(8, 6);
    s.ItemSpacing      = ImVec2(8, 6);
    s.ItemInnerSpacing = ImVec2(6, 4);

    // Colors — match hud::col palette in hud_panel.hpp
    auto rgba = [](float r, float g, float b, float a) { return ImVec4(r, g, b, a); };
    ImVec4* c = s.Colors;
    c[ImGuiCol_WindowBg]            = rgba(0.05f, 0.06f, 0.09f, 0.92f);
    c[ImGuiCol_PopupBg]             = rgba(0.05f, 0.06f, 0.09f, 0.96f);
    c[ImGuiCol_ChildBg]             = rgba(0.00f, 0.00f, 0.00f, 0.00f);
    c[ImGuiCol_Border]              = rgba(0.25f, 0.32f, 0.45f, 0.85f);
    c[ImGuiCol_BorderShadow]        = rgba(0.00f, 0.00f, 0.00f, 0.00f);
    c[ImGuiCol_TitleBg]             = rgba(0.10f, 0.13f, 0.20f, 1.00f);
    c[ImGuiCol_TitleBgActive]       = rgba(0.10f, 0.13f, 0.20f, 1.00f);
    c[ImGuiCol_TitleBgCollapsed]    = rgba(0.10f, 0.13f, 0.20f, 0.85f);
    c[ImGuiCol_Separator]           = rgba(0.25f, 0.32f, 0.45f, 0.55f);
    c[ImGuiCol_Text]                = rgba(0.92f, 0.95f, 1.00f, 1.00f);
    c[ImGuiCol_TextDisabled]        = rgba(0.55f, 0.62f, 0.73f, 1.00f);

    // Buttons / selectables
    c[ImGuiCol_Button]              = rgba(0.14f, 0.18f, 0.27f, 0.85f);
    c[ImGuiCol_ButtonHovered]       = rgba(0.20f, 0.27f, 0.40f, 0.95f);
    c[ImGuiCol_ButtonActive]        = rgba(0.32f, 0.55f, 0.85f, 0.95f);
    c[ImGuiCol_Header]              = rgba(0.32f, 0.55f, 0.85f, 0.55f);
    c[ImGuiCol_HeaderHovered]       = rgba(0.32f, 0.55f, 0.85f, 0.75f);
    c[ImGuiCol_HeaderActive]        = rgba(0.32f, 0.55f, 0.85f, 0.95f);
    c[ImGuiCol_FrameBg]             = rgba(0.10f, 0.13f, 0.20f, 0.85f);
    c[ImGuiCol_FrameBgHovered]      = rgba(0.18f, 0.24f, 0.36f, 0.95f);
    c[ImGuiCol_FrameBgActive]       = rgba(0.32f, 0.55f, 0.85f, 0.85f);

    // Scrollbar
    c[ImGuiCol_ScrollbarBg]         = rgba(0.05f, 0.06f, 0.09f, 0.50f);
    c[ImGuiCol_ScrollbarGrab]       = rgba(0.25f, 0.32f, 0.45f, 0.85f);
    c[ImGuiCol_ScrollbarGrabHovered]= rgba(0.32f, 0.55f, 0.85f, 0.85f);
    c[ImGuiCol_ScrollbarGrabActive] = rgba(0.45f, 0.70f, 1.00f, 0.95f);
}

// ────────────────────────────────────────────────────────────
// Light variant — off-white panels, dark text, blue accents.
// Call in place of applyAetherionStyle() when light mode is on.
// ────────────────────────────────────────────────────────────
inline void applyAetherionStyleLight() {
    ImGuiStyle& s = ImGui::GetStyle();

    // Geometry (identical to dark — keeps layout consistent)
    s.WindowRounding   = 6.0f;
    s.ChildRounding    = 6.0f;
    s.FrameRounding    = 4.0f;
    s.GrabRounding     = 4.0f;
    s.PopupRounding    = 6.0f;
    s.ScrollbarRounding= 6.0f;
    s.WindowBorderSize = 1.0f;
    s.FrameBorderSize  = 0.0f;
    s.WindowPadding    = ImVec2(14, 12);
    s.FramePadding     = ImVec2(8, 6);
    s.ItemSpacing      = ImVec2(8, 6);
    s.ItemInnerSpacing = ImVec2(6, 4);

    auto rgba = [](float r, float g, float b, float a) { return ImVec4(r, g, b, a); };
    ImVec4* c = s.Colors;
    c[ImGuiCol_WindowBg]            = rgba(0.96f, 0.96f, 0.98f, 0.94f);
    c[ImGuiCol_PopupBg]             = rgba(0.98f, 0.98f, 1.00f, 0.97f);
    c[ImGuiCol_ChildBg]             = rgba(0.00f, 0.00f, 0.00f, 0.00f);
    c[ImGuiCol_Border]              = rgba(0.70f, 0.72f, 0.80f, 0.80f);
    c[ImGuiCol_BorderShadow]        = rgba(0.00f, 0.00f, 0.00f, 0.00f);
    c[ImGuiCol_TitleBg]             = rgba(0.88f, 0.90f, 0.95f, 1.00f);
    c[ImGuiCol_TitleBgActive]       = rgba(0.82f, 0.86f, 0.94f, 1.00f);
    c[ImGuiCol_TitleBgCollapsed]    = rgba(0.88f, 0.90f, 0.95f, 0.85f);
    c[ImGuiCol_Separator]           = rgba(0.70f, 0.72f, 0.80f, 0.55f);
    c[ImGuiCol_Text]                = rgba(0.08f, 0.09f, 0.14f, 1.00f);
    c[ImGuiCol_TextDisabled]        = rgba(0.50f, 0.52f, 0.60f, 1.00f);

    c[ImGuiCol_Button]              = rgba(0.84f, 0.86f, 0.92f, 0.90f);
    c[ImGuiCol_ButtonHovered]       = rgba(0.72f, 0.78f, 0.92f, 0.95f);
    c[ImGuiCol_ButtonActive]        = rgba(0.24f, 0.48f, 0.82f, 0.95f);
    c[ImGuiCol_Header]              = rgba(0.24f, 0.48f, 0.82f, 0.40f);
    c[ImGuiCol_HeaderHovered]       = rgba(0.24f, 0.48f, 0.82f, 0.60f);
    c[ImGuiCol_HeaderActive]        = rgba(0.24f, 0.48f, 0.82f, 0.85f);
    c[ImGuiCol_FrameBg]             = rgba(0.90f, 0.91f, 0.95f, 0.90f);
    c[ImGuiCol_FrameBgHovered]      = rgba(0.80f, 0.84f, 0.95f, 0.95f);
    c[ImGuiCol_FrameBgActive]       = rgba(0.24f, 0.48f, 0.82f, 0.70f);

    c[ImGuiCol_ScrollbarBg]         = rgba(0.92f, 0.93f, 0.96f, 0.50f);
    c[ImGuiCol_ScrollbarGrab]       = rgba(0.65f, 0.68f, 0.78f, 0.80f);
    c[ImGuiCol_ScrollbarGrabHovered]= rgba(0.40f, 0.55f, 0.80f, 0.85f);
    c[ImGuiCol_ScrollbarGrabActive] = rgba(0.24f, 0.48f, 0.82f, 0.95f);
}

// ────────────────────────────────────────────────────────────
// Draw the preset selector. Drives the same `bh3d::State`
// machinery (presetMenu.{open,selected,scroll,t}) and calls
// `bh3d::switchToProfile` on confirmation.
//
// Should be called between ImGui::SFML::Update and
// ImGui::SFML::Render. Returns true if the menu rendered this
// frame (caller may use this to gate other input).
// ────────────────────────────────────────────────────────────
inline bool drawPresetMenu(bh3d::State& s) {
    auto& m = s.presetMenu;
    if (m.t <= 0.001f) return false;

    // Smoothly fade alpha with the same easing the legacy menu used.
    const float t     = std::min(std::max(m.t, 0.0f), 1.0f);
    const float inv   = 1.0f - t;
    const float eased = 1.0f - inv * inv * inv;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 320.0f;
    // Slide in from the right — at t=0 the panel sits off-screen.
    const float restX  = vp->WorkPos.x + vp->WorkSize.x - panelW - 16.0f;
    const float offX   = (1.0f - eased) * (panelW + 24.0f);
    const ImVec2 pos   = ImVec2(restX + offX, vp->WorkPos.y + vp->WorkSize.y * 0.5f);
    const ImVec2 size  = ImVec2(panelW, std::min(vp->WorkSize.y - 32.0f, 540.0f));

    ImGui::SetNextWindowPos(pos, ImGuiCond_Always, ImVec2(0.0f, 0.5f));
    ImGui::SetNextWindowSize(size, ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.92f * eased);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize        |
        ImGuiWindowFlags_NoMove          |
        ImGuiWindowFlags_NoCollapse      |
        ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoTitleBar      |
        ImGuiWindowFlags_NoFocusOnAppearing;

    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, eased);
    if (ImGui::Begin("##AetherionPresets", nullptr, flags)) {
        // Custom header bar
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.85f, 0.92f, 1.00f, 1.00f));
        ImGui::TextUnformatted("Available Presets");
        ImGui::PopStyleColor();
        ImGui::Separator();
        ImGui::Spacing();

        // Scrollable list
        const float footerH = ImGui::GetFrameHeightWithSpacing() + 80.0f;
        ImGui::BeginChild("##presetList",
                          ImVec2(0, -footerH),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_None);

        for (int i = 0; i < profiles::NUM_PROFILES; ++i) {
            const auto& p = s.profilesArr[i];
            const bool  isCur = (i == s.profileIdx);

            ImGui::PushID(i);
            char idxBuf[8];
            std::snprintf(idxBuf, sizeof(idxBuf), "%02d", i + 1);

            // Row layout: index dim, name accent if selected, dot if current
            char label[256];
            std::snprintf(label, sizeof(label), "%s  %s%s",
                          idxBuf, p.name.c_str(), isCur ? "  *" : "");

            const bool selected = (i == m.selected);
            if (ImGui::Selectable(label, selected,
                                  ImGuiSelectableFlags_AllowDoubleClick,
                                  ImVec2(0, 0))) {
                m.selected = i;
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    bh3d::switchToProfile(s, i);
                    m.open = false;
                }
            }
            // Auto-scroll to keep selection visible when navigating with keys.
            if (selected && ImGui::IsKeyPressed(ImGuiKey_UpArrow,  false)) {}
            if (selected && ImGui::IsKeyPressed(ImGuiKey_DownArrow,false)) {}

            ImGui::PopID();
        }
        ImGui::EndChild();

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Description block for the highlighted item
        if (m.selected >= 0 && m.selected < profiles::NUM_PROFILES) {
            const auto& p = s.profilesArr[m.selected];
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.85f, 0.92f, 1.00f, 1.00f));
            ImGui::TextUnformatted(p.name.c_str());
            ImGui::PopStyleColor();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("%s", p.description.c_str());
            ImGui::PopTextWrapPos();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("Up/Down  Enter  Esc");
    }
    ImGui::End();
    ImGui::PopStyleVar();

    // Keyboard navigation (only when actually visible and open).
    if (m.open) {
        if (ImGui::IsKeyPressed(ImGuiKey_UpArrow, true)) {
            m.selected = std::max(0, m.selected - 1);
        }
        if (ImGui::IsKeyPressed(ImGuiKey_DownArrow, true)) {
            m.selected = std::min(profiles::NUM_PROFILES - 1, m.selected + 1);
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Enter, false) ||
            ImGui::IsKeyPressed(ImGuiKey_KeypadEnter, false)) {
            bh3d::switchToProfile(s, m.selected);
            m.open = false;
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Escape, false)) {
            m.open = false;
        }
    }

    return true;
}

// ────────────────────────────────────────────────────────────
// Helpers for the data panels below.
// ────────────────────────────────────────────────────────────
namespace detail {
    inline ImVec4 accent()     { return ImVec4(0.32f, 0.65f, 1.00f, 1.00f); }
    inline ImVec4 valueOn()    { return ImVec4(0.40f, 0.90f, 0.55f, 1.00f); }
    inline ImVec4 valueOff()   { return ImVec4(0.85f, 0.42f, 0.45f, 1.00f); }
    inline ImVec4 textPrim()   { return ImVec4(0.92f, 0.95f, 1.00f, 1.00f); }
    inline ImVec4 textDim()    { return ImVec4(0.62f, 0.68f, 0.78f, 1.00f); }

    inline void kvRow(const char* key, const char* val, ImVec4 valColor) {
        const float rowW = ImGui::GetContentRegionAvail().x;
        ImGui::PushStyleColor(ImGuiCol_Text, textDim());
        ImGui::TextUnformatted(key);
        ImGui::PopStyleColor();
        ImGui::SameLine();
        const float vw = ImGui::CalcTextSize(val).x;
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + std::max(0.0f, rowW - vw - ImGui::CalcTextSize(key).x - 4.0f));
        ImGui::PushStyleColor(ImGuiCol_Text, valColor);
        ImGui::TextUnformatted(val);
        ImGui::PopStyleColor();
    }
    inline void onOffRow(const char* key, bool on) {
        kvRow(key, on ? "ON" : "OFF", on ? valueOn() : valueOff());
    }
    inline void sectionHeader(const char* label) {
        ImGui::PushStyleColor(ImGuiCol_Text, accent());
        ImGui::TextUnformatted(label);
        ImGui::PopStyleColor();
    }
} // namespace detail

// ────────────────────────────────────────────────────────────
// Status panel (top-left).
// ────────────────────────────────────────────────────────────
inline void drawStatusPanel(const bh3d::State& s) {
    using namespace detail;
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + 12.0f, vp->WorkPos.y + 12.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(280.0f, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.85f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize;

    if (ImGui::Begin("Status", nullptr, flags)) {
        const auto& snap = s.snap;
        const float rollDeg = snap.roll * (180.0f / 3.14159265f);

        sectionHeader("Profile");
        kvRow("Profile", snap.profileName.c_str(), textPrim());
        ImGui::Separator();

        sectionHeader("Camera");
        onOffRow("Freelook", snap.freelook);
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%d deg", int(rollDeg + 0.5f));
        kvRow("Tilt", buf, textPrim());
        ImGui::Separator();

        sectionHeader("Toggles");
        onOffRow("Jets",        snap.jetsEnabled);
        onOffRow("BLR",         snap.blrEnabled);
        onOffRow("Orb Bodies",  snap.orbBodyEnabled);
        onOffRow("Doppler",     snap.dopplerEnabled);
        onOffRow("Blueshift",   snap.blueshiftEnabled);
        onOffRow("Host Galaxy", snap.hostGalaxyEnabled);
        onOffRow("LAB",         snap.labEnabled);
        onOffRow("CGM",         snap.cgmEnabled);
        if (snap.animSpeed >= 1.0f) std::snprintf(buf, sizeof(buf), "%dx", int(snap.animSpeed + 0.5f));
        else                        std::snprintf(buf, sizeof(buf), "%.2fx", snap.animSpeed);
        kvRow("Anim Speed", buf, textPrim());
        ImGui::Separator();

        sectionHeader("Performance");
        kvRow("Render", snap.cinematicMode ? "CINEMATIC" : "FAST",
              snap.cinematicMode ? valueOn() : textPrim());
        std::snprintf(buf, sizeof(buf), "%d", int(snap.fps + 0.5f));
        kvRow("FPS", buf, textPrim());
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Debug panel (top-right).
// ────────────────────────────────────────────────────────────
inline void drawDebugPanel(const bh3d::State& s) {
    using namespace detail;
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 280.0f;
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - panelW - 12.0f,
                                   vp->WorkPos.y + 12.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(panelW, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.85f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize;

    if (ImGui::Begin("Debug", nullptr, flags)) {
        const auto& snap = s.snap;
        const float rollDeg = snap.roll * (180.0f / 3.14159265f);
        char buf[64];

        sectionHeader("Shader");
        kvRow("Photoreal", s.havePhotoreal ? "compiled" : "N/A",
              s.havePhotoreal ? valueOn() : valueOff());
        kvRow("Simple", s.haveSimple ? "compiled" : "N/A",
              s.haveSimple ? valueOn() : valueOff());
        ImGui::Separator();

        sectionHeader("Camera");
        std::snprintf(buf, sizeof(buf), "(%.1f, %.1f, %.1f)", snap.cameraPos.x, snap.cameraPos.y, snap.cameraPos.z);
        kvRow("Pos", buf, textPrim());
        std::snprintf(buf, sizeof(buf), "(%.2f, %.2f, %.2f)", snap.cameraDir.x, snap.cameraDir.y, snap.cameraDir.z);
        kvRow("Dir", buf, textPrim());
        std::snprintf(buf, sizeof(buf), "%.1f deg", rollDeg);
        kvRow("Tilt", buf, textPrim());
        ImGui::Separator();

        sectionHeader("Render");
        std::snprintf(buf, sizeof(buf), "%d", snap.maxSteps);
        kvRow("Max Steps", buf, textPrim());
        kvRow("Cinematic", snap.cinematicMode ? "ON" : "OFF",
              snap.cinematicMode ? valueOn() : valueOff());
        ImGui::Separator();

        sectionHeader("Animation");
        if (snap.animSpeed >= 1.0f) std::snprintf(buf, sizeof(buf), "%.0fx", snap.animSpeed);
        else                        std::snprintf(buf, sizeof(buf), "%.2fx", snap.animSpeed);
        kvRow("Speed", buf, textPrim());
        std::snprintf(buf, sizeof(buf), "%.1fs", snap.animTime);
        kvRow("Anim Time", buf, textPrim());
        std::snprintf(buf, sizeof(buf), "%.1fs", snap.totalTime);
        kvRow("Total Time", buf, textPrim());
        ImGui::Separator();

        sectionHeader("Window");
        std::snprintf(buf, sizeof(buf), "%dx%d  @ %.0f FPS", snap.windowW, snap.windowH, snap.fps);
        kvRow("Size", buf, textPrim());
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Controls hint strip (bottom-center).
// ────────────────────────────────────────────────────────────
inline void drawControlsHint() {
    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const char* txt =
        "[WASD] Move  [Q/E] Tilt  [Space/LCTRL] Up/Down  [LSHIFT] Fast  "
        "[RMB] Look  [N] Next Profile  [H] Hide HUD  [B] Debug  [M] Overlays";

    ImGuiStyle& st = ImGui::GetStyle();
    const ImVec2 ts = ImGui::CalcTextSize(txt);
    const float winW = ts.x + st.WindowPadding.x * 2.0f;
    const float winH = ts.y + st.WindowPadding.y * 2.0f;

    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + (vp->WorkSize.x - winW) * 0.5f,
                                   vp->WorkPos.y + vp->WorkSize.y - winH - 10.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(winW, winH), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.80f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar;

    if (ImGui::Begin("##controls", nullptr, flags)) {
        ImGui::TextUnformatted(txt);
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Overlays panel (bottom-left). Toggles are interactive.
// `spacetimeViz`, `gravRedshift`, `kerrShadowGuide` are stubs
// (the legacy panel marked them disabled); we render them
// disabled here too for parity.
// ────────────────────────────────────────────────────────────
inline void drawOverlaysPanel(bh3d::State& s) {
    if (!s.overlays.panelOpen) return;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 260.0f;
    const float estH   = 180.0f;
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + 12.0f,
                                   vp->WorkPos.y + vp->WorkSize.y - estH - 60.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(panelW, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.88f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing |
        ImGuiWindowFlags_AlwaysAutoResize;

    if (ImGui::Begin("Overlays", nullptr, flags)) {
        const float btnW = ImGui::GetContentRegionAvail().x;

        // Active: Label View
        if (ImGui::Selectable(s.overlays.labelView ? "Label View   [ON]" : "Label View   [OFF]",
                              s.overlays.labelView, 0, ImVec2(btnW, 24))) {
            s.overlays.labelView = !s.overlays.labelView;
        }

        // Stubbed togglables (kept disabled to match legacy parity)
        ImGui::BeginDisabled();
        bool dummy;
        dummy = s.overlays.spacetimeViz;
        ImGui::Selectable("Spacetime Visualization (TBD)", dummy, 0, ImVec2(btnW, 24));
        dummy = s.overlays.gravRedshift;
        ImGui::Selectable("Gravitational Redshift (TBD)", dummy, 0, ImVec2(btnW, 24));
        dummy = s.overlays.kerrShadowGuide;
        ImGui::Selectable("Kerr Shadow Guide (TBD)", dummy, 0, ImVec2(btnW, 24));
        ImGui::EndDisabled();

        ImGui::Spacing();
        ImGui::PushStyleColor(ImGuiCol_Text, detail::textDim());
        ImGui::TextUnformatted("Click to toggle.");
        ImGui::PopStyleColor();
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Convenience: draw the full ImGui HUD for the standalone host.
// Call between ImGui::SFML::Update and ImGui::SFML::Render.
// ────────────────────────────────────────────────────────────
inline void drawAll(bh3d::State& s) {
    if (s.showHUD) {
        drawStatusPanel(s);
        drawControlsHint();
        drawOverlaysPanel(s);
    }
    if (s.showDebugHUD) {
        drawDebugPanel(s);
    }
    drawPresetMenu(s);

    // ---- Tidal disruption event particle overlay ----
    if (s.snap.tdeActive) {
        auto* dl = ImGui::GetForegroundDrawList();
        const auto& snap = s.snap;
        int sw = snap.windowW, sh = snap.windowH;

        // Project disruption point for flash and label
        float fsx = 0.0f, fsy = 0.0f;
        bool flashVisible = hud::projectToScreen(
            snap.tdeEventPos, snap.cameraPos, snap.cameraDir, snap.cameraUp,
            snap.fov, sw, sh, fsx, fsy);

        if (flashVisible && snap.tdeFlashAlpha > 0.0f) {
            float a = snap.tdeFlashAlpha;
            float t = 1.0f - a;

            // Expanding shockwave ring
            float ringR = 18.0f + t * 80.0f;
            uint32_t ringCol = IM_COL32(255, 160, 60, (int)(a * 175));
            dl->AddCircle(ImVec2(fsx, fsy), ringR, ringCol, 48, 2.5f);

            // Inner bright core
            float coreR = a * 20.0f;
            uint32_t coreCol = IM_COL32(255, 235, 160, (int)(a * 210));
            dl->AddCircleFilled(ImVec2(fsx, fsy), coreR, coreCol, 32);

            // Label
            const char* lbl = "TIDAL DISRUPTION EVENT";
            uint32_t lblCol = IM_COL32(255, 140, 40, 220);
            dl->AddText(ImVec2(fsx - 105.0f, fsy + coreR + 8.0f), lblCol, lbl);
        }

        // Debris / gas-stream particles
        for (size_t i = 0; i < snap.tdeDebrisPos.size(); ++i) {
            float px = 0.0f, py = 0.0f;
            if (!hud::projectToScreen(snap.tdeDebrisPos[i],
                                       snap.cameraPos, snap.cameraDir, snap.cameraUp,
                                       snap.fov, sw, sh, px, py))
                continue;
            float lifeF = snap.tdeDebrisLifeF[i];
            float r = 3.5f * lifeF;
            if (r < 0.5f) continue;

            uint32_t col;
            if (snap.tdeDebrisIsFallback[i]) {
                col = IM_COL32(140, 195, 255, (int)(lifeF * 195));
            } else {
                int g = (int)(80 + lifeF * 120);
                col = IM_COL32(255, g, 30, (int)(lifeF * 215));
            }
            dl->AddCircleFilled(ImVec2(px, py), r, col, 8);
        }
    }
}

} // namespace hud_im
