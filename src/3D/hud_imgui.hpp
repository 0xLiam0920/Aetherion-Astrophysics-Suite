// ============================================================
// hud_imgui.hpp: ImGui-based HUD overlays for the standalone
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

    // Colors, match hud::col palette in hud_panel.hpp
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
// Light variant, off-white panels, dark text, blue accents.
// Call in place of applyAetherionStyle() when light mode is on.
// ────────────────────────────────────────────────────────────
inline void applyAetherionStyleLight() {
    ImGuiStyle& s = ImGui::GetStyle();

    // Geometry (identical to dark, keeps layout consistent)
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
    // Slide in from the right, at t=0 the panel sits off-screen.
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
// Draw the merger selector. Lists the curated secondary objects
// from presets.hpp and calls bh3d::startMerger3D on confirmation.
// Slides in from the left so it never overlaps the preset menu.
// ────────────────────────────────────────────────────────────
inline bool drawMergerMenu(bh3d::State& s) {
    auto& m = s.mergerMenu;
    if (m.t <= 0.001f) return false;

    const float t     = std::min(std::max(m.t, 0.0f), 1.0f);
    const float inv   = 1.0f - t;
    const float eased = 1.0f - inv * inv * inv;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 340.0f;
    // Slide in from the left.
    const float restX  = vp->WorkPos.x + 16.0f;
    const float offX   = (1.0f - eased) * (panelW + 24.0f);
    const ImVec2 pos   = ImVec2(restX - offX, vp->WorkPos.y + vp->WorkSize.y * 0.5f);
    const ImVec2 size  = ImVec2(panelW, std::min(vp->WorkSize.y - 32.0f, 560.0f));

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

    const int N = NUM_MERGER_SECONDARY_3D_PRESETS;
    if (m.selected < 0) m.selected = 0;
    if (m.selected >= N) m.selected = N - 1;

    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, eased);
    if (ImGui::Begin("##AetherionMerger", nullptr, flags)) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.00f, 0.86f, 0.55f, 1.00f));
        ImGui::TextUnformatted("Trigger Merger");
        ImGui::PopStyleColor();
        ImGui::TextDisabled("Send a secondary spiralling into the primary");
        ImGui::Separator();
        ImGui::Spacing();

        const float footerH = ImGui::GetFrameHeightWithSpacing() * 3.0f + 96.0f;
        ImGui::BeginChild("##mergerList",
                          ImVec2(0, -footerH),
                          ImGuiChildFlags_None,
                          ImGuiWindowFlags_None);

        for (int i = 0; i < N; ++i) {
            const auto& p = MERGER_SECONDARY_3D_PRESETS[i];
            ImGui::PushID(i);
            char label[256];
            std::snprintf(label, sizeof(label), "%02d  %s", i + 1, p.label);

            const bool selected = (i == m.selected);
            if (ImGui::Selectable(label, selected,
                                  ImGuiSelectableFlags_AllowDoubleClick,
                                  ImVec2(0, 0))) {
                m.selected = i;
                if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                    bh3d::startMerger3D(s, p.kind, p.massSolar, p.profile);
                    m.open = false;
                }
            }
            // Keep the keyboard-selected row visible in the scrollable list.
            if (selected && m.scrollToSelected) {
                ImGui::SetScrollHereY(0.5f);
            }
            ImGui::PopID();
        }
        m.scrollToSelected = false;
        ImGui::EndChild();

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        // Description for the highlighted secondary.
        if (m.selected >= 0 && m.selected < N) {
            const auto& p = MERGER_SECONDARY_3D_PRESETS[m.selected];
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.00f, 0.86f, 0.55f, 1.00f));
            ImGui::TextUnformatted(p.label);
            ImGui::PopStyleColor();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextDisabled("%s", p.blurb);
            ImGui::PopTextWrapPos();
        }

        ImGui::Spacing();
        // Cinematic time-scale selector.
        ImGui::TextDisabled("Pacing");
        auto scaleBtn = [&](const char* name, bh3d::State::MergerState3D::TimeScale ts) {
            const bool on = (s.merger.timeScale == ts);
            if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.32f, 0.55f, 0.85f, 0.95f));
            if (ImGui::Button(name)) s.merger.timeScale = ts;
            if (on) ImGui::PopStyleColor();
        };
        scaleBtn("Cinematic", bh3d::State::MergerState3D::TimeScale::Cinematic);
        ImGui::SameLine();
        scaleBtn("Default",   bh3d::State::MergerState3D::TimeScale::Default);
        ImGui::SameLine();
        scaleBtn("Realtime",  bh3d::State::MergerState3D::TimeScale::Realtime);

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("Up/Down  Enter=launch  Esc");
    }
    ImGui::End();
    ImGui::PopStyleVar();

    if (m.open) {
        if (ImGui::IsKeyPressed(ImGuiKey_UpArrow, true)) {
            m.selected = std::max(0, m.selected - 1);
            m.scrollToSelected = true;
        }
        if (ImGui::IsKeyPressed(ImGuiKey_DownArrow, true)) {
            m.selected = std::min(N - 1, m.selected + 1);
            m.scrollToSelected = true;
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Enter, false) ||
            ImGui::IsKeyPressed(ImGuiKey_KeypadEnter, false)) {
            const auto& p = MERGER_SECONDARY_3D_PRESETS[m.selected];
            bh3d::startMerger3D(s, p.kind, p.massSolar, p.profile);
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
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs;

    if (ImGui::Begin("Status", nullptr, flags)) {
        const auto& snap = s.snap;
        const float rollDeg = snap.roll * (180.0f / 3.14159265f);

        sectionHeader("Profile");
        kvRow("Profile", snap.profileName.c_str(), accent());
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
        onOffRow("RK4 Orbits",  s.physOverlay.orbitsEnabled);
        onOffRow("Photons",     s.physOverlay.photonsEnabled);
        onOffRow("Spacetime",   s.physOverlay.spacetimeEnabled);
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
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs;

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
// Physical scale panel (bottom-right, above the legend).
// Shows Rs in km, disk extent in AU or light-days, and a graphical
// scale bar that auto-switches units depending on BH mass.
// Also exposed via the Overlays panel as the true-scale toggle.
// ────────────────────────────────────────────────────────────
inline void drawPhysicalScalePanel(const bh3d::State& s) {
    const auto& snap = s.snap;
    if (snap.massSolar <= 0.0) return;

    // Physical constants
    const double Rs_km     = 2.953 * snap.massSolar;   // Rs_primary in km
    const double Rs_AU     = Rs_km / 1.496e8;           // Rs in AU
    const double Rs_ld     = Rs_km / 2.592e10;          // Rs in light-days

    // Choose display unit for scale bar based on BH mass
    // < 1e4 M☉ → km,  1e4–1e9 M☉ → AU,  > 1e9 M☉ → light-days
    const bool useKm = snap.massSolar < 1e4;
    const bool useAU = snap.massSolar >= 1e4 && snap.massSolar < 1e9;
    // (light-days for the SMBH end)

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 240.0f;
    ImGui::SetNextWindowPos(
        ImVec2(vp->WorkPos.x + vp->WorkSize.x - panelW - 12.0f,
               vp->WorkPos.y + vp->WorkSize.y - 180.0f),
        ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(panelW, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.82f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs;

    if (ImGui::Begin("##PhysScale", nullptr, flags)) {
        char buf[64];

        ImGui::TextColored(ImVec4(0.66f, 0.78f, 1.0f, 1.0f), "Physical Scale");
        ImGui::Separator();

        // Rs readout
        if (useKm) {
            if (Rs_km >= 1000.0)
                std::snprintf(buf, sizeof(buf), "%.0f km (%.2f Mm)", Rs_km, Rs_km / 1000.0);
            else
                std::snprintf(buf, sizeof(buf), "%.1f km", Rs_km);
        } else if (useAU) {
            std::snprintf(buf, sizeof(buf), "%.3g AU", Rs_AU);
        } else {
            std::snprintf(buf, sizeof(buf), "%.3g light-days", Rs_ld);
        }
        detail::kvRow("1 Rs  =", buf, detail::textPrim());

        // Disk outer extent
        if (snap.diskOuterRadius > 0.0f) {
            double extKm = snap.diskOuterRadius * Rs_km;
            if (useKm)
                std::snprintf(buf, sizeof(buf), "%.0f km", extKm);
            else if (useAU)
                std::snprintf(buf, sizeof(buf), "%.2f AU", extKm / 1.496e8);
            else
                std::snprintf(buf, sizeof(buf), "%.3g ld", extKm / 2.592e10);
            detail::kvRow("Disk outer  =", buf, detail::textPrim());
        }

        // Graphical scale bar: 1 Rs wide in screen-pixel analogy
        // We draw a coloured bar representing "1 Rs" to give a visual anchor.
        ImGui::Spacing();
        const float barW = ImGui::GetContentRegionAvail().x;
        ImVec2 p = ImGui::GetCursorScreenPos();
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float barH = 8.0f;
        // Gradient: dark teal → cyan (horizon colour)
        dl->AddRectFilledMultiColor(
            p, ImVec2(p.x + barW, p.y + barH),
            IM_COL32(20, 80, 120, 200),  IM_COL32(80, 200, 255, 220),
            IM_COL32(80, 200, 255, 220), IM_COL32(20, 80, 120, 200));
        ImGui::Dummy(ImVec2(barW, barH));
        ImGui::PushStyleColor(ImGuiCol_Text, detail::textDim());
        if (useKm)
            std::snprintf(buf, sizeof(buf), "\xe2\x86\x94 1 Rs = %.1f km", Rs_km);
        else if (useAU)
            std::snprintf(buf, sizeof(buf), "\xe2\x86\x94 1 Rs = %.3g AU", Rs_AU);
        else
            std::snprintf(buf, sizeof(buf), "\xe2\x86\x94 1 Rs = %.3g light-days", Rs_ld);
        ImGui::TextUnformatted(buf);
        ImGui::PopStyleColor();

        // True-scale mode indicator
        if (snap.trueScaleMode) {
            ImGui::Spacing();
            ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.2f, 1.0f),
                               "TRUE SCALE: bodies at physical size");
        }
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Gravitational-wave strain sparkline for mergers
// The purpose of this is to plot the scrolling h(t) ∝ (1/sep)·cos(2φ_orbit) so that
// frequency/amplitude sweep is visible as the bodies coalesce together. 
// ────────────────────────────────────────────────────────────
inline void drawMergerGWPlot(const bh3d::State& s) {
    const auto& snap = s.snap;
    if (!snap.mergerInspiral || snap.mergerGWWaveform.size() < 2) return;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    const float panelW = 300.0f, plotH = 70.0f;
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + (vp->WorkSize.x - panelW) * 0.5f,
                                   vp->WorkPos.y + 12.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(panelW, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.80f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoInputs;

    if (ImGui::Begin("##GWStrain", nullptr, flags)) {
        char ovl[64];
        std::snprintf(ovl, sizeof(ovl), "h(t)  sep %.1f Rs", snap.mergerSepRs);
        // Symmetric scale so the waveform stays centred as amplitude grows exponentially
        float amp = 0.0f;
        for (float v : snap.mergerGWWaveform) amp = std::max(amp, std::fabs(v));
        amp = std::max(amp, 1e-3f);
        ImGui::TextColored(ImVec4(0.66f, 0.78f, 1.0f, 1.0f), "Gravitational Wave Strain");
        ImGui::PlotLines("##gw", snap.mergerGWWaveform.data(),
                         (int)snap.mergerGWWaveform.size(), 0, ovl,
                         -amp, amp, ImVec2(panelW - 16.0f, plotH));
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
        ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoInputs;

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
    // Cap height so the panel scrolls instead of overflowing the viewport on
    // shorter windows. ~58% of viewport, hard-clamped to a sane range.
    const float maxH   = std::max(180.0f, std::min(vp->WorkSize.y * 0.58f, 420.0f));
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + 12.0f,
                                   vp->WorkPos.y + vp->WorkSize.y - maxH - 60.0f),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(panelW, maxH), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.88f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing;

    if (ImGui::Begin("Overlays", nullptr, flags)) {
        const float btnW = ImGui::GetContentRegionAvail().x;

        // Helper: a Selectable that paints itself green when active so the
        // audience can see at a glance which layers are live.
        auto toggleRow = [&](const char* labelOn, const char* labelOff,
                             bool active, ImVec2 sz) -> bool {
            const ImVec4 onCol(0.55f, 0.95f, 0.65f, 1.0f);
            if (active) ImGui::PushStyleColor(ImGuiCol_Text, onCol);
            bool clicked = ImGui::Selectable(active ? labelOn : labelOff,
                                             active, 0, sz);
            if (active) ImGui::PopStyleColor();
            return clicked;
        };

        ImGui::SeparatorText("Visuals");
        if (toggleRow("Label View   [ON]", "Label View   [OFF]",
                      s.overlays.labelView, ImVec2(btnW, 24))) {
            s.overlays.labelView = !s.overlays.labelView;
        }

        ImGui::SeparatorText("Physics");
        if (toggleRow("RK4 GR Orbits   [ON]", "RK4 GR Orbits   [OFF]",
                      s.physOverlay.orbitsEnabled, ImVec2(btnW, 24))) {
            s.physOverlay.orbitsEnabled = !s.physOverlay.orbitsEnabled;
            s.physOverlay.markDirty();
        }
        if (toggleRow("Photon Geodesics  [ON]", "Photon Geodesics  [OFF]",
                      s.physOverlay.photonsEnabled, ImVec2(btnW, 24))) {
            s.physOverlay.photonsEnabled = !s.physOverlay.photonsEnabled;
            s.physOverlay.markDirty();
        }
        if (toggleRow("Spacetime Curvature  [ON]", "Spacetime Curvature  [OFF]",
                      s.physOverlay.spacetimeEnabled, ImVec2(btnW, 24))) {
            s.physOverlay.spacetimeEnabled = !s.physOverlay.spacetimeEnabled;
            s.physOverlay.markDirty();
        }

        ImGui::SeparatorText("Scale");
        if (toggleRow("True Scale Bodies  [ON]", "True Scale Bodies  [OFF]",
                      s.trueScaleMode, ImVec2(btnW, 24))) {
            s.trueScaleMode = !s.trueScaleMode;
        }
        if (s.trueScaleMode) {
            ImGui::PushStyleColor(ImGuiCol_Text, detail::textDim());
            ImGui::TextWrapped("Bodies shown at physical Rs size (may be sub-pixel).");
            ImGui::PopStyleColor();
        }

        ImGui::SeparatorText("Disk Colour");
        if (toggleRow("Physical Disk Colour  [ON]", "Physical Disk Colour  [OFF]",
                      s.physicalDiskColor, ImVec2(btnW, 24))) {
            s.physicalDiskColor = !s.physicalDiskColor;
        }
        if (s.physicalDiskColor) {
            ImGui::PushStyleColor(ImGuiCol_Text, detail::textDim());
            ImGui::TextWrapped("T(r) \u221d r\u207b\u00be Novikov\u2013Thorne + Doppler/grav shift \u2192 CIE sRGB. "
                               "Changing spin or mass now visibly alters disc colour.");
            ImGui::PopStyleColor();
        }

        ImGui::SeparatorText("Experimental");
        // Stubbed togglables (kept disabled to match legacy parity)
        ImGui::BeginDisabled();
        bool dummy;
        dummy = s.overlays.spacetimeViz;
        ImGui::Selectable("Legacy Spacetime Stub (TBD)", dummy, 0, ImVec2(btnW, 24));
        dummy = s.overlays.gravRedshift;
        ImGui::Selectable("Gravitational Redshift (TBD)", dummy, 0, ImVec2(btnW, 24));
        dummy = s.overlays.kerrShadowGuide;
        ImGui::Selectable("Kerr Shadow Guide (TBD)", dummy, 0, ImVec2(btnW, 24));
        ImGui::EndDisabled();

        ImGui::Spacing();
        ImGui::PushStyleColor(ImGuiCol_Text, detail::textDim());
        ImGui::TextUnformatted("Click to toggle.  Green = active.");
        ImGui::PopStyleColor();
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Floating legend for the RK4/photon overlays. Anchored bottom-right.
// Auto-shows whenever either overlay is enabled.
// ────────────────────────────────────────────────────────────
inline void drawPhysOverlayLegend(bh3d::State& s) {
    if (!s.physOverlay.orbitsEnabled && !s.physOverlay.photonsEnabled &&
        !s.physOverlay.spacetimeEnabled) return;

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x - 16.0f,
                                   vp->WorkPos.y + vp->WorkSize.y - 16.0f),
                            ImGuiCond_Always, ImVec2(1.0f, 1.0f));
    ImGui::SetNextWindowBgAlpha(0.80f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoInputs;

    if (ImGui::Begin("##PhysLegend", nullptr, flags)) {
        auto swatch = [](ImVec4 c) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 p = ImGui::GetCursorScreenPos();
            float sz = ImGui::GetTextLineHeight();
            dl->AddRectFilled(p, ImVec2(p.x + sz, p.y + sz),
                              ImGui::ColorConvertFloat4ToU32(c));
            ImGui::Dummy(ImVec2(sz, sz));
            ImGui::SameLine();
        };
        if (s.physOverlay.photonsEnabled) {
            swatch(ImVec4(0.45f, 0.85f, 1.0f, 1.0f));
            ImGui::TextUnformatted("Photon escape (b > b_crit)");
            swatch(ImVec4(1.0f, 0.35f, 0.20f, 1.0f));
            ImGui::TextUnformatted("Photon capture (b < b_crit)");
        }
        if (s.physOverlay.orbitsEnabled) {
            swatch(ImVec4(0.75f, 0.85f, 1.0f, 1.0f));
            ImGui::TextUnformatted("RK4 GR orbit (rosette precession)");
        }
        if (s.physOverlay.spacetimeEnabled) {
            swatch(ImVec4(0.65f, 0.80f, 1.0f, 1.0f));
            ImGui::TextUnformatted("Spacetime curvature (Flamm)");
        }
        ImGui::Separator();
        ImGui::TextDisabled("K = orbits   L = photons   T = spacetime");
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Body inspector card.
// Anchored near the selected body's screen position (clamped to
// viewport edges).  Double-click the card to dismiss; Escape also
// works via onActionKey.  All data is already in the snapshot.
// ────────────────────────────────────────────────────────────
inline void drawBodyInspector(bh3d::State& s) {
    if (s.selectedBody == -1) return;
    const auto& snap = s.snap;
    const int   sw   = snap.windowW, sh = snap.windowH;

    // ── Resolve selected body data ─────────────────────────────────────────
    glm::vec3   worldPos;
    std::string bodyName, bodyType;
    float       visRadius  = 1.0f;
    float       semiMajorRs = 0.0f; // for Kepler period; 0 → skip

    const bool isMerger = (s.selectedBody == bh3d::State::BODY_MERGER_SECONDARY);
    if (isMerger) {
        if (!snap.mergerInspiral) { s.selectedBody = -1; return; }
        worldPos   = snap.mergerSecondaryPos;
        bodyName   = snap.mergerSecondaryLabel;
        bodyType   = "Merger secondary";
        visRadius  = s.merger.secRadius;
        // Semi-major: use current separation as proxy (no orbital config for merger)
        semiMajorRs = snap.mergerSepRs;
    } else {
        int idx = s.selectedBody;
        if (idx < 0 || idx >= (int)snap.orbBodyPositions.size()) {
            s.selectedBody = -1; return;
        }
        worldPos  = snap.orbBodyPositions[idx];
        bodyName  = (idx < (int)snap.orbBodyLabels.size()) ? snap.orbBodyLabels[idx] : "Body";
        visRadius = (idx < (int)snap.orbBodyRadii.size())  ? snap.orbBodyRadii[idx]  : 1.0f;
        int  btype = (idx < (int)snap.orbBodyTypes.size()) ? snap.orbBodyTypes[idx]  : 0;
        switch (btype) {
            case 0: bodyType = "Star";            break;
            case 1: bodyType = "Gas Cloud";       break;
            case 2: bodyType = "Stellar Cluster"; break;
            case 3: bodyType = "Dwarf Galaxy";    break;
            case 4: bodyType = "Neutron Star";    break;
            case 5: bodyType = "White Dwarf";     break;
            case 6: bodyType = "Companion Star";  break;
            default: bodyType = "Body";           break;
        }
        // Use semi-major axis from the OrbitalBody config if available
        if (idx < (int)s.orbBodies.size())
            semiMajorRs = s.orbBodies[idx].config().semiMajor;
    }

    // ── Derived physics ───────────────────────────────────────────────────
    const float Rs_km   = (float)(2.953 * std::max(1e-6, snap.massSolar));
    // Distance from BH in Rs (native sim unit)
    float rRs = glm::length(worldPos - snap.bhPosition);
    if (rRs < 1.001f) rRs = 1.001f; // clamp outside horizon

    // Physical distance from horizon: (r - 1) × Rs_km  [horizon = 1 Rs]
    float distKm  = (rRs - 1.0f) * Rs_km;
    // Adaptive unit formatter: km < 0.1 AU → AU < 0.1 ly → ly
    constexpr float km_per_AU  = 1.496e8f;
    constexpr float km_per_ly  = 9.461e12f;
    char distBuf[48];
    if (distKm < 0.1f * km_per_AU)
        std::snprintf(distBuf, sizeof(distBuf), "%.1f km", distKm);
    else if (distKm < 0.1f * km_per_ly)
        std::snprintf(distBuf, sizeof(distBuf), "%.3g AU", distKm / km_per_AU);
    else
        std::snprintf(distBuf, sizeof(distBuf), "%.3g ly", distKm / km_per_ly);

    // Gravitational time-dilation factor √(1 - 1/r)  (Schwarzschild, r in Rs)
    float tdil = std::sqrt(std::max(0.0f, 1.0f - 1.0f / rRs));

    // Gravitational redshift z = 1/√(1 - 1/r) - 1
    float gRedshift = (tdil > 1e-6f) ? (1.0f / tdil - 1.0f) : 1e6f;

    // Orbital velocity v/c (Newtonian approximation valid at r >> Rs)
    //   v = √(M / r)  in geometrised units (M = 0.5 Rs, r in Rs)
    float vOverC = std::sqrt(0.5f / rRs); // rough Keplerian

    // Keplerian period from snap: T = 2π √(a³ / M), M = Rs/2 (geom units)
    char periodBuf[48] = "—";
    if (semiMajorRs > 0.0f && snap.massSolar > 0.0) {
        // a in metres = semiMajorRs * Rs_km * 1000
        // M in metres = massSolar * 1476.25 (G M_sun / c²)
        double a_m = (double)semiMajorRs * (double)Rs_km * 1000.0;
        double M_m = snap.massSolar * 1476.25;  // geometric mass in metres
        double T_s = 2.0 * 3.14159265 * std::sqrt(a_m * a_m * a_m / M_m) / 3.0e8;
        if (T_s < 60.0)
            std::snprintf(periodBuf, sizeof(periodBuf), "%.1f s", T_s);
        else if (T_s < 3600.0)
            std::snprintf(periodBuf, sizeof(periodBuf), "%.2f min", T_s / 60.0);
        else if (T_s < 86400.0)
            std::snprintf(periodBuf, sizeof(periodBuf), "%.2f h", T_s / 3600.0);
        else if (T_s < 3.156e7)
            std::snprintf(periodBuf, sizeof(periodBuf), "%.2f days", T_s / 86400.0);
        else
            std::snprintf(periodBuf, sizeof(periodBuf), "%.3g yr", T_s / 3.156e7);
    }

    // ── Screen-space anchor (clamp to viewport) ───────────────────────────
    float bsx = 0.0f, bsy = 0.0f;
    const bool visible = hud::projectToScreen(worldPos,
        snap.cameraPos, snap.cameraDir, snap.cameraUp,
        snap.fov, sw, sh, bsx, bsy);
    // If the body went off-screen while followed, anchor to a corner
    if (!visible) { bsx = (float)sw - 280.0f; bsy = 200.0f; }

    const float cardW  = 260.0f;
    const float margin = 14.0f;
    float winX = bsx + 22.0f;
    float winY = bsy - 20.0f;
    // Clamp to screen
    winX = std::clamp(winX, margin, (float)sw - cardW - margin);
    winY = std::clamp(winY, margin, (float)sh - 300.0f);

    const ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + winX, vp->WorkPos.y + winY),
                            ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(cardW, 0.0f), ImGuiCond_Always);
    ImGui::SetNextWindowBgAlpha(0.90f);

    const ImGuiWindowFlags flags =
        ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings |
        ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav |
        ImGuiWindowFlags_AlwaysAutoResize;

    if (ImGui::Begin("##BodyInspector", nullptr, flags)) {
        // ── Header ──────────────────────────────────────────────────────────
        ImGui::TextColored(ImVec4(1.0f, 0.88f, 0.55f, 1.0f), "%s", bodyName.c_str());
        ImGui::SameLine();
        ImGui::TextDisabled("  %s", bodyType.c_str());
        ImGui::Separator();

        char buf[64];
        auto kvRow = [&](const char* key, const char* val) {
            ImGui::TextDisabled("%s", key);
            ImGui::SameLine(100.0f);
            ImGui::TextColored(detail::textPrim(), "%s", val);
        };

        // ── Position ──────────────────────────────────────────────────────────
        std::snprintf(buf, sizeof(buf), "%.2f Rs", rRs);
        kvRow("Dist (BH)", buf);
        kvRow("From horizon", distBuf);

        glm::vec3 delta = worldPos - snap.bhPosition;
        std::snprintf(buf, sizeof(buf), "(%.2f, %.2f, %.2f) Rs", delta.x, delta.y, delta.z);
        ImGui::TextDisabled("Vector ");
        ImGui::SameLine(100.0f);
        ImGui::TextColored(detail::textDim(), "%s", buf);

        ImGui::Separator();

        // ── Physics ────────────────────────────────────────────────────────────
        std::snprintf(buf, sizeof(buf), "%.4f", tdil);
        kvRow("Time dil √(1-1/r)", buf);

        std::snprintf(buf, sizeof(buf), "z = %.4f", gRedshift);
        kvRow("Grav. redshift", buf);

        std::snprintf(buf, sizeof(buf), "%.3f c", vOverC);
        kvRow("Orbital v ≈", buf);

        kvRow("Orbital period", periodBuf);

        // ── Dismiss hint ───────────────────────────────────────────────────────
        ImGui::Spacing();
        ImGui::TextDisabled("Double-click to close  |  Esc to deselect");

        // Double-click inside the card → dismiss
        if (ImGui::IsWindowHovered() && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
            bh3d::deselectBody(s);
        }
    }
    ImGui::End();
}

// ────────────────────────────────────────────────────────────
// Convenience: draw the full ImGui HUD for the standalone host.
// Call between ImGui::SFML::Update and ImGui::SFML::Render.
// ────────────────────────────────────────────────────────────
inline void drawAll(bh3d::State& s) {
    drawPhysOverlayLegend(s);
    if (s.showHUD) {
        drawStatusPanel(s);
        drawPhysicalScalePanel(s);
        drawControlsHint();
        drawOverlaysPanel(s);
    }
    if (s.showDebugHUD) {
        drawDebugPanel(s);
    }
    drawPresetMenu(s);
    drawMergerMenu(s);
    drawMergerGWPlot(s);
    drawBodyInspector(s);

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

    // ---- Visual black hole merger overlay ----
    if (s.snap.mergerActive) {
        auto* dl = ImGui::GetForegroundDrawList();
        const auto& snap = s.snap;
        const int sw = snap.windowW, sh = snap.windowH;

        // Death-spiral trail: project world points, draw fading segments.
        if (snap.mergerTrail.size() >= 2) {
            float px0 = 0.0f, py0 = 0.0f;
            bool  have0 = hud::projectToScreen(
                snap.mergerTrail[0], snap.cameraPos, snap.cameraDir, snap.cameraUp,
                snap.fov, sw, sh, px0, py0);
            const size_t n = snap.mergerTrail.size();
            for (size_t i = 1; i < n; ++i) {
                float px1 = 0.0f, py1 = 0.0f;
                bool have1 = hud::projectToScreen(
                    snap.mergerTrail[i], snap.cameraPos, snap.cameraDir, snap.cameraUp,
                    snap.fov, sw, sh, px1, py1);
                if (have0 && have1) {
                    float f = (float)i / (float)n;          // 0 old → 1 newest
                    int   a = (int)(f * 170.0f);
                    float w = 1.0f + f * 2.0f;
                    dl->AddLine(ImVec2(px0, py0), ImVec2(px1, py1),
                                IM_COL32(150, 200, 255, a), w);
                }
                px0 = px1; py0 = py1; have0 = have1;
            }
        }

        // Inspiralling secondary: marker + label + separation readout.
        if (snap.mergerInspiral) {
            float sx = 0.0f, sy = 0.0f;
            if (hud::projectToScreen(snap.mergerSecondaryPos,
                                     snap.cameraPos, snap.cameraDir, snap.cameraUp,
                                     snap.fov, sw, sh, sx, sy)) {
                dl->AddCircle(ImVec2(sx, sy), 14.0f, IM_COL32(255, 215, 140, 200), 24, 1.5f);
                char buf[96];
                std::snprintf(buf, sizeof(buf), "%s", snap.mergerSecondaryLabel.c_str());
                dl->AddText(ImVec2(sx + 18.0f, sy - 8.0f), IM_COL32(255, 225, 170, 230), buf);
                char sep[64];
                std::snprintf(sep, sizeof(sep), "sep %.1f Rs", snap.mergerSepRs);
                dl->AddText(ImVec2(sx + 18.0f, sy + 8.0f), IM_COL32(170, 200, 255, 200), sep);
            }
        }

        // Coalescence point: gravitational-wave shockwave ring + flash core.
        float ex = 0.0f, ey = 0.0f;
        const bool evtVisible = hud::projectToScreen(
            snap.mergerEventPos, snap.cameraPos, snap.cameraDir, snap.cameraUp,
            snap.fov, sw, sh, ex, ey);

        if (evtVisible && snap.mergerRingActive) {
            // Two concentric expanding rings read as a GW wavefront.
            float rt = std::min(std::max(snap.mergerRingT, 0.0f), 1.0f);
            float ease = 1.0f - (1.0f - rt) * (1.0f - rt);
            int   alpha = (int)((1.0f - rt) * 200.0f);
            float baseR = 24.0f + ease * (float)sw * 0.45f;
            dl->AddCircle(ImVec2(ex, ey), baseR,        IM_COL32(150, 200, 255, alpha), 64, 3.0f);
            dl->AddCircle(ImVec2(ex, ey), baseR * 0.7f, IM_COL32(200, 225, 255, alpha / 2), 64, 2.0f);
        }

        // Full-screen coalescence flash (fades from white).
        if (snap.mergerFlashAlpha > 0.0f) {
            int a = (int)(std::min(snap.mergerFlashAlpha, 1.0f) * 220.0f);
            dl->AddRectFilled(ImVec2(0, 0), ImVec2((float)sw, (float)sh),
                              IM_COL32(255, 250, 240, a));
            if (evtVisible) {
                float coreR = snap.mergerFlashAlpha * 40.0f;
                dl->AddCircleFilled(ImVec2(ex, ey), coreR,
                                    IM_COL32(255, 255, 255, (int)(snap.mergerFlashAlpha * 240)), 32);
            }
        }

        // Status label near the coalescence point.
        if (evtVisible) {
            if (snap.mergerFlashAlpha > 0.05f) {
                dl->AddText(ImVec2(ex - 34.0f, ey - 60.0f),
                            IM_COL32(255, 240, 180, 240), "COALESCENCE");
            } else if (snap.mergerRemnantAlpha > 0.0f) {
                int a = (int)(std::min(snap.mergerRemnantAlpha, 1.0f) * 230.0f);
                dl->AddText(ImVec2(ex - 28.0f, ey - 60.0f),
                            IM_COL32(180, 215, 255, a), "REMNANT");
            }
        }
    }
}

} // namespace hud_im
