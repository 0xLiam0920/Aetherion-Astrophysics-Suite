/*  pages/overview.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::overview {

    static void drawModeCard(const char* modeNum, const char* title,
                             const char* description,
                             ui::LaunchMode mode, ui::AppState& state) {
        ImVec2 avail = ImGui::GetContentRegionAvail();
        ImGui::PushStyleColor(ImGuiCol_FrameBg,  ImVec4(0.10f, 0.12f, 0.18f, 0.9f));
        ImGui::PushStyleColor(ImGuiCol_Border,   ImVec4(0.15f, 0.45f, 0.75f, 0.6f));
        ImGui::BeginChild(title, ImVec2((avail.x - 20) * 0.45f, 180), true);
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "%s", modeNum);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "%s", title);
        ImGui::Spacing();
        ImGui::PushTextWrapPos();
        ImGui::TextDisabled("%s", description);
        ImGui::PopTextWrapPos();
        ImGui::Spacing();
        if (ImGui::Button("Launch"))
            state.launchMode = mode;
        ImGui::EndChild();
        ImGui::PopStyleColor(2);
    }

    inline void draw(ui::AppState& state) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Overview");
        ImGui::TextDisabled("Start here, browse recent progress, or pick a mode.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "Select Mode");
        ImGui::Spacing();
        drawModeCard("MODE 01", "2D Research Mode",
            "Advanced calculations, data compilation, scientific analysis, BLR simulation",
            ui::LaunchMode::Sim2D, state);
        ImGui::SameLine();
        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + 20);
        drawModeCard("MODE 02", "3D Visualization Mode",
            "Real-time rendering, visual exploration, gravitational lensing preview",
            ui::LaunchMode::Sim3D, state);

        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "Recently used");
        ImGui::TextDisabled("Resume progress files from your last session.");
        ImGui::Spacing();

        const char* recentFiles[] = {
            "progress-2026-04-10.ahp",
            "workspace-sgrA-2d.ahp",
            "lens-study-3d.ahp"
        };
        for (int i = 0; i < 3; ++i) {
            ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0.10f, 0.15f, 0.22f, 0.9f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.18f, 0.25f, 0.35f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(0.12f, 0.20f, 0.30f, 1.0f));
            if (ImGui::Button(recentFiles[i], ImVec2(-1.0f, 0))) {
                state.statusText  = std::string("Restored ") + recentFiles[i];
                state.sessionState = 1;
            }
            ImGui::PopStyleColor(3);
        }
    }

} // namespace pages::overview
