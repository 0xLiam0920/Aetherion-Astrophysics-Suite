/*  pages/simulation.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::simulation {

    inline void draw(ui::AppState& state) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Simulation");
        ImGui::TextDisabled("Load a saved simulation or start a new 2D/3D run.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::InputText("Simulation file path", state.loadPath, sizeof(state.loadPath));
        if (ImGui::Button("Load simulation")) {
            state.statusText  = state.loadPath[0]
                ? std::string("Loaded ") + state.loadPath
                : "No file selected";
            if (state.loadPath[0]) state.sessionState = 1;
        }
        ImGui::SameLine();
        if (ImGui::Button("New 2D simulation"))
            state.launchMode = ui::LaunchMode::Sim2D;
        ImGui::SameLine();
        if (ImGui::Button("New 3D simulation"))
            state.launchMode = ui::LaunchMode::Sim3D;

        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "Preset previews");
        ImGui::Spacing();

        struct PreviewEntry { const char* title; const char* desc; ImVec4 color; };
        static constexpr PreviewEntry previews[] = {
            {"TON 618",   "Biggest confirmed quasar",         ImVec4(0.65f, 0.40f, 0.95f, 0.45f)},
            {"Sgr A*",    "Our local quasar",                 ImVec4(0.25f, 0.70f, 0.95f, 0.45f)},
            {"Custom BH", "Make your own black hole and galaxy", ImVec4(0.95f, 0.55f, 0.30f, 0.45f)}
        };
        ImVec2 avail = ImGui::GetContentRegionAvail();
        for (int i = 0; i < 3; ++i) {
            ImGui::PushStyleColor(ImGuiCol_Button,        previews[i].color);
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(
                previews[i].color.x + 0.1f, previews[i].color.y + 0.1f,
                previews[i].color.z + 0.1f, 0.6f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(
                previews[i].color.x + 0.05f, previews[i].color.y + 0.05f,
                previews[i].color.z + 0.05f, 0.7f));
            ImGui::BeginChild(previews[i].title,
                ImVec2((avail.x - 20.0f) / 3.0f, 140), true);
            ImGui::PushTextWrapPos();
            ImGui::TextColored(ImVec4(0.98f, 0.98f, 0.98f, 1.0f), "%s", previews[i].title);
            ImGui::TextDisabled("%s", previews[i].desc);
            ImGui::Spacing();
            if (ImGui::Button("Select preset", ImVec2(-1.0f, 0))) {
                state.statusText  = std::string("Preset selected: ") + previews[i].title;
                state.sessionState = 1;
            }
            ImGui::PopTextWrapPos();
            ImGui::EndChild();
            ImGui::PopStyleColor(3);
            if (i < 2) ImGui::SameLine();
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "Custom black hole");
        ImGui::InputText("Name", state.customBHName, sizeof(state.customBHName));
        ImGui::Checkbox("3D black hole", &state.customIs3D);
        ImGui::SliderFloat("Schwarzschild radius", &state.customBHRadius,
            0.1f, 50.0f, "%.2f Rs");
        if (ImGui::Button("Create custom black hole")) {
            state.statusText  = std::string("Created ") + state.customBHName;
            state.sessionState = 1;
        }
    }

} // namespace pages::simulation
