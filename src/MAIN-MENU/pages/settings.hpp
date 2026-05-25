/*  pages/settings.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::settings {

    inline void draw(ui::AppState& state) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Settings");
        ImGui::TextDisabled("Modify keybinds and common interface options.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::InputText("Move camera",       state.keybindMove,          sizeof(state.keybindMove));
        ImGui::InputText("Toggle preset",     state.keybindTogglePreset,  sizeof(state.keybindTogglePreset));
        ImGui::InputText("Export quick key",  state.keybindExport,        sizeof(state.keybindExport));
        ImGui::Checkbox("Show tooltips",      &state.showTooltips);
        ImGui::Spacing();
        if (ImGui::Button("Save settings"))
            state.statusText = "Settings updated";
    }

} // namespace pages::settings
