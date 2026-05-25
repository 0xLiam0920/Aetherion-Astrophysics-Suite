/*  pages/export.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::export_tab {

    inline void draw(ui::AppState& state) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Export");
        ImGui::TextDisabled("Choose export folder and file type for CSV/FITS output.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        ImGui::InputText("Export folder", state.exportFolder, sizeof(state.exportFolder));
        const char* types[] = {"CSV", "FITS", "Binary"};
        ImGui::Combo("File type", &state.exportTypeIndex, types, IM_ARRAYSIZE(types));
        ImGui::Spacing();
        ImGui::TextDisabled("All exported data will be stored under the selected folder by default.");
        ImGui::Spacing();
        if (ImGui::Button("Apply export settings")) {
            state.statusText = std::string("Export set to ")
                + state.exportFolder + " (" + types[state.exportTypeIndex] + ")";
        }
    }

} // namespace pages::export_tab
