/*  pages/data_analysis.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::data_analysis {

    inline void draw(ui::AppState& /*state*/) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Data Analysis");
        ImGui::TextDisabled("Coming soon, this space will host analysis controls and charts.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        ImGui::Dummy(ImVec2(0, 340));
        ImGui::TextDisabled("Blank workspace.");
    }

} // namespace pages::data_analysis
