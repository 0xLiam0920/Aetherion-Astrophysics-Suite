/*
A page dedicated to how our simulations work in various breakdowns.
First, the user will be prompted with a menu to choose the field in which they want a breakdown of,
then we ask for 2d explanation, 3d explanation, or both,
and then we show the breakdown in various forms, i.e, graphs, equations, logic, regular explanations ofc, etc
*/

/*  pages/breakdown.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"
namespace pages::selection {
    /*------------------------- FRONTEND OF SELECTION PAGE -------------------------*/
    inline void draw(ui::AppState& state) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Please select a field for breakdown:");
        ImGui::Spacing();
        ImGui::BeginChild("##selection", ImVec2(0, 600), true);
        
        ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), "Core Physics");
        if (ImGui::Button("General Relativity", ImVec2(-1, 0))) {
            state.current_breakdown_field = "General Relativity";
            state.current_page = "breakdown";
        }
        if (ImGui::Button("Orbital Mechanics", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Orbital Mechanics";
            state.current_page = "breakdown";
        }
        if (ImGui::Button("Accretion Physics", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Accretion Physics";
            state.current_page = "breakdown";
        }
        
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), "Radiation & Observations");
        if (ImGui::Button("Radiation & Optics", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Radiation & Optics";
            state.current_page = "breakdown";
        }
        if (ImGui::Button("Thermodynamics", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Thermodynamics";
            state.current_page = "breakdown";
        }
        
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), "Advanced Topics");
        if (ImGui::Button("Quantum Effects", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Quantum Effects";
            state.current_page = "breakdown";
        }
        if (ImGui::Button("Electromagnetism", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Electromagnetism";
            state.current_page = "breakdown";
        }
        
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), "Computational Methods");
        if (ImGui::Button("Numerical Methods", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Numerical Methods";
            state.current_page = "breakdown";
        }
        if (ImGui::Button("Signal Processing & Data", ImVec2(-1, 0))) {
            state.current_breakdown_field = "Signal Processing & Data";
            state.current_page = "breakdown";
        }
        
        ImGui::EndChild();
    }
    /*------------------------- BACKEND OF SELECTION PAGE -------------------------*/
}
namespace pages::breakdown {

    inline void draw(ui::AppState& /*state*/) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Breakdown");
        ImGui::TextDisabled("Coming soon, this space will host detailed breakdowns of our simulations.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        ImGui::Dummy(ImVec2(0, 340));
        ImGui::TextDisabled("Blank workspace.");
    }

} // namespace pages::breakdown