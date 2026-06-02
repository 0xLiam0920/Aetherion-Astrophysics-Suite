/*  pages/object_library.hpp  */
#pragma once
#include "../app_state.hpp"
#include "../libraries.hpp"

namespace pages::object_library {

    inline void draw(ui::AppState& /*state*/) {
        ImGui::SetCursorPosX(30);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Object Library");
        ImGui::TextDisabled("Browse custom presets for black holes and simulation objects.");
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        struct LibraryEntry { const char* title; const char* subtitle; const char* tag; ImVec4 color; };
        static constexpr LibraryEntry entries[] = {
            {"2D Micro BH",    "Thin-disk research preset",   "Black Hole", ImVec4(0.40f, 0.70f, 0.95f, 1.0f)},
            {"3D Rapid Spin",  "High-magnification disk",     "Black Hole", ImVec4(0.80f, 0.45f, 0.95f, 1.0f)},
            {"Pulsar Beacon",  "Fast-rotating neutron star",  "Object",     ImVec4(0.95f, 0.55f, 0.30f, 1.0f)},
            {"Nebula Cloud",   "Diffuse gas cluster",         "Object",     ImVec4(0.30f, 0.85f, 0.55f, 1.0f)},
            {"Asteroid Swarm", "Infalling debris field",      "Object",     ImVec4(0.95f, 0.75f, 0.35f, 1.0f)}
        };
        for (const auto& e : entries) {
            ImGui::PushStyleColor(ImGuiCol_Text, e.color);
            ImGui::Text("%s", e.title);
            ImGui::PopStyleColor();
            ImGui::TextDisabled("%s • %s", e.subtitle, e.tag);
            ImGui::Spacing();
        }
    }

} // namespace pages::object_library
