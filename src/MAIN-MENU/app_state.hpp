/*
    app_state.hpp
    Shared state types for the Aetherion main menu.
    All pages and the registry depend on this header — keep it minimal.
*/
#pragma once

#include <string>

namespace ui {

    enum class LaunchMode {
        None,
        Sim2D,
        Sim3D
    };

    // One entry per page. The registry uses these as identifiers.
    enum class Tab {
        Overview,
        Simulation,
        DataAnalysis,
        ObjectLibrary,
        Export,
        Settings,
        Breakdown
    };

    struct AppState {
        LaunchMode launchMode  = LaunchMode::None;
        bool       wantsExit   = false;
        std::string statusText = "ready";
        int  sessionState      = 0;   // 0 = no session, 1 = session loaded
        int  modelState        = 2;   // 0 = complete, 1/2 = pending
        Tab  currentTab        = Tab::Overview;

        // Simulation page
        char loadPath[256]  = "";
        char customBHName[64] = "Custom Black Hole";
        float customBHRadius  = 1.0f;
        bool  customIs3D      = false;

        // Export page
        char exportFolder[256] = "saves/";
        int  exportTypeIndex   = 0;

        // Settings page
        char keybindMove[32]         = "WASD";
        char keybindTogglePreset[32] = "T";
        char keybindExport[32]       = "X";
        bool showTooltips            = true;
    };

} // namespace ui
