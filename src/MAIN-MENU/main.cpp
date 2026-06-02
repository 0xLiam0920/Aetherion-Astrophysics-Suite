/*
    Author: Liam Nayyar
    Purpose:
    Main entry point for the Aetherion launcher.

    Adding a new page
    ─────────────────
    1. Create  src/MAIN-MENU/pages/my_page.hpp  with a  draw(AppState&)  function.
    2. Add an entry to the Tab enum in  app_state.hpp.
    3. #include and register it in  page_registry.hpp  (one line each).
    Done, sidebar and content update automatically.
*/

#include "libraries.hpp"
#include "launcher.hpp"
#include "page_registry.hpp"  // pulls in app_state + all page headers

namespace ui {

    // ── Theme ──────────────────────────────────────────────────────────────────
    void applyTheme() {
        ImGuiStyle& style = ImGui::GetStyle();
        style.WindowRounding = 0.0f;
        style.FrameRounding  = 4.0f;
        style.PopupRounding  = 4.0f;
        style.GrabRounding   = 4.0f;

        style.Colors[ImGuiCol_WindowBg]       = ImVec4(0.08f, 0.08f, 0.11f, 1.0f);
        style.Colors[ImGuiCol_FrameBg]        = ImVec4(0.12f, 0.12f, 0.16f, 0.8f);
        style.Colors[ImGuiCol_FrameBgHovered] = ImVec4(0.15f, 0.15f, 0.22f, 0.9f);
        style.Colors[ImGuiCol_FrameBgActive]  = ImVec4(0.18f, 0.18f, 0.25f, 1.0f);
        style.Colors[ImGuiCol_Button]         = ImVec4(0.12f, 0.35f, 0.60f, 0.85f);
        style.Colors[ImGuiCol_ButtonHovered]  = ImVec4(0.16f, 0.45f, 0.80f, 1.0f);
        style.Colors[ImGuiCol_ButtonActive]   = ImVec4(0.10f, 0.30f, 0.55f, 1.0f);
        style.Colors[ImGuiCol_Text]           = ImVec4(0.93f, 0.93f, 0.96f, 1.0f);
        style.Colors[ImGuiCol_TextDisabled]   = ImVec4(0.55f, 0.55f, 0.65f, 1.0f);
        style.Colors[ImGuiCol_Separator]      = ImVec4(0.25f, 0.25f, 0.35f, 0.5f);
        style.ItemSpacing      = ImVec2(8.0f, 8.0f);
        style.ItemInnerSpacing = ImVec2(6.0f, 6.0f);
        style.FramePadding     = ImVec2(8.0f, 6.0f);
        style.WindowPadding    = ImVec2(0.0f, 0.0f);
    }

    // ── Sidebar ────────────────────────────────────────────────────────────────
    // Iterates kPages; prints a section heading whenever page.section is set.
    void drawSidebar(AppState& state) {
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(280, ImGui::GetIO().DisplaySize.y));
        ImGui::Begin("##Sidebar", nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);

        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.08f, 1.0f));
        ImGui::Spacing();
        ImGui::SetCursorPosX(20);
        ImGui::TextColored(ImVec4(0.93f, 0.93f, 0.96f, 1.0f), "Aetherion");
        ImGui::SetCursorPosX(20);
        ImGui::TextDisabled("ASTROPHYSICS SUITE");
        ImGui::SetCursorPosX(20);
        ImGui::Separator();
        ImGui::Spacing();

        for (const auto& page : getPages()) {
            if (page.section != nullptr) {
                ImGui::SetCursorPosX(20);
                ImGui::Spacing();
                ImGui::SetCursorPosX(20);
                ImGui::Separator();
                ImGui::Spacing();
                ImGui::SetCursorPosX(20);
                ImGui::TextColored(ImVec4(0.55f, 0.55f, 0.65f, 1.0f), "%s", page.section);
            }
            ImGui::SetCursorPosX(30);
            bool selected = state.currentTab == page.id;
            if (ImGui::Selectable(page.label, selected, 0, ImVec2(0, 20)))
                state.currentTab = page.id;
        }

        ImGui::SetCursorPos(ImVec2(20, ImGui::GetWindowHeight() - 60));
        ImGui::Separator();
        ImGui::SetCursorPosX(20);
        ImGui::TextDisabled("v0.0.1 • open source");
        ImGui::PopStyleColor();
        ImGui::End();
    }

    // ── Main content ───────────────────────────────────────────────────────────
    // Finds the active page in kPages and calls its draw function.
    void drawMainContent(AppState& state) {
        ImGui::SetNextWindowPos(ImVec2(280, 0));
        ImGui::SetNextWindowSize(ImVec2(
            ImGui::GetIO().DisplaySize.x - 280,
            ImGui::GetIO().DisplaySize.y));
        ImGui::Begin("##MainContent", nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);

        for (const auto& page : getPages()) {
            if (page.id == state.currentTab) {
                page.draw(state);
                break;
            }
        }

        ImGui::End();
    }

    // ── Status bar ─────────────────────────────────────────────────────────────
    void drawStatusBar(AppState& state) {
        ImGui::SetNextWindowPos(ImVec2(0, ImGui::GetIO().DisplaySize.y - 30));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x, 30));
        ImGui::Begin("##StatusBar", nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
            ImGuiWindowFlags_NoScrollbar);

        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.06f, 0.06f, 0.09f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Border,   ImVec4(0.25f, 0.25f, 0.35f, 1.0f));

        ImGui::SetCursorPosX(10);
        ImGui::TextColored(ImVec4(0.30f, 0.90f, 0.30f, 1.0f),
            u8"● %s", state.statusText.c_str());

        ImGui::SameLine(150);
        ImGui::TextDisabled(state.sessionState == 1
            ? "| session loaded" : "| no session loaded");

        ImGui::SameLine(300);
        ImGui::TextDisabled("| blueshift model - %s",
            state.modelState == 0 ? "complete" : "pending");

        ImGui::PopStyleColor(2);
        ImGui::End();
    }

    void drawUI(AppState& state) {
        drawSidebar(state);
        drawMainContent(state);
        drawStatusBar(state);
    }

} // namespace ui

// ── Entry point ────────────────────────────────────────────────────────────────
int main() {
    constexpr unsigned int WIN_W = 1200;
    constexpr unsigned int WIN_H = 800;

    sf::RenderWindow window(
        sf::VideoMode(sf::Vector2u(WIN_W, WIN_H)),
        "Aetherion Astrophysics Suite",
        sf::Style::Close);
    window.setFramerateLimit(60);

    if (!ImGui::SFML::Init(window)) {
        std::cerr << "Error: Failed to initialize ImGui-SFML.\n";
        exit(1); 
    }

    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.0f;

    sf::Clock       deltaClock;
    ui::AppState    state;

    while (window.isOpen()) {
        while (auto evOpt = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *evOpt);
            if (evOpt->is<sf::Event::Closed>())
                window.close();
        }

        sf::Time dt = deltaClock.restart();
        ImGui::SFML::Update(window, dt);

        ui::applyTheme();
        ui::drawUI(state);

        // ── Handle launch request ──────────────────────────────────────────────
        if (state.launchMode != ui::LaunchMode::None) {
            const char* exec = nullptr;
            if (state.launchMode == ui::LaunchMode::Sim2D) exec = "blackhole-2D";
            if (state.launchMode == ui::LaunchMode::Sim3D) exec = "blackhole-3D";
            if (exec && !launcher::launchSimulator(exec))
                state.statusText = "launch failed";
            state.launchMode = ui::LaunchMode::None;
        }

        if (state.wantsExit)
            window.close();

        window.clear(sf::Color(13, 13, 18));
        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}