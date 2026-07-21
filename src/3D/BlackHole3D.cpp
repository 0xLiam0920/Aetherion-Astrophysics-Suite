// ============================================================
// BlackHole3D.cpp: Standalone window entry point.
// ============================================================
// All simulation, rendering, HUD, and input logic lives in
//   src/3D/bh3d_core.hpp
// shared with the Qt-embedded widget. This file only owns:
//   1. SFML window creation and main loop
//   2. CLI preset argument parsing
//   3. Routing window events into bh3d::* handlers
//   4. Mouse cursor grab/release (window-only concern)
// ============================================================

#include "platform.hpp"
#include <SFML/Window.hpp>
#include <SFML/Graphics/Image.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Window/Keyboard.hpp>
#include <SFML/Window/Mouse.hpp>

#include <cstring>
#include <iostream>

#include "bh3d_core.hpp"

#include "imgui.h"
#include "imgui-SFML.h"
#include "hud_imgui.hpp"

int main(int argc, char* argv[]) {
    // ── Profiles & CLI ────────────────────────────────────
    int profileIdx = 0; // default: TON 618
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--preset") == 0) {
            std::string id = argv[i + 1];
            if      (id == "ton618")   profileIdx = 0;
            else if (id == "sgra")     profileIdx = 1;
            else if (id == "3c273")    profileIdx = 2;
            else if (id == "j0529")    profileIdx = 3;
            else if (id == "phoenixa") profileIdx = 12;
            else std::cerr << "Unknown preset: " << id << ", using default.\n";
            break;
        }
    }

    // ── Build core state from selected profile ────────────
    auto bhProfilesProbe = profiles::allProfiles();                 // we have to call this to get the config for the initial camera setup, but we'll reload it properly in initProfiles 
    bh3d::State state(bhProfilesProbe[profileIdx].config.camera);   // right after. a bit wasteful but it only happens once and it's simpler than refactoring initProfiles to take a config arg or something.
    bh3d::initProfiles(state, profileIdx);

    // ── Window & GL context ───────────────────────────────
    sf::ContextSettings settings;
    settings.depthBits      = 24;
    settings.majorVersion   = 3;
    settings.minorVersion   = 3;
    settings.attributeFlags = sf::ContextSettings::Core;

    sf::RenderWindow window(
        sf::VideoMode({cfg::DEFAULT_WIDTH, cfg::DEFAULT_HEIGHT}),
        "GLSL Black Hole",
        sf::Style::Default,
        sf::State::Windowed,
        settings);
    window.setVerticalSyncEnabled(false);
    window.setKeyRepeatEnabled(false);
    if (!window.setActive(true)) {
        std::cerr << "Failed to activate OpenGL context.\n";
        return 1;
    }

    // ── Font ──────────────────────────────────────────────
    {
        sf::Font tmpFont;
        bh3d::initFont(state, tmpFont);
        if (!state.fontLoaded) std::cerr << "Warning: Could not load font for HUD.\n";
    }
    (void)window.setActive(true); // Re-activate the context: some platforms (notably macOS) need it active before loading textures, and initFont may have changed the active context.
    glViewport(0, 0, cfg::DEFAULT_WIDTH, cfg::DEFAULT_HEIGHT);

    // ── Resources / shaders / textures / bloom ────────────
    ResourceManager resources;
    if (!bh3d::initShaders(state, resources)) {
        std::cerr << "FATAL: No fragment shader compiled. Searched: "
                  << resources.paths() << "\n";
        sf::sleep(sf::seconds(3.0f));
        return 1;
    }
    state.quad = GLVertexArray::makeFullScreenQuad();
    bh3d::initTextures(state, resources);
    if (!state.bloom.init(cfg::DEFAULT_WIDTH, cfg::DEFAULT_HEIGHT)) {
        std::cerr << "Failed to initialize bloom pipeline!\n";
        return 1;
    }
    state.lastW = cfg::DEFAULT_WIDTH;
    state.lastH = cfg::DEFAULT_HEIGHT;
    state.actionKeys = bh3d::loadActionKeybinds();
    state.frameClock.restart();
    state.fpsClock.restart();

    // ── ImGui (SFML backend) ─────────────────────────────
    // Phase 1: only the standalone window uses ImGui. The preset selector
    // menu is drawn through ImGui; all other HUD elements still use the
    // GLBitmapFont path.
    if (!ImGui::SFML::Init(window)) {
        std::cerr << "Warning: ImGui::SFML::Init failed; HUD will fall back to "
                     "the bespoke renderer.\n";
    } else {
        state.useImGuiHud = true;
        hud_im::applyAetherionStyle();
        // Enable keyboard nav so arrow keys / enter work inside the preset menu
        // without the user clicking the window first.
        ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    }

    sf::Vector2i lastMousePos(0, 0);
    sf::Clock    imguiClock;

    // ════════════════════════════════════════════════════════
    //                       MAIN LOOP
    // ════════════════════════════════════════════════════════
    while (window.isOpen()) {
        const float dt = state.frameClock.restart().asSeconds();
        bh3d::tickFPS(state);

        // ── Events ────────────────────────────────────────
        while (auto evOpt = window.pollEvent()) {
            const auto& event = *evOpt;
            if (state.useImGuiHud) ImGui::SFML::ProcessEvent(window, event);

            const ImGuiIO& io = state.useImGuiHud ? ImGui::GetIO() : ImGui::GetIO();
            const bool imGuiKbd = state.useImGuiHud && io.WantCaptureKeyboard;
            const bool imGuiMouse = state.useImGuiHud && io.WantCaptureMouse;

            if (event.is<sf::Event::Closed>()) window.close();

            if (const auto* kp = event.getIf<sf::Event::KeyPressed>()) {
                // We always feed our own keystate (so WASD held during a menu
                // open keeps physics consistent on close), but we suppress
                // action-key dispatch when ImGui is consuming keyboard input -
                // except for the menu toggle (N) and release-mouse (Esc) keys
                // which the user expects to always work.
                state.keys.onKeyPressed(kp->code);
                const bool isToggle = (kp->code == state.actionKeys.nextProfile)
                                   || (kp->code == state.actionKeys.openMergerMenu)
                                   || (kp->code == state.actionKeys.releaseMouse);
                if (!imGuiKbd || isToggle) {
                    bh3d::onActionKey(state, kp->code);
                }
                if (kp->code == state.actionKeys.releaseMouse) {
                    window.setMouseCursorGrabbed(false);
                    window.setMouseCursorVisible(true);
                }
            }
            if (const auto* kr = event.getIf<sf::Event::KeyReleased>()) {
                state.keys.onKeyReleased(kr->code);
            }
            if (const auto* resized = event.getIf<sf::Event::Resized>()) {
                glViewport(0, 0, (GLsizei)resized->size.x, (GLsizei)resized->size.y);
            }
            if (const auto* mb = event.getIf<sf::Event::MouseButtonPressed>()) {
                if (mb->button == sf::Mouse::Button::Right) {
                    if (!imGuiMouse) {
                        state.looking = true;
                        window.setMouseCursorGrabbed(true);
                        window.setMouseCursorVisible(false);
                        lastMousePos = mb->position;
                    }
                } else if (mb->button == sf::Mouse::Button::Left) {
                    if (!imGuiMouse) {
                        const int sw = std::max(1, (int)window.getSize().x);
                        const int sh = std::max(1, (int)window.getSize().y);
                        // Body-pick takes priority; falls through to overlay buttons if nothing hit
                        if (!bh3d::onBodyClick(state, float(mb->position.x), float(mb->position.y), sw, sh))
                            bh3d::onLeftClickOverlays(state, float(mb->position.x), float(mb->position.y));
                    }
                }
            }
            if (const auto* mb = event.getIf<sf::Event::MouseButtonReleased>()) {
                if (mb->button == sf::Mouse::Button::Right) {
                    state.looking = false;
                    window.setMouseCursorGrabbed(false);
                    window.setMouseCursorVisible(true);
                }
            }
            if (const auto* mm = event.getIf<sf::Event::MouseMoved>()) {
                if (state.looking) {
                    sf::Vector2i delta = mm->position - lastMousePos;
                    lastMousePos = mm->position;
                    state.camera.onMouseDelta(float(delta.x), float(delta.y));
                }
            }
        }

        // ── Physics, snapshot, render ─────────────────────
        bh3d::tickPhysics(state, dt);
        const int w = std::max(1, (int)window.getSize().x);
        const int h = std::max(1, (int)window.getSize().y);
        bh3d::buildSnapshot(state, w, h, dt);
        bh3d::renderScene(state, w, h);
        bh3d::renderHUD(state, w, h);

        // ── ImGui frame (preset menu lives here) ─────────
        if (state.useImGuiHud) {
            ImGui::SFML::Update(window, imguiClock.restart());
            hud_im::drawAll(state);
            ImGui::SFML::Render(window);
        }

        window.display();
        while (glGetError() != GL_NO_ERROR) {} // clears any GL errors that may have been generated 
    }

    if (state.useImGuiHud) ImGui::SFML::Shutdown(window);

    std::cerr << "All GL resources cleaned automatically.\n";
    return 0;
}