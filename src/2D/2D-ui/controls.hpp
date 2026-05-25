#pragma once
#include <SFML/Graphics.hpp>
#include "../2D-simulation/simulation.hpp"
#include "../2D-utils/presets_2d.hpp"
#include "../2D-utils/key_config_2d.hpp"
#include "../2D-physics/units.hpp"
#include <algorithm>
#include <string>

// Runtime toggles controlled by the user.
struct UIState {
    bool   showRays          = true;
    bool   showPhotonSphere  = true;
    bool   showGalaxySystem  = true;
    bool   showInfluenceZones = true;
    bool   paused            = false;
    bool   presetActive      = false;
    bool   accuratePhysics   = true;
    bool   showTimeDilation  = false;
    int    presetIdx         = 0;
    double timeScale         = 1.0;
    double defaultPixelsPerM = 60.0;

    // Research mode toggles
    bool   showDataPanel     = false;
    bool   showTimeDilationMap = false;
    bool   showCaustics      = false;
    bool   showNumericalError = false;
    bool   showControlsPanel  = true;   // toggleable controls list
    bool   highResLensing    = false;   // log-spaced extra rays near b_crit
    std::string notification;
    int    notificationTimer = 0;

    static constexpr double presetHorizonPixelsTarget = 120.0;

    // ── Merger menu ──────────────────────────────────────────────────────────
    struct MergerMenuState {
        bool        open             = false;
        int         selectedIdx      = 0;       // 0..NUM_BH2D_PRESETS-1 = preset, NUM_BH2D_PRESETS = Custom
        bool        inputtingCustom  = false;   // user is typing a custom mass
        std::string customInput      = "10.0";  // raw digit string
        double      customMassSolar  = 10.0;    // last parsed custom mass
    } mergerMenu;
};

// Process a single SFML event and update UI / simulation state.
// Pass a KeyConfig2D to support remapped keybinds; defaults match original layout.
inline void handleInput(
    const sf::Event&  ev,
    Simulation&       sim,
    UIState&          ui,
    unsigned int      windowHeight,
    const KeyConfig2D& cfg = KeyConfig2D{})
{
    // ── Text entry for custom-mass field in merger menu ─────────────────────
    if (ui.mergerMenu.open && ui.mergerMenu.inputtingCustom) {
        if (auto* te = ev.getIf<sf::Event::TextEntered>()) {
            char32_t c = te->unicode;
            if (c >= '0' && c <= '9') {
                ui.mergerMenu.customInput += (char)c;
            } else if (c == '.') {
                if (ui.mergerMenu.customInput.find('.') == std::string::npos)
                    ui.mergerMenu.customInput += '.';
            } else if (c == 8 || c == 127) { // backspace / delete
                if (!ui.mergerMenu.customInput.empty())
                    ui.mergerMenu.customInput.pop_back();
            } else if (c == 13) { // Enter — confirm
                try {
                    double m = std::stod(ui.mergerMenu.customInput);
                    if (m > 0.0) {
                        ui.mergerMenu.customMassSolar = m;
                        sim.startMerger(m, windowHeight);
                        ui.mergerMenu.open            = false;
                        ui.mergerMenu.inputtingCustom = false;
                        ui.notification      = "Merger with custom BH initiated!";
                        ui.notificationTimer = 150;
                    }
                } catch (...) {
                    ui.notification      = "Invalid mass — enter a positive number";
                    ui.notificationTimer = 120;
                }
            }
        }
    }

    if (!ev.is<sf::Event::KeyPressed>()) return;

    const auto *key    = ev.getIf<sf::Event::KeyPressed>();
    const sf::Keyboard::Key code = key->code;
    const double bh_isco = sim.bh.metric.isco();

    // ── Merger menu navigation (intercepts all keys while open) ─────────────
    if (ui.mergerMenu.open) {
        if (code == sf::Keyboard::Key::Escape) {
            ui.mergerMenu.open            = false;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Up) {
            int total = NUM_BH2D_PRESETS + 1;
            ui.mergerMenu.selectedIdx = (ui.mergerMenu.selectedIdx - 1 + total) % total;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Down) {
            int total = NUM_BH2D_PRESETS + 1;
            ui.mergerMenu.selectedIdx = (ui.mergerMenu.selectedIdx + 1) % total;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Enter) {
            if (ui.mergerMenu.selectedIdx == NUM_BH2D_PRESETS) {
                // Custom option — begin text entry
                if (!ui.mergerMenu.inputtingCustom) {
                    ui.mergerMenu.inputtingCustom = true;
                    ui.mergerMenu.customInput     = "10.0";
                } else {
                    // Confirm what was typed
                    try {
                        double m = std::stod(ui.mergerMenu.customInput);
                        if (m > 0.0) {
                            ui.mergerMenu.customMassSolar = m;
                            sim.startMerger(m, windowHeight);
                            ui.mergerMenu.open            = false;
                            ui.mergerMenu.inputtingCustom = false;
                            ui.notification      = "Merger with custom BH initiated!";
                            ui.notificationTimer = 150;
                        }
                    } catch (...) {
                        ui.notification      = "Invalid mass — enter a positive number";
                        ui.notificationTimer = 120;
                    }
                }
            } else {
                double massSolar = BH2D_PRESETS[ui.mergerMenu.selectedIdx].massSolar;
                sim.startMerger(massSolar, windowHeight);
                ui.mergerMenu.open = false;
                ui.notification = std::string("Merger with ") +
                                  BH2D_PRESETS[ui.mergerMenu.selectedIdx].name + " initiated!";
                ui.notificationTimer = 150;
            }
        }
        return; // swallow all keypresses while menu is open
    }

    // ── Simulation ────────────────────────────────────────────────────────────
    if (code == cfg.pause) {
        ui.paused = !ui.paused;

    } else if (code == cfg.togglePreset) {
        ui.presetActive = !ui.presetActive;
        if (ui.presetActive) {
            sim.bh.metric.M       = units::solarMassToGeomMeters(BH2D_PRESETS[ui.presetIdx].massSolar);
            sim.params.pixelsPerM = UIState::presetHorizonPixelsTarget / (2.0 * sim.bh.metric.M);
            sim.tidalRadiusM      = BH2D_PRESETS[ui.presetIdx].zones.tidalDisruptionM;
            sim.reinitBodies();
            if (ui.showGalaxySystem && BH2D_PRESETS[ui.presetIdx].isGalacticCenter)
                sim.spawnGalaxySystem(ui.presetIdx);
        } else {
            sim.bh.metric.M       = 1.0;
            sim.params.pixelsPerM = ui.defaultPixelsPerM;
            sim.clearGalaxySystem();
            sim.reinitBodies();
            sim.rebuildPhotons(windowHeight);
        }

    } else if (code == cfg.nextPreset) {
        if (ui.presetActive) {
            ui.presetIdx = (ui.presetIdx + 1) % NUM_BH2D_PRESETS;
            sim.bh.metric.M       = units::solarMassToGeomMeters(BH2D_PRESETS[ui.presetIdx].massSolar);
            sim.params.pixelsPerM = UIState::presetHorizonPixelsTarget / (2.0 * sim.bh.metric.M);
            sim.tidalRadiusM      = BH2D_PRESETS[ui.presetIdx].zones.tidalDisruptionM;
            sim.reinitBodies();
            if (ui.showGalaxySystem && BH2D_PRESETS[ui.presetIdx].isGalacticCenter)
                sim.spawnGalaxySystem(ui.presetIdx);
            else
                sim.clearGalaxySystem();
            sim.rebuildPhotons(windowHeight);
        }

    } else if (code == cfg.prevPreset) {
        if (ui.presetActive) {
            ui.presetIdx = (ui.presetIdx - 1 + NUM_BH2D_PRESETS) % NUM_BH2D_PRESETS;
            sim.bh.metric.M       = units::solarMassToGeomMeters(BH2D_PRESETS[ui.presetIdx].massSolar);
            sim.params.pixelsPerM = UIState::presetHorizonPixelsTarget / (2.0 * sim.bh.metric.M);
            sim.tidalRadiusM      = BH2D_PRESETS[ui.presetIdx].zones.tidalDisruptionM;
            sim.reinitBodies();
            if (ui.showGalaxySystem && BH2D_PRESETS[ui.presetIdx].isGalacticCenter)
                sim.spawnGalaxySystem(ui.presetIdx);
            else
                sim.clearGalaxySystem();
            sim.rebuildPhotons(windowHeight);
        }

    } else if (code == cfg.reset) {
        sim.reset();
        ui.presetActive       = false;
        sim.params.pixelsPerM = ui.defaultPixelsPerM;
        sim.rebuildPhotons(windowHeight);

    // ── Orbital control ───────────────────────────────────────────────────
    } else if (code == cfg.zoomIn) {
        if (ui.presetActive) {
            sim.params.pixelsPerM *= 1.15;
        } else {
            sim.bh.metric.M *= 1.10;
            sim.reinitBodies();
            sim.rebuildPhotons(windowHeight);
        }

    } else if (code == cfg.zoomOut) {
        if (ui.presetActive) {
            sim.params.pixelsPerM = std::max(1e-20, sim.params.pixelsPerM / 1.15);
        } else {
            sim.bh.metric.M *= 0.9;
            sim.reinitBodies();
            sim.rebuildPhotons(windowHeight);
        }

    } else if (code == cfg.orbitIn) {
        if (!sim.bodies.empty()) {
            auto& b = sim.bodies[0];
            b.nominalA = std::max(bh_isco + 0.5 * sim.bh.metric.M,
                                  b.nominalA - 0.5 * sim.bh.metric.M);
            if (b.initFromKeplerian(sim.bh.metric, b.nominalA, b.nominalEcc)) {
                ui.notification = "Periapsis clamped to ISCO";
                ui.notificationTimer = 120;
            }
        }

    } else if (code == cfg.orbitOut) {
        if (!sim.bodies.empty()) {
            auto& b = sim.bodies[0];
            b.nominalA += 0.5 * sim.bh.metric.M;
            if (b.initFromKeplerian(sim.bh.metric, b.nominalA, b.nominalEcc)) {
                ui.notification = "Periapsis clamped to ISCO";
                ui.notificationTimer = 120;
            }
        }

    } else if (code == cfg.eccDecrease) {
        if (!sim.bodies.empty()) {
            auto& b = sim.bodies[0];
            b.nominalEcc = std::max(0.0, b.nominalEcc - 0.05);
            if (b.initFromKeplerian(sim.bh.metric, b.nominalA, b.nominalEcc)) {
                ui.notification = "Periapsis clamped to ISCO";
                ui.notificationTimer = 120;
            }
        }

    } else if (code == cfg.eccIncrease) {
        if (!sim.bodies.empty()) {
            auto& b = sim.bodies[0];
            b.nominalEcc = std::min(0.9, b.nominalEcc + 0.05);
            if (b.initFromKeplerian(sim.bh.metric, b.nominalA, b.nominalEcc)) {
                ui.notification = "Periapsis clamped to ISCO";
                ui.notificationTimer = 120;
            }
        }

    // ── Visualization ─────────────────────────────────────────────────────
    } else if (code == cfg.toggleRays) {
        ui.showRays = !ui.showRays;

    } else if (code == cfg.togglePhoton) {
        ui.showPhotonSphere = !ui.showPhotonSphere;

    } else if (code == cfg.toggleGalaxy) {
        ui.showGalaxySystem = !ui.showGalaxySystem;
        if (ui.presetActive) {
            if (ui.showGalaxySystem && BH2D_PRESETS[ui.presetIdx].isGalacticCenter)
                sim.spawnGalaxySystem(ui.presetIdx);
            else
                sim.clearGalaxySystem();
        }

    } else if (code == cfg.toggleInfluence) {
        ui.showInfluenceZones = !ui.showInfluenceZones;

    } else if (code == cfg.toggleHighResLensing) {
        ui.highResLensing = !ui.highResLensing;
        sim.rebuildPhotons(windowHeight, ui.highResLensing);
        ui.notification = ui.highResLensing ? "High-res lensing ON" : "High-res lensing OFF";
        ui.notificationTimer = 90;

    } else if (code == cfg.speedUp) {
        ui.timeScale = std::min(64.0, ui.timeScale * 2.0);
        {
            std::ostringstream os;
            if (ui.timeScale >= 1.0) os << (int)ui.timeScale;
            else os << std::fixed << std::setprecision(2) << ui.timeScale;
            ui.notification = "Speed: " + os.str() + "x";
        }
        ui.notificationTimer = 90;

    } else if (code == cfg.speedDown) {
        ui.timeScale = std::max(0.125, ui.timeScale / 2.0);
        {
            std::ostringstream os;
            if (ui.timeScale >= 1.0) os << (int)ui.timeScale;
            else os << std::fixed << std::setprecision(2) << ui.timeScale;
            ui.notification = "Speed: " + os.str() + "x";
        }
        ui.notificationTimer = 90;

    } else if (code == cfg.tempDown) {
        sim.gasTemperatureK = std::max(1e4, sim.gasTemperatureK / 2.0);

    } else if (code == cfg.tempUp) {
        sim.gasTemperatureK = std::min(1e12, sim.gasTemperatureK * 2.0);

    // ── Research ──────────────────────────────────────────────────────────
    } else if (code == cfg.toggleData) {
        ui.showDataPanel = !ui.showDataPanel;

    } else if (code == cfg.toggleDilation) {
        ui.showTimeDilationMap = !ui.showTimeDilationMap;

    } else if (code == cfg.toggleCaustics) {
        ui.showCaustics = !ui.showCaustics;

    } else if (code == cfg.toggleError) {
        ui.showNumericalError = !ui.showNumericalError;

    } else if (code == cfg.cycleBody) {
        sim.cycleSelectedBody();

    } else if (code == cfg.toggleControls) {
        ui.showControlsPanel = !ui.showControlsPanel;

    } else if (code == cfg.exportData) {
        if (cfg.exportFormat == "FITS")
            ui.notification = sim.exportAllFITS();
        else if (cfg.exportFormat == "Binary")
            ui.notification = sim.exportAllBinary();
        else
            ui.notification = sim.exportAllCSV();
        // Give extra time to read the path (or error message)
        ui.notificationTimer = 360;

    // ── Test scenarios ────────────────────────────────────────────────────
    } else if (code == cfg.testIsco) {
        sim.startISCOTest();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = "ISCO validation: 5M/6M/7M test";
        ui.notificationTimer = 120;

    } else if (code == cfg.testPhoton) {
        sim.startPhotonSphereTest();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = "Photon sphere test active";
        ui.notificationTimer = 120;

    } else if (code == cfg.testInfall) {
        sim.startRadialInfall();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = "Radial infall from 20M";
        ui.notificationTimer = 120;

    } else if (code == cfg.testTidal) {
        sim.startTidalDisruption();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = "Tidal disruption demo";
        ui.notificationTimer = 120;

    } else if (code == cfg.testPulsar) {
        sim.startPulsarOrbital();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = "Pulsar orbital sim — enable preset (T) for physical GW values";
        ui.notificationTimer = 150;

    } else if (code == cfg.mergerMenu) {
        ui.mergerMenu.open            = !ui.mergerMenu.open;
        ui.mergerMenu.inputtingCustom = false;
        ui.mergerMenu.selectedIdx     = 0;
    }
}
