#pragma once
#include <SFML/Graphics.hpp>
#include "../2D-simulation/simulation.hpp"
#include "../2D-utils/presets_2d.hpp"
#include "../2D-utils/key_config_2d.hpp"
#include "../2D-physics/units.hpp"
#include "../../platform.hpp"
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
    bool   lightMode         = false;   // light theme for the 2D HUD overlays
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

    // ── Custom-body creator menu ─────────────────────────────────────────────
    // Lets the user spawn a free-form body around the primary BH and save it
    // as a reusable preset. Fields 0..3 are input rows (name, type, sm, ecc).
    // focusedField >= kCustomBodyInputCount selects a saved preset row.
    static constexpr int kCustomBodyInputCount = 4;
    struct CustomBodyMenuState {
        bool        open            = false;
        int         focusedField    = 0;
        int         typeIdx         = 0;        // index into customBodyTypes()
        std::string nameInput;                  // optional; auto-named if empty
        std::string smInput         = "20.0";
        std::string eccInput        = "0.3";
    } customBodyMenu;

    // Fixed type list for the custom-body creator. Order is independent of the
    // GalaxyBodyType enum order so the menu can present a curated lineup.
    static constexpr int kCustomBodyTypeCount = 7;
    static const GalaxyBodyType* customBodyTypes() {
        static const GalaxyBodyType kTypes[kCustomBodyTypeCount] = {
            GalaxyBodyType::Star,
            GalaxyBodyType::GasCloud,
            GalaxyBodyType::StellarCluster,
            GalaxyBodyType::DwarfGalaxy,
            GalaxyBodyType::NeutronStar,
            GalaxyBodyType::WhiteDwarf,
            GalaxyBodyType::CompanionStar,
        };
        return kTypes;
    }
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
            } else if (c == 13) { // Enter, confirm
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
                    ui.notification      = "Invalid mass, enter a positive number";
                    ui.notificationTimer = 120;
                }
            }
        }
    }

    // ── Text entry for the custom-body creator (name + sm + ecc fields) ──
    if (ui.customBodyMenu.open) {
        const int f = ui.customBodyMenu.focusedField;
        if (f == 0 || f == 2 || f == 3) {
            if (auto* te = ev.getIf<sf::Event::TextEntered>()) {
                std::string* buf = nullptr;
                if      (f == 0) buf = &ui.customBodyMenu.nameInput;
                else if (f == 2) buf = &ui.customBodyMenu.smInput;
                else             buf = &ui.customBodyMenu.eccInput;
                char32_t c = te->unicode;
                if (f == 0) {
                    // Name field: any printable ASCII except tab/newline.
                    if (c >= 0x20 && c < 0x7f && c != '\t' && buf->size() < 40) {
                        *buf += (char)c;
                    } else if (c == 8 || c == 127) {
                        if (!buf->empty()) buf->pop_back();
                    }
                } else {
                    // Numeric fields.
                    if (c >= '0' && c <= '9') {
                        *buf += (char)c;
                    } else if (c == '.') {
                        if (buf->find('.') == std::string::npos) *buf += '.';
                    } else if (c == 8 || c == 127) {
                        if (!buf->empty()) buf->pop_back();
                    }
                }
            }
        }
    }

    if (!ev.is<sf::Event::KeyPressed>()) return;

    const auto *key    = ev.getIf<sf::Event::KeyPressed>();
    const sf::Keyboard::Key code = key->code;
    const double bh_isco = sim.bh.metric.isco();

    // ── Merger menu navigation (intercepts all keys while open) ─────────────
    // Row layout:
    //   [0 .. NUM_BH2D_PRESETS-1]            BH presets
    //   [NUM_BH2D_PRESETS .. +N_SEC-1]       Compact / stellar secondaries
    //   [last]                               Custom BH (mass typed by user)
    if (ui.mergerMenu.open) {
        const int kBHEnd     = NUM_BH2D_PRESETS;                          // exclusive
        const int kSecEnd    = kBHEnd + NUM_MERGER_SECONDARY_PRESETS;      // exclusive
        const int kCustomIdx = kSecEnd;                                    // single custom row
        const int kTotalRows = kCustomIdx + 1;
        if (code == sf::Keyboard::Key::Escape) {
            ui.mergerMenu.open            = false;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Up) {
            ui.mergerMenu.selectedIdx = (ui.mergerMenu.selectedIdx - 1 + kTotalRows) % kTotalRows;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Down) {
            ui.mergerMenu.selectedIdx = (ui.mergerMenu.selectedIdx + 1) % kTotalRows;
            ui.mergerMenu.inputtingCustom = false;
        } else if (code == sf::Keyboard::Key::Enter) {
            const int idx = ui.mergerMenu.selectedIdx;
            if (idx == kCustomIdx) {
                // Custom option, begin text entry
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
                        ui.notification      = "Invalid mass, enter a positive number";
                        ui.notificationTimer = 120;
                    }
                }
            } else if (idx >= kBHEnd && idx < kSecEnd) {
                // Compact / stellar secondary
                const auto& sec = MERGER_SECONDARY_PRESETS[idx - kBHEnd];
                auto kind = static_cast<Simulation::MergerSecondaryKind>(sec.kind);
                sim.startMerger(sec.massSolar, windowHeight, kind);
                ui.mergerMenu.open = false;
                ui.notification = std::string("Merger with ") + sec.name + " initiated!";
                ui.notificationTimer = 150;
            } else {
                double massSolar = BH2D_PRESETS[idx].massSolar;
                sim.startMerger(massSolar, windowHeight);
                ui.mergerMenu.open = false;
                ui.notification = std::string("Merger with ") +
                                  BH2D_PRESETS[idx].name + " initiated!";
                ui.notificationTimer = 150;
            }
        }
        return; // swallow all keypresses while menu is open
    }

    // ── Custom-body creator menu navigation ─────────────────────────────────
    if (ui.customBodyMenu.open) {
        auto& m = ui.customBodyMenu;
        const int nPresets = (int)sim.customPresets.size();
        const int nFields  = UIState::kCustomBodyInputCount;     // 4
        const int nRows    = nFields + nPresets;
        if (code == sf::Keyboard::Key::Escape) {
            m.open = false;
        } else if (code == sf::Keyboard::Key::Tab || code == sf::Keyboard::Key::Down) {
            m.focusedField = (m.focusedField + 1) % nRows;
        } else if (code == sf::Keyboard::Key::Up) {
            m.focusedField = (m.focusedField + nRows - 1) % nRows;
        } else if (m.focusedField == 1 &&
                   (code == sf::Keyboard::Key::Left || code == sf::Keyboard::Key::Right)) {
            int delta = (code == sf::Keyboard::Key::Right) ? 1 : -1;
            m.typeIdx = (m.typeIdx + delta + UIState::kCustomBodyTypeCount)
                        % UIState::kCustomBodyTypeCount;
        } else if (code == sf::Keyboard::Key::Delete && m.focusedField >= nFields) {
            int presetIdx = m.focusedField - nFields;
            std::string nm = sim.customPresets[presetIdx].name;
            sim.removeCustomPreset(presetIdx);
            int newRowCount = nFields + (int)sim.customPresets.size();
            if (m.focusedField >= newRowCount) m.focusedField = std::max(0, newRowCount - 1);
            ui.notification = "Deleted preset: " + nm;
            ui.notificationTimer = 120;
        } else if (code == sf::Keyboard::Key::Enter) {
            if (m.focusedField >= nFields) {
                const auto& p = sim.customPresets[m.focusedField - nFields];
                sim.addCustomBody(p.type, p.semiMajorM, p.ecc, p.name.c_str());
                m.open = false;
                ui.notification = "Spawned preset: " + p.name;
                ui.notificationTimer = 150;
            } else {
                double sm  = 20.0;
                double ecc = 0.3;
                try { sm  = std::stod(m.smInput);  } catch (...) {}
                try { ecc = std::stod(m.eccInput); } catch (...) {}
                if (sm < 2.5 || ecc < 0.0 || ecc >= 1.0) {
                    ui.notification = "Bad input: need sm >= 2.5 M and 0 <= e < 1";
                    ui.notificationTimer = 150;
                } else {
                    GalaxyBodyType t = UIState::customBodyTypes()[m.typeIdx];
                    std::string nm = sim.appendCustomPreset(m.nameInput, t, sm, ecc);
                    sim.addCustomBody(t, sm, ecc, nm.c_str());
                    m.open = false;
                    m.nameInput.clear();
                    ui.notification = "Spawned + saved: " + nm;
                    ui.notificationTimer = 150;
                }
            }
        }
        return;
    }

    // ── Simulation ────────────────────────────────────────────────────────────
    if (code == cfg.pause) {
        ui.paused = !ui.paused;

    } else if (code == cfg.togglePreset) {
        ui.presetActive = !ui.presetActive;
        sim.clearTransientEvents();
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
            sim.clearTransientEvents();
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
            sim.clearTransientEvents();
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

    } else if (code == cfg.resetView) {
        // Snap the camera back to a sensible default for the current scenario
        // without disturbing physics state (orbits, merger, etc.).
        if (ui.presetActive) {
            sim.params.pixelsPerM = UIState::presetHorizonPixelsTarget / (2.0 * sim.bh.metric.M);
        } else {
            sim.params.pixelsPerM = ui.defaultPixelsPerM;
        }
        ui.notification = "View reset";
        ui.notificationTimer = 90;

    } else if (code == cfg.toggleLightMode) {
        ui.lightMode = !ui.lightMode;
        ui.notification = ui.lightMode ? "Light mode" : "Dark mode";
        ui.notificationTimer = 90;

    } else if (code == cfg.learnMore) {
        // Open the active preset's Wikipedia / reference URL in the user's
        // default browser. Only built-in `BH2D_PRESETS` carry URLs; custom
        // presets and the merger menu are intentionally URL-less.
        if (ui.presetActive && ui.presetIdx >= 0 && ui.presetIdx < NUM_BH2D_PRESETS) {
            const char* url = BH2D_PRESETS[ui.presetIdx].learnMoreUrl;
            if (url && *url) {
                platformOpenUrl(url);
                ui.notification = std::string("Opening ") + BH2D_PRESETS[ui.presetIdx].name + " reference...";
            } else {
                ui.notification = "No reference URL for this preset";
            }
        } else {
            ui.notification = "Activate a preset (T) to open its reference";
        }
        ui.notificationTimer = 120;

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
        bool wasActive = sim.anyBodyWithScenarioTag(ResearchScenario::ISCOTest);
        sim.startISCOTest();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = wasActive
            ? "ISCO test bodies removed"
            : "ISCO validation: 5M/6M/7M test added";
        ui.notificationTimer = 120;

    } else if (code == cfg.testPhoton) {
        sim.startPhotonSphereTest();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = std::string("Photon sweep: ") + sim.photonSphereSweepName();
        ui.notificationTimer = 120;

    } else if (code == cfg.testInfall) {
        bool wasActive = sim.anyBodyWithScenarioTag(ResearchScenario::RadialInfall);
        sim.startRadialInfall();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = wasActive
            ? "Radial infall body removed"
            : "Radial infall from 20M added";
        ui.notificationTimer = 120;

    } else if (code == cfg.testTidal) {
        bool wasActive = sim.anyBodyWithScenarioTag(ResearchScenario::TidalDisruption);
        sim.startTidalDisruption();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = wasActive
            ? "Tidal test body removed"
            : "Tidal disruption demo added";
        ui.notificationTimer = 120;

    } else if (code == cfg.testPulsar) {
        bool wasActive = (sim.findPulsarIdx() >= 0);
        sim.togglePulsarOrbital();
        sim.rebuildPhotons(windowHeight);
        ui.notification      = wasActive
            ? "Pulsar removed"
            : "Pulsar added, enable preset (T) for physical GW values";
        ui.notificationTimer = 150;

    } else if (code == cfg.mergerMenu) {
        ui.mergerMenu.open            = !ui.mergerMenu.open;
        ui.mergerMenu.inputtingCustom = false;
        ui.mergerMenu.selectedIdx     = 0;

    } else if (code == cfg.customBodyMenu) {
        ui.customBodyMenu.open         = !ui.customBodyMenu.open;
        ui.customBodyMenu.focusedField = 0;

    } else if (code == cfg.cycleMergerSpeed) {
        int next = (static_cast<int>(sim.merger.timeScale) + 1) % 3;
        sim.merger.timeScale = static_cast<Simulation::MergerState::TimeScale>(next);
        ui.notification = std::string("Merger pacing: ") +
                          Simulation::MergerState::timeScaleName(sim.merger.timeScale);
        ui.notificationTimer = 120;
    }
}
