/*--------- BlackHole2D.cpp ---------*/
// Author: Liam N.
// Purpose: The main file for the 2D, more research & data-collection oriented half of the Aetherion Astrophysics Suite
// Status: In Beta development.


/*
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⢠⢀⡐⢄⢢⡐⢢⢁⠂⠄⠠⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀⠀⠀⠀⡄⣌⠰⣘⣆⢧⡜⣮⣱⣎⠷⣌⡞⣌⡒⠤⣈⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠒⠊⠀⠀⠀⠀⢀⠢⠱⡜⣞⣳⠝⣘⣭⣼⣾⣷⣶⣶⣮⣬⣥⣙⠲⢡⢂⠡⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠃⠀⠀⠀⠀⠀⠀⢀⠢⣑⢣⠝⣪⣵⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣯⣻⢦⣍⠢⢅⢂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠀⠀⢆⡱⠌⣡⢞⣵⣿⣿⣿⠿⠛⠛⠉⠉⠛⠛⠿⢷⣽⣻⣦⣎⢳⣌⠆⡱⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠂⠠⠌⢢⢃⡾⣱⣿⢿⡾⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⢻⣏⠻⣷⣬⡳⣤⡂⠜⢠⡀⣀⠀⠀⡀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⢀⠂⣌⢃⡾⢡⣿⢣⡏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⡇⡊⣿⣿⣾⣽⣛⠶⣶⣬⣭⣥⣙⣚⢷⣶⠦⡤⢀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⢁⠂⠰⡌⡼⠡⣼⢃⡿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢻⣿⣿⣿⣿⣿⣿⣿⣿⣾⡿⠿⣛⣯⡴⢏⠳⠁⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠠⠑⡌⠀⣉⣾⣩⣼⣿⣾⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣀⣠⣤⣤⣿⣿⣿⣿⡿⢛⣛⣯⣭⠶⣞⠻⣉⠒⠀⠂⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢀⣀⡶⢝⣢⣾⣿⣼⣿⣿⣿⣿⣿⣀⣼⣀⣀⣀⣤⣴⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⣿⠿⡛⠏⠍⠂⠁⢠⠁⠀⠀⠀⠀⠀⠀⠀
⠀⠠⢀⢥⣰⣾⣿⣯⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⣽⠟⣿⠐⠨⠑⡀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⡐⢢⣟⣾⣿⣿⣟⣛⣿⣿⣿⣿⢿⣝⠻⠿⢿⣯⣛⢿⣿⣿⣿⡛⠻⠿⣛⠻⠛⡛⠩⢁⣴⡾⢃⣾⠇⢀⠡⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠈⠁⠊⠙⠉⠩⠌⠉⠢⠉⠐⠈⠂⠈⠁⠉⠂⠐⠉⣻⣷⣭⠛⠿⣶⣦⣤⣤⣴⣴⡾⠟⣫⣾⣿⡏⠀⠂⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⢻⢿⢶⣤⣬⣉⣉⣭⣤⣴⣿⣿⡿⠃⠄⡈⠁⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠘⢊⠳⠭⡽⣿⠿⠿⠟⠛⠉⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠁⠈⠐⠀⠘⠀⠈⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

I know it has absolutely nothing to do with the program, but I just want to leave this beautiful art here <3 */
/*--------- Headers ---------*/
#include <SFML/Graphics.hpp>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "2D-physics/units.hpp"
#include "2D-simulation/simulation.hpp"
#include "2D-rendering/camera.hpp"
#include "2D-rendering/renderer.hpp"
#include "2D-visualization/ray_visualizer.hpp"
#include "2D-visualization/orbit_visualizer.hpp"
#include "2D-ui/controls.hpp"
#include "2D-utils/presets_2d.hpp"
#include <cstring>

/*--------- HUD text builder ---------*/
// this function is doing a lot of string formatting inline. I know it should probably be a separate
// class or at least split into smaller functions but every time I start refactoring it I end up
// breaking something subtle in the formatting, so it stays as one big function. for now.
static std::string buildHUD(const Simulation& sim, const UIState& ui) {
    auto pretty = [](double v) -> std::string {
        std::ostringstream os;
        if (std::fabs(v) >= 1e6 || (std::fabs(v) > 0 && std::fabs(v) < 1e-3))
            os << std::scientific << std::setprecision(3) << v;
        else
            os << std::fixed << std::setprecision(4) << v;
        return os.str();
    };

    std::ostringstream ss;
    if (ui.presetActive) {
        double M = sim.bh.metric.M;
        double r_h_m = sim.bh.metric.horizon();
        double r_h_km = r_h_m / 1000.0;
        double r_h_AU = r_h_m / 1.495978707e11;
        double r_h_light_days = r_h_m / (units::c_SI * 86400.0);

        const auto& preset = BH2D_PRESETS[ui.presetIdx];
        const double currentMsun = sim.currentMassSolar();
        ss << "Preset: " << preset.name
           << " (" << pretty(preset.massSolar) << " Msun";
        if (sim.merger.completed)
            ss << " → " << pretty(currentMsun) << " Msun post-merger";
        ss << ")\n";
        ss << preset.description << "\n";
        ss << "M (geom) = " << pretty(M) << " m\n";
        ss << "Event horizon r_h = " << pretty(r_h_m) << " m ("
           << pretty(r_h_km) << " km, " << pretty(r_h_AU) << " AU, "
           << pretty(r_h_light_days) << " light-days)\n";
        ss << "Display scale: horizon ~ "
           << static_cast<int>(2.0 * M * sim.params.pixelsPerM + 0.5)
           << " px (pixelsPerMeter=" << std::scientific << std::setprecision(3)
           << sim.params.pixelsPerM << ")\n";

        // Galaxy system body summary
        if (sim.galaxySystemActive) {
            int nStars = 0, nGas = 0, nCluster = 0, nDwarf = 0;
            int nNS = 0, nWD = 0, nComp = 0;
            for (const auto& b : sim.bodies) {
                if (!b.isGalaxyBody) continue;
                switch (b.bodyType) {
                    case GalaxyBodyType::Star:           nStars++;   break;
                    case GalaxyBodyType::GasCloud:       nGas++;     break;
                    case GalaxyBodyType::StellarCluster: nCluster++; break;
                    case GalaxyBodyType::DwarfGalaxy:    nDwarf++;   break;
                    case GalaxyBodyType::NeutronStar:    nNS++;      break;
                    case GalaxyBodyType::WhiteDwarf:     nWD++;      break;
                    case GalaxyBodyType::CompanionStar:  nComp++;    break;
                }
            }
            ss << "Galaxy system: ";
            if (nStars)   ss << nStars << " stars ";
            if (nGas)     ss << nGas << " gas clouds ";
            if (nCluster) ss << nCluster << " clusters ";
            if (nDwarf)   ss << nDwarf << " dwarfs ";
            if (nNS)      ss << nNS << " neutron stars ";
            if (nWD)      ss << nWD << " white dwarfs ";
            if (nComp)    ss << nComp << " companions ";
            ss << "\n";

            if (ui.showInfluenceZones && preset.isGalacticCenter) {
                ss << "Zones: tidal=" << pretty(preset.zones.tidalDisruptionM * M)
                   << "m  Bondi=" << pretty(preset.zones.bondiRadiusM * M)
                   << "m  SoI=" << pretty(preset.zones.sphereOfInfluenceM * M) << "m\n";
            }
        }

        if (!sim.bodies.empty()) {
            const auto& body = sim.bodies[0];
            ss << "r = " << pretty(body.radius()) << "  v = "
               << pretty(body.measurement.orbitalVelocity) << "c  1+z = "
               << pretty(body.measurement.redshift) << "\n";
        }
        ss << "Presets: [ / ] to switch | T toggle preset | Up/Down zoom\n";
    } else {
        double M = sim.bh.metric.M;
        // Pulsar orbital scenario: show GW summary instead of generic orbital state
        if (sim.activeScenario == ResearchScenario::PulsarOrbital) {
            const auto& pd = sim.pulsarData;
            auto gw = [](double v) -> std::string {
                std::ostringstream os;
                os << std::scientific << std::setprecision(2) << v;
                return os.str();
            };
            ss << "Pulsar orbital sim  [" << std::fixed << std::setprecision(2)
               << sim.pulsarState.massSolar << " Msun NS @ "
               << std::setprecision(1) << pd.semiMajorAxis_M << "M]"
               << (pd.disrupted ? "  DISRUPTED" : "") << "\n";
            ss << "Shapiro: " << gw(pd.shapiroDelay_us) << " us  "
               << "nu/nu0: " << gw(pd.gravRedshift) << "\n";
            ss << "h(1kpc): " << gw(pd.gwStrain)
               << "  f_GW: " << gw(pd.gwFreq_Hz) << " Hz\n";
            ss << "r/r_Roche: " << pretty(pd.rocheMargin)
               << "  e: " << pretty(pd.eccentricity) << "\n";
            ss << "D: data panel  5: restart scenario\n";
        } else {
            ss << "M = " << pretty(M) << "  (r_h=" << pretty(sim.bh.metric.horizon())
               << ", ISCO=" << pretty(sim.bh.metric.isco()) << ")\n";
        if (!sim.bodies.empty()) {
            const auto& body = sim.bodies[0];
            ss << "E = " << pretty(body.E)
               << "  L = " << pretty(body.L)
               << "  r = " << pretty(body.radius())
               << "  v = " << pretty(body.measurement.orbitalVelocity) << "c\n";
            ss << "1+z = " << pretty(body.measurement.redshift)
               << "  gamma = " << pretty(body.measurement.timeDilation)
               << "  (a~" << pretty(body.nominalA)
               << " e~" << pretty(body.nominalEcc) << ")\n";
        }
        ss << "Photon sphere (r=3M) shown: " << (ui.showPhotonSphere ? "YES" : "NO")
           << "    Rays: " << (ui.showRays ? "ON" : "OFF") << "\n";
        }
    }
    if (ui.timeScale != 1.0) {
        std::ostringstream ts;
        if (ui.timeScale >= 1.0) ts << (int)ui.timeScale;
        else ts << std::fixed << std::setprecision(2) << ui.timeScale;
        ss << "Speed: " << ts.str() << "x  |  ";
    }
    ss << "Press / for controls";
    return ss.str();
}

/*--------- Entry point ---------*/
int main(int argc, char* argv[]) {
    const unsigned int WIN_W = 1200, WIN_H = 800;
    sf::RenderWindow window(sf::VideoMode(sf::Vector2u(WIN_W, WIN_H)),
                            "Schwarzschild BH - 2D demo");
    window.setFramerateLimit(60);

    Simulation sim;
    sim.rebuildPhotons(WIN_H);
    sim.loadCustomPresets();

    Camera   camera{ sim.params.pixelsPerM, (float)WIN_W, (float)WIN_H };
    Renderer renderer(window);
    UIState  ui;

    // CLI: --preset <id>  (ton618, sgra, 3c273, j0529, m87, cygnusx1, gw150914, intermediate, primordial,
    //                       gaiabh1, gaiabh2, gaiabh3, v404cyg, a062000, groj165540, ngc1277, oj287)
    // useful for screenshots and demo, lets you launch directly into a specific black hole without
    // clicking through the menu. the string-to-index mapping below is clunky but there aren't that many
    // presets so I'm not going to build a proper lookup table for it
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--preset") == 0) {
            std::string id = argv[i + 1];
            int idx = -1;
            if      (id == "ton618")       idx = 0;
            else if (id == "sgra")         idx = 1;
            else if (id == "3c273")        idx = 2;
            else if (id == "j0529")        idx = 3;
            else if (id == "m87")          idx = 4;
            else if (id == "cygnusx1")     idx = 5;
            else if (id == "gw150914")     idx = 6;
            else if (id == "intermediate") idx = 7;
            else if (id == "primordial")   idx = 8;
            else if (id == "gaiabh1")      idx = 9;
            else if (id == "gaiabh2")      idx = 10;
            else if (id == "gaiabh3")      idx = 11;
            else if (id == "v404cyg")      idx = 12;
            else if (id == "a062000")      idx = 13;
            else if (id == "groj165540")   idx = 14;
            else if (id == "ngc1277")      idx = 15;
            else if (id == "oj287")        idx = 16;
            if (idx >= 0 && idx < NUM_BH2D_PRESETS) {
                ui.presetActive = true;
                ui.presetIdx    = idx;
                sim.bh.metric.M   = units::solarMassToGeomMeters(BH2D_PRESETS[idx].massSolar);
                sim.params.pixelsPerM = UIState::presetHorizonPixelsTarget / (2.0 * sim.bh.metric.M);
                sim.tidalRadiusM  = BH2D_PRESETS[idx].zones.tidalDisruptionM;
                sim.reinitBodies();
                if (ui.showGalaxySystem && BH2D_PRESETS[idx].isGalacticCenter)
                    sim.spawnGalaxySystem(idx);
                sim.rebuildPhotons(WIN_H);
            }
            break;
        }
    }

    // Set the initial view to match window size
    // this has to be done after the window is created and before the main loop.
    // if you forget it, SFML renders into a 1x1 pixel view or something equally absurd.
    // I forgot it once. it was a confusing 20 minutes.
    renderer.resetView(camera.viewWidth, camera.viewHeight);

    sf::Clock clock;
    sf::Clock bgClock;   // separate clock so background drift is independent of paused state
    std::vector<sf::Vertex> rayVertScratch; // reusable scratch buffer, before this existed we were allocating a new vector every frame for every ray.
                                             // at 120 rays × 60fps that was generating ~7200 heap allocations per second just for ray drawing.
                                             // embarrassing in retrospect but at least it's fixed now

    // HUD text caching to avoid expensive string formatting every frame
    std::string cachedHUD;
    int cachedPresetIdx = -1;
    bool cachedPresetActive = false;
    bool cachedGalaxySystemActive = false;
    size_t cachedBodyCount = 0;
    ResearchScenario cachedScenario = ResearchScenario::None;
    int cachedNotificationTimer = -1;
    double hudRefreshAccum = 0.0;
    constexpr double HUD_REFRESH_SECONDS = 0.2;

    auto hudNeedsRebuild = [&]() {
        return cachedPresetIdx != ui.presetIdx
            || cachedPresetActive != ui.presetActive
            || cachedGalaxySystemActive != sim.galaxySystemActive
            || cachedBodyCount != sim.bodies.size()
            || cachedScenario != sim.activeScenario
            || cachedNotificationTimer != ui.notificationTimer;
    };

    while (window.isOpen()) {
        /*--------- Input ---------*/
        while (auto evOpt = window.pollEvent()) {
            if (evOpt->is<sf::Event::Closed>()) { window.close(); continue; }

            // Handle window resize / fullscreen: update view + camera
            // SFML 3.x changed how resize events work vs 2.x, this tripped me up for a while.
            // the view needs to be explicitly reset or everything stays at the old dimensions
            // and you get stretched/clipped rendering. fun bug to track down at midnight.
            if (auto* resized = evOpt->getIf<sf::Event::Resized>()) {
                float w = (float)resized->size.x;
                float h = (float)resized->size.y;
                camera.updateSize(w, h);
                renderer.resetView(w, h);
            }

            handleInput(*evOpt, sim, ui, (unsigned int)camera.viewHeight);
        }

        /*--------- Update ---------*/
        double dt = clock.restart().asSeconds();
        if (!ui.paused)
            sim.update(dt * ui.timeScale);
        hudRefreshAccum += dt;

        // Auto-zoom during merger inspiral so both BHs stay framed even at
        // extreme mass ratios. Skip when finished or during the post-merger
        // flash so the remnant shot isn't suddenly rescaled.
        if (sim.merger.active && sim.merger.flashTimer <= 0.0 && sim.merger.r_M > 0.0) {
            const double sepM = sim.merger.r_M * sim.bh.metric.M;
            const double halfExtentM = sepM * 1.35 + 4.0 * sim.bh.metric.M;
            const double targetHalfPx = std::min(camera.viewWidth, camera.viewHeight) * 0.40;
            double targetPpm = targetHalfPx / std::max(halfExtentM, 1e-9);
            // Don't let auto-zoom blow the horizon up past ~250 px wide as r→0.
            const double ppmHorizonCap = 125.0 / std::max(sim.bh.metric.M, 1e-9);
            targetPpm = std::min(targetPpm, ppmHorizonCap);
            const double k = std::clamp(dt * 2.5, 0.0, 0.25);
            sim.params.pixelsPerM = sim.params.pixelsPerM + (targetPpm - sim.params.pixelsPerM) * k;
        }

        camera.pixelsPerM = sim.params.pixelsPerM;
        sf::Vector2f center = camera.center();
        double M = sim.bh.metric.M;

        // ── Barycentric binary mode (Gaia BH1 / BH2 / BH3) ────────────────────
        // For these systems the companion mass is significant enough that the BH
        // visibly wobbles around the shared center of mass, exactly the signal
        // Gaia measured astrometrically.  The Schwarzschild physics is still
        // integrated in BH-centred coords, but at render time we:
        //   • keep the *barycenter* at screen centre
        //   • draw the BH offset by  r_BH = -r_comp × q  (q = M_comp/(M_BH+M_comp))
        //   • draw the companion at   r_comp_bary = r_comp × (1-q)  (= r_comp × M_BH/(M_BH+M_comp))
        //   • show a "Center of Mass" marker at screen centre
        bool         barycentricMode = false;
        float        baryScale       = 1.0f;   // companion world-coord scale factor
        sf::Vector2f bhCenter        = center; // screen position of the BH (default: screen centre)

        if (ui.presetActive && ui.showGalaxySystem
            && BH2D_PRESETS[ui.presetIdx].isBinaryWithBarycenter
            && !sim.bodies.empty())
        {
            const auto&  preset  = BH2D_PRESETS[ui.presetIdx];
            double M_comp = preset.companionMassSolar;
            double M_bh   = preset.massSolar;
            double q      = M_comp / (M_bh + M_comp);   // fraction of separation = BH offset
            baryScale     = (float)(1.0 - q);            // companion scale: M_BH / (M_BH + M_comp)

            // Current companion position in BH-centred world coords
            double cx = sim.bodies[0].worldX();
            double cy = sim.bodies[0].worldY();

            // BH sits at -q × (companion vector) in barycenter-centred coords
            bhCenter = camera.worldToScreen((float)(-cx * q), (float)(-cy * q));
            barycentricMode = true;
        }

        // Merger barycentric mode: while the secondary is spiralling in, the
        // primary BH visibly orbits the common centre of mass too. We treat the
        // sim origin as the barycentre and offset the primary draw position by
        //   r_primary = -q * sep * [cos(phi), sin(phi)]
        // All primary visualisations (horizon, disk, photon sphere, bodies)
        // pick up this offset because they’re drawn relative to bhCenter.
        if (sim.merger.active && sim.merger.flashTimer <= 0.0) {
            const double q   = sim.merger.massRatioQ;
            // Honesty gate: the primary doesn’t actually move in the body
            // integrator (see TODO(physics-honesty) in simulation.hpp). Only
            // wobble the renderer when the lie is visually justified — i.e.
            // when q is high enough that a stationary primary would look
            // obviously wrong. Below q=0.05 the offset is sub-pixel anyway
            // and we get to be honest for free.
            if (q >= 0.05) {
                const double sep = sim.merger.r_M * M;
                const double bx  = q * sep * std::cos(sim.merger.phi);
                const double by  = q * sep * std::sin(sim.merger.phi);
                bhCenter = camera.worldToScreen((float)(-bx), (float)(-by));
            }
        }

        // Rebuild HUD text only when relevant simulation state changes
        if (hudRefreshAccum >= HUD_REFRESH_SECONDS || hudNeedsRebuild()) {
            cachedHUD = buildHUD(sim, ui);
            cachedPresetIdx = ui.presetIdx;
            cachedPresetActive = ui.presetActive;
            cachedGalaxySystemActive = sim.galaxySystemActive;
            cachedBodyCount = sim.bodies.size();
            cachedScenario = sim.activeScenario;
            cachedNotificationTimer = ui.notificationTimer;
            hudRefreshAccum = 0.0;
        }

        /*--------- Draw ---------*/
        renderer.setLightMode(ui.lightMode);
        renderer.beginFrame();
        renderer.drawStarfield(bgClock.getElapsedTime().asSeconds());

        // Influence zones (behind everything except starfield)
        if (ui.presetActive && ui.showInfluenceZones && ui.showGalaxySystem) {
            const auto& preset = BH2D_PRESETS[ui.presetIdx];
            if (preset.isGalacticCenter) {
                float ppm = (float)camera.pixelsPerM;
                float soiPx   = (float)(preset.zones.sphereOfInfluenceM * M * ppm);
                float bondiPx = (float)(preset.zones.bondiRadiusM * M * ppm);
                float tidalPx = (float)(preset.zones.tidalDisruptionM * M * ppm);

                // Influence zones are centred on the BH, not the barycenter
                renderer.drawInfluenceZone(bhCenter, soiPx,
                    sf::Color(100, 255, 100, 80), "Sphere of Influence");
                renderer.drawInfluenceZone(bhCenter, bondiPx,
                    sf::Color(255, 200, 50, 90), "Bondi Radius");
                renderer.drawInfluenceZone(bhCenter, tidalPx,
                    sf::Color(255, 80, 80, 120), "Tidal Disruption");

                if (preset.zones.hasJetCones) { // draws relativistic jets if the preset has them (most don't)
                    float jetLen = soiPx * 1.5f;
                    renderer.drawJetCones(bhCenter, jetLen);
                }
            }
        }

        renderer.drawHorizon(bhCenter, (float)(sim.bh.metric.horizon() * camera.pixelsPerM));

        if (ui.presetActive) { //  && BH2D_PRESETS[ui.presetIdx].hasAccretionDisk maybe?
            // the accretion disk is just a glowing ring drawn on screen. it's not physically simulated.
            // please do not ask me when I'm adding MHD accretion disk dynamics. the answer is never.
            float diskRadius = (float)(6.0 * M * camera.pixelsPerM);
            renderer.drawAccretionDisk(bhCenter, diskRadius);
        }

        if (ui.showPhotonSphere)
            renderer.drawPhotonSphere(bhCenter, (float)(sim.bh.metric.photonSphere() * camera.pixelsPerM));

        // Rays, visualization layer computes colors
        if (ui.showRays) {
            for (const auto& photon : sim.photons) {
                if (photon.captured) {
                    auto [v0, v1] = RayVisualizer::capturedRayLine(photon.impactParameter, camera);
                    renderer.drawCapturedRay(v0, v1);
                } else {
                    RayVisualizer::colorByRedshift(photon, sim.bh.metric, camera, rayVertScratch);
                    renderer.drawRayPath(rayVertScratch);
                }
            }
        }

        // Kerr equatorial overlay (key J) — thin amber lines for Kerr geodesics
        if (sim.kerrOverlayEnabled && !sim.kerrRays.empty()) {
            static constexpr sf::Color kKerrCol(255, 180, 50, 130);
            for (const auto& kr : sim.kerrRays) {
                const size_t nPts = kr.verts.size() / 2;
                if (nPts < 2) continue;
                sf::VertexArray va(sf::PrimitiveType::LineStrip, nPts);
                for (size_t i = 0; i < nPts; ++i) {
                    va[i].position = camera.worldToScreen(kr.verts[i*2], kr.verts[i*2+1]);
                    va[i].color    = kKerrCol;
                }
                renderer.drawOrbitPath(va);
            }
        }

        // Disk-emitter photons (key O) — thermal null geodesics from inner disk
        if (sim.diskEmitterEnabled && !sim.emittedPhotons.empty()) {
            for (const auto& photon : sim.emittedPhotons) {
                if (photon.captured) {
                    auto [v0, v1] = RayVisualizer::capturedRayLine(photon.impactParameter, camera);
                    renderer.drawCapturedRay(v0, v1);
                } else {
                    RayVisualizer::colorByRedshift(photon, sim.bh.metric, camera, rayVertScratch);
                    renderer.drawRayPath(rayVertScratch);
                }
            }
        }

        // Caustic highlights (where rays converge strongly)
        if (ui.showCaustics && !sim.lensingData.causticPoints.empty()) {
            for (const auto& [avgB, gradient] : sim.lensingData.causticPoints) {
                // Mark caustic regions on both sides of the BH
                sf::Vector2f pos1 = camera.worldToScreen(0.0f, (float)avgB);
                sf::Vector2f pos2 = camera.worldToScreen(0.0f, (float)-avgB);
                renderer.drawCausticMarker(pos1, gradient);
                renderer.drawCausticMarker(pos2, gradient);
            }
        }

        // Time dilation heatmap
        if (ui.showTimeDilationMap)
            renderer.drawTimeDilationMap(center, M, camera.pixelsPerM);

        // Barycenter marker, drawn before bodies so it sits behind them
        if (barycentricMode)
            renderer.drawBarycenterMarker(center);

        // Bodies + orbit paths
        for (size_t bodyIdx = 0; bodyIdx < sim.bodies.size(); ++bodyIdx) {
            const auto& body = sim.bodies[bodyIdx];
            // In barycentric mode the companion's world coords are scaled by baryScale
            // (= M_BH / (M_BH + M_comp)) so both objects orbit screen centre.
            float scale = (barycentricMode && body.isGalaxyBody) ? baryScale : 1.0f;
            auto bodyVis = OrbitVisualizer::computeBodyVisual(body, sim.bh.metric, camera, scale);
            if (body.isSecondaryBH) {
                // Mini SMBH visual (horizon + tinted disk). Sized cosmetically
                // so the binary nature is unmistakable even at preset zoom.
                const float secHorizonPx = 6.0f;
                const float secDiskPx    = 16.0f;
                renderer.drawMergerBH(bodyVis.screenPos, secHorizonPx,
                                      secDiskPx, 0.0f, body.label.c_str());
            } else {
                renderer.drawBody(bodyVis);
            }

            // Pulsar-specific visual decoration (jets, magnetic field lines, LC ring, flux tubes)
            if (sim.activeScenario == ResearchScenario::PulsarOrbital && body.isPulsar
                && !body.captured) {
                float lc_px = 0.0f;
                if (sim.pulsarState.lightCylRadius_m > 0.0)
                    lc_px = (float)(sim.pulsarState.lightCylRadius_m * camera.pixelsPerM);
                renderer.drawPulsarBody(
                    bodyVis.screenPos,
                    (float)sim.pulsarState.spinPhase,
                    (float)sim.pulsarState.precPhase,
                    lc_px,
                    sim.pulsarData.inLightCylinder,
                    bhCenter,                        // BH screen position
                    sim.pulsarData.magPower_ergs
                );
            }

            auto orbitPath = OrbitVisualizer::computeOrbitPath(body, sim.bh.metric, camera, scale);
            renderer.drawOrbitPath(orbitPath);

            // Label galaxy system bodies
            if (body.isGalaxyBody && ui.showGalaxySystem && !body.label.empty() && !body.isSecondaryBH) {
                sf::Color labelColor;
                switch (body.bodyType) {
                    case GalaxyBodyType::Star:           labelColor = sf::Color(255, 240, 200); break;
                    case GalaxyBodyType::GasCloud:       labelColor = sf::Color(255, 150, 100); break;
                    case GalaxyBodyType::StellarCluster:  labelColor = sf::Color(200, 220, 255); break;
                    case GalaxyBodyType::DwarfGalaxy:    labelColor = sf::Color(255, 200, 255); break;
                    case GalaxyBodyType::NeutronStar:    labelColor = sf::Color(160, 200, 255); break;
                    case GalaxyBodyType::WhiteDwarf:     labelColor = sf::Color(200, 215, 255); break;
                    case GalaxyBodyType::CompanionStar:  labelColor = sf::Color(255, 220, 160); break;
                    /*
                    case GalaxyBodyType::Blanet:         labelColor = sf::Color(100, 100, 100); break; // WILL COME BACK TO THESE AT ANOTHER TIME
                    case GalaxyBodyType::BrownDwarf:     labelColor = sf::Color(); break;
                    case GalaxyBodyType::G-Body:         labelColor = sf::Color(180, 255, 180); break;
                    case GalaxyBodyType::SecondaryBH:    labelColor = sf::Color(); break;
                    */ 
                }
                renderer.drawBodyLabel(bodyVis.screenPos, body.label, labelColor);
            }

            // Label research scenario bodies, makes things a lot less confusing
            if (!body.isGalaxyBody && !body.label.empty() &&
                sim.activeScenario != ResearchScenario::None) {
                renderer.drawBodyLabel(bodyVis.screenPos, body.label,
                    sf::Color(180, 220, 255));
            }
        }

        // ---- Tidal disruption event overlay ----
        if (sim.tidalEvent.active) {
            sf::Vector2f flashPos = {
                bhCenter.x + static_cast<float>(sim.tidalEvent.eventX * camera.pixelsPerM),
                bhCenter.y - static_cast<float>(sim.tidalEvent.eventY * camera.pixelsPerM)
            };
            float flashAlpha = (sim.tidalEvent.flashTimer > 0.0)
                ? static_cast<float>(sim.tidalEvent.flashTimer /
                                     Simulation::TidalEvent::FLASH_DURATION)
                : 0.0f;
            bool showLabel = flashAlpha > 0.0f || !sim.tidalEvent.particles.empty();
            std::vector<Renderer::TidalParticleVis> pvis;
            pvis.reserve(sim.tidalEvent.particles.size());
            for (const auto& p : sim.tidalEvent.particles) {
                float lifeF = static_cast<float>(p.lifetime / p.maxLifetime);
                pvis.push_back({
                    bhCenter.x + static_cast<float>(p.x * camera.pixelsPerM),
                    bhCenter.y - static_cast<float>(p.y * camera.pixelsPerM),
                    p.size * lifeF,
                    lifeF,
                    p.isFallback
                });
            }
            renderer.drawTidalEvent(flashPos, flashAlpha, pvis, showLabel,
                sim.tidalEvent.kind == Simulation::TidalEventKind::MergerShockwave ? "MERGER REMNANT"
              : sim.tidalEvent.kind == Simulation::TidalEventKind::XrayBurst       ? "X-RAY BURST"
              : "TIDAL DISRUPTION EVENT");
        }

        renderer.drawHUD(cachedHUD);

        // Controls panel (left side, toggleable with / or ? symbol)
        if (ui.showControlsPanel)
            renderer.drawControlsPanel(camera.viewHeight);

        // Data panel (right side)
        if (ui.showDataPanel)
            renderer.drawDataPanel(sim.formatDataPanel(),
                                   camera.viewWidth, camera.viewHeight);

        // Merger: draw death-spiral trails and incoming BH during inspiral.
        // Trails are stored in barycentre-centred world coords, so the draw
        // origin is the screen centre (= world origin = barycentre), not bhCenter.
        if (sim.merger.active && sim.merger.flashTimer <= 0.0) {
            const double q   = sim.merger.massRatioQ;
            const double sep = sim.merger.r_M * sim.bh.metric.M;

            // Secondary BH screen position: at +(1-q)*sep from the barycentre.
            const double r2   = sep * (1.0 - q);
            sf::Vector2f secPos = camera.worldToScreen(
                (float)( r2 * std::cos(sim.merger.phi)),
                (float)( r2 * std::sin(sim.merger.phi)));

            // Death-spiral trails. Both are drawn from the barycentre (= screen
            // centre), so primary + secondary trace symmetrically opposite arcs
            // around the centre of mass.
            renderer.drawMergerTrail(sim.merger.trail2, center,
                                     camera.pixelsPerM,
                                     sf::Color(180, 100, 255));   // secondary (purple)
            if (q > 0.01) {
                renderer.drawMergerTrail(sim.merger.trail1, center,
                                         camera.pixelsPerM,
                                         sf::Color(100, 220, 255)); // primary (cyan)
            }

            // Incoming secondary — black hole vs. compact / stellar object.
            // For BH secondaries: full horizon + disk (+ jets for AGN-class).
            // For NS/Pulsar/Star/WD: a small bright dot with kind colour.
            if (sim.merger.secondaryKind == Simulation::MergerSecondaryKind::BlackHole) {
                // Sized from its real mass so a TON 618-class secondary looks
                // visibly larger than a stellar BH. Disk radius mirrors the
                // primary’s 6M convention; jets are added for AGN-class
                // secondaries (>= 1e6 Msun) so visually it reads as an SMBH
                // falling in.
                double m2_geom = units::solarMassToGeomMeters(sim.merger.massSolar);
                float horizonPx2 = (float)(m2_geom * camera.pixelsPerM);
                float minPx = 4.0f;
                float diskPx2 = (float)(6.0 * m2_geom * camera.pixelsPerM);
                float jetPx2  = (sim.merger.massSolar >= 1.0e6)
                              ? std::max(40.0f, (float)(30.0 * m2_geom * camera.pixelsPerM))
                              : 0.0f;
                renderer.drawMergerBH(secPos, std::max(minPx, horizonPx2),
                                      diskPx2, jetPx2);
            } else {
                renderer.drawMergerCompactSecondary(
                    secPos, static_cast<int>(sim.merger.secondaryKind),
                    sim.merger.phi);
            }
            // The primary horizon/disk are already drawn at bhCenter, which has
            // been shifted to its barycentric position, so the wobble is visible
            // directly — no extra primary glow needed.
        }

        // Merger: white flash at coalescence
        if (sim.merger.active && sim.merger.flashTimer > 0.0) {
            float alpha = (float)(sim.merger.flashTimer / Simulation::MergerState::FLASH_DURATION);
            renderer.drawMergeFlash(alpha, camera.viewWidth, camera.viewHeight);
        }

        // Merger menu overlay (drawn on top of everything)
        if (ui.mergerMenu.open)
            renderer.drawMergerMenu(ui.mergerMenu, BH2D_PRESETS, NUM_BH2D_PRESETS,
                                    camera.viewWidth, camera.viewHeight,
                                    MERGER_SECONDARY_PRESETS, NUM_MERGER_SECONDARY_PRESETS);

        // Custom-body creator menu (also a modal, mutually exclusive with merger)
        if (ui.customBodyMenu.open) {
            GalaxyBodyType t = UIState::customBodyTypes()[ui.customBodyMenu.typeIdx];
            renderer.drawCustomBodyMenu(ui.customBodyMenu, bodyTypeName(t),
                                        sim.customPresets, bodyTypeName,
                                        camera.viewWidth, camera.viewHeight);
        }

        // Notification (bottom)
        if (ui.notificationTimer > 0) {
            renderer.drawNotification(ui.notification,
                                      camera.viewWidth, camera.viewHeight);
            ui.notificationTimer--;
        }

        renderer.endFrame();
    }

    // using exit(0) instead of return 0 because SFML's static cleanup on macOS can occasionally
    // hang waiting for threads that have already finished. exit() skips the destructors entirely.
    // ugly but functional, and I'd rather not track down a platform-specific cleanup hang right now.
    exit(0);
}
