#include "simulation_2d_widget.h"

// 2D simulation headers (header-only, so only included here)
#include "2D-physics/units.hpp"
#include "2D-simulation/simulation.hpp"
#include "2D-rendering/camera.hpp"
#include "2D-rendering/renderer.hpp"
#include "2D-visualization/ray_visualizer.hpp"
#include "2D-visualization/orbit_visualizer.hpp"
#include "2D-ui/controls.hpp"
#include "2D-utils/presets_2d.hpp"
#include "2D-utils/key_config_2d.hpp"

#include <QDir>
#include <QStandardPaths>
#include <QJsonObject>
#include <QSettings>

#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

// ─────────────────────────────────────────────────────────────
// HUD text builder (verbatim from BlackHole2D.cpp)
// ─────────────────────────────────────────────────────────────
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

        if (sim.galaxySystemActive) {
            int nStars = 0, nGas = 0, nCluster = 0, nDwarf = 0;
            for (const auto& b : sim.bodies) {
                if (!b.isGalaxyBody) continue;
                switch (b.bodyType) {
                    case GalaxyBodyType::Star:           nStars++;   break;
                    case GalaxyBodyType::GasCloud:       nGas++;     break;
                    case GalaxyBodyType::StellarCluster: nCluster++; break;
                    case GalaxyBodyType::DwarfGalaxy:    nDwarf++;   break;
                    case GalaxyBodyType::NeutronStar:    nStars++;   break;
                    case GalaxyBodyType::WhiteDwarf:     nStars++;   break;
                    case GalaxyBodyType::CompanionStar:  nStars++;   break;
                }
            }
            ss << "Galaxy system: ";
            if (nStars)   ss << nStars   << " stars ";
            if (nGas)     ss << nGas     << " gas clouds ";
            if (nCluster) ss << nCluster << " clusters ";
            if (nDwarf)   ss << nDwarf   << " dwarfs ";
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
    if (ui.timeScale != 1.0) {
        std::ostringstream ts;
        if (ui.timeScale >= 1.0) ts << (int)ui.timeScale;
        else ts << std::fixed << std::setprecision(2) << ui.timeScale;
        ss << "Speed: " << ts.str() << "x  |  ";
    }
    ss << "Press / for controls";
    return ss.str();
}

// ─────────────────────────────────────────────────────────────
// Lifecycle
// ─────────────────────────────────────────────────────────────

Simulation2DWidget::Simulation2DWidget(const QString &workspaceName, QWidget *parent)
    : QSFMLCanvas(parent, QSize(800, 600))
    , m_workspaceName(workspaceName)
{}

Simulation2DWidget::~Simulation2DWidget() = default;   // types are complete here

void Simulation2DWidget::setWorkspaceName(const QString &name) {
    m_workspaceName = name;
    if (sim_) sim_->exportName = name.toStdString();
}

void Simulation2DWidget::setLightMode(bool on) {
    if (ui_) ui_->lightMode = on;
    // Renderer is re-synced from ui_->lightMode every frame in onUpdate().
}

// ─────────────────────────────────────────────────────────────
// State serialization / deserialization
// ─────────────────────────────────────────────────────────────

QJsonObject Simulation2DWidget::getState() const
{
    if (!sim_ || !camera_ || !ui_) return {};

    QJsonObject o;
    // Black hole + sim params
    o["bhMass"]             = sim_->bh.metric.M;
    o["numRays"]            = sim_->params.numRays;
    o["pixelsPerM"]         = sim_->params.pixelsPerM;
    o["rMaxIntegrate"]      = sim_->params.rMaxIntegrate;
    o["galaxySystemActive"] = sim_->galaxySystemActive;
    o["activePresetIdx"]    = sim_->activePresetIdx;
    // UI state
    o["presetActive"]       = ui_->presetActive;
    o["presetIdx"]          = ui_->presetIdx;
    o["timeScale"]          = ui_->timeScale;
    o["showRays"]           = ui_->showRays;
    o["showPhotonSphere"]   = ui_->showPhotonSphere;
    o["showInfluenceZones"] = ui_->showInfluenceZones;
    o["showGalaxySystem"]   = ui_->showGalaxySystem;
    o["paused"]             = ui_->paused;
    // Primary body
    if (!sim_->bodies.empty()) {
        const auto &b = sim_->bodies[0];
        o["body0_a"]   = b.nominalA;
        o["body0_ecc"] = b.nominalEcc;
        o["body0_r"]   = b.r;
        o["body0_phi"] = b.phi;
        o["body0_vr"]  = b.vr;
        o["body0_E"]   = b.E;
        o["body0_L"]   = b.L;
    }
    return o;
}

void Simulation2DWidget::setPendingState(const QJsonObject &obj)
{
    pendingState_    = obj;
    hasPendingState_ = !obj.isEmpty();
}

// Applies a previously-stored state blob to the already-initialised objects.
static void applySimState(const QJsonObject &o,
                          Simulation &sim,
                          Camera     &camera,
                          UIState    &ui,
                          unsigned int windowHeight)
{
    auto dv = [&](const char *k, double def) -> double {
        return o.contains(k) ? o[k].toDouble(def) : def;
    };
    auto iv = [&](const char *k, int def) -> int {
        return o.contains(k) ? o[k].toInt(def) : def;
    };
    auto bv = [&](const char *k, bool def) -> bool {
        return o.contains(k) ? o[k].toBool(def) : def;
    };

    sim.bh.metric.M          = dv("bhMass",          sim.bh.metric.M);
    sim.params.numRays        = iv("numRays",         sim.params.numRays);
    sim.params.pixelsPerM     = dv("pixelsPerM",      sim.params.pixelsPerM);
    sim.params.rMaxIntegrate  = dv("rMaxIntegrate",   sim.params.rMaxIntegrate);
    sim.galaxySystemActive    = bv("galaxySystemActive", sim.galaxySystemActive);
    sim.activePresetIdx       = iv("activePresetIdx", sim.activePresetIdx);

    ui.presetActive    = bv("presetActive",    ui.presetActive);
    ui.presetIdx       = iv("presetIdx",       ui.presetIdx);
    ui.timeScale       = dv("timeScale",       ui.timeScale);
    ui.showRays        = bv("showRays",        ui.showRays);
    ui.showPhotonSphere= bv("showPhotonSphere",ui.showPhotonSphere);
    ui.showInfluenceZones = bv("showInfluenceZones", ui.showInfluenceZones);
    ui.showGalaxySystem   = bv("showGalaxySystem",   ui.showGalaxySystem);
    ui.paused          = bv("paused",          ui.paused);

    if (!sim.bodies.empty() && o.contains("body0_r")) {
        auto &b    = sim.bodies[0];
        b.nominalA   = dv("body0_a",   b.nominalA);
        b.nominalEcc = dv("body0_ecc", b.nominalEcc);
        b.r          = dv("body0_r",   b.r);
        b.phi        = dv("body0_phi", b.phi);
        b.vr         = dv("body0_vr",  b.vr);
        b.E          = dv("body0_E",   b.E);
        b.L          = dv("body0_L",   b.L);
        b.trail.clear();   // stale trail from a previous run
    }

    camera.pixelsPerM = sim.params.pixelsPerM;
    sim.rebuildPhotons(windowHeight);
}

void Simulation2DWidget::onInit()
{
    auto sz = getSize();
    const float w = static_cast<float>(std::max(sz.x, 1u));
    const float h = static_cast<float>(std::max(sz.y, 1u));

    sim_      = std::make_unique<Simulation>();
    sim_->exportName = m_workspaceName.toStdString();
    camera_   = std::make_unique<Camera>();
    camera_->pixelsPerM = sim_->params.pixelsPerM;
    camera_->viewWidth  = w;
    camera_->viewHeight = h;
    renderer_ = std::make_unique<Renderer>(static_cast<sf::RenderWindow&>(*this));
    ui_       = std::make_unique<UIState>();

    sim_->rebuildPhotons(static_cast<int>(h));
    sim_->loadCustomPresets();
    renderer_->resetView(w, h);
    clock_.restart();

    // Adopt the global light/dark setting at startup.
    ui_->lightMode = QSettings("Aetherion", "AetherionSuite")
                         .value("ui/lightMode", false).toBool();

    // Load saved 2D keybinds (falls back to defaults automatically)
#ifdef Q_OS_MACOS
    const std::string cfgPath =
        (QDir::homePath() + "/Library/Application Support/Aetherion/blackhole2d_keybinds.cfg")
        .toStdString();
#elif defined(Q_OS_WIN)
    const std::string cfgPath =
        (QStandardPaths::writableLocation(QStandardPaths::AppDataLocation)
         + "/blackhole2d_keybinds.cfg")
        .toStdString();
#else
    const std::string cfgPath =
        (QDir::homePath() + "/.local/share/Aetherion/blackhole2d_keybinds.cfg")
        .toStdString();
#endif
    keyConfig_ = loadKeyConfig2D(cfgPath);

    // Restore a previously-saved workspace state if one was supplied.
    if (hasPendingState_) {
        applySimState(pendingState_, *sim_, *camera_, *ui_,
                      static_cast<unsigned int>(h));
        hasPendingState_ = false;
    }
}

// ─────────────────────────────────────────────────────────────
// Render loop, mirrors BlackHole2D.cpp's draw block.
// QSFMLCanvas::paintEvent calls onUpdate() then display(), so
// we must NOT call renderer_->endFrame() (which also calls display()).
// ─────────────────────────────────────────────────────────────
void Simulation2DWidget::onUpdate()
{
    if (!sim_ || !camera_ || !renderer_ || !ui_) return;

    // Detect resize and update camera/view accordingly
    const auto sz = getSize();
    const float w = static_cast<float>(std::max(sz.x, 1u));
    const float h = static_cast<float>(std::max(sz.y, 1u));

    if (w != camera_->viewWidth || h != camera_->viewHeight) {
        camera_->updateSize(w, h);
        renderer_->resetView(w, h);
        sim_->rebuildPhotons(static_cast<int>(h));
    }

    // Physics update
    const double dt = clock_.restart().asSeconds();
    if (!ui_->paused)
        sim_->update(dt * ui_->timeScale);

    // Auto-zoom during merger inspiral so both BHs stay framed even at
    // extreme mass ratios (parity with standalone BlackHole2D).
    if (sim_->merger.active && sim_->merger.flashTimer <= 0.0 && sim_->merger.r_M > 0.0) {
        const double sepM = sim_->merger.r_M * sim_->bh.metric.M;
        const double halfExtentM = sepM * 1.35 + 4.0 * sim_->bh.metric.M;
        const double targetHalfPx = std::min(camera_->viewWidth, camera_->viewHeight) * 0.40;
        double targetPpm = targetHalfPx / std::max(halfExtentM, 1e-9);
        const double ppmHorizonCap = 125.0 / std::max(sim_->bh.metric.M, 1e-9);
        targetPpm = std::min(targetPpm, ppmHorizonCap);
        const double k = std::clamp(dt * 2.5, 0.0, 0.25);
        sim_->params.pixelsPerM = sim_->params.pixelsPerM + (targetPpm - sim_->params.pixelsPerM) * k;
    }

    camera_->pixelsPerM  = sim_->params.pixelsPerM;
    const sf::Vector2f barycentre = camera_->center();
    sf::Vector2f center = barycentre;
    const double M = sim_->bh.metric.M;

    // Barycentric binary mode (Gaia BH1/2/3, OJ 287): shift the primary off
    // screen-centre by  -q * r_comp  so both objects visibly orbit the shared
    // center of mass. The companion world coords are scaled by baryScale at
    // render time.
    bool  barycentricMode = false;
    float baryScale       = 1.0f;
    if (ui_->presetActive && ui_->showGalaxySystem
        && BH2D_PRESETS[ui_->presetIdx].isBinaryWithBarycenter
        && !sim_->bodies.empty()) {
        const auto&  preset = BH2D_PRESETS[ui_->presetIdx];
        const double M_comp = preset.companionMassSolar;
        const double M_bh   = preset.massSolar;
        const double q      = M_comp / (M_bh + M_comp);
        baryScale = static_cast<float>(1.0 - q);
        const double cx = sim_->bodies[0].worldX();
        const double cy = sim_->bodies[0].worldY();
        center = camera_->worldToScreen(static_cast<float>(-cx * q),
                                        static_cast<float>(-cy * q));
        barycentricMode = true;
    }

    // Merger barycentric mode: shift the primary draw position so both BHs
    // visibly orbit the common centre of mass during inspiral. All primary
    // visuals (horizon, disk, photon sphere, bodies) draw relative to `center`
    // so they pick this up automatically.
    if (sim_->merger.active && sim_->merger.flashTimer <= 0.0) {
        const double q   = sim_->merger.massRatioQ;
        // Skips the cosmetic wobble for extreme-mass-ratio mergers where the
        // primary really shouldn’t budge.
        if (q >= 0.05) {
            const double sep = sim_->merger.r_M * M;
            const double bx  = q * sep * std::cos(sim_->merger.phi);
            const double by  = q * sep * std::sin(sim_->merger.phi);
            center.x = barycentre.x - static_cast<float>(bx * camera_->pixelsPerM);
            center.y = barycentre.y + static_cast<float>(by * camera_->pixelsPerM);
        }
    }

    // ── Draw ──────────────────────────────────────────────────
    renderer_->setLightMode(ui_->lightMode);
    renderer_->beginFrame();
    renderer_->drawStarfield(bgClock_.getElapsedTime().asSeconds());

    // Influence zones
    if (ui_->presetActive && ui_->showInfluenceZones && ui_->showGalaxySystem) {
        const auto& preset = BH2D_PRESETS[ui_->presetIdx];
        if (preset.isGalacticCenter) {
            const float ppm    = static_cast<float>(camera_->pixelsPerM);
            const float soiPx  = static_cast<float>(preset.zones.sphereOfInfluenceM * M * ppm);
            const float bondiPx= static_cast<float>(preset.zones.bondiRadiusM        * M * ppm);
            const float tidalPx= static_cast<float>(preset.zones.tidalDisruptionM   * M * ppm);
            renderer_->drawInfluenceZone(center, soiPx,   sf::Color(100, 255, 100,  80), "Sphere of Influence");
            renderer_->drawInfluenceZone(center, bondiPx, sf::Color(255, 200,  50,  90), "Bondi Radius");
            renderer_->drawInfluenceZone(center, tidalPx, sf::Color(255,  80,  80, 120), "Tidal Disruption");
            if (preset.zones.hasJetCones)
                renderer_->drawJetCones(center, soiPx * 1.5f);
        }
    }

    renderer_->drawHorizon(center,
        static_cast<float>(sim_->bh.metric.horizon() * camera_->pixelsPerM));

    if (ui_->presetActive) {
        const float diskRadius = static_cast<float>(6.0 * M * camera_->pixelsPerM);
        renderer_->drawAccretionDisk(center, diskRadius);
    }

    if (ui_->showPhotonSphere)
        renderer_->drawPhotonSphere(center,
            static_cast<float>(sim_->bh.metric.photonSphere() * camera_->pixelsPerM));

    // Rays
    if (ui_->showRays) {
        for (const auto& photon : sim_->photons) {
            if (photon.captured) {
                auto [v0, v1] = RayVisualizer::capturedRayLine(photon.impactParameter, *camera_);
                renderer_->drawCapturedRay(v0, v1);
            } else {
                RayVisualizer::colorByRedshift(photon, sim_->bh.metric, *camera_, rayVertScratch_);
                renderer_->drawRayPath(rayVertScratch_);
            }
        }
    }

    // Caustic highlights
    if (ui_->showCaustics && !sim_->lensingData.causticPoints.empty()) {
        for (const auto& [avgB, gradient] : sim_->lensingData.causticPoints) {
            renderer_->drawCausticMarker(camera_->worldToScreen(0.0f,  static_cast<float>(avgB)), gradient);
            renderer_->drawCausticMarker(camera_->worldToScreen(0.0f, -static_cast<float>(avgB)), gradient);
        }
    }

    // Time-dilation heatmap
    if (ui_->showTimeDilationMap)
        renderer_->drawTimeDilationMap(center, M, camera_->pixelsPerM);

    // Barycenter marker (drawn before bodies so it sits behind them)
    if (barycentricMode)
        renderer_->drawBarycenterMarker(barycentre);

    // Bodies and orbit paths
    for (size_t bodyIdx = 0; bodyIdx < sim_->bodies.size(); ++bodyIdx) {
        const auto& body     = sim_->bodies[bodyIdx];
        const float scale    = (barycentricMode && body.isGalaxyBody) ? baryScale : 1.0f;
        const auto bodyVis   = OrbitVisualizer::computeBodyVisual(body, sim_->bh.metric, *camera_, scale);
        if (body.isSecondaryBH) {
            const float secHorizonPx = 6.0f;
            const float secDiskPx    = 16.0f;
            renderer_->drawMergerBH(bodyVis.screenPos, secHorizonPx,
                                    secDiskPx, 0.0f, body.label.c_str());
        } else {
            renderer_->drawBody(bodyVis);
        }

        // Pulsar-specific visual decoration (jets, magnetic field lines, LC ring, flux tubes)
        if (sim_->activeScenario == ResearchScenario::PulsarOrbital && body.isPulsar
            && !body.captured) {
            float lc_px = 0.0f;
            if (sim_->pulsarState.lightCylRadius_m > 0.0)
                lc_px = static_cast<float>(sim_->pulsarState.lightCylRadius_m * camera_->pixelsPerM);
            renderer_->drawPulsarBody(
                bodyVis.screenPos,
                static_cast<float>(sim_->pulsarState.spinPhase),
                static_cast<float>(sim_->pulsarState.precPhase),
                lc_px,
                sim_->pulsarData.inLightCylinder,
                center,
                sim_->pulsarData.magPower_ergs
            );
        }

        const auto orbitPath = OrbitVisualizer::computeOrbitPath(body, sim_->bh.metric, *camera_, scale);
        renderer_->drawOrbitPath(orbitPath);

        if (body.isGalaxyBody && ui_->showGalaxySystem && !body.label.empty() && !body.isSecondaryBH) {
            sf::Color labelColor;
            switch (body.bodyType) {
                case GalaxyBodyType::Star:           labelColor = sf::Color(255, 240, 200); break;
                case GalaxyBodyType::GasCloud:       labelColor = sf::Color(255, 150, 100); break;
                case GalaxyBodyType::StellarCluster: labelColor = sf::Color(200, 220, 255); break;
                case GalaxyBodyType::DwarfGalaxy:    labelColor = sf::Color(255, 200, 255); break;
                case GalaxyBodyType::NeutronStar:    labelColor = sf::Color(160, 200, 255); break;
                case GalaxyBodyType::WhiteDwarf:     labelColor = sf::Color(200, 215, 255); break;
                case GalaxyBodyType::CompanionStar:  labelColor = sf::Color(255, 200, 140); break;
            }
            renderer_->drawBodyLabel(bodyVis.screenPos, body.label, labelColor);
        }

        if (!body.isGalaxyBody && !body.label.empty() &&
            sim_->activeScenario != ResearchScenario::None) {
            renderer_->drawBodyLabel(bodyVis.screenPos, body.label, sf::Color(180, 220, 255));
        }
    }

    renderer_->drawHUD(buildHUD(*sim_, *ui_));

    // ---- Tidal disruption event overlay ----
    if (sim_->tidalEvent.active) {
        sf::Vector2f center = camera_->center();
        sf::Vector2f flashPos = {
            center.x + static_cast<float>(sim_->tidalEvent.eventX * camera_->pixelsPerM),
            center.y - static_cast<float>(sim_->tidalEvent.eventY * camera_->pixelsPerM)
        };
        float flashAlpha = (sim_->tidalEvent.flashTimer > 0.0)
            ? static_cast<float>(sim_->tidalEvent.flashTimer /
                                 Simulation::TidalEvent::FLASH_DURATION)
            : 0.0f;
        bool showLabel = flashAlpha > 0.0f || !sim_->tidalEvent.particles.empty();
        std::vector<Renderer::TidalParticleVis> pvis;
        pvis.reserve(sim_->tidalEvent.particles.size());
        for (const auto& p : sim_->tidalEvent.particles) {
            float lifeF = static_cast<float>(p.lifetime / p.maxLifetime);
            pvis.push_back({
                center.x + static_cast<float>(p.x * camera_->pixelsPerM),
                center.y - static_cast<float>(p.y * camera_->pixelsPerM),
                p.size * lifeF,
                lifeF,
                p.isFallback
            });
        }
        renderer_->drawTidalEvent(flashPos, flashAlpha, pvis, showLabel,
            sim_->tidalEvent.kind == Simulation::TidalEventKind::MergerShockwave ? "MERGER REMNANT"
          : sim_->tidalEvent.kind == Simulation::TidalEventKind::XrayBurst       ? "X-RAY BURST"
          : "TIDAL DISRUPTION EVENT");
    }

    if (ui_->showControlsPanel)
        renderer_->drawControlsPanel(camera_->viewHeight);

    if (ui_->showDataPanel)
        renderer_->drawDataPanel(sim_->formatDataPanel(), camera_->viewWidth, camera_->viewHeight);

    // Merger BH (inspiral body) + death-spiral trails
    if (sim_->merger.active && sim_->merger.flashTimer <= 0.0) {
        const double q   = sim_->merger.massRatioQ;
        const double sep = sim_->merger.r_M * sim_->bh.metric.M;

        // Death-spiral trails (barycentre-centred world coords → draw from screen centre)
        renderer_->drawMergerTrail(sim_->merger.trail2, barycentre,
                                   camera_->pixelsPerM,
                                   sf::Color(180, 100, 255));   // secondary (purple)
        if (q > 0.01) {
            renderer_->drawMergerTrail(sim_->merger.trail1, barycentre,
                                       camera_->pixelsPerM,
                                       sf::Color(100, 220, 255)); // primary (cyan)
        }

        // Secondary BH at +(1-q)*sep from the barycentre. Disk + jets scaled
        // from its actual mass so an SMBH secondary visibly outsizes a stellar
        // one and AGN-class secondaries get jets. NS/Pulsar/Star/WD render as
        // a small bright dot instead.
        const double r2 = sep * (1.0 - q);
        float sx = barycentre.x + static_cast<float>(r2 * std::cos(sim_->merger.phi) * camera_->pixelsPerM);
        float sy = barycentre.y - static_cast<float>(r2 * std::sin(sim_->merger.phi) * camera_->pixelsPerM);
        if (sim_->merger.secondaryKind == Simulation::MergerSecondaryKind::BlackHole) {
            double m2_geom = units::solarMassToGeomMeters(sim_->merger.massSolar);
            float horizonPx2 = static_cast<float>(m2_geom * camera_->pixelsPerM);
            float minPx = 4.0f;
            float diskPx2 = static_cast<float>(6.0 * m2_geom * camera_->pixelsPerM);
            float jetPx2  = (sim_->merger.massSolar >= 1.0e6)
                          ? std::max(40.0f, static_cast<float>(30.0 * m2_geom * camera_->pixelsPerM))
                          : 0.0f;
            renderer_->drawMergerBH({sx, sy}, std::max(minPx, horizonPx2),
                                    diskPx2, jetPx2);
        } else {
            renderer_->drawMergerCompactSecondary(
                {sx, sy}, static_cast<int>(sim_->merger.secondaryKind),
                sim_->merger.phi);
        }
        // Primary already drawn at the shifted `center` above.
    }

    // Merger white flash
    if (sim_->merger.active && sim_->merger.flashTimer > 0.0) {
        float alpha = (float)(sim_->merger.flashTimer / Simulation::MergerState::FLASH_DURATION);
        renderer_->drawMergeFlash(alpha, camera_->viewWidth, camera_->viewHeight);
    }

    // Merger menu overlay
    if (ui_->mergerMenu.open)
        renderer_->drawMergerMenu(ui_->mergerMenu, BH2D_PRESETS, NUM_BH2D_PRESETS,
                                   camera_->viewWidth, camera_->viewHeight,
                                   MERGER_SECONDARY_PRESETS, NUM_MERGER_SECONDARY_PRESETS);

    // Custom-body creator menu
    if (ui_->customBodyMenu.open) {
        GalaxyBodyType t = UIState::customBodyTypes()[ui_->customBodyMenu.typeIdx];
        renderer_->drawCustomBodyMenu(ui_->customBodyMenu, bodyTypeName(t),
                                       sim_->customPresets, bodyTypeName,
                                       camera_->viewWidth, camera_->viewHeight);
    }

    if (ui_->notificationTimer > 0) {
        renderer_->drawNotification(ui_->notification, camera_->viewWidth, camera_->viewHeight);
        ui_->notificationTimer--;
    }
    // paintEvent calls display() after onUpdate(), do NOT call endFrame() here
}

// ─────────────────────────────────────────────────────────────
// Input, forward to the 2D controls handler
// ─────────────────────────────────────────────────────────────

void Simulation2DWidget::onKeyPressed(sf::Keyboard::Key code)
{
    if (!sim_ || !ui_ || !camera_) return;
    // handleInput() takes sf::Event, construct a KeyPressed event wrapper
    sf::Event ev{ sf::Event::KeyPressed{ code, sf::Keyboard::Scancode::Unknown, false, false } };
    handleInput(ev, *sim_, *ui_, static_cast<unsigned int>(camera_->viewHeight), keyConfig_);
}

void Simulation2DWidget::onKeyReleased(sf::Keyboard::Key /*code*/)
{
    // 2D controls only act on KeyPressed; nothing to do here
}

void Simulation2DWidget::onTextEntered(char32_t unicode)
{
    if (!sim_ || !ui_ || !camera_) return;
    if (!ui_->mergerMenu.open || !ui_->mergerMenu.inputtingCustom) return;
    sf::Event ev{ sf::Event::TextEntered{ unicode } };
    handleInput(ev, *sim_, *ui_, static_cast<unsigned int>(camera_->viewHeight), keyConfig_);
}

void Simulation2DWidget::onMouseMoved(float /*x*/, float /*y*/)
{
    // 2D sim doesn't use mouse look
}

void Simulation2DWidget::onMousePressed(sf::Mouse::Button /*button*/, float /*x*/, float /*y*/)
{
    // Reserved for future 2D mouse-click interactions
}

void Simulation2DWidget::onMouseWheelScrolled(float delta, float /*x*/, float /*y*/)
{
    if (!sim_ || !camera_) return;
    // Zoom: replicate the Up/Down arrow key zoom from BlackHole2D
    const float factor = (delta > 0.0f) ? 1.15f : (1.0f / 1.15f);
    sim_->params.pixelsPerM *= factor;
    camera_->pixelsPerM      = sim_->params.pixelsPerM;
    sim_->rebuildPhotons(static_cast<int>(camera_->viewHeight));
}
