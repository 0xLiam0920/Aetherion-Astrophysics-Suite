#pragma once
#include "black_hole.hpp"
#include "photon.hpp"
#include "orbiting_body.hpp"
#include "research_data.hpp"
#include "csv_export.hpp"
#include "FITS_export.hpp"
#include "../2D-utils/presets_2d.hpp"
#include "../2D-physics/geodesic.hpp"
#include "../2D-physics/pulsar_orbital.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <filesystem>

// TODO: pack all viable headers into single header hpp files for 2d, 3d, and launcher
// hours wasted trying to fix linux package diff issues: 7

struct SimParams {
    int    numRays        = 120;
    double rMaxIntegrate  = 1e5;
    double pixelsPerM     = 60.0;
    float  timeScale      = 1.0f;   // user-controllable speed multiplier
};

/*--------- Research scenario types ---------*/
enum class ResearchScenario {
    None,
    ISCOTest,          // Particles at 5M, 6M, 7M
    PhotonSphereTest,  // Photon near critical b
    RadialInfall,      // Particle dropped from rest
    TidalDisruption,   // Body in tidal zone
    PulsarOrbital      // Neutron star inspiraling via GW emission
};

// Central simulation manager — orchestrates the black hole world.
class Simulation {
public:
    BlackHole                  bh;
    std::vector<Photon>        photons;
    std::vector<OrbitingBody>  bodies;
    SimParams                  params;
    bool                       galaxySystemActive = false;
    int                        activePresetIdx    = -1;

    // Research mode state
    LensingData                lensingData;
    std::vector<PhotonSphereTestResult> photonSphereResults;
    std::vector<ISCOTestResult>         iscoResults;
    ResearchScenario           activeScenario = ResearchScenario::None;
    int                        selectedBodyIdx = 0;
    double                     gasTemperatureK = 1e7;
    // Ambient gas density used by the Bondi-Hoyle accretion model.
    // Default ~1.7e-19 kg/m³ ≈ 100 protons/cm³ (typical hot ISM, Sgr A* environment).
    double                     gasDensity_kgm3 = 1.7e-19;
    // Running simulation time in geometric units of M; advanced in update(dt).
    double                     simTime_M       = 0.0;
    // Bondi-Hoyle-Lyttleton accretion + bolometric luminosity time series.
    BondiAccretionTracker      accretion;
    // Snapshot diagnostics rebuilt on export: analytic gas radial profile and
    // multi-color disk SED. Held here so they can also feed the data panel.
    GasRadialProfile           gasProfile;
    SEDSnapshot                sed;
    std::string                exportName      = "untitled_data"; ///< used as export subfolder name

    // Pulsar orbital scenario state
    PulsarState        pulsarState;
    PulsarOrbitalData  pulsarData;
    bool               tidalDisrupted = false; // set when tidal test body crosses disruption threshold

    // ── Tidal disruption event state ─────────────────────────────────────────
    struct TidalDebrisParticle {
        double x, y;           // world position (geometric units, BH at origin)
        double vx, vy;         // velocity (geometric units / wall-second)
        double lifetime;       // remaining lifetime (wall-seconds)
        double maxLifetime;
        float  size;           // base radius in pixels (at reference zoom)
        bool   isFallback;     // true = inward gas-stream particle
    };

    struct TidalEvent {
        bool   active         = false;
        double flashTimer     = 0.0;   // seconds remaining for the localised flash
        double eventX         = 0.0;   // disruption world X (geometric units)
        double eventY         = 0.0;   // disruption world Y (geometric units)
        std::vector<TidalDebrisParticle> particles;
        static constexpr double FLASH_DURATION  = 0.45;  // seconds
        static constexpr double DEBRIS_LIFETIME = 5.0;   // seconds
        static constexpr double STREAM_LIFETIME = 8.0;   // seconds (fallback gas)
    } tidalEvent;

    // Tidal disruption threshold in units of M.
    // Set from the active preset's zones.tidalDisruptionM; default matches the
    // startTidalDisruption() initial orbit so the scenario always triggers cleanly.
    double tidalRadiusM = 8.0;

    // ── Merger event state ───────────────────────────────────────────────────
    struct MergerState {
        bool   active          = false;
        bool   completed       = false;
        double r_M             = 40.0;   // separation in units of primary M (scales with M)
        double phi             = 0.0;    // orbital angle (radians)
        double massSolar       = 0.0;    // incoming BH mass in solar masses
        double secondaryMassGeom = 0.0;  // precomputed geometric mass of secondary
        double flashTimer      = 0.0;
        double preMergeMassSolar = 0.0;
        // Mass ratio q = M2/(M1+M2) ∈ [0,0.5].
        // q > 0.1 → both BHs visibly orbit the common barycenter ("death spiral").
        double massRatioQ      = 0.0;
        // Death-spiral trail: world-space positions of the SECONDARY in BH-centred
        // coords, accumulated at ~30 fps to draw the fading streak.
        static constexpr int MAX_TRAIL = 400;
        std::deque<std::pair<double,double>> trail2;  // secondary trail
        std::deque<std::pair<double,double>> trail1;  // primary barycentric wobble trail
        // Inspiral coefficient tuned for ~15s from r=40 to r=3.5.
        static constexpr double INSPIRAL_K    = 42670.0;
        static constexpr double OMEGA_K       = 12.0;
        static constexpr double FLASH_DURATION = 0.6;
        static constexpr double MERGE_RADIUS_M = 3.5;
    } merger;

    Simulation() {
        bodies.emplace_back(bh.metric, 10.0 * bh.metric.M, 0.3);
    }

    void update(double dt) {
        // Scale simulation time with M so animations run at similar visual
        // speed regardless of BH mass.  For the default M=1 case, dtSim=dt.
        // dt * M is the correct normalization: orbital period ∝ M, so this
        // keeps the fraction-of-orbit-per-frame mass-independent.
        double dtSim = dt * bh.metric.M * static_cast<double>(params.timeScale);

        for (auto& body : bodies)
            body.update(bh.metric, dtSim);

        // Advance simulation clock and update Bondi-Hoyle accretion model.
        // We feed the analytic BHL formula with the *current* selected body's
        // orbital v/c so Ṁ and L_bol modulate naturally with the orbit phase
        // (denominator (v²+c_s²)^(3/2) drops at apoapsis, peaks at periapsis).
        // simTime_M is in units of M (geometric); divide by M to keep the
        // FITS time column dimensionless across BH masses.
        simTime_M += dtSim / std::max(bh.metric.M, 1e-30);
        {
            const double M_BH_solar = currentMassSolar();
            const double cs_over_c  = Schwarzschild::soundSpeedFromTemp(gasTemperatureK);
            double v_over_c = 0.0;
            if (!bodies.empty() && selectedBodyIdx >= 0
                && selectedBodyIdx < (int)bodies.size())
            {
                v_over_c = bodies[selectedBodyIdx].measurement.orbitalVelocity;
            }
            accretion.update(simTime_M, M_BH_solar, cs_over_c,
                             gasDensity_kgm3, v_over_c);
        }

        // Update ISCO test results if active
        if (activeScenario == ResearchScenario::ISCOTest)
            updateISCOResults();

        // Update pulsar orbital scenario if active
        if (activeScenario == ResearchScenario::PulsarOrbital)
            updatePulsar(dtSim);

        // Update tidal disruption scenario if active
        if (activeScenario == ResearchScenario::TidalDisruption)
            updateTidalDisruption(dtSim);

        // Advance tidal event particle system (always, even after scenario completes)
        updateTidalEvent(dt);

        // Advance merger inspiral if active
        if (merger.active)
            updateMerger(dt); // wall-clock dt, not dtSim — so the animation plays at a fixed visual rate
    }

    // Black hole mass in solar masses, derived from the active preset.
    // Falls back to 0 (which the accretion tracker treats as "no data") when
    // no preset is selected, so we never feed garbage into the BHL formula.
    double currentMassSolar() const {
        if (activePresetIdx >= 0 && activePresetIdx < NUM_BH2D_PRESETS)
            return BH2D_PRESETS[activePresetIdx].massSolar;
        return 0.0;
    }

    void rebuildPhotons(unsigned int windowHeight, bool highRes = false) {
        photons.clear();
        double halfHeightM = ((double)windowHeight * 0.5) / params.pixelsPerM;

        // 120 rays uniformly spaced in impact parameter b across the visible screen height.
        // uniform-in-b is NOT the ideal sampling near b_crit, where deflection changes
        // really fast and you'd want to cluster more rays there. but it looks good on screen
        // and changing it would mess up the visual spacing, so here we are.
        for (int i = 0; i < params.numRays; ++i) {
            double t = (params.numRays <= 1)
                     ? 0.5
                     : (double)i / (double)(params.numRays - 1);
            double b = (2.0 * t - 1.0) * halfHeightM;

            Photon p;
            p.impactParameter = b;
            p.computePath(bh.metric, params.rMaxIntegrate);
            photons.push_back(std::move(p));
        }

        // High-res lensing mode: add extra rays log-spaced near b_crit on both sides.
        // These show the strong-lensing regime in much finer detail without disturbing
        // the uniform background grid.
        if (highRes) {
            double b_crit = bh.metric.criticalImpact();
            constexpr int N_HIGHRES = 20; // rays per side
            for (int sign : {-1, 1}) {
                for (int i = 0; i < N_HIGHRES; ++i) {
                    double t = (double)i / (double)(N_HIGHRES - 1);
                    // epsilon log-spaced from 1e-3 to 0.2 (0.1% to 20% offset from b_crit)
                    double eps = std::exp(std::log(1e-3) + t * (std::log(0.2) - std::log(1e-3)));
                    double b = b_crit * (1.0 + sign * eps);
                    Photon p;
                    p.impactParameter = b;
                    p.computePath(bh.metric, params.rMaxIntegrate, 0.001); // finer angular step near b_crit
                    photons.push_back(std::move(p));
                }
            }
        }

        // Rebuild lensing analysis
        buildLensingData();
    }

    /*--------- Black hole merger event ---------*/
    void startMerger(double massSolar, unsigned int /*windowHeight*/) {
        merger.active          = true;
        merger.completed       = false;
        merger.massSolar       = massSolar;
        merger.r_M             = 40.0;
        merger.phi             = M_PI / 4.0;  // start offset so it’s immediately visible
        merger.flashTimer      = 0.0;
        merger.preMergeMassSolar = 0.0;
        merger.secondaryMassGeom  = units::solarMassToGeomMeters(massSolar);
        merger.trail2.clear();
        merger.trail1.clear();
        // Mass ratio q = M2/(M1+M2). For q > 0.1 both BHs visibly orbit the
        // barycenter — the death spiral.
        const double m1Sol = currentMassSolar();
        const double totM  = m1Sol + massSolar;
        merger.massRatioQ = (totM > 0.0) ? massSolar / totM : 0.0;
    }

    // Apply gravitational perturbation from the in-spiraling secondary BH onto all
    // orbiting bodies. Uses Newtonian tidal acceleration (M2/d²) in geometric units,
    // applied as a proper-time-scaled impulse to vr and L, then E is recomputed.
    // The perturbation is physically negligible for large mass-ratio mergers (e.g.,
    // Sgr A* into TON 618) and dramatic for equal-mass mergers (TON 618 + TON 618).
    void applyMergerPerturbation(double dt) {
        const double M2 = merger.secondaryMassGeom;
        if (M2 < 1e-30) return;

        // Secondary BH world position (geometric units, BH at origin)
        const double r_BH = merger.r_M * bh.metric.M;
        const double bx   = r_BH * std::cos(merger.phi);
        const double by   = r_BH * std::sin(merger.phi);

        for (auto& body : bodies) {
            if (body.captured) continue;

            const double px = body.r * std::cos(body.phi);
            const double py = body.r * std::sin(body.phi);

            const double dx = px - bx;  // vector from secondary to body
            const double dy = py - by;
            const double d2 = dx*dx + dy*dy;
            const double d  = std::sqrt(d2);

            // Skip bodies too close to the secondary (they'd be "captured" by it)
            if (d < bh.metric.M * 0.5) {
                body.captured = true;
                continue;
            }

            // Gravitational acceleration toward secondary (geometric units, G=1)
            const double inv_d3 = 1.0 / (d2 * d);
            const double ax = -M2 * dx * inv_d3;  // attraction toward secondary
            const double ay = -M2 * dy * inv_d3;

            // Resolve into body's radial (outward) and tangential (CCW) components
            const double cos_p = std::cos(body.phi);
            const double sin_p = std::sin(body.phi);
            const double ar =  ax * cos_p + ay * sin_p;
            const double at = -ax * sin_p + ay * cos_p;

            // Approximate proper-time step for this body
            const double fr   = bh.metric.f(body.r);
            const double dtau = dt * std::max(fr, 1e-10) / std::max(body.E, 1e-10);

            // Apply impulse to geodesic state
            body.vr += ar * dtau;
            body.L  += at * body.r * dtau;

            // Recompute specific energy from updated (r, vr, L) in current metric.
            // E² = vr² + f(r)·(1 + L²/r²)
            const double new_E2 = body.vr * body.vr
                                 + fr * (1.0 + body.L * body.L / (body.r * body.r));
            if (new_E2 > 0.0)
                body.E = std::sqrt(new_E2);

            // Bodies perturbed into escape trajectories during inspiral fly off screen —
            // that's the point. Don't silently delete them, let them streak away.
            if (body.E > 1.2)
                body.markEjected();
        }
    }

    // Post-merger body update.
    // Galaxy bodies keep their current perturbed (r, phi, vr, L) and have E recomputed
    // under the new combined metric. Unbound bodies get ejected and streak off screen.
    // A GW recoil kick is applied to surviving bodies (the merged BH "jumps" due to
    // asymmetric GW emission; in the BH frame this appears as a sudden velocity shift
    // on all remaining matter). Finally a ring of particles is spawned to represent
    // the expanding EM / gravitational wave pulse.
    void reinitBodiesPostMerger() {
        const double newM    = bh.metric.M;
        const double newISCO = bh.metric.isco();

        // Reset the primary (non-galaxy) test body to a clean orbit
        for (auto& body : bodies) {
            if (body.isGalaxyBody) continue;
            body.initFromKeplerian(bh.metric, 10.0 * newM, body.nominalEcc);
            break;
        }

        // Galaxy bodies: recompute E under new metric, eject the unbound ones
        for (auto& body : bodies) {
            if (!body.isGalaxyBody || body.captured || body.ejected) continue;
            const double fr_new = bh.metric.f(body.r);
            const double new_E2 = body.vr * body.vr
                                 + fr_new * (1.0 + body.L * body.L / (body.r * body.r));
            body.E = (new_E2 > 0.0) ? std::sqrt(new_E2) : 1e-3;

            if (body.E >= 1.0) {
                // The instantaneous mass increase unbound this body — it’s being
                // flung into intergalactic space. Let it visibly streak away.
                body.markEjected();
            } else if (body.r < newISCO * 1.05 && body.vr >= 0.0) {
                body.vr = -0.01 * newM;  // nudge inward so geodesic integrator lets it plunge
            }
        }

        // GW recoil kick: the merged BH recoils due to asymmetric gravitational wave
        // emission. In the BH-centred frame all surviving bodies receive the opposite
        // velocity impulse (−v_kick). For equal-mass mergers this represents the
        // spin-induced “superkick”; for unequal masses it scales with the asymmetry.
        // Magnitude: ~1% of c at q=1 (equal mass), fading toward zero for extreme q.
        {
            const double m1 = merger.preMergeMassSolar;
            const double m2 = merger.massSolar;
            const double q  = (m1 > 0.0 && m2 > 0.0)
                            ? std::min(m1, m2) / std::max(m1, m2)
                            : 0.0;   // mass ratio [0,1]; 1 = equal mass
            // visual kick speed in geometric units; 0.01*newM puts it at ~1% c
            const double v_kick = 0.01 * q * newM;
            // deterministic random direction from mass parameters
            uint32_t seed = (static_cast<uint32_t>(m1) ^ 0xA5A5A5A5u)
                          + static_cast<uint32_t>(m2 * 7919.0);
            seed = seed * 1664525u + 1013904223u;
            const double kick_angle = 2.0 * M_PI * (seed >> 8) / double(1 << 24);
            const double kvx = v_kick * std::cos(kick_angle);
            const double kvy = v_kick * std::sin(kick_angle);

            for (auto& body : bodies) {
                if (body.captured || body.ejected || !body.isGalaxyBody) continue;
                const double cos_p = std::cos(body.phi);
                const double sin_p = std::sin(body.phi);
                body.vr += kvx * cos_p + kvy * sin_p;
                body.L  += body.r * (-kvx * sin_p + kvy * cos_p);
                const double fr2 = bh.metric.f(body.r);
                const double e2  = body.vr * body.vr
                                 + fr2 * (1.0 + body.L * body.L / (body.r * body.r));
                if (e2 > 0.0) body.E = std::sqrt(e2);
                if (body.E >= 1.0) body.markEjected();
            }
        }

        // GW / EM shockwave ring: spawn an outward ring of debris particles centered
        // on the merged BH. This represents the expanding gravitational wave and
        // quasar-level electromagnetic burst that a real ultramassive merger would
        // produce. Reuses the tidal event particle system so no extra rendering needed.
        {
            tidalEvent.active     = true;
            tidalEvent.flashTimer = TidalEvent::FLASH_DURATION * 2.0;  // longer flash for a merger
            tidalEvent.eventX     = 0.0;
            tidalEvent.eventY     = 0.0;

            uint32_t rng = 0xDEADBEEFu ^ static_cast<uint32_t>(newM * 1e-6);
            auto rand01 = [&]() -> double {
                rng = rng * 1664525u + 1013904223u;
                return (rng >> 8) / double(1 << 24);
            };

            // Fast outward ring — the “shockwave” expanding at a large fraction of c
            const double V_RING = 0.5 * newM;
            for (int i = 0; i < 72; ++i) {
                const double angle = 2.0 * M_PI * i / 72.0 + rand01() * 0.08;
                const double speed = V_RING * (0.9 + 0.2 * rand01());
                TidalDebrisParticle p;
                p.x          = 0.0;  p.y = 0.0;
                p.vx         = std::cos(angle) * speed;
                p.vy         = std::sin(angle) * speed;
                p.maxLifetime = 3.5 + 2.0 * rand01();
                p.lifetime    = p.maxLifetime;
                p.size        = 3.5f + (float)(2.0 * rand01());
                p.isFallback  = false;
                tidalEvent.particles.push_back(p);
            }
            // Slower inward-then-outward fallback debris — accretion disk splashback
            for (int i = 0; i < 24; ++i) {
                const double angle = 2.0 * M_PI * rand01();
                const double speed = V_RING * (0.15 + 0.25 * rand01());
                TidalDebrisParticle p;
                p.x          = newM * (rand01() - 0.5) * 4.0;
                p.y          = newM * (rand01() - 0.5) * 4.0;
                p.vx         = std::cos(angle) * speed;
                p.vy         = std::sin(angle) * speed;
                p.maxLifetime = 5.0 + 3.0 * rand01();
                p.lifetime    = p.maxLifetime;
                p.size        = 2.0f + (float)(1.5 * rand01());
                p.isFallback  = true;
                tidalEvent.particles.push_back(p);
            }
        }
    }

    void updateMerger(double dt) {
        if (!merger.active) return;

        // Flash phase: wait for the white flash to finish, then finalize
        if (merger.flashTimer > 0.0) {
            merger.flashTimer -= dt;
            if (merger.flashTimer <= 0.0) {
                merger.flashTimer = 0.0;
                merger.completed  = true;
                merger.active     = false;
            }
            return;
        }

        applyMergerPerturbation(dt);

        // GW inspiral: Peters-formula shape — separation shrinks as r^{-3}
        double r3 = merger.r_M * merger.r_M * merger.r_M;
        if (r3 < 1e-30) r3 = 1e-30;
        merger.r_M -= (MergerState::INSPIRAL_K * dt) / r3;

        // Angular advance — spin up as separation closes (Ω ∝ r^{-3/2})
        double r_safe = std::max(merger.r_M, 0.1);
        merger.phi += (MergerState::OMEGA_K * dt) / std::pow(r_safe, 1.5);

        // Record trail positions in BH-centred geometric units.
        // For similar-mass mergers (q > 0.1), apply barycentric correction:
        //   secondary orbits at  r2 = r * (1 - q)  from the CoM
        //   primary   orbits at  r1 = r * q         from the CoM (opposite side)
        // For extreme mass ratios (q ≈ 0) this reduces to the old behaviour.
        {
            const double q   = merger.massRatioQ;   // M2/(M1+M2) in [0,0.5]
            const double sep = merger.r_M * bh.metric.M;  // geometric separation
            // Secondary position (from BH origin = total-mass centroid)
            const double r2 = sep * (1.0 - q);
            const double sx  =  r2 * std::cos(merger.phi);
            const double sy  =  r2 * std::sin(merger.phi);
            merger.trail2.emplace_back(sx, sy);
            if ((int)merger.trail2.size() > MergerState::MAX_TRAIL)
                merger.trail2.pop_front();

            // Primary barycentric wobble (only visible for q > 0.1)
            if (q > 0.01) {
                const double r1 = sep * q;
                const double px  = -r1 * std::cos(merger.phi);  // opposite side
                const double py  = -r1 * std::sin(merger.phi);
                merger.trail1.emplace_back(px, py);
                if ((int)merger.trail1.size() > MergerState::MAX_TRAIL)
                    merger.trail1.pop_front();
            }
        }

        // Merge trigger
        if (merger.r_M <= MergerState::MERGE_RADIUS_M) {
            double m1Solar = bh.metric.M / units::solarMassToGeomMeters(1.0);
            double mFinalSolar = m1Solar + merger.massSolar;
            merger.preMergeMassSolar = m1Solar;

            double oldM = bh.metric.M;
            bh.metric.M = units::solarMassToGeomMeters(mFinalSolar);
            params.pixelsPerM *= oldM / bh.metric.M;

            activeScenario = ResearchScenario::None;
            reinitBodiesPostMerger();

            merger.flashTimer = MergerState::FLASH_DURATION;
        }
    }

    void reinitBodies() {        double M = bh.metric.M;
        if (!bodies.empty()) {
            auto& b = bodies[0];
            b.initFromKeplerian(bh.metric, 10.0 * M, b.nominalEcc);
        }
        bodies.erase(
            std::remove_if(bodies.begin(), bodies.end(),
                [](const OrbitingBody& b){ return b.isGalaxyBody; }),
            bodies.end());
        if (galaxySystemActive && activePresetIdx >= 0)
            spawnGalaxySystem(activePresetIdx);
    }

    void spawnGalaxySystem(int presetIdx) {
        // Remove all bodies (including the default test body) — it doesn't belong
        // in a galaxy preset and would just orbit with no physical context.
        bodies.clear();

        if (presetIdx < 0 || presetIdx >= NUM_BH2D_PRESETS) return;
        const auto& preset = BH2D_PRESETS[presetIdx];
        if (!preset.isGalacticCenter || preset.numGalaxyBodies == 0) return;

        double M = bh.metric.M;
        for (int i = 0; i < preset.numGalaxyBodies; ++i) {
            const auto& gsb = preset.galaxyBodies[i];
            double aSim = gsb.semiMajorM * M;
            OrbitingBody ob(bh.metric, aSim, gsb.eccentricity);
            ob.label = gsb.label;
            ob.bodyType = gsb.type;
            ob.isGalaxyBody = true;
            ob.phi = (2.0 * M_PI * i) / preset.numGalaxyBodies; // stagger starting angles evenly so bodies don't all pile up at phi=0 which looks terrible
            ob.trail.clear();
            ob.trail.emplace_back(ob.worldX(), ob.worldY());
            bodies.push_back(std::move(ob));
        }
        activePresetIdx = presetIdx;
        galaxySystemActive = true;
    }

    void clearGalaxySystem() {
        bodies.clear();
        // Restore the default single test body now that we're back in free-orbit mode.
        bodies.emplace_back(bh.metric, 10.0 * bh.metric.M, 0.3);
        galaxySystemActive = false;
    }

    void reset() {
        bh.metric.M = 1.0;
        bodies.clear();
        bodies.emplace_back(bh.metric, 10.0 * bh.metric.M, 0.3);
        galaxySystemActive = false;
        activePresetIdx = -1;
        activeScenario = ResearchScenario::None;
        iscoResults.clear();
        photonSphereResults.clear();
        selectedBodyIdx = 0;
        // note: this does NOT reset the exportName. you're welcome to argue that it should,
        // but I've already accidentally lost data once by having reset() wipe a name I forgot to save.
        // so it stays. fight me.
    }

    /*--------- Lensing analysis ---------*/
    void buildLensingData() {
        lensingData.criticalImpactParam = bh.metric.criticalImpact();
        lensingData.deflectionTable.clear();
        lensingData.causticPoints.clear();

        for (const auto& p : photons) {
            PhotonDeflection d;
            d.impactParameter = p.impactParameter;
            d.captured = p.captured;
            d.deflectionAngle = p.deflectionAngle;
            d.deflectionDeg   = p.deflectionDeg;
            lensingData.deflectionTable.push_back(d);
        }

        // Caustic detection: find where adjacent non-captured rays cross
        // by checking if deflection angles cause path inversion
        for (size_t i = 1; i < lensingData.deflectionTable.size(); ++i) {
            const auto& prev = lensingData.deflectionTable[i-1];
            const auto& curr = lensingData.deflectionTable[i];
            if (prev.captured || curr.captured) continue;

            // Large deflection gradient → ray convergence (caustic region)
            double db = curr.impactParameter - prev.impactParameter;
            if (std::abs(db) < 1e-15) continue;
            double dDeflect = curr.deflectionAngle - prev.deflectionAngle;
            double gradient = std::abs(dDeflect / db);

            // threshold of 0.5 rad/M is completely empirical — I tuned it until the caustic markers
            // appeared in roughly the right places visually. I cannot rigorously justify this number.
            // if caustics are appearing in weird places, this is probably why.
            if (gradient > 0.5) {
                float avgB = (float)(0.5 * (prev.impactParameter + curr.impactParameter));
                lensingData.causticPoints.emplace_back(avgB, (float)gradient);
            }
        }
    }

    /*--------- ISCO validation test ---------*/
    void startISCOTest() {
        activeScenario = ResearchScenario::ISCOTest;
        // Remove non-galaxy bodies (except first default)
        bodies.clear();
        iscoResults.clear();

        double M = bh.metric.M;
        double testRadii[] = { 5.0 * M, 6.0 * M, 7.0 * M };
        const char* labels[] = { "5M (Unstable)", "6M (Critical)", "7M (Stable)" };

        for (int i = 0; i < 3; ++i) {
            OrbitingBody ob(bh.metric, testRadii[i], 0.0);
            ob.initCircularOrbit(bh.metric, testRadii[i]);
            ob.label = labels[i];
            // tiny radial kick to break the exact circular orbit so we can actually see stability play out.
            // 1e-6 is small enough not to visually disturb the orbit but large enough to diverge over time.
            // the sign matters: inward kick for the unstable case pushes it toward the horizon faster
            ob.vr = 1e-6 * (i == 0 ? -1.0 : 1.0);
            bodies.push_back(std::move(ob));

            ISCOTestResult res;
            res.testRadius_M = testRadii[i] / M;
            res.classification = bh.metric.stabilityClassification(testRadii[i]);
            iscoResults.push_back(res);
        }
        selectedBodyIdx = 0;
    }

    void updateISCOResults() {
        for (size_t i = 0; i < iscoResults.size() && i < bodies.size(); ++i) {
            double M = bh.metric.M;
            iscoResults[i].captured = bodies[i].captured;
            iscoResults[i].radiusDrift = (bodies[i].r - iscoResults[i].testRadius_M * M) / M;
            iscoResults[i].survivalTime = bodies[i].measurement.properTime;
        }
    }

    /*--------- Photon sphere test ---------*/
    void startPhotonSphereTest() {
        activeScenario = ResearchScenario::PhotonSphereTest;
        photonSphereResults.clear();

        double b_crit = bh.metric.criticalImpact();
        // Test photons near the critical impact parameter
        double epsilons[] = { 0.001, 0.01, 0.05, -0.001, -0.01, -0.05 };

        for (double eps : epsilons) {
            double b = b_crit * (1.0 + eps);

            Photon p;
            p.impactParameter = b;
            p.computePath(bh.metric, params.rMaxIntegrate, 0.001);

            PhotonSphereTestResult res;
            res.impactParameter = b;
            res.captured = p.captured;
            res.escaped  = !p.captured;

            if (!p.captured) {
                // Count orbits: total angle / (2π)
                res.stabilityAngle = p.deflectionAngle;
                res.orbitsBeforeDecay = std::abs(p.deflectionAngle) / (2.0 * M_PI);
            } else {
                // For captured photons, estimate from path length
                if (p.path.size() >= 2) {
                    double totalAngle = 0.0;
                    for (size_t i = 1; i < p.path.size(); ++i) {
                        double a1 = std::atan2(p.path[i-1].y, p.path[i-1].x);
                        double a2 = std::atan2(p.path[i].y, p.path[i].x);
                        double da = a2 - a1;
                        if (da > M_PI) da -= 2.0 * M_PI;
                        if (da < -M_PI) da += 2.0 * M_PI;
                        totalAngle += std::abs(da);
                    }
                    res.orbitsBeforeDecay = totalAngle / (2.0 * M_PI);
                    res.stabilityAngle = totalAngle;
                }
            }
            photonSphereResults.push_back(res);
        }
    }

    /*--------- Radial infall test ---------*/
    void startRadialInfall() {
        activeScenario = ResearchScenario::RadialInfall;
        bodies.clear();

        double M = bh.metric.M;
        double r0 = 20.0 * M; // start at 20M because anything closer and the thing hits the horizon
                               // before the user has time to read the data panel. tried 10M once. not great.

        OrbitingBody ob(bh.metric, r0, 0.0);
        ob.initRadialInfall(bh.metric, r0);
        ob.label = "Radial infall from 20M";
        bodies.push_back(std::move(ob));
        selectedBodyIdx = 0;
    }

    /*--------- Pulsar orbital simulation ---------*/
    void startPulsarOrbital() {
        activeScenario = ResearchScenario::PulsarOrbital;
        bodies.clear();

        double M   = bh.metric.M;
        double a0  = 20.0 * M;  // start at 20M — comfortably outside ISCO at 6M
        double e0  = 0.3;       // moderate eccentricity for visible variation

        OrbitingBody ob(bh.metric, a0, e0);
        ob.label = "Pulsar (1.4 Msun NS)";
        bodies.push_back(std::move(ob));
        selectedBodyIdx = 0;

        // Initialise pulsar physics state
        pulsarState.init(a0, e0);

        // Reset data struct and set normalisation references
        pulsarData = PulsarOrbitalData{};
        pulsarData.a0 = a0 / M;   // in units of M_BH for history normalisation

        double M_NS  = pulsarState.massGeom;
        double M_tot = M + M_NS;
        double T0_g  = pulsarOrbitalPeriodGeom(a0, M_tot);
        pulsarData.T0 = T0_g / units::c_SI;   // initial period in seconds

        double h0 = pulsarGWStrain(a0, e0, M, M_NS);
        pulsarData.h0 = (h0 > 1e-50) ? h0 : 1.0;
    }

    // Adiabatic GW-decay update called every sim frame when scenario is PulsarOrbital.
    void updatePulsar(double dtSim) {
        if (bodies.empty() || !pulsarState.active) return;

        auto& body = bodies[0];
        if (body.captured) { pulsarData.disrupted = true; return; }

        double M_BH  = bh.metric.M;
        double M_NS  = pulsarState.massGeom;
        double R_NS  = pulsarState.radiusM;
        double M_tot = M_BH + M_NS;

        // ── Adiabatic GW + magnetic orbital decay ────────────────────────────────
        // da/dt and de/dt are in geometric units (m per m of coordinate time).
        // min(1, M_BH/M_NS) suppresses GW decay for the unphysical M=1 default case.
        double mass_ratio_clamp = std::min(1.0, M_BH / std::max(M_NS, 1e-30));
        double vis = PulsarState::GW_VIS_SCALE * mass_ratio_clamp;

        double da_gw  = pulsarDadt(pulsarState.a, pulsarState.e, M_BH, M_NS) * dtSim * vis;
        double de_raw = pulsarDedt(pulsarState.a, pulsarState.e, M_BH, M_NS) * dtSim * vis;

        // Magnetic orbital decay (unipolar inductor, amplified for data-panel display).
        // Cap at 5 % of the GW contribution to avoid dominating the inspiral.
        double da_mag_real = pulsarMagDadt(pulsarState.a, M_BH, M_NS,
                                           pulsarState.magField, R_NS);
        double da_mag_vis  = da_mag_real * dtSim * PulsarState::MAG_VIS_SCALE;
        double da_mag_cap  = std::max(da_mag_vis, -std::abs(da_gw) * 0.05);
        double da_raw      = da_gw + da_mag_cap;

        // Safety clamp: limit to max 3 % of current a per step to prevent
        // numerical blow-up in the M=1 default-unit case.
        double max_da = -0.03 * pulsarState.a;
        double da     = std::max(da_raw, max_da);
        double de     = de_raw;  // de/dt is always negative (orbit circularises)

        pulsarState.a = std::max(pulsarState.a + da, bh.metric.isco() * 1.05);
        pulsarState.e = std::clamp(pulsarState.e + de, 0.0, 0.98);

        // Update body's conserved E/L to match the decayed (a, e).
        // The body continues from its current (r, phi, vr) but now evolves
        // on the new, slightly smaller Schwarzschild orbit.
        double r_peri = pulsarState.a * (1.0 - pulsarState.e);
        double r_apo  = pulsarState.a * (1.0 + pulsarState.e);
        r_peri = std::max(r_peri, bh.metric.isco() + 0.05 * M_BH);
        r_apo  = std::max(r_apo,  r_peri + 0.05 * M_BH);
        auto newEL  = bh.metric.boundOrbitEL(r_peri, r_apo);
        body.E       = newEL.E;
        body.L       = newEL.L;
        body.nominalA   = pulsarState.a;
        body.nominalEcc = pulsarState.e;

        // Merge detection
        if (pulsarState.a <= bh.metric.isco() * 1.06) {
            pulsarData.disrupted = true;
        }

        // ── Category 1: Pulsar Timing ──────────────────────────────────────────
        double r   = body.r;
        double phi = body.phi;

        pulsarData.shapiroDelay_us = pulsarShapiroDelayUs(phi, M_BH);
        pulsarData.shapiroMax_us   = pulsarShapiroMaxUs(M_BH);
        pulsarData.gravRedshift    = pulsarGravRedshift(r, M_BH);
        pulsarData.dopplerFactor   = pulsarDopplerFactor(r, phi, M_BH);

        // EMA baseline for timing residual (alpha = 0.015 → slow follower)
        const double ALPHA = 0.015;
        pulsarData.shapiro_ema_us = (1.0 - ALPHA) * pulsarData.shapiro_ema_us
                                  + ALPHA * pulsarData.shapiroDelay_us;
        pulsarData.timingResidual_us = pulsarData.shapiroDelay_us
                                     - pulsarData.shapiro_ema_us;

        // ── Category 2: Orbital Evolution ──────────────────────────────────────
        pulsarData.semiMajorAxis_M = pulsarState.a / M_BH;
        pulsarData.eccentricity    = pulsarState.e;
        double T_geom              = pulsarOrbitalPeriodGeom(pulsarState.a, M_tot);
        pulsarData.period_s        = T_geom / units::c_SI;
        pulsarData.dadt_real       = pulsarDadt(pulsarState.a, pulsarState.e, M_BH, M_NS);
        pulsarData.timeToMerger_s  = pulsarMergerTimeS(pulsarState.a, pulsarState.e,
                                                        M_BH, M_NS);

        // ── Category 3: Gravitational Waves ────────────────────────────────────
        pulsarData.gwStrain        = pulsarGWStrain(pulsarState.a, pulsarState.e,
                                                     M_BH, M_NS);
        pulsarData.gwFreq_Hz       = pulsarGWFreqHz(pulsarState.a, M_tot);
        pulsarData.gwPower_ergs    = pulsarGWPowerErgs(pulsarState.a, pulsarState.e,
                                                        M_BH, M_NS);
        pulsarData.chirpMass_solar = pulsarChirpMassSolar(M_BH, M_NS);

        // ── Category 4: Structural Integrity ────────────────────────────────────
        double r_roche             = pulsarRocheLimit(M_BH, M_NS, R_NS);
        pulsarData.rocheLimit_M    = r_roche / M_BH;
        pulsarData.currentRadius_M = r / M_BH;
        pulsarData.tidalForce_geom = 2.0 * M_BH / (r * r * r);
        pulsarData.bindingEnergy_J = pulsarBindingEnergyJ(M_NS, R_NS);
        pulsarData.rocheMargin     = (r_roche > 1e-30) ? r / r_roche : 1e9;
        pulsarData.rocheViolated   = (r < r_roche && !body.captured);
        if (pulsarData.rocheViolated) pulsarData.disrupted = true;

        // ── Category 5: Spin & Magnetosphere ──────────────────────────────────
        // Spin-down: target period to double in ~100 orbital periods, which gives a
        // visible slowdown of the jets over the same timescale as the GW inspiral.
        {
            double T_orb_vis = pulsarOrbitalPeriodGeom(pulsarState.a, M_tot);
            double frac      = 0.005 * (dtSim / std::max(T_orb_vis, 1e-30));
            pulsarState.spinOmega *= (1.0 - frac);
            double omega0          = 2.0 * M_PI / pulsarState.spinPeriod;
            pulsarState.spinOmega  = std::max(pulsarState.spinOmega, 0.02 * omega0);
            pulsarState.lightCylRadius_m = pulsarLightCylRadius_m(pulsarState.spinOmega);

            // Visual spin phase: ~1 rotation per orbit so jets sweep visibly but
            // don't blur into a smear. Slows naturally as spinOmega decays.
            double vis_rate = (2.0 * M_PI / std::max(T_orb_vis, 1e-30))
                            * (pulsarState.spinOmega / omega0);
            pulsarState.spinPhase += vis_rate * dtSim;

            // Geodetic precession (visually amplified to match inspiral timescale).
            double prec_geom = 1.5 * std::pow(M_BH, 1.5) / std::pow(pulsarState.a, 2.5);
            pulsarState.precPhase += prec_geom * dtSim * PulsarState::GW_VIS_SCALE;

            // Real (unscaled) observables for the data panel
            double Omega_real = 2.0 * M_PI / pulsarState.spinPeriod;
            double dO_dt_real = pulsarSpinDownRate(pulsarState.magField, R_NS, Omega_real);
            double R_LC_m     = pulsarState.lightCylRadius_m;

            pulsarData.spinPeriod_s     = 2.0 * M_PI / pulsarState.spinOmega;
            pulsarData.Pdot             = -(2.0 * M_PI / (Omega_real * Omega_real)) * dO_dt_real;
            pulsarData.spinDownLum_ergs = pulsarSpinDownLum_ergs(pulsarState.magField, R_NS,
                                                                   pulsarState.spinOmega);
            pulsarData.lightCylRadius_M = R_LC_m / std::max(M_BH, 1e-30);
            pulsarData.inLightCylinder  = (r < R_LC_m);
            pulsarData.geodeticRate_deg =
                pulsarGeodeticPrecRate_radps(pulsarState.a, M_BH)
                * (180.0 / M_PI) * (365.25 * 86400.0);
        }

        // ── Category 6: Magnetic Coupling (unipolar inductor) ─────────────────
        {
            double P_uni_W         = pulsarUnipolePower_W(r, M_BH, pulsarState.magField, R_NS);
            pulsarData.magPower_ergs = P_uni_W * 1.0e7;  // W → erg/s
            pulsarData.magDadt_real  = pulsarMagDadt(pulsarState.a, M_BH, M_NS,
                                                      pulsarState.magField, R_NS);
        }

        // ── History ─────────────────────────────────────────────────────────────
        pulsarData.coordTime += dtSim;
        pulsarData.pushHistory();
    }

    /*--------- Tidal disruption demo ---------*/
    void startTidalDisruption() {
        activeScenario  = ResearchScenario::TidalDisruption;
        tidalDisrupted  = false;
        bodies.clear();

        double M = bh.metric.M;

        // eccentricity 0.85 — this is deliberately extreme so the body actually reaches tidal disruption range.
        // if you set it lower the body just orbits happily and nothing dramatic happens, which defeats the point.
        // if it crosses ISCO it'll plunge in, which is also fine — that's kind of what tidal disruption events do.
        OrbitingBody ob(bh.metric, 8.0 * M, 0.85);
        ob.label = "Tidal test body";
        bodies.push_back(std::move(ob));
        selectedBodyIdx = 0;
    }

    void updateTidalDisruption(double /*dtSim*/) {
        if (bodies.empty() || tidalDisrupted) return;
        const auto& body = bodies[0];
        if (body.captured) {
            tidalDisrupted = true;
            return;
        }
        double r_tidal = tidalRadiusM * bh.metric.M;
        if (body.r < r_tidal) {
            tidalDisrupted = true;
            triggerTidalDisruption(0);
        }
    }

    // Spawn a visible tidal disruption event for the body at bodyIdx, then
    // remove it from the bodies list.
    // Velocities are in geometric units / wall-second, scaled by M so that
    // visual drift rate (≈ M * 0.3 / ppm) stays ~18 px/s across BH masses
    // because default ppm = 60/M.
    void triggerTidalDisruption(int bodyIdx) {
        if (bodyIdx < 0 || bodyIdx >= (int)bodies.size()) return;
        const auto& body = bodies[bodyIdx];

        const double M  = bh.metric.M;
        const double cx = body.r * std::cos(body.phi);
        const double cy = body.r * std::sin(body.phi);

        tidalEvent.active     = true;
        tidalEvent.flashTimer = TidalEvent::FLASH_DURATION;
        tidalEvent.eventX     = cx;
        tidalEvent.eventY     = cy;
        tidalEvent.particles.clear();

        // Prograde tangential direction and inward radial direction
        const double tang_x   = -std::sin(body.phi);
        const double tang_y   =  std::cos(body.phi);
        const double inward_x = -std::cos(body.phi);
        const double inward_y = -std::sin(body.phi);

        // Simple deterministic LCG seeded from the disruption position
        uint32_t rng = 0x12345678u ^ static_cast<uint32_t>(std::abs(cx + cy) * 1000.0);
        auto rand01 = [&]() -> double {
            rng = rng * 1664525u + 1013904223u;
            return (rng >> 8) / double(1 << 24);
        };

        const double V_OUT  = 0.30 * M;  // geometric units / wall-second
        const double V_FALL = 0.14 * M;

        // 20 outward debris fragments
        for (int i = 0; i < 20; ++i) {
            const double angle = 2.0 * M_PI * rand01();
            const double speed = V_OUT * (0.7 + 0.6 * rand01());
            TidalDebrisParticle p;
            p.x   = cx;  p.y   = cy;
            p.vx  = std::cos(angle) * speed + tang_x * V_OUT * 0.4;
            p.vy  = std::sin(angle) * speed + tang_y * V_OUT * 0.4;
            p.maxLifetime = TidalEvent::DEBRIS_LIFETIME * (0.6 + 0.4 * rand01());
            p.lifetime    = p.maxLifetime;
            p.size        = 2.5f + (float)(1.5 * rand01());
            p.isFallback  = false;
            tidalEvent.particles.push_back(p);
        }

        // 10 inward fallback gas-stream particles
        for (int i = 0; i < 10; ++i) {
            const double jitter = (rand01() - 0.5) * 0.8;
            const double speed  = V_FALL * (0.8 + 0.4 * rand01());
            TidalDebrisParticle p;
            p.x   = cx + tang_x * M * 0.5 * (rand01() - 0.5);
            p.y   = cy + tang_y * M * 0.5 * (rand01() - 0.5);
            p.vx  = (inward_x + tang_x * jitter) * speed;
            p.vy  = (inward_y + tang_y * jitter) * speed;
            p.maxLifetime = TidalEvent::STREAM_LIFETIME * (0.7 + 0.3 * rand01());
            p.lifetime    = p.maxLifetime;
            p.size        = 2.0f + (float)(rand01());
            p.isFallback  = true;
            tidalEvent.particles.push_back(p);
        }

        bodies.erase(bodies.begin() + bodyIdx);
    }

    // Advance the tidal particle system using wall-clock dt.
    void updateTidalEvent(double dt) {
        if (!tidalEvent.active) return;

        if (tidalEvent.flashTimer > 0.0)
            tidalEvent.flashTimer = std::max(0.0, tidalEvent.flashTimer - dt);

        for (auto& p : tidalEvent.particles) {
            p.x += p.vx * dt;
            p.y += p.vy * dt;
            p.lifetime -= dt;
        }

        tidalEvent.particles.erase(
            std::remove_if(tidalEvent.particles.begin(), tidalEvent.particles.end(),
                [](const TidalDebrisParticle& p){ return p.lifetime <= 0.0; }),
            tidalEvent.particles.end());

        if (tidalEvent.particles.empty() && tidalEvent.flashTimer <= 0.0)
            tidalEvent.active = false;
    }

    /*--------- Data panel info builder ---------*/
    DataPanelInfo buildDataPanel() const {
        DataPanelInfo info;
        if (bodies.empty() || selectedBodyIdx >= (int)bodies.size())
            return info;

        const auto& body = bodies[selectedBodyIdx];
        double M = bh.metric.M;

        info.radius_M = body.r / M;
        info.velocity_c = body.measurement.orbitalVelocity;
        info.timeDilation = body.measurement.timeDilation;
        info.redshift_z = body.measurement.redshift - 1.0;
        info.energy_E = body.E;
        info.angularMomentum_L = body.L;
        info.properTime = body.measurement.properTime;
        info.coordinateTime = body.measurement.coordinateTime;
        info.escapeVelocity_c = bh.metric.escapeVelocity(body.r);
        info.isBound = (body.E < 1.0);
        info.stabilityClass = bh.metric.stabilityClassification(body.r);

        // Tidal
        info.tidalForce = bh.metric.tidalForce(body.r);
        info.tidalStress1m = bh.metric.tidalStress(body.r, 1.0);

        // Conservation
        info.energyDrift = body.conservation.relativeDrift();
        info.maxEnergyDrift = body.conservation.maxAbsDrift;

        // Numerical error
        info.rk4LastError = body.errorTracker.lastError;
        info.rk4MaxError  = body.errorTracker.maxError;
        info.rk4AvgError  = body.errorTracker.avgError;

        // Precession
        info.precessionDeg = body.precession.lastPrecessionDeg();
        info.theoreticalPrecessionDeg =
            bh.metric.theoreticalPrecession(body.nominalA, body.nominalEcc) * 180.0 / M_PI;
        info.orbitsCompleted = body.precession.orbitsCompleted();

        // Bondi
        double cs = Schwarzschild::soundSpeedFromTemp(gasTemperatureK);
        info.gasTemperatureK = gasTemperatureK;
        info.soundSpeed_c = cs;
        info.bondiRadius_M = bh.metric.bondiRadius(cs) / M;
        info.gasDensity_kgm3   = gasDensity_kgm3;
        info.mdotBondi_MsunYr  = accretion.last_mdotBondi_kgs
                                  * BondiAccretionTracker::YEAR_S
                                  / BondiAccretionTracker::MSUN_KG;
        info.mdotBHL_MsunYr    = accretion.last_mdotBHL_MsunYr;
        info.Lbol_W            = accretion.last_Lbol_W;
        info.Lbol_Lsun         = accretion.last_Lbol_Lsun;

        return info;
    }

    /*--------- Format data panel as string ---------*/
    std::string formatDataPanel() const {
        auto info = buildDataPanel();
        auto sci = [](double v, int prec = 3) -> std::string {
            std::ostringstream os;
            if (std::fabs(v) >= 1e5 || (std::fabs(v) > 0 && std::fabs(v) < 1e-3))
                os << std::scientific << std::setprecision(prec) << v;
            else
                os << std::fixed << std::setprecision(prec + 1) << v;
            return os.str();
        };

        std::ostringstream ss;
        ss << "=== RESEARCH DATA ===\n";
        ss << "Body [" << selectedBodyIdx << "/"
           << (int)bodies.size()-1 << "]";
        if (selectedBodyIdx < (int)bodies.size() && !bodies[selectedBodyIdx].label.empty())
            ss << " " << bodies[selectedBodyIdx].label;
        ss << "\n";

        ss << "--- STATE ---\n";
        ss << "Radius:   " << sci(info.radius_M) << "M\n";
        ss << "Velocity: " << sci(info.velocity_c) << "c\n";
        ss << "Escape v: " << sci(info.escapeVelocity_c) << "c\n";
        ss << "Bound:    " << (info.isBound ? "YES" : "NO") << "\n";
        ss << "Stability:" << info.stabilityClass << "\n";

        ss << "--- RELATIVISTIC ---\n";
        ss << "dt/dtau:  " << sci(info.timeDilation) << "\n";
        ss << "Proper t: " << sci(info.properTime) << "\n";
        ss << "Coord t:  " << sci(info.coordinateTime) << "\n";
        ss << "z:        " << sci(info.redshift_z) << "\n";

        ss << "--- CONSERVED ---\n";
        ss << "E:        " << sci(info.energy_E) << "\n";
        ss << "L:        " << sci(info.angularMomentum_L) << "\n";
        ss << "E drift:  " << sci(info.energyDrift * 100.0) << "%\n";
        ss << "Max drift:" << sci(info.maxEnergyDrift * 100.0) << "%\n";

        ss << "--- PRECESSION ---\n";
        ss << "Orbits:   " << info.orbitsCompleted << "\n";
        ss << "Last dPhi:" << sci(info.precessionDeg) << " deg/orb\n";
        ss << "Theory:   " << sci(info.theoreticalPrecessionDeg) << " deg/orb\n";

        ss << "--- TIDAL ---\n";
        ss << "F_tidal:  " << sci(info.tidalForce) << " /M^2\n";
        ss << "Stress 1m:" << sci(info.tidalStress1m) << "\n";

        ss << "--- BONDI ---\n";
        ss << "T_gas:    " << sci(info.gasTemperatureK, 1) << " K\n";
        ss << "c_s:      " << sci(info.soundSpeed_c) << " c\n";
        ss << "r_Bondi:  " << sci(info.bondiRadius_M) << "M\n";
        ss << "rho_gas:  " << sci(info.gasDensity_kgm3) << " kg/m^3\n";

        ss << "--- ACCRETION ---\n";
        ss << "Mdot_B:   " << sci(info.mdotBondi_MsunYr) << " Msun/yr\n";
        ss << "Mdot_BHL: " << sci(info.mdotBHL_MsunYr)   << " Msun/yr\n";
        ss << "L_bol:    " << sci(info.Lbol_W)    << " W\n";
        ss << "L_bol:    " << sci(info.Lbol_Lsun) << " Lsun\n";

        ss << "--- NUMERICAL ---\n";
        ss << "RK4 last: " << sci(info.rk4LastError) << "\n";
        ss << "RK4 max:  " << sci(info.rk4MaxError) << "\n";
        ss << "RK4 avg:  " << sci(info.rk4AvgError) << "\n";

        // ISCO test results
        if (activeScenario == ResearchScenario::ISCOTest && !iscoResults.empty()) {
            ss << "--- ISCO TEST ---\n";
            for (const auto& r : iscoResults) {
                ss << sci(r.testRadius_M) << "M: "
                   << r.classification
                   << " drift=" << sci(r.radiusDrift)
                   << (r.captured ? " CAPTURED" : "") << "\n";
            }
        }

        // Photon sphere test
        if (activeScenario == ResearchScenario::PhotonSphereTest && !photonSphereResults.empty()) {
            ss << "--- PHOTON SPHERE ---\n";
            ss << "b_crit=" << sci(bh.metric.criticalImpact()) << "\n";
            for (const auto& r : photonSphereResults) {
                ss << "b=" << sci(r.impactParameter)
                   << " orbits=" << sci(r.orbitsBeforeDecay, 1)
                   << (r.escaped ? " ESC" : " CAP") << "\n";
            }
        }

        // Lensing
        if (!lensingData.deflectionTable.empty()) {
            ss << "--- LENSING ---\n";
            ss << "b_crit: " << sci(lensingData.criticalImpactParam) << "\n";
            ss << "Caustic pts: " << lensingData.causticPoints.size() << "\n";
        }

        // Pulsar orbital scenario
        if (activeScenario == ResearchScenario::PulsarOrbital) {
            const auto& pd = pulsarData;
            ss << "--- PULSAR ORBITAL ---\n";
            ss << "NS: " << std::fixed << std::setprecision(2)
               << pulsarState.massSolar << " Msun";
            if (pd.disrupted) ss << "  !! DISRUPTED !!";
            ss << "\n";

            ss << "--- 1. TIMING ---\n";
            ss << "Shapiro:  " << sci(pd.shapiroDelay_us) << " us\n";
            ss << "Max Shap: " << sci(pd.shapiroMax_us)   << " us\n";
            ss << "Residual: " << sci(pd.timingResidual_us) << " us\n";
            ss << "nu_obs/0: " << sci(pd.gravRedshift)    << "\n";
            ss << "Doppler:  " << sci(pd.dopplerFactor)   << "\n";

            ss << "--- 2. ORBIT DECAY ---\n";
            ss << "a/M_BH:  " << sci(pd.semiMajorAxis_M) << "\n";
            ss << "e:       " << sci(pd.eccentricity)     << "\n";
            ss << "Period:  " << sci(pd.period_s)         << " s\n";
            ss << "da/dt:   " << sci(pd.dadt_real)        << " m/m\n";

            // Merger time: display in years if large, seconds if small
            double t_merge = pd.timeToMerger_s;
            if (t_merge > 3.156e7)
                ss << "Merger:  ~" << sci(t_merge / 3.156e7) << " yr\n";
            else
                ss << "Merger:  ~" << sci(t_merge) << " s\n";

            ss << "--- 3. GRAV WAVES ---\n";
            ss << "h(1kpc): " << sci(pd.gwStrain)        << "\n";
            ss << "f_GW:    " << sci(pd.gwFreq_Hz)       << " Hz\n";
            ss << "P_GW:    " << sci(pd.gwPower_ergs)    << " erg/s\n";
            ss << "M_chirp: " << sci(pd.chirpMass_solar) << " Msun\n";

            ss << "--- 4. INTEGRITY ---\n";
            ss << "r_Roche: " << sci(pd.rocheLimit_M)    << "M\n";
            ss << "r_curr:  " << sci(pd.currentRadius_M) << "M\n";
            ss << "r/r_R:   " << sci(pd.rocheMargin)     << "\n";
            ss << "F_tidal: " << sci(pd.tidalForce_geom) << " /M^2\n";
            ss << "E_bind:  " << sci(pd.bindingEnergy_J) << " J\n";
            if (pd.rocheViolated) ss << "!! ROCHE LIMIT BREACHED !!\n";

            ss << "--- 5. SPIN & MAGNETOSPHERE ---\n";
            ss << "P_spin: " << sci(pd.spinPeriod_s)     << " s\n";
            ss << "P_dot:  " << sci(pd.Pdot)             << " s/s\n";
            ss << "L_sd:   " << sci(pd.spinDownLum_ergs) << " erg/s\n";
            ss << "R_LC:   " << sci(pd.lightCylRadius_M) << " M\n";
            ss << "In LC:  " << (pd.inLightCylinder ? "YES (magnetically coupled)" : "NO") << "\n";
            ss << "Prec:   " << sci(pd.geodeticRate_deg) << " deg/yr\n";

            ss << "--- 6. MAGNETIC COUPLING ---\n";
            ss << "P_uni:  " << sci(pd.magPower_ergs)  << " erg/s\n";
            ss << "da/dt_B:" << sci(pd.magDadt_real)   << " m/m\n";
            ss << "(BH unipolar inductor: field lines connect NS to BH)\n";
        }

        ss << "\nH: cycle body  X: export data";
        return ss.str();
    }

    /*--------- Helpers ---------*/

    // Strip characters that are invalid or problematic in directory names.
    static std::string sanitizeName(const std::string& s) {
        std::string out;
        out.reserve(s.size());
        for (unsigned char c : s) {
            if (std::isalnum(c) || c == '_' || c == '-' || c == '.' || c == ' ')
                out += static_cast<char>(c);
            else
                out += '_';
        }
        // Trim leading/trailing whitespace
        size_t start = out.find_first_not_of(' ');
        if (start == std::string::npos) return "untitled";
        size_t end = out.find_last_not_of(' ');
        return out.substr(start, end - start + 1);
    }

    // Returns the workspace-level export directory, or "" on failure.
    std::string makeExportDir() const {
        const char *home = getenv("HOME");
        if (!home) home = "/tmp";
        std::string safeName = sanitizeName(exportName.empty() ? "untitled_data" : exportName);
        std::string dir;
#ifdef __APPLE__
        dir = std::string(home) + "/Library/Application Support/Aetherion/exports/" + safeName + "/";
#elif defined(__linux__)
        dir = std::string(home) + "/.local/share/Aetherion/exports/" + safeName + "/";
#else
        dir = std::string(home) + "/Aetherion/exports/" + safeName + "/";
#endif
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        if (ec) return "";   // Signal failure to callers
        return dir;
    }

    // Returns a per-body subdirectory (body_<idx>[_<label>]/) so multiple
    // galaxy bodies never overwrite each other.  Returns "" on failure.
    std::string makeBodyExportDir() const {
        std::string base = makeExportDir();
        if (base.empty()) return "";

        std::string subdir = "body_" + std::to_string(selectedBodyIdx);
        if (selectedBodyIdx < (int)bodies.size() && !bodies[selectedBodyIdx].label.empty())
            subdir += "_" + sanitizeName(bodies[selectedBodyIdx].label);

        std::string dir = base + subdir + "/";
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        if (ec) return "";
        return dir;
    }

    // Abbreviates the home directory as "~" so messages fit in the HUD.
    static std::string abbreviateHome(const std::string& path) {
        const char *home = getenv("HOME");
        if (home && path.substr(0, std::strlen(home)) == home)
            return "~" + path.substr(std::strlen(home));
        return path;
    }

    /*--------- CSV export ---------*/
    std::string exportAllCSV() const {
        std::string msg;
        std::string displayDir;

        if (!bodies.empty() && selectedBodyIdx < (int)bodies.size()) {
            if (bodies[selectedBodyIdx].trail.size() <= 1)
                return "No orbit data yet — let the simulation run first";
            std::string dir = makeBodyExportDir();
            if (dir.empty()) return "Export failed: could not create directory";
            displayDir = dir;
            const auto& body = bodies[selectedBodyIdx];
            if (CSVExport::exportOrbitData(dir + "orbit_data.csv", body.trail))
                msg += "orbit_data.csv ";
            if (CSVExport::exportConservationHistory(dir + "conservation.csv", body.conservation))
                msg += "conservation.csv ";
            if (CSVExport::exportPrecessionData(dir + "precession.csv", body.precession))
                msg += "precession.csv ";
        }
        if (!lensingData.deflectionTable.empty()) {
            std::string base = makeExportDir();
            if (!base.empty()) {
                if (CSVExport::exportDeflectionTable(base + "deflection.csv", lensingData.deflectionTable))
                    msg += "deflection.csv ";
                if (displayDir.empty()) displayDir = base;
            }
        }
        if (msg.empty()) return "No data to export";
        return "Saved: " + msg + "| " + abbreviateHome(displayDir);
    }

    /*--------- FITS export ---------*/
    std::string exportAllFITS() const {
        std::string msg;
        std::string displayDir;

        // Build the simulation metadata block that gets baked into every
        // FITS primary HDU we write below. Research-grade FITS requires the
        // simulation parameters to live in the file header itself so that
        // downstream tools can recover units and physical context without
        // any out-of-band sidecar metadata.
        FITSExport::SimulationMetadata meta;
        meta.M_BH    = (activePresetIdx >= 0 && activePresetIdx < NUM_BH2D_PRESETS)
                         ? BH2D_PRESETS[activePresetIdx].massSolar
                         : 0.0;
        meta.c_s     = Schwarzschild::soundSpeedFromTemp(gasTemperatureK);
        meta.r_Bondi = (meta.c_s > 0.0)
                         ? bh.metric.bondiRadius(meta.c_s) / bh.metric.M
                         : 0.0;
        meta.beta    = (!bodies.empty() && selectedBodyIdx >= 0 && selectedBodyIdx < (int)bodies.size())
                         ? bodies[selectedBodyIdx].measurement.orbitalVelocity
                         : 0.0;
        meta.valid   = true;

        if (!bodies.empty() && selectedBodyIdx < (int)bodies.size()) {
            if (bodies[selectedBodyIdx].trail.size() <= 1)
                return "No orbit data yet — let the simulation run first";
            std::string dir = makeBodyExportDir();
            if (dir.empty()) return "Export failed: could not create directory";
            displayDir = dir;
            const auto& body = bodies[selectedBodyIdx];
            if (FITSExport::exportOrbitData(dir + "orbit_data.fits", body.trail, meta))
                msg += "orbit_data.fits ";
            if (FITSExport::exportConservationHistory(dir + "conservation.fits", body.conservation, meta))
                msg += "conservation.fits ";
            if (FITSExport::exportPrecessionData(dir + "precession.fits", body.precession, meta))
                msg += "precession.fits ";
        }
        // Accretion history is global (not per-body) — export it whenever we
        // have any sampled entries, regardless of which body is selected.
        if (!accretion.history.empty()) {
            std::string base = makeExportDir();
            if (!base.empty()) {
                if (FITSExport::exportAccretionHistory(base + "accretion.fits", accretion, meta))
                    msg += "accretion.fits ";
                if (displayDir.empty()) displayDir = base;
            }
        }
        // Gas radial profile + SED snapshots — recomputed at export time from
        // current sim state so the FITS file always reflects the latest values.
        {
            const double M_BH_solar = currentMassSolar();
            const double cs_over_c  = Schwarzschild::soundSpeedFromTemp(gasTemperatureK);
            const double cs_ms      = cs_over_c * 2.99792458e8;
            if (M_BH_solar > 0.0 && cs_ms > 0.0 && gasDensity_kgm3 > 0.0) {
                GasRadialProfile prof;
                prof.snapshot(M_BH_solar, gasDensity_kgm3, cs_ms);
                SEDSnapshot sedSnap;
                // Use the latest BHL rate (orbit-modulated). Fall back to pure
                // Bondi if the orbit hasn't accumulated yet.
                const double mdot = (accretion.last_mdotBHL_kgs > 0.0)
                                  ? accretion.last_mdotBHL_kgs
                                  : accretion.last_mdotBondi_kgs;
                if (mdot > 0.0) sedSnap.snapshot(M_BH_solar, mdot);

                if (prof.valid || sedSnap.valid) {
                    std::string base = makeExportDir();
                    if (!base.empty()) {
                        if (prof.valid &&
                            FITSExport::exportGasRadialProfile(base + "gas_profile.fits", prof, meta))
                            msg += "gas_profile.fits ";
                        if (sedSnap.valid &&
                            FITSExport::exportSED(base + "sed.fits", sedSnap, meta))
                            msg += "sed.fits ";
                        if (displayDir.empty()) displayDir = base;
                    }
                }
            }
        }
        if (!lensingData.deflectionTable.empty()) {
            std::string base = makeExportDir();
            if (!base.empty()) {
                if (FITSExport::exportDeflectionTable(base + "deflection.fits", lensingData.deflectionTable, meta))
                    msg += "deflection.fits ";
                if (displayDir.empty()) displayDir = base;
            }
        }
        if (msg.empty()) return "No data to export";
        return "Saved: " + msg + "| " + abbreviateHome(displayDir);
    }

    /*--------- Binary export ---------*/
    std::string exportAllBinary() const {
        std::string msg;
        std::string displayDir;

        if (!bodies.empty() && selectedBodyIdx < (int)bodies.size()) {
            if (bodies[selectedBodyIdx].trail.size() <= 1)
                return "No orbit data yet — let the simulation run first";
            std::string dir = makeBodyExportDir();
            if (dir.empty()) return "Export failed: could not create directory";
            displayDir = dir;
            const auto& body = bodies[selectedBodyIdx];
            // Orbit trail: count + (x,y) pairs as doubles
            {
                std::ofstream f(dir + "orbit_data.bin", std::ios::binary);
                if (f) {
                    uint64_t n = body.trail.size();
                    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
                    for (const auto& [x, y] : body.trail) {
                        f.write(reinterpret_cast<const char*>(&x), sizeof(x));
                        f.write(reinterpret_cast<const char*>(&y), sizeof(y));
                    }
                    msg += "orbit_data.bin ";
                }
            }
            // Conservation drift: count + (time, drift) pairs
            {
                std::ofstream f(dir + "conservation.bin", std::ios::binary);
                if (f) {
                    uint64_t n = body.conservation.driftHistory.size();
                    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
                    for (const auto& [t, d] : body.conservation.driftHistory) {
                        f.write(reinterpret_cast<const char*>(&t), sizeof(t));
                        f.write(reinterpret_cast<const char*>(&d), sizeof(d));
                    }
                    msg += "conservation.bin ";
                }
            }
            // Precession: count + precession_per_orbit doubles
            {
                std::ofstream f(dir + "precession.bin", std::ios::binary);
                if (f) {
                    uint64_t n = body.precession.precessionPerOrbit.size();
                    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
                    for (double p : body.precession.precessionPerOrbit)
                        f.write(reinterpret_cast<const char*>(&p), sizeof(p));
                    msg += "precession.bin ";
                }
            }
        }
        if (!lensingData.deflectionTable.empty()) {
            std::string base = makeExportDir();
            if (!base.empty()) {
                std::ofstream f(base + "deflection.bin", std::ios::binary);
                if (f) {
                    uint64_t n = lensingData.deflectionTable.size();
                    f.write(reinterpret_cast<const char*>(&n), sizeof(n));
                    for (const auto& d : lensingData.deflectionTable) {
                        f.write(reinterpret_cast<const char*>(&d.impactParameter),  sizeof(double));
                        f.write(reinterpret_cast<const char*>(&d.deflectionAngle),  sizeof(double));
                        f.write(reinterpret_cast<const char*>(&d.deflectionDeg),    sizeof(double));
                        uint8_t cap = d.captured ? 1 : 0;
                        f.write(reinterpret_cast<const char*>(&cap), sizeof(cap));
                    }
                    msg += "deflection.bin ";
                    if (displayDir.empty()) displayDir = base;
                }
            }
        }
        if (msg.empty()) return "No data to export";
        return "Saved: " + msg + "| " + abbreviateHome(displayDir);
    }

    // Select next body for data panel
    void cycleSelectedBody() {
        if (bodies.empty()) return;
        selectedBodyIdx = (selectedBodyIdx + 1) % (int)bodies.size();
    }
};
