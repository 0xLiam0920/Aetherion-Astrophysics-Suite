#pragma once
#include "black_hole.hpp"
#include "photon.hpp"
#include "orbiting_body.hpp"
#include "research_data.hpp"
#include "csv_export.hpp"
#include "FITS_export.hpp"
#include "../2D-utils/presets_2d.hpp"
#include "../2D-physics/geodesic.hpp"
#include "../2D-physics/kerr.hpp"
#include "../2D-physics/pulsar_orbital.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
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
// (ResearchScenario enum is defined in research_data.hpp so it is also
//  visible to OrbitingBody, which tags itself with the scenario that
//  spawned it.)

// Central simulation manager, orchestrates the black hole world.
class Simulation {
public:
    BlackHole                  bh;
    std::vector<Photon>        photons;
    // Photons fired from on-body hot-spot emitters (see OrbitingBody::emitsPhotons).
    // Kept separate from `photons` so the global lensing/deflection table
    // built by buildLensingData() — which assumes rays come from infinity —
    // is not contaminated by emitter-spawned rays.
    std::vector<Photon>        emittedPhotons;
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

    // ── Disk-annulus photon emitter (batch 2 #2) ────────────────────────────────────────
    // When enabled, the simulation periodically spawns null geodesics from
    // uniformly random points around the inner accretion disk annulus
    // r ∈ [r_ISCO, 12M]. Each ray is tagged with a thermal-disk wavelength
    // (λ ∝ r^(3/4) following T ∝ r^(-3/4)), so the spectrum panel shows the
    // multi-colour disk SED character even from a Monte-Carlo sample.
    bool   diskEmitterEnabled = false;
    double diskEmitterInterval= 0.35;   // wall-seconds between batches
    double diskEmitterTimer   = 0.0;
    int    diskEmitterRays    = 6;      // photons per batch

    // ── Kerr equatorial overlay (batch 2 #4) ───────────────────────────────────────────
    // Spin is in units of M; 0 = Schwarzschild, 0.998 = near-extremal Kerr.
    bool   kerrOverlayEnabled = false;
    double kerrSpin           = 0.7;    // a / M
    std::vector<KerrPathResult> kerrRays;

    // Polarization-tick overlay (batch 2 #3): off by default to keep the
    // canvas uncluttered. Renderer reads this to decide whether to draw
    // little ticks at each emitted-photon endpoint.
    bool   polarizationOverlay = false;

    // ── Tidal disruption event state ─────────────────────────────────────────
    struct TidalDebrisParticle {
        double x, y;           // world position (geometric units, BH at origin)
        double vx, vy;         // velocity (geometric units / wall-second)
        double lifetime;       // remaining lifetime (wall-seconds)
        double maxLifetime;
        float  size;           // base radius in pixels (at reference zoom)
        bool   isFallback;     // true = inward gas-stream particle
        // Specific orbital energy at spawn, in units where v² ~ 2M/r ≈ 1.
        // Used by the renderer to colour the stream by binding (red = bound /
        // returns, blue = unbound / escapes), and lets the simulation skip
        // the gravity solve for particles that have been flagged as escaping
        // straight-line debris.
        double specificEnergy = 0.0;
        bool   unbound        = false;
    };

    // Distinguishes which physical event the shared particle/flash system is
    // representing, so renderers can pick the right label/colour.
    //   TidalDisruption  — a star/WD shredded by the BH (red/orange stream)
    //   MergerShockwave  — GW + EM ring from a BH–BH coalescence (white-purple)
    //   XrayBurst        — NS/pulsar swallowed by the BH (white-blue burst)
    enum class TidalEventKind { TidalDisruption, MergerShockwave, XrayBurst };

    struct TidalEvent {
        bool   active         = false;
        TidalEventKind kind   = TidalEventKind::TidalDisruption;
        double flashTimer     = 0.0;   // seconds remaining for the localised flash
        double eventX         = 0.0;   // disruption world X (geometric units)
        double eventY         = 0.0;   // disruption world Y (geometric units)
        std::vector<TidalDebrisParticle> particles;
        static constexpr double FLASH_DURATION  = 0.45;  // seconds
        static constexpr double DEBRIS_LIFETIME = 5.0;   // seconds
        static constexpr double STREAM_LIFETIME = 18.0;  // seconds (fallback gas — long enough for bound debris to return to periapsis)
    } tidalEvent;

    // Tidal disruption threshold in units of M.
    // Set from the active preset's zones.tidalDisruptionM; default matches the
    // startTidalDisruption() initial orbit so the scenario always triggers cleanly.
    double tidalRadiusM = 8.0;

    // ── Merger event state ───────────────────────────────────────────────────
    // What kind of object is the BH eating? Branches both visuals and the
    // post-merger particle event (BH→GW ring, NS/Pulsar→X-ray burst, Star/WD→
    // tidal-disruption fallback stream).
    enum class MergerSecondaryKind {
        BlackHole = 0,
        NeutronStar = 1,
        Pulsar = 2,
        Star = 3,
        WhiteDwarf = 4
    };

    struct MergerState {
        bool   active          = false;
        bool   completed       = false;
        // User-facing pacing dial. Default keeps the existing ~15s wall-clock
        // pacing; Cinematic is a slow-mo for screen recording; Realtime nudges
        // toward the Peters Mtot-aware speed-up (true wall-clock for SMBH+SMBH
        // would be unwatchable, so this is an approximate factor only).
        enum class TimeScale { Cinematic = 0, Default = 1, Realtime = 2 };
        TimeScale timeScale = TimeScale::Default;
        static double timeScaleFactor(TimeScale t) {
            switch (t) {
                case TimeScale::Cinematic: return 0.35;
                case TimeScale::Default:   return 1.0;
                case TimeScale::Realtime:  return 3.0;
            }
            return 1.0;
        }
        static const char* timeScaleName(TimeScale t) {
            switch (t) {
                case TimeScale::Cinematic: return "Cinematic";
                case TimeScale::Default:   return "Default";
                case TimeScale::Realtime:  return "Realtime";
            }
            return "";
        }
        double r_M             = 40.0;   // separation in units of primary M (scales with M)
        double phi             = 0.0;    // orbital angle (radians)
        double massSolar       = 0.0;    // incoming BH mass in solar masses
        double secondaryMassGeom = 0.0;  // precomputed geometric mass of secondary
        double flashTimer      = 0.0;
        double preMergeMassSolar = 0.0;
        MergerSecondaryKind secondaryKind = MergerSecondaryKind::BlackHole;
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

        // Safety: if some other action (preset switch, scenario reset, …)
        // removed the pulsar body, retire the pulsar state so the HUD/data
        // panel and scenario flag don't get stranded.
        if (pulsarState.active && findPulsarIdx() < 0) {
            pulsarState = PulsarState{};
            pulsarData  = PulsarOrbitalData{};
            if (activeScenario == ResearchScenario::PulsarOrbital)
                activeScenario = ResearchScenario::None;
        }

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
        if (anyBodyWithScenarioTag(ResearchScenario::ISCOTest))
            updateISCOResults();

        // Update pulsar orbital scenario if active
        if (pulsarState.active && anyBodyWithScenarioTag(ResearchScenario::PulsarOrbital))
            updatePulsar(dtSim);

        // Update tidal disruption scenario if active
        if (anyBodyWithScenarioTag(ResearchScenario::TidalDisruption))
            updateTidalDisruption(dtSim);

        // Advance tidal event particle system (always, even after scenario completes)
        updateTidalEvent(dt);

        // Hot-spot photon emitters (per-body pulsed null-geodesic emission)
        updateHotspotEmission(dt);

        // Disk-annulus emitter (random thermal photons from inner disk)
        if (diskEmitterEnabled) updateDiskEmitter(dt);

        // Rebuild the Kerr equatorial overlay when enabled. The ray count is
        // kept low and the integration cheap, so re-sweeping per frame is
        // fine and means the user can drag the spin slider live.
        if (kerrOverlayEnabled) updateKerrOverlay();

        // Advance merger inspiral if active
        if (merger.active)
            updateMerger(dt); // wall-clock dt, not dtSim, so the animation plays at a fixed visual rate
    }

    // Per-body hot-spot photon emission. Each OrbitingBody with
    // emitsPhotons=true periodically fires a fan of N null geodesics outward
    // from its instantaneous position, letting the user watch how light from
    // a compact emitter (a flare, a magnetic reconnection event near the
    // ISCO, an orbiting AGN cloudlet, etc.) is lensed and partially captured
    // by the black hole. The emitted rays are appended to `photons`, but
    // tagged with `fromEmitter = true` so the global lensing analysis
    // (buildLensingData) ignores them.
    void updateHotspotEmission(double dt) {
        for (auto& b : bodies) {
            if (!b.emitsPhotons || b.captured || b.ejected) continue;
            b.emitTimer -= dt;
            if (b.emitTimer > 0.0) continue;
            b.emitTimer = b.emitInterval;
            const int N = std::max(2, b.emitFanCount);
            // Fan α from a small angle away from the radial-inward direction
            // (α = π) to a small angle away from radial-outward (α = 0) so
            // we always have some rays falling in and some escaping.
            for (int i = 0; i < N; ++i) {
                // α ∈ (0.05, π − 0.05) — avoid the exact radial limits, which
                // computeEmissionPath treats as the straight-line case.
                double alpha = 0.05 + (M_PI - 0.10) * (double(i) / double(N - 1));
                Photon p;
                p.wavelength_nm     = b.emitWavelength_nm;
                p.wavelength_obs_nm = b.emitWavelength_nm;
                p.computeEmissionPath(bh.metric, b.r, b.phi, alpha);
                p.fromEmitter = true;
                // Spectrum bin contribution (only for rays that escape).
                if (!p.captured) addSpectrumSample(p.wavelength_obs_nm);
                emittedPhotons.push_back(std::move(p));
            }
        }
        // Cap emitted photons total so long sessions don't grow without
        // bound. Keep the most recent ones.
        constexpr size_t EMITTED_CAP = 240;
        if (emittedPhotons.size() > EMITTED_CAP) {
            emittedPhotons.erase(emittedPhotons.begin(),
                emittedPhotons.begin() + (emittedPhotons.size() - EMITTED_CAP));
        }
    }

    // ── Spectrum binning helper ─────────────────────────────────────────────
    void addSpectrumSample(double lambda_nm) {
        if (!(lambda_nm > 0.0) || !std::isfinite(lambda_nm)) return;
        constexpr double lo = LensingData::SPECTRUM_MIN_NM;
        constexpr double hi = LensingData::SPECTRUM_MAX_NM;
        if (lambda_nm < lo || lambda_nm >= hi) return;
        const int N  = LensingData::SPECTRUM_BINS;
        const int b  = std::min(N - 1, static_cast<int>((lambda_nm - lo) / (hi - lo) * N));
        ++lensingData.spectrumCounts[b];
        ++lensingData.spectrumTotal;
    }

    // ── Disk-annulus photon emitter ─────────────────────────────────────────
    // Spawns a small Monte-Carlo batch of photons from random points around
    // the inner accretion disk r ∈ [r_ISCO, 12M], each tagged with a thermal
    // disk wavelength λ(r) ∝ r^(3/4) (since T ∝ r^(-3/4) and λ_peak ∝ 1/T).
    void updateDiskEmitter(double dt) {
        diskEmitterTimer -= dt;
        if (diskEmitterTimer > 0.0) return;
        diskEmitterTimer = diskEmitterInterval;

        // Cheap LCG so we don't depend on <random> here.
        static uint32_t rng = 0xC0FFEEu;
        auto frand = [&]() -> double {
            rng = rng * 1664525u + 1013904223u;
            return (rng >> 8) / double(1u << 24);   // [0,1)
        };

        const double rIsco = bh.metric.isco();
        const double rOut  = 12.0 * bh.metric.M;
        // λ_peak at ISCO ≈ 480 nm (blue), at 12M ≈ 480·(12/6)^(3/4) ≈ 808 nm (deep red).
        const double lambdaIsco_nm = 480.0;

        for (int i = 0; i < diskEmitterRays; ++i) {
            const double r0   = rIsco + (rOut - rIsco) * frand();
            const double phi0 = 2.0 * M_PI * frand();
            // emission angle α relative to local radial-outward; pick from a
            // broad bell so most rays escape but some plunge.
            const double alpha = 0.2 + 0.6 * (frand() - 0.5);   // ≈ [-0.1, 0.5]·… → keep π/2-ish
            const double alphaFinal = M_PI * 0.5 + 0.9 * (frand() - 0.5);
            (void)alpha;   // (kept the broader-band line above for readability)

            // λ(r) = λ_ISCO · (r / r_ISCO)^(3/4)  (thermal multi-colour disk)
            const double lambda_rest = lambdaIsco_nm * std::pow(r0 / rIsco, 0.75);

            Photon p;
            p.wavelength_nm     = lambda_rest;
            p.wavelength_obs_nm = lambda_rest;
            p.computeEmissionPath(bh.metric, r0, phi0, alphaFinal);
            p.fromEmitter = true;
            if (!p.captured) addSpectrumSample(p.wavelength_obs_nm);
            emittedPhotons.push_back(std::move(p));
        }

        constexpr size_t EMITTED_CAP = 240;
        if (emittedPhotons.size() > EMITTED_CAP) {
            emittedPhotons.erase(emittedPhotons.begin(),
                emittedPhotons.begin() + (emittedPhotons.size() - EMITTED_CAP));
        }
    }

    // ── Kerr equatorial overlay rebuild ─────────────────────────────────────
    void updateKerrOverlay() {
        Kerr k;
        k.M = bh.metric.M;
        k.a = std::clamp(kerrSpin, 0.0, 0.998) * bh.metric.M;
        const double bMax = 8.0 * bh.metric.M;
        kerrRays = sweepKerrEquatorial(k, 24, bMax, 60.0 * bh.metric.M);
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
    void startMerger(double massSolar, unsigned int windowHeight) {
        startMerger(massSolar, windowHeight, MergerSecondaryKind::BlackHole);
    }

    void startMerger(double massSolar, unsigned int /*windowHeight*/,
                     MergerSecondaryKind kind) {
        merger.active          = true;
        merger.completed       = false;
        merger.secondaryKind   = kind;
        merger.massSolar       = massSolar;
        merger.r_M             = 40.0;
        merger.phi             = M_PI / 4.0;  // start offset so it’s immediately visible
        merger.flashTimer      = 0.0;
        merger.preMergeMassSolar = 0.0;
        merger.secondaryMassGeom  = units::solarMassToGeomMeters(massSolar);
        merger.trail2.clear();
        merger.trail1.clear();
        // Mass ratio q = M2/(M1+M2). For q > 0.1 both BHs visibly orbit the
        // barycenter, the death spiral.
        const double m1Sol = currentMassSolar();
        const double totM  = m1Sol + massSolar;
        merger.massRatioQ = (totM > 0.0) ? massSolar / totM : 0.0;
    }

    // Apply gravitational perturbation from the in-spiraling secondary BH onto all
    // orbiting bodies. Uses Newtonian tidal acceleration (M2/d²) in geometric units,
    // applied as a proper-time-scaled impulse to vr and L, then E is recomputed.
    // The perturbation is physically negligible for large mass-ratio mergers (e.g.,
    // Sgr A* into TON 618) and dramatic for equal-mass mergers (TON 618 + TON 618).
    //
    // TODO(physics-honesty): the primary BH itself never moves in the body
    // integrator — geodesics are integrated around a metric pinned to the
    // origin while the renderer wobbles the primary for show. So the bodies
    // feel a kick from where the secondary actually is in barycentric
    // coordinates, but the primary's gravitational well stays politely
    // stationary. To fix properly: integrate in the barycentric frame and
    // shift the metric source to the primary's true position each step. Until
    // then, this is an honest approximation for q<<1 and a cinematic fib for
    // q~0.5. Don't tell the reviewers.
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

            // Bodies perturbed into escape trajectories during inspiral fly off screen -
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
                // The instantaneous mass increase unbound this body, it’s being
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
        // Only meaningful for BH–BH coalescence; NS/star captures emit far less GW
        // energy at merger and don’t produce a measurable recoil.
        if (merger.secondaryKind == MergerSecondaryKind::BlackHole)
        {
            const double m1 = merger.preMergeMassSolar;
            const double m2 = merger.massSolar;
            const double q  = (m1 > 0.0 && m2 > 0.0)
                            ? std::min(m1, m2) / std::max(m1, m2)
                            : 0.0;   // mass ratio [0,1]; 1 = equal mass
            // visual kick speed as a fraction of c (G=c=1 geometric units, so
            // velocities are dimensionless). 0.01 ≈ 3000 km/s at q=1 — the
            // upper end of measured GW super-kicks. Previously this was
            // multiplied by newM, which is a length in geometric metres, so
            // for SMBH mergers the kick became enormously superluminal and
            // instantly ejected every surviving body.
            const double v_kick = 0.01 * q;
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
        // BH–BH only — NS and star captures get their own flavour below.
        if (merger.secondaryKind == MergerSecondaryKind::BlackHole)
        {
            tidalEvent.active     = true;
            tidalEvent.kind       = TidalEventKind::MergerShockwave;
            tidalEvent.flashTimer = TidalEvent::FLASH_DURATION * 2.0;  // longer flash for a merger
            tidalEvent.eventX     = 0.0;
            tidalEvent.eventY     = 0.0;

            uint32_t rng = 0xDEADBEEFu ^ static_cast<uint32_t>(newM * 1e-6);
            auto rand01 = [&]() -> double {
                rng = rng * 1664525u + 1013904223u;
                return (rng >> 8) / double(1 << 24);
            };

            // Fast outward ring, the “shockwave” expanding at a large fraction of c.
            // Geometric units: c=1, so velocities are dimensionless. 0.5 = 0.5c.
            const double V_RING = 0.5;
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
            // Slower inward-then-outward fallback debris, accretion disk splashback
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

        // NS / pulsar capture — short, hot X-ray burst with no slow fallback.
        // The neutron star is shredded into nuclear-density debris that flashes
        // briefly in X-rays as it crosses the horizon; no GW ring at our visual
        // scale, no measurable BH recoil.
        if (merger.secondaryKind == MergerSecondaryKind::NeutronStar ||
            merger.secondaryKind == MergerSecondaryKind::Pulsar)
        {
            tidalEvent.active     = true;
            tidalEvent.kind       = TidalEventKind::XrayBurst;
            tidalEvent.flashTimer = TidalEvent::FLASH_DURATION * 1.5;
            tidalEvent.eventX     = 0.0;
            tidalEvent.eventY     = 0.0;
            tidalEvent.particles.clear();

            uint32_t rng = 0xC0FFEE11u ^ static_cast<uint32_t>(merger.massSolar * 7919.0);
            auto rand01 = [&]() -> double {
                rng = rng * 1664525u + 1013904223u;
                return (rng >> 8) / double(1 << 24);
            };

            // Fast, sparse burst — fewer particles, higher velocity, short life.
            const double V_BURST = 0.7;
            for (int i = 0; i < 36; ++i) {
                const double angle = 2.0 * M_PI * i / 36.0 + rand01() * 0.15;
                const double speed = V_BURST * (0.85 + 0.3 * rand01());
                TidalDebrisParticle p;
                p.x = 0.0; p.y = 0.0;
                p.vx = std::cos(angle) * speed;
                p.vy = std::sin(angle) * speed;
                p.maxLifetime = 1.2 + 0.8 * rand01();
                p.lifetime    = p.maxLifetime;
                p.size        = 2.0f + (float)(1.2 * rand01());
                p.isFallback  = false;
                tidalEvent.particles.push_back(p);
            }
        }

        // Star / white-dwarf capture — the BH just ate a fluffy object. The
        // realistic picture is a tidal disruption event: a long, curving
        // stream of stellar debris re-accreting over many dynamical times.
        // Reuse the existing TidalDisruption particle path at the origin.
        if (merger.secondaryKind == MergerSecondaryKind::Star ||
            merger.secondaryKind == MergerSecondaryKind::WhiteDwarf)
        {
            tidalEvent.active     = true;
            tidalEvent.kind       = TidalEventKind::TidalDisruption;
            tidalEvent.flashTimer = TidalEvent::FLASH_DURATION;
            tidalEvent.eventX     = 0.0;
            tidalEvent.eventY     = 0.0;
            tidalEvent.particles.clear();

            uint32_t rng = 0xFEEDFACEu ^ static_cast<uint32_t>(merger.massSolar * 1234.5);
            auto rand01 = [&]() -> double {
                rng = rng * 1664525u + 1013904223u;
                return (rng >> 8) / double(1 << 24);
            };

            // Bound + unbound branches of the fallback stream. WDs are denser
            // than MS stars so their disruption stream is shorter-lived; we
            // shave the lifetime a touch for the WD case.
            const bool isWD = (merger.secondaryKind == MergerSecondaryKind::WhiteDwarf);
            const double lifeScale = isWD ? 0.7 : 1.0;
            const double V_TIDAL = 0.18;  // dimensionless: ~0.18c, near-periapsis speed
            for (int i = 0; i < 48; ++i) {
                const double angle = 2.0 * M_PI * rand01();
                const double speed = V_TIDAL * (0.4 + 1.2 * rand01());
                const bool unbound = (rand01() > 0.5);
                TidalDebrisParticle p;
                p.x = newM * (rand01() - 0.5) * 2.0;
                p.y = newM * (rand01() - 0.5) * 2.0;
                p.vx = std::cos(angle) * speed * (unbound ? 1.0 : 0.6);
                p.vy = std::sin(angle) * speed * (unbound ? 1.0 : 0.6);
                p.maxLifetime = lifeScale * (unbound ? 4.0 + 2.0 * rand01()
                                                     : 8.0 + 4.0 * rand01());
                p.lifetime    = p.maxLifetime;
                p.size        = (isWD ? 1.6f : 2.4f) + (float)(1.0 * rand01());
                p.isFallback  = !unbound;
                p.unbound     = unbound;
                tidalEvent.particles.push_back(p);
            }
        }
    }

    void updateMerger(double dt) {
        if (!merger.active) return;

        // Apply the user-selected pacing dial. Flash phase keeps wall-clock
        // timing so the burn-in doesn't stretch out at Cinematic.
        const double inspiralDt = dt * MergerState::timeScaleFactor(merger.timeScale);

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

        applyMergerPerturbation(inspiralDt);

        // GW inspiral — Peters (1964) shape with the correct mass-ratio scaling.
        //   dr/dt ∝ -M1 * M2 * Mtot / r^3
        // In our normalised coordinate x = r / M1 this becomes
        //   dx/dt ∝ -eta / x^3   with eta = M1*M2/(M1+M2)^2 = q*(1-q)
        // (the symmetric mass ratio, peaks at 0.25 for equal masses, → 0 for
        // extreme q). Normalising by 0.25 means an equal-mass merger keeps
        // the old pacing and extreme-q mergers spiral noticeably slower —
        // a TON 618 swallowing a stellar BH now takes its sweet time, as it
        // should.
        //
        // TODO(realism): true Peters wall-clock also scales with Mtot itself
        // (heavier binaries shrink faster in metres/second). Honouring that
        // would make a TON 618+TON 618 merger play out ~10^9× faster than a
        // stellar one, which is correct but unwatchable. Pacing is currently
        // mass-normalised on purpose. If you ever add a "realtime / slow-mo /
        // cinematic" speed dial, plug the Mtot factor in there and apologise
        // to anyone who picks "realtime" for an SMBH merger.
        const double q_sym = merger.massRatioQ;            // M2/(M1+M2) ∈ [0,1]
        const double eta   = q_sym * (1.0 - q_sym);        // symmetric mass ratio
        const double etaScale = (eta > 0.0) ? (eta / 0.25) : 0.0;

        double r3 = merger.r_M * merger.r_M * merger.r_M;
        if (r3 < 1e-30) r3 = 1e-30;
        merger.r_M -= (MergerState::INSPIRAL_K * etaScale * inspiralDt) / r3;

        // Keplerian angular advance — Ω ∝ r^{-3/2} in units of 1/M, which is
        // exactly the shape below. Mass cancels out in these normalised
        // coordinates so OMEGA_K stays constant; for inertial wall-clock
        // accuracy see the TODO on dr/dt above.
        double r_safe = std::max(merger.r_M, 0.1);
        merger.phi += (MergerState::OMEGA_K * inspiralDt) / std::pow(r_safe, 1.5);

        // Record trail positions in BARYCENTER-centred geometric units, so the
        // renderer can draw both BHs spiralling toward the common centre of
        // mass. The renderer uses the screen centre (= world origin = barycentre)
        // as the trail draw origin.
        //   secondary orbits at  r2 = +(1-q) * sep   along the unit vector
        //   primary   orbits at  r1 = -q     * sep   along the unit vector
        {
            const double q    = merger.massRatioQ;   // M2/(M1+M2) in [0,0.5]
            const double sep  = merger.r_M * bh.metric.M;  // geometric separation
            const double cosP = std::cos(merger.phi);
            const double sinP = std::sin(merger.phi);

            const double r2 = sep * (1.0 - q);
            merger.trail2.emplace_back(r2 * cosP, r2 * sinP);
            if ((int)merger.trail2.size() > MergerState::MAX_TRAIL)
                merger.trail2.pop_front();

            // Primary barycentric wobble (only meaningful for q > 0.01)
            if (q > 0.01) {
                const double r1 = sep * q;
                merger.trail1.emplace_back(-r1 * cosP, -r1 * sinP);
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
        // Remove all bodies (including the default test body), it doesn't belong
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

    /*--------- Custom-body presets (persist to disk) ---------*/
    struct CustomBodyPreset {
        std::string    name;
        GalaxyBodyType type;
        double         semiMajorM;
        double         ecc;
    };
    std::vector<CustomBodyPreset> customPresets;

    static const char* bodyTypeKey(GalaxyBodyType t) {
        switch (t) {
            case GalaxyBodyType::Star:           return "Star";
            case GalaxyBodyType::GasCloud:       return "GasCloud";
            case GalaxyBodyType::StellarCluster: return "StellarCluster";
            case GalaxyBodyType::DwarfGalaxy:    return "DwarfGalaxy";
            case GalaxyBodyType::NeutronStar:    return "NeutronStar";
            case GalaxyBodyType::WhiteDwarf:     return "WhiteDwarf";
            case GalaxyBodyType::CompanionStar:  return "CompanionStar";
        }
        return "Star";
    }
    static bool bodyTypeFromKey(const std::string& s, GalaxyBodyType& out) {
        if (s == "Star")           { out = GalaxyBodyType::Star;           return true; }
        if (s == "GasCloud")       { out = GalaxyBodyType::GasCloud;       return true; }
        if (s == "StellarCluster") { out = GalaxyBodyType::StellarCluster; return true; }
        if (s == "DwarfGalaxy")    { out = GalaxyBodyType::DwarfGalaxy;    return true; }
        if (s == "NeutronStar")    { out = GalaxyBodyType::NeutronStar;    return true; }
        if (s == "WhiteDwarf")     { out = GalaxyBodyType::WhiteDwarf;     return true; }
        if (s == "CompanionStar")  { out = GalaxyBodyType::CompanionStar;  return true; }
        return false;
    }

    static std::string customPresetsPath() {
        const char* home = getenv("HOME");
        if (!home) home = "/tmp";
        std::string dir;
#ifdef __APPLE__
        dir = std::string(home) + "/Library/Application Support/Aetherion/saves/";
#elif defined(__linux__)
        dir = std::string(home) + "/.local/share/Aetherion/saves/";
#else
        dir = std::string(home) + "/Aetherion/saves/";
#endif
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
        return dir + "custom_presets.tsv";
    }

    void loadCustomPresets() {
        customPresets.clear();
        std::ifstream f(customPresetsPath());
        if (!f) return;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            // Format: name<TAB>typeKey<TAB>semiMajorM<TAB>ecc
            size_t t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            size_t t2 = line.find('\t', t1 + 1);
            if (t2 == std::string::npos) continue;
            size_t t3 = line.find('\t', t2 + 1);
            if (t3 == std::string::npos) continue;
            CustomBodyPreset p;
            p.name = line.substr(0, t1);
            std::string typeKey = line.substr(t1 + 1, t2 - t1 - 1);
            if (!bodyTypeFromKey(typeKey, p.type)) continue;
            try {
                p.semiMajorM = std::stod(line.substr(t2 + 1, t3 - t2 - 1));
                p.ecc        = std::stod(line.substr(t3 + 1));
            } catch (...) { continue; }
            customPresets.push_back(std::move(p));
        }
    }

    bool saveCustomPresets() const {
        std::ofstream f(customPresetsPath(), std::ios::trunc);
        if (!f) return false;
        f << "# Aetherion 2D custom body presets v1\n";
        f << "# name<TAB>type<TAB>semiMajor_M<TAB>eccentricity\n";
        for (const auto& p : customPresets) {
            // Strip tabs/newlines from name to keep the format unambiguous.
            std::string safe = p.name;
            for (char& c : safe) if (c == '\t' || c == '\n' || c == '\r') c = ' ';
            f << safe << '\t' << bodyTypeKey(p.type) << '\t'
              << std::fixed << std::setprecision(4) << p.semiMajorM << '\t'
              << std::fixed << std::setprecision(4) << p.ecc << '\n';
        }
        return f.good();
    }

    // Add+persist; if name is empty, synthesise one. Returns the stored name.
    std::string appendCustomPreset(std::string name, GalaxyBodyType type,
                                   double semiMajorM, double ecc) {
        if (ecc < 0.0) ecc = 0.0;
        if (ecc > 0.99) ecc = 0.99;
        if (semiMajorM < 2.5) semiMajorM = 2.5;
        if (name.empty()) {
            std::ostringstream os;
            os << bodyTypeName(type) << " a=" << std::fixed << std::setprecision(1)
               << semiMajorM << "M e=" << std::setprecision(2) << ecc;
            name = os.str();
        }
        customPresets.push_back({name, type, semiMajorM, ecc});
        saveCustomPresets();
        return name;
    }

    void removeCustomPreset(int idx) {
        if (idx < 0 || idx >= (int)customPresets.size()) return;
        customPresets.erase(customPresets.begin() + idx);
        saveCustomPresets();
    }

    // User-defined body added via the custom-body creator menu. Appended to the
    // existing `bodies` vector (does not clear other bodies). Tagged as a galaxy
    // body so the existing render / event paths pick it up unchanged.
    // semiMajorM is in units of primary M; ecc is clamped to [0, 0.99].
    void addCustomBody(GalaxyBodyType type, double semiMajorM, double ecc,
                       const char* label = nullptr) {
        if (ecc < 0.0) ecc = 0.0;
        if (ecc > 0.99) ecc = 0.99;
        if (semiMajorM < 2.5) semiMajorM = 2.5;
        const double M = bh.metric.M;
        OrbitingBody ob(bh.metric, semiMajorM * M, ecc);
        ob.label = label ? label : bodyTypeName(type);
        ob.bodyType = type;
        ob.isGalaxyBody = true;
        ob.phi = static_cast<double>(bodies.size()) * 0.7;
        ob.trail.clear();
        ob.trail.emplace_back(ob.worldX(), ob.worldY());
        bodies.push_back(std::move(ob));
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
        clearTransientEvents();
        // note: this does NOT reset the exportName. you're welcome to argue that it should,
        // but I've already accidentally lost data once by having reset() wipe a name I forgot to save.
        // so it stays. fight me.
    }

    // Clear in-flight tidal/merger event state so stale "MERGER REMNANT" /
    // "X-RAY BURST" labels and inspiral trails don't bleed into a freshly
    // switched preset. Safe to call any time.
    void clearTransientEvents() {
        tidalEvent.active     = false;
        tidalEvent.flashTimer = 0.0;
        tidalEvent.particles.clear();
        merger.active     = false;
        merger.completed  = false;
        merger.flashTimer = 0.0;
        merger.trail1.clear();
        merger.trail2.clear();
    }

    /*--------- Lensing analysis ---------*/
    void buildLensingData() {
        lensingData.criticalImpactParam = bh.metric.criticalImpact();
        lensingData.deflectionTable.clear();
        lensingData.causticPoints.clear();
        lensingData.einsteinRings.clear();

        for (const auto& p : photons) {
            // Skip emitter-spawned photons in the global deflection table:
            // they don't represent the "ray from infinity" geometry that
            // the lensing analytics assume.
            if (p.fromEmitter) continue;
            PhotonDeflection d;
            d.impactParameter = p.impactParameter;
            d.captured = p.captured;
            d.deflectionAngle = p.deflectionAngle;
            d.deflectionDeg   = p.deflectionDeg;
            d.shapiroDelay_M  = p.shapiroDelay;
            d.coordTime_M     = p.coordTime;
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

            // threshold of 0.5 rad/M is completely empirical, I tuned it until the caustic markers
            // appeared in roughly the right places visually. I cannot rigorously justify this number.
            // if caustics are appearing in weird places, this is probably why.
            if (gradient > 0.5) {
                float avgB = (float)(0.5 * (prev.impactParameter + curr.impactParameter));
                lensingData.causticPoints.emplace_back(avgB, (float)gradient);
            }
        }

        // Einstein-ring / image-multiplicity finder.
        // A photon that picks up Δ = (2n−1)π of deflection comes back to the
        // source line on the opposite side ⇒ produces an image. n=1 is the
        // primary Einstein ring, n=2,3,… are higher-order relativistic rings.
        // We sweep the deflection table separately on the b>0 and b<0 halves
        // (which correspond to images above / below the lens), monotone in
        // |Δ| between b_crit and b → ∞, so a sign-change in (Δ − (2n−1)π)
        // between adjacent samples brackets exactly one image of order n.
        auto sweepHalf = [&](bool positive) {
            int sgn = positive ? +1 : -1;
            // Collect non-captured entries with this sign of b, sorted by |b| ascending
            std::vector<const PhotonDeflection*> halfTable;
            halfTable.reserve(lensingData.deflectionTable.size());
            for (const auto& d : lensingData.deflectionTable) {
                if (d.captured) continue;
                if ((d.impactParameter > 0.0) != positive) continue;
                halfTable.push_back(&d);
            }
            std::sort(halfTable.begin(), halfTable.end(),
                [](const PhotonDeflection* a, const PhotonDeflection* b){
                    return std::abs(a->impactParameter) < std::abs(b->impactParameter);
                });

            // Look for crossings up through n = 3 (primary + 2 relativistic rings).
            // Higher orders need rays within ~exp(−2π) of b_crit, which the
            // uniform-in-b sampler does not provide — high-res lensing mode adds
            // log-spaced rays right where they're needed.
            for (int n = 1; n <= 3; ++n) {
                double target = (2.0 * n - 1.0) * M_PI; // |Δ| target
                for (size_t i = 1; i < halfTable.size(); ++i) {
                    double da = std::abs(halfTable[i-1]->deflectionAngle);
                    double db = std::abs(halfTable[i  ]->deflectionAngle);
                    double lo = std::min(da, db), hi = std::max(da, db);
                    if (target < lo || target > hi) continue;
                    // Linear interp on |b| vs |Δ|
                    double b1 = std::abs(halfTable[i-1]->impactParameter);
                    double b2 = std::abs(halfTable[i  ]->impactParameter);
                    double denom = (db - da);
                    if (std::abs(denom) < 1e-18) continue;
                    double t = (target - da) / denom;
                    double bRing = b1 + t * (b2 - b1);
                    EinsteinRingImage img;
                    img.order = n;
                    img.impactParameter = sgn * bRing;
                    img.deflectionDeg   = target * 180.0 / M_PI;
                    img.sign = sgn;
                    lensingData.einsteinRings.push_back(img);
                    break; // one image of this order per side
                }
            }
        };
        sweepHalf(true);
        sweepHalf(false);
    }

    /*--------- ISCO validation test ---------*/
    // Returns true if any scenario-tagged bodies were removed (toggle-off).
    bool removeBodiesByScenarioTag(ResearchScenario tag) {
        size_t before = bodies.size();
        bodies.erase(
            std::remove_if(bodies.begin(), bodies.end(),
                [tag](const OrbitingBody& b){ return b.scenarioTag == tag; }),
            bodies.end());
        if (selectedBodyIdx >= (int)bodies.size())
            selectedBodyIdx = bodies.empty() ? 0 : (int)bodies.size() - 1;
        return bodies.size() != before;
    }

    bool anyBodyWithScenarioTag(ResearchScenario tag) const {
        for (const auto& b : bodies)
            if (b.scenarioTag == tag) return true;
        return false;
    }

    void startISCOTest() {
        // Toggle-off path: if ISCO test bodies are already present, remove
        // them and leave the rest of the scene untouched.
        if (anyBodyWithScenarioTag(ResearchScenario::ISCOTest)) {
            removeBodiesByScenarioTag(ResearchScenario::ISCOTest);
            iscoResults.clear();
            if (activeScenario == ResearchScenario::ISCOTest)
                activeScenario = ResearchScenario::None;
            return;
        }

        activeScenario = ResearchScenario::ISCOTest;
        iscoResults.clear();

        double M = bh.metric.M;
        double testRadii[] = { 5.0 * M, 6.0 * M, 7.0 * M };
        const char* labels[] = { "5M (Unstable)", "6M (Critical)", "7M (Stable)" };

        for (int i = 0; i < 3; ++i) {
            OrbitingBody ob(bh.metric, testRadii[i], 0.0);
            ob.initCircularOrbit(bh.metric, testRadii[i]);
            ob.label = labels[i];
            ob.scenarioTag = ResearchScenario::ISCOTest;
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
        selectedBodyIdx = (int)bodies.size() - 1;
    }

    void updateISCOResults() {
        // Locate ISCO test bodies (the order in which they were inserted is
        // preserved, so we can still pair them with iscoResults entries).
        std::vector<size_t> iscoIdx;
        iscoIdx.reserve(iscoResults.size());
        for (size_t i = 0; i < bodies.size(); ++i)
            if (bodies[i].scenarioTag == ResearchScenario::ISCOTest)
                iscoIdx.push_back(i);

        for (size_t i = 0; i < iscoResults.size() && i < iscoIdx.size(); ++i) {
            const auto& b = bodies[iscoIdx[i]];
            double M = bh.metric.M;
            iscoResults[i].captured = b.captured;
            iscoResults[i].radiusDrift = (b.r - iscoResults[i].testRadius_M * M) / M;
            iscoResults[i].survivalTime = b.measurement.properTime;
        }
    }

    /*--------- Photon sphere test ---------*/
    // Sweep preset cycles through three epsilon windows. Pressing the
    // "test photon sphere" key while the test is already active steps to
    // the next preset and re-runs the sweep, so the user can drill in from
    // the broad-band view down to ε ~ 1e−4 of the photon sphere.
    int photonSphereSweepPreset = 0;
    static constexpr int N_PHOTON_SPHERE_PRESETS = 3;
    const char* photonSphereSweepName() const {
        switch (photonSphereSweepPreset) {
            case 0: return "broad  (eps=1e-3..0.05)";
            case 1: return "narrow (eps=1e-4..1e-2)";
            case 2: return "ultra  (eps=1e-5..1e-3)";
            default: return "broad";
        }
    }

    void startPhotonSphereTest() {
        // Toggle-off when already active and we're back at preset 0 after a
        // full cycle: clear results and exit the scenario. Otherwise cycle
        // the sweep preset and re-run.
        if (activeScenario == ResearchScenario::PhotonSphereTest) {
            photonSphereSweepPreset = (photonSphereSweepPreset + 1) % N_PHOTON_SPHERE_PRESETS;
        } else {
            photonSphereSweepPreset = 0;
        }

        activeScenario = ResearchScenario::PhotonSphereTest;
        photonSphereResults.clear();

        const double b_crit = bh.metric.criticalImpact();

        // Build the epsilon list for this preset. Each preset symmetrically
        // samples both sides of b_crit (b = b_crit · (1 ± ε)) so the data
        // panel shows both grazing-escape and just-barely-captured rays.
        std::vector<double> magnitudes;
        switch (photonSphereSweepPreset) {
            case 0: magnitudes = { 0.001, 0.01,  0.025, 0.05 };          break;
            case 1: magnitudes = { 1e-4, 5e-4, 1e-3, 5e-3, 1e-2 };       break;
            case 2: magnitudes = { 1e-5, 5e-5, 1e-4, 5e-4, 1e-3 };       break;
        }
        std::vector<double> epsilons;
        epsilons.reserve(magnitudes.size() * 2);
        for (double m : magnitudes) { epsilons.push_back(+m); epsilons.push_back(-m); }

        for (double eps : epsilons) {
            double b = b_crit * (1.0 + eps);

            Photon p;
            p.impactParameter = b;
            p.computePath(bh.metric, params.rMaxIntegrate, 0.001);

            PhotonSphereTestResult res;
            res.impactParameter = b;
            res.epsilon         = eps;
            res.captured        = p.captured;
            res.escaped         = !p.captured;
            res.rMin_M          = p.rMin_M;
            res.deflectionDeg   = p.deflectionDeg;
            res.shapiroDelay_M  = p.shapiroDelay;
            res.coordTime_M     = p.coordTime;

            if (!p.captured) {
                res.stabilityAngle    = p.deflectionAngle;
                res.orbitsBeforeDecay = std::abs(p.deflectionAngle) / (2.0 * M_PI);
            } else {
                // For captured photons, estimate from path length.
                if (p.path.size() >= 2) {
                    double totalAngle = 0.0;
                    for (size_t i = 1; i < p.path.size(); ++i) {
                        double a1 = std::atan2(p.path[i-1].y, p.path[i-1].x);
                        double a2 = std::atan2(p.path[i  ].y, p.path[i  ].x);
                        double da = a2 - a1;
                        if (da >  M_PI) da -= 2.0 * M_PI;
                        if (da < -M_PI) da += 2.0 * M_PI;
                        totalAngle += std::abs(da);
                    }
                    res.orbitsBeforeDecay = totalAngle / (2.0 * M_PI);
                    res.stabilityAngle    = totalAngle;
                }
            }
            photonSphereResults.push_back(res);
        }

        // Sort by |epsilon| ascending so the panel reads from "closest to
        // b_crit" outward, regardless of insertion order.
        std::sort(photonSphereResults.begin(), photonSphereResults.end(),
            [](const PhotonSphereTestResult& a, const PhotonSphereTestResult& b){
                return std::abs(a.epsilon) < std::abs(b.epsilon);
            });
    }

    /*--------- Radial infall test ---------*/
    void startRadialInfall() {
        // Toggle-off
        if (anyBodyWithScenarioTag(ResearchScenario::RadialInfall)) {
            removeBodiesByScenarioTag(ResearchScenario::RadialInfall);
            if (activeScenario == ResearchScenario::RadialInfall)
                activeScenario = ResearchScenario::None;
            return;
        }

        activeScenario = ResearchScenario::RadialInfall;

        double M = bh.metric.M;
        double r0 = 20.0 * M; // start at 20M because anything closer and the thing hits the horizon
                               // before the user has time to read the data panel. tried 10M once. not great.

        OrbitingBody ob(bh.metric, r0, 0.0);
        ob.initRadialInfall(bh.metric, r0);
        ob.label = "Radial infall from 20M";
        ob.scenarioTag = ResearchScenario::RadialInfall;
        bodies.push_back(std::move(ob));
        selectedBodyIdx = (int)bodies.size() - 1;
    }

    /*--------- Pulsar orbital simulation ---------*/
    // Find the pulsar body, if one currently exists. Returns -1 when absent.
    int findPulsarIdx() const {
        for (size_t i = 0; i < bodies.size(); ++i)
            if (bodies[i].isPulsar) return static_cast<int>(i);
        return -1;
    }

    // Toggle a pulsar companion in/out of the active scene WITHOUT disturbing
    // any other orbiting bodies (galaxy stars, default test body, etc).
    //   - If a pulsar is already present: remove it, clear pulsar state.
    //   - Otherwise: append a fresh pulsar at a0 = 20M, e0 = 0.3.
    // Always re-initialises pulsar state to the canonical starting orbit on
    // re-enable, so repeated presses never accumulate inspiral drift.
    void togglePulsarOrbital() {
        // Already active → toggle OFF
        if (findPulsarIdx() >= 0) {
            removeBodiesByScenarioTag(ResearchScenario::PulsarOrbital);
            pulsarState = PulsarState{};      // fully reset (active=false)
            pulsarData  = PulsarOrbitalData{};
            if (activeScenario == ResearchScenario::PulsarOrbital)
                activeScenario = ResearchScenario::None;
            return;
        }

        // Toggle ON, append pulsar without clearing existing bodies.
        activeScenario = ResearchScenario::PulsarOrbital;

        double M   = bh.metric.M;
        double a0  = 20.0 * M;  // start at 20M, comfortably outside ISCO at 6M
        double e0  = 0.3;       // moderate eccentricity for visible variation

        OrbitingBody ob(bh.metric, a0, e0);
        ob.label       = "Pulsar (1.4 Msun NS)";
        ob.isPulsar    = true;
        ob.scenarioTag = ResearchScenario::PulsarOrbital;
        bodies.push_back(std::move(ob));
        selectedBodyIdx = (int)bodies.size() - 1;

        // Initialise pulsar physics state (fully reset every enable)
        pulsarState = PulsarState{};
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

    // Backwards-compatible alias (kept so callers that say "start" still work).
    void startPulsarOrbital() { togglePulsarOrbital(); }

    // Adiabatic GW-decay update called every sim frame when scenario is PulsarOrbital.
    void updatePulsar(double dtSim) {
        if (!pulsarState.active) return;
        int pidx = findPulsarIdx();
        if (pidx < 0) return;

        auto& body = bodies[pidx];
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
        // Toggle-off
        if (anyBodyWithScenarioTag(ResearchScenario::TidalDisruption)) {
            removeBodiesByScenarioTag(ResearchScenario::TidalDisruption);
            tidalDisrupted = false;
            if (activeScenario == ResearchScenario::TidalDisruption)
                activeScenario = ResearchScenario::None;
            return;
        }

        activeScenario  = ResearchScenario::TidalDisruption;
        tidalDisrupted  = false;

        double M = bh.metric.M;

        // eccentricity 0.85, this is deliberately extreme so the body actually reaches tidal disruption range.
        // if you set it lower the body just orbits happily and nothing dramatic happens, which defeats the point.
        // if it crosses ISCO it'll plunge in, which is also fine, that's kind of what tidal disruption events do.
        OrbitingBody ob(bh.metric, 8.0 * M, 0.85);
        ob.label       = "Tidal test body";
        ob.scenarioTag = ResearchScenario::TidalDisruption;
        bodies.push_back(std::move(ob));
        selectedBodyIdx = (int)bodies.size() - 1;
    }

    void updateTidalDisruption(double /*dtSim*/) {
        if (tidalDisrupted) return;
        // Find the tagged tidal test body.
        int tidx = -1;
        for (size_t i = 0; i < bodies.size(); ++i) {
            if (bodies[i].scenarioTag == ResearchScenario::TidalDisruption) {
                tidx = (int)i;
                break;
            }
        }
        if (tidx < 0) return;
        const auto& body = bodies[tidx];
        if (body.captured) {
            tidalDisrupted = true;
            return;
        }
        double r_tidal = tidalRadiusM * bh.metric.M;
        if (body.r < r_tidal) {
            tidalDisrupted = true;
            triggerTidalDisruption(tidx);
        }
    }

    // Spawn a visible tidal disruption event for the body at bodyIdx, then
    // remove it from the bodies list.
    // Physics: the disrupted body is treated as having a small radius R⋆
    // (≪ r_tidal) so each fragment inherits the body's bulk velocity plus a
    // frozen-in energy spread dE = ± M·R⋆/r² (Lacy/Carter '82, Rees '88).
    // Fragments with E < 0 are bound and return to periapsis on Keplerian
    // ellipses → the classic t⁻⁵ᐟ³ fallback stream. Unbound fragments
    // (E > 0) escape on hyperbolae. We then integrate every particle under
    // Newtonian gravity in updateTidalEvent() so the stream actually
    // *looks* like a stream instead of straight-line debris.
    void triggerTidalDisruption(int bodyIdx) {
        if (bodyIdx < 0 || bodyIdx >= (int)bodies.size()) return;
        const auto& body = bodies[bodyIdx];

        const double M  = bh.metric.M;
        const double cx = body.r * std::cos(body.phi);
        const double cy = body.r * std::sin(body.phi);

        tidalEvent.active     = true;
        tidalEvent.kind       = TidalEventKind::TidalDisruption;
        tidalEvent.flashTimer = TidalEvent::FLASH_DURATION;
        tidalEvent.eventX     = cx;
        tidalEvent.eventY     = cy;
        tidalEvent.particles.clear();

        // Bulk velocity of the disrupted body in world (Cartesian) coordinates.
        // Same conversion markEjected() uses.
        const double cos_p = std::cos(body.phi);
        const double sin_p = std::sin(body.phi);
        const double v_tan = body.L / std::max(1e-30, body.r);
        const double Vx    = body.vr * cos_p - v_tan * sin_p;
        const double Vy    = body.vr * sin_p + v_tan * cos_p;

        // Spread direction along the tidal axis (= radial), perpendicular = tangential.
        const double rad_x = cos_p, rad_y = sin_p;
        const double tan_x = -sin_p, tan_y = cos_p;

        // Characteristic spread: at r_tidal the body of radius R⋆ has fragments
        // gaining/losing velocity ~ sqrt(2 M R⋆ / r²). We don't have R⋆ from
        // the body data so we pick R⋆ such that the spread is a small but
        // visible fraction of the orbital velocity at this radius.
        const double v_circ = std::sqrt(M / std::max(1e-30, body.r));
        const double dV     = 0.45 * v_circ; // half-width of the energy spread

        // Deterministic RNG seeded from the disruption position
        uint32_t rng = 0x12345678u ^ static_cast<uint32_t>(std::abs(cx + cy) * 1000.0 + body.r * 17.0);
        auto rand01 = [&]() -> double {
            rng = rng * 1664525u + 1013904223u;
            return (rng >> 8) / double(1 << 24);
        };
        auto randSym = [&]() { return 2.0 * rand01() - 1.0; };

        // 48 main stream particles: linear in the −1..+1 spread coordinate
        // s ∈ [−1, +1]; s < 0 ⇒ slowed ⇒ bound ⇒ returns; s > 0 ⇒ unbound.
        const int N_MAIN = 48;
        for (int i = 0; i < N_MAIN; ++i) {
            double s = (i + rand01() * 0.6 - 0.3) / double(N_MAIN - 1) * 2.0 - 1.0; // jittered linear sweep
            double dv = s * dV;
            // Velocity offset is radial (along the tidal axis), with a tiny tangential jitter for visual width.
            double jitter = randSym() * 0.04 * v_circ;
            double vx = Vx + rad_x * dv + tan_x * jitter;
            double vy = Vy + rad_y * dv + tan_y * jitter;
            // Position offset: small smear along the same axis so we don't all start at one pixel
            double pos_smear = randSym() * 0.6 * M;
            double px = cx + rad_x * pos_smear + tan_x * randSym() * 0.4 * M;
            double py = cy + rad_y * pos_smear + tan_y * randSym() * 0.4 * M;
            double rSq = px*px + py*py;
            double Esp = 0.5 * (vx*vx + vy*vy) - M / std::sqrt(std::max(1e-30, rSq));
            TidalDebrisParticle p;
            p.x = px; p.y = py; p.vx = vx; p.vy = vy;
            p.maxLifetime  = TidalEvent::STREAM_LIFETIME * (0.8 + 0.4 * rand01());
            p.lifetime     = p.maxLifetime;
            p.size         = 1.8f + (float)(1.2 * rand01());
            p.isFallback   = (Esp < 0.0);
            p.specificEnergy = Esp;
            p.unbound      = (Esp >= 0.0);
            tidalEvent.particles.push_back(p);
        }

        // 16 extra "flash" sprites — small isotropic outburst at the disruption point,
        // so the moment of disruption reads visually even before the stream develops.
        for (int i = 0; i < 16; ++i) {
            double angle = 2.0 * M_PI * rand01();
            double speed = 0.35 * v_circ * (0.5 + rand01());
            TidalDebrisParticle p;
            p.x = cx; p.y = cy;
            p.vx = Vx + std::cos(angle) * speed;
            p.vy = Vy + std::sin(angle) * speed;
            p.maxLifetime = TidalEvent::DEBRIS_LIFETIME * (0.6 + 0.4 * rand01());
            p.lifetime    = p.maxLifetime;
            p.size        = 2.5f + (float)(1.5 * rand01());
            p.isFallback  = false;
            double rSq = p.x*p.x + p.y*p.y;
            p.specificEnergy = 0.5 * (p.vx*p.vx + p.vy*p.vy) - M / std::sqrt(std::max(1e-30, rSq));
            p.unbound     = (p.specificEnergy >= 0.0);
            tidalEvent.particles.push_back(p);
        }

        bodies.erase(bodies.begin() + bodyIdx);
    }

    // Advance the tidal particle system using wall-clock dt.
    // Bound particles feel Newtonian gravity from the BH (M/r²) so they
    // produce the familiar fallback stream curving back toward periapsis.
    // Unbound particles already escape; we still apply gravity for visual
    // coherence (orbits curve as they fly away). Particles captured by the
    // horizon are removed.
    void updateTidalEvent(double dt) {
        if (!tidalEvent.active) return;

        if (tidalEvent.flashTimer > 0.0)
            tidalEvent.flashTimer = std::max(0.0, tidalEvent.flashTimer - dt);

        const double M  = bh.metric.M;
        const double rH = bh.metric.horizon();
        const double rH2 = rH * rH;

        // Symplectic Euler (kick-drift): cheap, stable for visualisation,
        // doesn't conserve energy long-term but neither did the old straight-
        // line code, and over ~5–8s lifetimes the drift is invisible.
        for (auto& p : tidalEvent.particles) {
            double r2 = p.x * p.x + p.y * p.y;
            if (r2 < rH2) {
                // Captured — mark for removal via expired lifetime.
                p.lifetime = 0.0;
                continue;
            }
            double r  = std::sqrt(r2);
            double aMag = -M / r2;
            double ax = aMag * (p.x / r);
            double ay = aMag * (p.y / r);
            p.vx += ax * dt;
            p.vy += ay * dt;
            p.x  += p.vx * dt;
            p.y  += p.vy * dt;
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
            ss << "Sweep: " << photonSphereSweepName() << "\n";
            ss << "b_crit = " << sci(bh.metric.criticalImpact()) << " M\n";
            ss << "r_ph   = " << sci(bh.metric.photonSphere())  << " M\n";

            // Aggregate stats over this sweep.
            int escCount = 0, capCount = 0;
            double rMinMin   =  std::numeric_limits<double>::infinity();
            double defMax    = 0.0;
            double orbitsMax = 0.0;
            double shapMax   = 0.0;
            for (const auto& r : photonSphereResults) {
                if (r.escaped) ++escCount; else ++capCount;
                rMinMin   = std::min(rMinMin,   r.rMin_M);
                defMax    = std::max(defMax,    std::abs(r.deflectionDeg));
                orbitsMax = std::max(orbitsMax, r.orbitsBeforeDecay);
                shapMax   = std::max(shapMax,   std::abs(r.shapiroDelay_M));
            }
            ss << "N=" << photonSphereResults.size()
               << "  escaped=" << escCount
               << "  captured=" << capCount << "\n";
            ss << "r_min(min) = " << sci(rMinMin) << " M\n";
            ss << "|Defl|(max)= " << sci(defMax)  << " deg\n";
            ss << "orbits(max)= " << sci(orbitsMax) << "\n";
            ss << "Shap(max)  = " << sci(shapMax) << " M\n";

            ss << "eps      b/M     rmin    deg     orb    shap   st\n";
            for (const auto& r : photonSphereResults) {
                char sign = (r.epsilon >= 0) ? '+' : '-';
                ss << sign << sci(std::abs(r.epsilon))
                   << " " << sci(r.impactParameter)
                   << " " << sci(r.rMin_M)
                   << " " << sci(r.deflectionDeg)
                   << " " << sci(r.orbitsBeforeDecay, 1)
                   << " " << sci(r.shapiroDelay_M)
                   << (r.escaped ? " ESC" : " CAP") << "\n";
            }
            ss << "(re-press test key to cycle sweep preset)\n";
        }

        // Lensing
        if (!lensingData.deflectionTable.empty()) {
            ss << "--- LENSING ---\n";
            ss << "b_crit: " << sci(lensingData.criticalImpactParam) << " M\n";
            ss << "Caustic pts: " << lensingData.causticPoints.size() << "\n";

            // Aggregate deflection + Shapiro stats over the current sweep.
            int    nEsc = 0, nCap = 0;
            double defMin = std::numeric_limits<double>::infinity();
            double defMax = 0.0, defSum = 0.0;
            double shapMin = std::numeric_limits<double>::infinity();
            double shapMax = 0.0, shapSum = 0.0;
            int    nDef = 0, nShap = 0;
            for (const auto& e : lensingData.deflectionTable) {
                if (e.captured) { ++nCap; continue; }
                ++nEsc;
                double d = std::abs(e.deflectionDeg);
                defMin = std::min(defMin, d);
                defMax = std::max(defMax, d);
                defSum += d;
                ++nDef;
                double s = std::abs(e.shapiroDelay_M);
                shapMin = std::min(shapMin, s);
                shapMax = std::max(shapMax, s);
                shapSum += s;
                ++nShap;
            }
            ss << "Rays: esc=" << nEsc << " cap=" << nCap << "\n";
            if (nDef > 0) {
                ss << "|Defl| min/avg/max: "
                   << sci(defMin) << " / "
                   << sci(defSum / nDef) << " / "
                   << sci(defMax) << " deg\n";
            }
            if (nShap > 0) {
                ss << "Shapiro min/avg/max: "
                   << sci(shapMin) << " / "
                   << sci(shapSum / nShap) << " / "
                   << sci(shapMax) << " M\n";
            }
            // Einstein rings (batch 1).
            if (!lensingData.einsteinRings.empty()) {
                ss << "Einstein rings (n=" << lensingData.einsteinRings.size() << "):\n";
                size_t cap = std::min<size_t>(lensingData.einsteinRings.size(), 4);
                for (size_t i = 0; i < cap; ++i) {
                    const auto& er = lensingData.einsteinRings[i];
                    ss << "  b=" << sci(er.impactParameter)
                       << "  defl=" << sci(er.deflectionDeg) << " deg\n";
                }
                if (cap < lensingData.einsteinRings.size())
                    ss << "  ... (+" << (lensingData.einsteinRings.size() - cap) << " more)\n";
            }
            // Spectrum (batch 2).
            if (lensingData.spectrumTotal > 0) {
                int peakBin = 0;
                int peakCount = 0;
                for (int i = 0; i < LensingData::SPECTRUM_BINS; ++i) {
                    if (lensingData.spectrumCounts[i] > peakCount) {
                        peakCount = lensingData.spectrumCounts[i];
                        peakBin = i;
                    }
                }
                double binW = (LensingData::SPECTRUM_MAX_NM - LensingData::SPECTRUM_MIN_NM)
                              / double(LensingData::SPECTRUM_BINS);
                double peakLambda = LensingData::SPECTRUM_MIN_NM + (peakBin + 0.5) * binW;
                ss << "Spectrum: N=" << lensingData.spectrumTotal
                   << "  peak~" << sci(peakLambda, 1) << " nm"
                   << " (count=" << peakCount << ")\n";
            }
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
                return "No orbit data yet, let the simulation run first";
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
                return "No orbit data yet, let the simulation run first";
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
        // Accretion history is global (not per-body), export it whenever we
        // have any sampled entries, regardless of which body is selected.
        if (!accretion.history.empty()) {
            std::string base = makeExportDir();
            if (!base.empty()) {
                if (FITSExport::exportAccretionHistory(base + "accretion.fits", accretion, meta))
                    msg += "accretion.fits ";
                if (displayDir.empty()) displayDir = base;
            }
        }
        // Gas radial profile + SED snapshots, recomputed at export time from
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
                return "No orbit data yet, let the simulation run first";
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
