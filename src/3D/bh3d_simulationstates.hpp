#pragma once
// ============================================================
// bh3d_simulationstates.hpp  (header still says simulation_state,
// never got around to renaming it): plain-data physics snapshot.
//
// It's the hand-off point between the physics side (camera, config)
// and the render side (shaders, bloom, HUD), so neither one has to
// include the other's headers.
// ============================================================

#include <glm/glm.hpp>
#include <string>
#include <vector>

struct PhysicsSnapshot {
    /*--------- Camera ---------*/
    glm::vec3 cameraPos;
    glm::vec3 cameraDir;
    glm::vec3 cameraUp;
    float     fov;
    float     roll;          // radians
    bool      freelook;

    /*--------- Time ---------*/
    float totalTime;
    float animTime;
    float animSpeed;
    float dt;

    /*--------- Toggles / feature flags ---------*/
    bool jetsEnabled;
    bool blrEnabled;
    bool dopplerEnabled;
    bool blueshiftEnabled;
    bool cinematicMode;

    /*--------- Black hole parameters (from config) ---------*/
    float     bhRadius;
    glm::vec3 bhPosition;
    float     bhSpin;

    /*--------- Disk ---------*/
    float diskInnerRadius;
    float diskOuterRadius;
    float diskHalfThickness;
    float diskPeakTemp;
    float diskDisplayTempInner;
    float diskDisplayTempOuter;
    float diskSatBoostInner;
    float diskSatBoostOuter;

    /*--------- Jets ---------*/
    float     jetRadius;
    float     jetLength;
    glm::vec3 jetColor;

    /*--------- BLR ---------*/
    float blrInnerRadius;
    float blrOuterRadius;
    float blrThickness;
    float blrStrength;   // [0,1] cloud opacity scaler

    /*--------- Orbital body ---------*/
    std::vector<glm::vec3> orbBodyPositions;
    std::vector<float>     orbBodyRadii;
    std::vector<glm::vec3> orbBodyColors;
    std::vector<int>       orbBodyTypes;     // GalaxyBody3DType enum value per body
    std::vector<std::string> orbBodyLabels;  // Display name per body (label-view mode)
    bool                   orbBodyEnabled;

    /*--------- Large-scale structures ---------*/
    bool hostGalaxyEnabled; // the host galaxy itself
    bool labEnabled;        // Lyman-alpha Blob overlay
    bool cgmEnabled;        // circumgalactic medium (diffuse gas halo around the galaxy)

    /*--------- Render settings (more steps = higher quality, slower) ---------*/
    int maxSteps;

    /*--------- Profile metadata (for HUD) ---------*/
    std::string profileName;
    double      massSolar = 0.0;  // this stores primary BH mass [M☉] for physical scale readout

    /*--------- True-scale mode ---------*/
    // When this is true, orbBodyRadii are physical (Rs) rather than visual.  Toggled by
    // the user from the Overlays panel; the contrast with the visual default is
    // the educational point.
    bool trueScaleMode = false;
    // Per-body physical radius in Rs (populated when trueScaleMode is on).
    std::vector<float> orbBodyPhysRadii;

    /*--------- Barycentric binary mode (Gaia BH1/2/3) ---------*/
    // When true, bhPosition is offset from origin and the barycenter sits at
    // world origin (= screen centre).  The HUD draws a "Center of Mass" marker.
    bool      barycentricMode = false;
    glm::vec3 barycenterPos   = glm::vec3(0.0f); // always world origin when active

    /*--------- Window ---------*/
    int windowW;
    int windowH;

    /*--------- Performance ---------*/
    float fps; // TODO: expose frame time (ms) alongside FPS for more useful profiling

    /*--------- Tidal disruption event (3D overlay) ---------*/
    bool                   tdeActive      = false;
    float                  tdeFlashAlpha  = 0.0f;   // 0..1
    glm::vec3              tdeEventPos    = glm::vec3(0.0f); // world position
    std::vector<glm::vec3> tdeDebrisPos;             // world positions
    std::vector<float>     tdeDebrisLifeF;           // lifetime fraction [0,1]
    std::vector<bool>      tdeDebrisIsFallback;      // gas-stream vs debris

    /*--------- Black hole merger (3D visual) ---------*/
    // Visual-only merger mirroring the 2D simulator: a secondary object
    // spirals into the primary, coalesces in a flash, and leaves a grown
    // remnant + gravitational-wave shockwave ring. The secondary itself is
    // injected into the orbBody* arrays above (so the shader renders it);
    // these fields drive the HUD overlays (trail, flash, ring, labels).
    bool                   mergerActive   = false;  // inspiral or aftermath running
    bool                   mergerInspiral = false;  // secondary still incoming
    float                  mergerFlashAlpha = 0.0f; // 0..1 full-screen coalescence flash
    glm::vec3              mergerEventPos = glm::vec3(0.0f); // world coalescence point
    bool                   mergerRingActive = false; // GW shockwave ring visible
    float                  mergerRingT    = 0.0f;   // 0..1 ring expansion progress
    int                    mergerSecondaryKind = 0; // MergerSecondaryKind3D
    std::string            mergerSecondaryLabel;    // "Black Hole", "Pulsar", ...
    glm::vec3              mergerSecondaryPos = glm::vec3(0.0f); // world position of incoming body
    float                  mergerSepRs    = 0.0f;   // current separation [Rs] (HUD readout)
    float                  mergerRemnantAlpha = 0.0f; // 0..1 "remnant" label fade
    std::vector<glm::vec3> mergerTrail;             // secondary death-spiral trail (world, oldest→newest)

    // Scrolling gravitational-wave strain record h(t) ∝ (1/sep)·cos(2φ_orbit),
    // sampled per-frame during inspiral. Mirrors the 2D pulsar hHistory; the
    // frequency/amplitude sweep-up is drawn as an ImGui::PlotLines sparkline.
    static constexpr int   MERGER_GW_SAMPLES = 256;
    std::vector<float>     mergerGWWaveform;


    // Disk / spin / jet identity of the inspiralling secondary black hole, so it
    // looks the same as it would if it were the standalone primary. The sizes
    // here are ratios of the secondary's shader radius (orbBodyRadius), not
    // absolute units. Only read when secBHActive is set.
    // (secBH = secondary black hole; kept short so the uniform names line up.)
    bool                   secBHActive     = false;
    glm::vec3              secBHDiskNormal = glm::vec3(0.0f, 1.0f, 0.0f);
    float                  secBHSpin       = 0.0f;
    float                  secBHDiskInner  = 2.2f;   // × secondary radius
    float                  secBHDiskOuter  = 12.0f;  // × secondary radius
    float                  secBHDiskStrength = 0.0f; // 0 = no disk
    glm::vec3              secBHColorInner = glm::vec3(1.0f, 0.85f, 0.55f);
    glm::vec3              secBHColorOuter = glm::vec3(0.9f, 0.45f, 0.20f);
    bool                   secBHShowJets   = false;
    glm::vec3              secBHJetColor   = glm::vec3(0.4f, 0.7f, 1.0f);
    float                  secBHJetRadius  = 0.15f;  // × secondary radius
    float                  secBHJetLength  = 10.0f;  // × secondary radius
};
