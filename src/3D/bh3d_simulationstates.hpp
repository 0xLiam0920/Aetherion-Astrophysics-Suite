#pragma once
// ============================================================
// simulation_state.hpp: Pure-data snapshot of physics state
//
// This struct is the ONLY bridge between the physics layer
// (camera_controller, config) and the rendering layer (shaders,
// bloom, HUD).  Neither side includes the other's headers.
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
    bool hostGalaxyEnabled; // the galaxy itself, wowie :))))))))
    bool labEnabled;        // broad-line region label overlay, not a literal lab (unfortunately)
    bool cgmEnabled;        // circumgalactic medium, basically the galaxy's personal bubble of gas nobody talks about

    /*--------- Render settings (more steps = prettier, slower, warmer laptop) ---------*/
    int maxSteps;

    /*--------- Profile metadata (for HUD) ---------*/
    std::string profileName;

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
    // pro tip: if this drops below 10, you probably enabled every feature at once. Best of luck with that.

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
};
