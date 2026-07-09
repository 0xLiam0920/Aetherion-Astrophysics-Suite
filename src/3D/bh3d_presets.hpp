#pragma once
// ============================================================
// presets.hpp: Black hole profile presets with full SimConfig
// ============================================================

#include "config.hpp"
#include <string>
#include <array>
#include <vector>

// Body archetypes for systems orbiting the black hole.
// Mirrors the 2D simulator's GalaxyBodyType so that the same astrophysical
// systems (companion stars in HMXBs, S-cluster stars around Sgr A*, etc.)
// can be visualized in 3D.
enum class GalaxyBody3DType {
    Star,           // Generic star
    GasCloud,       // Gas/accretion clump
    StellarCluster, // Dense star cluster
    DwarfGalaxy,    // Dwarf-galaxy / satellite remnant
    NeutronStar,    // Compact remnant
    WhiteDwarf,     // Compact remnant
    CompanionStar   // Stellar binary companion
};

// A single orbiting member, in 3D-shader units (Rs).
struct GalaxyBody3D {
    GalaxyBody3DType type;
    float            semiMajorRs;   // Semi-major axis [Rs]
    float            eccentricity;
    float            inclination;   // Radians
    std::string      label;         // Optional display name for label-view mode (empty = derive from type)
};

// A complete black hole profile: config + metadata + default feature flags.
struct BlackHoleProfile {
    std::string       name;
    std::string       description;
    double            massSolar;      // Solar masses (informational)
    cfg::SimConfig    config;
    bool              defaultJets;    // Whether jets are visible by default
    bool              defaultBLR;     // Whether BLR is visible by default
    bool              defaultOrbBody; // Whether orbiting body is visible by default
    bool              defaultDoppler; // Whether Doppler beaming is on by default
    bool              defaultHostGalaxy; // Whether host galaxy is visible by default
    bool              defaultLAB;    // Whether Lyman-alpha Blob is visible by default
    bool              defaultCGM;    // Whether Circumgalactic Medium is visible by default

    // Orbiting bodies (stars / clouds / clusters / companions). When non-empty,
    // these REPLACE the legacy single config.orbital body for this profile.
    std::vector<GalaxyBody3D> galaxyBodies;

    // Barycentric binary: true when the companion mass is significant enough
    // that both the BH and companion visibly orbit the shared centre of mass.
    // Set for Gaia BH1/BH2/BH3, the astrometric wobble these systems were
    // discovered through IS the BH's displacement from the barycenter.
    bool   isBinaryWithBarycenter = false;
    double companionMassSolar     = 0.0;   // companion mass [Msun] for mass-ratio calc
};

// ============================================================
// Visual black hole merger: the kinds of secondary object that
// can be sent spiralling into the active primary, mirroring the
// 2D simulator's merger system (but rendered, not physically
// integrated). Order matters: indices are persisted in the menu
// and mapped to shader body types in bh3d_core.hpp.
// ============================================================
enum class MergerSecondaryKind3D {
    BlackHole = 0,   // dark horizon + photon ring (BODY_BLACKHOLE)
    NeutronStar,     // compact magenta remnant
    Pulsar,          // neutron star with lighthouse beam
    Star,            // main-sequence star (tidal-disruption aftermath)
    WhiteDwarf       // pale compact remnant
};

// A selectable secondary preset shown in the merger menu.
struct MergerSecondary3DPreset {
    const char*           label;       // Display name
    MergerSecondaryKind3D kind;        // Visual archetype
    double                massSolar;   // Secondary mass [Msun] (drives mass ratio + growth)
    const char*           blurb;       // One-line description for the menu
    // Name of the full BlackHoleProfile this secondary should visually mirror
    // (must match a profiles::allProfiles() name). When set, the inspiralling
    // secondary renders with that profile's accretion disk, spin, jets and
    // colours so it looks 1:1 with its standalone counterpart. nullptr → a
    // generic mass-derived disk (archetypes with no dedicated profile).
    const char*           profile = nullptr;
};

// Curated secondary objects, ordered light → heavy within each class.
// Named black holes mirror the 2D simulator's BH2D_PRESETS list so every
// object that can be merged in 2D is also a valid, working secondary in 3D.
inline constexpr MergerSecondary3DPreset MERGER_SECONDARY_3D_PRESETS[] = {
    // -- Generic black-hole archetypes --
    { "Stellar-mass BH",  MergerSecondaryKind3D::BlackHole,   10.0,    "A 10 M\u2609 stellar black hole spirals in." },
    { "Intermediate BH",  MergerSecondaryKind3D::BlackHole,   1.0e3,   "A 1000 M\u2609 IMBH, a near-equal merger for stellar primaries." },
    { "Supermassive BH",  MergerSecondaryKind3D::BlackHole,   1.0e8,   "A 10^8 M\u2609 SMBH, a galactic-scale coalescence." },

    // -- Compact / stellar secondaries --
    { "Neutron Star",     MergerSecondaryKind3D::NeutronStar, 1.4,     "A 1.4 M\u2609 neutron star plunges to the horizon." },
    { "Pulsar",           MergerSecondaryKind3D::Pulsar,      1.4,     "A spinning pulsar, lighthouse beam sweeping as it falls." },
    { "Main-sequence Star", MergerSecondaryKind3D::Star,      1.0,     "A Sun-like star, shredded into a tidal stream." },
    { "White Dwarf",      MergerSecondaryKind3D::WhiteDwarf,  0.6,     "A 0.6 M\u2609 white dwarf, disrupted on final approach." },

    // -─ Named stellar-mass black holes (ported  diectly from 2D) --
    { "GRO J1655-40",     MergerSecondaryKind3D::BlackHole,   6.3,     "Microquasar with relativistic jets, ~6.3 M\u2609, spirals in.", "GRO J1655-40" },
    { "A0620-00",         MergerSecondaryKind3D::BlackHole,   6.61,    "First confirmed stellar BH, quiescent ~6.6 M\u2609 binary.", "A0620-00" },
    { "V404 Cygni",       MergerSecondaryKind3D::BlackHole,   9.0,     "X-ray nova microquasar, ~9 M\u2609, plunges to the horizon.", "V404 Cygni" },
    { "Gaia BH2",         MergerSecondaryKind3D::BlackHole,   8.94,    "Dormant BH with red-giant companion, ~8.9 M\u2609.", "Gaia BH2" },
    { "Gaia BH1",         MergerSecondaryKind3D::BlackHole,   9.62,    "Nearest dormant BH, ~9.6 M\u2609, spirals into the primary.", "Gaia BH1" },
    { "Cygnus X-1",       MergerSecondaryKind3D::BlackHole,   21.2,    "Famous HMXB stellar black hole, ~21 M\u2609." },
    { "Gaia BH3",         MergerSecondaryKind3D::BlackHole,   32.7,    "Most massive nearby dormant BH, ~33 M\u2609.", "Gaia BH3" },
    { "LIGO GW150914",    MergerSecondaryKind3D::BlackHole,   62.0,    "First detected merger remnant, ~62 M\u2609, coalesces again." },

    // -- Exotic / intermediate-mass --
    { "Primordial BH",    MergerSecondaryKind3D::BlackHole,   1.0e-5,  "A tiny ~10^-5 M\u2609 primordial black hole grazes the horizon." },
    { "Intermediate-mass BH", MergerSecondaryKind3D::BlackHole, 1.0e4, "A hypothetical ~10^4 M\u2609 intermediate-mass black hole." },

    // -- Named supermassive black holes (ported from 2D) --
    // A secondary this massive shreds nearby stars and gas from a great distance.
    { "Sgr A* SMBH",      MergerSecondaryKind3D::BlackHole,   4.3e6,   "The Milky Way's 4.3\u00d710^6 M\u2609 nucleus spirals in.", "Sgr A*" },
    { "3C 273",           MergerSecondaryKind3D::BlackHole,   8.86e8,  "The first quasar identified, ~8.9\u00d710^8 M\u2609.", "3C 273" },
    { "M87* SMBH",        MergerSecondaryKind3D::BlackHole,   6.5e9,   "Virgo A's 6.5\u00d710^9 M\u2609 giant, the first BH ever imaged.", "M87*" },
    { "J0529-4351",       MergerSecondaryKind3D::BlackHole,   1.7e10,  "The most luminous quasar known, ~1.7\u00d710^10 M\u2609.", "J0529-4351" },
    { "NGC 1277",         MergerSecondaryKind3D::BlackHole,   1.7e10,  "Overmassive SMBH in a compact elliptical, ~1.7\u00d710^10 M\u2609.", "NGC 1277" },
    { "OJ 287",           MergerSecondaryKind3D::BlackHole,   1.8e10,  "SMBH binary system with orbital outbursts, ~1.8\u00d710^10 M\u2609.", "OJ 287" },
    { "TON 618 SMBH",     MergerSecondaryKind3D::BlackHole,   6.6e10,  "The 6.6\u00d710^10 M\u2609 ultramassive quasar, a titanic merger.", "TON 618" },
};
inline constexpr int NUM_MERGER_SECONDARY_3D_PRESETS =
    (int)(sizeof(MERGER_SECONDARY_3D_PRESETS) / sizeof(MERGER_SECONDARY_3D_PRESETS[0]));

namespace profiles {

/*--------- Type → visual archetype mapping ---------*/
// Returns a fully-populated OrbitalConfig with type-appropriate radius/color
// for use with the OrbitalBody Keplerian integrator.
inline cfg::OrbitalConfig makeOrbitalBody(const GalaxyBody3D& b) {
    cfg::OrbitalConfig oc;
    oc.semiMajor    = b.semiMajorRs;
    oc.eccentricity = b.eccentricity;
    oc.inclination  = b.inclination;
    oc.bodyType     = static_cast<int>(b.type);
    switch (b.type) {
        case GalaxyBody3DType::Star:
            oc.bodyRadius = 0.8f;
            oc.bodyColor  = glm::vec3(1.00f, 0.93f, 0.78f);  // Yellow-white star
            break;
        case GalaxyBody3DType::CompanionStar:
            oc.bodyRadius = 1.4f;
            oc.bodyColor  = glm::vec3(1.00f, 0.85f, 0.55f);  // Warm K/G companion
            break;
        case GalaxyBody3DType::GasCloud:
            oc.bodyRadius = 1.2f;
            oc.bodyColor  = glm::vec3(0.95f, 0.55f, 0.30f);  // Hot ionised gas
            break;
        case GalaxyBody3DType::StellarCluster:
            oc.bodyRadius = 2.0f;
            oc.bodyColor  = glm::vec3(0.90f, 0.80f, 0.65f);  // Mixed stellar pop
            break;
        case GalaxyBody3DType::DwarfGalaxy:
            oc.bodyRadius = 2.5f;
            oc.bodyColor  = glm::vec3(0.70f, 0.72f, 0.95f);  // Distant bluish smudge
            break;
        case GalaxyBody3DType::NeutronStar:
            oc.bodyRadius = 0.40f;
            // Vivid magenta, physically a pulsar is blue-white, but here we
            // tint it for at-a-glance identification against nearby S-stars.
            oc.bodyColor  = glm::vec3(1.00f, 0.35f, 0.90f);
            break;
        case GalaxyBody3DType::WhiteDwarf:
            oc.bodyRadius = 0.50f;
            oc.bodyColor  = glm::vec3(1.00f, 1.00f, 0.95f);  // Pale white
            break;
    }
    return oc;
}


/*--------- TON 618, Most massive known quasar ---------*/
inline BlackHoleProfile ton618() {
    BlackHoleProfile p;
    p.name        = "TON 618";
    p.description = "Most massive known quasar (~6.6e10 Msun)";
    p.massSolar   = 6.6e10;
    p.config      = cfg::defaultConfig();
    // TON 618: hot blue-white inner disk fading to warm orange outer
    p.config.disk.peakTemp           = 30000.0f;
    p.config.disk.displayTempInner   = 5500.0f;  // Warm white
    p.config.disk.displayTempOuter   = 2200.0f;  // Deep orange
    p.config.disk.saturationBoostInner = 2.0f;
    p.config.disk.saturationBoostOuter = 1.1f;
    // Orbiting body at larger distance for quasar scale
    p.config.orbital.semiMajor = 100.0f;  // Further out for quasar
    p.config.orbital.bodyRadius = 1.5f;   // Larger body for visibility
    p.config.blr.strength = 1.0f;  // TON 618: strongest observed BLR
    p.defaultJets    = true;
    p.defaultBLR     = true;
    p.defaultOrbBody = true;
    p.defaultDoppler = true;
    p.defaultHostGalaxy = true;
    p.defaultLAB = true;
    p.defaultCGM = false;  // Diffuse fog off by default
    // 6-body galactic-nuclear environment (mirrors 2D TON618_BODIES)
    p.galaxyBodies = {
        {GalaxyBody3DType::GasCloud,        8.0f,  0.15f, -0.40f},
        {GalaxyBody3DType::Star,           15.0f,  0.55f, -0.22f},
        {GalaxyBody3DType::GasCloud,       25.0f,  0.30f, -0.04f},
        {GalaxyBody3DType::Star,           12.0f,  0.80f,  0.14f},
        {GalaxyBody3DType::Star,           50.0f,  0.20f,  0.32f},
        {GalaxyBody3DType::StellarCluster, 80.0f,  0.10f,  0.50f}
    };
    return p;
}

/*--------- Sagittarius A*, Milky Way center ---------*/
// Low-luminosity, radiatively inefficient accretion flow (RIAF).
// Weak/absent jets, no broad-line region, small dim disk.
inline BlackHoleProfile sgrAstar() {
    BlackHoleProfile p;
    p.name        = "Sgr A*";
    p.description = "Milky Way center (~4.3e6 Msun), dim, quiet";
    p.massSolar   = 4.3e6;

    cfg::SimConfig c;
    // Moderate spin, current EHT constraints suggest ~0.5–0.9
    c.blackHole.spinParameter = 0.5f;
    c.blackHole.radius        = 1.0f;

    // Small, dim accretion flow, RIAF geometry
    c.disk.innerRadius   = 3.0f;   // Larger ISCO (lower spin)
    c.disk.outerRadius   = 10.0f;  // Compact accretion flow
    c.disk.halfThickness = 0.04f;  // Thicker, puffy RIAF

    // Weak jets (barely visible if toggled on)
    c.jet.radius = 0.15f;
    c.jet.length = 10.0f;
    c.jet.color  = glm::vec3(0.15f, 0.5f, 0.8f);  // Faint blue

    // No significant BLR
    c.blr.innerRadius = 8.0f;
    c.blr.outerRadius = 14.0f;
    c.blr.thickness   = 2.0f;
    c.blr.strength    = 0.0f;  // Sgr A*: RIAF, no detectable BLR

    // Camera closer, smaller object
    c.camera.initialPos = glm::vec3(0.0f, 4.0f, 18.0f);

    // Subtle bloom (dim source)
    c.bloom.threshold  = 1.5f;
    c.bloom.intensity  = 0.35f;
    c.bloom.exposure   = 0.9f;

    // Sgr A*: reddish-orange/yellow based on EHT images
    c.disk.peakTemp           = 8000.0f;   // Lower temp RIAF
    c.disk.displayTempInner   = 3200.0f;   // Reddish-orange inner
    c.disk.displayTempOuter   = 1600.0f;   // Deep red outer
    c.disk.saturationBoostInner = 2.8f;    // Strong warm saturation
    c.disk.saturationBoostOuter = 1.8f;    // Keep outer vivid

    p.config = c;
    p.defaultJets    = false;  // Sgr A* has no prominent jets
    p.defaultBLR     = false;  // No BLR in a low-luminosity AGN
    p.defaultOrbBody = true;   // S-cluster + circumnuclear gas orbit Sgr A*
    p.defaultDoppler = true;
    p.defaultHostGalaxy = false;
    p.defaultLAB = false;
    p.defaultCGM = false;
    // S-cluster stars + circumnuclear gas (mirrors 2D SGRA_BODIES) + the
    // 8.19 ms pulsar discovered in tight orbit around Sgr A*.
    p.galaxyBodies = {
        {GalaxyBody3DType::Star,        20.0f, 0.88f, -0.35f, "S2"},
        {GalaxyBody3DType::Star,        10.0f, 0.95f, -0.15f, "S14"},
        {GalaxyBody3DType::Star,        40.0f, 0.30f,  0.05f, "IRS 16"},
        {GalaxyBody3DType::GasCloud,    60.0f, 0.15f,  0.25f, "Circumnuclear Gas"},
        {GalaxyBody3DType::Star,        30.0f, 0.50f,  0.45f, "S-cluster"},
        // BLPSR, 8.19 ms pulsar in tight relativistic orbit around Sgr A*.
        // Compact, fast-spinning recycled neutron star; emits a coherent
        // radio/X-ray beam every 8.19 ms (~122 Hz). Coloured magenta in-sim
        // purely as a visual identifier against the surrounding S-cluster.
        {GalaxyBody3DType::NeutronStar,  6.5f, 0.72f, -0.55f, "BLPSR"}
    };
    return p;
}

/*--------- 3C 273, First identified quasar ---------*/
// Powerful relativistic jet, bright accretion disk, active BLR.
inline BlackHoleProfile qso3c273() {
    BlackHoleProfile p;
    p.name        = "3C 273";
    p.description = "First quasar identified (~8.9e8 Msun), bright jet";
    p.massSolar   = 8.86e8;

    cfg::SimConfig c;
    // High spin, jet-producing quasar
    c.blackHole.spinParameter = 0.9f;
    c.blackHole.radius        = 1.0f;

    // Bright, extended accretion disk
    c.disk.innerRadius   = 1.8f;   // Close to ISCO at high spin
    c.disk.outerRadius   = 22.0f;
    c.disk.halfThickness = 0.015f;

    // Powerful relativistic jet (3C 273 is famous for its jet)
    c.jet.radius = 0.45f;
    c.jet.length = 55.0f;  // Enormous jet
    c.jet.color  = glm::vec3(0.3f, 0.85f, 1.5f);  // Bright synchrotron blue

    // Active broad-line region
    c.blr.innerRadius = 16.0f;
    c.blr.outerRadius = 38.0f;
    c.blr.thickness   = 7.0f;
    c.blr.strength    = 0.75f;  // 3C 273: bright quasar BLR

    c.camera.initialPos = glm::vec3(0.0f, 7.0f, 28.0f);

    // Bright bloom, luminous quasar
    c.bloom.threshold  = 1.0f;
    c.bloom.intensity  = 0.60f;
    c.bloom.exposure   = 1.1f;

    // 3C 273: bright blue-white quasar disk
    c.disk.peakTemp           = 35000.0f;
    c.disk.displayTempInner   = 6500.0f;   // Hot blue-white core
    c.disk.displayTempOuter   = 2500.0f;   // Warm yellow outer
    c.disk.saturationBoostInner = 1.8f;
    c.disk.saturationBoostOuter = 1.2f;

    p.config = c;
    p.defaultJets    = true;   // Famous jet
    p.defaultBLR     = true;
    p.defaultOrbBody = true;
    p.defaultDoppler = true;
    p.defaultHostGalaxy = true;
    p.defaultLAB = true;
    p.defaultCGM = false;
    // Quasar environment (mirrors 2D QSO3C273_BODIES)
    p.galaxyBodies = {
        {GalaxyBody3DType::GasCloud,    10.0f, 0.20f, -0.30f},  // Inner jet-base cloud
        {GalaxyBody3DType::GasCloud,    20.0f, 0.40f, -0.10f},  // BLR cloud
        {GalaxyBody3DType::DwarfGalaxy, 70.0f, 0.25f,  0.10f},  // Stripped dwarf remnant
        {GalaxyBody3DType::Star,        35.0f, 0.60f,  0.30f}   // Close stellar orbit
    };
    return p;
}

/*--------- J0529-4351, Most luminous quasar known ---------*/
// Extreme accretion rate, enormous disk, blazing luminosity.
inline BlackHoleProfile j0529() {
    BlackHoleProfile p;
    p.name        = "J0529-4351";
    p.description = "Most luminous quasar (~1.7e10 Msun), extreme disk";
    p.massSolar   = 1.7e10;

    cfg::SimConfig c;
    // Very high spin, extreme accretion
    c.blackHole.spinParameter = 0.95f;
    c.blackHole.radius        = 1.0f;

    // Massive, bright accretion disk, largest known
    c.disk.innerRadius   = 1.5f;   // Very close ISCO at near-maximal spin
    c.disk.outerRadius   = 30.0f;  // Enormous disk
    c.disk.halfThickness = 0.01f;  // Thin, efficient disk

    // Strong jets
    c.jet.radius = 0.5f;
    c.jet.length = 50.0f;
    c.jet.color  = glm::vec3(0.35f, 0.95f, 1.6f);  // Intense synchrotron

    // Extensive BLR
    c.blr.innerRadius = 20.0f;
    c.blr.outerRadius = 42.0f;
    c.blr.thickness   = 8.0f;
    c.blr.strength    = 0.90f;  // J0529-4351: hyperluminous, near-maximal BLR

    // Pull camera back, big object
    c.camera.initialPos = glm::vec3(0.0f, 8.0f, 32.0f);

    // Intense bloom, most luminous object in the universe
    c.bloom.threshold  = 0.8f;
    c.bloom.intensity  = 0.75f;
    c.bloom.exposure   = 1.2f;
    c.bloom.softKnee   = 0.7f;

    // J0529-4351: blazing white-blue disk (extreme luminosity)
    c.disk.peakTemp           = 40000.0f;
    c.disk.displayTempInner   = 7000.0f;   // Very hot blue-white
    c.disk.displayTempOuter   = 2800.0f;   // Yellow-orange fringe
    c.disk.saturationBoostInner = 1.6f;
    c.disk.saturationBoostOuter = 1.3f;

    p.config = c;
    p.defaultJets    = true;
    p.defaultBLR     = true;
    p.defaultOrbBody = true;
    p.defaultDoppler = true;
    p.defaultHostGalaxy = true;
    p.defaultLAB = true;
    p.defaultCGM = false;
    // Hyperluminous quasar environment (mirrors 2D J0529_BODIES)
    p.galaxyBodies = {
        {GalaxyBody3DType::GasCloud,        9.0f, 0.10f, -0.40f},  // Fast accretion blob
        {GalaxyBody3DType::GasCloud,       18.0f, 0.35f, -0.20f},  // UV-bright clump
        {GalaxyBody3DType::Star,           14.0f, 0.70f,  0.00f},  // Tidally disrupting star
        {GalaxyBody3DType::GasCloud,       45.0f, 0.20f,  0.20f},  // Outer gas stream
        {GalaxyBody3DType::StellarCluster, 65.0f, 0.15f,  0.40f}   // Infalling cluster
    };
    return p;
}

/*--------- Gaia BH1, Nearest known dormant black hole ---------*/
// Discovered via astrometric wobble of its G-dwarf companion (Gaia DR3).
// Dormant: NO detectable accretion disk. Disk params represent theoretical
// Bondi-Hoyle capture of stellar wind at ~10^-8 L_Edd, far below any
// observable threshold. The companion star is the entire visual story here.
inline BlackHoleProfile gaiaBH1() {
    BlackHoleProfile p;
    p.name        = "Gaia BH1";
    p.description = "Nearest dormant BH (~9.6 Msun), discovered via astrometry, no disk";
    p.massSolar   = 9.62;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.3f;
    c.blackHole.radius        = 1.0f;
    // No accretion disk, dormant black hole, discovered by astrometry only.
    // outerRadius = 0 makes the shader condition (r < diskOuterRadius) always
    // false, suppressing all disk rendering passes.
    c.disk.innerRadius   = 3.0f;
    c.disk.outerRadius   = 0.0f;    // Zero, disables disk entirely in shader
    c.disk.halfThickness = 0.0f;
    c.disk.peakTemp             = 1000.0f;
    c.disk.displayTempInner     = 1000.0f;
    c.disk.displayTempOuter     = 1000.0f;
    c.disk.saturationBoostInner = 0.0f;
    c.disk.saturationBoostOuter = 0.0f;
    c.jet.radius = 0.08f; c.jet.length = 4.0f;
    c.jet.color  = glm::vec3(0.1f, 0.3f, 0.6f);
    c.orbital.semiMajor  = 28.0f;
    c.orbital.bodyRadius = 1.2f;    // G-dwarf companion, the visual focus
    c.camera.initialPos  = glm::vec3(0.0f, 3.5f, 15.0f);
    c.bloom.threshold    = 4.0f;    // Nothing blooms at this luminosity
    c.bloom.intensity    = 0.04f;
    c.bloom.exposure     = 0.40f;
    p.config = c;
    p.defaultJets        = false;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // G-dwarf binary companion. Real measured eccentricity e ≈ 0.451
    // (Gaia DR3), the visual focus of this dormant BH system.
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 28.0f, 0.45f, 0.10f}
    };
    p.isBinaryWithBarycenter = true;
    p.companionMassSolar     = 0.93;   // G-dwarf, Gaia DR3
    return p;
}

/*--------- Gaia BH2, Dormant BH with red giant companion ---------*/
// Discovered via astrometric wobble of a red giant companion ~1.16 kpc away.
// Dormant: NO detectable accretion disk. The red giant overfills its Roche
// lobe slightly, but the resulting RIAF luminosity is still undetectable.
// The oversized companion body is the sole visual identifier for this system.
inline BlackHoleProfile gaiaBH2() {
    BlackHoleProfile p;
    p.name        = "Gaia BH2";
    p.description = "Dormant BH (~8.9 Msun), discovered via astrometry, red giant companion";
    p.massSolar   = 8.94;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.2f;
    c.blackHole.radius        = 1.0f;
    // No accretion disk, dormant black hole, discovered by astrometry only.
    // outerRadius = 0 makes the shader condition (r < diskOuterRadius) always
    // false, suppressing all disk rendering passes.
    c.disk.innerRadius   = 3.0f;
    c.disk.outerRadius   = 0.0f;    // Zero, disables disk entirely in shader
    c.disk.halfThickness = 0.0f;
    c.disk.peakTemp             = 1000.0f;
    c.disk.displayTempInner     = 1000.0f;
    c.disk.displayTempOuter     = 1000.0f;
    c.disk.saturationBoostInner = 0.0f;
    c.disk.saturationBoostOuter = 0.0f;
    c.jet.radius = 0.08f; c.jet.length = 4.0f;
    c.jet.color  = glm::vec3(0.1f, 0.3f, 0.6f);
    c.orbital.semiMajor  = 30.0f;
    c.orbital.bodyRadius = 2.4f;    // Red giant is noticeably larger, the visual focus
    c.camera.initialPos  = glm::vec3(0.0f, 3.5f, 15.0f);
    c.bloom.threshold    = 4.0f;
    c.bloom.intensity    = 0.04f;
    c.bloom.exposure     = 0.40f;
    p.config = c;
    p.defaultJets        = false;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // Red giant companion. Real measured eccentricity e ≈ 0.518 (Gaia DR3).
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 30.0f, 0.52f, 0.12f}
    };
    p.isBinaryWithBarycenter = true;
    p.companionMassSolar     = 1.07;   // Red giant, Gaia DR3
    return p;
}

/*--------- Gaia BH3, Most massive nearby dormant BH ---------*/
// Discovered via astrometric wobble of a metal-poor giant companion.
// ~33 Msun, well above typical stellar BH mass, yet completely dormant.
// No detectable X-ray or optical emission from any accretion disk.
// Disk params: theoretical RIAF trace only; renders invisible by design.
inline BlackHoleProfile gaiaBH3() {
    BlackHoleProfile p;
    p.name        = "Gaia BH3";
    p.description = "Most massive nearby dormant BH (~33 Msun), discovered via astrometry, no disk";
    p.massSolar   = 32.7;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.45f;
    c.blackHole.radius        = 1.0f;
    // No accretion disk, dormant black hole, discovered by astrometry only.
    // outerRadius = 0 makes the shader condition (r < diskOuterRadius) always
    // false, suppressing all disk rendering passes.
    c.disk.innerRadius   = 2.8f;
    c.disk.outerRadius   = 0.0f;    // Zero, disables disk entirely in shader
    c.disk.halfThickness = 0.0f;
    c.disk.peakTemp             = 1000.0f;
    c.disk.displayTempInner     = 1000.0f;
    c.disk.displayTempOuter     = 1000.0f;
    c.disk.saturationBoostInner = 0.0f;
    c.disk.saturationBoostOuter = 0.0f;
    c.jet.radius = 0.1f; c.jet.length = 5.0f;
    c.jet.color  = glm::vec3(0.1f, 0.35f, 0.7f);
    c.orbital.semiMajor  = 22.0f;
    c.orbital.bodyRadius = 1.6f;    // Metal-poor giant companion
    c.camera.initialPos  = glm::vec3(0.0f, 4.0f, 22.0f);
    c.bloom.threshold    = 4.0f;
    c.bloom.intensity    = 0.04f;
    c.bloom.exposure     = 0.40f;
    p.config = c;
    p.defaultJets        = false;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // Metal-poor giant companion. Real measured eccentricity e ≈ 0.729
    // (Gaia DR3, 2024), a dramatic wide-binary swing from periapsis to apoapsis.
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 38.0f, 0.73f, 0.08f}
    };
    p.isBinaryWithBarycenter = true;
    p.companionMassSolar     = 0.76;   // Metal-poor giant, Gaia DR3 2024
    return p;
}

/*--------- V404 Cygni, X-ray nova microquasar ---------*/
// Represents the 2015 outburst: brightest stellar BH transient in decades.
// Hot bright disk. kelvinToRGB(7500) → near-white with yellow tinge.
// Well-distinguished from quiescent Gaia BHs by temperature and brightness.
inline BlackHoleProfile v404Cyg() {
    BlackHoleProfile p;
    p.name        = "V404 Cygni";
    p.description = "Triple system: ~9 Msun BH + K-giant donor + distant tertiary (V404 Cyg C)";
    p.massSolar   = 9.0;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.5f;
    c.blackHole.radius        = 1.0f;
    c.disk.innerRadius   = 2.5f;
    c.disk.outerRadius   = 12.0f;
    c.disk.halfThickness = 0.03f;   // Thin, efficient accretion (outburst)
    // kelvinToRGB(7500) → bright pale-yellow/almost-white inner core.
    // kelvinToRGB(3800) → orange middle ring, classic X-ray nova gradient.
    c.disk.peakTemp             = 30000.0f;
    c.disk.displayTempInner     = 7500.0f;   // Hot pale-yellow inner ring
    c.disk.displayTempOuter     = 3800.0f;   // Orange outer disk
    c.disk.saturationBoostInner = 2.0f;
    c.disk.saturationBoostOuter = 3.0f;      // Punchy orange outer
    c.jet.radius = 0.20f;
    c.jet.length = 14.0f;
    c.jet.color  = glm::vec3(0.2f, 0.65f, 1.3f);
    c.orbital.semiMajor  = 26.0f;
    c.orbital.bodyRadius = 1.4f;
    c.camera.initialPos  = glm::vec3(0.0f, 3.5f, 16.0f);
    c.bloom.threshold    = 1.35f;
    c.bloom.intensity    = 0.45f;
    c.bloom.exposure     = 1.0f;
    p.config = c;
    p.defaultJets        = true;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // Triple system:
    //   • Donor (K-giant, ~0.7 Msun) on the close 6.47-day orbit being stripped
    //     into the disk, modelled as the inner CompanionStar.
    //   • V404 Cyg C, distant tertiary (~3500 AU, ~70 000-yr period, evolved
    //     star / future red giant) discovered 2024. Compressed to a far orbit
    //     in simulation units so it stays visible alongside the inner pair.
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 26.0f,  0.034f,  0.05f, "K-giant donor"},
        {GalaxyBody3DType::Star,          95.0f,  0.30f,  -0.18f, "V404 Cyg C"}
    };
    return p;
}

/*--------- A0620-00, First confirmed stellar BH ---------*/
// Quiescent LMXB, archetypal dormant system; K-dwarf donor, very low accretion.
// Faintest disk of any preset. kelvinToRGB(3000) → deep red-orange ember glow.
// Distinctive as the "dimmest", users can see it's barely accreting.
inline BlackHoleProfile a062000() {
    BlackHoleProfile p;
    p.name        = "A0620-00";
    p.description = "First confirmed stellar BH (~6.6 Msun), quiescent K-dwarf";
    p.massSolar   = 6.61;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.12f;  // Low spin, lowest Novikov-Thorne flux
    c.blackHole.radius        = 1.0f;
    c.disk.innerRadius   = 3.2f;
    c.disk.outerRadius   = 8.0f;
    c.disk.halfThickness = 0.065f;  // Puffy ADAF
    // kelvinToRGB(3000) → vivid deep red-orange (the ember look).
    // kelvinToRGB(2100) → smouldering red edge.
    c.disk.peakTemp             = 8000.0f;
    c.disk.displayTempInner     = 3000.0f;   // Deep red-orange ember
    c.disk.displayTempOuter     = 2100.0f;   // Smouldering red outer
    c.disk.saturationBoostInner = 3.6f;      // Max saturation to make the ember glow pop
    c.disk.saturationBoostOuter = 2.8f;
    c.jet.radius = 0.08f; c.jet.length = 4.0f;
    c.jet.color  = glm::vec3(0.1f, 0.25f, 0.5f);
    c.orbital.semiMajor  = 24.0f;
    c.orbital.bodyRadius = 1.1f;
    c.camera.initialPos  = glm::vec3(0.0f, 3.0f, 13.0f);
    c.bloom.threshold    = 1.95f;
    c.bloom.intensity    = 0.22f;
    c.bloom.exposure     = 0.72f;
    p.config = c;
    p.defaultJets        = false;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // K-dwarf companion (circular orbit)
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 24.0f, 0.0f, 0.04f}
    };
    return p;
}

/*--------- GRO J1655-40, Microquasar with relativistic jets ---------*/
// Best-measured spin in the stellar-mass regime (a* ~ 0.7).
// Hot, efficient thin disk during outburst.
// kelvinToRGB(9500) → white with a perceptible blue-violet tinge, the hottest
// stellar BH in the set. The jet color reinforces the blue identity.
inline BlackHoleProfile groJ165540() {
    BlackHoleProfile p;
    p.name        = "GRO J1655-40";
    p.description = "Microquasar with relativistic jets (~6.3 Msun), F-star companion";
    p.massSolar   = 6.3;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.7f;
    c.blackHole.radius        = 1.0f;
    c.disk.innerRadius   = 1.9f;    // Smaller ISCO at a*=0.7
    c.disk.outerRadius   = 10.0f;
    c.disk.halfThickness = 0.025f;  // Very thin, efficient disk
    // kelvinToRGB(9500) → cold white with faint blue cast, clearly hotter than V404.
    // kelvinToRGB(4500) → warm yellow outer; keeps the gradient interesting.
    c.disk.peakTemp             = 40000.0f;
    c.disk.displayTempInner     = 9500.0f;   // Near-white with blue cast
    c.disk.displayTempOuter     = 4500.0f;   // Yellow-orange middle/outer
    c.disk.saturationBoostInner = 1.5f;      // Low, blue-white desaturates naturally
    c.disk.saturationBoostOuter = 2.6f;      // Orange outer pops against the white core
    c.jet.radius = 0.24f;
    c.jet.length = 18.0f;
    c.jet.color  = glm::vec3(0.3f, 0.8f, 1.6f);   // Vivid synchrotron blue
    c.orbital.semiMajor  = 22.0f;
    c.orbital.bodyRadius = 1.3f;
    c.camera.initialPos  = glm::vec3(0.0f, 3.5f, 16.0f);
    c.bloom.threshold    = 1.25f;
    c.bloom.intensity    = 0.50f;
    c.bloom.exposure     = 1.05f;
    p.config = c;
    p.defaultJets        = true;
    p.defaultBLR         = false;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = false;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // F-star companion (circular orbit)
    p.galaxyBodies = {
        {GalaxyBody3DType::CompanionStar, 22.0f, 0.0f, 0.06f}
    };
    return p;
}

/*--------- NGC 1277, Overmassive SMBH in compact elliptical ---------*/
// ~14% of host-galaxy bulge mass, 10× above the M-σ relation.
// AGN-class accretion; should visually resemble a compact 3C 273 but with
// a hotter, more saturated disk to distinguish it.
// kelvinToRGB(6800) → warm-white with slight gold tinge.
inline BlackHoleProfile ngc1277() {
    BlackHoleProfile p;
    p.name        = "NGC 1277";
    p.description = "Overmassive SMBH in compact elliptical (~1.7e10 Msun)";
    p.massSolar   = 1.7e10;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.65f;
    c.blackHole.radius        = 1.0f;
    c.disk.innerRadius   = 1.9f;
    c.disk.outerRadius   = 20.0f;
    c.disk.halfThickness = 0.014f;
    // kelvinToRGB(6800) → gold-white inner. kelvinToRGB(3000) → orange outer.
    // The gold-to-orange gradient makes it feel distinctly different from
    // the bluer OJ 287 or the warm-white quasars.
    c.disk.peakTemp             = 28000.0f;
    c.disk.displayTempInner     = 6800.0f;   // Gold-white inner AGN disk
    c.disk.displayTempOuter     = 3000.0f;   // Orange outer, vivid contrast
    c.disk.saturationBoostInner = 2.2f;
    c.disk.saturationBoostOuter = 3.2f;      // Make the orange outer edge striking
    c.jet.radius = 0.30f;
    c.jet.length = 28.0f;
    c.jet.color  = glm::vec3(0.22f, 0.70f, 1.35f);
    c.blr.innerRadius = 13.0f;
    c.blr.outerRadius = 28.0f;
    c.blr.thickness   = 5.0f;
    c.blr.strength    = 0.55f;  // NGC 1277: compact SMBH, moderate BLR
    c.orbital.semiMajor  = 80.0f;
    c.orbital.bodyRadius = 1.5f;
    c.camera.initialPos  = glm::vec3(0.0f, 6.0f, 24.0f);
    c.bloom.threshold    = 1.05f;
    c.bloom.intensity    = 0.58f;
    c.bloom.exposure     = 1.10f;
    p.config = c;
    p.defaultJets        = false;
    p.defaultBLR         = true;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = true;
    p.defaultLAB         = false;
    p.defaultCGM         = false;
    // NGC 1277 nuclear environment (mirrors 2D NGC1277_BODIES)
    p.galaxyBodies = {
        {GalaxyBody3DType::Star,           12.0f, 0.30f, -0.25f},  // Inner stellar orbit A
        {GalaxyBody3DType::Star,           20.0f, 0.45f, -0.05f},  // Inner stellar orbit B
        {GalaxyBody3DType::GasCloud,       30.0f, 0.20f,  0.15f},  // Gas cloud
        {GalaxyBody3DType::StellarCluster, 75.0f, 0.10f,  0.35f}   // Stellar cluster
    };
    return p;
}

/*--------- OJ 287, SMBH binary with optical outbursts ---------*/
// Primary BH ~1.8×10¹⁰ Msun; secondary (~1.5×10⁸ Msun) punches through the
// disk every ~12 years, triggering UV/optical flares.
// Blazar-class: jet pointed near line of sight → Doppler-boosted bright jet.
// kelvinToRGB(11000) → distinctly icy-blue-white; hottest and brightest of the set.
inline BlackHoleProfile oj287() {
    BlackHoleProfile p;
    p.name        = "OJ 287";
    p.description = "SMBH binary system (~1.8e10 Msun primary), blazar, optical outbursts";
    p.massSolar   = 1.8e10;

    cfg::SimConfig c;
    c.blackHole.spinParameter = 0.82f;
    c.blackHole.radius        = 1.0f;
    c.disk.innerRadius   = 1.75f;   // Very close ISCO at a*=0.82
    c.disk.outerRadius   = 26.0f;
    c.disk.halfThickness = 0.011f;  // Extremely thin Shakura-Sunyaev disk
    // kelvinToRGB(11000) → icy blue-white, immediately distinguishable from all others.
    // kelvinToRGB(4200) → warm yellow-orange outer, strong thermal gradient.
    c.disk.peakTemp             = 38000.0f;
    c.disk.displayTempInner     = 11000.0f;  // Icy blue-white blazar inner disk
    c.disk.displayTempOuter     = 4200.0f;   // Yellow-orange outer, broad gradient
    c.disk.saturationBoostInner = 1.4f;      // Cool, blue-white is naturally pale
    c.disk.saturationBoostOuter = 2.8f;      // Pop the yellow outer ring
    c.jet.radius = 0.45f;
    c.jet.length = 48.0f;
    c.jet.color  = glm::vec3(0.35f, 0.85f, 1.65f);   // Intense electric-blue jet
    c.blr.innerRadius = 17.0f;
    c.blr.outerRadius = 36.0f;
    c.blr.thickness   = 7.0f;
    c.blr.strength    = 0.65f;  // OJ 287: blazar-class BLR, variable
    c.orbital.semiMajor  = 14.0f;    // Secondary SMBH in tight orbit
    c.orbital.bodyRadius = 2.2f;
    c.camera.initialPos  = glm::vec3(0.0f, 7.0f, 27.0f);
    c.bloom.threshold    = 0.88f;
    c.bloom.intensity    = 0.68f;
    c.bloom.exposure     = 1.15f;
    c.bloom.softKnee     = 0.62f;
    p.config = c;
    p.defaultJets        = true;
    p.defaultBLR         = true;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = true;
    p.defaultLAB         = true;
    p.defaultCGM         = false;
    // OJ 287 binary system (mirrors 2D OJ287_BODIES)
    // Secondary SMBH represented as a bright star, the famous ~12yr periodic
    // disk-piercing companion that triggers optical outbursts.
    p.galaxyBodies = {
        {GalaxyBody3DType::Star,     15.0f, 0.66f, -0.20f},  // Secondary SMBH (~1.5e8 Msun)
        {GalaxyBody3DType::GasCloud, 22.0f, 0.30f,  0.05f},  // Accretion blob
        {GalaxyBody3DType::GasCloud, 40.0f, 0.20f,  0.30f}   // Outer gas cloud
    };
    return p;
}

/*--------- Phoenix A, Largest known black hole ---------*/
// Ultramassive BH in the BCG of the Phoenix Cluster (SPT-CLJ2344-4243), z = 0.597.
// Mass ~1.0×10¹¹ Msun, the largest known black hole, ~1.52× TON 618.
// Schwarzschild radius Rs ≈ 1970 AU ≈ 6.4 light-days (~50× Pluto's orbital radius).
// Spin unconstrained; ultramassive BHs grown partly via mergers of randomly-aligned
// spins tend toward lower net spin → conservative a* = 0.35.
// The Phoenix Cluster hosts one of the highest BCG star-formation rates ever observed
// (~740 Msun/yr), driven by a prodigious cooling flow that also feeds the AGN.
// Jet cavities ~200 kpc across in the cluster ICM confirm past extreme jet activity.
// Disk is slightly cooler than TON 618: T_inner ∝ M^{-1/4} at comparable Eddington
// fraction → amber-gold (4800 K) vs TON 618's warm-white (5500 K).
inline BlackHoleProfile phoenixA() {
    BlackHoleProfile p;
    p.name        = "Phoenix A";
    p.description = "Largest known BH (~1.0e11 Msun), Rs~1970 AU, Phoenix Cluster BCG";
    p.massSolar   = 1.0e11;

    cfg::SimConfig c;
    // Suppressed spin: mergers of randomly-oriented progenitors tend to cancel.
    c.blackHole.spinParameter = 0.35f;
    c.blackHole.radius        = 1.0f;

    // At a*=0.35, ISCO ≈ 5.3 Rs. Larger outer disk than TON 618 to convey extra scale.
    c.disk.innerRadius   = 5.0f;
    c.disk.outerRadius   = 30.0f;
    c.disk.halfThickness = 0.010f;  // Thin efficient Shakura-Sunyaev disk
    // kelvinToRGB(4800) → amber-gold (cooler than TON 618's 5500 K warm-white).
    // kelvinToRGB(1900) → deep crimson-orange outer rim, visually distinct from all others.
    c.disk.peakTemp             = 26000.0f;
    c.disk.displayTempInner     = 4800.0f;   // Amber-gold AGN inner disk
    c.disk.displayTempOuter     = 1900.0f;   // Deep crimson-orange outer rim
    c.disk.saturationBoostInner = 2.4f;      // Rich amber saturation
    c.disk.saturationBoostOuter = 1.6f;      // Warm crimson glow at the rim

    // Powerful relativistic jets, X-ray cavities ~200 kpc wide confirmed.
    // Longer than TON 618 jets; wide radius befitting the extreme jet power.
    c.jet.radius = 0.50f;
    c.jet.length = 55.0f;
    c.jet.color  = glm::vec3(0.30f, 0.72f, 1.50f);   // Synchrotron blue-white

    // Broad-line region: powerful cooling-flow-fed AGN
    c.blr.innerRadius = 16.0f;
    c.blr.outerRadius = 35.0f;
    c.blr.thickness   = 7.0f;
    c.blr.strength    = 0.70f;  // Phoenix A: ultramassive, luminous BLR

    // Legacy single orbital body (overridden by galaxyBodies below)
    c.orbital.semiMajor  = 110.0f;
    c.orbital.bodyRadius = 1.5f;

    // Camera pulled back to encompass the larger disk
    c.camera.initialPos = glm::vec3(0.0f, 8.0f, 34.0f);

    // Strong bloom, luminous cooling-flow AGN
    c.bloom.threshold = 0.92f;
    c.bloom.intensity = 0.70f;
    c.bloom.exposure  = 1.15f;
    c.bloom.softKnee  = 0.58f;

    p.config = c;
    p.defaultJets        = true;
    p.defaultBLR         = true;
    p.defaultOrbBody     = true;
    p.defaultDoppler     = true;
    p.defaultHostGalaxy  = true;
    p.defaultLAB         = true;
    p.defaultCGM         = false;

    // Phoenix Cluster BCG nuclear environment:
    //   Inner condensations: cooling-flow filaments that feed the AGN
    //   Mid-distance gas:    ionised nebula in the BCG core
    //   Nuclear star:        analogue to S-stars near a galactic nucleus
    //   Inner/outer clusters: trace the massive BCG stellar envelope
    //   Infalling satellite: BCG cannibalising a satellite galaxy (common in clusters)
    // Orbital radii are deliberately large to reinforce the extreme physical scale.
    p.galaxyBodies = {
        {GalaxyBody3DType::GasCloud,        9.0f,  0.10f, -0.30f},  // Cooling condensation A
        {GalaxyBody3DType::GasCloud,       14.0f,  0.20f, -0.08f},  // Cooling condensation B
        {GalaxyBody3DType::GasCloud,       22.0f,  0.35f,  0.12f},  // Ionised gas filament
        {GalaxyBody3DType::Star,           18.0f,  0.60f,  0.28f},  // Nuclear star
        {GalaxyBody3DType::StellarCluster, 55.0f,  0.15f,  0.40f},  // Inner BCG cluster
        {GalaxyBody3DType::StellarCluster, 90.0f,  0.08f, -0.20f},  // Outer BCG cluster
        {GalaxyBody3DType::DwarfGalaxy,   130.0f,  0.25f,  0.55f}   // Infalling satellite
    };
    return p;
}

/*--------- M87* (Virgo A), first black hole ever imaged ---------*/
// The Event Horizon Telescope target: a 6.5e9 Msun SMBH in the giant
// elliptical M87, famous for its bright asymmetric photon ring and a
// ~5000 ly relativistic jet (HST-1 knot). LINER-class nucleus (modest BLR)
// but an extraordinarily powerful, well-collimated jet.
inline BlackHoleProfile m87() {
    BlackHoleProfile p;
    p.name        = "M87*";
    p.description = "First imaged BH (~6.5e9 Msun), EHT ring + relativistic jet";
    p.massSolar   = 6.5e9;

    cfg::SimConfig c;
    // High spin powers the Blandford-Znajek jet; EHT favours a*~0.9.
    c.blackHole.spinParameter = 0.9f;
    c.blackHole.radius        = 1.0f;

    // Bright, fairly compact ring with a hot inner edge (close ISCO at a*=0.9).
    c.disk.innerRadius   = 2.0f;
    c.disk.outerRadius   = 18.0f;
    c.disk.halfThickness = 0.05f;   // Geometrically thick RIAF-like flow

    // The signature M87 jet: long, well-collimated, brilliant synchrotron blue.
    c.jet.radius = 0.40f;
    c.jet.length = 60.0f;           // ~5000 ly one-sided jet
    c.jet.color  = glm::vec3(0.35f, 0.80f, 1.50f);

    // LINER nucleus: weak broad-line region.
    c.blr.innerRadius = 12.0f;
    c.blr.outerRadius = 26.0f;
    c.blr.thickness   = 5.0f;
    c.blr.strength    = 0.15f;      // M87: low-luminosity AGN, faint BLR

    c.orbital.semiMajor  = 60.0f;
    c.orbital.bodyRadius = 1.5f;

    c.camera.initialPos = glm::vec3(0.0f, 7.0f, 26.0f);

    // Moderate bloom; the ring is bright but the nucleus is sub-Eddington.
    c.bloom.threshold = 1.1f;
    c.bloom.intensity = 0.55f;
    c.bloom.exposure  = 1.05f;

    // EHT image palette: warm orange-gold asymmetric ring.
    c.disk.peakTemp             = 12000.0f;
    c.disk.displayTempInner     = 4200.0f;   // Warm gold inner ring
    c.disk.displayTempOuter     = 2000.0f;   // Deep orange outer
    c.disk.saturationBoostInner = 2.6f;
    c.disk.saturationBoostOuter = 1.7f;

    p.config = c;
    p.defaultJets       = true;    // M87's defining feature
    p.defaultBLR        = false;   // LINER: no prominent BLR
    p.defaultOrbBody    = true;    // nuclear gas disk, stars, globular clusters
    p.defaultDoppler    = true;    // strong jet Doppler asymmetry
    p.defaultHostGalaxy = true;    // embedded in the giant elliptical Virgo A
    p.defaultLAB        = false;
    p.defaultCGM        = false;

    // M87 nuclear environment: the famous ionised gas disk (Ford/Harms 1994),
    // a few nuclear stars, and a couple of the galaxy's ~12,000 globular
    // clusters threading the inner halo. These are the bodies a merging
    // secondary will tug and tidally shred.
    p.galaxyBodies = {
        {GalaxyBody3DType::GasCloud,        8.0f,  0.12f, -0.22f, "Ionised Gas Disk"},
        {GalaxyBody3DType::GasCloud,       13.0f,  0.22f,  0.06f, "Nuclear Filament"},
        {GalaxyBody3DType::Star,           17.0f,  0.55f,  0.30f, "Nuclear Star"},
        {GalaxyBody3DType::Star,           24.0f,  0.40f, -0.18f, "Inner S-star"},
        {GalaxyBody3DType::StellarCluster, 70.0f,  0.12f,  0.42f, "Globular Cluster"},
        {GalaxyBody3DType::DwarfGalaxy,   120.0f,  0.22f, -0.30f, "Infalling Satellite"}
    };
    return p;
}

/*--------- All profiles in order ---------*/
inline std::array<BlackHoleProfile, 14> allProfiles() {
    return { ton618(), sgrAstar(), m87(), qso3c273(), j0529(),
             gaiaBH1(), gaiaBH2(), gaiaBH3(),
             v404Cyg(), a062000(), groJ165540(),
             ngc1277(), oj287(), phoenixA() };
}

constexpr int NUM_PROFILES = 14;

// Look up a full profile by its display name (matches BlackHoleProfile::name).
// Returns nullptr if no profile carries that name. Used to give an inspiralling
// merger secondary the accretion-disk / spin / jet identity of its standalone
// counterpart. The returned pointer refers to a function-local static, valid
// for the lifetime of the program.
inline const BlackHoleProfile* findProfileByName(const char* name) {
    if (!name) return nullptr;
    static const std::array<BlackHoleProfile, 14> kProfiles = allProfiles();
    for (const auto& p : kProfiles) {
        if (p.name == name) return &p;
    }
    return nullptr;
}

} // namespace profiles

/*--------- Legacy compatibility for BlackHole2D ---------*/
// The 2D simulator uses a flat name/mass/description struct.
struct BlackHolePreset {
    std::string name;
    double massSolar;
    std::string description;
};

const inline BlackHolePreset BH_PRESETS[] = {
    {"TON 618",           6.6e10, "Most massive known quasar BH (~6.6e10 Msun)"},
    {"Sgr A*",            4.3e6,  "Milky Way center (~4.3e6 Msun)"},
    {"3C 273",            8.86e8, "First quasar identified (~8.9e8 Msun)"},
    {"J0529-4351",        1.7e10, "Most luminous quasar (~1.7e10 Msun)"},
    {"M87",               6.5e9,  "M87 galaxy BH (~6.5e9 Msun, EHT image)"},
    {"Cygnus X-1",        21.2,   "Famous stellar-mass HMXB black hole (~21 Msun)"},
    {"LIGO GW150914",     62.0,   "First gravitational wave merger remnant (~62 Msun)"},
    {"Intermediate-mass", 1e4,    "Hypothetical intermediate-mass BH (~1e4 Msun)"},
    {"Primordial",        1e-5,   "Tiny primordial black hole (~1e-5 Msun)"},
    {"Gaia BH1",          9.62,   "Nearest dormant BH, G-dwarf companion (~9.6 Msun)"},
    {"Gaia BH2",          8.94,   "Dormant BH with red giant companion (~8.9 Msun)"},
    {"Gaia BH3",          32.7,   "Most massive nearby dormant BH (~33 Msun)"},
    {"V404 Cygni",        9.0,    "Triple system: ~9 Msun BH + K-giant donor + distant tertiary (V404 Cyg C)"},
    {"A0620-00",          6.61,   "First confirmed stellar BH, quiescent K-dwarf (~6.6 Msun)"},
    {"GRO J1655-40",      6.3,    "Microquasar with relativistic jets (~6.3 Msun)"},
    {"NGC 1277",          1.7e10, "Overmassive SMBH in compact elliptical (~1.7e10 Msun)"},
    {"OJ 287",            1.8e10, "SMBH binary, optical outburst source (~1.8e10 Msun)"},
    {"Phoenix A",         1.0e11, "Largest known BH, Phoenix Cluster BCG (~1.0e11 Msun)"}
};
constexpr int NUM_BH_PRESETS = sizeof(BH_PRESETS) / sizeof(BH_PRESETS[0]);
