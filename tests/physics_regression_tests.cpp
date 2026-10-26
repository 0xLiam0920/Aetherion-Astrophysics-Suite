#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include "src/2D/2D-physics/geodesic.hpp"
#include "src/2D/2D-physics/schwarzschild.hpp"
#include "src/2D/2D-physics/units.hpp"
#include "src/2D/2D-simulation/photon.hpp"
#include "src/3D/bh3d_blackbody.hpp"
#include "src/3D/bh3d_shaderinject.hpp"

// Directory holding the .frag / .glsl sources. CMake passes the absolute path;
// a bare `c++ -I.` CI compile run from the repo root uses the relative default.
#ifndef AETHERION_SHADER_DIR
#  define AETHERION_SHADER_DIR "src/3D"
#endif

// The JSON-catalog round-trip test needs glm (config.hpp) and nlohmann/json.
// The lightweight CI job (a bare `c++ -I. …` invocation) has neither on its
// include path, so we compile this block in only when both are available —
// which is exactly the CMake build, where the catalog dir is passed via
// AETHERION_CATALOG_DIR. CI keeps running the five GUI-free physics tests.
#if defined(__has_include)
#  if __has_include(<glm/glm.hpp>) && __has_include(<nlohmann/json.hpp>)
#    define AETHERION_HAVE_CATALOG 1
#    include "src/common/catalog.hpp"
#    include "src/3D/bh3d_presets.hpp"
#    include "src/3D/bh3d_catalog_adapter.hpp"
#    include "src/3D/orbital_body.hpp"
#    include <vector>
#  endif
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace {

bool testSolarGrazingDeflection() {
    // Reference: weak-field GR prediction for a grazing solar ray is ~1.75 arcsec.
    // In geometric units (G=c=1): alpha ~= 4M / b.
    Schwarzschild bh;
    bh.M = units::solarMassToGeomMeters(1.0);

    constexpr double solarRadiusMeters = 6.9634e8;
    constexpr double expectedArcsec = 1.75;
    constexpr double expectedDeg = expectedArcsec / 3600.0;

    Photon photon;
    photon.impactParameter = solarRadiusMeters;
    photon.computePath(bh, /*rMax=*/1.0e12, /*dphi=*/0.0015);

    if (photon.captured) {
        std::cerr << "[FAIL] solar deflection: photon was unexpectedly captured\n";
        return false;
    }

    const double measuredDeg = std::abs(photon.deflectionDeg);
    const double absErrDeg = std::abs(measuredDeg - expectedDeg);

    // Conservative beta tolerance: 0.12 arcsec allows discretization + finite-rMax effects.
    constexpr double toleranceDeg = 0.12 / 3600.0;
    if (absErrDeg > toleranceDeg) {
        std::cerr << "[FAIL] solar deflection: expected ~" << expectedDeg
                  << " deg (" << expectedArcsec << " arcsec), got " << measuredDeg
                  << " deg, abs err " << absErrDeg << " deg\n";
        return false;
    }

    std::cout << "[PASS] solar deflection: " << measuredDeg << " deg (target "
              << expectedDeg << " deg)\n";
    return true;
}

bool testTimelikeConservationLongRun() {
    Schwarzschild bh;
    bh.M = 1.0;

    // Moderately eccentric bound orbit, safely outside strong-instability region.
    const double rPeri = 12.0 * bh.M;
    const double rApo = 28.0 * bh.M;
    auto params = bh.boundOrbitEL(rPeri, rApo);

    TimelikeState state{rPeri, 0.0, 0.0};
    const double E0 = bh.computeEnergy(state.r, state.vr, params.L);
    const double L0 = params.L;

    constexpr double dtau = 0.02;
    constexpr int steps = 60000;

    double maxRelEnergyDrift = 0.0;
    double maxRelLDrift = 0.0;

    TimelikeState prev = state;
    for (int i = 0; i < steps; ++i) {
        state = stepTimelikeGeodesic(bh, state, params.L, dtau);

        if (!std::isfinite(state.r) || !std::isfinite(state.phi) || !std::isfinite(state.vr)) {
            std::cerr << "[FAIL] conservation: non-finite state at step " << i << "\n";
            return false;
        }

        const double E = bh.computeEnergy(state.r, state.vr, params.L);
        const double relEDrift = std::abs(E - E0) / std::max(1e-15, std::abs(E0));
        maxRelEnergyDrift = std::max(maxRelEnergyDrift, relEDrift);

        // Infer L numerically from dphi/dtau to check angular momentum consistency.
        const double dphi = state.phi - prev.phi;
        const double rAvg = 0.5 * (state.r + prev.r);
        const double Linferred = (dphi / dtau) * rAvg * rAvg;
        const double relLDrift = std::abs(Linferred - L0) / std::max(1e-15, std::abs(L0));
        maxRelLDrift = std::max(maxRelLDrift, relLDrift);

        prev = state;
    }

    // Beta gate: generous enough for RK4 long runs, strict enough to catch regressions.
    constexpr double maxEnergyDriftAllowed = 5e-3;
    constexpr double maxLDriftAllowed = 2e-2;

    if (maxRelEnergyDrift > maxEnergyDriftAllowed) {
        std::cerr << "[FAIL] conservation: max relative energy drift " << maxRelEnergyDrift
                  << " exceeds " << maxEnergyDriftAllowed << "\n";
        return false;
    }
    if (maxRelLDrift > maxLDriftAllowed) {
        std::cerr << "[FAIL] conservation: max relative angular momentum drift " << maxRelLDrift
                  << " exceeds " << maxLDriftAllowed << "\n";
        return false;
    }

    std::cout << "[PASS] conservation: max dE/E=" << maxRelEnergyDrift
              << ", max dL/L=" << maxRelLDrift << "\n";
    return true;
}

bool testPhotonCaptureThreshold() {
    // Photons with impact parameter b < b_crit = 3*sqrt(3)*M are always
    // captured; those with b > b_crit always escape. Check both sides of the
    // critical value to make sure the capture logic matches the analytic
    // photon-sphere result.
    Schwarzschild bh;
    bh.M = 1.0;
    const double bCrit = bh.criticalImpact();

    Photon below;
    below.impactParameter = 0.98 * bCrit; // just inside -> must be captured
    below.computePath(bh, /*rMax=*/1.0e6, /*dphi=*/0.0015);
    if (!below.captured) {
        std::cerr << "[FAIL] capture: b=" << below.impactParameter
                  << " (< b_crit=" << bCrit << ") escaped but should be captured\n";
        return false;
    }

    Photon above;
    above.impactParameter = 1.02 * bCrit; // just outside -> must escape
    above.computePath(bh, /*rMax=*/1.0e6, /*dphi=*/0.0015);
    if (above.captured) {
        std::cerr << "[FAIL] capture: b=" << above.impactParameter
                  << " (> b_crit=" << bCrit << ") captured but should escape\n";
        return false;
    }

    std::cout << "[PASS] capture threshold: b_crit=" << bCrit
              << " (captured below, escaped above)\n";
    return true;
}

bool testPeriapsisTurningPoint() {
    // findPeriapsis(b) must return the radius r where the null turning-point
    // condition r^2 / f(r) = b^2 holds. Verify the residual is ~0 across a
    // range of impact parameters, and that sub-critical rays report no
    // turning point (-1).
    Schwarzschild bh;
    bh.M = 1.0;

    bool ok = true;
    double worstResidual = 0.0;
    for (double b : {6.0, 8.0, 12.0, 25.0, 100.0}) {
        const double r = bh.findPeriapsis(b);
        if (r <= 0.0) {
            std::cerr << "[FAIL] periapsis: b=" << b
                      << " (> b_crit) returned no turning point\n";
            ok = false;
            continue;
        }
        // Turning-point residual normalised by b^2.
        const double residual = std::abs((r * r) / bh.f(r) - b * b) / (b * b);
        worstResidual = std::max(worstResidual, residual);
        if (residual > 1e-6) {
            std::cerr << "[FAIL] periapsis: b=" << b << " r=" << r
                      << " turning-point residual " << residual << " too large\n";
            ok = false;
        }
    }

    // Sub-critical impact parameter has no real turning point.
    if (bh.findPeriapsis(0.9 * bh.criticalImpact()) > 0.0) {
        std::cerr << "[FAIL] periapsis: sub-critical b reported a turning point\n";
        ok = false;
    }

    if (ok) {
        std::cout << "[PASS] periapsis turning point: worst residual "
                  << worstResidual << "\n";
    }
    return ok;
}

bool testStrongDeflectionGrowth() {
    // The deflection angle grows without bound (logarithmically) as b -> b_crit
    // from above, far exceeding the weak-field value at large b. Check that a
    // near-critical ray deflects substantially more than one full turn, and
    // much more than a wide ray.
    Schwarzschild bh;
    bh.M = 1.0;
    const double bCrit = bh.criticalImpact();

    Photon wide;
    wide.impactParameter = 50.0 * bh.M; // weak field: alpha ~ 4M/b ~ 0.08 rad
    wide.computePath(bh, /*rMax=*/1.0e7, /*dphi=*/0.0015);

    Photon strong;
    strong.impactParameter = 1.01 * bCrit; // strong field: many radians
    strong.computePath(bh, /*rMax=*/1.0e7, /*dphi=*/0.0005);

    if (wide.captured || strong.captured) {
        std::cerr << "[FAIL] strong deflection: a test ray was captured\n";
        return false;
    }

    const double wideDefl = std::abs(wide.deflectionAngle);
    const double strongDefl = std::abs(strong.deflectionAngle);

    // Weak-field ray should be a small fraction of a radian.
    if (wideDefl > 0.2) {
        std::cerr << "[FAIL] strong deflection: wide ray deflection " << wideDefl
                  << " rad unexpectedly large\n";
        return false;
    }
    // Near-critical ray should exceed a full half-turn (pi) and dwarf the wide ray.
    if (strongDefl < M_PI) {
        std::cerr << "[FAIL] strong deflection: near-critical ray deflection "
                  << strongDefl << " rad did not reach pi\n";
        return false;
    }
    if (strongDefl < 10.0 * wideDefl) {
        std::cerr << "[FAIL] strong deflection: near-critical/wide ratio "
                  << (strongDefl / wideDefl) << " too small\n";
        return false;
    }

    std::cout << "[PASS] strong deflection: wide=" << wideDefl
              << " rad, near-critical=" << strongDefl << " rad\n";
    return true;
}

// ────────────────────────────────────────────────────────────────────────────
// ACCRETION-DISK COLOUR PHYSICS
//
// These validate the physical colour pipeline shared with the GPU shader
// through src/3D/bh3d_blackbody.hpp (the same code that builds the blackbody
// LUT the shader samples). Since that header is shared, these tests guarantee
// the disk always shows a physically-correct temperature-to-colour gradient
// instead of quietly regressing back to a flat white sheet.
// ────────────────────────────────────────────────────────────────────────────

// Works out the "blueness" of a colour. Positive means it leans blue,
// negative means it leans red.
double blueness(const std::array<double, 3>& rgb) { return rgb[2] - rgb[0]; }

bool testBlackbodyMonotonicBlueness() {
    // Cooler blackbodies have to be redder and hotter ones bluer, with no
    // exceptions along the way, so this should increase strictly monotonically.
    const double temps[] = {1500, 2500, 3500, 5000, 6600, 9000, 15000, 25000, 40000};
    double prev = -2.0;
    for (double T : temps) {
        const double bl = blueness(bh3d::physics::blackbodyRGB(T));
        if (bl <= prev) {
            std::cerr << "[FAIL] blackbody blueness not monotincreasing at T=" << T
                      << " K (blueness " << bl << " <= previous " << prev << ")\n";
            return false;
        }
        prev = bl;
    }
    std::cout << "[PASS] blackbody blueness increases monotonically with temperature\n";
    return true;
}

bool testBlackbodyReferencePoints() {
    using bh3d::physics::blackbodyRGB;

    // Every entry needs to be a valid, chromaticity-normalised sRGB triple.
    for (double T = 1000.0; T <= 40000.0; T += 250.0) {
        const auto c = blackbodyRGB(T);
        const double mx = std::max({c[0], c[1], c[2]});
        for (double ch : c) {
            if (ch < -1e-9 || ch > 1.0 + 1e-9) {
                std::cerr << "[FAIL] blackbody channel out of [0,1] at T=" << T
                          << " K: " << ch << "\n";
                return false;
            }
        }
        if (mx < 0.98) {  // normalised so max=1, gamma(1)=1
            std::cerr << "[FAIL] blackbody not normalised at T=" << T
                      << " K (max channel " << mx << ")\n";
            return false;
        }
    }

    // A cool star around 3000 K should read warm and orange, so red ends up
    // the strongest channel and blue stays weak.
    const auto warm = blackbodyRGB(3000.0);
    if (!(warm[0] > warm[1] && warm[1] > warm[2] && warm[2] < 0.6)) {
        std::cerr << "[FAIL] 3000 K should be orange (R>G>B, weak blue); got "
                  << warm[0] << "," << warm[1] << "," << warm[2] << "\n";
        return false;
    }

    // Something Sun-like around 6600 K should look close to white, so red and
    // blue ought to be fairly balanced.
    const auto white = blackbodyRGB(6600.0);
    if (std::abs(white[0] - white[2]) > 0.12) {
        std::cerr << "[FAIL] 6600 K should be near-white (|R-B|<=0.12); got R="
                  << white[0] << " B=" << white[2] << "\n";
        return false;
    }

    // A hot O-star around 25000 K should read blue-white, so blue needs to be
    // the dominant channel here.
    const auto hot = blackbodyRGB(25000.0);
    if (!(hot[2] >= hot[0] && hot[2] >= hot[1])) {
        std::cerr << "[FAIL] 25000 K should be blue-dominant; got "
                  << hot[0] << "," << hot[1] << "," << hot[2] << "\n";
        return false;
    }

    std::cout << "[PASS] blackbody reference points (orange / white / blue) correct\n";
    return true;
}

bool testNovikovThorneProfile() {
    using bh3d::physics::novikovThorneTemperature;
    using bh3d::physics::computeISCO;

    // Check the relativistic profile across three spins spanning the range.
    for (double a : {0.0, 0.5, 0.9}) {
        const double rISCO = computeISCO(a);

        // Zero-torque inner boundary: temperature vanishes at the ISCO.
        if (novikovThorneTemperature(rISCO, a) > 1e-6) {
            std::cerr << "[FAIL] NT profile should be ~0 at the ISCO (a=" << a
                      << "); got " << novikovThorneTemperature(rISCO, a) << "\n";
            return false;
        }

        // Scan outward: locate the peak, confirm it normalises to ~1, and that
        // the profile stays bounded within [0,1] the whole way.
        double maxT = 0.0, rAtMax = rISCO;
        for (double r = rISCO; r <= 40.0; r += 0.05) {
            const double t = novikovThorneTemperature(r, a);
            if (t < -1e-9 || t > 1.0 + 1e-6) {
                std::cerr << "[FAIL] NT profile out of [0,1] at r=" << r
                          << " (a=" << a << "): " << t << "\n";
                return false;
            }
            if (t > maxT) { maxT = t; rAtMax = r; }
        }
        if (std::abs(maxT - 1.0) > 5e-3) {
            std::cerr << "[FAIL] NT profile peak should normalise to 1.0 (a=" << a
                      << "); got " << maxT << " at r=" << rAtMax << "\n";
            return false;
        }

        // Beyond the peak it must decrease monotonically outward.
        double prev = 2.0;
        for (double r = rAtMax; r <= 40.0; r += 0.5) {
            const double t = novikovThorneTemperature(r, a);
            if (t > prev + 1e-6) {
                std::cerr << "[FAIL] NT profile not decreasing outward at r=" << r
                          << " (a=" << a << "): " << t << " > " << prev << "\n";
                return false;
            }
            prev = t;
        }
    }

    std::cout << "[PASS] Novikov-Thorne profile (Page-Thorne): zero at ISCO, "
                 "peak=1, monotone falloff across spins\n";
    return true;
}

bool testComputeISCO() {
    using bh3d::physics::computeISCO;

    // Schwarzschild: exactly 3 Rs (r_ms = 6M).
    if (std::abs(computeISCO(0.0) - 3.0) > 1e-9) {
        std::cerr << "[FAIL] Schwarzschild ISCO should be 3 Rs; got "
                  << computeISCO(0.0) << "\n";
        return false;
    }
    // Concrete Kerr checkpoint: a*=0.9 prograde ISCO ≈ 1.161 Rs (r_ms ≈ 2.321 M).
    if (std::abs(computeISCO(0.9) - 1.1605) > 2e-3) {
        std::cerr << "[FAIL] Kerr ISCO at a=0.9 should be ~1.161 Rs; got "
                  << computeISCO(0.9) << "\n";
        return false;
    }
    // Prograde ISCO shrinks monotonically with spin toward the 0.5 Rs extremal
    // limit, and never dips below it.
    double prev = 4.0;
    for (double a = 0.0; a <= 0.999; a += 0.02) {
        const double r = computeISCO(a);
        if (r > prev + 1e-9) {
            std::cerr << "[FAIL] ISCO should decrease with spin; a=" << a
                      << " gave " << r << " > " << prev << "\n";
            return false;
        }
        if (r < 0.5 - 1e-6) {
            std::cerr << "[FAIL] ISCO fell below extremal 0.5 Rs at a=" << a
                      << ": " << r << "\n";
            return false;
        }
        prev = r;
    }
    std::cout << "[PASS] ISCO(a*): 3 Rs at a=0, 1.16 Rs at a=0.9, monotone to 0.5 Rs\n";
    return true;
}

bool testDiskColourGradient() {
    using bh3d::physics::blackbodyRGB;
    using bh3d::physics::computeISCO;
    using bh3d::physics::diskColorTemperature;

     // Since every black-hole preset's real disk config are mirrored from
    // src/3D/bh3d_presets.hpp ({name, displayTempInner [K], spin a*,
    // innerRadius [Rs], outerRadius [Rs]}), we have to keep in sync with the presets. Given that
    // bh3d_presets.hpp itself isn't included here because it pulls in glm, 
    // it's largely a GUI-free CI job intentionally doesn't depend on.
    // outerRadius == 0 marks presets that disable the disk entirely, 
    // for now as the Gaia series of BlackHoles given they have zero accretion disks since none are active.
    struct Preset {
        const char* name;
        double inner;   // displayTempInner [K]
        double spin;
        double rIn;
        double rOut;
    };
    const Preset presets[] = {
        {"TON 618",       5500.0,  0.80,  2.0,   20.0},
        {"Sgr A*",        3200.0,  0.5,   3.0,   10.0},
        {"M87*",          4200.0,  0.9,   2.0,   18.0},
        {"3C 273",        6500.0,  0.90,  1.8,   22.0},
        {"J0529-4351",    7000.0,  0.95,  1.5,   30.0},
        {"Gaia BH1",      1000.0,  0.3,   3.0,   0.0},   // disk disabled
        {"Gaia BH2",      1000.0,  0.2,   3.0,   0.0},   // disk disabled
        {"Gaia BH3",      1000.0,  0.45,  2.8,   0.0},   // disk disabled
        {"V404 Cygni",    7500.0,  0.5,   2.5,   12.0},
        {"A0620-00",      3000.0,  0.12,  3.2,   8.0},
        {"GRO J1655-40",  9500.0,  0.7,   1.9,   10.0},
        {"NGC 1277",      6800.0,  0.65,  1.9,   20.0},
        {"OJ 287",       11000.0,  0.82,  1.75,  26.0},
        {"Phoenix A",     4800.0,  0.35,  5.0,   30.0},
    };

    bool ok = true;
    int diskCount = 0, disabledCount = 0;
    for (const auto& p : presets) {
        if (p.rOut <= p.rIn) {
            // The disk is intentionally disabled here (detached, dormant Gaia
            // BHs fall into this bucket), so there's no colour to validate.
            // We still count it, just to confirm it's genuinely inert.
            ++disabledCount;
            continue;
        }
        ++diskCount;

        const double rISCO = computeISCO(p.spin);
        const double rPeak = std::clamp((49.0 / 36.0) * rISCO, p.rIn, p.rOut);

        const double tInner = diskColorTemperature(p.inner, rPeak, p.spin);
        const double tOuter = diskColorTemperature(p.inner, p.rOut, p.spin);

        if (p.inner < 1000.0 || p.inner > 40000.0) {
            std::cerr << "[FAIL] " << p.name << " displayTempInner " << p.inner
                      << " K outside LUT range [1000,40000]\n";
            ok = false;
            continue;
        }
        // The disk has to actually cool outward, since that's the whole point
        // of having a gradient in the first place.
        if (!(tInner > tOuter + 200.0)) {
            std::cerr << "[FAIL] " << p.name << " disk not cooling outward: inner "
                      << tInner << " K vs outer " << tOuter << " K (rISCO="
                      << rISCO << ", rPeak=" << rPeak << ")\n";
            ok = false;
            continue;
        }
        // And that temperature drop needs to actually show up as colour too,
        // with the inner disk reading bluer than the outer disk.
        const double blInner = blueness(blackbodyRGB(tInner));
        const double blOuter = blueness(blackbodyRGB(tOuter));
        if (!(blInner > blOuter + 0.05)) {
            std::cerr << "[FAIL] " << p.name << " colour gradient too flat: inner "
                      << "blueness " << blInner << " (" << tInner << " K) vs outer "
                      << blOuter << " (" << tOuter << " K), so this disk risks "
                      << "rendering as a flat colour\n";
            ok = false;
        }
    }

    if (ok) {
        std::cout << "[PASS] disk colour gradient: " << diskCount
                  << " presets cool outward with a visible colour gradient ("
                  << disabledCount << " disk-disabled presets skipped)\n";
    }
    return ok;
}

// ===========================================================================
// Step 10: CPU reference geodesic integrator (GL-free, matches bh3d_physics.glsl)
// ---------------------------------------------------------------------------
// Units: Rs = 1, c = 1, GM = 0.5. A photon's path obeys
//   d2x/dlambda2 = -1.5 * h2 * x / r^5,   h = |x cross v|  (conserved).
// This mirrors geoAccel()/geoAdvance() in the shared shader so the shadow
// edge, weak deflection and low-quality path are all checked in CI.
// ===========================================================================
struct GVec3 { double x, y, z; };
inline GVec3 operator+(GVec3 a, GVec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline GVec3 operator-(GVec3 a, GVec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline GVec3 operator*(GVec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double gdot(GVec3 a, GVec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline GVec3 gcross(GVec3 a, GVec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double glength(GVec3 a) { return std::sqrt(gdot(a, a)); }
inline GVec3 gnormalize(GVec3 a) { double l = glength(a); return l > 0.0 ? a * (1.0 / l) : a; }

inline GVec3 geoAccel(GVec3 x, double h2) {
    const double r2 = gdot(x, x);
    return x * (-1.5 * h2 / (r2 * r2 * std::sqrt(r2)));
}

// RK4 advance (photoreal quality): 4 acceleration evals per step.
inline void geoAdvanceRK4(GVec3& x, GVec3& v, double h2, double dt) {
    const GVec3 k1x = v,                    k1v = geoAccel(x, h2);
    const GVec3 k2x = v + k1v * (0.5 * dt), k2v = geoAccel(x + k1x * (0.5 * dt), h2);
    const GVec3 k3x = v + k2v * (0.5 * dt), k3v = geoAccel(x + k2x * (0.5 * dt), h2);
    const GVec3 k4x = v + k3v * dt,         k4v = geoAccel(x + k3x * dt, h2);
    x = x + (k1x + k2x * 2.0 + k3x * 2.0 + k4x) * (dt / 6.0);
    v = v + (k1v + k2v * 2.0 + k3v * 2.0 + k4v) * (dt / 6.0);
}

// Velocity Verlet advance (low quality): 1 acceleration eval per step,
// carrying `a` between calls exactly as the QUALITY==0 shader path does.
inline void geoAdvanceVerlet(GVec3& x, GVec3& v, GVec3& a, double h2, double dt) {
    v = v + a * (0.5 * dt);
    x = x + v * dt;
    a = geoAccel(x, h2);
    v = v + a * (0.5 * dt);
}

inline double gsmoothstep(double e0, double e1, double x) {
    const double t = std::clamp((x - e0) / (e1 - e0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}
// Low-quality step length proportional to r (shader QUALITY==0: FINE/COARSE).
inline double geoStepLow(double r) {
    const double FINE = 0.05, COARSE = 0.14;
    return (FINE + (COARSE - FINE) * gsmoothstep(1.5, 6.0, r)) * r;
}

// Ray from (-D, b, 0) heading +X. Returns true if captured (r < 1);
// on escape, *vOut is the final (normalised) direction.
inline bool traceRK4(double b, double D, GVec3* vOut = nullptr) {
    GVec3 x{-D, b, 0.0}, v{1.0, 0.0, 0.0};
    const GVec3 h = gcross(x, v);
    const double h2 = gdot(h, h);
    for (int i = 0; i < 400000; ++i) {
        const double r = glength(x);
        if (r < 1.0) return true;
        if (r > 2.0 * D && gdot(x, v) > 0.0) { if (vOut) *vOut = gnormalize(v); return false; }
        geoAdvanceRK4(x, v, h2, std::min(0.01 * r, 5.0));
    }
    return false;
}

// Same trace, but using the low-quality Verlet integrator + geoStepLow.
inline bool traceVerletLow(double b, double D, GVec3* vOut = nullptr) {
    GVec3 x{-D, b, 0.0}, v{1.0, 0.0, 0.0};
    const GVec3 h = gcross(x, v);
    const double h2 = gdot(h, h);
    GVec3 a = geoAccel(x, h2);
    for (int i = 0; i < 400000; ++i) {
        const double r = glength(x);
        if (r < 1.0) return true;
        if (r > 2.0 * D && gdot(x, v) > 0.0) { if (vOut) *vOut = gnormalize(v); return false; }
        geoAdvanceVerlet(x, v, a, h2, geoStepLow(r));
    }
    return false;
}

// Bisect the capture/escape boundary (the shadow edge) in b.
template <typename TraceFn>
double shadowEdge(TraceFn trace, double D, int iters = 40) {
    double lo = 2.0, hi = 3.5;   // captured below b_crit, escapes above
    for (int i = 0; i < iters; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (trace(mid, D)) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

bool testGeodesicShadowEdge() {
    // The RK4 integrator's capture boundary must land on the analytic
    // Schwarzschild critical impact parameter b_crit = 3*sqrt(3)/2 Rs.
    const double bCrit = 3.0 * std::sqrt(3.0) / 2.0;   // 2.598076...
    const double b = shadowEdge([](double bb, double D) { return traceRK4(bb, D); }, 1000.0);
    if (std::abs(b - bCrit) > 0.005) {
        std::cerr << "[FAIL] geodesic shadow edge: expected b_crit=" << bCrit
                  << ", got " << b << " (err " << std::abs(b - bCrit) << ")\n";
        return false;
    }
    std::cout << "[PASS] geodesic shadow edge: b_crit=" << b
              << " (analytic " << bCrit << ")\n";
    return true;
}

bool testGeodesicWeakDeflection() {
    // Far from the hole the bending is the weak-field value alpha = 2 Rs / b
    // (= 4GM/(c^2 b) with GM = 0.5). b = 200 => alpha = 0.0100 rad.
    const double b = 200.0;
    GVec3 vOut{1.0, 0.0, 0.0};
    if (traceRK4(b, 1.0e4, &vOut)) {
        std::cerr << "[FAIL] geodesic weak deflection: ray at b=200 was captured\n";
        return false;
    }
    const double alpha = std::acos(std::clamp(vOut.x, -1.0, 1.0));
    const double expected = 2.0 / b;   // 0.0100
    if (std::abs(alpha - expected) > 0.03 * expected) {
        std::cerr << "[FAIL] geodesic weak deflection: expected ~" << expected
                  << " rad, got " << alpha << " rad\n";
        return false;
    }
    std::cout << "[PASS] geodesic weak deflection: alpha=" << alpha
              << " rad vs 2/b=" << expected << "\n";
    return true;
}

bool testGeodesicMatches2D() {
    // The 3D Cartesian RK4 integrator must reproduce the independent 2D
    // angular integrator (Photon::computePath) for the same M = 0.5 hole.
    Schwarzschild bh;
    bh.M = 0.5;   // Rs = 2M = 1, matching the 3D units

    bool ok = true;
    for (double b : {5.0, 8.0, 12.0, 20.0}) {
        Photon ph;
        ph.impactParameter = b;
        ph.computePath(bh, /*rMax=*/1.0e5, /*dphi=*/0.0005);
        if (ph.captured) {
            std::cerr << "[FAIL] geodesic-2D match: 2D ray captured at b=" << b << "\n";
            ok = false;
            continue;
        }
        const double d2d = std::abs(ph.deflectionDeg);

        GVec3 vOut{1.0, 0.0, 0.0};
        if (traceRK4(b, 1.0e4, &vOut)) {
            std::cerr << "[FAIL] geodesic-2D match: 3D ray captured at b=" << b << "\n";
            ok = false;
            continue;
        }
        const double d3d = std::acos(std::clamp(vOut.x, -1.0, 1.0)) * 180.0 / M_PI;

        // Two different integrators (angular vs Cartesian RK4); allow a small
        // relative gap plus an absolute floor for the tiny-angle cases.
        const double tol = std::max(0.03 * d2d, 0.05);
        if (std::abs(d3d - d2d) > tol) {
            std::cerr << "[FAIL] geodesic-2D match at b=" << b << ": 2D=" << d2d
                      << " deg vs 3D=" << d3d << " deg (tol " << tol << ")\n";
            ok = false;
        }
    }
    if (ok)
        std::cout << "[PASS] geodesic 3D RK4 matches 2D Photon::computePath across b\n";
    return ok;
}

bool testLowQualityShadow() {
    // The low-quality Verlet path (i3 GPUs) must still resolve the shadow edge
    // to within 1% of the analytic b_crit, guarding the cheap integrator.
    const double bCrit = 3.0 * std::sqrt(3.0) / 2.0;
    const double b = shadowEdge([](double bb, double D) { return traceVerletLow(bb, D); }, 1000.0);
    if (std::abs(b - bCrit) > 0.01 * bCrit) {
        std::cerr << "[FAIL] low-quality shadow edge: expected ~" << bCrit
                  << ", got " << b << " (err " << std::abs(b - bCrit) << ")\n";
        return false;
    }
    std::cout << "[PASS] low-quality (Verlet) shadow edge: b=" << b
              << " within 1% of " << bCrit << "\n";
    return true;
}

bool testDiskRedshift() {
    using bh3d::physics::diskRedshift;

    // Face-on at the a=0 ISCO (r=3, Lz=0): pure gravitational + transverse
    // shift g = sqrt(1/2) = 0.70711.
    const double g0 = diskRedshift(3.0, 0.0, 0.0, /*doppler=*/true);
    if (std::abs(g0 - 0.70711) > 1e-3) {
        std::cerr << "[FAIL] disk redshift anchor: r=3,a=0,Lz=0 expected 0.70711, got "
                  << g0 << "\n";
        return false;
    }

    // Doppler asymmetry: the approaching side (Lz<0) is blueshifted (g>1) and
    // the receding side (Lz>0) is redshifted (g<1).
    const double gApproach = diskRedshift(3.0, 0.0, -5.0, true);
    const double gRecede   = diskRedshift(3.0, 0.0,  5.0, true);
    if (!(gApproach > 1.0)) {
        std::cerr << "[FAIL] disk redshift: approaching side should be blue (g>1); got "
                  << gApproach << "\n";
        return false;
    }
    if (!(gRecede < 1.0)) {
        std::cerr << "[FAIL] disk redshift: receding side should be red (g<1); got "
                  << gRecede << "\n";
        return false;
    }
    // g must fall monotonically as Lz increases.
    if (!(gApproach > g0 && g0 > gRecede)) {
        std::cerr << "[FAIL] disk redshift not monotone in Lz: " << gApproach << " / "
                  << g0 << " / " << gRecede << "\n";
        return false;
    }

    // Doppler off: Lz is ignored, so any Lz returns the gravitational anchor.
    const double gOff = diskRedshift(3.0, 0.0, -5.0, /*doppler=*/false);
    if (std::abs(gOff - g0) > 1e-9) {
        std::cerr << "[FAIL] disk redshift: Doppler-off should ignore Lz; got "
                  << gOff << " vs " << g0 << "\n";
        return false;
    }

    std::cout << "[PASS] disk redshift: 0.7071 anchor, blue approach / red recede, "
                 "Doppler-off gravitational only\n";
    return true;
}

bool testBlackbodyLinearRadiance() {
    using bh3d::physics::blackbodyLinear;
    auto luminance = [](const std::array<double, 3>& c) {
        return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
    };

    // Normalisation: luminance Y(6500 K) == 1.
    const double y65 = luminance(blackbodyLinear(6500.0));
    if (std::abs(y65 - 1.0) > 0.02) {
        std::cerr << "[FAIL] blackbodyLinear: Y(6500 K) should be 1.0; got " << y65 << "\n";
        return false;
    }

    // Unlike blackbodyRGB, this keeps absolute brightness: hotter is brighter.
    const double y3k = luminance(blackbodyLinear(3000.0));
    const double y10k = luminance(blackbodyLinear(10000.0));
    if (!(y10k > y65 && y65 > y3k)) {
        std::cerr << "[FAIL] blackbodyLinear: luminance should grow with T; got "
                  << y3k << " / " << y65 << " / " << y10k << "\n";
        return false;
    }

    // Hotter is bluer: the blue/red ratio must increase with temperature.
    auto blueRedRatio = [](const std::array<double, 3>& c) {
        return c[2] / std::max(c[0], 1e-12);
    };
    if (!(blueRedRatio(blackbodyLinear(10000.0)) > blueRedRatio(blackbodyLinear(3000.0)))) {
        std::cerr << "[FAIL] blackbodyLinear: hotter should be bluer\n";
        return false;
    }

    std::cout << "[PASS] blackbody linear radiance: Y(6500 K)=" << y65
              << ", brightness grows with T, hotter is bluer\n";
    return true;
}

// ===========================================================================
// Step 0: shared-GLSL injection
// ===========================================================================
bool testShaderInjection() {
    // A leading comment before #version is legal GLSL and is what the real
    // BlackHole3D.frag looks like, so #version is NOT necessarily line 1.
    const std::string src =
        "// header comment\n"
        "#version 330 core // trailing note\n"
        "uniform float x;\n"
        "void main() {}\n";
    const std::string common = "uniform float shared1;\nfloat helper() { return 1.0; }\n";

    const std::string out = bh3d::injectCommon(src, common, 1);
    const size_t vPos = out.find("#version 330 core");
    const size_t qPos = out.find("#define QUALITY 1\n");
    const size_t cPos = out.find("uniform float shared1;");
    const size_t lPos = out.find("#line 3\n");   // #version is line 2 -> next line is 3
    const size_t bPos = out.find("uniform float x;");

    if (vPos == std::string::npos || qPos == std::string::npos ||
        cPos == std::string::npos || lPos == std::string::npos || bPos == std::string::npos) {
        std::cerr << "[FAIL] shader injection: missing a piece in:\n" << out << "\n";
        return false;
    }
    if (!(vPos < qPos && qPos < cPos && cPos < lPos && lPos < bPos)) {
        std::cerr << "[FAIL] shader injection: pieces are out of order:\n" << out << "\n";
        return false;
    }
    size_t occurrences = 0;
    for (size_t p = out.find("shared1"); p != std::string::npos; p = out.find("shared1", p + 1))
        ++occurrences;
    if (occurrences != 1) {   // the shared block must be spliced in exactly once
        std::cerr << "[FAIL] shader injection: shared code appears " << occurrences
                  << " times, expected 1\n";
        return false;
    }

    // Quality 0 flips the define; the original body must be preserved verbatim.
    const std::string low = bh3d::injectCommon(src, common, 0);
    if (low.find("#define QUALITY 0\n") == std::string::npos ||
        low.find("void main() {}\n") == std::string::npos) {
        std::cerr << "[FAIL] shader injection: quality 0 / body not preserved\n";
        return false;
    }

    // Common without a trailing newline must not glue onto the #line directive.
    const std::string noNl = bh3d::injectCommon(src, "float z;", 1);
    if (noNl.find("float z;\n#line 3\n") == std::string::npos) {
        std::cerr << "[FAIL] shader injection: missing newline after shared block\n";
        return false;
    }

    // Failure modes return empty so callers treat them like a failed read.
    if (!bh3d::injectCommon("void main() {}\n", common, 1).empty()) {
        std::cerr << "[FAIL] shader injection: source without #version should fail\n";
        return false;
    }
    if (!bh3d::injectCommon(src, "", 1).empty()) {
        std::cerr << "[FAIL] shader injection: empty shared code should fail\n";
        return false;
    }

    std::cout << "[PASS] shader injection: order, QUALITY define, #line, failure modes\n";
    return true;
}

namespace {
std::string readTextFile(const std::string& path) {
    std::ifstream f(path);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Collects every top-level uniform name and function name declared in `src`.
std::set<std::string> declaredNames(const std::string& src) {
    static const std::regex uniformRe(R"(^\s*uniform\s+\w+\s+(\w+))");
    static const std::regex funcRe(
        R"(^(?:float|int|bool|void|vec[234]|ivec[234]|mat[234])\s+(\w+)\s*\()");
    std::set<std::string> names;
    std::istringstream in(src);
    std::string line;
    while (std::getline(in, line)) {
        std::smatch m;
        if (std::regex_search(line, m, uniformRe)) names.insert("uniform " + m[1].str());
        else if (std::regex_search(line, m, funcRe)) names.insert("func " + m[1].str());
    }
    return names;
}
}  // namespace

bool testSharedShaderNoDuplicates() {
    // The shared GLSL is injected into both fragment shaders, so any symbol
    // declared in both places is a GLSL redefinition error that only shows up
    // on a machine with a GPU. Catch it here instead.
    const std::string dir = AETHERION_SHADER_DIR;
    const std::string common = readTextFile(dir + "/bh3d_physics.glsl");
    if (common.empty()) {
        std::cerr << "[FAIL] shared shader: could not read " << dir << "/bh3d_physics.glsl\n";
        return false;
    }
    const auto sharedNames = declaredNames(common);
    if (sharedNames.count("func getRayDir") == 0 || sharedNames.count("func hash12") == 0 ||
        sharedNames.count("uniform spinParameter") == 0) {
        std::cerr << "[FAIL] shared shader: expected getRayDir/hash12/spinParameter in "
                     "bh3d_physics.glsl\n";
        return false;
    }

    bool ok = true;
    for (const char* name : {"BlackHole3D.frag", "BlackHole3D_PhotorealDisk.frag"}) {
        const std::string frag = readTextFile(dir + "/" + name);
        if (frag.empty()) {
            std::cerr << "[FAIL] shared shader: could not read " << dir << "/" << name << "\n";
            ok = false;
            continue;
        }
        for (const auto& n : declaredNames(frag)) {
            if (sharedNames.count(n)) {
                std::cerr << "[FAIL] shared shader: '" << n << "' is declared in both "
                          << name << " and bh3d_physics.glsl\n";
                ok = false;
            }
        }
        // And the shader must still be injectable (has a #version line).
        if (bh3d::injectCommon(frag, common, 1).empty()) {
            std::cerr << "[FAIL] shared shader: " << name << " cannot be injected\n";
            ok = false;
        }
    }
    if (ok)
        std::cout << "[PASS] shared shader: no symbol declared in both bh3d_physics.glsl and "
                     "either .frag (" << sharedNames.size() << " shared symbols)\n";
    return ok;
}

// ===========================================================================
// Steps 2 + 3: exact disk physics (CPU twins of the shader code)
// ===========================================================================
bool testNTSetupMatchesReference() {
    using namespace bh3d::physics;
    for (double a : {0.0, 0.38, 0.5, 0.9, 0.998}) {
        const NTSetup s = ntSetup(a);
        if (std::abs(s.isco - computeISCO(a)) > 1e-12 ||
            std::abs(s.x0 - std::sqrt(2.0 * s.isco)) > 1e-12) {
            std::cerr << "[FAIL] ntSetup: isco/x0 wrong at a=" << a << "\n";
            return false;
        }
        double maxT = 0.0;
        for (double r = s.isco * 1.0001; r <= 40.0; r += 0.05) {
            const double t = ntProfileFromSetup(r, s);
            const double ref = novikovThorneTemperature(r, a);
            // Same math, but the reference clamps to [0,1] and normalises with its own scan.
            if (std::abs(std::min(t, 1.0) - ref) > 1e-9) {
                std::cerr << "[FAIL] ntSetup profile mismatch at a=" << a << " r=" << r
                          << ": setup " << t << " vs reference " << ref << "\n";
                return false;
            }
            maxT = std::max(maxT, t);
        }
        if (std::abs(maxT - 1.0) > 5e-3) {
            std::cerr << "[FAIL] ntSetup: peak should normalise to 1; got " << maxT
                      << " at a=" << a << "\n";
            return false;
        }
        if (ntProfileFromSetup(s.isco, s) > 1e-6) {
            std::cerr << "[FAIL] ntSetup: profile must be zero at the ISCO (a=" << a << ")\n";
            return false;
        }
    }
    std::cout << "[PASS] NT setup: per-spin constants reproduce the Page-Thorne profile "
                 "(a = 0..0.998), zero at ISCO, peak = 1\n";
    return true;
}

namespace {
// Mirrors blackbodyRadiance() in bh3d_physics.glsl: log-mapped lookup with
// linear texture filtering and clamp-to-edge.
std::array<double, 3> sampleRadianceLut(const std::vector<float>& tab, double T) {
    using namespace bh3d::physics;
    const int n = int(tab.size() / 3);
    const double t = std::log(std::clamp(T, kRadianceLutTmin, kRadianceLutTmax) / kRadianceLutTmin) /
                     std::log(kRadianceLutTmax / kRadianceLutTmin);
    const double f = t * (n - 1);
    const int i0 = std::min(int(f), n - 2);
    const double w = f - i0;
    std::array<double, 3> out{};
    for (int k = 0; k < 3; ++k)
        out[k] = tab[std::size_t(i0) * 3 + k] * (1.0 - w) + tab[std::size_t(i0 + 1) * 3 + k] * w;
    return out;
}
}  // namespace

bool testRadianceLut() {
    using namespace bh3d::physics;
    const std::vector<float> tab = radianceLutTable(1024);
    if (tab.size() != 1024u * 3u) {
        std::cerr << "[FAIL] radiance LUT: wrong size " << tab.size() << "\n";
        return false;
    }
    for (double T : {400.0, 1000.0, 6500.0, 12345.0, 38000.0, 2.0e5, 1.0e7}) {
        const auto lut = sampleRadianceLut(tab, T);
        const auto ref = blackbodyLinear(T);
        const double yl = linearLuminance(lut), yr = linearLuminance(ref);
        if (std::abs(yl - yr) > 0.01 * yr + 1e-12) {
            std::cerr << "[FAIL] radiance LUT at " << T << " K: luminance " << yl
                      << " vs direct " << yr << "\n";
            return false;
        }
    }
    // Exposure maps the unshifted hottest annulus to luminance 1.
    for (double peak : {5000.0, 38000.0, 1.0e7}) {
        const double y = linearLuminance(sampleRadianceLut(tab, peak)) * sciExposure(peak);
        if (std::abs(y - 1.0) > 0.01) {
            std::cerr << "[FAIL] sciExposure: peak " << peak << " K maps to " << y << ", expected 1\n";
            return false;
        }
    }
    // The shader must use the same range and size as the CPU table.
    const std::string glsl = readTextFile(std::string(AETHERION_SHADER_DIR) + "/bh3d_physics.glsl");
    if (glsl.find("clamp(T, 300.0, 1e8)") == std::string::npos ||
        glsl.find("log(1e8 / 300.0)") == std::string::npos ||
        glsl.find("1023.0") == std::string::npos || glsl.find("1024.0") == std::string::npos) {
        std::cerr << "[FAIL] radiance LUT: bh3d_physics.glsl no longer matches the 300 K..1e8 K, "
                     "1024-sample mapping\n";
        return false;
    }
    std::cout << "[PASS] radiance LUT: log-mapped lookup within 1% of direct blackbody, "
                 "exposure normalises peak to 1, shader mapping matches\n";
    return true;
}

bool testScientificDiskEmission() {
    using namespace bh3d::physics;
    // Rebuilds diskEmissionPhysical() on the CPU and checks the physical signatures.
    const double spin = 0.38, peak = 38000.0;
    const NTSetup s = ntSetup(spin);
    const std::vector<float> tab = radianceLutTable(1024);
    const double expo = sciExposure(peak);
    auto emit = [&](double r, double Lz, bool doppler) {
        const double T = peak * ntProfileFromSetup(r, s);
        const double g = diskRedshift(r, spin, Lz, doppler);
        auto c = sampleRadianceLut(tab, g * T);
        for (double& v : c) v *= expo;
        return c;
    };
    const double r = 4.0;
    const auto face = emit(r, 0.0, true);
    const auto appr = emit(r, -4.0, true);   // approaching side: blueshifted + beamed
    const auto rece = emit(r, +4.0, true);
    if (!(linearLuminance(appr) > linearLuminance(face) && linearLuminance(face) > linearLuminance(rece))) {
        std::cerr << "[FAIL] scientific emission: beaming should make approach > face-on > recede\n";
        return false;
    }
    if (!(appr[2] / appr[0] > rece[2] / rece[0])) {
        std::cerr << "[FAIL] scientific emission: approaching side should be bluer\n";
        return false;
    }
    const auto off = emit(r, -4.0, false);
    if (std::abs(linearLuminance(off) - linearLuminance(face)) > 1e-6) {
        std::cerr << "[FAIL] scientific emission: Doppler off must ignore Lz\n";
        return false;
    }
    // The disk is dark at and inside the ISCO and brightest well outside it.
    if (linearLuminance(emit(s.isco * 0.999, 0.0, true)) > 1e-6 ||
        !(linearLuminance(emit(s.isco * 1.4, 0.0, true)) > linearLuminance(emit(30.0, 0.0, true)))) {
        std::cerr << "[FAIL] scientific emission: radial brightness profile wrong\n";
        return false;
    }
    std::cout << "[PASS] scientific disk emission: beamed/blue approach side, dark at ISCO, "
                 "Doppler-off ignores Lz\n";
    return true;
}

#ifdef AETHERION_HAVE_CATALOG
// Round-trip: the JSON catalog must reconstruct the compiled-in profiles
// exactly. Keeps the two in lockstep until the C++ profiles are retired.
bool testCatalogRoundTrip() {
#ifndef AETHERION_CATALOG_DIR
#  define AETHERION_CATALOG_DIR "src/common/resources/catalog"
#endif
    const std::string dir = AETHERION_CATALOG_DIR;

    auto feq = [](float a, float b) {
        return std::fabs(a - b) <= 1e-3 * std::max(1.0f, std::fabs(a));
    };

    auto lr = catalog::loadDir(dir);
    if (!lr.ok()) {
        std::cerr << "[FAIL] catalog round-trip: no entries loaded from " << dir << "\n";
        return false;
    }
    for (const auto& err : lr.errors)
        std::cerr << "[WARN] catalog round-trip: skipped " << err << "\n";

    const auto builtin = profiles::allProfiles();
    if (lr.entries.size() != builtin.size()) {
        std::cerr << "[FAIL] catalog round-trip: " << lr.entries.size()
                  << " JSON entries vs " << builtin.size() << " built-in profiles\n";
        return false;
    }

    bool ok = true;
    auto fail = [&](size_t i, const std::string& what) {
        std::cerr << "[FAIL] catalog round-trip [" << i << "] "
                  << builtin[i].name << ": " << what << "\n";
        ok = false;
    };

    for (size_t i = 0; i < builtin.size(); ++i) {
        const BlackHoleProfile got = bh3d_catalog::toProfile(lr.entries[i]);
        const BlackHoleProfile& exp = builtin[i];

        if (got.name != exp.name) { fail(i, "name mismatch: " + got.name); continue; }
        if (got.description != exp.description) fail(i, "description mismatch");
        if (std::abs(got.massSolar - exp.massSolar) > 1e-3 * std::max(1.0, std::abs(exp.massSolar)))
            fail(i, "massSolar mismatch");

        const auto& gc = got.config; const auto& ec = exp.config;
        if (!feq(gc.blackHole.spinParameter, ec.blackHole.spinParameter)) fail(i, "spin mismatch");
        if (!feq(gc.disk.innerRadius, ec.disk.innerRadius)) fail(i, "disk.innerRadius");
        if (!feq(gc.disk.outerRadius, ec.disk.outerRadius)) fail(i, "disk.outerRadius");
        if (!feq(gc.disk.halfThickness, ec.disk.halfThickness)) fail(i, "disk.halfThickness");
        if (!feq(gc.disk.peakTemp, ec.disk.peakTemp)) fail(i, "disk.peakTemp");
        if (!feq(gc.disk.displayTempInner, ec.disk.displayTempInner)) fail(i, "disk.displayTempInner");
        if (!feq(gc.disk.displayTempOuter, ec.disk.displayTempOuter)) fail(i, "disk.displayTempOuter");
        if (!feq(gc.disk.saturationBoostInner, ec.disk.saturationBoostInner)) fail(i, "disk.satInner");
        if (!feq(gc.disk.saturationBoostOuter, ec.disk.saturationBoostOuter)) fail(i, "disk.satOuter");
        if (!feq(gc.jet.radius, ec.jet.radius)) fail(i, "jet.radius");
        if (!feq(gc.jet.length, ec.jet.length)) fail(i, "jet.length");
        if (!feq(gc.jet.color.r, ec.jet.color.r) || !feq(gc.jet.color.g, ec.jet.color.g) ||
            !feq(gc.jet.color.b, ec.jet.color.b)) fail(i, "jet.color");
        if (!feq(gc.blr.innerRadius, ec.blr.innerRadius)) fail(i, "blr.innerRadius");
        if (!feq(gc.blr.outerRadius, ec.blr.outerRadius)) fail(i, "blr.outerRadius");
        if (!feq(gc.blr.thickness, ec.blr.thickness)) fail(i, "blr.thickness");
        if (!feq(gc.blr.strength, ec.blr.strength)) fail(i, "blr.strength");
        if (!feq(gc.orbital.semiMajor, ec.orbital.semiMajor)) fail(i, "orbital.semiMajor");
        if (!feq(gc.orbital.bodyRadius, ec.orbital.bodyRadius)) fail(i, "orbital.bodyRadius");
        if (!feq(gc.bloom.threshold, ec.bloom.threshold)) fail(i, "bloom.threshold");
        if (!feq(gc.bloom.softKnee, ec.bloom.softKnee)) fail(i, "bloom.softKnee");
        if (!feq(gc.bloom.intensity, ec.bloom.intensity)) fail(i, "bloom.intensity");
        if (!feq(gc.bloom.exposure, ec.bloom.exposure)) fail(i, "bloom.exposure");
        if (!feq(gc.camera.initialPos.x, ec.camera.initialPos.x) ||
            !feq(gc.camera.initialPos.y, ec.camera.initialPos.y) ||
            !feq(gc.camera.initialPos.z, ec.camera.initialPos.z)) fail(i, "camera.initialPos");

        if (got.defaultJets       != exp.defaultJets)       fail(i, "defaultJets");
        if (got.defaultBLR        != exp.defaultBLR)        fail(i, "defaultBLR");
        if (got.defaultOrbBody    != exp.defaultOrbBody)    fail(i, "defaultOrbBody");
        if (got.defaultDoppler    != exp.defaultDoppler)    fail(i, "defaultDoppler");
        if (got.defaultHostGalaxy != exp.defaultHostGalaxy) fail(i, "defaultHostGalaxy");
        if (got.defaultLAB        != exp.defaultLAB)        fail(i, "defaultLAB");
        if (got.defaultCGM        != exp.defaultCGM)        fail(i, "defaultCGM");

        if (got.isBinaryWithBarycenter != exp.isBinaryWithBarycenter) fail(i, "barycentric");
        if (std::abs(got.companionMassSolar - exp.companionMassSolar) > 1e-6) fail(i, "companionMassSolar");

        if (got.galaxyBodies.size() != exp.galaxyBodies.size()) {
            fail(i, "galaxyBodies count");
            continue;
        }
        for (size_t b = 0; b < exp.galaxyBodies.size(); ++b) {
            const auto& gb = got.galaxyBodies[b];
            const auto& eb = exp.galaxyBodies[b];
            if (gb.type != eb.type) fail(i, "body type");
            if (!feq(gb.semiMajorRs, eb.semiMajorRs)) fail(i, "body semiMajorRs");
            if (!feq(gb.eccentricity, eb.eccentricity)) fail(i, "body ecc");
            if (!feq(gb.inclination, eb.inclination)) fail(i, "body inclination");
            if (gb.label != eb.label) fail(i, "body label");
            if (!feq(gb.massSolar, eb.massSolar)) fail(i, "body massSolar");
        }
        if (got.hasSecondaryBH != exp.hasSecondaryBH) fail(i, "hasSecondaryBH");
        if (exp.hasSecondaryBH &&
            !feq(got.secondaryConfig.disk.innerRadius, exp.secondaryConfig.disk.innerRadius))
            fail(i, "secondary disk.innerRadius");
    }

    if (ok)
        std::cout << "[PASS] catalog round-trip: " << builtin.size()
                  << " JSON profiles match built-in profiles field-for-field\n";
    return ok;
}

// Step 8: the OJ 287 orbit (a = 24.2 Rs, e = 0.65) must advance its periapsis
// by ~38.6 deg per orbit through the built-in 1PN precession formula.
bool testOrbitPrecession() {
    cfg::OrbitalConfig oc;
    oc.semiMajor    = 24.2f;
    oc.eccentricity = 0.65f;
    oc.inclination  = 0.0f;   // keep the orbit in the XZ plane

    OrbitalBody body(oc);
    body.setTimeScale(1.0f);  // one unit of sim-time per dt

    // Period in these units (GM = 0.5): T = 2*pi*sqrt(2 a^3).
    const double a = oc.semiMajor;
    const double period = 2.0 * M_PI * std::sqrt(2.0 * a * a * a);

    // Measure the periapsis direction at two successive periapsis passages. The
    // azimuthal angle at periapsis equals the accumulated precession (true
    // anomaly is 0 there), so the gap between consecutive periapses is exactly
    // the precession per orbit -- robust to any period rounding, unlike reading
    // the angle after a fixed number of steps (fragile at high eccentricity).
    const int stepsPerOrbit = 50000;
    const float dt = (float)(period / stepsPerOrbit);
    const int totalSteps = stepsPerOrbit * 22 / 10;   // ~2.2 orbits

    double rPrev = 1e300, rCur = 1e300, anglePrev = 0.0;
    double peri1 = 1e300, peri2 = 1e300;
    int found = 0;
    for (int i = 0; i < totalSteps && found < 2; ++i) {
        body.update(dt);
        const glm::vec3 p = body.position();
        const double r = std::sqrt((double)p.x * p.x + (double)p.z * p.z);
        const double ang = std::atan2((double)p.z, (double)p.x);
        // rPrev (i-2), rCur (i-1), r (i): a local minimum sits at i-1. Require
        // rCur < a so only true periapses count -- apoapsis is nearly flat and
        // float jitter there would otherwise register spurious minima.
        if (i >= 2 && rCur < rPrev && rCur < r && rCur < (double)oc.semiMajor) {
            if (found == 0) peri1 = anglePrev;
            else            peri2 = anglePrev;
            ++found;
        }
        rPrev = rCur; rCur = r; anglePrev = ang;
    }

    if (found < 2) {
        std::cerr << "[FAIL] orbit precession: only found " << found
                  << " periapsis passages in ~2.2 orbits\n";
        return false;
    }

    double deltaDeg = (peri2 - peri1) * 180.0 / M_PI;
    while (deltaDeg < 0.0)   deltaDeg += 360.0;
    while (deltaDeg >= 360.0) deltaDeg -= 360.0;

    if (std::abs(deltaDeg - 38.6) > 0.5) {
        std::cerr << "[FAIL] orbit precession: expected ~38.6 deg/orbit at a=24.2, e=0.65; got "
                  << deltaDeg << " deg\n";
        return false;
    }
    std::cout << "[PASS] orbit precession: " << deltaDeg
              << " deg/orbit at a=24.2, e=0.65 (target 38.6)\n";
    return true;
}

// Step 5: loading oj_287.json must auto-resolve the disk inner radius to the
// prograde Kerr ISCO for the primary spin (a*=0.38 -> ~2.344 Rs).
bool testISCOAuto() {
    const std::string dir = AETHERION_CATALOG_DIR;
    auto lr = catalog::loadDir(dir);
    if (!lr.ok()) {
        std::cerr << "[FAIL] ISCO auto: no catalog entries loaded from " << dir << "\n";
        return false;
    }

    const catalog::Entry* oj = nullptr;
    for (const auto& e : lr.entries)
        if (e.name == "OJ 287") { oj = &e; break; }
    if (!oj) {
        std::cerr << "[FAIL] ISCO auto: OJ 287 entry not found in catalog\n";
        return false;
    }

    const BlackHoleProfile p = bh3d_catalog::toProfile(*oj);
    const double expected = bh3d::physics::computeISCO(0.38);   // ~2.3431
    if (std::abs(p.config.disk.innerRadius - expected) > 5e-3) {
        std::cerr << "[FAIL] ISCO auto: OJ 287 inner radius expected ~" << expected
                  << " Rs, got " << p.config.disk.innerRadius << "\n";
        return false;
    }
    std::cout << "[PASS] ISCO auto: OJ 287 inner radius = " << p.config.disk.innerRadius
              << " Rs (computeISCO(0.38))\n";

    // Secondary resolves "isco" from its OWN spin (a*=0.5 -> ~2.117 Rs).
    const double expectedSec = bh3d::physics::computeISCO(0.5);
    if (!p.hasSecondaryBH ||
        std::abs(p.secondaryConfig.disk.innerRadius - expectedSec) > 5e-3 ||
        std::abs(expectedSec - 2.117) > 5e-3) {
        std::cerr << "[FAIL] ISCO auto: OJ 287 secondary inner radius expected ~" << expectedSec
                  << " Rs, got " << p.secondaryConfig.disk.innerRadius << "\n";
        return false;
    }
    std::cout << "[PASS] ISCO auto: OJ 287 secondary inner radius = "
              << p.secondaryConfig.disk.innerRadius << " Rs (own spin 0.5)\n";

    // Step 7c: the secondary BH body carries its mass for true-scale horizons.
    bool haveMass = false;
    for (const auto& b : p.galaxyBodies)
        if (b.type == GalaxyBody3DType::BlackHole && std::abs(b.massSolar - 1.5e8f) < 1e3f)
            haveMass = true;
    if (!haveMass) {
        std::cerr << "[FAIL] OJ 287 secondary BlackHole body is missing massSolar=1.5e8\n";
        return false;
    }
    std::cout << "[PASS] OJ 287 secondary body massSolar = 1.5e8 (true-scale ratio "
              << 1.5e8 / p.massSolar << ")\n";

    // The mass must survive into the OrbitalBody the snapshot reads from.
    for (const auto& b : p.galaxyBodies) {
        if (b.type != GalaxyBody3DType::BlackHole) continue;
        const OrbitalBody ob(profiles::makeOrbitalBody(b));
        if (ob.bodyType() != static_cast<int>(GalaxyBody3DType::BlackHole) ||
            std::abs(ob.bodyMass() - 1.5e8f) > 1e3f) {
            std::cerr << "[FAIL] makeOrbitalBody dropped the black-hole mass: " << ob.bodyMass() << "\n";
            return false;
        }
    }
    std::cout << "[PASS] black-hole mass carried through makeOrbitalBody -> OrbitalBody\n";
    return true;
}
#endif  // AETHERION_HAVE_CATALOG

}  // namespace

int main() {
    bool ok = true;

    ok = testSolarGrazingDeflection() && ok;
    ok = testTimelikeConservationLongRun() && ok;
    ok = testPhotonCaptureThreshold() && ok;
    ok = testPeriapsisTurningPoint() && ok;
    ok = testStrongDeflectionGrowth() && ok;

    // Accretion-disk colour physics, shared with the GPU shader's LUT.
    ok = testBlackbodyMonotonicBlueness() && ok;
    ok = testBlackbodyReferencePoints() && ok;
    ok = testNovikovThorneProfile() && ok;
    ok = testComputeISCO() && ok;
    ok = testDiskColourGradient() && ok;

    // Beta 3 engine: exact geodesic integrator + disk physics (shared shader math).
    ok = testGeodesicShadowEdge() && ok;
    ok = testGeodesicWeakDeflection() && ok;
    ok = testGeodesicMatches2D() && ok;
    ok = testLowQualityShadow() && ok;
    ok = testDiskRedshift() && ok;
    ok = testBlackbodyLinearRadiance() && ok;

    // Step 0: shared GLSL injection + no duplicate declarations across shaders.
    ok = testShaderInjection() && ok;
    ok = testSharedShaderNoDuplicates() && ok;

    // Steps 2 + 3: exact Kerr redshift, NT constants, radiance LUT.
    ok = testNTSetupMatchesReference() && ok;
    ok = testRadianceLut() && ok;
    ok = testScientificDiskEmission() && ok;

#ifdef AETHERION_HAVE_CATALOG
    // JSON catalog must reconstruct the compiled-in profiles exactly.
    ok = testCatalogRoundTrip() && ok;
    // Beta 3 engine: auto-ISCO inner radius and GR orbital precession.
    ok = testISCOAuto() && ok;
    ok = testOrbitPrecession() && ok;
#endif

    if (!ok) {
        std::cerr << "Physics regression tests FAILED\n";
        return EXIT_FAILURE;
    }

    std::cout << "Physics regression tests PASSED\n";
    return EXIT_SUCCESS;
}
