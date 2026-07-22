#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "src/2D/2D-physics/geodesic.hpp"
#include "src/2D/2D-physics/schwarzschild.hpp"
#include "src/2D/2D-physics/units.hpp"
#include "src/2D/2D-simulation/photon.hpp"
#include "src/3D/bh3d_blackbody.hpp"

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
    const double rISCO = 2.0;

    // The peak should land at r = (49/36) times rISCO, and the profile is
    // normalised so it hits exactly 1.0 right there.
    const double rPeak = (49.0 / 36.0) * rISCO;
    const double peak = novikovThorneTemperature(rPeak, rISCO);
    if (std::abs(peak - 1.0) > 1e-3) {
        std::cerr << "[FAIL] NT profile peak should be 1.0 at r=" << rPeak
                  << "; got " << peak << "\n";
        return false;
    }

    // Since the ISCO is a zero-torque inner boundary, the temperature should
    // basically vanish right there.
    if (novikovThorneTemperature(rISCO, rISCO) > 1e-6) {
        std::cerr << "[FAIL] NT profile should be ~0 at the ISCO; got "
                  << novikovThorneTemperature(rISCO, rISCO) << "\n";
        return false;
    }

    // Beyond the peak it should decrease monotonically outward, and stay
    // bounded within [0,1] the whole way.
    double prev = 2.0;
    for (double r = rPeak; r <= 40.0; r += 0.5) {
        const double t = novikovThorneTemperature(r, rISCO);
        if (t < -1e-9 || t > 1.0 + 1e-6) {
            std::cerr << "[FAIL] NT profile out of [0,1] at r=" << r << ": " << t << "\n";
            return false;
        }
        if (t > prev + 1e-6) {
            std::cerr << "[FAIL] NT profile not decreasing outward at r=" << r
                      << " (" << t << " > " << prev << ")\n";
            return false;
        }
        prev = t;
    }

    std::cout << "[PASS] Novikov-Thorne profile: peak=1, zero at ISCO, monotone falloff\n";
    return true;
}

bool testComputeISCO() {
    using bh3d::physics::computeISCO;

    if (std::abs(computeISCO(0.0) - 3.0) > 1e-9) {
        std::cerr << "[FAIL] Schwarzschild ISCO should be 3 Rs; got "
                  << computeISCO(0.0) << "\n";
        return false;
    }
    // The prograde ISCO should keep shrinking as spin increases.
    double prev = 4.0;
    for (double a = 0.0; a <= 0.99; a += 0.05) {
        const double r = computeISCO(a);
        if (r > prev + 1e-9) {
            std::cerr << "[FAIL] ISCO should decrease with spin; a=" << a
                      << " gave " << r << " > " << prev << "\n";
            return false;
        }
        if (r < 0.6 - 1e-9) {
            std::cerr << "[FAIL] ISCO fell below physical floor at a=" << a
                      << ": " << r << "\n";
            return false;
        }
        prev = r;
    }
    std::cout << "[PASS] ISCO(a*): 3 Rs at a=0, monotone shrink with spin\n";
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

        const double tInner = diskColorTemperature(p.inner, rPeak, rISCO);
        const double tOuter = diskColorTemperature(p.inner, p.rOut, rISCO);

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

    if (!ok) {
        std::cerr << "Physics regression tests FAILED\n";
        return EXIT_FAILURE;
    }

    std::cout << "Physics regression tests PASSED\n";
    return EXIT_SUCCESS;
}
