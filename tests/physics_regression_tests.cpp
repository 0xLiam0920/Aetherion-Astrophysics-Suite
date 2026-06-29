#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "src/2D/2D-physics/geodesic.hpp"
#include "src/2D/2D-physics/schwarzschild.hpp"
#include "src/2D/2D-physics/units.hpp"
#include "src/2D/2D-simulation/photon.hpp"

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

} // namespace

int main() {
    bool ok = true;

    ok = testSolarGrazingDeflection() && ok;
    ok = testTimelikeConservationLongRun() && ok;
    ok = testPhotonCaptureThreshold() && ok;
    ok = testPeriapsisTurningPoint() && ok;
    ok = testStrongDeflectionGrowth() && ok;

    if (!ok) {
        std::cerr << "Physics regression tests FAILED\n";
        return EXIT_FAILURE;
    }

    std::cout << "Physics regression tests PASSED\n";
    return EXIT_SUCCESS;
}
