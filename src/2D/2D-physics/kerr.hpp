#pragma once
// ─────────────────────────────────────────────────────────────────────────────
// kerr.hpp — equatorial Kerr null geodesic integrator (batch 2 #4).
//
// Restricted to the equatorial plane (θ = π/2), so the Carter constant Q = 0
// and only the (E, L_z) conserved pair remains. We integrate in Boyer-Lindquist
// (t, r, φ) with affine parameter λ. Geometric units G = c = M = 1 are used
// internally; the caller passes a Schwarzschild handle so the mass scale lines
// up with the rest of the simulation, and an `a` parameter for the spin in
// units of M (a ∈ [0, 1)).
//
// Equations of motion for equatorial null geodesics (see Bardeen 1973,
// Chandrasekhar 1983 §63):
//
//   Σ = r²            (in the equatorial plane)
//   Δ = r² − 2Mr + a²
//
//   dt/dλ  = (1/Δ) · [ (r² + a² + 2Ma²/r) · E − (2Ma/r) · L ]
//   dφ/dλ  = (1/Δ) · [ (2Ma/r) · E + (1 − 2M/r) · L ]
//   (dr/dλ)² = E² + (2M/r)(a·E − L)² / r² + (a²·E² − L²) / r²
//            = R(r) / r⁴                          (Carter form, with Q=0)
//
// We integrate r with a sign flag for dr/dλ that flips at turning points
// (R(r) = 0), which is the standard trick for getting orbits through pericenter
// without dropping accuracy. φ and t advance monotonically.
// ─────────────────────────────────────────────────────────────────────────────

#include "schwarzschild.hpp"
#include "integrator.hpp"
#include <cmath>
#include <algorithm>

struct Kerr {
    double M = 1.0;
    double a = 0.0;   // spin in geometric units (|a| < M)

    double horizonOuter() const {
        const double disc = std::max(0.0, M * M - a * a);
        return M + std::sqrt(disc);
    }
    double horizonInner() const {
        const double disc = std::max(0.0, M * M - a * a);
        return M - std::sqrt(disc);
    }
    // Equatorial ergosphere outer radius: r_ergo(θ=π/2) = 2M (independent of spin).
    double ergosphereEquatorial() const { return 2.0 * M; }
    double delta(double r) const { return r * r - 2.0 * M * r + a * a; }

    // R(r) for an equatorial null geodesic with conserved (E, L). Sign of R
    // gives the allowed region for r-motion; zero = turning point.
    double R_radial(double r, double E, double L) const {
        const double term1 = E * E * (r * r * r * r + a * a * r * r + 2.0 * M * a * a * r);
        const double term2 = -4.0 * M * a * E * L * r;
        const double term3 = -(r * r - 2.0 * M * r) * L * L;
        return term1 + term2 + term3;
    }
};

// State along an affine-parameterised equatorial null geodesic.
struct KerrNullState {
    double r;
    double phi;
    double t;          // coordinate time (geometric units)
    double pr_sign;    // +1 outgoing, −1 ingoing (flips at turning points)
};

// Integrate one equatorial null geodesic from (r0, phi0) with conserved
// (E, L). The φ direction is set by the sign of L (Bardeen convention:
// L > 0 = prograde with respect to BH spin). Stops on horizon capture
// (r ≤ r+·1.001) or on r > rMax.
//
// Returns sample points (x, y) in equatorial-plane Cartesian (x = r·cosφ,
// y = r·sinφ) suitable for a polyline renderer.
struct KerrPathResult {
    std::vector<float> verts;    // alternating x,y pairs
    bool   captured = false;
    int    steps    = 0;
    double finalR   = 0.0;
};

inline KerrPathResult integrateKerrEquatorial(
    const Kerr&  kerr,
    double r0, double phi0,
    double E,  double L,
    int    initial_pr_sign = -1,         // -1 = ingoing (most rays from infinity)
    double rMax    = 1.0e3,
    double dlam    = 0.25,
    int    maxStep = 4000)
{
    KerrPathResult out;
    out.verts.reserve(static_cast<size_t>(maxStep) * 2);
    KerrNullState s {r0, phi0, 0.0, double(initial_pr_sign < 0 ? -1 : 1)};

    const double rH = kerr.horizonOuter() * 1.001;

    auto pushVert = [&](double r, double phi) {
        out.verts.push_back(static_cast<float>(r * std::cos(phi)));
        out.verts.push_back(static_cast<float>(r * std::sin(phi)));
    };
    pushVert(s.r, s.phi);

    for (int i = 0; i < maxStep; ++i) {
        if (s.r <= rH) { out.captured = true; break; }
        if (s.r > rMax) break;

        const double r  = s.r;
        const double Δ  = std::max(1e-9, kerr.delta(r));
        const double r2 = r * r;

        // dt/dλ and dφ/dλ from the standard equatorial Kerr null geodesic.
        const double a  = kerr.a;
        const double Mm = kerr.M;
        const double dphi_dl = ( (2.0 * Mm * a / r) * E + (1.0 - 2.0 * Mm / r) * L ) / Δ;
        // |dr/dλ|² = R(r) / r⁴
        const double R = kerr.R_radial(r, E, L);
        double dr_abs = (R > 0.0) ? std::sqrt(R) / r2 : 0.0;
        if (R <= 0.0) {
            // turning point reached → flip the radial direction
            s.pr_sign = -s.pr_sign;
            dr_abs = 1e-6;   // nudge so we don't get stuck
        }
        const double dr_dl = s.pr_sign * dr_abs;

        s.r   += dr_dl   * dlam;
        s.phi += dphi_dl * dlam;
        // (we don't render t, but advance it so callers could compute it later)
        const double dt_dl = ( (r2 + a*a + 2.0*Mm*a*a/r) * E - (2.0*Mm*a/r) * L ) / Δ;
        s.t   += dt_dl   * dlam;

        pushVert(s.r, s.phi);
        ++out.steps;
    }
    out.finalR = s.r;
    return out;
}

// Build a fan of Kerr equatorial null geodesics sampled across impact
// parameters b ∈ [−bMax, +bMax] coming in from large r0, matched to the
// canonical Schwarzschild conventions used elsewhere in the simulation:
//   E = 1, L = b · E    (for an ingoing photon at infinity)
//   initial pr_sign = −1
// This is the Kerr analogue of `Simulation::rebuildPhotons()` and is intended
// for visual comparison with the Schwarzschild case.
inline std::vector<KerrPathResult> sweepKerrEquatorial(
    const Kerr& kerr, int nRays, double bMax, double r0 = 60.0)
{
    std::vector<KerrPathResult> out;
    out.reserve(static_cast<size_t>(nRays));
    for (int i = 0; i < nRays; ++i) {
        const double t = (nRays <= 1) ? 0.5 : double(i) / double(nRays - 1);
        const double b = (2.0 * t - 1.0) * bMax;
        const double E = 1.0;
        const double L = b * E;
        // Aim the photon roughly along −x at +y = b at r0.
        const double phi0 = std::atan2(b, -std::sqrt(std::max(0.0, r0 * r0 - b * b)));
        out.push_back(integrateKerrEquatorial(kerr, r0, phi0, E, L, -1));
    }
    return out;
}
