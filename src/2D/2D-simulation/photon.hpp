/*--------------- HEADERS ---------------*/
#pragma once
#include "../2D-physics/schwarzschild.hpp"
#include "../2D-physics/geodesic.hpp"
#include "../2D-core/types.hpp"
#include <vector>
#include <cmath>
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*---------------- main struct for photons ------------------*/
struct Photon {
    double impactParameter = 0.0;
    bool captured = false;
    std::vector<Vec2> path;      // downsampled for display
    std::vector<Vec2> fullPath;  // full-resolution for export

    // deflection tracking logic 
    double deflectionAngle = 0.0;   // total deflection in radians
    double deflectionDeg   = 0.0;   // total deflection in degrees
    double phiInfinity     = 0.0;   // angle swept periapsis as photon approaches infinity (used for deflection calculation)

    // Shapiro / coordinate-time bookkeeping (set by computePath).
    // coordTime  = ∫ dt along the integrated half-path (one side of periapsis)
    //              in geometric units (M); multiplied by 2 to cover the full
    //              symmetric trajectory. Diverges as r → ∞ (light still has
    //              finite path length but time is unbounded for unbounded r),
    //              so this is really only meaningful as a *relative* number
    //              between rays integrated to the same rMax.
    // shapiroDelay = ∫ (1/f(r) − 1) · r²/|b| dφ across the full path; this is
    //                the *excess* coordinate time over a hypothetical flat
    //                traversal of the same φ sweep. Converges (integrand → 0
    //                as r → ∞), so it is the physically meaningful quantity
    //                to display per-ray. Geometric units (M).
    double coordTime    = 0.0;
    double shapiroDelay = 0.0;

    // Closest approach radius along the integrated path, in geometric units
    // (M). For non-captured rays this equals the analytic periapsis from
    // findPeriapsis(b); for captured rays it is the smallest r reached
    // before the horizon-crossing termination. Surfaced in the data-analysis
    // panel so photon-sphere sweeps can show how close each ray got.
    double rMin_M = 0.0;

    // Hot-spot / disk-source emission flag. When true the path was built by
    // computeEmissionPath() starting at (r_emit, phi_emit) with a local
    // emission angle, rather than from infinity via the impact parameter.
    // Affects nothing in the integrator, but lets the renderer style the ray
    // differently and lets the deflection/Shapiro analyses skip it.
    bool   fromEmitter   = false;
    double emitRadius_M  = 0.0;  // r_emit / M at spawn (display only)

    // ── Wavelength tagging (batch 2 #1) ─────────────────────────────────────
    // Rest-frame emitted wavelength (set by the caller); the observed value
    // after gravitational redshift along the path is filled in by
    // computeEmissionPath() (for emitter rays) or left equal to the rest value
    // for background sweeps. Stored in nanometres so the spectrum panel can
    // bin them directly. Default 550nm = green / visible reference.
    double wavelength_nm     = 550.0;
    double wavelength_obs_nm = 550.0;

    // ── Linear polarization (batch 2 #5) ────────────────────────────────────
    // Stokes I/Q/U in the photon's local frame. V is omitted (no circular
    // polarization handling). The polarizationAngle field stores the angle
    // of the linear-polarization axis (radians, [0, π)) at the observer end
    // of the ray, after parallel-transport along the null geodesic.
    // In Schwarzschild, in the equatorial plane, the Walker-Penrose constant
    // simplifies to: the angle between the polarization vector and the local
    // (e_r, e_φ) frame is parallel-transported, which (for a planar geodesic)
    // means the observed polarization-angle rotation equals half the total
    // azimuthal sweep relative to the emission direction.
    double stokesI = 1.0;
    double stokesQ = 1.0;   // default = fully linearly polarized along emission radial
    double stokesU = 0.0;
    double polarizationAngle = 0.0;   // radians, [0, π)

    void computePath(const Schwarzschild& bh, double rMax = 1e5, double dphi = 0.0025) { // logic here is a little scuffed, 
        path.clear(); // but basically we integrate the photon's path in small angular steps until it either escapes (r → ∞) or is captured (crosses the event horizon). 
        fullPath.clear();
        captured = false; // We track the path points for rendering, and we also keep track of the deflection angle and the angle at which the photon escapes to infinity for analysis. 
        // The integration uses the null geodesic equations in Schwarzschild spacetime, and we have to be careful to handle the asymptotic behavior correctly to get accurate deflection angles.

        // Scale the integration cutoff with M so the loop doesn't bail on step 0 at preset scales (M ~ 1e9 m makes the default 1e5 dwarfed by typical periapsis).
        rMax = std::max(rMax, 200.0 * bh.M);

        double b = std::abs(impactParameter);
        if (b <= bh.criticalImpact()) { captured = true; rMin_M = bh.horizon() / bh.M; return; }

        double rmin = bh.findPeriapsis(b);
        if (rmin <= 0.0) { captured = true; rMin_M = bh.horizon() / bh.M; return; }

        // Record closest approach for the data panel. For non-captured rays
        // this is exactly the analytic periapsis; for captured rays we
        // overwrite below as the photon spirals inward.
        rMin_M = rmin / bh.M;

        // Initial conditions at periapsis (phi=0, du/dphi=0)
        GeodesicState state = { 1.0 / rmin, 0.0 };
        double phi = 0.0;

        // Integrate outgoing half from periapsis
        std::vector<std::pair<double,double>> half;  // (r, phi)
        half.emplace_back(1.0 / state.u, phi);

        // Track last valid outgoing state for asymptotic extrapolation.
        // This is always the most recent step where u > 0 and the photon
        // is heading outward (du/dφ < 0). 
        double ext_phi = phi;
        double ext_u   = state.u;
        double ext_du  = state.du_dphi;

        const double horizonU = 1.0 / (bh.horizon() * 1.001); // capture threshold
        const int maxSteps = 200000; // Call it overenginering, call it paranoia, but I want to avoid infinite loops, thank you very much.
        auto startTime = std::chrono::steady_clock::now();
        constexpr int BUDGET_CHECK_INTERVAL = 512;
        constexpr auto TIME_BUDGET = std::chrono::milliseconds(50);

        // Shapiro / coord-time accumulators (one-sided; doubled below for the full trajectory).
        // dt/dφ along a null geodesic in Schwarzschild = r² / (|b| · f(r)).
        // Subtracting the flat-space r²/|b| leaves the Shapiro excess, which converges as r → ∞.
        double t_half       = 0.0;
        double shapiro_half = 0.0;
        double r_prev       = 1.0 / state.u; // = rmin at φ=0
        double f_prev       = bh.f(r_prev);
        double abs_b        = std::max(1e-30, b);

        for (int i = 0; i < maxSteps; ++i) {
            // Check time budget periodically to avoid stalling the frame
            if ((i & (BUDGET_CHECK_INTERVAL - 1)) == 0 && i > 0) {
                auto now = std::chrono::steady_clock::now();
                if (now - startTime > TIME_BUDGET) break;
            }
            state = stepNullGeodesic(bh, state, dphi);
            phi += dphi;

            // NOTE: Escape: u → 0 means that r → ∞
            if (state.u <= 0.0) break;

            double r = 1.0 / state.u;

            // Trapezoidal accumulation of t and Shapiro-excess between (r_prev, r)
            // Skip while still inside / right at horizon to avoid the f→0 singularity blowing up.
            double f_curr = bh.f(r);
            if (f_prev > 1e-6 && f_curr > 1e-6) {
                double dt_curr = (r * r) / (abs_b * f_curr);
                double dt_prev = (r_prev * r_prev) / (abs_b * f_prev);
                t_half += 0.5 * (dt_prev + dt_curr) * dphi;

                double sh_curr = (1.0 / f_curr - 1.0) * (r * r) / abs_b;
                double sh_prev = (1.0 / f_prev - 1.0) * (r_prev * r_prev) / abs_b;
                shapiro_half += 0.5 * (sh_prev + sh_curr) * dphi;
            }
            r_prev = r;
            f_prev = f_curr;

            // NOTE: Capture: if the photon crosses the horizon (with a small safety margin) we consider it captured. We check this inside 
            // the loop because the photon could be on an unstable orbit and hover near the horizon for a while before finally crossing it, 
            // and we want to capture that behavior accurately.
            if (state.u >= horizonU) { captured = true; break; }

            // Save outgoing state for extrapolation (u > 0 guaranteed here,
            // du < 0 once photon has passed periapsis and is heading out).
            if (state.du_dphi < 0.0) {
                ext_phi = phi;
                ext_u   = state.u;
                ext_du  = state.du_dphi;
            }

            half.emplace_back(r, phi);

            // Escaped to max distance
            if (r >= rMax) break;
        }

        if (captured) return;

        // Extrapolate φ to r = ∞ using the asymptotic flat-space solution.
        // Far from the black hole the Binet equation reduces to d²u/dφ² + u ≈ 0
        // of which the solution is u = A cos(φ − φ₀).  The remaining angle from an
        // outgoing state (u, du/dφ) to the zero-crossing u = 0 (i.e. r → ∞) is
        // atan2(u, −du/dφ).  Using the last saved outgoing state makes this
        // sub-step accurate regardless of whether the loop terminated because
        // u crossed zero or because r exceeded rMax.


        double phi_inf = ext_phi + std::atan2(ext_u, -ext_du);

        if (!std::isfinite(phi_inf)) {
        // Asymptotic extrapolation failed (degenerate state near critical b).
        // Mark as captured so this ray is excluded from lensing analytics
        // and FITS export instead of poisoning them with NaN.
            captured = true;
            return;
        }
        double phi_offset = M_PI + phi_inf;

        // Deflection = 2·φ_∞ - π  (straight line sweeps π)
        phiInfinity     = phi_inf;
        deflectionAngle = 2.0 * phi_inf - M_PI;
        deflectionDeg   = deflectionAngle * 180.0 / M_PI;

        // Symmetric trajectory: outgoing half == incoming half (time-reversal symmetry),
        // so the full-path totals are just 2× the one-sided accumulators.

        coordTime    = std::isfinite(t_half)       ? 2.0 * t_half       : 0.0;
        shapiroDelay = std::isfinite(shapiro_half) ? 2.0 * shapiro_half : 0.0;

        // Build full-resolution path for export BEFORE downsampling.
        // Cap to MAX_EXPORT_HALF to prevent unbounded allocation on near-critical-impact photons.
        // Previously had no cap, so at 120 rays × 200K steps this would cause a memory leak nearing ~500 MB.
        // MAX_EXPORT_HALF gives 5× more detail than display while keeping memory bounded.
        // [FIXED 2026-04-24: added stride-sampled cap, was reserve(2*half.size()) with no limit]
        static constexpr int MAX_EXPORT_HALF = 5000;
        const std::vector<std::pair<double,double>>* exportSrc = &half;
        std::vector<std::pair<double,double>> exportSampled;
        if ((int)half.size() > MAX_EXPORT_HALF) {
            exportSampled.reserve(MAX_EXPORT_HALF);
            float stride = float(half.size() - 1) / float(MAX_EXPORT_HALF - 1);
            for (int s = 0; s < MAX_EXPORT_HALF; ++s)
                exportSampled.push_back(half[size_t(s * stride + 0.5f)]);
            exportSrc = &exportSampled;
        }
        fullPath.reserve(2 * exportSrc->size());
        for (int i = (int)exportSrc->size() - 1; i >= 0; --i) {
            double r = (*exportSrc)[i].first;
            double phi_local = -(*exportSrc)[i].second + phi_offset;
            fullPath.push_back({
                (float)(r * std::cos(phi_local)),
                (float)(r * std::sin(phi_local))
            });
        }
        for (size_t i = 0; i < exportSrc->size(); ++i) {
            double r = (*exportSrc)[i].first;
            double phi_local = (*exportSrc)[i].second + phi_offset;
            fullPath.push_back({
                (float)(r * std::cos(phi_local)),
                (float)(r * std::sin(phi_local))
            });
        }

        // Downsample half-path for display: cap to MAX_DISPLAY_HALF points to
        // prevent huge allocations for strongly-lensed near-critical-impact
        // photons (which can accumulate up to 200K integration steps).
        // I learned this the hard way by running it on Sgr A* and watching the process balloon.
        // do not remove this limit.
        static constexpr int MAX_DISPLAY_HALF = 1000; 
        if ((int)half.size() > MAX_DISPLAY_HALF) {
            std::vector<std::pair<double,double>> sampled;
            sampled.reserve(MAX_DISPLAY_HALF);
            float stride = float(half.size() - 1) / float(MAX_DISPLAY_HALF - 1);
            for (int s = 0; s < MAX_DISPLAY_HALF; ++s)
                sampled.push_back(half[size_t(s * stride + 0.5f)]);
            half = std::move(sampled);
        }
        path.reserve(2 * half.size()); // 

        // Build full path: incoming (mirror of outgoing) then outgoing
        for (int i = (int)half.size() - 1; i >= 0; --i) {
            double r = half[i].first;
            double phi_local = -half[i].second + phi_offset;
            path.push_back({
                (float)(r * std::cos(phi_local)),
                (float)(r * std::sin(phi_local))
            });
        }
        for (size_t i = 0; i < half.size(); ++i) { // for the size of the half vector, push back the full path of the photon, which is the mirror of the outgoing path followed by the outgoing path itself
            double r = half[i].first;
            double phi_local = half[i].second + phi_offset;
            path.push_back({
                (float)(r * std::cos(phi_local)),
                (float)(r * std::sin(phi_local))
            });
        }
    }

    /*------------- Hot-spot / disk-source emission --------------------*/
    // Build a null geodesic that starts at (r0, phi0) heading at local
    // emission angle alpha measured from the *outward radial* direction
    // (alpha=0 → straight out, alpha=π/2 → prograde tangent, alpha=−π/2 →
    // retrograde tangent). Integrates forward until the ray escapes to
    // r >= rMax or falls through the horizon. No mirroring/symmetry is
    // applied — this is a one-shot emitter, not a passing flyby.
    //
    // The conserved impact parameter for the resulting ray is
    //     b = r0 · sin(alpha) / sqrt(f(r0))
    // (derived from the local-tetrad components of the null tangent in
    //  the Schwarzschild static frame).  Storing it makes the existing
    //  capturedRayLine / colorByRedshift visualizers Just Work.
    void computeEmissionPath(const Schwarzschild& bh,
                             double r0, double phi0, double alpha,
                             double rMax = 1e5, double dphi = 0.002)
    {
        path.clear();
        fullPath.clear();
        captured     = false;
        fromEmitter  = true;
        emitRadius_M = r0 / std::max(1e-30, bh.M);

        // Same M-scaling as computePath: required so disk-emitter rays (O) render at preset scales.
        rMax = std::max(rMax, 200.0 * bh.M);

        // If we're emitting from inside the horizon, the photon is captured by definition.
        const double rH = bh.horizon();
        if (r0 <= rH * 1.001) { captured = true; return; }

        // Initial Binet state. du/dφ = -u · cot(alpha); sign of dφ direction tracked separately.
        const double s = std::sin(alpha);
        const double c = std::cos(alpha);
        const double u0 = 1.0 / r0;
        // Photon going straight outward (sin=0) has du/dφ undefined; integrate trivially
        // outward by skipping the geodesic stepper and just emitting a radial line.
        if (std::abs(s) < 1e-6) {
            // Radial null geodesic: pure straight line outward (or inward if c<0).
            const double dir = (c >= 0.0) ? 1.0 : -1.0;
            if (dir < 0.0) { captured = true; return; }
            const double cos_p = std::cos(phi0), sin_p = std::sin(phi0);
            const int N = 64;
            for (int i = 0; i <= N; ++i) {
                double r = r0 + (rMax - r0) * (double(i) / N);
                path.push_back({(float)(r * cos_p), (float)(r * sin_p)});
            }
            fullPath = path;
            // Coordinate-time integral along a radial null geodesic:
            //   dt/dr = 1/f(r) ⇒ ∫ = (rMax−r0) + 2M ln((rMax−2M)/(r0−2M))
            coordTime    = (rMax - r0) + 2.0 * bh.M *
                           std::log((rMax - 2.0*bh.M) / std::max(1e-30, r0 - 2.0*bh.M));
            shapiroDelay = 2.0 * bh.M *
                           std::log((rMax - 2.0*bh.M) / std::max(1e-30, r0 - 2.0*bh.M));
            impactParameter   = 0.0;
            deflectionAngle   = 0.0;
            deflectionDeg     = 0.0;
            return;
        }

        // Conserved b (positive). Direction (prograde/retrograde) handled via the sign of dφ.
        const double abs_b = r0 * std::abs(s) / std::sqrt(std::max(1e-12, bh.f(r0)));
        const double phiSign = (s > 0.0) ? 1.0 : -1.0;
        impactParameter = phiSign * abs_b;

        // du/dφ at emission: dr/dφ = r·cot(alpha) (with the photon moving outward when c>0)
        //   du/dφ = -u² · dr/dφ = -u · cot(alpha)
        // We always step in the direction of increasing φ (positive dφ for prograde,
        // negative for retrograde) — we do this by flipping the sign of the φ increment.
        GeodesicState state = { u0, -u0 * (c / s) };
        double phi_local = 0.0;
        double dphi_signed = phiSign * dphi;

        const double horizonU = 1.0 / (rH * 1.001);

        std::vector<std::pair<double,double>> traj; // (r, phi_local)
        traj.emplace_back(r0, 0.0);

        const int maxSteps = 200000;
        auto startTime = std::chrono::steady_clock::now();
        constexpr int BUDGET_CHECK_INTERVAL = 512;
        constexpr auto TIME_BUDGET = std::chrono::milliseconds(40);

        double r_prev = r0;
        double f_prev = bh.f(r0);
        double t_acc  = 0.0;
        double sh_acc = 0.0;

        for (int i = 0; i < maxSteps; ++i) {
            if ((i & (BUDGET_CHECK_INTERVAL - 1)) == 0 && i > 0) {
                auto now = std::chrono::steady_clock::now();
                if (now - startTime > TIME_BUDGET) break;
            }
            state = stepNullGeodesic(bh, state, dphi_signed);
            phi_local += dphi_signed;

            if (state.u <= 0.0) break; // escaped to r = ∞
            if (state.u >= horizonU) { captured = true; break; }

            double r = 1.0 / state.u;
            double f_curr = bh.f(r);
            if (f_prev > 1e-6 && f_curr > 1e-6) {
                double dt_curr = (r * r)         / (abs_b * f_curr);
                double dt_prev = (r_prev * r_prev) / (abs_b * f_prev);
                t_acc  += 0.5 * (dt_prev + dt_curr) * std::abs(dphi_signed);
                sh_acc += 0.5 *
                          ((1.0/f_prev - 1.0) * (r_prev*r_prev) / abs_b +
                           (1.0/f_curr - 1.0) * (r*r)           / abs_b)
                          * std::abs(dphi_signed);
            }
            r_prev = r;
            f_prev = f_curr;

            traj.emplace_back(r, phi_local);
            if (r >= rMax) break;
        }

        coordTime    = t_acc;
        shapiroDelay = sh_acc;

        // Decimate for display, copy to fullPath at higher density.
        static constexpr int MAX_DISPLAY = 800;
        static constexpr int MAX_EXPORT  = 4000;
        auto stamp = [&](std::vector<Vec2>& dst, int cap) {
            if ((int)traj.size() <= cap) {
                dst.reserve(traj.size());
                for (auto& p : traj) {
                    double phi_w = phi0 + p.second;
                    dst.push_back({(float)(p.first * std::cos(phi_w)),
                                   (float)(p.first * std::sin(phi_w))});
                }
            } else {
                dst.reserve(cap);
                float stride = float(traj.size() - 1) / float(cap - 1);
                for (int i = 0; i < cap; ++i) {
                    auto& p = traj[size_t(i * stride + 0.5f)];
                    double phi_w = phi0 + p.second;
                    dst.push_back({(float)(p.first * std::cos(phi_w)),
                                   (float)(p.first * std::sin(phi_w))});
                }
            }
        };
        stamp(path,     MAX_DISPLAY);
        stamp(fullPath, MAX_EXPORT);

        // Deflection is the change in φ between start and the asymptotic
        // direction; for an emitter we just report |Δφ| along the integrated
        // arc, since there is no "incoming" leg to mirror.
        deflectionAngle = phi_local;
        deflectionDeg   = phi_local * 180.0 / M_PI;

        // Gravitational redshift of the rest-frame wavelength:
        //   λ_obs / λ_emit = sqrt(f(r_obs)) / sqrt(f(r_emit))
        // For an observer at r → ∞ this reduces to 1/sqrt(f(r_emit)).
        // Sky-mapped emitter rays are by construction observed at infinity, so:
        const double fEmit = std::max(1e-12, bh.f(r0));
        wavelength_obs_nm = wavelength_nm / std::sqrt(fEmit);

        // Parallel-transported linear polarization (Walker-Penrose).
        // For an equatorial null geodesic in Schwarzschild, Q+iU rotates by
        // 2·Δψ where Δψ = ½·(swept azimuth) relative to the local radial.
        // Take the swept azimuth as |phi_local| (one-sided emission); the
        // observed polarization angle is the emission angle plus Δψ.
        const double psiEmit = 0.5 * std::atan2(stokesU, stokesQ);
        const double psiObs  = psiEmit + 0.5 * std::abs(phi_local);
        // Normalise into [0, π).
        double pa = std::fmod(psiObs, M_PI);
        if (pa < 0) pa += M_PI;
        polarizationAngle = pa;
        const double P = std::sqrt(stokesQ*stokesQ + stokesU*stokesU);
        stokesQ = P * std::cos(2.0 * polarizationAngle);
        stokesU = P * std::sin(2.0 * polarizationAngle);
    }
};  // fun fact! I spent 10 minutes trying to figure out why my code wouldn't compile, including rewriting a quarter of this file, only to find that everything
    // was caused by forgetting this semicolon! 
