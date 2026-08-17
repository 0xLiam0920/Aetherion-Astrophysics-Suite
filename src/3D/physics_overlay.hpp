#pragma once
// ============================================================
// physics_overlay.hpp
// ============================================================
// World-space line overlays for the 3D viewer. Two of them:
//
//   1. RK4 timelike-geodesic orbits.
//      Integrates the Schwarzschild geodesic for each orbiting body from its
//      (a, e, incl) -> (E, L) and draws the precessing rosette, i.e. the real
//      GR path rather than the Keplerian one the body markers actually follow.
//      (So this and the markers drift apart near periapsis - that's expected.)
//
//   2. Null-geodesic photon rays.
//      A fan of photons at different impact parameters in the equatorial plane;
//      Binet equation (d²u/dφ² + u = 3Mu²) integrated with RK4. Captured rays
//      go red, escaping ones cyan, fading out along the path.
//
// Kept self-contained (no src/2D includes) so it can build into both
// blackhole-3D and, through bh3d_core.hpp, blackhole-sim.
// ============================================================

#include "platform.hpp"
#include "gl_types.hpp"
#include "bh3d_simulationstates.hpp"
#include "orbital_body.hpp"
#include "config.hpp"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/constants.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace bh3d {

// ============================================================
//  Tiny self-contained Schwarzschild + RK4 helpers
//  (mirrors src/2D/2D-physics/{schwarzschild,geodesic,integrator}.hpp
//   but kept local to avoid cross-module include dependencies)
// ============================================================
namespace overlay_phys {

struct Sch {
    double M;
    double horizon()        const { return 2.0 * M; }
    double criticalImpact() const { return 3.0 * std::sqrt(3.0) * M; }
    double f(double r)      const { return 1.0 - 2.0 * M / r; }
    double radialAccel(double r, double L) const {
        double r2 = r*r, r3 = r2*r, r4 = r3*r;
        return -M/r2 + L*L/r3 - 3.0*M*L*L/r4;
    }
    // E, L for bound orbit between periapsis r_p and apoapsis r_a.
    // Falls back to a sane default if the inputs degenerate.
    bool boundOrbitEL(double rp, double ra, double& Eout, double& Lout) const {
        const double rMin = 3.0 * M * 1.01;
        if (rp < rMin) rp = rMin;
        if (ra < rp)   ra = rp * 1.0001;
        double fp = f(rp), fa = f(ra);
        double rp2 = rp*rp, ra2 = ra*ra;
        double denom = fp/rp2 - fa/ra2;
        if (std::abs(denom) < 1e-30) return false;
        double L2 = (fa - fp) / denom;
        if (L2 <= 0.0) return false;
        double E2 = fp * (1.0 + L2/rp2);
        if (E2 <= 0.0) return false;
        Eout = std::sqrt(E2);
        Lout = std::sqrt(L2);
        return true;
    }
};

struct Timelike { double r, phi, vr; };
inline Timelike operator+(const Timelike& a, const Timelike& b) { return {a.r+b.r, a.phi+b.phi, a.vr+b.vr}; }
inline Timelike operator*(const Timelike& a, double s)          { return {a.r*s,   a.phi*s,   a.vr*s  }; }

struct NullState { double u, du; };
inline NullState operator+(const NullState& a, const NullState& b) { return {a.u+b.u, a.du+b.du}; }
inline NullState operator*(const NullState& a, double s)           { return {a.u*s,   a.du*s  }; }

template<typename S, typename F>
inline S rk4(const S& y, double h, F f) {
    S k1 = f(y);
    S k2 = f(y + k1 * (h*0.5));
    S k3 = f(y + k2 * (h*0.5));
    S k4 = f(y + k3 * h);
    return y + (k1 + k2*2.0 + k3*2.0 + k4) * (h/6.0);
}

inline Timelike stepTimelike(const Sch& s, const Timelike& y, double L, double dtau) {
    const double rMin = s.horizon() * 1.01;
    Timelike c = y;
    if (c.r < rMin) c.r = rMin;
    Timelike out = rk4(c, dtau, [&](const Timelike& q) -> Timelike {
        double r = std::max(q.r, rMin);
        return { q.vr, L/(r*r), s.radialAccel(r, L) };
    });
    if (!std::isfinite(out.r) || !std::isfinite(out.phi) || !std::isfinite(out.vr)) return y;
    return out;
}

inline NullState stepNull(const Sch& s, const NullState& y, double dphi) {
    return rk4(y, dphi, [&](const NullState& q) -> NullState {
        return { q.du, -q.u + 3.0 * s.M * q.u * q.u };
    });
}


// ============================================================
//  Kerr characteristic radii (all returned in units of M).
//  Spin `a` here is the dimensionless a* = J/(GM²/c) in [0, 1).
// ============================================================

// Bardeen (1972) marginally-stable circular orbit (ISCO). The prograde and
// retrograde roots collapse onto the single Schwarzschild value 6M when a = 0.
inline void kerrIsco(double a, double& rPro, double& rRetro) {
    const double Z1 = 1.0 + std::cbrt(1.0 - a*a)
                          * (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    const double Z2 = std::sqrt(3.0*a*a + Z1*Z1);
    const double s  = std::sqrt(std::max(0.0, (3.0 - Z1) * (3.0 + Z1 + 2.0*Z2)));
    rPro   = 3.0 + Z2 - s;   // co-rotating photons truncate the disk here
    rRetro = 3.0 + Z2 + s;
}

// Equatorial unstable photon-orbit radii. Both branches give 3M (= 1.5 Rs)
// at a = 0, which is the familiar Schwarzschild photon sphere.
inline void kerrPhotonEquatorial(double a, double& rPro, double& rRetro) {
    rPro   = 2.0 * (1.0 + std::cos((2.0/3.0) * std::acos(-a)));
    rRetro = 2.0 * (1.0 + std::cos((2.0/3.0) * std::acos( a)));
}

} // namespace overlay_phys

// ============================================================
//  Overlay shader sources
// ============================================================
inline const char* kOverlayVS = R"(#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec4 aColor;
uniform mat4 uMVP;
out vec4 vColor;
void main() {
    gl_Position = uMVP * vec4(aPos, 1.0);
    vColor = aColor;
}
)";

inline const char* kOverlayFS = R"(#version 330 core
in vec4 vColor;
out vec4 fragColor;
uniform float uAlpha;
void main() {
    fragColor = vec4(vColor.rgb, vColor.a * uAlpha);
}
)";

// ============================================================
//  PhysicsOverlay: owns GL resources + vertex caches for the two overlays.
// ============================================================
class PhysicsOverlay {
public:
    struct Vertex { glm::vec3 pos; glm::vec4 color; };

    // User-facing toggles + parameters
    bool  orbitsEnabled    = false;
    bool  photonsEnabled   = false;
    bool  spacetimeEnabled = false;
    int   orbitSamples     = 1200;   // vertices per body's RK4 trail
    int   orbitPeriods     = 2;      // how many Keplerian periods to integrate
    int   photonCount      = 40;     // photons per fan
    float photonBMinMul    = 0.45f;  // impact parameter / b_crit, min
    float photonBMaxMul    = 3.00f;  // impact parameter / b_crit, max
    int   spacetimeRings   = 18;     // concentric rings on Flamm's paraboloid
    int   spacetimeSpokes  = 24;     // radial spokes
    int   spacetimeRingSeg = 96;     // segments per ring (smoothness)
    float spacetimeOuterMul= 18.0f;  // outer radius in units of M

    // Kerr characteristic-geometric overlay. Basically, the spin comes from the snapshot's
    // bhSpin each rebuild, so all of these follow the spin slider live.
    bool  iscoEnabled          = false;  // marginally-stable orbit rings
    bool  photonSphereEnabled  = false;  // unstable photon-orbit shell
    bool  ergosphereEnabled    = false;  // static-limit surface (Kerr only)
    bool  shadowContourEnabled = false;  // Bardeen apparent shadow edge
    bool  bfieldEnabled        = false;  // Blandford-Znajek field-line sketch
    int   ringSegments         = 128;    // smoothness for the characteristic rings. 

    bool init() {
        if (inited_) return true;
        if (!prog_.build(kOverlayVS, kOverlayFS)) return false;
        glGenVertexArrays(1, &vao_);
        glGenBuffers(1, &vbo_);
        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                              (void*)offsetof(Vertex, pos));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                              (void*)offsetof(Vertex, color));
        glEnableVertexAttribArray(1);
        glBindVertexArray(0);
        inited_ = true;
        return true;
    }

    ~PhysicsOverlay() {
        if (vbo_) glDeleteBuffers(1, &vbo_);
        if (vao_) glDeleteVertexArrays(1, &vao_);
    }

    // Rebuild the cached vertex buffer from current snapshot + bodies.
    // Cheap enough to call when the user toggles overlays on, switches
    // profile, or resizes; do NOT call every frame.
    void rebuild(const PhysicsSnapshot& snap,
                 const std::vector<OrbitalBody>& bodies)
    {
        if (!inited_) return;
        vertices_.clear();
        strips_.clear();
        spacetimeStart_ = -1;
        spacetimeCount_ = 0;

        if (orbitsEnabled)    buildOrbits(snap, bodies);
        if (photonsEnabled)   buildPhotons(snap);
        if (spacetimeEnabled) {
            spacetimeStart_ = (int)strips_.size();
            buildSpacetime(snap);
            spacetimeCount_ = (int)strips_.size() - spacetimeStart_;
        }

        if (iscoEnabled)          buildISCO(snap);
        if (photonSphereEnabled)  buildPhotonSphere(snap);
        if (ergosphereEnabled)    buildErgosphere(snap);
        if (shadowContourEnabled) buildShadowContour(snap);
        if (bfieldEnabled)        buildBField(snap);

        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(vertices_.size() * sizeof(Vertex)),
                     vertices_.empty() ? nullptr : vertices_.data(),
                     GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        dirty_ = false;
    }

    void markDirty() { dirty_ = true; }
    bool dirty() const { return dirty_; }

    // Marks dirty if any BH/disk scale parameter has changed since the last
    // rebuild. Without this, sliders that resize the BH or disk at runtime
    // would leave the spacetime well and photon fan stuck at the old scale.
    // Spin is folded in too so the Kerr overlays follow the spin slider.
    void notifyScale(float bhRadius, float diskOuterRadius, float bhSpin = -1.0f) {
        bool changed = false;
        if (std::abs(bhRadius - lastBHRadius_) > 1e-5f)         { lastBHRadius_  = bhRadius;        changed = true; }
        if (std::abs(diskOuterRadius - lastDiskOuter_) > 1e-4f) { lastDiskOuter_ = diskOuterRadius; changed = true; }
        if (bhSpin >= 0.0f && std::abs(bhSpin - lastBHSpin_) > 1e-5f) { lastBHSpin_ = bhSpin;       changed = true; }
        if (changed) dirty_ = true;
    }

    // Draw the cached strips. Builds view/proj from camera snapshot.
    void draw(const PhysicsSnapshot& snap) {
        if (!inited_ || strips_.empty()) return;

        const float aspect = (snap.windowH > 0)
            ? (float)snap.windowW / (float)snap.windowH : 1.0f;
        const float fovRad = glm::radians(snap.fov);
        // Far plane chosen large enough to contain integrated orbits and
        // photon trails even for the biggest presets.
        glm::mat4 proj = glm::perspective(fovRad, aspect, 0.1f, 100000.0f);
        glm::vec3 center = snap.cameraPos + snap.cameraDir;
        glm::mat4 view = glm::lookAt(snap.cameraPos, center, snap.cameraUp);
        glm::mat4 mvp  = proj * view;

        // Lines render on top of the screen-space ray-marched scene; no
        // depth buffer is shared so we just disable depth test/write and
        // rely on additive-ish alpha blending for the soft glow look.
        GLboolean prevDepth = glIsEnabled(GL_DEPTH_TEST);
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);   // additive
        glLineWidth(1.5f);

        prog_.use();
        glUniformMatrix4fv(prog_.uniform("uMVP"), 1, GL_FALSE, &mvp[0][0]);
        prog_.set1f("uAlpha", 1.0f);

        glBindVertexArray(vao_);

        // Pass 1: orbits + photons with additive blend + thin lines (glow).
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);
        glLineWidth(1.5f);
        const int spaceEnd = (spacetimeStart_ >= 0)
            ? spacetimeStart_ + spacetimeCount_ : (int)strips_.size();
        for (int i = 0; i < (int)strips_.size(); ++i) {
            if (spacetimeStart_ >= 0 && i >= spacetimeStart_ && i < spaceEnd) continue;
            glDrawArrays(GL_LINE_STRIP, strips_[i].first, strips_[i].count);
        }

        // Pass 2: spacetime grid with standard alpha + thicker lines so it
        // reads against both the bright disk and dark sky.
        if (spacetimeStart_ >= 0 && spacetimeCount_ > 0) {
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glLineWidth(2.5f);
            for (int i = spacetimeStart_; i < spaceEnd; ++i) {
                glDrawArrays(GL_LINE_STRIP, strips_[i].first, strips_[i].count);
            }
        }

        glBindVertexArray(0);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        if (prevDepth) glEnable(GL_DEPTH_TEST);
    }

private:
    struct Strip { GLint first; GLsizei count; };

    void appendStrip(const std::vector<Vertex>& v) {
        if (v.size() < 2) return;
        Strip s;
        s.first = (GLint)vertices_.size();
        s.count = (GLsizei)v.size();
        vertices_.insert(vertices_.end(), v.begin(), v.end());
        strips_.push_back(s);
    }

    void buildOrbits(const PhysicsSnapshot& snap,
                     const std::vector<OrbitalBody>& bodies)
    {
        if (bodies.empty() || snap.bhRadius <= 0.0f) return;

        overlay_phys::Sch sch{ (double)snap.bhRadius * 0.5 };  // M = Rs / 2
        const glm::vec3 origin = snap.bhPosition;

        for (size_t bi = 0; bi < bodies.size(); ++bi) {
            const auto& bc = bodies[bi].config();
            // Convert "units of Rs" → world units. semiMajor is in Rs and the
            // OrbitalBody analytical path uses the raw float as if it were in
            // world units, so the GR overlay matches by treating them the
            // same way. (We just need a self-consistent M to compute E, L.)
            const double a   = (double)bc.semiMajor;
            const double e   = std::clamp((double)bc.eccentricity, 0.0, 0.97);
            const double rp  = a * (1.0 - e);
            const double ra  = a * (1.0 + e);
            const double incl = (double)bc.inclination;

            double E, L;
            if (!sch.boundOrbitEL(rp, ra, E, L)) continue;

            // Integrate from periapsis with vr=0; pick a proper-time step
            // sized to the analytical period so orbitPeriods controls how
            // many laps we draw without runaway sample counts.
            const double period = 2.0 * glm::pi<double>() * std::sqrt(a*a*a);
            const int    nSteps = std::max(64, orbitSamples);
            const double tEnd   = period * std::max(1, orbitPeriods);
            const double dtau   = tEnd / (double)nSteps;

            overlay_phys::Timelike st{ rp, 0.0, 0.0 };
            std::vector<Vertex> verts;
            verts.reserve(nSteps + 1);

            const glm::vec4 baseColor = glm::vec4(
                bc.bodyColor.r * 0.9f + 0.15f,
                bc.bodyColor.g * 0.9f + 0.15f,
                bc.bodyColor.b * 0.9f + 0.40f,
                0.85f);

            for (int i = 0; i <= nSteps; ++i) {
                double r = std::max(st.r, sch.horizon() * 1.01);
                // Position in orbital plane, then rotate by inclination
                // (around X axis) to match OrbitalBody's transform.
                double xo = r * std::cos(st.phi);
                double zo = r * std::sin(st.phi);
                glm::vec3 local((float)xo,
                                (float)(zo * std::sin(incl)),
                                (float)(zo * std::cos(incl)));
                Vertex v;
                v.pos   = origin + local;
                // Subtle pulse along the trail so motion direction reads.
                float t = (float)i / (float)nSteps;
                float gloss = 0.55f + 0.45f * std::sin(t * 24.0f);
                v.color = glm::vec4(baseColor.r * gloss,
                                    baseColor.g * gloss,
                                    baseColor.b * gloss,
                                    baseColor.a);
                verts.push_back(v);
                st = overlay_phys::stepTimelike(sch, st, L, dtau);
                if (st.r <= sch.horizon() * 1.05) break; // plunged in
                if (st.r > a * 10.0) break;              // escaped, safety
            }
            appendStrip(verts);
        }
    }

    void buildPhotons(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f || photonCount < 2) return;

        overlay_phys::Sch sch{ (double)snap.bhRadius * 0.5 };
        const glm::vec3 origin = snap.bhPosition;
        const double bCrit = sch.criticalImpact();

        // Equatorial plane fan: photons in the (x, z) plane around the BH.
        // Each photon comes from "infinity" (φ = π) at impact parameter b.
        const int N = photonCount;
        for (int i = 0; i < N; ++i) {
            const double t = (N == 1) ? 0.5 : (double)i / (double)(N - 1);
            const double bMul = (double)photonBMinMul
                              + t * ((double)photonBMaxMul - (double)photonBMinMul);
            const double b = bMul * bCrit;
            const bool captured = (bMul < 1.0);

            std::vector<Vertex> verts;
            verts.reserve(800);

            // Start far out along the negative x axis, integrate inward
            // toward periapsis and back out (or into the horizon).
            //
            // For escape (b > b_crit) use periapsis = positive real root of
            // 1/b² = u²(1 − 2Mu). Approximate by Newton iteration; for
            // capture (b ≤ b_crit) start at the photon sphere u₀ = 1/(3M)
            // and integrate inward.
            double u0;
            if (captured) {
                u0 = 1.0 / (3.0 * sch.M);                 // photon sphere
            } else {
                // Newton on g(u) = 1/b² − u²(1 − 2Mu) starting at u = 1/(3M·b/bCrit)
                u0 = 1.0 / (3.5 * sch.M);
                for (int k = 0; k < 24; ++k) {
                    double g  = 1.0/(b*b) - u0*u0 * (1.0 - 2.0*sch.M*u0);
                    double gp = -2.0*u0 + 6.0*sch.M*u0*u0;
                    if (std::abs(gp) < 1e-30) break;
                    u0 -= g / gp;
                    if (u0 <= 0.0) { u0 = 1e-6; break; }
                }
            }

            overlay_phys::NullState ns{ u0, 0.0 };
            const double horizonU = 1.0 / (sch.horizon() * 1.001);
            const double dphi     = 0.01;
            const int    maxSteps = 1800;

            // Trace one branch in polar coords (r, phi) starting at periapsis.
            // We keep polar samples here and rotate the whole trajectory at
            // the end so all incoming rays share a common asymptotic
            // direction (φ_in = π, i.e. coming from the left side of the
            // BH) — that's the textbook "parallel beam deflection" picture.
            struct PolarSample { double r, phi; int k; };
            auto traceHalfPolar = [&](double phiStep) {
                overlay_phys::NullState s = ns;
                double phi = 0.0;
                std::vector<PolarSample> out;
                out.reserve(maxSteps);
                for (int k = 0; k < maxSteps; ++k) {
                    if (!std::isfinite(s.u) || s.u <= 0.0) break;
                    double r = 1.0 / s.u;
                    if (r > 250.0 * sch.M) break;
                    if (s.u >= horizonU) break;
                    out.push_back({ r, phi, k });
                    s = overlay_phys::stepNull(sch, s, phiStep);
                    phi += phiStep;
                }
                return out;
            };

            std::vector<PolarSample> outHalf = traceHalfPolar(+dphi);
            std::vector<PolarSample> inHalf;
            if (!captured) inHalf = traceHalfPolar(-dphi);

            // Rotation: line up the most-distant incoming sample (escape)
            // or the most-distant single-branch sample (capture) with the
            // world-space direction φ = π. Result: every ray's "source at
            // infinity" lies on the same side of the BH.
            double rot = 0.0;
            if (captured) {
                // Captured rays go inward; align their starting point with π.
                if (!outHalf.empty()) {
                    rot = glm::pi<double>() - outHalf.front().phi;
                }
            } else {
                // The incoming branch's last sample is the farthest point on
                // that side — treat that as "from infinity".
                if (!inHalf.empty()) {
                    rot = glm::pi<double>() - inHalf.back().phi;
                }
            }

            auto pushPoint = [&](const PolarSample& p, std::vector<Vertex>& dst) {
                double phi = p.phi + rot;
                glm::vec3 local((float)(p.r * std::cos(phi)),
                                0.0f,
                                (float)(p.r * std::sin(phi)));
                Vertex v;
                v.pos = origin + local;
                float a = std::max(0.15f, 1.0f - (float)p.k / (float)maxSteps);
                if (captured) v.color = glm::vec4(1.00f, 0.35f, 0.20f, 0.7f * a);
                else          v.color = glm::vec4(0.45f, 0.85f, 1.00f, 0.7f * a);
                dst.push_back(v);
            };

            std::vector<Vertex> full;
            full.reserve(outHalf.size() + inHalf.size());
            if (captured) {
                for (const auto& p : outHalf) pushPoint(p, full);
            } else {
                // Walk incoming from infinity (its last sample) back to
                // periapsis, then outgoing away to infinity.
                for (auto it = inHalf.rbegin(); it != inHalf.rend(); ++it) pushPoint(*it, full);
                for (const auto& p : outHalf) pushPoint(p, full);
            }
            appendStrip(full);
        }
    }

    // Flamm's paraboloid: the embedding diagram of an equatorial slice of
    // Schwarzschild spacetime into flat 3D. z(r) = -2 sqrt(2M (r - 2M)) for
    // r >= 2M. Rendered as a ring + spoke grid so the well's geometry reads
    // clearly under the equatorial plane.
    void buildSpacetime(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;

        const double M     = (double)snap.bhRadius * 0.5; // Rs = 2M
        const double r2M   = 2.0 * M;
        const double rIn   = r2M * 1.01;                  // start just outside horizon

        // Outer extent: extend well beyond the disk to convey the BH's
        // long-range gravitational reach (sphere of influence). Use ~5x
        // the visible disk radius, floored to 30 Schwarzschild radii so
        // it never collapses to a tiny puddle when the disk is small.
        double rOut = (snap.diskOuterRadius > 0.0f)
            ? (double)snap.diskOuterRadius * 5.0
            : std::max((double)spacetimeOuterMul, 6.0) * M;
        rOut = std::max(rOut, r2M * 30.0);
        if (rOut <= rIn) return;

        // Vertical centering: the Flamm embedding has y=0 at the horizon
        // and rises to y_max at the outer rim. Without a shift the whole
        // funnel sits above the BH; offset by -y_max/2 so the well is
        // visually centered on the BH at the equatorial plane.
        const double yMax  = 2.0 * std::sqrt(2.0 * M * (rOut - r2M));
        const double yBias = -0.5 * yMax;

        const glm::vec3 origin = snap.bhPosition;

        auto embed = [&](double r) -> double {
            // Flamm paraboloid, shifted so the well is centered vertically
            // on the BH instead of sitting entirely above the equator.
            return 2.0 * std::sqrt(std::max(0.0, 2.0 * M * (r - r2M))) + yBias;
        };

        auto colorFor = [&](double r) {
            // Warm cream — slight brightness ramp from throat to rim so the
            // funnel still reads as 3D even with a single hue.
            double t = (r - rIn) / (rOut - rIn);
            t = std::clamp(t, 0.0, 1.0);
            float b = 0.85f + 0.15f * (float)t;
            return glm::vec4(1.00f * b, 0.96f * b, 0.82f * b, 1.0f);
        };

        // ── Rings (concentric, equatorial) ──────────────────────────────
        const int   nRings = std::max(4, spacetimeRings);
        const int   nSeg   = std::max(16, spacetimeRingSeg);
        for (int i = 0; i < nRings; ++i) {
            // Log-spaced radii so the well is sampled finely near the horizon.
            double t = (double)i / (double)(nRings - 1);
            double r = rIn * std::pow(rOut / rIn, t);
            double y = embed(r);
            glm::vec4 col = colorFor(r);

            std::vector<Vertex> verts;
            verts.reserve(nSeg + 1);
            for (int s = 0; s <= nSeg; ++s) {
                double phi = (2.0 * glm::pi<double>()) * (double)s / (double)nSeg;
                glm::vec3 local((float)(r * std::cos(phi)),
                                (float)y,
                                (float)(r * std::sin(phi)));
                verts.push_back({ origin + local, col });
            }
            appendStrip(verts);
        }

        // ── Spokes (radial profile curves) ──────────────────────────────
        const int nSpokes  = std::max(4, spacetimeSpokes);
        const int nSpokeSm = std::max(32, spacetimeRingSeg / 2);
        for (int j = 0; j < nSpokes; ++j) {
            double phi = (2.0 * glm::pi<double>()) * (double)j / (double)nSpokes;
            double cphi = std::cos(phi), sphi = std::sin(phi);

            std::vector<Vertex> verts;
            verts.reserve(nSpokeSm + 1);
            for (int k = 0; k <= nSpokeSm; ++k) {
                double t = (double)k / (double)nSpokeSm;
                double r = rIn * std::pow(rOut / rIn, t);
                double y = embed(r);
                glm::vec3 local((float)(r * cphi),
                                (float)y,
                                (float)(r * sphi));
                verts.push_back({ origin + local, colorFor(r) });
            }
            appendStrip(verts);
        }
    }

    // ============================================================
    //  Kerr characteristic-geometry builders (items 14-18)
    // ============================================================

    // Shared helper: a full circle in the equatorial (XZ) plane. The disk in
    // this engine lies in XZ with the spin axis along +Y, so every ISCO /
    // photon ring is just a horizontal circle around the hole.
    void appendEquatorialRing(const glm::vec3& origin, double r, const glm::vec4& col) {
        if (r <= 0.0) return;
        const int nSeg = std::max(24, ringSegments);
        std::vector<Vertex> verts;
        verts.reserve(nSeg + 1);
        for (int s = 0; s <= nSeg; ++s) {
            double phi = (2.0 * glm::pi<double>()) * (double)s / (double)nSeg;
            glm::vec3 local((float)(r * std::cos(phi)), 0.0f, (float)(r * std::sin(phi)));
            verts.push_back({ origin + local, col });
        }
        appendStrip(verts);
    }

    // Shared helper: a lat/long wireframe sphere around the +Y spin axis.
    void appendWireSphere(const glm::vec3& origin, double r, const glm::vec4& col) {
        if (r <= 0.0) return;
        const int nLat = 9;
        const int nSeg = std::max(32, ringSegments / 2);
        for (int i = 1; i < nLat; ++i) {                 // latitude rings
            double theta = glm::pi<double>() * (double)i / (double)nLat;
            double y   = r * std::cos(theta);
            double rho = r * std::sin(theta);
            std::vector<Vertex> verts;
            verts.reserve(nSeg + 1);
            for (int s = 0; s <= nSeg; ++s) {
                double phi = 2.0 * glm::pi<double>() * (double)s / (double)nSeg;
                verts.push_back({ origin + glm::vec3((float)(rho * std::cos(phi)),
                                                     (float)y,
                                                     (float)(rho * std::sin(phi))), col });
            }
            appendStrip(verts);
        }
        const int nMer = 8;                              // meridian arcs (pole to pole)
        const int nArc = std::max(24, nSeg / 2);
        for (int j = 0; j < nMer; ++j) {
            double phi = glm::pi<double>() * (double)j / (double)nMer;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            std::vector<Vertex> verts;
            verts.reserve(nArc + 1);
            for (int k = 0; k <= nArc; ++k) {
                double theta = glm::pi<double>() * (double)k / (double)nArc;
                double y   = r * std::cos(theta);
                double rho = r * std::sin(theta);
                verts.push_back({ origin + glm::vec3((float)(rho * cphi),
                                                     (float)y,
                                                     (float)(rho * sphi)), col });
            }
            appendStrip(verts);
        }
    }

    // ISCO rings. The prograde radius (gold) is the physically relevant one —
    // that's where a co-rotating accretion disk actually ends. The retrograde
    // radius (dim red) is drawn alongside so the spin dependence is obvious:
    // crank a* up and watch the gold ring shrink toward the horizon while the
    // red one grows. At a* = 0 they sit on top of each other at 3 Rs.
    void buildISCO(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;
        const double M = (double)snap.bhRadius * 0.5;
        const double a = std::clamp((double)snap.bhSpin, 0.0, 0.9999);
        double rPro, rRetro;
        overlay_phys::kerrIsco(a, rPro, rRetro);
        appendEquatorialRing(snap.bhPosition, rPro   * M, glm::vec4(1.00f, 0.82f, 0.28f, 0.90f));
        appendEquatorialRing(snap.bhPosition, rRetro * M, glm::vec4(0.95f, 0.28f, 0.22f, 0.50f));
    }

    // Photon-orbit shell. Kerr's photon region isn't actually a sphere (the
    // orbit radius depends on inclination), so this is a readable stand-in:
    // a single wireframe shell at the prograde equatorial photon radius, which
    // is exactly the 1.5 Rs Schwarzschild photon sphere when a* = 0.
    void buildPhotonSphere(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;
        const double M = (double)snap.bhRadius * 0.5;
        const double a = std::clamp((double)snap.bhSpin, 0.0, 0.9999);
        double rPro, rRetro;
        overlay_phys::kerrPhotonEquatorial(a, rPro, rRetro);
        appendWireSphere(snap.bhPosition, rPro * M, glm::vec4(0.65f, 0.80f, 1.00f, 0.30f));
    }

    // Ergosphere / static-limit surface. r_ergo(θ)/M = 1 + sqrt(1 - a²cos²θ),
    // with the polar angle θ measured from the +Y spin axis. It touches the
    // horizon at the poles and bulges out to 2M (= Rs) at the equator, so it
    // only exists for a spinning hole; a* = 0 gets skipped.
    void buildErgosphere(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;
        const double a = std::clamp((double)snap.bhSpin, 0.0, 0.9999);
        if (a < 1e-3) return;
        const double M = (double)snap.bhRadius * 0.5;
        const glm::vec3 origin = snap.bhPosition;
        const glm::vec4 col(0.78f, 0.45f, 0.98f, 0.32f);

        auto rErgo = [&](double theta) {
            double ct = std::cos(theta);
            return (1.0 + std::sqrt(std::max(0.0, 1.0 - a*a*ct*ct))) * M;
        };

        const int nLat = 14;
        const int nSeg = std::max(32, ringSegments / 2);
        for (int i = 1; i < nLat; ++i) {                 // latitude rings
            double theta = glm::pi<double>() * (double)i / (double)nLat;
            double rE  = rErgo(theta);
            double y   = rE * std::cos(theta);
            double rho = rE * std::sin(theta);
            std::vector<Vertex> verts;
            verts.reserve(nSeg + 1);
            for (int s = 0; s <= nSeg; ++s) {
                double phi = 2.0 * glm::pi<double>() * (double)s / (double)nSeg;
                verts.push_back({ origin + glm::vec3((float)(rho * std::cos(phi)),
                                                     (float)y,
                                                     (float)(rho * std::sin(phi))), col });
            }
            appendStrip(verts);
        }
        const int nMer = 12;                             // meridian arcs
        const int nArc = std::max(24, ringSegments / 4);
        for (int j = 0; j < nMer; ++j) {
            double phi = 2.0 * glm::pi<double>() * (double)j / (double)nMer;
            double cphi = std::cos(phi), sphi = std::sin(phi);
            std::vector<Vertex> verts;
            verts.reserve(nArc + 1);
            for (int k = 0; k <= nArc; ++k) {
                double theta = glm::pi<double>() * (double)k / (double)nArc;
                double rE  = rErgo(theta);
                double y   = rE * std::cos(theta);
                double rho = rE * std::sin(theta);
                verts.push_back({ origin + glm::vec3((float)(rho * cphi),
                                                     (float)y,
                                                     (float)(rho * sphi)), col });
            }
            appendStrip(verts);
        }
    }

    // Apparent shadow edge, Bardeen (1972) celestial coordinates (α, β). The
    // shadow is what the EHT actually images, and for a* > 0 it's flattened on
    // one side because photons co-rotating with the hole fall in more easily.
    //
    // We compute the critical curve in the (α, β) image plane, then billboard
    // it in front of the hole: β runs along the projected spin axis, α runs
    // perpendicular to it. That keeps the flattened edge pointing the right way
    // no matter where the camera is. It does mean the curve has to rebuild as
    // the camera moves, which the draw gate handles by marking it dirty.
    void buildShadowContour(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;
        const double M = (double)snap.bhRadius * 0.5;
        const double a = std::clamp((double)snap.bhSpin, 0.0, 0.9999);

        glm::vec3 toBH = snap.bhPosition - snap.cameraPos;
        double len = (double)glm::length(toBH);
        if (len < 1e-4) return;
        glm::vec3 viewDir = toBH / (float)len;           // camera -> hole
        glm::vec3 obsDir  = -viewDir;                    // hole -> observer

        // Inclination is the angle between the spin axis (+Y) and the observer.
        double cosI = std::clamp((double)obsDir.y, -1.0, 1.0);
        double sinI = std::sqrt(std::max(1e-6, 1.0 - cosI * cosI));

        // Billboard basis: "up" is the spin axis projected onto the image plane
        // (so β aligns with it), "right" is perpendicular for α.
        glm::vec3 yAxis(0.0f, 1.0f, 0.0f);
        glm::vec3 up = yAxis - viewDir * glm::dot(yAxis, viewDir);
        if (glm::length(up) < 1e-4f) up = snap.cameraUp; // near pole-on, pick any stable up
        up = glm::normalize(up);
        glm::vec3 right = glm::normalize(glm::cross(viewDir, up));

        const glm::vec4 col(1.00f, 0.95f, 0.72f, 0.90f);
        std::vector<Vertex> verts;

        auto worldOf = [&](double alpha, double beta) {
            // α, β are in units of M; place the point on the billboard plane.
            return snap.bhPosition + (float)M * ((float)alpha * right + (float)beta * up);
        };

        if (a < 1e-3) {
            // Schwarzschild: a perfect circle of radius 3√3 M (= √27).
            const double rSh = std::sqrt(27.0);
            const int nSeg = std::max(64, ringSegments);
            verts.reserve(nSeg + 1);
            for (int s = 0; s <= nSeg; ++s) {
                double ang = 2.0 * glm::pi<double>() * (double)s / (double)nSeg;
                verts.push_back({ worldOf(rSh * std::cos(ang), rSh * std::sin(ang)), col });
            }
            appendStrip(verts);
            return;
        }

        // Conserved quantities of the spherical photon orbits, parametrised by
        // Boyer-Lindquist radius r (units M). These are the standard Bardeen
        // expressions; α = -λ/sinθ, β² = η + a²cos²θ - λ²cot²θ.
        auto lambdaEta = [&](double r, double& lam, double& eta) {
            double rm1 = r - 1.0;
            lam = (r*r*(3.0 - r) - a*a*(r + 1.0)) / (a * rm1);
            eta = (r*r*r*(4.0*a*a - r*(r - 3.0)*(r - 3.0))) / (a*a * rm1*rm1);
        };

        double r1, r2;
        overlay_phys::kerrPhotonEquatorial(a, r1, r2);   // prograde .. retrograde
        const double cot2 = (cosI * cosI) / (sinI * sinI);
        const int nR = 160;
        verts.reserve(2 * nR + 2);

        auto addPoint = [&](double r, double sign) {
            double lam, eta;
            lambdaEta(r, lam, eta);
            double beta2 = eta + a*a*cosI*cosI - lam*lam*cot2;
            if (beta2 < 0.0) beta2 = 0.0;               // clamp the tiny endpoints
            double alpha = -lam / sinI;
            double beta  = sign * std::sqrt(beta2);
            verts.push_back({ worldOf(alpha, beta), col });
        };

        for (int i = 0; i <= nR; ++i) { double r = r1 + (r2 - r1) * (double)i / (double)nR; addPoint(r, +1.0); }
        for (int i = nR; i >= 0; --i) { double r = r1 + (r2 - r1) * (double)i / (double)nR; addPoint(r, -1.0); }
        appendStrip(verts);
    }

    // Magnetic-field overlay. It
    // gets drawn depends on what we actually know about the system, because a
    // real field topology is almost never possible to observe, at least with current technology. it also doesn't
    //  help that we essentially have a great wall of china for interstellar gas and dust.
    //   1. Jet-launching hole  -> Blandford-Znajek funnel (needs spin + a real
    //      large-scale field to tap, so we only draw it when the sim is actively
    //      showing jets and the hole spins).
    //   2. Wind-fed binary     -> a weak, unresolved stellar-wind field streaming
    //      off the companion and being gravitationally focused past the hole. This
    //      is the Gaia BH1 case: hypothetical, faint, and it tracks the companion's
    //      live orbital position rather than being baked to a preset.
    //   3. Isolated + quiescent -> we honestly don't know the field, so only a
    //      barely-there dipole cage is hinted.
    //
    // Everything scales off snapshot properties (spin, Rs, companion position and
    // radius, jet reach, activity) - no per-preset hand tuning.
    void buildBField(const PhysicsSnapshot& snap) {
        if (snap.bhRadius <= 0.0f) return;
        const double a = std::clamp((double)snap.bhSpin, 0.0, 0.9999);
        // BZ needs both spin and an accretion-powered field; "jets on" is our
        // proxy for an active, jet-launching state.
        const bool jetActive = snap.jetsEnabled && a > 0.08;
        const int  compIdx   = findWindCompanion(snap);

        if (jetActive)        buildBZFunnel(snap, a);
        if (compIdx >= 0)     buildWindField(snap, compIdx, jetActive);
        if (!jetActive && compIdx < 0) buildQuiescentHint(snap);
    }

    // Nearest wind-bearing companion (a real star, not a compact remnant).
    // Returns the value '-1' when the hole is isolated or the companion is a compact object (i.e, like a neutron star or secondary black hole, as seen somewhere like OJ 287).
    int findWindCompanion(const PhysicsSnapshot& snap) const {
        if (!snap.orbBodyEnabled) return -1;
        int best = -1;
        double bestD = 1e30;
        const glm::vec3 B = snap.bhPosition;
        for (size_t i = 0; i < snap.orbBodyPositions.size(); ++i) {
            const int t = (i < snap.orbBodyTypes.size()) ? snap.orbBodyTypes[i] : 0;
            if (t != 0 && t != 6) continue;              // Star or CompanionStar only
            const double d = (double)glm::length(snap.orbBodyPositions[i] - B);
            if (d < bestD) { bestD = d; best = (int)i; }
        }
        return best;
    }

    // This is essentially just a topology solution.
    // solution. Higher spin when both winds the field up harder (twist) and collimates
    // it more tightly to induce some form of flare, so the shape actually responds to a*.
    void buildBZFunnel(const PhysicsSnapshot& snap, double a) {
        const double Rs = (double)snap.bhRadius;
        const glm::vec3 origin = snap.bhPosition;
        const glm::vec4 baseCol(0.10f, 0.20f, 0.85f, 0.75f);
        const glm::vec4 tipCol (0.25f, 0.90f, 1.45f, 0.60f);

        const double jetLen = std::max((double)snap.jetLength, Rs * 6.0);
        const double footR  = 1.15 * Rs;                 // footpoint just outside the horizon
        const double twist  = a * 6.0;                   // frame-drag winding grows with spin
        const double flare  = 1.30 - 0.60 * a;           // high spin = tighter funnel
        const int nLines = 12;
        const int nSeg   = 60;

        for (int L = 0; L < nLines; ++L) {
            double phi0 = 2.0 * glm::pi<double>() * (double)L / (double)nLines;
            for (int hemi = 0; hemi < 2; ++hemi) {
                double dir = (hemi == 0) ? 1.0 : -1.0;   // north / south funnel
                std::vector<Vertex> verts;
                verts.reserve(nSeg + 1);
                for (int k = 0; k <= nSeg; ++k) {
                    double t = (double)k / (double)nSeg;
                    double h   = dir * jetLen * t;
                    double rho = footR + flare * (0.6 * Rs * std::sqrt(t) + 0.15 * jetLen * t);
                    double phi = phi0 + twist * t;
                    glm::vec3 p(origin.x + (float)(rho * std::cos(phi)),
                                origin.y + (float)h,
                                origin.z + (float)(rho * std::sin(phi)));
                    verts.push_back({ p, glm::mix(baseCol, tipCol, (float)t) });
                }
                appendStrip(verts);
            }
        }
    }

    // Weak stellar-wind field from the companion, gravitationally
    // focused as it streams past the hole. Kept faint and unassertive on purpose
    // this is a model of an unresolved wind, and for a
    // quiescent (non-accreting) system the bending stays gentle.
    void buildWindField(const PhysicsSnapshot& snap, int idx, bool accreting) {
        const double Rs = (double)snap.bhRadius;
        const glm::vec3 B = snap.bhPosition;
        const glm::vec3 C = snap.orbBodyPositions[idx];
        const double Rc = (idx < (int)snap.orbBodyRadii.size())
                        ? std::max((double)snap.orbBodyRadii[idx], Rs * 0.4) : Rs;

        glm::vec3 CB = B - C;
        double sep = (double)glm::length(CB);
        if (sep < 1e-4) return;
        glm::vec3 axis = CB / (float)sep;                // companion -> hole
        glm::vec3 up   = (std::abs(axis.y) < 0.9f) ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0);
        glm::vec3 e1   = glm::normalize(glm::cross(up, axis));
        glm::vec3 e2   = glm::cross(axis, e1);


        const float aMax = accreting ? 0.42f : 0.26f;
        const glm::vec4 nearCol(0.55f, 0.70f, 0.95f, aMax);
        const glm::vec4 farCol (0.30f, 0.55f, 0.90f, aMax * 0.45f);

        const double ds      = std::max(Rs * 0.35, sep * 0.012);
        const int    maxSteps = 240;
        const double focusK  = accreting ? 0.60 : 0.32; // how hard the wind bends toward the hole
        const int    nLines  = 16;

        for (int L = 0; L < nLines; ++L) {
            double ang    = 2.0 * glm::pi<double>() * (double)L / (double)nLines;
            double spread = 0.15 + 0.75 * (double)((L % 4)) / 3.0;   // unfold the launch cone out
            glm::vec3 dirLaunch = glm::normalize(
                axis + (float)(spread * std::cos(ang)) * e1
                     + (float)(spread * std::sin(ang)) * e2);
            glm::vec3 p = C + dirLaunch * (float)Rc;     // start at the stellar surface
            glm::vec3 v = dirLaunch;                     // wind leaves the star roughly radial

            std::vector<Vertex> verts;
            verts.reserve(maxSteps);
            for (int k = 0; k < maxSteps; ++k) {
                double rB = (double)glm::length(B - p);
                if (rB < 1.4 * Rs) break;                // captured near the hole
                if ((double)glm::length(p - C) > sep * 1.6) break;   // swept past / escaped
                float frac = (float)std::min(1.0, (double)glm::length(p - C) / sep);
                verts.push_back({ p, glm::mix(nearCol, farCol, frac) });

                glm::vec3 toB    = glm::normalize(B - p);
                glm::vec3 radial = glm::normalize(p - C);
                double bend = focusK * std::min(1.0, (Rs * Rs) / (rB * rB)); // ~ (Rs/r)^2, capped
                v = glm::normalize(v + (float)bend * toB + 0.04f * radial);
                p += v * (float)ds;
            }
            if (verts.size() >= 2) appendStrip(verts);
        }

        buildCoronalLoops(C, Rc, axis, e1, e2, aMax);
    }

    // Small closed dipole loops anchored on the companion, echoing the classic
    // solar-type coronal field. Purely illustrative of the wind's stellar origin.
    void buildCoronalLoops(const glm::vec3& C, double Rc, const glm::vec3& axis,
                           const glm::vec3& e1, const glm::vec3& e2, float aMax) {
        const glm::vec4 col(0.65f, 0.75f, 1.0f, aMax * 0.8f);
        const int nLoops = 6, nArc = 24;
        glm::vec3 anti = -axis;                          // anchor on the near/side face
        for (int l = 0; l < nLoops; ++l) {
            double ang = 2.0 * glm::pi<double>() * (double)l / (double)nLoops;
            glm::vec3 base = glm::normalize(anti * 0.5f
                             + (float)std::cos(ang) * e1
                             + (float)std::sin(ang) * e2);
            glm::vec3 perp = glm::normalize(glm::cross(base, axis));
            const double sepAng = 0.35;                  // angular footpoint separation
            glm::vec3 f1 = glm::normalize(base * (float)std::cos(sepAng) + perp * (float)std::sin(sepAng));
            glm::vec3 f2 = glm::normalize(base * (float)std::cos(sepAng) - perp * (float)std::sin(sepAng));
            std::vector<Vertex> verts;
            verts.reserve(nArc + 1);
            for (int k = 0; k <= nArc; ++k) {
                double u = (double)k / (double)nArc;
                glm::vec3 dirp = glm::normalize(glm::mix(f1, f2, (float)u));
                double bulge = 1.0 + 0.6 * std::sin(glm::pi<double>() * u); // loop rises above the surface
                verts.push_back({ C + dirp * (float)(Rc * bulge), col });
            }
            appendStrip(verts);
        }
    }

    // We don't know the field here, so this
    // is intentionally almost invisible. Call it pseudo science, but a faint dipole cage just so the toggle
    // shows the hole isn't claimed to be field-free.
    void buildQuiescentHint(const PhysicsSnapshot& snap) {
        const double Rs = (double)snap.bhRadius;
        const glm::vec3 O = snap.bhPosition;
        const glm::vec4 col(0.35f, 0.45f, 0.78f, 0.15f);
        const int nLines = 8, nSeg = 60;
        const double shells[3] = { 4.0 * Rs, 7.0 * Rs, 11.0 * Rs };
        for (double Lmax : shells) {
            for (int L = 0; L < nLines; ++L) {
                double phi = 2.0 * glm::pi<double>() * (double)L / (double)nLines;
                double cphi = std::cos(phi), sphi = std::sin(phi);
                std::vector<Vertex> verts;
                verts.reserve(nSeg + 1);
                for (int k = 0; k <= nSeg; ++k) {
                    // Dipole field line r = Lmax sin^2(theta), poles trimmed off.
                    double th = glm::pi<double>() * (0.04 + 0.92 * (double)k / (double)nSeg);
                    double r  = Lmax * std::sin(th) * std::sin(th);
                    double y  = r * std::cos(th);
                    double rho = r * std::sin(th);
                    verts.push_back({ O + glm::vec3((float)(rho * cphi), (float)y, (float)(rho * sphi)), col });
                }
                appendStrip(verts);
            }
        }
    }

    GLProgram prog_;
    GLuint    vao_ = 0;
    GLuint    vbo_ = 0;
    std::vector<Vertex> vertices_;
    std::vector<Strip>  strips_;
    int  spacetimeStart_ = -1;
    int  spacetimeCount_ = 0;
    bool inited_ = false;
    bool dirty_  = true;
    float lastBHRadius_  = -1.0f;
    float lastDiskOuter_ = -1.0f;
    float lastBHSpin_    = -1.0f;
};

} // namespace bh3d
