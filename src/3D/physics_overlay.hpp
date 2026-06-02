#pragma once
// ============================================================
// physics_overlay.hpp
// ============================================================
// World-space line overlays for the 3D black-hole viewer:
//
//   1. RK4 timelike-geodesic orbits
//      Numerically integrates the Schwarzschild geodesic equation for each
//      orbital body using its semi-major axis / eccentricity / inclination to
//      derive (E, L). Draws the resulting rosette (apsidal-precession) curve
//      that shows the *real* GR orbit, not the analytical Keplerian path the
//      body markers follow.
//
//   2. Null-geodesic photon rays
//      Casts a fan of photons at varying impact parameters in the equatorial
//      plane, integrates the Binet equation (d²u/dφ² + u = 3Mu²) with RK4,
//      colors captured rays red and escaping rays cyan with brightness
//      falloff along the path.
//
// Self-contained: no dependency on src/2D headers. Compiles into blackhole-3D
// and (via bh3d_core.hpp) into blackhole-sim.
// ============================================================

#include "platform.hpp"
#include "gl_types.hpp"
#include "simulation_state.hpp"
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

// ────────────────────────────────────────────────────────────
//  Tiny self-contained Schwarzschild + RK4 helpers
//  (mirrors src/2D/2D-physics/{schwarzschild,geodesic,integrator}.hpp
//   but kept local to avoid cross-module include dependencies)
// ────────────────────────────────────────────────────────────
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

} // namespace overlay_phys

// ────────────────────────────────────────────────────────────
//  Overlay shader sources
// ────────────────────────────────────────────────────────────
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

// ────────────────────────────────────────────────────────────
//  PhysicsOverlay: owns GL resources + vertex caches for the two overlays.
// ────────────────────────────────────────────────────────────
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
    // rebuild — without this, sliders that resize the BH or disk at runtime
    // would leave the spacetime well and photon fan stuck at the old scale.
    void notifyScale(float bhRadius, float diskOuterRadius) {
        if (std::abs(bhRadius - lastBHRadius_)         > 1e-5f ||
            std::abs(diskOuterRadius - lastDiskOuter_) > 1e-4f) {
            lastBHRadius_  = bhRadius;
            lastDiskOuter_ = diskOuterRadius;
            dirty_ = true;
        }
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
};

} // namespace bh3d
