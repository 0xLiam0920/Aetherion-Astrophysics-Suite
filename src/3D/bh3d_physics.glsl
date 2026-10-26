// ============================================================================
// Author: 0xLiam
// description: bh3d_physics.glsl: the (new) shared GLSL used by BOTH BlackHole3D.frag (simple/LQP) and
// BlackHole3D_PhotorealDisk.frag (photoreal/realistic).
//
// initShaders() (bh3d_core.hpp) injects this file
// into each fragment shader right after its #version line, via injectCommon()
// in bh3d_shaderinject.hpp. It is preceded by:
//     #define QUALITY 0   // simple shader   (Verlet, cheap visuals)
//     #define QUALITY 1   // photoreal shader (RK4, full visuals)
// Low quality changes step counts and detail depth
//
// Units applied: Rs = 1, c = 1, so GM = 0.5.
// ============================================================================

// ---- Uniforms used by both shaders ----------------------------------------
uniform vec2  resolution;          // Render-target size in pixels
uniform float spinParameter;       // Dimensionless spin a* (0 = Schwarzschild, 1 = extremal Kerr)
uniform float blackHoleRadius;     // Schwarzschild radius (event horizon)
uniform vec3  blackHolePos;        // Black hole position in world space
uniform float diskInnerRadius;     // Disk inner edge (world units)
uniform float diskOuterRadius;     // Disk outer edge (world units)
uniform int   showDoppler;         // Toggle Doppler beaming/shift (V key)
uniform float uTime;               // Elapsed time for animation

// ---- Disk physics uniforms (Steps 2, 3, 5) --------------------------------
uniform float diskPeakTemp;        // Physical peak disk temperature [K]
uniform int   scientificMode;      // 1 = exact physics, no artistic layers; 0 = cinematic
uniform sampler2D blackbodyRadianceLUT; // 1024x1 linear RGB radiance, log-spaced 300 K..1e8 K
uniform float sciExposure;         // 1 / luminance of the LUT at diskPeakTemp (CPU)
// Novikov-Thorne per-spin constants, computed on the CPU when spin changes.
uniform float ntX0;                // sqrt(2 * isco)
uniform float ntA;                 // clamped spin
uniform float ntPeakInv;           // 1 / peak of the flux shape
uniform float isco;                // prograde ISCO [Rs]
uniform vec3  ntX;                 // roots of x^3 - 3x + 2a
uniform vec3  ntC;                 // log coefficients

// Disk inner edge: the ISCO in scientific mode, the preset value otherwise.
float diskInner() { return scientificMode == 1 ? isco : diskInnerRadius; }

// ---- Camera ray -------------------------------------------------------------
// Convert screen UV to a world-space ray direction (standard FPS-style camera).
vec3 getRayDir(vec2 uv, vec3 camPos, vec3 camDir, vec3 camUp, float fovDeg) {
    float aspect = resolution.x / resolution.y;
    float px = (uv.x - 0.5) * 2.0 * aspect;
    float py = (uv.y - 0.5) * 2.0;
    float angle = radians(fovDeg * 0.5);
    float dz = 1.0 / tan(angle);
    vec3 right = normalize(cross(camDir, camUp));
    vec3 up = normalize(camUp);
    vec3 forward = normalize(camDir);
    return normalize(px * right + py * up + dz * forward);
}

// ---- Hash -------------------------------------------------------------------
// 2D -> 1D hash. Used by the noise code and by the per-pixel first-step jitter
// of the geodesic march (Step 1b).
float hash12(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

// ---- Novikov-Thorne disk temp. theorem ---------------------------------------
// Unnormalised Page-Thorne flux shape at x = sqrt(r_M), given the cubic roots
// of x^3 - 3x + 2a and their log coefficients (Page & Thorne 1974). The
// relativistic Kerr generalisation of the Newtonian [1 - sqrt(rISCO/r)] / r^3.
float ptFluxShape(float x, float x0, float a,
                  float x1, float x2, float x3,
                  float c1, float c2, float c3) {
    if (x <= x0) return 0.0;
    float Q = x - x0 - 1.5 * a * log(x / x0)
            - c1 * log((x - x1) / (x0 - x1))
            - c2 * log((x - x2) / (x0 - x2))
            - c3 * log((x - x3) / (x0 - x3));
    float denom = x * x * x * x * (x * x * x - 3.0 * x + 2.0 * a);
    return max(Q / denom, 0.0);
}

// T(r) / T_peak in [0,1] from the uniforms precomputed by ntSetup() on the CPU.
float ntProfile(float rRs) {
    float x = sqrt(2.0 * rRs);
    float F = ptFluxShape(x, ntX0, ntA, ntX.x, ntX.y, ntX.z, ntC.x, ntC.y, ntC.z);
    return pow(max(F * ntPeakInv, 0.0), 0.25);
}

// ---- Exact disk redshift (Step 2) -------------------------------------------
// g = nu_obs / nu_emit for gas on a prograde circular Kerr orbit, seen by a
// distant static observer. r in Rs. Lz is the conserved photon L/E (Rs units)
// of the traced ray; see photonLz().
float diskRedshift(float r, float Lz) {
    float a  = clamp(spinParameter, 0.0, 0.998);
    float rM = 2.0 * r;                    // Rs -> M units
    float s  = rM * sqrt(rM);              // rM^1.5
    float Om = 1.0 / (s + a);              // gas angular velocity, M units
    float ut = (s + a) / (sqrt(s) * sqrt(max(s - 3.0 * sqrt(rM) + 2.0 * a, 1e-5)));
    float LzM = (showDoppler == 0) ? 0.0 : 2.0 * Lz;   // Doppler off keeps gravitational + transverse
    return 1.0 / (ut * (1.0 + Om * LzM));
}

// Photon angular momentum about the disk axis (Y), in units of Rs. Conserved
// along the whole geodesic, so the camera ray is enough: no need to recompute
// it at the disk hit. Disk gas rotates with velocity (-z, 0, x), so its angular
// momentum points along -Y, hence the minus sign. Camera at +Z looking at the
// hole: gas at +X approaches and comes out blue (g > 1).
float photonLz(vec3 camPos, vec3 rayDir) {
    return -cross(camPos - blackHolePos, rayDir).y / max(blackHoleRadius, 1e-6);
}

// Static observer at finite distance sits inside the well too: blueshift.
float observerShift(vec3 camPos) {
    float rObs = length(camPos - blackHolePos) / max(blackHoleRadius, 1e-6);
    return 1.0 / sqrt(max(1.0 - 1.0 / rObs, 1e-3));
}

// ---- Blackbody radiance ---------------------------------------------
// A blackbody seen with frequency shift g looks
// exactly like a blackbody at g*T (I_nu / nu^3 is invariant), so if an observer looks up the disk, the observed radiance is that of a blackbody at g*T.
// Thank god for my teacher. You're a godsend mr J.
vec3 blackbodyRadiance(float T) {
    float t = log(clamp(T, 300.0, 1e8) / 300.0) / log(1e8 / 300.0);
    return texture(blackbodyRadianceLUT, vec2((0.5 + t * 1023.0) / 1024.0, 0.5)).rgb;
}

// Exact disk emission at radius r (Rs) for a ray with photon angular momentum Lz.
vec3 diskEmissionPhysical(float r, float Lz, vec3 camPos) {
    float T = diskPeakTemp * ntProfile(r);           // Novikov-Thorne
    if (T < 300.0) return vec3(0.0);                 // ISCO edge / below LUT range
    float g = diskRedshift(r, Lz) * observerShift(camPos);
    return blackbodyRadiance(g * T) * sciExposure;
}
