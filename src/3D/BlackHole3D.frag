// BlackHole3D.frag
// GLSL Fragment Shader for a Ray-Traced Black Hole Effect
// Author: Copilot (GPT-4.1)
// Usage: Render a full-screen quad with this shader. Provide a background texture and set uniforms as needed.

#version 330 core // Fun fact: This was released back in 2010. Should I use it? maybe not, but for compatibility sake we will make do.

in vec2 fragUV; // UV coordinates from vertex shader (0,0) to (1,1)
out vec4 FragColor;

uniform sampler2D backgroundTex; // Background (star field, accretion disk, etc.)
uniform sampler2D diskTex;       // Accretion disk RGBA (polar mapped)
uniform vec2 resolution;         // Screen resolution
uniform vec3 cameraPos;          // Camera position in world space
uniform vec3 cameraDir;          // Camera forward direction
uniform vec3 cameraUp;           // Camera up vector
uniform float fov;               // Field of view in degrees
uniform float blackHoleRadius;   // Schwarzschild radius (event horizon)
uniform vec3 blackHolePos;       // Black hole position in world space

uniform float diskInnerRadius;   // world units
uniform float diskOuterRadius;   // world units
uniform float diskHalfThickness; // world units (half thickness)

uniform float jetRadius;         // world units
uniform float jetLength;         // world units (each direction)
uniform vec3  jetColor;          // emissive color

uniform int   showDoppler;       // Toggle Doppler shift (V key)
uniform int   showBlueshift;     // Toggle gravitational blueshift of background (U key)
uniform float spinParameter;     // Kerr spin parameter (0-1)

uniform float uTime;             // Elapsed time for animation
uniform int   showJets;          // Toggle jet emission (J key)     [ADDED 2026-04-24: was unconditional]
uniform int   showBLR;           // Toggle Broad Line Region          [ADDED 2026-04-24: uniform existed on CPU, never declared here]
uniform float blrInnerRadius;    // BLR inner edge (world units)     [ADDED 2026-04-24]
uniform float blrOuterRadius;    // BLR outer edge (world units)     [ADDED 2026-04-24]
uniform float blrThickness;      // BLR vertical half-thickness      [ADDED 2026-04-24]
uniform float blrStrength;       // [0,1] cloud opacity scaler per black hole
uniform int   maxStepsOverride;  // Runtime step limit [ADDED 2026-05-25: was previously dead]

// ── Orbiting bodies [ADDED 2026: ported from the photoreal shader so the ──
//    low-quality path shows orbiting stars/clouds, tidal-disruption victims
//    and the inspiralling secondary black hole during a merger.
const int MAX_ORB_BODIES = 10;
uniform vec3  orbBodyPos[MAX_ORB_BODIES];    // world-space centres
uniform float orbBodyRadius[MAX_ORB_BODIES]; // world-unit radii
uniform vec3  orbBodyColor[MAX_ORB_BODIES];  // base colours
uniform int   orbBodyType[MAX_ORB_BODIES];   // BODY_* archetype
uniform int   numOrbBodies;                  // active count (<= MAX_ORB_BODIES)
uniform int   showOrbBody;                   // master toggle

// Body archetypes — must mirror profiles::GalaxyBody3DType, with type 7 the
// transient secondary black hole injected during a merger.
#define BODY_STAR      0
#define BODY_GASCLOUD  1
#define BODY_CLUSTER   2
#define BODY_DGAL      3
#define BODY_NEUTRON   4
#define BODY_WDWARF    5
#define BODY_COMPANION 6
#define BODY_BLACKHOLE 7

// ============================================================================
// NOISE UTILITIES
// ============================================================================
// Value hash, maps a vec3 lattice point to a pseudo-random float [0,1].
// Classic permutation-free hash from Dave Hoskins.
float hash31(vec3 p) {
    p = fract(p * vec3(443.897, 441.423, 437.195));
    p += dot(p, p.yzx + 19.19);
    return fract((p.x + p.y) * p.z);
}

float hash11(float p) {
    return fract(sin(p * 127.1) * 43758.5453);
}

// Smooth 3-D value noise, trilinear interpolation of random lattice values.
float valueNoise3(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    vec3 u = f * f * (3.0 - 2.0 * f); // smoothstep curve

    float v000 = hash31(i + vec3(0,0,0));
    float v100 = hash31(i + vec3(1,0,0));
    float v010 = hash31(i + vec3(0,1,0));
    float v110 = hash31(i + vec3(1,1,0));
    float v001 = hash31(i + vec3(0,0,1));
    float v101 = hash31(i + vec3(1,0,1));
    float v011 = hash31(i + vec3(0,1,1));
    float v111 = hash31(i + vec3(1,1,1));

    return mix(
        mix(mix(v000, v100, u.x), mix(v010, v110, u.x), u.y),
        mix(mix(v001, v101, u.x), mix(v011, v111, u.x), u.y),
        u.z);
}

// Fractional Brownian Motion, 4 octaves of value noise.
// Gives the clumpy, self-similar structure of turbulent astrophysical gas.
float fbm(vec3 p) {
    float v = 0.0;
    float a = 0.5;
    vec3  shift = vec3(100.0);
    for (int o = 0; o < 4; ++o) {
        v += a * valueNoise3(p);
        p  = p * 2.1 + shift;
        a *= 0.5;
    }
    return v;
}

// Constants
// [FIXED 2026-05-25] MAX_STEPS is now the hard upper bound for the GLSL loop; the actual
// per-frame iteration count is clamped from maxStepsOverride (CPU side sends 200 for FAST,
// 300 for CINEMATIC). Early-break keeps perf identical to old 140 when override is low.
const int MAX_STEPS = 320;
const float EPSILON = 1e-4;
const float PI = 3.14159265359;

// Utility: Convert screen UV to world ray direction - standard FPS type camera behavior.
vec3 getRayDir(vec2 uv, vec3 camPos, vec3 camDir, vec3 camUp, float fov) {
    float aspect = resolution.x / resolution.y;
    float px = (uv.x - 0.5) * 2.0 * aspect;
    float py = (uv.y - 0.5) * 2.0;
    float angle = radians(fov * 0.5);
    float dz = 1.0 / tan(angle);
    vec3 right = normalize(cross(camDir, camUp));
    vec3 up = normalize(camUp);
    vec3 forward = normalize(camDir);
    return normalize(px * right + py * up + dz * forward);
}

bool rayHitsSphere(vec3 ro, vec3 rd, vec3 center, float radius) {
    vec3 oc = ro - center;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - radius * radius;
    float h = b * b - c;
    if (h < 0.0) return false;
    float t = -b - sqrt(h);
    return t > 0.0;
}

// Simple lensing: iteratively bend the ray toward the black hole.
// Not physically exact, but stable and produces a convincing distortion.
vec3 sampleBackground(vec3 samplePos, vec3 bhPos) {
    vec3 rel = normalize(samplePos - bhPos);
    float u = 0.5 + atan(rel.z, rel.x) / (2.0 * PI);
    float v = 0.5 - asin(rel.y) / PI;
    vec3 bg = texture(backgroundTex, vec2(u, v)).rgb;
    // Boost stars without lifting the dark sky: brighten highlights, keep blacks
    float lum = dot(bg, vec3(0.2126, 0.7152, 0.0722));
    float boost = smoothstep(0.05, 0.5, lum) * 1.8 + 1.0;
    bg *= boost;
    return bg;
}

// Disk is assumed to lie in the XZ plane at bhPos (normal = +Y).
vec4 sampleDiskRGBA(vec3 pos, vec3 bhPos) {
    vec3 rel = pos - bhPos;
    float r = length(rel.xz);
    if (r < diskInnerRadius || r > diskOuterRadius) return vec4(0.0);
    if (abs(rel.y) > diskHalfThickness) return vec4(0.0);

    // Keplerian differential rotation: inner disk orbits faster
    float omega = 1.5 * pow(max(r, diskInnerRadius), -1.5);
    float rotAngle = uTime * omega;
    float cosR = cos(rotAngle);
    float sinR = sin(rotAngle);
    vec2 rotXZ = vec2(rel.x * cosR - rel.z * sinR,
                      rel.x * sinR + rel.z * cosR);

    // diskTex is a planar RGBA ring image centered in the texture.
    vec2 uv = rotXZ / (diskOuterRadius * 2.0) + 0.5;
    vec4 d = texture(diskTex, uv);
    // Fade by radius so the hard edge of the texture doesn't look cut-out.
    float radFade = smoothstep(diskOuterRadius, diskOuterRadius * 0.92, r) *
                    smoothstep(diskInnerRadius, diskInnerRadius * 1.10, r);
    d.a *= radFade;

    // ── Magnetic instability / MRI turbulence knots ────────────────────────
    // The magnetorotational instability (MRI) creates clumpy density structures
    // in the disk. Model them as animated fBm noise co-rotating with the disk.
    // Noise is evaluated in a frame that rotates with the local Keplerian speed
    // so knots appear stationary in the disk (not in the lab frame).
    float phi_rot = atan(rotXZ.y, rotXZ.x); // azimuth in co-rotating frame
    vec3 noiseCoord = vec3(
        r * cos(phi_rot + uTime * 0.03),
        r * sin(phi_rot + uTime * 0.03),
        uTime * 0.07            // slow vertical shimmer
    ) * 0.55;                   // scale: controls knot spatial frequency
    float turb = fbm(noiseCoord);
    // Shape: bright knots (turb > 0.55) pop out against a dimmer average disk.
    // The 2x factor keeps average brightness neutral; clamp prevents blow-out.
    float knotBright = clamp(turb * 2.0 - 0.5, 0.5, 2.5);
    // Knots are hotter, bias color toward blue-white at peak intensity
    vec3 knotTint = mix(vec3(1.0, 0.85, 0.6), vec3(1.1, 1.05, 1.3),
                        smoothstep(0.55, 1.0, turb));
    d.rgb *= knotBright * knotTint;

    return d;
}

// ============================================================================
// SPECTRAL FREQUENCY SHIFT
// ============================================================================
// Redistributes RGB energy to approximate a spectral wavelength shift.
// freqRatio > 1: blueshift (light gained energy), < 1: redshift (light lost energy).
vec3 applyFrequencyShift(vec3 color, float freqRatio) {
    float shift = freqRatio - 1.0;
    vec3 shifted = color;

    if (shift > 0.0) {
        // Blueshift: energy migrates R→G→B (shorter wavelengths)
        // Phase 1 (mild, freqRatio 1–3): visible color shift toward blue
        float s1 = clamp(shift, 0.0, 2.0);
        shifted.r *= (1.0 - s1 * 0.45);
        shifted.g = color.g * (1.0 - s1 * 0.15) + color.r * s1 * 0.25;
        shifted.b = color.b + color.g * s1 * 0.40 + color.r * s1 * 0.10;

        // Phase 2 (strong, freqRatio 3–8): light shifts into UV.
        // Visible light fades to blue-violet, then desaturates toward white
        // as the entire Planck curve floods all visible channels.
        float s2 = clamp(shift - 2.0, 0.0, 5.0) / 5.0;
        vec3 uvWhite = vec3(0.7, 0.75, 1.0);
        shifted = mix(shifted, uvWhite * dot(shifted, vec3(0.33)), s2);

        // Phase 3 (extreme, freqRatio >8): light shifts beyond visible entirely.
        // Observer sees diminishing violet glow, then near-black as all photon
        // energy has been boosted past the visible band into X-ray/gamma.
        float s3 = clamp(shift - 7.0, 0.0, 8.0) / 8.0;
        float fadeToInvisible = 1.0 - s3 * s3;
        shifted *= fadeToInvisible;
    } else {
        // Redshift: energy migrates B→G→R (longer wavelengths)
        float s1 = clamp(-shift, 0.0, 2.0);
        shifted.b *= (1.0 - s1 * 0.45);
        shifted.g = color.g * (1.0 - s1 * 0.15) + color.b * s1 * 0.25;
        shifted.r = color.r + color.g * s1 * 0.40 + color.b * s1 * 0.10;

        // Deep redshift: light fades into infrared, disappearing from visible
        float s2 = clamp(-shift - 2.0, 0.0, 5.0) / 5.0;
        shifted *= (1.0 - s2 * 0.85);
    }

    // Relativistic intensity: I_obs ∝ D^3 (Liouville's theorem)
    // Cap raised to 20 [FIXED 2026-04-24: was lowered to 8 which under-lit the disk;
    // 50 was the original (blew out large-disk configs), 20 is the balance]
    float brightFactor = pow(clamp(freqRatio, 0.05, 12.0), 2.5);
    brightFactor = min(brightFactor, 20.0);
    return max(shifted * brightFactor, vec3(0.0));
}

// Physically motivated Doppler + gravitational redshift for the accretion disk.
// Computes Keplerian orbital velocity, relativistic Doppler factor, and
// gravitational redshift corrected for observer position.
vec3 applyDoppler(vec3 color, vec3 rel, vec3 viewDir) {
    float r = length(rel.xz);
    float rs = blackHoleRadius;

    // Keplerian orbital velocity: v/c = √(Rs / (2r))
    float v = sqrt(rs / (2.0 * max(r, rs * 2.0)));
    v = clamp(v, 0.0, 0.7);

    // Tangential velocity direction (prograde, counter-clockwise in XZ plane)
    vec3 velDir = normalize(vec3(-rel.z, 0.0, rel.x));
    float cosTheta = dot(velDir, -viewDir);

    // Special relativistic Doppler factor: D = 1/(γ(1 - β·cosθ))
    float gamma = 1.0 / sqrt(max(1.0 - v * v, 0.01));
    float D = 1.0 / (gamma * (1.0 - v * cosTheta));
    D = clamp(D, 0.15, 5.0);

    // Toggle Doppler on/off
    if (showDoppler == 0) D = 1.0;

    // Gravitational redshift: f_obs/f_emit = √((1-Rs/r_emit)/(1-Rs/r_obs))
    // Accounts for observer position in the gravitational well
    float rObs = length(cameraPos - blackHolePos);
    float gEmit = sqrt(max(1.0 - rs / max(r, rs * 1.01), 0.01));
    float gObs  = sqrt(max(1.0 - rs / rObs, 0.01));
    float grav  = gEmit / max(gObs, 0.05);

    // Combined frequency shift: Doppler × gravitational
    float freqShift = D * grav;

    // Apply spectral shift and relativistic beaming to the color
    return applyFrequencyShift(color, freqShift);
}

// ============================================================================
// JET EMISSION
// ============================================================================
// Volumetric relativistic jets along ±Y with:
//   • Noise knots (MHD shocks) animated outward at ~0.9c
//   • Relativistic Doppler beaming: jet aimed toward viewer appears longer/brighter
//
// Doppler beaming derivation:
//   β = jet speed / c (0.9)
//   cosθ = dot(jetDir, -viewDir)   (positive when jet points toward viewer)
//   D = 1 / (γ(1 − βcosθ))        (relativistic Doppler factor)
//   I_obs ∝ D^3                    (Liouville's theorem for a continuous jet)
// The counter-jet (pointing away) has cosθ < 0 → D < 1 → dimmer / shorter apparent length.
//
vec3 jetEmission(vec3 pos, vec3 bhPos, vec3 viewDir) {
    vec3 rel = pos - bhPos;
    float y  = rel.y;
    float ay = abs(y);
    if (ay < diskHalfThickness) return vec3(0.0);
    if (ay > jetLength)         return vec3(0.0);

    // ── Relativistic Doppler beaming ──────────────────────────────────────
    // Use D^2 (brightness modulation) rather than D^3 (specific intensity)
    // because in this ray-marcher each step already accumulates depth; D^3
    // per-step would compound and blow out entirely when the jet faces us.
    // D^2 still gives a dramatic approaching/receding asymmetry while keeping
    // structure visible.  Cap is tight (6×) so knot detail isn't washed out.
    const float BETA    = 0.92;
    float gamma_rel     = 1.0 / sqrt(1.0 - BETA * BETA);   // ≈ 2.55
    vec3  jetDir        = vec3(0.0, sign(y), 0.0);
    float cosTheta      = dot(jetDir, -viewDir);
    float D             = 1.0 / max(gamma_rel * (1.0 - BETA * cosTheta), 0.05);
    float beaming       = clamp(D * D, 0.05, 6.0);
    float effectiveLength = jetLength * clamp(D * 0.5 + 0.5, 0.2, 2.0);
    if (ay > effectiveLength) return vec3(0.0);

    // ── Episodic accretion bursts ─────────────────────────────────────────
    // Accretion onto a SMBH like TON 618 is wildly variable, the disk is
    // magnetically unstable and dumps clumps of material onto the jet base
    // on timescales of minutes to hours (compressed here for visibility).
    // Model as a smooth low-frequency noise (≈0.07 Hz) so the whole jet
    // occasionally surges then dims, with no two cycles identical.
    float burstNoise = valueNoise3(vec3(uTime * 0.07, 3.7, 1.2));
    // Map [0,1] → [0.25, 1.8] so dim troughs are still visible but peaks blaze
    float burstPulse = 0.25 + 1.55 * burstNoise;

    // ── Kelvin-Helmholtz helical wobble ───────────────────────────────────
    // Shear between the jet and surrounding gas drives KHI, which manifests
    // as a slow helical/sinusoidal oscillation of the jet spine.
    // Amplitude grows with distance from the BH (jets widen and destabilise).
    // KHI helical wobble, amplitude 40% of jet radius, grows toward tip
    float wobbleAmp   = 0.40 * jetRadius * (ay / jetLength);
    float wobbleFreq  = 4.0 / max(jetLength, EPSILON);
    float wobbleSpeed = 0.35;
    float wobblePhase = ay * wobbleFreq - uTime * wobbleSpeed;
    vec2  wobble      = wobbleAmp * vec2(sin(wobblePhase), cos(wobblePhase * 0.71));

    vec2 coreOffset = rel.xz - wobble;
    float radial    = length(coreOffset);
    if (radial > jetRadius * 2.0) return vec3(0.0);

    float core  = exp(-radial * radial / max(jetRadius * jetRadius * 0.25, EPSILON));
    float along = smoothstep(jetLength, jetLength * 0.2, ay);

    // ── Shock knots: advected outward at β ~ 0.9c ────────────────────────
    const float BETA_KNOT = 0.90;
    float travelOffset    = uTime * BETA_KNOT;

    float scaleCoarse = 2.5  / max(jetLength, EPSILON);
    float scaleFine   = 13.0 / max(jetLength, EPSILON);

    // Stochastic knot injection: the rate at which new knots are launched
    // varies over time. We discretise time into ~4-second "episodes" and
    // assign each a random injection strength, then smoothly interpolate.
    // This creates clusters of knots followed by quiet spells.
    float episodeLen = 4.0;
    float epIdx      = floor(uTime / episodeLen);
    float epFrac     = fract(uTime / episodeLen);
    float injCurrent = hash11(epIdx);
    float injNext    = hash11(epIdx + 1.0);
    // Smoothstep so injection ramps up/down rather than snapping
    float injRate    = mix(injCurrent, injNext, smoothstep(0.0, 1.0, epFrac));
    // Low-injection episodes: faint, sparse knots. High-injection: dense bright chain.
    float knotScale  = 0.3 + injRate * 1.4;

    vec3 knotPosCoarse = vec3(coreOffset * scaleCoarse,
                              (ay - travelOffset) * scaleCoarse);
    vec3 knotPosFine   = vec3(coreOffset * scaleFine,
                              (ay - travelOffset * 1.15) * scaleFine);

    float knotCoarse = valueNoise3(knotPosCoarse);
    float knotFine   = fbm(knotPosFine);
    float knot       = (knotCoarse * 0.65 + knotFine * 0.35) * knotScale;

    // High contrast: dark lanes (knot < 0.4) are near-black, bright shocks pop
    float shockIntensity  = smoothstep(0.40, 0.68, knot);
    // baseEmission is low so inter-knot regions are genuinely dim
    float baseEmission    = 0.08 + 0.12 * knot;
    float emissionProfile = (baseEmission + shockIntensity * 2.2) * burstPulse * 0.55;

    vec3 knotColor = mix(jetColor,
                         jetColor * vec3(0.8, 0.92, 1.5),
                         shockIntensity * 0.8);

    // ── Bright collimated spine ───────────────────────────────────────────
    // A narrow hot channel down the jet axis keeps the beam reading as a
    // focused, relativistic ejection rather than a diffuse cone.
    float spine    = exp(-radial * radial / max(jetRadius * jetRadius * 0.06, EPSILON));
    vec3  spineCol = mix(knotColor, vec3(0.85, 0.95, 1.45), 0.6);   // hot blue-white core

    // ── Launch / collimation glow at the jet base ─────────────────────────
    // Material is brightest where it is accelerated and collimated just above
    // the disc; this anchors the jet visually to the black hole.
    float baseFall = exp(-ay / max(jetLength * 0.16, EPSILON));
    float baseGlow = baseFall * exp(-radial * radial / max(jetRadius * jetRadius * 0.5, EPSILON));

    vec3 emission = knotColor * core     * emissionProfile
                  + spineCol  * spine    * (0.5 + 0.9 * shockIntensity) * burstPulse * 0.6
                  + spineCol  * baseGlow * burstPulse * 0.9;

    return emission * along * beaming;
}

void traceBentRay(
    in vec3 ro,
    in vec3 rd,
    in vec3 bhPos,
    in float bhRadius,
    out vec3 outPos,
    out vec3 outDir,
    out bool absorbed,
    out vec3 outJet,
    out vec4 outDisk,
    out vec3 outBLR)
{
    vec3 pos = ro;
    vec3 dir = normalize(rd);
    absorbed = false;
    outJet = vec3(0.0);
    outDisk = vec4(0.0);
    outBLR = vec3(0.0);

    // Scene radius scales with the largest user-configured feature so the full disk/jets are always reached.
    // [FIXED 2026-04-24: was hardcoded MAX_DIST=140.0, which clipped outer disk/jets on large configs]
    float sceneRadius = max(max(diskOuterRadius, jetLength) * 1.5, 140.0);
    // [FIXED 2026-05-25] Step count now driven by CPU-side maxStepsOverride (200=fast, 300=cinematic).
    // Clamped to [80, MAX_STEPS] for safety against unset/garbage uniform values.
    int   nSteps   = clamp(maxStepsOverride, 80, MAX_STEPS);
    float baseStep = sceneRadius / float(nSteps);
    // [FIXED 2026-05-25] Deflection strength bumped 1.0→1.5 so grazing rays actually wrap the
    // shadow instead of producing a Newtonian-looking ring with a giant empty gap.
    float k = bhRadius * bhRadius * 1.5;

    float lastDiskSide = (pos - bhPos).y;
    vec3  prevPos      = pos;

    for (int i = 0; i < MAX_STEPS; ++i) {
        if (i >= nSteps) break;
        vec3 toBH = bhPos - pos;
        float r = length(toBH);

        // [FIXED 2026-05-25] Adaptive step size: shrink to ~35% of baseStep within 3 Rs of the
        // event horizon, ramping back to full baseStep by 20 Rs. This gives the integrator
        // enough resolution near the BH to bend grazing rays around the shadow (Interstellar
        // top-arc lensing) without paying the cost in the empty regions far from the hole.
        float stepLen = baseStep * mix(0.35, 1.0,
                                       smoothstep(3.0 * bhRadius, 20.0 * bhRadius, r));

        // Exact EH absorption: check if this step segment intersects the event horizon sphere.
        // [FIXED 2026-04-24: replaced the old 'r < bhRadius + EPSILON' point-check, which missed
        //  near-grazing rays when stepLen ≥ bhRadius (the step would skip over the horizon),
        //  producing the bright seam / 'hole' at the shadow edge. Sphere intersection is exact
        //  regardless of step size, so there is no need to change stepLen near the BH.]
        {
            float tca = dot(toBH, dir);     // signed distance to closest approach
            float d2  = dot(toBH, toBH) - tca * tca;  // squared miss distance
            if (d2 < bhRadius * bhRadius) {
                float thc = sqrt(bhRadius * bhRadius - d2);
                float tHit = tca - thc;     // entry intersection along ray
                if (tHit >= 0.0 && tHit <= stepLen) {
                    absorbed = true;
                    break;
                }
            }
        }

        // Jet emission, gated by showJets toggle.
        // Pass the original ray direction (rd) for Doppler beaming calculations;
        // the bent 'dir' at this sample point would give incorrect beaming angles.
        if (showJets != 0)
            outJet += jetEmission(pos, bhPos, rd) * stepLen;

        // BLR cloud field, discrete gas blobs orbiting at high speed.
        // Controlled by blrStrength [0,1]: TON 618 densest, non-AGN zero.
        if (showBLR != 0 && blrStrength > 0.0) {
            float blrAbsY = abs((pos - bhPos).y);
            if (r > blrInnerRadius && r < blrOuterRadius && blrAbsY < blrThickness) {
                float radFade  = smoothstep(blrInnerRadius, blrInnerRadius * 1.2, r)
                               * smoothstep(blrOuterRadius, blrOuterRadius * 0.85, r);
                float vertFade = smoothstep(blrThickness, blrThickness * 0.4, blrAbsY);

                // Clouds are fast-moving: use angle + time to scatter them.
                // azimuth φ changes with 1/r^0.5 (Keplerian tangential speed)
                vec3 rel = pos - bhPos;
                float phi = atan(rel.z, rel.x) + uTime * 0.18 / sqrt(max(r, 0.1));

                // Three octaves of 3-D value noise to break up the ribbon.
                vec3 cloudP = vec3(phi * 3.5, r * 0.22, rel.y * 0.35 + uTime * 0.04);
                float cloud = valueNoise3(cloudP * 1.0)
                            + 0.50 * valueNoise3(cloudP * 2.1)
                            + 0.25 * valueNoise3(cloudP * 4.3);
                cloud = cloud / 1.75;          // normalise to [0,1]
                cloud = pow(cloud, 1.4);        // mild sharpening, keeps dense filling
                cloud *= radFade * vertFade * blrStrength;

                // BLR line emission: Hα red core → blue-shifted CIV outer wing
                vec3 blrColor = mix(vec3(1.0, 0.45, 0.2),   // bright Hα orange-red
                                    vec3(0.6, 0.82, 1.0),    // blue CIV wing
                                    smoothstep(blrInnerRadius, blrOuterRadius, r));
                outBLR += blrColor * 0.22 * cloud * stepLen;
            }
        }

        // Detect crossing of the disk plane and sample if within thickness.
        float diskSide = (pos - bhPos).y;
        if (sign(diskSide) != sign(lastDiskSide)) {
            // Approximate intersection point at the plane.
            float t = lastDiskSide / max(lastDiskSide - diskSide, EPSILON);
            // [FIXED 2026-05-25] Use prevPos (true previous sample) instead of pos - dir*baseStep,
            // since stepLen is now adaptive and dir has been updated by the bend at this iteration.
            vec3 hitPos = mix(prevPos, pos, clamp(t, 0.0, 1.0));
            vec4 d = sampleDiskRGBA(hitPos, bhPos);
            if (d.a > 0.001) {
                outDisk = d;
                outPos = hitPos;
                outDir = dir;
                return;
            }
        }
        lastDiskSide = diskSide;

        vec3 pullDir = toBH / max(r, EPSILON);
        // Inverse-square-ish pull; softened for stability
        float strength = k / (r * r + bhRadius * bhRadius);
        dir = normalize(dir + pullDir * strength * stepLen);
        prevPos = pos;
        pos += dir * stepLen;
    }

    outPos = pos;
    outDir = dir;
}

// ============================================================================
// ORBITING BODIES  [ADDED 2026: low-quality parity with the photoreal shader]
// ============================================================================
// Surface shading + additive glow for the bodies that orbit the black hole
// (stars, gas clouds, clusters, remnants) and for the inspiralling secondary
// black hole during a merger (BODY_BLACKHOLE). Simplified ports of the
// photoreal shader's shadeOrbBody / orbBodyEmission so the fast render path
// shows the same objects, tidal-disruption victims and merger companions.

// Nearest positive ray-sphere intersection distance, or -1.0 on a miss.
float intersectSphereT(vec3 ro, vec3 rd, vec3 center, float radius) {
    vec3  oc = ro - center;
    float b  = dot(oc, rd);
    float c  = dot(oc, oc) - radius * radius;
    float disc = b * b - c;
    if (disc < 0.0) return -1.0;
    float s  = sqrt(disc);
    float t0 = -b - s;
    if (t0 > EPSILON) return t0;
    float t1 = -b + s;
    if (t1 > EPSILON) return t1;
    return -1.0;
}

// Which body types present a hard, light-blocking surface.
bool isOpaqueBody(int t) {
    return t == BODY_STAR || t == BODY_COMPANION
        || t == BODY_NEUTRON || t == BODY_WDWARF
        || t == BODY_BLACKHOLE;
}

// Hard-surface shading for opaque bodies. The secondary black hole returns
// black (its photon ring is drawn additively in bodyEmissionSimple).
vec3 shadeBodySimple(int type, vec3 hitPos, vec3 normal, vec3 center, vec3 col) {
    if (type == BODY_BLACKHOLE) return vec3(0.0);

    vec3  view = normalize(hitPos - cameraPos);
    float mu   = max(dot(normal, -view), 0.0);     // 1 at disc centre, 0 at limb
    float limb = 0.4 + 0.6 * mu;                    // Eddington limb darkening

    // Gravitational dimming as the body approaches the primary's horizon.
    float rBH  = length(hitPos - blackHolePos);
    float grav = sqrt(max(1.0 - blackHoleRadius / max(rBH, blackHoleRadius * 1.01), 0.05));

    if (type == BODY_NEUTRON) {
        vec3 hot = mix(col, vec3(1.0), pow(mu, 2.0) * 0.6);
        return hot * (0.55 + 0.45 * mu) * 12.0 * grav;
    }
    if (type == BODY_WDWARF) {
        vec3 hot = mix(col, vec3(1.0), 0.7);
        return hot * (0.5 + 0.5 * mu) * 9.0 * grav;
    }

    // Star / companion: warm limb-darkened photosphere with faint granulation
    // and tidal reddening when it strays close to the black hole.
    vec3  surfP   = (hitPos - center) * (type == BODY_COMPANION ? 3.0 : 4.0);
    float granMod = 0.85 + 0.30 * fbm(surfP);
    vec3  hot     = mix(col, vec3(1.0, 0.97, 0.90), 0.40);
    vec3  limbCol = mix(col, vec3(1.0, 0.55, 0.20), 0.35);
    vec3  photo   = mix(limbCol, hot, mu);
    float tidal   = smoothstep(3.0, 1.5, rBH / blackHoleRadius);
    photo = mix(photo, vec3(1.0, 0.30, 0.10), tidal * 0.6);
    return photo * limb * granMod * (type == BODY_COMPANION ? 5.0 : 6.0) * grav;
}

// Additive emission for every body: stellar coronae, volumetric clouds, and
// the gold photon ring of the inspiralling secondary black hole.
vec3 bodyEmissionSimple(int idx, vec3 ro, vec3 rd) {
    int   type = orbBodyType[idx];
    vec3  c    = orbBodyPos[idx];
    float r    = orbBodyRadius[idx];
    vec3  col  = orbBodyColor[idx];

    // Closest approach of the ray to the body centre.
    float sClosest = max(dot(c - ro, rd), 0.0);
    vec3  closestP = ro + rd * sClosest;
    float dPerp    = length(closestP - c);

    if (type == BODY_BLACKHOLE) {
        // Photon ring + warm glow + cool halo, gated to the outside of the disc.
        float outside = smoothstep(r * 0.82, r * 1.0, dPerp);
        float dRing   = dPerp - r * 1.10;
        float ring    = exp(-dRing * dRing / max(r * r * 0.0324, EPSILON));   // (0.18r)^2
        float glow    = exp(-max(dPerp - r, 0.0) / (r * 0.85));
        float halo    = exp(-max(dPerp - r, 0.0) / (r * 3.0)) * 0.25;
        vec3  ringCol = vec3(1.0, 0.84, 0.52);
        vec3  haloCol = vec3(0.45, 0.62, 1.0);
        return (ringCol * ring * 3.4 + ringCol * glow * 0.45 + haloCol * halo) * outside;
    }
    if (type == BODY_GASCLOUD) {
        float sigma = 1.5 * r;
        float env   = exp(-dPerp * dPerp / (2.0 * sigma * sigma));
        float turb  = 0.55 + 0.9 * fbm(closestP * 0.6);
        vec3  tint  = mix(col, vec3(1.0, 0.85, 0.7), 0.25);
        return tint * env * turb * 2.0;
    }
    if (type == BODY_CLUSTER) {
        float sigma = 1.2 * r;
        float halo  = exp(-dPerp * dPerp / (2.0 * sigma * sigma));
        return mix(col, vec3(1.0), 0.4) * halo * 1.4;
    }
    if (type == BODY_DGAL) {
        float sigma = 1.4 * r;
        float env   = exp(-dPerp * dPerp / (2.0 * sigma * sigma));
        float core  = exp(-dPerp * dPerp / (2.0 * (0.35 * r) * (0.35 * r)));
        return mix(col, vec3(0.85, 0.9, 1.0), 0.35) * (env * 0.7 + core * 1.4);
    }

    // Opaque stellar bodies: inverse-square corona halo outside the surface.
    float halo = r / max(dPerp, 0.5 * r);
    halo = halo * halo;
    halo *= smoothstep(6.0 * r, 1.0 * r, dPerp);
    vec3  coronaCol = col;
    float bright    = 0.6;
    if      (type == BODY_NEUTRON) { coronaCol = mix(col, vec3(0.6, 0.8, 1.0), 0.5); bright = 1.1; }
    else if (type == BODY_WDWARF)  { coronaCol = mix(col, vec3(0.8, 0.9, 1.0), 0.4); bright = 0.7; }
    return coronaCol * halo * bright;
}

void main() {
    // Get ray direction from camera through this pixel
    vec3 rayDir = getRayDir(fragUV, cameraPos, cameraDir, cameraUp, fov);
    vec3 rayOrigin = cameraPos;

    // ── Orbiting bodies (straight-ray approximation) ──────────────────────
    // Bodies are tested against the un-bent ray: lensing of the bodies
    // themselves is a second-order effect at their orbital distance, and the
    // straight ray fixes their screen-space position. Track the nearest
    // opaque surface hit (for occlusion) and additively accumulate every
    // body's glow / corona / photon ring.
    float orbBodyDist  = -1.0;
    int   hitBodyIndex = -1;
    vec3  coronaAcc    = vec3(0.0);
    if (showOrbBody != 0) {
        for (int i = 0; i < numOrbBodies; ++i) {
            if (isOpaqueBody(orbBodyType[i])) {
                float t = intersectSphereT(rayOrigin, rayDir, orbBodyPos[i], orbBodyRadius[i]);
                if (t > 0.0 && (orbBodyDist < 0.0 || t < orbBodyDist)) {
                    orbBodyDist  = t;
                    hitBodyIndex = i;
                }
            }
            coronaAcc += bodyEmissionSimple(i, rayOrigin, rayDir);
        }
    }
    vec3 orbBodyHitColor = vec3(0.0);
    if (hitBodyIndex >= 0) {
        vec3 hp  = rayOrigin + rayDir * orbBodyDist;
        vec3 nrm = normalize(hp - orbBodyPos[hitBodyIndex]);
        orbBodyHitColor = shadeBodySimple(orbBodyType[hitBodyIndex], hp, nrm,
                                          orbBodyPos[hitBodyIndex],
                                          orbBodyColor[hitBodyIndex]);
    }

    // Bend the ray (gravitational lensing approximation)
    vec3 bentPos;
    vec3 bentDir;
    bool absorbed;
    vec3 jetAcc;
    vec4 diskHit;
    vec3 blrAcc;
    traceBentRay(rayOrigin, rayDir, blackHolePos, blackHoleRadius, bentPos, bentDir, absorbed, jetAcc, diskHit, blrAcc);

    // ── Base scene colour + its depth along the ray ───────────────────────
    vec3  color;
    float sceneDist;
    if (absorbed) {
        // Ray fell through the event horizon: black silhouette. Depth is the
        // straight-ray distance to the horizon so a body in front still wins.
        color = vec3(0.0);
        float eh = intersectSphereT(rayOrigin, rayDir, blackHolePos, blackHoleRadius);
        sceneDist = eh > 0.0 ? eh : length(bentPos - rayOrigin);
    } else {
        // Gravitational blueshift factor for background starlight.
        float rCamera = length(cameraPos - blackHolePos);
        float gCamera = sqrt(max(1.0 - blackHoleRadius / rCamera, 0.001));
        float bgFreqShift = (showBlueshift != 0 && gCamera < 0.98) ? (1.0 / max(gCamera, 0.02)) : 1.0;
        float sceneRadius = max(max(diskOuterRadius, jetLength) * 1.5, 140.0);

        if (diskHit.a > 0.001) {
            vec3 rel = (bentPos - blackHolePos);
            vec3 diskCol = applyDoppler(diskHit.rgb, rel, bentDir);
            vec3 bg = sampleBackground(bentPos + bentDir * max(diskOuterRadius * 0.5, 40.0), blackHolePos);
            if (bgFreqShift > 1.01) bg = applyFrequencyShift(bg, bgFreqShift);
            color = mix(bg, diskCol, diskHit.a);
            sceneDist = length(bentPos - rayOrigin);
        } else {
            color = sampleBackground(bentPos + bentDir * (sceneRadius * 0.5), blackHolePos);
            if (bgFreqShift > 1.01) color = applyFrequencyShift(color, bgFreqShift);
            sceneDist = 1e9;
        }
    }

    // ── Opaque orbiting body occludes the scene when it is nearer ─────────
    if (orbBodyDist > 0.0 && orbBodyDist < sceneDist) {
        color = orbBodyHitColor;
    }

    // ── Emissive layers on top (jets, BLR, body coronae / photon rings) ───
    color += jetAcc + blrAcc + coronaAcc;

    FragColor = vec4(color, 1.0);
}
