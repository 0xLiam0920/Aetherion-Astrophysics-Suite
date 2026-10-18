#ifndef BH3D_BLACKBODY_HPP
#define BH3D_BLACKBODY_HPP
// ============================================================================
// bh3d_blackbody.hpp, file for GL-free physical colour + disk-temperature helpers.
// NOTE: This file was separately prototyped on a different directory, before being implemented and refined over 1.5 weeks
//
// The reason this exists as its own header is so there's exactly one place
// where the accretion-disk colour physics lives. Since the GPU LUT builder
// and the CPU test suite both need the same numbers, we couldn't have two
// copies of this math drifting apart over time, so everything got pulled out
// here instead.
// NOTE: THIS HEADER HAS ZERO OpenGL / SFML / glm dependency, so it can be shared by the following:
//1. the GPU LUT builder  (texture_utils.cpp :: createBlackbodyLUT), and
//2. the physics regression tests (tests/physics_regression_tests.cpp),
// which means the values the shader samples are guaranteed to be exactly the values CI checks.
//
// Physics:
//   blackbodyRGB(T)  Planck B(λ,T) × CIE-1931 2° CMF over 380-780 nm, converted
//                    to XYZ then D65 sRGB, gamma-encoded to display sRGB in
//                    [0,1]. This is precisely what createBlackbodyLUT() stores
//                    (times 255) and the shader reads back via kelvinToRGBLUT().
//   computeISCO(a)               Exact Kerr prograde ISCO (Bardeen 1972).
//   pageThorneFlux(r,a)          Relativistic Kerr disk flux (Page-Thorne 1974).
//   novikovThorneTemperature()   Normalised spin-dependent T(r;a) (shader-matched).
//   diskColorTemperature()       Visible-band colour anchor used by the
//                                physical disk-colour path.
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>

namespace bh3d {
namespace physics {
// CIE 1931 2° standard-observer colour-matching functions,
// sampled 380–780 nm at 10 nm steps (41 samples).
inline constexpr int kCMFSamples = 41;

inline constexpr double kCMF_x[kCMFSamples] = {
    0.001368, 0.004243, 0.014310, 0.043510, 0.134380, 0.283900, 0.348280,
    0.336200, 0.290800, 0.195360, 0.095640, 0.032010, 0.004900, 0.009300,
    0.063270, 0.165500, 0.290400, 0.433450, 0.594500, 0.762100, 0.916300,
    1.026300, 1.062200, 1.002600, 0.854450, 0.642400, 0.447900, 0.283500,
    0.164900, 0.087400, 0.046770, 0.022700, 0.011359, 0.005790, 0.002899,
    0.001440, 0.000690, 0.000332, 0.000166, 0.000083, 0.000041
};
inline constexpr double kCMF_y[kCMFSamples] = {
    0.000039, 0.000120, 0.000396, 0.001210, 0.004000, 0.011600, 0.023000, // The exact values of this are based on a (presumably human) observer. 
    0.038000, 0.060000, 0.090980, 0.139020, 0.208020, 0.323000, 0.503000,
    0.710000, 0.862000, 0.954000, 0.994950, 0.995000, 0.952000, 0.870000,
    0.757000, 0.631000, 0.503000, 0.381000, 0.265000, 0.175000, 0.107000,
    0.061000, 0.032000, 0.017000, 0.008210, 0.004102, 0.002091, 0.001047,
    0.000520, 0.000249, 0.000120, 0.000060, 0.000030, 0.000015
};
inline constexpr double kCMF_z[kCMFSamples] = {
    0.006450, 0.020050, 0.067850, 0.207400, 0.645600, 1.385600, 1.747060, 
    1.772110, 1.669200, 1.287640, 0.812950, 0.465180, 0.272000, 0.158200,
    0.078250, 0.042160, 0.020300, 0.008750, 0.003900, 0.002100, 0.001650,
    0.001100, 0.000800, 0.000340, 0.000190, 0.000070, 0.000020, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000
};

// hc/k = 0.014388 m·K, the Planck second radiation constant.
inline constexpr double kHC_K = 0.014388;

// sRGB (IEC 61966-2-1) opto-electronic transfer function. Basically the
// gamma curve that turns linear light values into the sRGB values a monitor
// actually expects.
inline double srgbGamma(double c) {
    return c <= 0.0031308 ? 12.92 * c : 1.055 * std::pow(c, 1.0 / 2.4) - 0.055;
}

// Works out the physical blackbody colour at temperature `kelvin`.
// Returns a gamma-encoded, chromaticity-normalised display sRGB triple in [0,1].
inline std::array<double, 3> blackbodyRGB(double kelvin) {
    const double T = std::clamp(kelvin, 1000.0, 40000.0);

    double X = 0.0, Y = 0.0, Z = 0.0; 
    for (int j = 0; j < kCMFSamples; ++j) {
        const double lam = (380.0 + j * 10.0) * 1e-9;   // metres
        const double e   = kHC_K / (lam * T);
        // Planck spectral radiance, unnormalised for now. We guard the exp
        // call since it can blow up for very cool temperatures at short wavelengths.
        const double B = (e < 709.0)
                       ? 1.0 / (std::pow(lam, 5.0) * (std::exp(e) - 1.0))
                       : 0.0;
        X += B * kCMF_x[j];
        Y += B * kCMF_y[j];
        Z += B * kCMF_z[j];
    }

    // Convert XYZ to linear sRGB using the D65 primaries.
    double r =  3.2406 * X - 1.5372 * Y - 0.4986 * Z;
    double g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
    double b =  0.0557 * X - 0.2040 * Y + 1.0570 * Z;

    r = std::max(r, 0.0);
    g = std::max(g, 0.0);
    b = std::max(b, 0.0);

    // Since we only care about the hue here and not the absolute brightness,
    // we normalise so the brightest channel is 1. The shader's own flux and
    // beaming math handles how bright things actually get.
    const double mx = std::max({r, g, b});
    if (mx > 1e-12) { r /= mx; g /= mx; b /= mx; }

    return { srgbGamma(r), srgbGamma(g), srgbGamma(b) };
}

// Linear sRGB-primary blackbody radiance: NOT gamma-encoded and NOT
// chromaticity-normalised, so it carries absolute brightness as well as hue.
// Scaled so the luminance Y(6500 K) == 1. This is the quantity the Step 3
// "radiance" LUT stores (GL_RGB32F, log-spaced), and because the frequency
// shift g on a blackbody maps B(T) -> B(gT) exactly (I_nu/nu^3 is invariant),
// a single lookup at g*T gives colour, beaming and brightness together.
inline std::array<double, 3> blackbodyLinear(double kelvin) {
    auto cieXYZ = [](double T, double& X, double& Y, double& Z) {
        X = Y = Z = 0.0;
        for (int j = 0; j < kCMFSamples; ++j) {
            const double lam = (380.0 + j * 10.0) * 1e-9;   // metres
            const double e   = kHC_K / (lam * T);
            const double B   = (e < 709.0)
                             ? 1.0 / (std::pow(lam, 5.0) * std::expm1(e))
                             : 0.0;
            X += B * kCMF_x[j];
            Y += B * kCMF_y[j];
            Z += B * kCMF_z[j];
        }
    };
    // Reference luminance at 6500 K, computed once.
    static const double Yref = [&] {
        double X, Y, Z;
        cieXYZ(6500.0, X, Y, Z);
        return Y;
    }();

    double X, Y, Z;
    cieXYZ(std::max(kelvin, 1.0), X, Y, Z);
    X /= Yref; Y /= Yref; Z /= Yref;

    return { std::max(0.0,  3.2406 * X - 1.5372 * Y - 0.4986 * Z),
             std::max(0.0, -0.9689 * X + 1.8758 * Y + 0.0415 * Z),
             std::max(0.0,  0.0557 * X - 0.2040 * Y + 1.0570 * Z) };
}

/*───────Novikov-Thorne / Page & Thorne 1974)───────*/

// The three functions below are the standard relativistic Kerr 
// temperature properties/profile, referenced by
// BlackHole3D_PhotorealDisk.frag so the shader's disk colour is exactly what CI
// checks here.
//
// Description: this is the
// Novikov-Thorne disk.  That last part is the glaring assumption: real
// inner disks can carry magnetic stress across the ISCO, which lifts the edge
// temperature instead of pinning it to zero. We take the classical zero-torque
// edge because that's what "Novikov-Thorne" means. We also ignore returning
// radiation and self-irradiation. The profile is spin-dependent through both
// the ISCO location and the relativistic flux correction — that spin term is
// exactly what a plain Shakura-Sunyaev r^-3/4·[1-√(rISCO/r)]^1/4 profile drops,
// and dropping it was the old (mislabelled) behaviour.

inline constexpr double kPi = 3.14159265358979323846;

// Exact Kerr prograde ISCO radius in units of Rs (Bardeen, Press & Teukolsky
// 1972, ApJ 178, 347). a is the dimensionless spin a* in [0,1). Returns 3 Rs at
// a=0 (Schwarzschild) and shrinks toward 0.5 Rs as a -> 1 (prograde extremal).
inline double computeISCO(double a) {
    a = std::clamp(a, 0.0, 0.9999);
    const double a2 = a * a;
    const double Z1 = 1.0 + std::cbrt(1.0 - a2) *
                            (std::cbrt(1.0 + a) + std::cbrt(1.0 - a));
    const double Z2 = std::sqrt(3.0 * a2 + Z1 * Z1);
    const double rMsOverM = 3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2));
    return 0.5 * rMsOverM;   // Rs = 2M, so r/Rs = (r/M)/2
}

// Unnormalised Page-Thorne radiative flux shape F(r) for an equatorial Kerr
// disk (r in Rs, converted to M internally). This is the relativistic
// generalisation of the Newtonian [1-√(rISCO/r)]/r^3 law: it adds the
// spin-dependent correction that Shakura-Sunyaev omits. Valid for r > rISCO;
// returns 0 at or inside the ISCO. Page & Thorne 1974, ApJ 191, 499 (eq. 15n).
inline double pageThorneFlux(double rRs, double a) {
    a = std::clamp(a, 0.0, 0.9999);
    const double rIscoRs = computeISCO(a);
    if (rRs <= rIscoRs) return 0.0;

    // Work in M-units 
    // r_M = 2·r_Rs, x = √(r_M).
    const double x  = std::sqrt(2.0 * rRs);
    const double x0 = std::sqrt(2.0 * rIscoRs);

    // Roots of x^3 - 3x + 2a = 0.
    const double ac = std::acos(a);
    const double x1 =  2.0 * std::cos(ac / 3.0 - kPi / 3.0);
    const double x2 =  2.0 * std::cos(ac / 3.0 + kPi / 3.0);
    const double x3 = -2.0 * std::cos(ac / 3.0);

    // The x2 coefficient is a 0/0 as a -> 0 (x2 -> 0); the limit is 0.
    const double c1 = 3.0 * (x1 - a) * (x1 - a) / (x1 * (x1 - x2) * (x1 - x3));
    const double c2 = (std::abs(x2) < 1e-9)
                    ? 0.0
                    : 3.0 * (x2 - a) * (x2 - a) / (x2 * (x2 - x1) * (x2 - x3));
    const double c3 = 3.0 * (x3 - a) * (x3 - a) / (x3 * (x3 - x1) * (x3 - x2));

    const double Q = x - x0 - 1.5 * a * std::log(x / x0)
                   - c1 * std::log((x - x1) / (x0 - x1))
                   - c2 * std::log((x - x2) / (x0 - x2))
                   - c3 * std::log((x - x3) / (x0 - x3));

    const double denom = x * x * x * x * (x * x * x - 3.0 * x + 2.0 * a);
    const double F = Q / denom;
    return F > 0.0 ? F : 0.0;
}

// Peak of the unnormalised flux shape at spin a (the hottest annulus sits
// within ~1.5·rISCO, so a scan out to 3·rISCO comfortably brackets it).
inline double ntPeakFlux(double a) {
    const double rin = computeISCO(a);
    double best = 1e-300;
    for (int i = 0; i <= 256; ++i) {
        const double r = rin * (1.0 + 2.0 * double(i) / 256.0);
        best = std::max(best, pageThorneFlux(r, a));
    }
    return best;
}

// Normalised relativistic disk temperature in [0,1], peaking at 1.0 at the
// hottest annulus. T ∝ F^(1/4) (Stefan-Boltzmann). a is the spin a*.
inline double novikovThorneTemperature(double rRs, double a) {
    const double F = pageThorneFlux(rRs, a);
    if (F <= 0.0) return 0.0;
    const double t = std::pow(F / ntPeakFlux(a), 0.25);
    return std::clamp(t, 0.0, 1.0);
}

// Visible-band disk colour temperature at radius r, matching the physical
// disk-colour path in BlackHole3D_PhotorealDisk.frag:
//   T_colour(r) = displayTempInner × NT_profile(r; a)
// diskPeakTemp is left out on purpose here since that's the bolometric/UV
// effective temperature and it drives luminosity, not hue.
inline double diskColorTemperature(double displayTempInner, double rRs, double a) {
    return displayTempInner * novikovThorneTemperature(rRs, a);
}

// Exact frequency-shift factor g = nu_obs / nu_emit for gas on a prograde
// circular Kerr orbit at radius r (Rs), as seen by a distant static observer.
// This is the GL-free twin of diskRedshift() in the shader; keeping it here
// means the shader's disk colour and the CI checks sample the same formula.
//   a   : dimensionless spin a* (clamped to [0, 0.998]).
//   Lz  : conserved photon L/E of the traced ray (0 => pure gravitational +
//         transverse shift, i.e. the "Doppler off" case).
//   doppler : when false, Lz is ignored (showDoppler == 0 in the shader).
// Reference checkpoint: r=3, a=0, Lz=0 gives g = sqrt(1/2) = 0.70711.
inline double diskRedshift(double rRs, double a, double Lz, bool doppler) {
    a = std::clamp(a, 0.0, 0.998);
    const double rM = 2.0 * rRs;              // Rs -> M units
    const double s  = rM * std::sqrt(rM);     // rM^1.5
    const double Om = 1.0 / (s + a);          // gas angular velocity (M units)
    const double ut = (s + a) /
        (std::sqrt(s) * std::sqrt(std::max(s - 3.0 * std::sqrt(rM) + 2.0 * a, 1e-5)));
    const double LzM = doppler ? 2.0 * Lz : 0.0;
    return 1.0 / (ut * (1.0 + Om * LzM));
}

}  // namespace physics
}  // namespace bh3d

#endif  // BH3D_BLACKBODY_HPP

