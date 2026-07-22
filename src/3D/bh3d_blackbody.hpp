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
//   computeISCO(a)               Kerr prograde ISCO fit (mirrors shader).
//   novikovThorneTemperature()   Normalised T(r) ∝ r^-3/4 profile (shader).
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
    0.000039, 0.000120, 0.000396, 0.001210, 0.004000, 0.011600, 0.023000,
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

// Works out the Kerr prograde ISCO radius (in Rs) for a given spin a*.
// It's a smooth fit: ISCO(0) = 3 (Schwarzschild), shrinking down to about 0.6
// near extremal spin. Mirrors the same formula used in the shader.
inline double computeISCO(double a) {
    return std::max(3.0 - 2.0 * a * (1.0 + 0.27 * a), 0.6);
}

// Normalised Novikov-Thorne thin-disk temperature profile.
// T(r) ∝ r^-3/4 · [1 - sqrt(rISCO/r)]^1/4, normalised so the peak comes out to
// 1.0 (that peak happens to land at r = (49/36)·rISCO). Comes out to 0 right
// at or inside the ISCO, since that's the zero-torque inner boundary.
inline double novikovThorneTemperature(double r, double rISCO) {
    const double x = r / rISCO;
    const double noTorque = std::max(1.0 - std::sqrt(1.0 / x), 0.0);
    const double profile  = std::pow(1.0 / x, 0.75) * std::pow(noTorque, 0.25);

    const double xPeak   = 49.0 / 36.0;
    const double peakVal = std::pow(1.0 / xPeak, 0.75) *
                           std::pow(1.0 - std::sqrt(1.0 / xPeak), 0.25);

    return profile / std::max(peakVal, 0.001);
}

// Visible-band disk colour temperature at radius r, matching the physical
// disk-colour path in BlackHole3D_PhotorealDisk.frag:
//   T_colour(r) = displayTempInner times NT_profile(r)
// diskPeakTemp is left out on purpose here since that's the bolometric/UV
// effective temperature and it drives luminosity, not hue.
inline double diskColorTemperature(double displayTempInner, double r, double rISCO) {
    return displayTempInner * novikovThorneTemperature(r, rISCO);
}

}  // namespace physics
}  // namespace bh3d

#endif  // BH3D_BLACKBODY_HPP
