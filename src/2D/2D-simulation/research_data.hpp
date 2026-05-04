/*
TODO: DWARF GALAXY & generic body SIMULATION
=====================================================

PHASE 1: CORE PHYSICS VALIDATION
--------------------------------
- Implement custom FITS unit parsing for 'rg', 'M' in astropy (for export compatibility)
- Export simulation parameters to PRIMARY HDU: M_BH, c_s, r_Bondi/rg, beta
- Add geodesic convergence tests (ISCO radius, frame-dragging angles)
- Compute/export ENERGY_DRIFT per timestep for all integrators

PHASE 2: DYNAMICAL DIAGNOSTICS  
-----------------------------
- Track stream morphology: WIDTH(rg), LENGTH(rg), DENSITY profile
- Loss-cone statistics: THETA_LC, FEED_RATE(stars/s), J_SCATTER
- Angular momentum histogram: J_MIN, J_ISCO per particle
- Precession module: PERIAPSIS_PHI evolution (fix empty table)

PHASE 3: ACCRETION & OBSERVABLES
-------------------------------
- Bondi-Hoyle accretion: MDOT_BONDI(t), MDOT_HOLE(t/M) 
- Gas dynamics grid: RHO_GAS(r,phi), VRAD(r), CS_PROFILE(r)
- Bolometric light curve: L_BOL(t), peak time/flux
- Spectral energy distribution mockup (SED)

PHASE 4: NUMERICAL CONVERGENCE  
-----------------------------
- Adaptive timestep convergence: DT_MIN history, CFL_NUMBER
- Spatial resolution sweep: N_PARTICLES=1e3→1e5
- Integrator comparison: RK4 vs GEOKON vs BULIRSCH-STOER  
- Particle noise analysis: shot noise vs dynamical friction

PHASE 5: PUBLICATION EXPORTS
---------------------------
- Full FITS suite with WCS headers (astropy-compliant)
- HDF5 particle dump for community codes (GADGET/AREPO)
- LaTeX table generator for paper (parameters, convergence)
- JWST/ELT mock images of stream (surface brightness)

SIMULATION BOUNDS (Sgr A* baseline):
- M_BH = 4e6 Msun → rg = 1.2e10 cm
- r_Bondi ≈ 1e4 rg (c_s=10 km/s)
- r_SOI ≈ 1e6 rg (3 pc)
- t_final = 10^5 M (~1 day)

Maybe a bug if it still exists: ENERGY_DRIFT < 1e-10, J_CONSERVATION < 1e-12


              _-o#&&*''''?d:>b\_
          _o/"`''  '',, dMF9MMMMMHo_
       .o&#'        `"MbHMMMMMMMMMMMHo.
     .o"" '         vodM*$&&HMMMMMMMMMM?.
    ,'              $M&ood,~'`(&##MMMMMMH\
   /               ,MMMMMMM#b?#bobMMMMHMMML
  &              ?MMMMMMMMMMMMMMMMM7MMM$R*Hk
 ?$.            :MMMMMMMMMMMMMMMMMMM/HMMM|`*L
|               |MMMMMMMMMMMMMMMMMMMMbMH'   T,
$H#:            `*MMMMMMMMMMMMMMMMMMMMb#}'  `?
]MMH#             ""*""""*#MMMMMMMMMMMM'     -
MMMMMb_                   |MMMMMMMMMMMP'     :
HMMMMMMMHo                 `MMMMMMMMMT       .
?MMMMMMMMP                  9MMMMMMMM}       -
-?MMMMMMM                  |MMMMMMMMM?,d-    '
 :|MMMMMM-                 `MMMMMMMT .M|.   :
  .9MMM[                    &MMMMM*' `'    .
   :9MMk                    `MMM#"        -
     &M}                     `          .-
      `&.                             .
        `~,   .                     ./
            . _                  .-
              '`--._,dd###pp=""'

*/

/*---------- Header files ---------*/
#pragma once
#include <vector>
#include <deque>
#include <array>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846 // Yes, this is a nitpick. We need 20 or so significant figures since precession and photon angles can be extremely small. something something radians as well
#endif

/*--------- Precession tracking for orbiting bodies ---------*/
struct PrecessionTracker {
    std::vector<double> periapsisAngles;
    std::vector<double> precessionPerOrbit;
    // 3-sample sliding window of r (and its phi) for local-minimum detection.
    double r_prev  = -1.0;
    double r_prev2 = -1.0;
    double phi_prev = 0.0;

    // periapsis detection via local minimum of r over a 3-sample window.
    // This is robust to inspiraling orbits where vr never flips sign positive
    // (GW decay, drag, capture trajectories) — the original vr-sign-flip
    // detector silently produced zero entries for those cases.
    // We register the periapsis at the *middle* sample's phi.
    void detectPeriapsis(double r, double /*vr*/, double phi) {
        if (r_prev > 0.0 && r_prev2 > 0.0) {
            // Require strict descent then ascent-or-equal across the window.
            // The relative tolerance suppresses jitter when r is near-stationary.
            const double tol = 1e-9 * r_prev;
            const bool isMin = (r_prev2 - r_prev) > tol && (r - r_prev) >= -tol;
            if (isMin) {
                periapsisAngles.push_back(phi_prev);
                if (periapsisAngles.size() >= 2) {
                    size_t n = periapsisAngles.size();
                    double delta = periapsisAngles[n-1] - periapsisAngles[n-2] - 2.0 * M_PI;
                    precessionPerOrbit.push_back(delta);
                }
            }
        }
        r_prev2 = r_prev;
        r_prev  = r;
        phi_prev = phi;
    }

    double lastPrecessionRad() const {
        if (precessionPerOrbit.empty()) return 0.0;
        return precessionPerOrbit.back();
    }

    double lastPrecessionDeg() const {
        return lastPrecessionRad() * 180.0 / M_PI;
    }

    double averagePrecessionDeg() const {
        if (precessionPerOrbit.empty()) return 0.0;
        double sum = 0.0;
        for (double p : precessionPerOrbit) sum += p;
        return (sum / (double)precessionPerOrbit.size()) * 180.0 / M_PI;
    }

    // "Orbits completed" = number of detected periapsis passages. The UI uses
    // this for live progress; precession entries are always one less than this
    // (we need two periapsides to measure one Δφ).
    int orbitsCompleted() const {
        return (int)periapsisAngles.size();
    }

    void reset() {
        periapsisAngles.clear();
        precessionPerOrbit.clear();
        r_prev  = -1.0;
        r_prev2 = -1.0;
        phi_prev = 0.0;
    }
};

/*--------- Energy conservation tracking ---------*/
struct ConservationTracker {
    double E_initial = 0.0;
    double E_current = 0.0;
    double maxAbsDrift = 0.0;
    std::deque<std::pair<double, double>> driftHistory; // deque for O(1) pop_front(); was std::vector [FIXED 2026-04-24]
    static constexpr size_t MAX_HISTORY = 2000;
    int sampleCounter = 0;

    // NOTE: E_initial is set at orbit init, not at first step — so any numerical transient
    // at the very first step counts against drift. this is intentional: we want to catch
    // initialization errors too, not just integration drift.
    void init(double E0) {
        E_initial = E0;
        E_current = E0;
        maxAbsDrift = 0.0;
        driftHistory.clear();
        sampleCounter = 0;
    }

    void update(double E_now, double properTime) {
        E_current = E_now;
        double drift = (E_initial > 1e-15)
            ? (E_now - E_initial) / E_initial : 0.0;
        maxAbsDrift = std::max(maxAbsDrift, std::abs(drift));

        // only record every 50 steps, at 4000 sub-steps/frame we'd generate thousands of
        // history points per second otherwise, which would fill MAX_HISTORY almost instantly.
        // 50 is a rough balance between resolution and not eating all your RAM.
        sampleCounter++;
        if (sampleCounter % 50 == 0) {
            driftHistory.emplace_back(properTime, drift);
            if (driftHistory.size() > MAX_HISTORY)      // [FIXED 2026-04-24: was erase(begin()), O(n) on vector — now pop_front() on deque, O(1)]
                driftHistory.pop_front();
        }
    }

    double relativeDrift() const {
        if (E_initial < 1e-15) return 0.0;
        return (E_current - E_initial) / E_initial;
    }

    void reset() {
        E_initial = E_current = 0.0;
        maxAbsDrift = 0.0;
        driftHistory.clear();
        sampleCounter = 0;
    }
};

/*--------- RK4 error estimation (step-doubling) ---------*/
// step-doubling: take one big step, then two half-steps, compare. error ≈ |big - double_half| / 15.
// the 15 is from Richardson extrapolation for a 4th-order method (2^4 - 1 = 15).
// we don't do this every step because it triples the computation cost.
struct NumericalErrorTracker {
    double lastError = 0.0;
    double maxError  = 0.0;
    double avgError  = 0.0;
    int stepCount = 0;
    int checkInterval = 100;
    int stepsSinceCheck = 0;

    void recordError(double error) {
        lastError = error;
        maxError = std::max(maxError, error);
        avgError = (avgError * stepCount + error) / (stepCount + 1);
        stepCount++;
    }

    bool shouldCheck() {
        stepsSinceCheck++;
        if (stepsSinceCheck >= checkInterval) {
            stepsSinceCheck = 0;
            return true;
        }
        return false;
    }

    void reset() {
        lastError = maxError = avgError = 0.0;
        stepCount = 0;
        stepsSinceCheck = 0;
    }
};

/*--------- Bondi-Hoyle accretion tracker ---------*/
// Implements the analytic Bondi-Hoyle-Lyttleton (BHL) spherical accretion model:
//
//     Ṁ_BHL = 4π λ G² M² ρ∞ / (v² + c_s²)^(3/2)        [kg/s]
//
// with λ = 0.25 for an adiabatic γ=5/3 gas (the textbook value; e.g. Bondi 1952,
// Edgar 2004 review). The "pure Bondi" limit (v=0) is reported separately so the
// effect of the accretor's relative motion through the gas is visible.
//
// Bolometric luminosity assumes thin-disk radiative efficiency at the
// Schwarzschild ISCO:
//
//     η = 1 - sqrt(8/9) ≈ 0.0572
//     L_bol = η Ṁ c²                                    [W]
//
// This is a strictly post-processing diagnostic — none of the gas is actually
// simulated, we just feed M_BH, c_s, and the selected body's orbital v through
// the analytic formula at every update and record the time series. It still
// gives a research-grade observable curve once the orbit has had time to evolve
// (the v-dependent (v²+c_s²)^(3/2) denominator makes Ṁ_BHL dip at periapsis).
struct BondiAccretionEntry {
    double t_M             = 0.0;   // simulation time [units of M]
    double mdotBondi_kgs   = 0.0;   // pure-Bondi rate (v=0)
    double mdotBHL_kgs     = 0.0;   // Bondi-Hoyle-Lyttleton rate (with v_rel)
    double L_bol_W         = 0.0;   // bolometric luminosity from MDOT_BHL [W]
};

struct BondiAccretionTracker {
    static constexpr double LAMBDA_ADIABATIC  = 0.25;        // λ for γ=5/3
    static constexpr double ETA_SCHWARZSCHILD = 0.0572;      // 1 - sqrt(8/9)
    static constexpr double G_SI    = 6.6743e-11;
    static constexpr double C_SI    = 2.99792458e8;
    static constexpr double MSUN_KG = 1.98892e30;
    static constexpr double YEAR_S  = 3.15576e7;
    static constexpr size_t MAX_HISTORY    = 2000;
    static constexpr int    SAMPLE_INTERVAL = 4;             // record every Nth update tick

    std::deque<BondiAccretionEntry> history;
    int sampleCounter = 0;

    // Last evaluated values — exposed for the data panel.
    double last_mdotBondi_kgs    = 0.0;
    double last_mdotBHL_kgs      = 0.0;
    double last_Lbol_W           = 0.0;
    double last_mdotBHL_MsunYr   = 0.0;   // convenience: BHL rate in M_sun/yr (the unit astronomers actually quote)
    double last_Lbol_Lsun        = 0.0;   // convenience: L_bol in solar luminosities (3.828e26 W)

    // Recompute Ṁ_Bondi, Ṁ_BHL, L_bol from current sim state.
    //   t_M        : simulation time in geometric units of M (any monotonic clock)
    //   M_BH_solar : black hole mass [M_sun]
    //   cs_over_c  : sound speed [c]
    //   rho_kgm3   : ambient gas density [kg/m^3]
    //   v_over_c   : relative velocity of accretor through the gas [c]
    void update(double t_M, double M_BH_solar, double cs_over_c,
                double rho_kgm3, double v_over_c) {
        if (M_BH_solar <= 0.0 || cs_over_c <= 0.0 || rho_kgm3 <= 0.0) return;

        const double M_kg = M_BH_solar * MSUN_KG;
        const double cs   = cs_over_c * C_SI;
        const double v    = v_over_c  * C_SI;

        // 4π λ G² M² ρ∞ — common prefactor. Pulling this out of the 1/denom
        // form keeps the two rates differing only in the denominator, which
        // is what we want for a clean BHL-vs-Bondi comparison.
        const double prefactor = 4.0 * M_PI * LAMBDA_ADIABATIC
                               * G_SI * G_SI * M_kg * M_kg * rho_kgm3;

        const double cs3      = cs * cs * cs;
        const double mdotBondi = (cs3 > 0.0) ? prefactor / cs3 : 0.0;

        const double denom    = std::pow(v * v + cs * cs, 1.5);
        const double mdotBHL  = (denom > 0.0) ? prefactor / denom : 0.0;

        // L_bol = η Ṁ c² — Schwarzschild ISCO efficiency.
        const double L_bol_W  = ETA_SCHWARZSCHILD * mdotBHL * C_SI * C_SI;

        last_mdotBondi_kgs   = mdotBondi;
        last_mdotBHL_kgs     = mdotBHL;
        last_Lbol_W          = L_bol_W;
        last_mdotBHL_MsunYr  = mdotBHL * YEAR_S / MSUN_KG;
        last_Lbol_Lsun       = L_bol_W / 3.828e26;

        // Sub-sample history to bound memory. At ~60 update calls/sec this gives
        // O(15) recorded entries per second, capped at MAX_HISTORY ≈ 2.2 minutes
        // of wall-clock playback before the oldest entries get popped.
        sampleCounter++;
        if (sampleCounter % SAMPLE_INTERVAL == 0) {
            history.push_back({t_M, mdotBondi, mdotBHL, L_bol_W});
            if (history.size() > MAX_HISTORY) history.pop_front();
        }
    }

    void reset() {
        history.clear();
        sampleCounter = 0;
        last_mdotBondi_kgs = last_mdotBHL_kgs = last_Lbol_W = 0.0;
        last_mdotBHL_MsunYr = last_Lbol_Lsun = 0.0;
    }
};

/*--------- Gas radial profile (analytic Bondi-limit snapshot) ---------*/
//
// Pure post-processing diagnostic: there is no real hydro in the simulation,
// so we evaluate the steady-state spherical Bondi γ=5/3 solution at a set of
// radii and export it as a snapshot table. Inside the Bondi radius
// r_B = GM / c_s∞² the flow becomes free-fall:
//     v_rad(r) ≈ -√(2GM/r)
//     ρ(r)    ≈ ρ_∞ (r_B / r)^(3/2)         (continuity + free-fall)
//     c_s(r)  ≈ c_s∞ √(r_B / r)             (adiabatic γ=5/3)
// Outside r_B the flow asymptotes to the unperturbed ISM (ρ_∞, c_s∞, v≈0).
//
// snapshot() fills the arrays from current simulation state. The struct holds
// no time history — it's a one-shot picture rebuilt on each export.
struct GasRadialProfile {
    static constexpr int    N_SAMPLES = 96;
    static constexpr double R_MIN_M   = 2.0;        // inside ISCO (in M units)
    static constexpr double R_MAX_M   = 1.0e5;      // far past r_Bondi for typical Sgr A*
    static constexpr double G_SI      = 6.6743e-11;
    static constexpr double C_SI      = 2.99792458e8;
    static constexpr double MSUN_KG   = 1.98892e30;

    std::array<double, N_SAMPLES> r_M{};            // sample radii [M]
    std::array<double, N_SAMPLES> rho_kgm3{};       // density       [kg/m³]
    std::array<double, N_SAMPLES> vRad_ms{};        // radial vel    [m/s] (negative = inflow)
    std::array<double, N_SAMPLES> cs_ms{};          // sound speed   [m/s]
    bool   valid       = false;
    double M_BH_solar  = 0.0;
    double rho_inf     = 0.0;
    double cs_inf      = 0.0;
    double r_Bondi_M   = 0.0;

    // Compute snapshot. r_g = GM/c² is the geometric length unit for r_M[].
    void snapshot(double M_BH_solar_, double rho_inf_kgm3, double cs_inf_ms) {
        if (M_BH_solar_ <= 0.0 || rho_inf_kgm3 <= 0.0 || cs_inf_ms <= 0.0) {
            valid = false;
            return;
        }
        M_BH_solar = M_BH_solar_;
        rho_inf    = rho_inf_kgm3;
        cs_inf     = cs_inf_ms;

        const double M_kg = M_BH_solar * MSUN_KG;
        const double GM   = G_SI * M_kg;
        const double rg_m = GM / (C_SI * C_SI);     // 1 M in metres
        const double rB_m = GM / (cs_inf * cs_inf); // Bondi radius in metres
        r_Bondi_M = rB_m / rg_m;

        // Log-spaced sample grid in r/M.
        const double logMin = std::log10(R_MIN_M);
        const double logMax = std::log10(R_MAX_M);
        for (int i = 0; i < N_SAMPLES; ++i) {
            const double frac = (N_SAMPLES > 1) ? double(i) / (N_SAMPLES - 1) : 0.0;
            const double rM   = std::pow(10.0, logMin + frac * (logMax - logMin));
            const double r_m  = rM * rg_m;
            r_M[i] = rM;

            // Smooth interpolation from free-fall (r ≪ r_B) to ISM (r ≫ r_B).
            // Use blending factor ξ = r_B / (r_B + r): ξ→1 at r→0, ξ→0 at r→∞.
            const double xi = rB_m / (rB_m + r_m);
            const double xi32 = std::pow(xi, 1.5);

            // Density: ρ_∞ × (1 + (r_B/r)^(3/2)) gives the right small-r limit
            // ρ ∝ r^(-3/2) and the right far-field limit ρ → ρ_∞.
            const double r_ratio = rB_m / r_m;
            const double rho     = rho_inf * (1.0 + std::pow(r_ratio, 1.5));
            rho_kgm3[i] = rho;

            // Free-fall velocity, smoothly damped to zero outside r_B.
            const double vff = std::sqrt(2.0 * GM / r_m);
            vRad_ms[i] = -vff * xi32;

            // Adiabatic c_s: c_s² ∝ ρ^(γ-1) with γ=5/3 → c_s ∝ ρ^(1/3)
            cs_ms[i] = cs_inf * std::cbrt(rho / rho_inf);
        }
        valid = true;
    }

    void reset() { valid = false; }
};

/*--------- Spectral energy distribution (multi-color disk) ---------*/
//
// Standard Shakura-Sunyaev geometrically-thin optically-thick disk:
//     T(r) = [ 3 G M Ṁ / (8 π σ_SB r³) × (1 - √(r_in/r)) ]^(1/4)
// Each annulus radiates as a blackbody at T(r); total ν L_ν is the
// integral of 4π² r B_ν(T(r)) dr from r_in (ISCO = 6M) to r_out.
//
// This is a snapshot — recomputed when exported. A toy SED in W per ν.
struct SEDSnapshot {
    static constexpr int    N_FREQ      = 200;
    static constexpr double NU_MIN_HZ   = 1.0e10;     // 10 GHz (radio)
    static constexpr double NU_MAX_HZ   = 1.0e20;     // 100 EHz (gamma)
    static constexpr int    N_RADIAL    = 256;
    static constexpr double R_OUT_FACTOR = 1000.0;    // r_out = 1000 × r_in
    static constexpr double G_SI        = 6.6743e-11;
    static constexpr double C_SI        = 2.99792458e8;
    static constexpr double H_PLANCK    = 6.62607015e-34;
    static constexpr double K_BOLTZ     = 1.380649e-23;
    static constexpr double SIGMA_SB    = 5.670374419e-8;
    static constexpr double MSUN_KG     = 1.98892e30;

    std::array<double, N_FREQ> nu_Hz{};
    std::array<double, N_FREQ> nuLnu_W{};      // ν L_ν  [W]
    std::array<double, N_FREQ> Lnu_WHz{};      // L_ν    [W/Hz]
    bool   valid     = false;
    double T_peak_K  = 0.0;     // hottest annulus temperature
    double L_int_W   = 0.0;     // ∫ L_ν dν  (sanity check vs. η Ṁ c²)

    // Compute SED snapshot. mdot_kgs is the BHL rate (we feed the hole at the
    // current orbit-modulated value for an "instantaneous" SED).
    void snapshot(double M_BH_solar, double mdot_kgs) {
        if (M_BH_solar <= 0.0 || mdot_kgs <= 0.0) {
            valid = false;
            return;
        }
        const double M_kg = M_BH_solar * MSUN_KG;
        const double GM   = G_SI * M_kg;
        const double r_in = 6.0 * GM / (C_SI * C_SI);             // ISCO (Schwarzschild)
        const double r_out = R_OUT_FACTOR * r_in;

        // Pre-tabulate T(r) on a log grid.
        std::array<double, N_RADIAL> r_m{}, T_K{};
        const double logRin  = std::log10(r_in);
        const double logRout = std::log10(r_out);
        T_peak_K = 0.0;
        for (int j = 0; j < N_RADIAL; ++j) {
            const double frac = double(j) / (N_RADIAL - 1);
            r_m[j] = std::pow(10.0, logRin + frac * (logRout - logRin));
            const double r3 = r_m[j] * r_m[j] * r_m[j];
            const double bracket = (3.0 * GM * mdot_kgs) / (8.0 * M_PI * SIGMA_SB * r3)
                                 * (1.0 - std::sqrt(r_in / r_m[j]));
            T_K[j] = (bracket > 0.0) ? std::pow(bracket, 0.25) : 0.0;
            if (T_K[j] > T_peak_K) T_peak_K = T_K[j];
        }

        // Integrate B_ν(T(r)) × 2πr over both faces (factor 2) at each ν.
        const double logNuMin = std::log10(NU_MIN_HZ);
        const double logNuMax = std::log10(NU_MAX_HZ);
        L_int_W = 0.0;
        double prev_nuLnu = 0.0, prev_nu = 0.0;
        for (int i = 0; i < N_FREQ; ++i) {
            const double frac = double(i) / (N_FREQ - 1);
            const double nu = std::pow(10.0, logNuMin + frac * (logNuMax - logNuMin));
            nu_Hz[i] = nu;

            // L_ν = 4π² ∫_{r_in}^{r_out} r B_ν(T(r)) dr  (both faces, isotropic emission per face)
            // Trapezoidal in r, with dr varying (log grid).
            double Lnu = 0.0;
            for (int j = 0; j < N_RADIAL - 1; ++j) {
                if (T_K[j] <= 0.0 && T_K[j+1] <= 0.0) continue;
                auto Bnu = [&](double T) -> double {
                    if (T <= 0.0) return 0.0;
                    const double x = (H_PLANCK * nu) / (K_BOLTZ * T);
                    if (x > 700.0) return 0.0;          // exp overflow guard
                    const double denom = std::expm1(x);  // exp(x) - 1
                    if (denom <= 0.0) return 0.0;
                    return (2.0 * H_PLANCK * nu * nu * nu) / (C_SI * C_SI * denom);
                };
                const double r0 = r_m[j],   r1 = r_m[j+1];
                const double f0 = r0 * Bnu(T_K[j]);
                const double f1 = r1 * Bnu(T_K[j+1]);
                Lnu += 0.5 * (f0 + f1) * (r1 - r0);
            }
            Lnu *= 4.0 * M_PI * M_PI;
            Lnu_WHz[i] = Lnu;
            const double nuLnu = nu * Lnu;
            nuLnu_W[i] = nuLnu;

            // Trapezoidal integration of L_ν over ν for sanity check.
            if (i > 0) L_int_W += 0.5 * (Lnu + prev_nuLnu / std::max(prev_nu, 1.0)) * (nu - prev_nu);
            prev_nuLnu = nuLnu;
            prev_nu = nu;
        }
        valid = true;
    }

    void reset() { valid = false; }
};

/*--------- Photon deflection data ---------*/
struct PhotonDeflection {
    double impactParameter = 0.0;
    double deflectionAngle = 0.0;
    double deflectionDeg   = 0.0;
    bool   captured = false;
};

/*--------- Lensing analysis ---------*/
struct LensingData {
    double criticalImpactParam = 0.0;
    std::vector<PhotonDeflection> deflectionTable;
    std::vector<std::pair<float, float>> causticPoints;
};

/*--------- Photon sphere test ---------*/
struct PhotonSphereTestResult {
    double impactParameter    = 0.0;
    double orbitsBeforeDecay  = 0.0;
    bool   escaped  = false;
    bool   captured = false;
    double stabilityAngle = 0.0;   // total φ swept before escape or capture
    // NOTE: for b exactly at b_crit, the theoretical stabilityAngle is infinite (photon orbits forever).
    // in practice we always get a finite number because floating point isn't exact and the
    // unstable orbit eventually tips one way or the other. the nearer to b_crit, the bigger the number.
};

/*--------- ISCO test ---------*/
struct ISCOTestResult {
    double testRadius_M = 0.0;
    std::string classification;     // "Stable", "Critical", "Unstable"
    double radiusDrift = 0.0;       // how much r has drifted from starting radius (in units of M)
    bool   captured = false;
    double survivalTime = 0.0;
    // if radiusDrift is growing monotonically — the orbit is unstable, as expected.
    // if it oscillates around zero — stable orbit, just normal orbital motion.
    // in practice even the "stable" 7M orbit shows a tiny secular drift due to RK4 not being
    // a symplectic integrator. it's small enough that it doesn't matter for demo purposes, but
    // if you run it for long enough it will eventually drift. this is a known RK4 limitation.
};

/*--------- Data panel display info ---------*/
struct DataPanelInfo {
    // Body state
    double radius_M     = 0.0;
    double velocity_c   = 0.0;
    double timeDilation  = 0.0;
    double redshift_z    = 0.0;
    double energy_E      = 0.0;
    double angularMomentum_L = 0.0;
    double properTime    = 0.0;
    double coordinateTime = 0.0;
    double escapeVelocity_c = 0.0;
    bool   isBound = true;
    std::string stabilityClass;

    // Tidal
    double tidalForce    = 0.0;
    double tidalStress1m = 0.0;

    // Conservation / error
    double energyDrift    = 0.0;
    double maxEnergyDrift = 0.0;
    double rk4LastError = 0.0;
    double rk4MaxError  = 0.0;
    double rk4AvgError  = 0.0;

    // Precession
    double precessionDeg = 0.0;
    double theoreticalPrecessionDeg = 0.0;
    int    orbitsCompleted = 0;

    // Bondi
    double bondiRadius_M     = 0.0;
    double gasTemperatureK   = 1e7;
    double soundSpeed_c      = 0.0;
    double gasDensity_kgm3   = 0.0;       // ambient gas density used for Bondi rate
    double mdotBondi_MsunYr  = 0.0;       // pure-Bondi (v=0) rate [M_sun/yr]
    double mdotBHL_MsunYr    = 0.0;       // Bondi-Hoyle-Lyttleton rate [M_sun/yr]
    double Lbol_W            = 0.0;       // bolometric luminosity [W]
    double Lbol_Lsun         = 0.0;       // bolometric luminosity [L_sun]

    // Lensing
    double photonImpactParam = 0.0;
    double photonDeflectionDeg = 0.0;
    // ok this should be enough. Lensing is a prototype for now, we can always add more fields later
};
