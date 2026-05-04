/*------------------------------------ FITS Export Utilities ------------------------------------*/
#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include "research_data.hpp"

namespace FITSExport {

/*------------------------------------ FITS Export Utilities ------------------------------------*/

// FITS headers are 80-char cards, padded with spaces, in 2880-byte blocks
static void writeCard(std::ofstream& f, const std::string& card) {
    char buf[81] = {};
    std::memset(buf, ' ', 80);
    size_t len = std::min(card.size(), size_t(80));
    std::memcpy(buf, card.c_str(), len);
    f.write(buf, 80);
}

static void writePadBlock(std::ofstream& f, int cardsWritten) {
    // Pad header to 2880-byte boundary with spaces
    int remainder = cardsWritten % 36; // 36 cards per block
    if (remainder != 0) {
        char blank[80];
        std::memset(blank, ' ', 80);
        for (int i = remainder; i < 36; ++i)
            f.write(blank, 80);
    }
}

// Swap byte order of a double (for big-endian FITS format)
static double swapDouble(double v) {
    uint64_t tmp;
    std::memcpy(&tmp, &v, 8);
    tmp = ((tmp & 0x00000000000000FFULL) << 56) |
          ((tmp & 0x000000000000FF00ULL) << 40) | // This is the least significant byte, which goes to the most significant position
          ((tmp & 0x0000000000FF0000ULL) << 24) | // Please, PLEASE be careful with these bitmasks and shifts. It's easy to mess up the order.
          ((tmp & 0x00000000FF000000ULL) <<  8) |
          ((tmp & 0x000000FF00000000ULL) >>  8) |
          ((tmp & 0x0000FF0000000000ULL) >> 24) |
          ((tmp & 0x00FF000000000000ULL) >> 40) |
          ((tmp & 0xFF00000000000000ULL) >> 56);
    std::memcpy(&v, &tmp, 8); // copies the swapped bytes back to the original var (double)
    return v; 
}

static int32_t swapInt32(int32_t v) { // This is literally just the 32 bit equivalent of the above, but for 4 instead of 8 bytes. I could probably templatize this but... do I really want to?
    uint32_t tmp;
    std::memcpy(&tmp, &v, 4);
    tmp = ((tmp & 0x000000FF) << 24) | // Again, be very careful with these bitmasks and shifts. Triple-check them if you fork or make a pull request
          ((tmp & 0x0000FF00) <<  8) |
          ((tmp & 0x00FF0000) >>  8) |
          ((tmp & 0xFF000000) >> 24);
    std::memcpy(&v, &tmp, 4); // copies the swapped bytes back to the original var (int32_t)
    return v;
}

// Format a fixed-width FITS header card
static std::string keyVal(const std::string& key, int val, const std::string& comment = "") {
    char buf[81];
    std::memset(buf, ' ', 80); // Fill with spaces first
    buf[80] = '\0';            // Null terminate only for snprintf's safety
    
    if (comment.empty())
        snprintf(buf, 81, "%-8s= %20d", key.c_str(), val);
    else
        snprintf(buf, 81, "%-8s= %20d / %s", key.c_str(), val, comment.c_str());

    // bug found - turns out writing 81 null terminators for FITS headers causes... problems. Who would have thought? 
    // So we need to overwrite that null terminator (and any chars after it) back to spaces to ensure the entire 80-char card is valid.
    std::string result(buf, 80);
    for (size_t i = result.find('\0'); i < 80; ++i) {
        result[i] = ' ';
    }
    return result;
}
// Tested as of Wednesday, April 8, 2026 - this file functions mostly correctly.

// Overload for string values (e.g. EXTNAME), which need to be quoted in FITS headers
static std::string keyVal(const std::string& key, const std::string& val, const std::string& comment = "") {
    char buf[81] = {};
    std::string quoted = "'" + val + "'";
    if (comment.empty())
        snprintf(buf, 81, "%-8s= %-20s", key.c_str(), quoted.c_str());
    else
        snprintf(buf, 81, "%-8s= %-20s / %s", key.c_str(), quoted.c_str(), comment.c_str());
    // CRITICAL: snprintf put a \0 at the end of the text. 
    // We must overwrite that \0 (and everything after it) back to spaces.
    std::string result(buf, 80);
    for (size_t i = result.find('\0'); i < 80; ++i) {
        result[i] = ' ';
    }
    return result;
}

// Overload for double values — fixed-width FITS floating-point card.
// FITS uses ASCII E-notation, right-justified in columns 11–30.
static std::string keyVal(const std::string& key, double val, const std::string& comment = "") {
    char buf[81] = {};
    if (comment.empty())
        snprintf(buf, 81, "%-8s= %20.12E", key.c_str(), val);
    else
        snprintf(buf, 81, "%-8s= %20.12E / %s", key.c_str(), val, comment.c_str());
    std::string result(buf, 80);
    for (size_t i = result.find('\0'); i < 80; ++i) {
        result[i] = ' ';
    }
    return result;
}

/*--------- Simulation metadata baked into PRIMARY HDU ---------*/
// Research-grade FITS files need the simulation parameters that produced
// the data attached to the primary HDU header so downstream analysis tools
// can reconstruct units and physical context without out-of-band metadata.
struct SimulationMetadata {
    double M_BH      = 0.0;   // black hole geometric mass [M_sun]
    double c_s       = 0.0;   // sound speed [units of c]
    double r_Bondi   = 0.0;   // Bondi radius [units of M = r_g/2 ; equivalently the value reported by Schwarzschild::bondiRadius / M]
    double beta      = 0.0;   // dimensionless v/c of the exported body (orbital velocity)
    bool   valid     = false; // when false, only the standard SIMPLE/BITPIX/NAXIS/EXTEND cards are written
};

/*--------- Primary HDU (mandatory empty image for pure-table FITS) ---------*/
static int writePrimaryHDU(std::ofstream& f, const SimulationMetadata& meta = {}) {
    int n = 0;
    writeCard(f, "SIMPLE  =                    T / Aetherion Astrophysics Suite"); n++;
    writeCard(f, "BITPIX  =                    8"); n++;
    writeCard(f, "NAXIS   =                    0"); n++;
    writeCard(f, "EXTEND  =                    T"); n++;

    // Bake simulation parameters into the primary header (HIERARCH-style
    // 8-char keys keep the standard FITS card layout). These are the four
    // core Bondi-accretion parameters required for research reproducibility.
    if (meta.valid) {
        writeCard(f, keyVal("M_BH",    meta.M_BH,    "[Msun] black hole mass"));                 n++;
        writeCard(f, keyVal("C_S",     meta.c_s,     "[c] gas sound speed"));                    n++;
        writeCard(f, keyVal("R_BONDI", meta.r_Bondi, "[M] Bondi accretion radius"));             n++;
        writeCard(f, keyVal("BETA",    meta.beta,    "[c] orbital velocity v/c of exported body")); n++;
        writeCard(f, "COMMENT Parameters above describe the simulation that produced this HDU.");  n++;
        writeCard(f, "COMMENT Units: M_BH solar masses; C_S, BETA in units of c; R_BONDI in M.");  n++;
    }

    writeCard(f, "END"); n++;
    writePadBlock(f, n);
    return n;
}

/*------------------------------------  Binary Table HDU Writer ------------------------------------*/
// colNames:   column labels
// colFormats: FITS TFORM strings e.g. "1D" (double), "1J" (int32), "1L" (logical)
// colUnits:   physical units per column (empty string = none)
// nRows:      number of data rows to write
// writeData:  lambda/callback that writes the binary row data
template<typename DataWriter> // will this cause a memory leak? Who knows! it's a gamble :)
static void writeBinTableHDU(
    std::ofstream& f,
    const std::string& extname,
    const std::vector<std::string>& colNames,
    const std::vector<std::string>& colFormats,
    const std::vector<std::string>& colUnits,
    int nRows, 
    int rowBytes,
    DataWriter writeData)
{
    int nCols = (int)colNames.size();

    // Build header cards
    std::vector<std::string> cards;
    cards.push_back("XTENSION= 'BINTABLE'           / Binary table extension");
    cards.push_back(keyVal("BITPIX",   8));
    cards.push_back(keyVal("NAXIS",    2));
    cards.push_back(keyVal("NAXIS1",   rowBytes, "bytes per row"));
    cards.push_back(keyVal("NAXIS2",   nRows,    "number of rows"));
    cards.push_back(keyVal("PCOUNT",   0));
    cards.push_back(keyVal("GCOUNT",   1));
    cards.push_back(keyVal("TFIELDS",  nCols,    "number of columns"));
    cards.push_back(keyVal("EXTNAME",  extname));

    for (int i = 0; i < nCols; ++i) {
        char key[16];
        snprintf(key, 16, "TTYPE%d", i + 1);
        cards.push_back(keyVal(key, colNames[i]));
        snprintf(key, 16, "TFORM%d", i + 1);
        cards.push_back(keyVal(key, colFormats[i]));
        if (!colUnits[i].empty()) {
            snprintf(key, 16, "TUNIT%d", i + 1);
            cards.push_back(keyVal(key, colUnits[i]));
        }
    }
    cards.push_back("END");

    for (auto& c : cards) {
        // Pad each card to exactly 80 chars
        c.resize(80, ' ');
        f.write(c.c_str(), 80);
    }
    writePadBlock(f, (int)cards.size());

    // Write binary data
    writeData(f);

    // Pads the data block to a 2880-byte boundary with zeros, otherwise the FITS reader is going to throw a fit (insert laughtrack here).
        long pos = (long)f.tellp();
    long rem = pos % 2880L;
    if (rem != 0) {
        long padSize = 2880L - rem;
        std::vector<char> pad(padSize, 0);
        f.write(pad.data(), padSize);
    }
    }

/*--------- Public export functions ---------*/

template<typename TrailContainer>
inline bool exportOrbitData(
    const std::string& filename,
    const TrailContainer& trail,
    const SimulationMetadata& meta = {})
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)trail.size();
    // 4 doubles per row = 32 bytes
    writeBinTableHDU(f, "ORBIT",
        {"WORLD_X", "WORLD_Y", "RADIUS", "PHI"},
        {"1D",      "1D",      "1D",     "1D"},
        {"rg",      "rg",      "rg",     "rad"},
        nRows, 32,
        [&](std::ofstream& out) {
            for (auto& [x, y] : trail) {
                double r = std::hypot(x, y);
                double p = std::atan2(y, x);
                double vals[4] = {x, y, r, p};
                for (double v : vals) out.write(
                    reinterpret_cast<char*>(&(v = swapDouble(v))), 8);
            }
        });
    return true;
}

inline bool exportDeflectionTable(
    const std::string& filename,
    const std::vector<PhotonDeflection>& table,
    const SimulationMetadata& meta = {})
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)table.size();
    // 3 doubles + 1 int32 = 28 bytes per row
    writeBinTableHDU(f, "DEFLECTION",
        {"IMPACT_B",         "DEFL_RAD",      "DEFL_DEG",      "CAPTURED"},
        {"1D",                "1D",             "1D",             "1J"},
        {"rg",                "rad",            "deg",            ""},
        nRows, 28,
        [&](std::ofstream& out) {
            for (auto& d : table) {
                double b   = swapDouble(d.impactParameter);
                double dr  = swapDouble(d.deflectionAngle);
                double ddg = swapDouble(d.deflectionDeg);
                int32_t cap = swapInt32(d.captured ? 1 : 0);
                out.write(reinterpret_cast<char*>(&b),   8);
                out.write(reinterpret_cast<char*>(&dr),  8);
                out.write(reinterpret_cast<char*>(&ddg), 8);
                out.write(reinterpret_cast<char*>(&cap), 4);
            }
        });
    return true;
}

inline bool exportConservationHistory(
    const std::string& filename,
    const ConservationTracker& tracker,
    const SimulationMetadata& meta = {})
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)tracker.driftHistory.size();
    // 2 doubles = 16 bytes
    writeBinTableHDU(f, "CONSERVATION",
        {"PROPER_TIME", "ENERGY_DRIFT"},
        {"1D",          "1D"},
        {"M",           ""},
        nRows, 16,
        [&](std::ofstream& out) {
            for (auto& [t, d] : tracker.driftHistory) {
                double tv = swapDouble(t);
                double dv = swapDouble(d);
                out.write(reinterpret_cast<char*>(&tv), 8);
                out.write(reinterpret_cast<char*>(&dv), 8);
            }
        });
    return true;
}

inline bool exportPrecessionData( // ok so this is just the same as the above but with different columns. I could probably templatize this whole thing but... do I really want to?
    const std::string& filename,
    const PrecessionTracker& tracker,
    const SimulationMetadata& meta = {})
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta); // writes the MPHDU with the standard header cards. This is not the function you're looking for.

    int nRows = (int)tracker.precessionPerOrbit.size(); // number of orbits tracked, which determines the number of rows in the table. This is just the size of the precessionPerOrbit vector, which should have one entry per orbit completed.
    // 1 int32 + 3 doubles = 28 bytes
    writeBinTableHDU(f, "PRECESSION",
        {"ORBIT_NUM", "PERIAPSIS_PHI", "PRECESS_RAD", "PRECESS_DEG"},
        {"1J",        "1D",             "1D",           "1D"},
        {"",          "rad",            "rad",           "deg"},
        nRows, 28,
        [&](std::ofstream& out) {
            for (size_t i = 0; i < tracker.precessionPerOrbit.size(); ++i) {
                int32_t orb = swapInt32((int32_t)(i + 1));
                double phi  = swapDouble( // the phi of the periapsis for this orbit, which is the angle at which the orbit comes closest to the black hole. This is tracked in the PrecessionTracker and updated each time an orbit completes. We need to swap it to big-endian because FITS is a big baby.
                    (i + 1 < tracker.periapsisAngles.size())
                    ? tracker.periapsisAngles[i + 1] : 0.0);
                double pr   = swapDouble(tracker.precessionPerOrbit[i]);
                double pd   = swapDouble(tracker.precessionPerOrbit[i] * 180.0 / M_PI);
                out.write(reinterpret_cast<char*>(&orb), 4); // start of orbit data: write the orbit number as a 32-bit integer
                out.write(reinterpret_cast<char*>(&phi), 8);
                out.write(reinterpret_cast<char*>(&pr),  8);
                out.write(reinterpret_cast<char*>(&pd),  8);
            }
        });
    return true;
}

/*--------- Bondi-Hoyle accretion history (Ṁ_Bondi, Ṁ_BHL, L_bol) ---------*/
inline bool exportAccretionHistory(
    const std::string& filename,
    const BondiAccretionTracker& tracker,
    const SimulationMetadata& meta = {})
{
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)tracker.history.size();
    // 4 doubles = 32 bytes per row: TIME [M], MDOT_BND [kg/s], MDOT_BHL [kg/s], L_BOL [W]
    writeBinTableHDU(f, "ACCRETION",
        {"TIME",   "MDOT_BND",  "MDOT_BHL",  "L_BOL"},
        {"1D",     "1D",         "1D",         "1D"},
        {"M",      "kg/s",       "kg/s",       "W"},
        nRows, 32,
        [&](std::ofstream& out) {
            for (const auto& e : tracker.history) {
                double t  = swapDouble(e.t_M);
                double m1 = swapDouble(e.mdotBondi_kgs);
                double m2 = swapDouble(e.mdotBHL_kgs);
                double L  = swapDouble(e.L_bol_W);
                out.write(reinterpret_cast<char*>(&t),  8);
                out.write(reinterpret_cast<char*>(&m1), 8);
                out.write(reinterpret_cast<char*>(&m2), 8);
                out.write(reinterpret_cast<char*>(&L),  8);
            }
        });
    return true;
}

/*--------- Gas radial profile snapshot (Bondi-limit ρ, v_rad, c_s vs r) ---------*/
inline bool exportGasRadialProfile(
    const std::string& filename,
    const GasRadialProfile& prof,
    const SimulationMetadata& meta = {})
{
    if (!prof.valid) return false;
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)GasRadialProfile::N_SAMPLES;
    // 4 doubles = 32 bytes per row: R [M], RHO [kg/m^3], V_RAD [m/s], C_S [m/s]
    writeBinTableHDU(f, "GAS_PROF",
        {"R",      "RHO",       "V_RAD",   "C_S"},
        {"1D",     "1D",         "1D",      "1D"},
        {"M",      "kg/m3",      "m/s",     "m/s"},
        nRows, 32,
        [&](std::ofstream& out) {
            for (int i = 0; i < nRows; ++i) {
                double r  = swapDouble(prof.r_M[i]);
                double rh = swapDouble(prof.rho_kgm3[i]);
                double vr = swapDouble(prof.vRad_ms[i]);
                double cs = swapDouble(prof.cs_ms[i]);
                out.write(reinterpret_cast<char*>(&r),  8);
                out.write(reinterpret_cast<char*>(&rh), 8);
                out.write(reinterpret_cast<char*>(&vr), 8);
                out.write(reinterpret_cast<char*>(&cs), 8);
            }
        });
    return true;
}

/*--------- Spectral energy distribution snapshot (multi-color disk) ---------*/
inline bool exportSED(
    const std::string& filename,
    const SEDSnapshot& sed,
    const SimulationMetadata& meta = {})
{
    if (!sed.valid) return false;
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    writePrimaryHDU(f, meta);

    int nRows = (int)SEDSnapshot::N_FREQ;
    // 3 doubles = 24 bytes per row: NU [Hz], LNU [W/Hz], NU_LNU [W]
    writeBinTableHDU(f, "SED",
        {"NU",     "LNU",        "NU_LNU"},
        {"1D",     "1D",          "1D"},
        {"Hz",     "W/Hz",        "W"},
        nRows, 24,
        [&](std::ofstream& out) {
            for (int i = 0; i < nRows; ++i) {
                double nu  = swapDouble(sed.nu_Hz[i]);
                double Lnu = swapDouble(sed.Lnu_WHz[i]);
                double nL  = swapDouble(sed.nuLnu_W[i]);
                out.write(reinterpret_cast<char*>(&nu),  8);
                out.write(reinterpret_cast<char*>(&Lnu), 8);
                out.write(reinterpret_cast<char*>(&nL),  8);
            }
        });
    return true;
}

} // namespace FITSExport


//   absolute cinema
//         \O/
//          |
//         / \ 