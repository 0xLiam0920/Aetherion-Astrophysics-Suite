#pragma once
// ============================================================
// catalog.hpp
// PURPOSE: shared black hole catalog (module-agnostic).
// Physical form data only. 3D visual tuning lives in each entry's
// "visual3d" JSON blob, applied by bh3d_catalog_adapter.hpp;
// "research2d" is reserved for the 2D subprogram.
//
// This header is deliberately Qt-free, glm-free, and free of any
// 3D-module type (cfg::SimConfig, BlackHoleProfile, etc etc). It carries
// physical parameters plus opaque JSON blobs that each module's
// adapter interprets for itself.
// ============================================================

#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

namespace catalog {

// One orbiting member. Type is a string here; each module maps it
// onto its own enum (GalaxyBodyType / GalaxyBody3DType).
struct BodyEntry {
    std::string type = "Star";
    double semiMajorRs  = 10.0;
    double eccentricity = 0.0;
    double inclination  = 0.0;   // stored in radians; 2D program ignores this though.
    std::string label;           // optional display name
};

struct Entry {
    int         schema = 1;
    std::string name;
    std::string description;
    double      massSolar = 10.0;
    double      spin      = 0.0;      // dimensionless a*
    std::string learnMoreUrl;
    int         sortKey   = 1000;     // menu ordering (filesystem order is not stable)
    bool        isGalacticCenter  = false;
    bool        barycentric       = false;
    double      companionMassSolar = 0.0;
    std::vector<BodyEntry> bodies;
    nlohmann::json visual3d;          // passed through untouched
    nlohmann::json research2d;        // passed through untouched
    std::string sourceFile;           // diagnostics
};

struct LoadResult {
    std::vector<Entry>       entries;
    std::vector<std::string> errors;  // one string per unreadable file
    bool ok() const { return !entries.empty(); }
};

inline Entry parseEntry(const nlohmann::json& j) {
    Entry e;
    e.schema       = j.value("schema", 1);
    e.name         = j.at("name").get<std::string>();      // required — throws if absent
    e.massSolar    = j.at("massSolar").get<double>();      // required
    e.description  = j.value("description", "");
    e.spin         = j.value("spin", 0.0);
    e.learnMoreUrl = j.value("learnMoreUrl", "");
    e.sortKey      = j.value("sortKey", 1000);
    e.isGalacticCenter = j.value("isGalacticCenter", false);
    if (j.contains("binary")) {
        e.barycentric        = j["binary"].value("barycentric", false);
        e.companionMassSolar = j["binary"].value("companionMassSolar", 0.0);
    }
    for (const auto& jb : j.value("galaxyBodies", nlohmann::json::array())) {
        BodyEntry b;
        b.type         = jb.value("type", "Star");
        b.semiMajorRs  = jb.value("semiMajorRs", 10.0);
        b.eccentricity = jb.value("ecc", 0.0);
        b.inclination  = jb.value("inclination", 0.0);
        b.label        = jb.value("label", "");
        e.bodies.push_back(std::move(b));
    }
    e.visual3d   = j.value("visual3d",   nlohmann::json::object());
    e.research2d = j.value("research2d", nlohmann::json::object());
    return e;
}

// Load every *.json in `dir` (non-recursive). A bad file is skipped and
// reported; it never poisons the rest of the catalog.
inline LoadResult loadDir(const std::filesystem::path& dir) {
    LoadResult out;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return out;
    for (const auto& de : std::filesystem::directory_iterator(dir, ec)) {
        if (de.path().extension() != ".json") continue;
        std::ifstream f(de.path());
        try {
            // ignore_comments=true: user-edited files get to keep // notes
            auto j = nlohmann::json::parse(f, nullptr, true, true);
            Entry e = parseEntry(j);
            e.sourceFile = de.path().string();
            out.entries.push_back(std::move(e));
        } catch (const std::exception& ex) {
            out.errors.push_back(de.path().filename().string() + ": " + ex.what());
        }
    }
    std::sort(out.entries.begin(), out.entries.end(),
              [](const Entry& a, const Entry& b) {
                  return a.sortKey != b.sortKey ? a.sortKey < b.sortKey
                                                : a.name    < b.name;
              });
    return out;
}

} // namespace catalog
