#pragma once
// ============================================================
// bh3d_catalog_adapter.hpp.
// PURPOSE:Converts shared catalog entries into the 3D module's
// BlackHoleProfile, applying the "visual3d" section as sparse
// overrides on cfg::defaultConfig(). Absent keys keep the
// compiled-in default, exactly like the C++ profile functions
// only ever set deltas.
// ============================================================

#include "../common/catalog.hpp"
#include "bh3d_presets.hpp"     // BlackHoleProfile, GalaxyBody3D(Type) (global namespace)
#include "config.hpp"
#include "resource_manager.hpp"
#include <iostream>
#include <vector>

namespace bh3d_catalog {
using nlohmann::json;

inline GalaxyBody3DType bodyType(const std::string& s) {
    if (s == "GasCloud")       return GalaxyBody3DType::GasCloud;
    if (s == "StellarCluster") return GalaxyBody3DType::StellarCluster;
    if (s == "DwarfGalaxy")    return GalaxyBody3DType::DwarfGalaxy;
    if (s == "NeutronStar")    return GalaxyBody3DType::NeutronStar;
    if (s == "WhiteDwarf")     return GalaxyBody3DType::WhiteDwarf;
    if (s == "CompanionStar")  return GalaxyBody3DType::CompanionStar;
    return GalaxyBody3DType::Star;
}

inline glm::vec3 vec3Of(const json& j, glm::vec3 cur) {
    return (j.is_array() && j.size() == 3)
        ? glm::vec3(j[0].get<float>(), j[1].get<float>(), j[2].get<float>())
        : cur;
}

inline void applyDisk(const json& j, cfg::DiskConfig& d) {
    d.innerRadius          = j.value("innerRadius",          d.innerRadius);
    d.outerRadius          = j.value("outerRadius",          d.outerRadius);
    d.halfThickness        = j.value("halfThickness",        d.halfThickness);
    d.peakTemp             = j.value("peakTemp",             d.peakTemp);
    d.displayTempInner     = j.value("displayTempInner",     d.displayTempInner);
    d.displayTempOuter     = j.value("displayTempOuter",     d.displayTempOuter);
    d.saturationBoostInner = j.value("saturationBoostInner", d.saturationBoostInner);
    d.saturationBoostOuter = j.value("saturationBoostOuter", d.saturationBoostOuter);
    d.physicalDiskColor    = j.value("physicalDiskColor",    d.physicalDiskColor);
}

inline void applyJet(const json& j, cfg::JetConfig& c) {
    c.radius = j.value("radius", c.radius);
    c.length = j.value("length", c.length);
    if (j.contains("color")) c.color = vec3Of(j["color"], c.color);
}

inline void applyBLR(const json& j, cfg::BLRConfig& c) {
    c.innerRadius = j.value("innerRadius", c.innerRadius);
    c.outerRadius = j.value("outerRadius", c.outerRadius);
    c.thickness   = j.value("thickness",   c.thickness);
    c.strength    = j.value("strength",    c.strength);
}

inline void applyOrbital(const json& j, cfg::OrbitalConfig& c) {
    c.semiMajor    = j.value("semiMajor",    c.semiMajor);
    c.eccentricity = j.value("eccentricity", c.eccentricity);
    c.inclination  = j.value("inclination",  c.inclination);
    c.bodyRadius   = j.value("bodyRadius",   c.bodyRadius);
    if (j.contains("bodyColor")) c.bodyColor = vec3Of(j["bodyColor"], c.bodyColor);
}

inline void applyBloom(const json& j, cfg::BloomConfig& c) {
    c.threshold = j.value("threshold", c.threshold);
    c.softKnee  = j.value("softKnee",  c.softKnee);
    c.intensity = j.value("intensity", c.intensity);
    c.exposure  = j.value("exposure",  c.exposure);
    if (j.contains("mipWeights") && j["mipWeights"].is_array() &&
        j["mipWeights"].size() == 4) {
        for (int i = 0; i < 4; ++i) c.mipWeights[i] = j["mipWeights"][i].get<float>();
    }
}

inline void applyCamera(const json& j, cfg::CameraConfig& c) {
    if (j.contains("initialPos")) c.initialPos = vec3Of(j["initialPos"], c.initialPos);
    c.initialYaw   = j.value("initialYaw",   c.initialYaw);
    c.initialPitch = j.value("initialPitch", c.initialPitch);
    c.fov          = j.value("fov",          c.fov);
}

inline BlackHoleProfile toProfile(const catalog::Entry& e) {
    BlackHoleProfile p;
    p.name        = e.name;
    p.description = e.description;
    p.massSolar   = e.massSolar;
    p.config      = cfg::defaultConfig();
    p.config.blackHole.spinParameter = (float)e.spin;

    const json& v = e.visual3d;
    if (v.contains("disk"))    applyDisk   (v["disk"],    p.config.disk);
    if (v.contains("jet"))     applyJet    (v["jet"],     p.config.jet);
    if (v.contains("blr"))     applyBLR    (v["blr"],     p.config.blr);
    if (v.contains("orbital")) applyOrbital(v["orbital"], p.config.orbital);
    if (v.contains("bloom"))   applyBloom  (v["bloom"],   p.config.bloom);
    if (v.contains("camera"))  applyCamera (v["camera"],  p.config.camera);

    const json d = v.value("defaults", json::object());
    p.defaultJets       = d.value("jets",       false);
    p.defaultBLR        = d.value("blr",        false);
    p.defaultOrbBody    = d.value("orbBody",    false);
    p.defaultDoppler    = d.value("doppler",    false);
    p.defaultHostGalaxy = d.value("hostGalaxy", false);
    p.defaultLAB        = d.value("lab",        false);
    p.defaultCGM        = d.value("cgm",        false);

    p.isBinaryWithBarycenter = e.barycentric;
    p.companionMassSolar     = e.companionMassSolar;
    for (const auto& b : e.bodies)
        p.galaxyBodies.push_back({ bodyType(b.type), (float)b.semiMajorRs,
                                   (float)b.eccentricity, (float)b.inclination,
                                   b.label });
    return p;
}

// Entry point: JSON catalog if found, compiled-in profiles otherwise.
inline std::vector<BlackHoleProfile> loadProfiles(const ResourceManager& res) {
    std::string dir = res.findDir("catalog");
    if (!dir.empty()) {
        auto lr = catalog::loadDir(dir);
        for (const auto& err : lr.errors)
            std::cerr << "Catalog: skipped " << err << "\n";
        if (lr.ok()) {
            std::vector<BlackHoleProfile> out;
            out.reserve(lr.entries.size());
            for (const auto& e : lr.entries) out.push_back(toProfile(e));
            std::cerr << "Catalog: " << out.size() << " profiles from " << dir << "\n";
            return out;
        }
    }
    std::cerr << "Catalog: no JSON found (searched: " << res.paths()
              << "); using built-in profiles\n";
    auto arr = profiles::allProfiles();
    return { arr.begin(), arr.end() };
}

} // namespace bh3d_catalog
