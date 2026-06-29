#pragma once
// ============================================================
// orbital_body.hpp: Keplerian orbit with GR apsidal precession
// ============================================================

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <cmath>

#include "config.hpp"

class OrbitalBody {
public:
    explicit OrbitalBody(const cfg::OrbitalConfig& c = {})
        : cfg_(c)
    {
        // Compute derived constants once
        const float a = cfg_.semiMajor;
        const float e = cfg_.eccentricity;
        period_ = 2.0f * glm::pi<float>() * std::sqrt(a * a * a); // T = 2π√(a³/M), M≈1
        precessionPerOrbit_ = 3.0f * glm::pi<float>() / (a * (1.0f - e * e));
        precessionRate_ = precessionPerOrbit_ / period_;
        meanMotion_ = 2.0f * glm::pi<float>() / period_;
    }

    /*--------- Per-frame update ---------*/
    void update(float dt) {
        const float scaledDt = dt * timeScale_;

        // Advance mean anomaly
        meanAnomaly_ += meanMotion_ * scaledDt;
        meanAnomaly_ = std::fmod(meanAnomaly_, 2.0f * glm::pi<float>());

        // Accumulate GR precession
        precession_ += precessionRate_ * scaledDt;

        // Solve Kepler's equation: M = E − e sin(E) via Newton–Raphson
        // with convergence early-out.
        // Kepler himself solved this by hand in 1609. We have a computer and it's still annoying.
        const float e = cfg_.eccentricity;
        float E = meanAnomaly_; // initial guess
        for (int k = 0; k < cfg_.keplerIterMax; ++k) {
            float dE = (E - e * std::sin(E) - meanAnomaly_)
                     / (1.0f - e * std::cos(E));
            E -= dE;
            if (std::abs(dE) < cfg_.keplerTolerance) break; // converged
        }

        // True anomaly from eccentric anomaly
        float trueAnomaly = 2.0f * std::atan2(
            std::sqrt(1.0f + e) * std::sin(E * 0.5f),
            std::sqrt(1.0f - e) * std::cos(E * 0.5f));

        // Orbital radius from conic equation: r = a(1−e²) / (1 + e cos ν)
        // four hundred years of orbital mechanics distilled into one line
        float r = cfg_.semiMajor * (1.0f - e * e)
                / (1.0f + e * std::cos(trueAnomaly));

        // Position in orbital plane (with precession rotating the ellipse)
        float totalAngle = trueAnomaly + precession_;
        float xOrb = r * std::cos(totalAngle);
        float zOrb = r * std::sin(totalAngle);

        // Rotate orbital plane by inclination around X-axis
        float incl = cfg_.inclination;
        position_ = glm::vec3(xOrb,
                               zOrb * std::sin(incl),
                               zOrb * std::cos(incl));
    }

    /*--------- External perturbation (merger secondary) ---------*/
    // Accumulate a gravitational tug from an external mass (the in-spiralling
    // merger secondary). The displacement is layered on top of the Keplerian
    // position and persists, so the orbit visibly bends. A light velocity drag
    // keeps it bounded; the offset is clamped to maxOffset so a body that is
    // flung outward streaks to the edge of its region instead of vanishing.
    void applyExternalAccel(const glm::vec3& accel, float dt, float maxOffset) {
        perturbVel_ += accel * dt;
        perturbVel_ *= 0.992f;                 // light drag bounds runaway
        perturb_    += perturbVel_ * dt;
        float len = glm::length(perturb_);
        if (len > maxOffset && len > 0.0f) perturb_ *= (maxOffset / len);
    }
    void resetPerturbation() {
        perturb_    = glm::vec3(0.0f);
        perturbVel_ = glm::vec3(0.0f);
        disrupted_  = false;
    }
    void setDisrupted()    { disrupted_ = true; }
    bool disrupted() const { return disrupted_; }
    // Unperturbed Keplerian position (reference orbit, before external tugs).
    glm::vec3 keplerPosition() const { return position_; }

    /*--------- Phase control ---------*/
    // Set the initial orbital phase (mean anomaly in radians). Used to
    // stagger multiple bodies so they don't all start at periapsis on the
    // same frame, which is especially noticeable for high-eccentricity
    // orbits where periapsis lies very close to (or inside) the BH.
    void setInitialPhase(float meanAnomalyRad) {
        meanAnomaly_ = std::fmod(meanAnomalyRad, 2.0f * glm::pi<float>());
        update(0.0f); // refresh position_ so the first frame already shows the staggered phase
    }

    /*--------- Speed control ---------*/
    void setTimeScale(float s) { timeScale_ = s; }
    float timeScale() const     { return timeScale_; }
    // TODO: expose a smooth lerp between normal and fast timescale instead of snapping
    void toggleFast() {
        if (timeScale_ < cfg_.fastTimeScale)
            timeScale_ = cfg_.fastTimeScale;
        else
            timeScale_ = cfg_.timeScale;
    }
    bool isFast() const { return timeScale_ >= cfg_.fastTimeScale; }

    /*--------- Accessors ---------*/
    // World position including any accumulated external perturbation.
    glm::vec3 position()   const { return position_ + perturb_; }
    float     bodyRadius() const { return cfg_.bodyRadius; }
    glm::vec3 bodyColor()  const { return cfg_.bodyColor; }
    int       bodyType()   const { return cfg_.bodyType; }
    const cfg::OrbitalConfig& config() const { return cfg_; }

private:
    cfg::OrbitalConfig cfg_;

    float period_;
    float precessionPerOrbit_;
    float precessionRate_;
    float meanMotion_;

    float meanAnomaly_  = 0.0f;
    float precession_   = 0.0f;
    float timeScale_    = 0.5f; // default from config
    glm::vec3 position_ = glm::vec3(0.0f);

    // External-perturbation state (merger secondary tug + tidal disruption).
    glm::vec3 perturb_    = glm::vec3(0.0f);
    glm::vec3 perturbVel_ = glm::vec3(0.0f);
    bool      disrupted_  = false;
};
