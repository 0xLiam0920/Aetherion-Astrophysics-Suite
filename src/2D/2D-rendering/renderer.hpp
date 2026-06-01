/* ------------------------------------------------ *\
   Renderer.hpp, SFML-based 2D renderer for black hole visualization.
   All visual data is pre-computed by the visualization layer; this class
   just draws shapes and sprites based on that data.

   Part of Aetherion Suite: https://github.com/0xLiam0920/Aetherion-Astrophysics-Suite
\* ------------------------------------------------ */
#pragma once
#include <SFML/Graphics.hpp>
#include "../2D-core/body_visual.hpp"
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <cstdint>
#include <algorithm>

// Renderer, keeps it dumb.  NEVER computes physics.
// All visual data is pre-computed by the visualization layer.
class Renderer {
    sf::RenderWindow& window_;
    sf::CircleShape   horizonShape_;
    sf::CircleShape   photonRing_;
    sf::Font          font_;
    sf::Text          infoText_;
    sf::Texture       diskTexture_;
    sf::Sprite        diskSprite_;
    bool              diskTextureLoaded_;
    bool              lightMode_ = false;

public:
    // Toggle light/dark HUD palette. Affects clear colour, starfield, HUD
    // text, and the data/controls panels. Does not retint physics visuals
    // (horizon, jets, accretion disk) which stay scientifically meaningful.
    void setLightMode(bool on) {
        lightMode_ = on;
        infoText_.setFillColor(on ? sf::Color(20, 22, 36)
                                  : sf::Color::White);
    }
    bool lightMode() const { return lightMode_; }

    explicit Renderer(sf::RenderWindow& win)
        : window_(win)
        , horizonShape_(10.0f)
        , photonRing_(15.0f)
        , infoText_(font_, "", 14)
        , diskSprite_(diskTexture_)
        , diskTextureLoaded_(false)
    {
        // Pre-check existence before calling openFromFile so SFML never prints
        // a "Failed to load font" error for paths that simply don't exist on
        // the current platform (e.g. Linux paths on macOS and vice versa).
        // Flatpak runtime: /usr/share/fonts/dejavu/  (no truetype/ subdir)
        // Host system:     /usr/share/fonts/truetype/dejavu/
        // Host via flatpak bridge: /run/host/fonts/truetype/dejavu/
        auto tryFont = [&](const std::string& path) -> bool {
            std::ifstream probe(path);
            return probe.good() && font_.openFromFile(path);
        }; /* ----- Try these common system font paths, in order ----- */
        bool fontLoaded =
            tryFont("/usr/share/fonts/dejavu/DejaVuSans.ttf")                       ||
            tryFont("/usr/share/fonts/liberation-fonts/LiberationSans-Regular.ttf") ||
            tryFont("/run/host/fonts/truetype/dejavu/DejaVuSans.ttf")               ||
            tryFont("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")              ||
            tryFont("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf") ||
            tryFont("/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf")               ||
            tryFont("/System/Library/Fonts/Helvetica.ttc")                          ||
            tryFont("C:/Windows/Fonts/segoeui.ttf")                                 ||
            tryFont("C:/Windows/Fonts/arial.ttf")                                   ||
            tryFont("C:/Windows/Fonts/tahoma.ttf");
        if (!fontLoaded)
            std::cerr << "Warning: could not load system font; text may not display.\n";
        infoText_ = sf::Text(font_, "", 14);
        infoText_.setFillColor(sf::Color::White);
        infoText_.setPosition({8.0f, 8.0f});

        // Build a list of candidate paths for the disk texture,
        // covering: dev builds (relative), flatpak (/app/share), XDG installs.
        // Pre-check existence with ifstream so SFML doesn't log spurious errors.
        auto tryDiskTexture = [&](const std::string& p) {
            if (!diskTextureLoaded_) {
                std::ifstream probe(p);
                if (probe.good())
                    diskTextureLoaded_ = diskTexture_.loadFromFile(p);
            }
        };
        tryDiskTexture("disk_texture.png");
        tryDiskTexture("src/disk_texture.png");
        tryDiskTexture("../src/disk_texture.png");
        tryDiskTexture("/app/share/aetherionsuite/disk_texture.png");
        tryDiskTexture("/usr/local/share/aetherionsuite/disk_texture.png");
        tryDiskTexture("/usr/share/aetherionsuite/disk_texture.png");
        // Also walk XDG_DATA_DIRS (set by Flatpak and system installs)
        if (!diskTextureLoaded_) {
            const char* xdg = std::getenv("XDG_DATA_DIRS");
            if (xdg && xdg[0]) {
                std::string dirs(xdg);
                std::string::size_type start = 0, end;
                while (!diskTextureLoaded_ && (end = dirs.find(':', start)) != std::string::npos) {
                    tryDiskTexture(dirs.substr(start, end - start) + "/aetherionsuite/disk_texture.png");
                    start = end + 1;
                }
                if (!diskTextureLoaded_)
                    tryDiskTexture(dirs.substr(start) + "/aetherionsuite/disk_texture.png");
            }
        }
        if (diskTextureLoaded_) {
            diskSprite_ = sf::Sprite(diskTexture_);
            diskSprite_.setOrigin(sf::Vector2f(
                diskTexture_.getSize().x / 2.f,
                diskTexture_.getSize().y / 2.f));
        }

        horizonShape_.setFillColor(sf::Color::Black);
        horizonShape_.setOutlineThickness(0.0f);

        photonRing_.setFillColor(sf::Color::Transparent);
        photonRing_.setOutlineThickness(2.0f);
        photonRing_.setOutlineColor(sf::Color(180, 180, 255, 200));
    }

    /*--------- Frame lifecycle ---------*/
    void beginFrame() {
        window_.clear(lightMode_ ? sf::Color(232, 234, 244)
                                 : sf::Color(10, 10, 30));
    }
    void endFrame()   { window_.display(); }

    /*--------- Background ---------*/
    // Three parallax star layers with per-star twinkle and slow horizontal
    // drift. `tSec` is total elapsed wall-clock time in seconds; pass 0 for a
    // static field. Drift wraps modulo the viewport so there's no popping.
    void drawStarfield(float tSec = 0.0f) {
        auto winSz = window_.getSize();
        const float W = (float)winSz.x;
        const float H = (float)winSz.y;

        // Three layers: (count, base alpha, drift px/sec, radius, twinkle hz).
        struct Layer { int n; float a; float vx; float r; float tw; };
        const Layer layers[3] = {
            { 120, 1.00f, 1.5f, 1.0f, 0.6f },  // foreground sparkle
            {  90, 0.70f, 0.6f, 0.8f, 0.4f },  // mid
            {  60, 0.45f, 0.25f, 0.6f, 0.25f } // distant
        };

        for (int L = 0; L < 3; ++L) {
            const Layer& lay = layers[L];
            for (int i = 0; i < lay.n; ++i) {
                const float seed = (float)(i + L * 911);
                float bx = std::fmod(seed * 37.1234f + 300.5f, W);
                float by = std::fmod(seed * 71.4321f + 100.5f, H);
                float sx = std::fmod(bx + tSec * lay.vx + W, W);
                float sy = by;

                // Twinkle: per-star phase, smooth sine in [0.55, 1.0].
                float phase = seed * 1.9173f;
                float tw = 0.775f + 0.225f * std::sin(tSec * lay.tw * 6.2831853f + phase);

                sf::Color base = lightMode_ ? sf::Color(120, 130, 170)
                                            : sf::Color(220, 220, 255);
                base.a = (std::uint8_t)std::clamp(int(lay.a * tw * (lightMode_ ? 110.0f : 80.0f)),
                                                  0, 255);

                sf::CircleShape s(lay.r);
                s.setOrigin({lay.r, lay.r});
                s.setPosition({sx, sy});
                s.setFillColor(base);
                window_.draw(s);
            }
        }
    }

    /*--------- Black hole visuals (pre-computed pixel values) ---------*/
    void drawHorizon(sf::Vector2f center, float radiusPx) {
        horizonShape_.setRadius(radiusPx);
        horizonShape_.setOrigin({radiusPx, radiusPx});
        horizonShape_.setPosition(center);
        window_.draw(horizonShape_);
    }

    void drawPhotonSphere(sf::Vector2f center, float radiusPx) {
        photonRing_.setRadius(radiusPx);
        photonRing_.setOrigin({radiusPx, radiusPx});
        photonRing_.setPosition(center);
        window_.draw(photonRing_);
    }

    void drawAccretionDisk(sf::Vector2f center, float radiusPx) {
        drawAccretionDisk(center, radiusPx, sf::Color::White);
    }

    // Tinted variant: applies an SFML sprite colour multiplier so the disk can
    // be re-coloured (e.g. purple for the inspiralling secondary BH so equal-
    // mass mergers stay visually distinguishable).
    void drawAccretionDisk(sf::Vector2f center, float radiusPx, sf::Color tint) {
        if (!diskTextureLoaded_) return;
        float texRadius = diskTexture_.getSize().x / 2.0f;
        float scale = radiusPx / texRadius;
        diskSprite_.setScale({scale, scale});
        diskSprite_.setPosition(center);
        diskSprite_.setColor(tint);
        window_.draw(diskSprite_);
        diskSprite_.setColor(sf::Color::White);  // restore default for next caller
    }

    /*--------- Rays (vertices pre-computed by RayVisualizer) ---------*/
    void drawRayPath(const std::vector<sf::Vertex>& verts) {
        if (verts.size() >= 2)
            window_.draw(verts.data(), verts.size(), sf::PrimitiveType::LineStrip);
    }

    void drawCapturedRay(const sf::Vertex& v0, const sf::Vertex& v1) {
        sf::Vertex line[] = { v0, v1 };
        window_.draw(line, 2, sf::PrimitiveType::Lines);
    }

    /*--------- Bodies (BodyVisual pre-computed by OrbitVisualizer) ---------*/
    void drawBody(const BodyVisual& vis) {
        sf::RectangleShape rect(sf::Vector2f(vis.width, vis.height));
        rect.setOrigin(rect.getSize() * 0.5f);
        rect.setPosition(vis.screenPos);
        rect.setRotation(sf::degrees(vis.angleDeg));
        rect.setFillColor(vis.color);
        window_.draw(rect);
    }

    // ── Orbit path (VertexArray pre-computed by OrbitVisualizer)
    void drawOrbitPath(const sf::VertexArray& path) {
        if (path.getVertexCount() > 1) window_.draw(path);
    }

    /*--------- Influence zones (labeled radial rings) ---------*/
    void drawInfluenceZone(sf::Vector2f center, float radiusPx,
                           sf::Color color, const std::string& label) {
        // Dashed-style ring
        sf::CircleShape ring(radiusPx);
        ring.setOrigin({radiusPx, radiusPx});
        ring.setPosition(center);
        ring.setFillColor(sf::Color::Transparent);
        ring.setOutlineThickness(1.5f);
        ring.setOutlineColor(color);
        window_.draw(ring);

        // Label at top of ring
        sf::Text lbl(font_, label, 11);
        lbl.setFillColor(color);
        lbl.setPosition({center.x - (float)label.size() * 3.0f,
                         center.y - radiusPx - 14.0f});
        window_.draw(lbl);
    }

    /*--------- Relativistic jet cones ---------*/
    void drawJetCones(sf::Vector2f center, float lengthPx) {
        drawJetCones(center, lengthPx,
                     sf::Color(120, 180, 255, 50),
                     sf::Color(150, 200, 255, 100));
    }

    // Colour-customised jets, used by the merger renderer so the secondary BH
    // gets purple jets matching its identity colour.
    void drawJetCones(sf::Vector2f center, float lengthPx,
                      sf::Color fillColor, sf::Color outlineColor) {
        // Upward jet
        sf::ConvexShape jetUp;
        jetUp.setPointCount(3);
        jetUp.setPoint(0, {center.x - 4.0f, center.y});
        jetUp.setPoint(1, {center.x + 4.0f, center.y});
        jetUp.setPoint(2, {center.x, center.y - lengthPx});
        jetUp.setFillColor(fillColor);
        jetUp.setOutlineThickness(1.0f);
        jetUp.setOutlineColor(outlineColor);
        window_.draw(jetUp);

        // Downward jet
        sf::ConvexShape jetDown;
        jetDown.setPointCount(3);
        jetDown.setPoint(0, {center.x - 4.0f, center.y});
        jetDown.setPoint(1, {center.x + 4.0f, center.y});
        jetDown.setPoint(2, {center.x, center.y + lengthPx});
        jetDown.setFillColor(fillColor);
        jetDown.setOutlineThickness(1.0f);
        jetDown.setOutlineColor(outlineColor);
        window_.draw(jetDown);
    }

    /*--------- Body label text ---------*/
    void drawBodyLabel(sf::Vector2f screenPos, const std::string& text, sf::Color color) {
        sf::Text lbl(font_, text, 10);
        lbl.setFillColor(color);
        lbl.setPosition({screenPos.x + 8.0f, screenPos.y - 6.0f});
        window_.draw(lbl);
    }

    /*--------- Barycenter / center-of-mass marker ---------*/
    // Draws a small cross-hair and "Center of Mass" label at the given screen position.
    // Used for binary presets (Gaia BH1/2/3) where both objects orbit the barycenter.
    void drawBarycenterMarker(sf::Vector2f pos) {
        constexpr float ARM = 7.0f;
        sf::Color col(200, 255, 200, 200);

        sf::RectangleShape h(sf::Vector2f(ARM * 2.0f, 1.5f));
        h.setOrigin({ARM, 0.75f});
        h.setPosition(pos);
        h.setFillColor(col);
        window_.draw(h);

        sf::RectangleShape v(sf::Vector2f(1.5f, ARM * 2.0f));
        v.setOrigin({0.75f, ARM});
        v.setPosition(pos);
        v.setFillColor(col);
        window_.draw(v);

        sf::Text lbl(font_, "Center of Mass", 10);
        lbl.setFillColor(col);
        lbl.setPosition({pos.x + ARM + 4.0f, pos.y - 7.0f});
        window_.draw(lbl);
    }

    /*--------- HUD ---------*/
    void drawHUD(const std::string& s) {
        infoText_.setString(s);
        window_.draw(infoText_);
    }

    /*--------- Data panel (right side, uses view coords) ---------*/
    void drawDataPanel(const std::string& text, float viewW, float viewH) {
        float panelW = 260.0f;
        float panelX = viewW - panelW - 4.0f;

        const sf::Color bgCol      = lightMode_ ? sf::Color(248, 249, 252, 235)
                                                : sf::Color(0, 0, 0, 180);
        const sf::Color outlineCol = lightMode_ ? sf::Color(160, 170, 200, 200)
                                                : sf::Color(80, 120, 200, 150);
        const sf::Color textCol    = lightMode_ ? sf::Color(20, 22, 36)
                                                : sf::Color(200, 220, 255);

        sf::RectangleShape bg(sf::Vector2f(panelW, viewH - 8.0f));
        bg.setPosition({panelX, 4.0f});
        bg.setFillColor(bgCol);
        bg.setOutlineThickness(1.0f);
        bg.setOutlineColor(outlineCol);
        window_.draw(bg);

        sf::Text panelText(font_, text, 11);
        panelText.setFillColor(textCol);
        panelText.setPosition({panelX + 6.0f, 8.0f});
        window_.draw(panelText);
    }

    /*--------- Controls panel (left side, toggleable) ---------*/
    void drawControlsPanel(float viewH) {
        float panelW = 280.0f;
        float panelX = 4.0f;
        float panelY = 4.0f;

        const sf::Color bgCol      = lightMode_ ? sf::Color(248, 249, 252, 235)
                                                : sf::Color(0, 0, 0, 190);
        const sf::Color outlineCol = lightMode_ ? sf::Color(150, 160, 190, 200)
                                                : sf::Color(120, 140, 180, 150);
        const sf::Color hdrCol     = lightMode_ ? sf::Color(180, 80, 30)
                                                : sf::Color(255, 220, 120);
        const sf::Color sectCol    = lightMode_ ? sf::Color(40, 80, 160, 230)
                                                : sf::Color(160, 200, 255, 220);
        const sf::Color lineCol    = lightMode_ ? sf::Color(30, 32, 52, 230)
                                                : sf::Color(200, 200, 210, 200);
        const sf::Color dividerCol = lightMode_ ? sf::Color(180, 120, 60, 200)
                                                : sf::Color(180, 120, 60, 140);
        const sf::Color advHdrCol  = lightMode_ ? sf::Color(180, 90, 20)
                                                : sf::Color(220, 160, 80, 220);
        const sf::Color advLineCol = lightMode_ ? sf::Color(60, 62, 80, 220)
                                                : sf::Color(170, 170, 185, 180);

        sf::RectangleShape bg(sf::Vector2f(panelW, viewH - 8.0f));
        bg.setPosition({panelX, panelY});
        bg.setFillColor(bgCol);
        bg.setOutlineThickness(1.0f);
        bg.setOutlineColor(outlineCol);
        window_.draw(bg);

        // Header
        sf::Text header(font_, "CONTROLS  (? to hide)", 13);
        header.setFillColor(hdrCol);
        header.setPosition({panelX + 8.0f, panelY + 6.0f});
        window_.draw(header);

        // Primary controls
        const char* primaryLines[] = {
            "",
            "--- General ---",
            "Space     Pause / resume",
            "R         Reset simulation",
            "V         Reset view (camera)",
            "B         Toggle light / dark",
            "L         Open preset reference (web)",
            "- / =     Slow down / speed up",
            "?         Toggle this panel",
            "",
            "--- View ---",
            "Up/Down   Zoom (preset) / Mass",
            "S         Toggle ray display",
            "P         Toggle photon sphere",
            "I         Toggle influence zones",
            "",
            "--- Orbit ---",
            "Left/Right  Semi-major axis",
            "Q / E       Eccentricity +/-",
            "",
            "--- Presets ---",
            "T         Toggle preset mode",
            "[ / ]     Cycle presets",
            "G         Toggle galaxy system",
        };
        constexpr int nPrimary = sizeof(primaryLines) / sizeof(primaryLines[0]);

        float y = panelY + 26.0f;
        for (int i = 0; i < nPrimary; ++i) {
            const char* line = primaryLines[i];
            bool isSectionHeader = (line[0] == '-');
            sf::Text lt(font_, line, 11);
            lt.setFillColor(isSectionHeader ? sectCol : lineCol);
            lt.setPosition({panelX + 10.0f, y});
            window_.draw(lt);
            y += isSectionHeader ? 16.0f : 14.0f;
        }

        // Divider before advanced section
        y += 4.0f;
        sf::RectangleShape divider(sf::Vector2f(panelW - 20.0f, 1.0f));
        divider.setPosition({panelX + 10.0f, y});
        divider.setFillColor(dividerCol);
        window_.draw(divider);
        y += 5.0f;

        sf::Text advHdr(font_, "--- Research / Debug ---", 11);
        advHdr.setFillColor(advHdrCol);
        advHdr.setPosition({panelX + 10.0f, y});
        window_.draw(advHdr);
        y += 16.0f;

        const char* advLines[] = {
            "F         Time-dilation map",
            "C         Caustic highlights",
            "D         Data panel",
            "H         Cycle selected body",
            "X         Export CSV data",
            ", / .     Gas temp down/up",
            "",
            "--- Scenarios ---",
            "1         ISCO validation test",
            "2         Photon sphere test",
            "3         Radial infall test",
            "4         Tidal disruption demo",
            "5         Pulsar orbital sim",
            "6         Merger event menu",
            "7         Custom-body creator",
            "M         Merger pacing",
        };
        constexpr int nAdv = sizeof(advLines) / sizeof(advLines[0]);
        for (int i = 0; i < nAdv; ++i) {
            const char* line = advLines[i];
            bool isSectionHeader = (line[0] == '-');
            sf::Text lt(font_, line, 11);
            lt.setFillColor(isSectionHeader ? sectCol : advLineCol);
            lt.setPosition({panelX + 10.0f, y});
            window_.draw(lt);
            y += isSectionHeader ? 16.0f : 14.0f;
        }
    }

    /*--------- Pulsar body (clean compact-secondary look) ---------*/
    // Mirrors drawMergerCompactSecondary's pulsar branch — kind-coloured glow +
    // bright core + two opposed sweeping beams. Signature preserved so callers
    // don't break; the magnetosphere / flux-tube parameters are no longer used
    // visually (their physics is still surfaced in the HUD / data panel).
    void drawPulsarBody(sf::Vector2f screenPos,
                        float spinPhase_rad,
                        float precPhase_rad,
                        [[maybe_unused]] float lcRadius_px,
                        [[maybe_unused]] bool  inLC,
                        [[maybe_unused]] sf::Vector2f bhPos,
                        [[maybe_unused]] double magPower_ergs)
    {
        constexpr float CORE_R    = 5.0f;
        constexpr float GLOW_R    = 14.0f;
        constexpr float BEAM_LEN  = 80.0f;
        constexpr float PI        = 3.14159265f;

        // Outer soft glow
        sf::CircleShape g(GLOW_R);
        g.setOrigin({GLOW_R, GLOW_R});
        g.setPosition(screenPos);
        g.setFillColor(sf::Color(150, 200, 255, 110));
        window_.draw(g);

        // Bright core
        sf::CircleShape c(CORE_R);
        c.setOrigin({CORE_R, CORE_R});
        c.setPosition(screenPos);
        c.setFillColor(sf::Color(235, 245, 255));
        window_.draw(c);

        // Two opposed sweeping beams. Spin phase drives the sweep; precession
        // adds a slow drift so the beam orientation isn't dead-aligned forever.
        const float beamAng = spinPhase_rad + 0.25f * precPhase_rad;
        for (int s = 0; s < 2; ++s) {
            const float a = beamAng + (s ? PI : 0.0f);
            sf::Vertex line[2];
            line[0].position = screenPos;
            line[0].color    = sf::Color(220, 235, 255, 230);
            line[1].position = {screenPos.x + std::cos(a) * BEAM_LEN,
                                screenPos.y + std::sin(a) * BEAM_LEN};
            line[1].color    = sf::Color(120, 180, 255, 0);
            window_.draw(line, 2, sf::PrimitiveType::Lines);
        }

        // Small label so the pulsar is unmistakable.
        sf::Text lbl(font_, "Pulsar", 10);
        lbl.setFillColor(sf::Color(220, 230, 255));
        lbl.setPosition({screenPos.x + GLOW_R + 3.0f, screenPos.y - 7.0f});
        window_.draw(lbl);
    }


    /*--------- Time dilation heatmap (concentric rings) ---------*/
    void drawTimeDilationMap(sf::Vector2f center, double M, double pixelsPerM) {
        const int nRings = 30;
        double rMax = 20.0 * M;
        double rMin = 2.1 * M;

        for (int i = 0; i < nRings; ++i) {
            double t = (double)i / (double)nRings;
            double r = rMax * (1.0 - t) + rMin * t;

            double fr = 1.0 - 2.0 * M / r;
            if (fr <= 0.0) continue;
            double gamma = 1.0 / std::sqrt(fr);

            // Map gamma [1, ~4] to color: blue → yellow → red
            double norm = std::min((gamma - 1.0) / 3.0, 1.0);
            uint8_t rc = (uint8_t)(norm * 255);
            uint8_t gc = (uint8_t)(std::max(0.0, 1.0 - std::abs(norm - 0.4) * 3.0) * 200);
            uint8_t bc = (uint8_t)((1.0 - norm) * 200);

            float radiusPx = (float)(r * pixelsPerM);
            sf::CircleShape ring(radiusPx);
            ring.setOrigin({radiusPx, radiusPx});
            ring.setPosition(center);
            ring.setFillColor(sf::Color(rc, gc, bc, 25));
            ring.setOutlineThickness(0.0f);
            window_.draw(ring);
        }

        // Labels at key radii
        auto drawLabel = [&](double r, const std::string& lbl) {
            float px = (float)(r * pixelsPerM);
            sf::Text t(font_, lbl, 9);
            t.setFillColor(sf::Color(200, 200, 255, 160));
            t.setPosition({center.x + px + 3.0f, center.y - 8.0f});
            window_.draw(t);
        };
        drawLabel(6.0 * M, "ISCO 6M");
        drawLabel(3.0 * M, "Photon 3M");
        drawLabel(2.0 * M, "Horizon 2M");
    }

    /*--------- Caustic point highlighting ---------*/
    void drawCausticMarker(sf::Vector2f pos, float intensity) {
        float radius = 3.0f + intensity * 2.0f;
        sf::CircleShape marker(radius);
        marker.setOrigin({radius, radius});
        marker.setPosition(pos);
        uint8_t alpha = (uint8_t)std::min(255.0f, 80.0f + intensity * 40.0f);
        marker.setFillColor(sf::Color(255, 255, 100, alpha));
        window_.draw(marker);
    }

    /*--------- Export status notification ---------*/
    void drawNotification(const std::string& text, float viewW, float viewH) {
        sf::Text notif(font_, text, 14);
        notif.setFillColor(sf::Color(100, 255, 100));
        // Anchor at the left margin so long export-path strings aren't clipped
        // by the right edge or the HUD panel. If the text is still wider than
        // the view (very long path + many filenames), shift left so the tail
        // remains visible.
        const float margin = 10.0f;
        const float textW = notif.getLocalBounds().size.x;
        float x = margin;
        if (textW + margin * 2.0f > viewW) {
            x = viewW - textW - margin;
            if (x < margin) x = margin;
        }
        notif.setPosition({x, viewH - 30.0f});
        window_.draw(notif);
    }

    /*--------- Merger: incoming BH ---------*/
    // Draws the secondary BH at `pos`. Renders (in order) accretion disk (if
    // diskPx > 0), purple jets (if jetLenPx > 0), purple outer glow, the dark
    // horizon core with a purple rim, and a label. The disk and jets are
    // tinted purple so the secondary stays visually distinct from the primary
    // even in an equal-mass merger.
    void drawMergerBH(sf::Vector2f pos, float radiusPx,
                      float diskPx = 0.0f, float jetLenPx = 0.0f,
                      const char* labelText = "Incoming BH") {
        // Accretion disk first (drawn beneath the horizon)
        if (diskPx > radiusPx + 0.5f) {
            drawAccretionDisk(pos, diskPx, sf::Color(210, 140, 255, 220));
        }

        // Relativistic jets for AGN-class secondaries
        if (jetLenPx > 0.0f) {
            drawJetCones(pos, jetLenPx,
                         sf::Color(200, 120, 255, 55),
                         sf::Color(220, 160, 255, 130));
        }

        // Outer glow (always)
        float glowR = radiusPx + 6.0f;
        sf::CircleShape glow(glowR);
        glow.setOrigin({glowR, glowR});
        glow.setPosition(pos);
        glow.setFillColor(sf::Color(180, 120, 255, 50));
        window_.draw(glow);

        // Core
        sf::CircleShape core(radiusPx);
        core.setOrigin({radiusPx, radiusPx});
        core.setPosition(pos);
        core.setFillColor(sf::Color(10, 10, 10));
        core.setOutlineThickness(2.0f);
        core.setOutlineColor(sf::Color(180, 100, 255, 200));
        window_.draw(core);

        // Label
        sf::Text lbl(font_, labelText, 10);
        lbl.setFillColor(sf::Color(200, 160, 255));
        lbl.setPosition({pos.x + radiusPx + 4.0f, pos.y - 7.0f});
        window_.draw(lbl);
    }

    // Compact / stellar secondary during inspiral. Not a black hole — just a
    // bright point with a kind-appropriate colour and label.
    // kindCode matches Simulation::MergerSecondaryKind: 1=NS, 2=Pulsar,
    // 3=Star, 4=WhiteDwarf (caller passes int to avoid header coupling).
    void drawMergerCompactSecondary(sf::Vector2f pos, int kindCode, double phase) {
        sf::Color core, glow;
        const char* label = "Secondary";
        float coreR = 3.5f;
        float glowR = 9.0f;
        switch (kindCode) {
            case 1: core = sf::Color(220, 235, 255); glow = sf::Color(120, 180, 255, 90);
                    label = "Neutron star"; break;
            case 2: core = sf::Color(235, 245, 255); glow = sf::Color(150, 200, 255, 110);
                    label = "Pulsar";        coreR = 4.0f; break;
            case 3: core = sf::Color(255, 230, 160); glow = sf::Color(255, 200, 100, 110);
                    label = "Star";          coreR = 4.5f; glowR = 12.0f; break;
            case 4: core = sf::Color(245, 245, 255); glow = sf::Color(200, 220, 255, 100);
                    label = "White dwarf";   coreR = 3.0f; break;
            default: core = sf::Color::White; glow = sf::Color(200, 200, 255, 80); break;
        }

        sf::CircleShape g(glowR);
        g.setOrigin({glowR, glowR});
        g.setPosition(pos);
        g.setFillColor(glow);
        window_.draw(g);

        sf::CircleShape c(coreR);
        c.setOrigin({coreR, coreR});
        c.setPosition(pos);
        c.setFillColor(core);
        window_.draw(c);

        // Pulsar beam — two thin opposed rays that sweep with the orbital phase.
        if (kindCode == 2) {
            const float beamLen = 22.0f;
            const float beamAng = static_cast<float>(phase * 4.0); // sweeps ~4x orbit rate
            for (int s = 0; s < 2; ++s) {
                const float a = beamAng + (s ? 3.14159265f : 0.0f);
                sf::Vertex line[2];
                line[0].position = pos;
                line[0].color    = sf::Color(220, 235, 255, 220);
                line[1].position = {pos.x + std::cos(a) * beamLen,
                                    pos.y + std::sin(a) * beamLen};
                line[1].color    = sf::Color(120, 180, 255, 0);
                window_.draw(line, 2, sf::PrimitiveType::Lines);
            }
        }

        sf::Text lbl(font_, label, 10);
        lbl.setFillColor(sf::Color(220, 230, 255));
        lbl.setPosition({pos.x + glowR + 3.0f, pos.y - 7.0f});
        window_.draw(lbl);
    }

    // Draw the fading death-spiral trail for one of the two BHs.
    // trailColor: base RGB of this BH's streak (purple for secondary, cyan for primary).
    void drawMergerTrail(
        const std::deque<std::pair<double,double>>& trail,
        sf::Vector2f bhCenter,
        float pixelsPerM,
        sf::Color baseColor)
    {
        if (trail.size() < 2) return;
        sf::VertexArray line(sf::PrimitiveType::LineStrip, trail.size());
        for (size_t i = 0; i < trail.size(); ++i) {
            float wx = (float)trail[i].first;
            float wy = (float)trail[i].second;
            line[i].position = {
                bhCenter.x + wx * pixelsPerM,
                bhCenter.y - wy * pixelsPerM   // y-flip: world +y = screen up
            };
            // Fade from transparent (old) to fully opaque (recent)
            float t = (float)i / (float)(trail.size() - 1);
            line[i].color = sf::Color(
                baseColor.r,
                baseColor.g,
                baseColor.b,
                (uint8_t)(t * 200.0f)
            );
        }
        window_.draw(line);
    }

    /*--------- Merger: white flash at the moment of coalescence ---------*/
    // alpha: 0.0 = transparent, 1.0 = full white
    void drawMergeFlash(float alpha, float viewW, float viewH) {
        sf::RectangleShape flash(sf::Vector2f(viewW, viewH));
        flash.setPosition({0.f, 0.f});
        uint8_t a = (uint8_t)(std::min(1.0f, alpha) * 220.0f);
        flash.setFillColor(sf::Color(255, 255, 255, a));
        window_.draw(flash);
    }

    /*--------- Merger: selection menu overlay ---------*/
    // presets / numPresets: the BH preset table.
    // secondaries / numSecondaries: compact / stellar secondaries (NS / Pulsar /
    // Star / WD), rendered between the BH list and the trailing “Custom” row.
    // The compact list is optional — pass nullptr/0 to render only BH presets +
    // Custom (legacy behaviour, used by 3D / other callers).
    template<typename MergerMenuState, typename BHPreset, typename SecondaryPreset = int>
    void drawMergerMenu(const MergerMenuState& state,
                        const BHPreset* presets, int numPresets,
                        float viewW, float viewH,
                        const SecondaryPreset* secondaries = nullptr,
                        int numSecondaries = 0)
    {
        // Theme palette
        const bool L = lightMode_;
        const sf::Color BG_FILL    = L ? sf::Color(248, 249, 252, 235) : sf::Color(8, 8, 24, 230);
        const sf::Color BG_OUTLINE = L ? sf::Color(150, 90, 200, 220) : sf::Color(160, 80, 255, 200);
        const sf::Color TITLE_COL  = L ? sf::Color(120, 40, 160)      : sf::Color(220, 160, 255);
        const sf::Color SUB_COL    = L ? sf::Color(80, 80, 120, 220)  : sf::Color(160, 160, 200, 200);
        const sf::Color DIV_COL    = L ? sf::Color(150, 90, 200, 110) : sf::Color(160, 80, 255, 120);
        const sf::Color SEL_FILL   = L ? sf::Color(150, 90, 200, 60)  : sf::Color(160, 60, 255, 80);
        const sf::Color ARROW_COL  = L ? sf::Color(120, 40, 160)      : sf::Color(220, 160, 255);
        const sf::Color BH_NAME_SEL = L ? sf::Color(60, 20, 90)       : sf::Color(255, 220, 255);
        const sf::Color BH_NAME     = L ? sf::Color(30, 32, 48)       : sf::Color(210, 210, 230);
        const sf::Color BH_MASS     = L ? sf::Color(40, 80, 140, 200) : sf::Color(160, 200, 255, 180);
        const sf::Color SEC_NAME_SEL = L ? sf::Color(20, 80, 50)      : sf::Color(220, 255, 230);
        const sf::Color SEC_NAME     = L ? sf::Color(40, 100, 70)     : sf::Color(190, 230, 210);
        const sf::Color SEC_KIND     = L ? sf::Color(40, 120, 100, 200) : sf::Color(160, 220, 200, 180);
        const sf::Color CUST_SEL    = L ? sf::Color(140, 100, 20)     : sf::Color(255, 255, 180);
        const sf::Color CUST_NORM   = L ? sf::Color(120, 90, 40)      : sf::Color(210, 210, 160);
        const sf::Color DIV2_COL    = L ? sf::Color(150, 90, 200, 80) : sf::Color(160, 80, 255, 80);
        const sf::Color HINT_COL    = L ? sf::Color(90, 90, 120, 200) : sf::Color(140, 140, 170, 180);

        constexpr float PANEL_W  = 440.0f;
        constexpr float ROW_H    = 18.0f;
        constexpr float HEADER_H = 48.0f;
        constexpr float FOOTER_H = 36.0f;
        // Row layout: BH presets, compact secondaries, then a single Custom row.
        int totalRows = numPresets + numSecondaries + 1;
        float panelH  = HEADER_H + totalRows * ROW_H + FOOTER_H + 8.0f;
        float panelX  = (viewW - PANEL_W) * 0.5f;
        float panelY  = (viewH - panelH)  * 0.5f;

        // Background
        sf::RectangleShape bg(sf::Vector2f(PANEL_W, panelH));
        bg.setPosition({panelX, panelY});
        bg.setFillColor(BG_FILL);
        bg.setOutlineThickness(2.0f);
        bg.setOutlineColor(BG_OUTLINE);
        window_.draw(bg);

        // Title
        sf::Text title(font_, "MERGER EVENT  (Enter: select  Esc: cancel)", 13);
        title.setFillColor(TITLE_COL);
        title.setPosition({panelX + 10.0f, panelY + 8.0f});
        window_.draw(title);

        sf::Text subtitle(font_, "Choose the merging black hole:", 11);
        subtitle.setFillColor(SUB_COL);
        subtitle.setPosition({panelX + 10.0f, panelY + 26.0f});
        window_.draw(subtitle);

        // Divider
        sf::RectangleShape div(sf::Vector2f(PANEL_W - 20.0f, 1.0f));
        div.setPosition({panelX + 10.0f, panelY + HEADER_H - 4.0f});
        div.setFillColor(DIV_COL);
        window_.draw(div);

        // Preset rows
        float rowY = panelY + HEADER_H;
        auto formatMass = [](double msun) -> std::string {
            std::ostringstream os;
            if (msun >= 1e9)       os << std::scientific << std::setprecision(2) << msun << " Msun";
            else if (msun >= 1e6)  os << std::fixed << std::setprecision(2) << (msun / 1e6) << "M Msun";
            else if (msun >= 1e3)  os << std::fixed << std::setprecision(0) << msun << " Msun";
            else                   os << std::fixed << std::setprecision(2) << msun << " Msun";
            return os.str();
        };

        for (int i = 0; i < totalRows; ++i) {
            bool isSelected = (state.selectedIdx == i);
            bool isBH       = (i < numPresets);
            bool isSecondary = (!isBH && i < numPresets + numSecondaries);
            bool isCustom   = (!isBH && !isSecondary);

            // Highlight bar
            if (isSelected) {
                sf::RectangleShape sel(sf::Vector2f(PANEL_W - 12.0f, ROW_H - 1.0f));
                sel.setPosition({panelX + 6.0f, rowY + 1.0f});
                sel.setFillColor(SEL_FILL);
                window_.draw(sel);
            }

            // Cursor arrow
            sf::Text arrow(font_, isSelected ? ">" : " ", 11);
            arrow.setFillColor(ARROW_COL);
            arrow.setPosition({panelX + 8.0f, rowY + 2.0f});
            window_.draw(arrow);

            if (isBH) {
                // Name
                sf::Text name(font_, presets[i].name, 11);
                name.setFillColor(isSelected ? BH_NAME_SEL : BH_NAME);
                name.setPosition({panelX + 22.0f, rowY + 2.0f});
                window_.draw(name);

                // Mass (right-aligned)
                std::string massStr = formatMass(presets[i].massSolar);
                sf::Text massText(font_, massStr, 10);
                massText.setFillColor(BH_MASS);
                float massX = panelX + PANEL_W - massText.getLocalBounds().size.x - 14.0f;
                massText.setPosition({massX, rowY + 3.0f});
                window_.draw(massText);
            } else if (isSecondary) {
                const auto& sec = secondaries[i - numPresets];
                sf::Text name(font_, sec.name, 11);
                name.setFillColor(isSelected ? SEC_NAME_SEL : SEC_NAME);
                name.setPosition({panelX + 22.0f, rowY + 2.0f});
                window_.draw(name);

                std::ostringstream rightOs;
                rightOs << sec.kindLabel << "  " << formatMass(sec.massSolar);
                sf::Text rightText(font_, rightOs.str(), 10);
                rightText.setFillColor(SEC_KIND);
                float rightX = panelX + PANEL_W - rightText.getLocalBounds().size.x - 14.0f;
                rightText.setPosition({rightX, rowY + 3.0f});
                window_.draw(rightText);
            } else {
                // Custom row
                std::string customLabel = "Custom BH";
                if (state.inputtingCustom) {
                    customLabel = "Custom BH  mass: " + state.customInput + "_";
                } else {
                    customLabel = "Custom BH  (press Enter, then type mass in solar masses)";
                }
                sf::Text cust(font_, customLabel, 11);
                cust.setFillColor(isSelected ? CUST_SEL : CUST_NORM);
                cust.setPosition({panelX + 22.0f, rowY + 2.0f});
                window_.draw(cust);
            }

            rowY += ROW_H;
        }

        // Footer hint
        sf::RectangleShape div2(sf::Vector2f(PANEL_W - 20.0f, 1.0f));
        div2.setPosition({panelX + 10.0f, rowY + 2.0f});
        div2.setFillColor(DIV2_COL);
        window_.draw(div2);

        sf::Text hint(font_, "Up/Down: navigate   Enter: confirm   Esc: cancel", 10);
        hint.setFillColor(HINT_COL);
        hint.setPosition({panelX + 10.0f, rowY + 6.0f});
        window_.draw(hint);
    }

    /*--------- Custom-body creator menu overlay ---------*/
    // Fields: 0=name, 1=type (cycled with Left/Right), 2=semi-major (M), 3=ecc.
    // focusedField >= 4 indexes into the saved-presets list.
    // PresetList must be iterable and each element must expose .name (string-like),
    // .semiMajorM (double), .ecc (double), .type (GalaxyBodyType-compatible).
    template<typename CustomBodyMenuState, typename PresetList, typename TypeNameFn>
    void drawCustomBodyMenu(const CustomBodyMenuState& state,
                            const char* typeLabel,
                            const PresetList& presets,
                            TypeNameFn typeNameFn,
                            float viewW, float viewH)
    {
        // Theme palette
        const bool L = lightMode_;
        const sf::Color BG_FILL    = L ? sf::Color(248, 250, 249, 235) : sf::Color(10, 12, 24, 235);
        const sf::Color BG_OUTLINE = L ? sf::Color(40, 150, 110, 220)  : sf::Color(80, 200, 160, 220);
        const sf::Color TITLE_COL  = L ? sf::Color(20, 90, 60)         : sf::Color(200, 255, 220);
        const sf::Color SUB_COL    = L ? sf::Color(60, 110, 90, 220)   : sf::Color(160, 200, 180, 200);
        const sf::Color DIV_COL    = L ? sf::Color(40, 150, 110, 110)  : sf::Color(80, 200, 160, 120);
        const sf::Color ROW_SEL    = L ? sf::Color(40, 170, 120, 65)   : sf::Color(60, 200, 160, 70);
        const sf::Color LAB_SEL    = L ? sf::Color(10, 60, 40)         : sf::Color(220, 255, 235);
        const sf::Color LAB_NORM   = L ? sf::Color(30, 50, 40)         : sf::Color(200, 220, 210);
        const sf::Color VAL_SEL    = L ? sf::Color(140, 100, 20)       : sf::Color(255, 255, 200);
        const sf::Color VAL_NORM   = L ? sf::Color(40, 80, 60)         : sf::Color(200, 230, 220);
        const sf::Color LIST_DIV   = L ? sf::Color(40, 150, 110, 90)   : sf::Color(80, 200, 160, 80);
        const sf::Color HDR_COL    = L ? sf::Color(30, 90, 70)         : sf::Color(180, 230, 210);
        const sf::Color LIST_SEL   = L ? sf::Color(40, 170, 120, 85)   : sf::Color(60, 200, 160, 90);
        const sf::Color NAME_NORM  = L ? sf::Color(30, 60, 50)         : sf::Color(210, 230, 220);
        const sf::Color RHS_NORM   = L ? sf::Color(60, 110, 90)        : sf::Color(170, 200, 185);
        const sf::Color RHS_SEL    = L ? sf::Color(120, 90, 20)        : sf::Color(220, 245, 230);
        const sf::Color HINT_COL   = L ? sf::Color(70, 110, 90, 220)   : sf::Color(140, 180, 160, 200);

        constexpr float PANEL_W      = 460.0f;
        constexpr float ROW_H        = 24.0f;
        constexpr float HEADER_H     = 46.0f;
        constexpr float FOOTER_H     = 40.0f;
        constexpr float PRESETS_HDR  = 22.0f;
        const int   nInput = 4;
        const int   nPresets = (int)presets.size();
        const int   listRows = std::min(nPresets, 8);
        const float listH    = listRows > 0 ? (PRESETS_HDR + listRows * ROW_H + 4.0f) : 0.0f;
        const float panelH   = HEADER_H + nInput * ROW_H + listH + FOOTER_H + 10.0f;
        const float panelX   = (viewW - PANEL_W) * 0.5f;
        const float panelY   = (viewH - panelH) * 0.5f;

        sf::RectangleShape bg(sf::Vector2f(PANEL_W, panelH));
        bg.setPosition({panelX, panelY});
        bg.setFillColor(BG_FILL);
        bg.setOutlineThickness(1.5f);
        bg.setOutlineColor(BG_OUTLINE);
        window_.draw(bg);

        sf::Text title(font_, "Create custom body", 14);
        title.setStyle(sf::Text::Bold);
        title.setFillColor(TITLE_COL);
        title.setPosition({panelX + 12.0f, panelY + 10.0f});
        window_.draw(title);

        sf::Text sub(font_, "Spawned around the primary BH (saved to disk)", 10);
        sub.setFillColor(SUB_COL);
        sub.setPosition({panelX + 12.0f, panelY + 28.0f});
        window_.draw(sub);

        sf::RectangleShape div(sf::Vector2f(PANEL_W - 20.0f, 1.0f));
        div.setPosition({panelX + 10.0f, panelY + HEADER_H - 2.0f});
        div.setFillColor(DIV_COL);
        window_.draw(div);

        auto drawRow = [&](float y, int idx, const char* label, const std::string& value,
                           bool isType, bool isText) {
            const bool focused = (state.focusedField == idx);
            if (focused) {
                sf::RectangleShape sel(sf::Vector2f(PANEL_W - 12.0f, ROW_H - 4.0f));
                sel.setPosition({panelX + 6.0f, y});
                sel.setFillColor(ROW_SEL);
                window_.draw(sel);
            }
            sf::Text lab(font_, label, 11);
            lab.setFillColor(focused ? LAB_SEL : LAB_NORM);
            lab.setPosition({panelX + 14.0f, y + 4.0f});
            window_.draw(lab);

            std::string shown = value;
            if (focused && isText) shown += "_";
            if (focused && isType) shown = std::string("< ") + value + " >";
            if (shown.empty() && isText && !focused) shown = "(auto)";
            sf::Text val(font_, shown, 11);
            val.setFillColor(focused ? VAL_SEL : VAL_NORM);
            float vx = panelX + PANEL_W - val.getLocalBounds().size.x - 16.0f;
            val.setPosition({vx, y + 4.0f});
            window_.draw(val);
        };

        const float inputBaseY = panelY + HEADER_H;
        drawRow(inputBaseY + 0 * ROW_H + 2.0f, 0, "Name",                state.nameInput, false, true);
        drawRow(inputBaseY + 1 * ROW_H + 2.0f, 1, "Type",                typeLabel,       true,  false);
        drawRow(inputBaseY + 2 * ROW_H + 2.0f, 2, "Semi-major axis (M)", state.smInput,   false, true);
        drawRow(inputBaseY + 3 * ROW_H + 2.0f, 3, "Eccentricity",        state.eccInput,  false, true);

        float listY = inputBaseY + nInput * ROW_H + 6.0f;
        if (listRows > 0) {
            sf::RectangleShape divL(sf::Vector2f(PANEL_W - 20.0f, 1.0f));
            divL.setPosition({panelX + 10.0f, listY});
            divL.setFillColor(LIST_DIV);
            window_.draw(divL);

            std::ostringstream hdr;
            hdr << "Saved presets (" << nPresets << ")";
            sf::Text hdrT(font_, hdr.str(), 11);
            hdrT.setStyle(sf::Text::Bold);
            hdrT.setFillColor(HDR_COL);
            hdrT.setPosition({panelX + 12.0f, listY + 4.0f});
            window_.draw(hdrT);

            // Scroll window: anchor so focused row stays visible.
            int focusedPreset = state.focusedField - nInput;
            int firstShown = 0;
            if (focusedPreset >= 0 && nPresets > listRows) {
                firstShown = std::clamp(focusedPreset - listRows / 2,
                                        0, nPresets - listRows);
            }

            for (int i = 0; i < listRows; ++i) {
                int presetIdx = firstShown + i;
                if (presetIdx >= nPresets) break;
                const auto& p = presets[presetIdx];
                const float y = listY + PRESETS_HDR + i * ROW_H;
                int rowField = nInput + presetIdx;
                bool focused = (state.focusedField == rowField);
                if (focused) {
                    sf::RectangleShape sel(sf::Vector2f(PANEL_W - 12.0f, ROW_H - 4.0f));
                    sel.setPosition({panelX + 6.0f, y});
                    sel.setFillColor(LIST_SEL);
                    window_.draw(sel);
                }
                sf::Text nameT(font_, p.name, 11);
                nameT.setFillColor(focused ? VAL_SEL : NAME_NORM);
                nameT.setPosition({panelX + 14.0f, y + 4.0f});
                window_.draw(nameT);

                std::ostringstream rhs;
                rhs << typeNameFn(p.type)
                    << "  a=" << std::fixed << std::setprecision(1) << p.semiMajorM
                    << "M  e=" << std::setprecision(2) << p.ecc;
                sf::Text rhsT(font_, rhs.str(), 10);
                rhsT.setFillColor(focused ? RHS_SEL : RHS_NORM);
                float vx = panelX + PANEL_W - rhsT.getLocalBounds().size.x - 16.0f;
                rhsT.setPosition({vx, y + 5.0f});
                window_.draw(rhsT);
            }
        }

        float footY = panelY + panelH - FOOTER_H + 6.0f;
        sf::RectangleShape div2(sf::Vector2f(PANEL_W - 20.0f, 1.0f));
        div2.setPosition({panelX + 10.0f, footY});
        div2.setFillColor(LIST_DIV);
        window_.draw(div2);

        const bool onPreset = (state.focusedField >= nInput);
        const char* hintTxt = onPreset
            ? "Up/Down: row   Enter: spawn preset   Del: remove   Esc: cancel"
            : "Tab/Up/Down: field   Left/Right: cycle type   Enter: spawn + save   Esc: cancel";
        sf::Text hint(font_, hintTxt, 10);
        hint.setFillColor(HINT_COL);
        hint.setPosition({panelX + 10.0f, footY + 6.0f});
        window_.draw(hint);
    }

    /*--------- Tidal disruption event overlay ---------*/
    // Caller pre-projects world positions to screen coords.
    struct TidalParticleVis {
        float sx, sy;       // screen position
        float size;         // dot radius in pixels (already lifetime-scaled)
        float lifeF;        // lifetime fraction [0,1]
        bool  isFallback;   // true = gas-stream (blue-white), false = debris (orange)
    };

    // flashPos: screen coords of disruption point
    // flashAlpha: 0 (no flash) to 1 (peak flash), derived from flashTimer/FLASH_DURATION
    // label: short overlay text (e.g. "TIDAL DISRUPTION EVENT", "MERGER REMNANT").
    //        If null, no label is drawn.
    void drawTidalEvent(sf::Vector2f flashPos, float flashAlpha,
                        const std::vector<TidalParticleVis>& particles,
                        bool showLabel,
                        const char* label = "TIDAL DISRUPTION EVENT") {
        // Localised expanding ring flash
        if (flashAlpha > 0.0f) {
            float t = 1.0f - flashAlpha;  // 0 = just triggered, 1 = faded
            // Inner bright core
            float coreR = flashAlpha * 22.0f;
            sf::CircleShape core(coreR);
            core.setOrigin({coreR, coreR});
            core.setPosition(flashPos);
            core.setFillColor(sf::Color(255, 230, 160,
                                        static_cast<uint8_t>(flashAlpha * 210)));
            window_.draw(core);

            // Expanding shockwave ring
            float ringR = 20.0f + t * 90.0f;
            sf::CircleShape ring(ringR);
            ring.setOrigin({ringR, ringR});
            ring.setPosition(flashPos);
            ring.setFillColor(sf::Color::Transparent);
            ring.setOutlineThickness(2.5f);
            ring.setOutlineColor(sf::Color(255, 160, 60,
                                           static_cast<uint8_t>(flashAlpha * 180)));
            window_.draw(ring);
        }

        // Particle cloud
        for (const auto& p : particles) {
            if (p.size < 0.5f) continue;
            sf::CircleShape dot(p.size);
            dot.setOrigin({p.size, p.size});
            dot.setPosition({p.sx, p.sy});

            if (p.isFallback) {
                // Gas stream: blue-white
                dot.setFillColor(sf::Color(140, 195, 255,
                                           static_cast<uint8_t>(p.lifeF * 195)));
            } else {
                // Debris: orange → yellow fade as it cools
                uint8_t g = static_cast<uint8_t>(80 + p.lifeF * 120);
                dot.setFillColor(sf::Color(255, g, 30,
                                           static_cast<uint8_t>(p.lifeF * 215)));
            }
            window_.draw(dot);
        }

        // Notification label near the disruption point
        if (showLabel && label && label[0]) {
            sf::Text lbl(font_, label, 13);
            lbl.setFillColor(sf::Color(255, 140, 40, 220));
            lbl.setPosition({flashPos.x - 110.0f, flashPos.y + 32.0f});
            window_.draw(lbl);
        }
    }

    /*--------- Reset view after window resize ---------*/
    void resetView(float viewW, float viewH) {
        sf::View view(sf::FloatRect({0.f, 0.f}, {viewW, viewH}));
        window_.setView(view);
    }
};
