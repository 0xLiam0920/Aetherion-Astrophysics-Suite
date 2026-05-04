// ============================================================
// simulation_3d_widget.cpp — Qt-embedded 3D black-hole widget.
// ============================================================
// Thin glue around the shared core (src/3D/bh3d_core.hpp). The
// QSFMLCanvas owns the GL context / event source; everything
// simulation- or render-related is delegated to bh3d::*.
// ============================================================

#include "simulation_3d_widget.h"
#include "custom_bh_dialog.h"

// All GL / SFML / simulation headers; never include from the .h
#include "platform.hpp"
#include <SFML/Graphics/Font.hpp>

#include <algorithm>
#include <iostream>
#include <memory>

#include "bh3d_core.hpp"

#ifdef AETHERION_QT_HOST
#include "hud_imgui.hpp"
#include "imgui_qt_adapter.hpp"
#endif

// ─────────────────────────────────────────────────────────────
// Pimpl: Sim3DState wraps the shared bh3d::State plus widget-only
// extras (resource manager + mouse-look state for Qt event source).
// ─────────────────────────────────────────────────────────────
struct Sim3DState {
    bh3d::State     core;
    ResourceManager resources;

    explicit Sim3DState(const cfg::CameraConfig& camCfg) : core(camCfg) {}
};

// ─────────────────────────────────────────────────────────────
// Widget lifecycle
// ─────────────────────────────────────────────────────────────

Simulation3DWidget::Simulation3DWidget(QWidget *parent)
    : QSFMLCanvas(parent, QSize(1000, 700),
                  []{ sf::ContextSettings s;
                      s.depthBits = 24; s.stencilBits = 8;
                      s.majorVersion = 3; s.minorVersion = 3;
                      s.attributeFlags = sf::ContextSettings::Core; return s; }())
{}

Simulation3DWidget::~Simulation3DWidget()
{
#ifdef AETHERION_QT_HOST
    if (imguiReady_) {
        (void)setActive(true);
        imgui_qt::shutdown();
        imguiReady_ = false;
    }
#endif
}

void Simulation3DWidget::setPendingConfig(const CustomBH3DConfig &cfg)
{
    pendingConfig_    = cfg;
    hasPendingConfig_ = true;
}

void Simulation3DWidget::onInit()
{
    (void)setActive(true);
    platformInitGL();

    // Build pimpl with the camera config from the default profile so the
    // CameraController can be constructed; bh3d::initProfiles will reset it.
    {
        auto probe = profiles::allProfiles();
        state_ = std::make_unique<Sim3DState>(probe[0].config.camera);
    }
    auto& s = state_->core;

    bh3d::initProfiles(s, /*initialIdx=*/0); // start on TON 618

    // Apply custom config if the caller pushed one before show()
    if (hasPendingConfig_) {
        const auto& c = pendingConfig_;
        s.config.blackHole.spinParameter = c.spinParameter;
        s.config.disk.innerRadius        = c.diskInnerRadius;
        s.config.disk.outerRadius        = c.diskOuterRadius;
        s.config.disk.halfThickness      = c.diskHalfThickness;
        s.config.disk.peakTemp           = c.diskPeakTemp;
        s.config.jet.length              = c.jetLength;
        s.jetsEnabled                    = c.jetsEnabled;
        s.blrEnabled                     = c.blrEnabled;
        s.dopplerEnabled                 = c.dopplerEnabled;
        s.hostGalaxyEnabled              = c.hostGalaxyEnabled;
        hasPendingConfig_ = false;
    }

    // Viewport sized to the actual widget
    const auto sz = getSize();
    const int w = std::max((int)sz.x, 1);
    const int h = std::max((int)sz.y, 1);
    glViewport(0, 0, w, h);

    // Font
    {
        sf::Font tmpFont;
        bh3d::initFont(s, tmpFont);
        if (!s.fontLoaded)
            std::cerr << "[3D widget] Warning: Could not load font for HUD.\n";
    }
    (void)setActive(true);

    // Shaders / textures / bloom
    if (!bh3d::initShaders(s, state_->resources)) {
        std::cerr << "[3D widget] FATAL: no fragment shader compiled.\n";
        return;
    }
    s.quad = GLVertexArray::makeFullScreenQuad();
    bh3d::initTextures(s, state_->resources);

    if (!s.bloom.init(w, h))
        std::cerr << "[3D widget] bloom init failed\n";
    s.lastW = w; s.lastH = h;

    s.actionKeys = bh3d::loadActionKeybinds();
    s.frameClock.restart();
    s.fpsClock.restart();

#ifdef AETHERION_QT_HOST
    imgui_qt::init(w, h);
    s.useImGuiHud = true;
    imguiReady_   = true;
#endif
}

// ─────────────────────────────────────────────────────────────
// Per-frame update (paintEvent calls onUpdate, then display())
// ─────────────────────────────────────────────────────────────
void Simulation3DWidget::onUpdate()
{
    if (!state_ || !state_->core.activeProgram) return;
    auto& s = state_->core;
    (void)setActive(true);

    const float dt = s.frameClock.restart().asSeconds();
    bh3d::tickFPS(s);
    bh3d::tickPhysics(s, dt);

    const auto sz = getSize();
    const int w = std::max((int)sz.x, 1);
    const int h = std::max((int)sz.y, 1);

    bh3d::buildSnapshot(s, w, h, dt);
    bh3d::renderScene(s, w, h);

#ifdef AETHERION_QT_HOST
    if (imguiReady_) {
        bh3d::renderHUD(s, w, h);            // legacy bitmap panels (gated by !useImGuiHud)
        imgui_qt::newFrame(w, h, dt);
        hud_im::drawAll(s);
        imgui_qt::render();
    } else {
        bh3d::renderHUD(s, w, h);
    }
#else
    bh3d::renderHUD(s, w, h);
#endif

    while (glGetError() != GL_NO_ERROR) {}
    // Note: paintEvent calls display() after onUpdate(); do not call here.
}

// ─────────────────────────────────────────────────────────────
// Keyboard
// ─────────────────────────────────────────────────────────────
void Simulation3DWidget::onKeyPressed(sf::Keyboard::Key code)
{
    if (!state_) return;
    auto& s = state_->core;
    s.keys.onKeyPressed(code);
    bh3d::onActionKey(s, code);
}

void Simulation3DWidget::onKeyReleased(sf::Keyboard::Key code)
{
    if (!state_) return;
    state_->core.keys.onKeyReleased(code);
}

// ─────────────────────────────────────────────────────────────
// Mouse
// ─────────────────────────────────────────────────────────────
void Simulation3DWidget::onMouseMoved(float x, float y)
{
    if (!state_) return;
    auto& s = state_->core;
    if (s.looking) {
        const float dx = x - s.lastMouseX;
        const float dy = y - s.lastMouseY;
        s.camera.onMouseDelta(dx, dy);
    }
    s.lastMouseX = x;
    s.lastMouseY = y;
}

void Simulation3DWidget::onMousePressed(sf::Mouse::Button button, float x, float y)
{
    if (!state_) return;
    auto& s = state_->core;
    if (button == sf::Mouse::Button::Right) {
        s.looking    = true;
        s.lastMouseX = x;
        s.lastMouseY = y;
    } else if (button == sf::Mouse::Button::Left) {
        bh3d::onLeftClickOverlays(s, x, y);
    }
}

void Simulation3DWidget::onMouseReleased(sf::Mouse::Button button, float /*x*/, float /*y*/)
{
    if (!state_) return;
    if (button == sf::Mouse::Button::Right)
        state_->core.looking = false;
}

void Simulation3DWidget::onMouseWheelScrolled(float /*delta*/, float /*x*/, float /*y*/)
{
    // (Intentionally no-op: legacy widget left zoom unimplemented.)
}

// ─────────────────────────────────────────────────────────────
// Qt-level event overrides (AETHERION_QT_HOST)
// Feed raw Qt events into ImGui before the QSFMLCanvas base class
// converts them to SFML types and calls the onXxx() virtuals above.
// ─────────────────────────────────────────────────────────────
#ifdef AETHERION_QT_HOST

void Simulation3DWidget::keyPressEvent(QKeyEvent* ev)
{
    if (imguiReady_) imgui_qt::feedKeyPress(ev);
    // Forward to simulation only when ImGui is not consuming keyboard input.
    // (WantCaptureKeyboard reflects the previous frame — standard ImGui practice.)
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureKeyboard)
        QSFMLCanvas::keyPressEvent(ev);
}

void Simulation3DWidget::keyReleaseEvent(QKeyEvent* ev)
{
    if (imguiReady_) imgui_qt::feedKeyRelease(ev);
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureKeyboard)
        QSFMLCanvas::keyReleaseEvent(ev);
}

void Simulation3DWidget::mouseMoveEvent(QMouseEvent* ev)
{
    if (imguiReady_) imgui_qt::feedMouseMove(ev);
    // Always forward mouse moves — camera look uses drag (right-button held),
    // and WantCaptureMouse only covers ImGui windows.
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureMouse)
        QSFMLCanvas::mouseMoveEvent(ev);
    else
        QWidget::mouseMoveEvent(ev);   // skip SFML delta-look, but let Qt propagate
}

void Simulation3DWidget::mousePressEvent(QMouseEvent* ev)
{
    if (imguiReady_) imgui_qt::feedMousePress(ev);
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureMouse)
        QSFMLCanvas::mousePressEvent(ev);
    else
        QWidget::mousePressEvent(ev);
}

void Simulation3DWidget::mouseReleaseEvent(QMouseEvent* ev)
{
    if (imguiReady_) imgui_qt::feedMouseRelease(ev);
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureMouse)
        QSFMLCanvas::mouseReleaseEvent(ev);
    else
        QWidget::mouseReleaseEvent(ev);
}

void Simulation3DWidget::wheelEvent(QWheelEvent* ev)
{
    if (imguiReady_) imgui_qt::feedWheel(ev);
    if (!imguiReady_ || !ImGui::GetIO().WantCaptureMouse)
        QSFMLCanvas::wheelEvent(ev);
    else
        QWidget::wheelEvent(ev);
}

#endif // AETHERION_QT_HOST
