// imgui_qt_adapter.hpp: Qt event → ImGuiIO bridge + imgui_impl_opengl3 lifecycle.
//
// This header is only compiled into translation units that define AETHERION_QT_HOST.
// It provides a thin `imgui_qt` namespace that:
//   - initialises / shuts down an ImGui context using the raw OpenGL3 backend
//     (no SFML involved: the Qt widget owns the GL context directly)
//   - translates raw Qt key/mouse/wheel events into ImGui AddXxxEvent() calls
//     before they are converted to SFML keys by QSFMLCanvas
//
// Usage pattern in Simulation3DWidget:
//   onInit()    → imgui_qt::init(w, h)
//   onUpdate()  → imgui_qt::newFrame(w, h, dt)
//              → hud_im::drawAll(state)
//              → imgui_qt::render()
//   destructor  → imgui_qt::shutdown()
//   Qt events   → imgui_qt::feedKeyPress/Release/MouseMove/Press/Release/Wheel(ev)
//                 (call BEFORE forwarding to the base QSFMLCanvas handler so that
//                  WantCaptureKeyboard / WantCaptureMouse reflect the previous frame)
//
// NOTE: ImGui::GetIO().WantCaptureKeyboard / WantCaptureMouse from the **previous**
// frame is used to gate whether Qt events are forwarded to the simulation.
#pragma once
#ifdef AETHERION_QT_HOST

#include <QKeyEvent>
#include <QMouseEvent>
#include <QWheelEvent>

#include "imgui.h"
#include "imgui_impl_opengl3.h"
#include "hud_imgui.hpp"   // hud_im::applyAetherionStyle()

namespace imgui_qt {

// ─────────────────────────────────────────────────────────────────────────────
// Qt key → ImGuiKey mapping
// ─────────────────────────────────────────────────────────────────────────────
inline ImGuiKey qtKeyToImGui(int k)
{
    // Letter keys A–Z
    if (k >= Qt::Key_A && k <= Qt::Key_Z)
        return static_cast<ImGuiKey>(ImGuiKey_A + (k - Qt::Key_A));
    // Digit keys 0–9
    if (k >= Qt::Key_0 && k <= Qt::Key_9)
        return static_cast<ImGuiKey>(ImGuiKey_0 + (k - Qt::Key_0));
    // Function keys F1–F12
    if (k >= Qt::Key_F1 && k <= Qt::Key_F12)
        return static_cast<ImGuiKey>(ImGuiKey_F1 + (k - Qt::Key_F1));

    switch (k) {
        case Qt::Key_Return:       return ImGuiKey_Enter;
        case Qt::Key_Enter:        return ImGuiKey_KeypadEnter;
        case Qt::Key_Escape:       return ImGuiKey_Escape;
        case Qt::Key_Backspace:    return ImGuiKey_Backspace;
        case Qt::Key_Delete:       return ImGuiKey_Delete;
        case Qt::Key_Tab:          return ImGuiKey_Tab;
        case Qt::Key_Space:        return ImGuiKey_Space;
        case Qt::Key_Left:         return ImGuiKey_LeftArrow;
        case Qt::Key_Right:        return ImGuiKey_RightArrow;
        case Qt::Key_Up:           return ImGuiKey_UpArrow;
        case Qt::Key_Down:         return ImGuiKey_DownArrow;
        case Qt::Key_Home:         return ImGuiKey_Home;
        case Qt::Key_End:          return ImGuiKey_End;
        case Qt::Key_PageUp:       return ImGuiKey_PageUp;
        case Qt::Key_PageDown:     return ImGuiKey_PageDown;
        case Qt::Key_Insert:       return ImGuiKey_Insert;
        case Qt::Key_Shift:
        case Qt::Key_CapsLock:     return ImGuiKey_LeftShift;
        case Qt::Key_Control:      return ImGuiKey_LeftCtrl;
        case Qt::Key_Alt:          return ImGuiKey_LeftAlt;
        case Qt::Key_Meta:         return ImGuiKey_LeftSuper;
        case Qt::Key_BracketLeft:  return ImGuiKey_LeftBracket;
        case Qt::Key_BracketRight: return ImGuiKey_RightBracket;
        case Qt::Key_Semicolon:    return ImGuiKey_Semicolon;
        case Qt::Key_Apostrophe:   return ImGuiKey_Apostrophe;
        case Qt::Key_Comma:        return ImGuiKey_Comma;
        case Qt::Key_Period:       return ImGuiKey_Period;
        case Qt::Key_Slash:        return ImGuiKey_Slash;
        case Qt::Key_Question:     return ImGuiKey_Slash;  // Shift+/ on most layouts
        case Qt::Key_Backslash:    return ImGuiKey_Backslash;
        case Qt::Key_Minus:        return ImGuiKey_Minus;
        case Qt::Key_Equal:        return ImGuiKey_Equal;
        case Qt::Key_QuoteLeft:    return ImGuiKey_GraveAccent;
        default:                   return ImGuiKey_None;
    }
}

inline void applyModifiers(Qt::KeyboardModifiers mods)
{
    ImGuiIO& io = ImGui::GetIO();
    io.AddKeyEvent(ImGuiMod_Ctrl,  (mods & Qt::ControlModifier) != 0);
    io.AddKeyEvent(ImGuiMod_Shift, (mods & Qt::ShiftModifier)   != 0);
    io.AddKeyEvent(ImGuiMod_Alt,   (mods & Qt::AltModifier)     != 0);
    io.AddKeyEvent(ImGuiMod_Super, (mods & Qt::MetaModifier)    != 0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Lifecycle, call inside the widget's GL context
// ─────────────────────────────────────────────────────────────────────────────

/// Call once after the GL context is active (e.g. end of onInit()).
/// @param w,h  Initial widget pixel dimensions.
inline void init(int w, int h)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.DisplaySize  = ImVec2(static_cast<float>(w), static_cast<float>(h));
    io.DeltaTime    = 1.0f / 60.0f;
    // Disable the default imgui.ini file, this widget is embedded in Qt.
    io.IniFilename  = nullptr;
    hud_im::applyAetherionStyle();
    ImGui_ImplOpenGL3_Init("#version 330 core");
}

/// Call with the GL context active (e.g. widget destructor after setActive(true)).
inline void shutdown()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui::DestroyContext();
}

// ─────────────────────────────────────────────────────────────────────────────
// Per-frame, call inside onUpdate() with the GL context active
// ─────────────────────────────────────────────────────────────────────────────

/// Begin a new ImGui frame.  Call before any hud_im::drawXxx calls.
/// @param w,h  Current widget pixel dimensions (re-read each frame).
/// @param dt   Elapsed seconds since last frame (from the simulation clock).
inline void newFrame(int w, int h, float dt)
{
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(static_cast<float>(w), static_cast<float>(h));
    io.DeltaTime   = (dt > 0.0f) ? dt : (1.0f / 60.0f);
    ImGui_ImplOpenGL3_NewFrame();
    ImGui::NewFrame();
}

/// Finalise and render the ImGui draw data.  Call after all hud_im::drawXxx calls.
inline void render()
{
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// ─────────────────────────────────────────────────────────────────────────────
// Event feeds, call from the widget's Qt event overrides
// ─────────────────────────────────────────────────────────────────────────────

inline void feedKeyPress(QKeyEvent* ev)
{
    if (!ev) return;
    ImGuiIO& io = ImGui::GetIO();
    applyModifiers(ev->modifiers());
    const ImGuiKey key = qtKeyToImGui(ev->key());
    if (key != ImGuiKey_None)
        io.AddKeyEvent(key, true);
    // Feed text input for printable characters (enables text fields, search boxes, etc.)
    const QString text = ev->text();
    for (const QChar& ch : text)
        if (ch.isPrint())
            io.AddInputCharacter(static_cast<unsigned int>(ch.unicode()));
}

inline void feedKeyRelease(QKeyEvent* ev)
{
    if (!ev) return;
    applyModifiers(ev->modifiers());
    const ImGuiKey key = qtKeyToImGui(ev->key());
    if (key != ImGuiKey_None)
        ImGui::GetIO().AddKeyEvent(key, false);
}

inline void feedMouseMove(QMouseEvent* ev)
{
    if (!ev) return;
    ImGui::GetIO().AddMousePosEvent(
        static_cast<float>(ev->pos().x()),
        static_cast<float>(ev->pos().y()));
}

inline void feedMousePress(QMouseEvent* ev)
{
    if (!ev) return;
    ImGuiMouseButton btn = ImGuiMouseButton_Left;
    if (ev->button() == Qt::RightButton)  btn = ImGuiMouseButton_Right;
    if (ev->button() == Qt::MiddleButton) btn = ImGuiMouseButton_Middle;
    ImGui::GetIO().AddMouseButtonEvent(btn, true);
}

inline void feedMouseRelease(QMouseEvent* ev)
{
    if (!ev) return;
    ImGuiMouseButton btn = ImGuiMouseButton_Left;
    if (ev->button() == Qt::RightButton)  btn = ImGuiMouseButton_Right;
    if (ev->button() == Qt::MiddleButton) btn = ImGuiMouseButton_Middle;
    ImGui::GetIO().AddMouseButtonEvent(btn, false);
}

/// @param ev  Qt wheel event (angleDelta().y() in 1/8° units; 120 = one notch).
inline void feedWheel(QWheelEvent* ev)
{
    if (!ev) return;
    ImGui::GetIO().AddMouseWheelEvent(
        static_cast<float>(ev->angleDelta().x()) / 120.0f,
        static_cast<float>(ev->angleDelta().y()) / 120.0f);
}

} // namespace imgui_qt
#endif // AETHERION_QT_HOST
