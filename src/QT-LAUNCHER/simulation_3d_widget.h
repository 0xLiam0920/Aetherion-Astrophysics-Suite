#ifndef SIMULATION_3D_WIDGET_H
#define SIMULATION_3D_WIDGET_H

#include "sfml_canvas.h"
#include "custom_bh_dialog.h"  // for CustomBH3DConfig
#include <memory>

#ifdef AETHERION_QT_HOST
#include <QKeyEvent>
#include <QMouseEvent>
#include <QWheelEvent>
#endif

// Pimpl: all GL / simulation state lives in this struct, defined in the .cpp.
// This keeps all GL headers out of the Qt header.
struct Sim3DState;

/**
 * Simulation3DWidget — embeds the full 3D ray-marching black hole simulation
 * inside a Qt tab. Uses the QSFMLCanvas OpenGL context; all GL objects are
 * created in onInit() after the context is active.
 */
class Simulation3DWidget : public QSFMLCanvas
{
    Q_OBJECT

public:
    explicit Simulation3DWidget(QWidget *parent = nullptr);
    ~Simulation3DWidget();   // Non-inline: defined in .cpp where Sim3DState is complete

    /// Store a custom config to be applied once onInit() runs (call before the widget is shown).
    void setPendingConfig(const CustomBH3DConfig &cfg);

protected:
    void onInit() override;
    void onUpdate() override;
    void onKeyPressed(sf::Keyboard::Key code) override;
    void onKeyReleased(sf::Keyboard::Key code) override;
    void onMouseMoved(float x, float y) override;
    void onMousePressed(sf::Mouse::Button button, float x, float y) override;
    void onMouseReleased(sf::Mouse::Button button, float x, float y) override;
    void onMouseWheelScrolled(float delta, float x, float y) override;

#ifdef AETHERION_QT_HOST
    // Override Qt-level events to feed raw input into ImGui before the
    // QSFMLCanvas base converts them to SFML key/button codes.
    void keyPressEvent(QKeyEvent* ev) override;
    void keyReleaseEvent(QKeyEvent* ev) override;
    void mouseMoveEvent(QMouseEvent* ev) override;
    void mousePressEvent(QMouseEvent* ev) override;
    void mouseReleaseEvent(QMouseEvent* ev) override;
    void wheelEvent(QWheelEvent* ev) override;
#endif

private:
    std::unique_ptr<Sim3DState> state_;
    bool              hasPendingConfig_ = false;
    CustomBH3DConfig  pendingConfig_;

#ifdef AETHERION_QT_HOST
    bool imguiReady_ = false;   ///< true once imgui_qt::init() succeeded in onInit()
#endif
};

#endif // SIMULATION_3D_WIDGET_H
