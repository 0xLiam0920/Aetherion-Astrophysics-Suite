#include "bh_preview_widget.h"

#include <QPainter>
#include <QPainterPath>
#include <QRadialGradient>
#include <QLinearGradient>
#include <cmath>
#include <random>

static constexpr int   TIMER_MS = 30;
static constexpr float PI       = 3.14159265f;

BlackHolePreviewWidget::BlackHolePreviewWidget(QWidget *parent)
    : QWidget(parent)
{
    setMinimumHeight(180);
    setMaximumHeight(220);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    setStyleSheet("background: transparent;");
    seedStars();
    m_timer = new QTimer(this);
    connect(m_timer, &QTimer::timeout, this, &BlackHolePreviewWidget::tick);
}

void BlackHolePreviewWidget::setStyle(PreviewStyle s)
{
    m_style = s;
    if (s == PreviewStyle::NONE) {
        m_timer->stop();
        setVisible(false);
    } else {
        setVisible(true);
        m_timer->start(TIMER_MS);
    }
    update();
}

void BlackHolePreviewWidget::setStyle(const QString &key)
{
    if      (key == "sgra")             setStyle(PreviewStyle::SGR_A);
    else if (key == "m87")              setStyle(PreviewStyle::M87);
    else if (key == "ton618")           setStyle(PreviewStyle::TON618);
    else if (key == "3c273")            setStyle(PreviewStyle::C3_273);
    else if (key == "j0529")            setStyle(PreviewStyle::J0529);
    else if (key == "cygnusx1")         setStyle(PreviewStyle::CYGNUS_X1);
    else if (key == "gw150914")         setStyle(PreviewStyle::GW150914);
    else if (key == "gaiabh1")          setStyle(PreviewStyle::GAIA_BH1);
    else if (key == "gaiabh2")          setStyle(PreviewStyle::GAIA_BH2);
    else if (key == "gaiabh3")          setStyle(PreviewStyle::GAIA_BH3);
    else if (key == "v404cygni")        setStyle(PreviewStyle::V404_CYGNI);
    else if (key == "phoenixa")         setStyle(PreviewStyle::PHOENIX_A);
    // Research Scenarios
    else if (key == "isco_unstable")    setStyle(PreviewStyle::ISCO_UNSTABLE);
    else if (key == "isco_critical")    setStyle(PreviewStyle::ISCO_CRITICAL);
    else if (key == "isco_stable")      setStyle(PreviewStyle::ISCO_STABLE);
    else if (key == "photon_sphere")    setStyle(PreviewStyle::PHOTON_SPHERE);
    else if (key == "radial_infall")    setStyle(PreviewStyle::RADIAL_INFALL);
    else if (key == "tidal_disrupt")    setStyle(PreviewStyle::TIDAL_DISRUPTION);
    else if (key == "pulsar_inspiral")  setStyle(PreviewStyle::PULSAR_INSPIRAL);
    // Sgr A* system
    else if (key == "sgra_s2")          setStyle(PreviewStyle::SGRA_S2);
    else if (key == "sgra_s14")         setStyle(PreviewStyle::SGRA_S14);
    else if (key == "sgra_irs16")       setStyle(PreviewStyle::SGRA_IRS16);
    else if (key == "sgra_gasclump")    setStyle(PreviewStyle::SGRA_GASCLUMP);
    else if (key == "sgra_smember")     setStyle(PreviewStyle::SGRA_SMEMBER);
    // TON 618 system
    else if (key == "ton618_blr")       setStyle(PreviewStyle::TON618_BLR);
    else if (key == "ton618_tdstar")    setStyle(PreviewStyle::TON618_TDSTAR);
    else if (key == "ton618_gasfilm")   setStyle(PreviewStyle::TON618_GASFILM);
    else if (key == "ton618_plunge")    setStyle(PreviewStyle::TON618_PLUNGE);
    else if (key == "ton618_outerstar") setStyle(PreviewStyle::TON618_OUTERSTAR);
    else if (key == "ton618_cluster")   setStyle(PreviewStyle::TON618_CLUSTER);
    // 3C 273 system
    else if (key == "3c273_jetblob")    setStyle(PreviewStyle::C3273_JETBLOB);
    else if (key == "3c273_blr")        setStyle(PreviewStyle::C3273_BLR);
    else if (key == "3c273_dwarf")      setStyle(PreviewStyle::C3273_DWARF);
    else if (key == "3c273_closestar")  setStyle(PreviewStyle::C3273_CLOSESTAR);
    // J0529-4351 system
    else if (key == "j0529_fastblob")   setStyle(PreviewStyle::J0529_FASTBLOB);
    else if (key == "j0529_uvclump")    setStyle(PreviewStyle::J0529_UVCLUMP);
    else if (key == "j0529_tdestar")    setStyle(PreviewStyle::J0529_TDESTAR);
    else if (key == "j0529_gasstream")  setStyle(PreviewStyle::J0529_GASSTREAM);
    else if (key == "j0529_cluster")    setStyle(PreviewStyle::J0529_CLUSTER);
    // M87 system
    else if (key == "m87_jetknot")      setStyle(PreviewStyle::M87_JETKNOT);
    else if (key == "m87_innerstar")    setStyle(PreviewStyle::M87_INNERSTAR);
    else if (key == "m87_hotgas")       setStyle(PreviewStyle::M87_HOTGAS);
    else if (key == "m87_globular")     setStyle(PreviewStyle::M87_GLOBULAR);
    else if (key == "m87_dwarf")        setStyle(PreviewStyle::M87_DWARF);
    else                                setStyle(PreviewStyle::NONE);
}

void BlackHolePreviewWidget::tick()
{
    m_angle += 0.4f;
    if (m_angle >= 360.f) m_angle -= 360.f;
    m_phase += 0.04f;
    ++m_frame;
    update();
}

void BlackHolePreviewWidget::seedStars()
{
    auto fill = [](QVector<Star> &v, int n, unsigned seed, float maxSz) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<float> pos(0.f, 1.f);
        std::uniform_real_distribution<float> sz (0.4f, maxSz);
        std::uniform_real_distribution<float> br (0.25f, 1.f);
        v.resize(n);
        for (auto &s : v) {
            s.x = pos(rng); s.y = pos(rng);
            s.r = sz(rng);  s.brightness = br(rng);
        }
    };
    fill(m_stars,   80,  42u, 1.8f);
    fill(m_gcStars, 160, 99u, 1.4f);
}

void BlackHolePreviewWidget::paintEvent(QPaintEvent *)
{
    if (m_style == PreviewStyle::NONE) return;
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    switch (m_style) {
    case PreviewStyle::SGR_A:             paintSgrA(p);            break;
    case PreviewStyle::M87:               paintM87(p);             break;
    case PreviewStyle::TON618:            paintTon618(p);          break;
    case PreviewStyle::C3_273:            paint3c273(p);           break;
    case PreviewStyle::J0529:             paintJ0529(p);           break;
    case PreviewStyle::CYGNUS_X1:         paintCygX1(p);           break;
    case PreviewStyle::GW150914:          paintGW(p);              break;
    case PreviewStyle::GAIA_BH1:          paintGaiaBH1(p);         break;
    case PreviewStyle::GAIA_BH2:          paintGaiaBH2(p);         break;
    case PreviewStyle::GAIA_BH3:          paintGaiaBH3(p);         break;
    case PreviewStyle::V404_CYGNI:        paintV404Cygni(p);       break;
    case PreviewStyle::PHOENIX_A:         paintPhoenixA(p);        break;
    // Research Scenarios
    case PreviewStyle::ISCO_UNSTABLE:     paintIscoUnstable(p);   break;
    case PreviewStyle::ISCO_CRITICAL:     paintIscoCritical(p);   break;
    case PreviewStyle::ISCO_STABLE:       paintIscoStable(p);     break;
    case PreviewStyle::PHOTON_SPHERE:     paintPhotonSphere(p);   break;
    case PreviewStyle::RADIAL_INFALL:     paintRadialInfall(p);   break;
    case PreviewStyle::TIDAL_DISRUPTION:  paintTidalDisrupt(p);   break;
    case PreviewStyle::PULSAR_INSPIRAL:   paintPulsarInspiral(p); break;
    // Sgr A* system
    case PreviewStyle::SGRA_S2:           paintSgraS2(p);         break;
    case PreviewStyle::SGRA_S14:          paintSgraS14(p);        break;
    case PreviewStyle::SGRA_IRS16:        paintSgraIrs16(p);      break;
    case PreviewStyle::SGRA_GASCLUMP:     paintSgraGasClump(p);   break;
    case PreviewStyle::SGRA_SMEMBER:      paintSgraSMember(p);    break;
    // TON 618 system
    case PreviewStyle::TON618_BLR:        paintTon618Blr(p);      break;
    case PreviewStyle::TON618_TDSTAR:     paintTon618TdStar(p);   break;
    case PreviewStyle::TON618_GASFILM:    paintTon618GasFilm(p);  break;
    case PreviewStyle::TON618_PLUNGE:     paintTon618Plunge(p);   break;
    case PreviewStyle::TON618_OUTERSTAR:  paintTon618OuterStar(p);break;
    case PreviewStyle::TON618_CLUSTER:    paintTon618Cluster(p);  break;
    // 3C 273 system
    case PreviewStyle::C3273_JETBLOB:     paint3c273JetBlob(p);   break;
    case PreviewStyle::C3273_BLR:         paint3c273Blr(p);       break;
    case PreviewStyle::C3273_DWARF:       paint3c273Dwarf(p);     break;
    case PreviewStyle::C3273_CLOSESTAR:   paint3c273CloseStar(p); break;
    // J0529-4351 system
    case PreviewStyle::J0529_FASTBLOB:    paintJ0529FastBlob(p);  break;
    case PreviewStyle::J0529_UVCLUMP:     paintJ0529UvClump(p);   break;
    case PreviewStyle::J0529_TDESTAR:     paintJ0529TdeStar(p);   break;
    case PreviewStyle::J0529_GASSTREAM:   paintJ0529GasStream(p); break;
    case PreviewStyle::J0529_CLUSTER:     paintJ0529Cluster(p);   break;
    // M87 system
    case PreviewStyle::M87_JETKNOT:       paintM87JetKnot(p);     break;
    case PreviewStyle::M87_INNERSTAR:     paintM87InnerStar(p);   break;
    case PreviewStyle::M87_HOTGAS:        paintM87HotGas(p);      break;
    case PreviewStyle::M87_GLOBULAR:      paintM87Globular(p);    break;
    case PreviewStyle::M87_DWARF:         paintM87Dwarf(p);       break;
    default: break;
    }
}

// ---------------------------------------------------------------------------
// Shared utilities
// ---------------------------------------------------------------------------

void BlackHolePreviewWidget::drawBackground(QPainter &p, QColor nearColor, QColor farColor)
{
    QRectF rect(0, 0, width(), height());
    QRadialGradient g(rect.center(), rect.height() * 0.9f);
    g.setColorAt(0.0, nearColor);
    g.setColorAt(1.0, farColor);
    p.fillRect(rect, g);
}

void BlackHolePreviewWidget::drawStarfield(QPainter &p, bool dense)
{
    const auto  &stars = dense ? m_gcStars : m_stars;
    const float  W = float(width()), H = float(height());
    for (const auto &s : stars) {
        float tw = 0.75f + 0.25f * std::sin(float(m_frame) * 0.05f + s.x * 37.f);
        int   a  = int(s.brightness * tw * (dense ? 155.f : 200.f));
        QColor c(255, 255, 255, qMin(a, 255));
        p.setPen(c); p.setBrush(c);
        p.drawEllipse(QPointF(s.x * W, s.y * H), s.r, s.r);
    }
}

// Draw one half of the 3D accretion disc.
// Clip is set in screen coordinates BEFORE the local translate/rotate so the
// disc appears to physically wrap around the black-hole sphere.
//   front=false -> upper clip (far side of disc, passes behind BH, dimmer)
//   front=true  -> lower clip (near side of disc, passes in front of BH, brighter)
void BlackHolePreviewWidget::drawDiscHalf(QPainter &p,
                                           QPointF c, float bhR, float outerR, float scaleY,
                                           QColor inner, QColor mid, QColor outer,
                                           float alphaScale, bool front)
{
    p.save();
    p.setClipRect(front
                  ? QRect(0, int(c.y()), width(), height())
                  : QRect(0, 0,          width(), int(c.y())));
    p.translate(c);
    p.rotate(m_angle * 0.35f);

    const int   N    = 28;
    const float iR   = bhR * 1.1f;
    const float span = outerR - iR;

    for (int i = N; i >= 1; --i) {
        float t  = float(i) / float(N);
        float ri = iR + span * t;
        float b  = std::pow(1.f - t, 1.6f);

        QColor col;
        if (t < 0.30f) {
            float lt = t / 0.30f;
            col = QColor(
                int(float(inner.red())   * (1.f-lt) + float(mid.red())   * lt),
                int(float(inner.green()) * (1.f-lt) + float(mid.green()) * lt),
                int(float(inner.blue())  * (1.f-lt) + float(mid.blue())  * lt),
                int(float(inner.alpha()) * (1.f-lt) + float(mid.alpha()) * lt));
        } else {
            float lt = (t - 0.30f) / 0.70f;
            col = QColor(
                int(float(mid.red())   * (1.f-lt) + float(outer.red())   * lt),
                int(float(mid.green()) * (1.f-lt) + float(outer.green()) * lt),
                int(float(mid.blue())  * (1.f-lt) + float(outer.blue())  * lt),
                int(float(mid.alpha()) * (1.f-lt) + float(outer.alpha()) * lt));
        }

        float backDim = front ? 1.0f : 0.60f;
        col.setAlpha(qMax(0, qMin(255, int(float(col.alpha()) * b * alphaScale * backDim))));
        p.setPen(Qt::NoPen);
        p.setBrush(col);
        p.drawEllipse(QPointF(0.f, 0.f), ri, ri * scaleY);
    }
    p.restore();
}

void BlackHolePreviewWidget::drawBHSphere(QPainter &p, QPointF c, float r, QColor rim)
{
    QRadialGradient shadow(c, r * 2.5f);
    shadow.setColorAt(0.0f, QColor(0,0,0,200));
    shadow.setColorAt(0.5f, QColor(0,0,0,110));
    shadow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(shadow); p.setPen(Qt::NoPen);
    p.drawEllipse(c, r * 2.5f, r * 2.5f);

    QPointF lit(c.x() - r * 0.28f, c.y() - r * 0.28f);
    QRadialGradient sphere(lit, r * 1.5f);
    sphere.setColorAt(0.0f, QColor(14,11,24));
    sphere.setColorAt(0.6f, QColor(5, 4, 12));
    sphere.setColorAt(1.0f, Qt::black);
    p.setBrush(sphere);
    p.drawEllipse(c, r, r);

    QPointF rimPt(c.x(), c.y() + r * 0.62f);
    QRadialGradient rimGrad(rimPt, r * 1.05f);
    rimGrad.setColorAt(0.0f, QColor(rim.red(), rim.green(), rim.blue(), 85));
    rimGrad.setColorAt(1.0f, Qt::transparent);
    p.setBrush(rimGrad);
    p.drawEllipse(c, r, r);
}

void BlackHolePreviewWidget::drawPhotonRing(QPainter &p, QPointF c,
                                             float r, float scaleY, QColor col)
{
    float rp = r * 1.52f;
    for (int i = 3; i >= 1; --i) {
        float rr = rp + float(i) * r * 0.11f;
        QPen glow(QColor(col.red(), col.green(), col.blue(), 14*(4-i)), r * 0.10f * float(i));
        p.setPen(glow); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, rr, rr * scaleY);
    }
    QPen sharp(col, r * 0.07f);
    p.setPen(sharp); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, rp, rp * scaleY);
}

void BlackHolePreviewWidget::drawJet(QPainter &p, QPointF c, float bhR,
                                      float angleDeg, float length, float baseW,
                                      QColor col, float alpha)
{
    float ar = angleDeg * PI / 180.f;
    QPointF dir( std::sin(ar), -std::cos(ar));
    QPointF perp(-dir.y(), dir.x());
    QPointF tip (c.x() + dir.x() * length,     c.y() + dir.y() * length);
    QPointF base(c.x() + dir.x() * bhR * 0.5f, c.y() + dir.y() * bhR * 0.5f);
    QPointF bL = base + perp * baseW;
    QPointF bR = base - perp * baseW;

    float pulse = 0.68f + 0.32f * std::sin(m_phase + angleDeg * 0.05f);
    QLinearGradient grad(base, tip);
    grad.setColorAt(0.0f, QColor(col.red(), col.green(), col.blue(),
                                 qMin(255, int(alpha * pulse * 210.f))));
    grad.setColorAt(0.45f, QColor(col.red(), col.green(), col.blue(),
                                  qMin(255, int(alpha * pulse * 100.f))));
    grad.setColorAt(1.0f, Qt::transparent);

    QPainterPath cone;
    cone.moveTo(tip); cone.lineTo(bL); cone.lineTo(bR); cone.closeSubpath();
    p.setBrush(grad); p.setPen(Qt::NoPen);
    p.drawPath(cone);
}

void BlackHolePreviewWidget::drawCompanion(QPainter &p, QPointF c,
                                            float orbitR, float starR,
                                            float orbitAngleDeg, QColor col,
                                            bool addStream, QPointF bhPos)
{
    float ar = orbitAngleDeg * PI / 180.f;
    QPointF sp(c.x() + orbitR * std::cos(ar),
               c.y() + orbitR * std::sin(ar) * 0.35f);

    QRadialGradient glow(sp, starR * 3.5f);
    glow.setColorAt(0.0f, QColor(col.red(), col.green(), col.blue(), 95));
    glow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(glow); p.setPen(Qt::NoPen);
    p.drawEllipse(sp, starR * 3.5f, starR * 3.5f);

    QPointF litPt(sp.x() - starR * 0.25f, sp.y() - starR * 0.25f);
    QRadialGradient star(litPt, starR * 1.4f);
    star.setColorAt(0.0f, Qt::white);
    star.setColorAt(0.35f, col.lighter(115));
    star.setColorAt(1.0f, col);
    p.setBrush(star);
    p.drawEllipse(sp, starR, starR);

    if (addStream) {
        QPointF mid((sp.x() + bhPos.x()) * 0.5f,
                    (sp.y() + bhPos.y()) * 0.5f - 7.f);
        QPainterPath stream;
        stream.moveTo(sp); stream.quadTo(mid, bhPos);
        QPen wp(QColor(col.red(), col.green(), col.blue(), 55), 1.2f, Qt::DashLine);
        p.setPen(wp); p.setBrush(Qt::NoBrush);
        p.drawPath(stream);
    }
}

// ---------------------------------------------------------------------------
// Per-object paint methods
// ---------------------------------------------------------------------------

// Sgr A*: sub-Eddington RIAF, dim wispy disc, dense galactic star field,
// dusty haze, no jets.
void BlackHolePreviewWidget::paintSgrA(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.110f, oR=H*0.280f, sY=0.20f;

    drawBackground(p, QColor(18,11,8), QColor(5,3,6));

    QRadialGradient haze(c, H * 0.42f);
    haze.setColorAt(0.0f, QColor(90,42,12,28)); haze.setColorAt(1.0f, Qt::transparent);
    p.setBrush(haze); p.setPen(Qt::NoPen);
    p.drawEllipse(c, H * 0.42f, H * 0.32f);

    drawStarfield(p, true);

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(195,95,35,115), QColor(145,55,18,75), QColor(75,22,8,38), 0.70f, false);
    drawBHSphere(p, c, bR, QColor(195,95,35));
    drawPhotonRing(p, c, bR, sY, QColor(215,135,75,150));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(195,95,35,125), QColor(145,55,18,85), QColor(75,22,8,42), 0.80f, true);
}

// M87*: one-sided relativistic jet (Doppler boosted upper-left), faint
// counter-jet lower-right, bright orange disc, prominent photon ring.
void BlackHolePreviewWidget::paintM87(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.130f, oR=H*0.360f, sY=0.20f;

    drawBackground(p, QColor(12,8,18), QColor(3,3,10));
    drawStarfield(p);

    drawJet(p, c, bR, 140.f, H*0.34f, bR*0.95f, QColor(70,90,200),  0.30f);
    drawJet(p, c, bR, -42.f, H*0.46f, bR*1.10f, QColor(130,170,255),0.92f);

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,215,140,215), QColor(255,135,38,185), QColor(195,65,14,125), 1.00f, false);
    drawBHSphere(p, c, bR, QColor(255,165,55));
    drawPhotonRing(p, c, bR, sY, QColor(255,225,155,225));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,215,140,232), QColor(255,135,38,205), QColor(195,65,14,145), 1.18f, true);
}

// TON 618: hyperluminous, 6.6e10 Msun, blazing white-yellow disc almost
// overexposed at centre, powerful symmetric twin jets, faint host-galaxy glow.
void BlackHolePreviewWidget::paintTon618(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.148f, oR=H*0.440f, sY=0.18f;

    drawBackground(p, QColor(14,10,20), QColor(5,4,12));

    QRadialGradient galaxy(c, H * 0.50f);
    galaxy.setColorAt(0.0f, QColor(60,50,100,22)); galaxy.setColorAt(1.0f, Qt::transparent);
    p.setBrush(galaxy); p.setPen(Qt::NoPen);
    p.drawEllipse(c, H * 0.50f, H * 0.36f);

    drawStarfield(p);

    drawJet(p, c, bR,   0.f, H*0.47f, bR*1.30f, QColor(155,195,255), 0.95f);
    drawJet(p, c, bR, 180.f, H*0.47f, bR*1.30f, QColor(155,195,255), 0.85f);

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,255,215,225), QColor(255,205,75,205), QColor(255,115,18,155), 1.10f, false);
    drawBHSphere(p, c, bR, QColor(255,238,155));
    drawPhotonRing(p, c, bR, sY, QColor(255,248,195,238));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,255,215,240), QColor(255,205,75,225), QColor(255,115,18,175), 1.22f, true);

    QRadialGradient bloom(c, bR * 2.1f);
    bloom.setColorAt(0.0f, QColor(255,252,235,95));
    bloom.setColorAt(0.5f, QColor(255,240,195,40));
    bloom.setColorAt(1.0f, Qt::transparent);
    p.setBrush(bloom); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.1f, bR * 2.1f);
}

// 3C 273: first quasar identified (1963), yellow-white disc, one prominent
// angled optical jet upper-right, faint counter-jet lower-left.
void BlackHolePreviewWidget::paint3c273(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.118f, oR=H*0.342f, sY=0.20f;

    drawBackground(p, QColor(10,8,16), QColor(3,3,10));
    drawStarfield(p);

    drawJet(p, c, bR, 215.f, H*0.27f, bR*0.80f, QColor(110,125,205), 0.28f);
    drawJet(p, c, bR,  35.f, H*0.44f, bR*1.02f, QColor(210,215,255), 0.85f);

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,242,175,215), QColor(255,178,55,180), QColor(215,95,18,130), 1.00f, false);
    drawBHSphere(p, c, bR, QColor(255,195,75));
    drawPhotonRing(p, c, bR, sY, QColor(255,232,165,218));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(255,242,175,230), QColor(255,178,55,198), QColor(215,95,18,148), 1.12f, true);
}

// J0529-4351: most luminous quasar (2024), super-Eddington ~400 Msun/yr.
// Greenish-yellow disc from extreme UV ionisation; disc-wind arcs expand
// outward beyond the disc edge, distinguishing it visually from TON 618.
void BlackHolePreviewWidget::paintJ0529(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.138f, oR=H*0.400f, sY=0.19f;

    drawBackground(p, QColor(7,11,9), QColor(3,5,4));
    drawStarfield(p);

    // Animated disc-wind arcs
    {
        float bp = std::fmod(m_phase * 0.28f, 1.0f);
        for (int w = 0; w < 3; ++w) {
            float wp = std::fmod(bp + float(w) * 0.333f, 1.0f);
            float wr = oR * (1.08f + wp * 0.52f);
            int   a  = int(72.f * std::pow(1.f - wp, 1.4f));
            QPen wPen(QColor(175,255,155,a), 1.4f);
            p.setPen(wPen); p.setBrush(Qt::NoBrush);
            p.drawEllipse(c, wr, wr * sY * 1.5f);
        }
    }

    drawJet(p, c, bR,   0.f, H*0.44f, bR*1.22f, QColor(165,225,135), 0.90f);
    drawJet(p, c, bR, 180.f, H*0.44f, bR*1.22f, QColor(165,225,135), 0.80f);

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(215,255,155,225), QColor(255,218,45,205), QColor(255,125,18,155), 1.10f, false);
    drawBHSphere(p, c, bR, QColor(215,255,95));
    drawPhotonRing(p, c, bR, sY, QColor(215,255,155,222));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(215,255,155,240), QColor(255,218,45,225), QColor(255,125,18,175), 1.22f, true);
}

// Cygnus X-1: first confirmed BH (1972). Compact hot X-ray disc (blue-white
// inner region). Large blue OB supergiant companion with wind-accretion streams.
void BlackHolePreviewWidget::paintCygX1(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.092f, oR=H*0.265f, sY=0.22f;

    drawBackground(p, QColor(8,10,20), QColor(3,4,12));
    drawStarfield(p);

    float compAngle = m_angle * 0.85f;
    float orbitR    = float(width()) * 0.290f;
    float starR     = H * 0.072f;
    drawCompanion(p, c, orbitR, starR, compAngle, QColor(115,155,255), true, c);

    // Wind accretion streams
    {
        float ar = compAngle * PI / 180.f;
        QPointF sp(c.x() + orbitR * std::cos(ar),
                   c.y() + orbitR * std::sin(ar) * 0.35f);
        for (int i = 0; i < 3; ++i) {
            float off = float(i - 1) * 7.f;
            QPointF mid1((sp.x() * 0.58f + c.x() * 0.42f) + off,
                         (sp.y() + c.y()) * 0.50f - 10.f);
            QPainterPath stream;
            stream.moveTo(sp); stream.quadTo(mid1, c);
            QPen wp(QColor(100,128,255, 32 - i*7), 0.9f + float(i)*0.3f);
            p.setPen(wp); p.setBrush(Qt::NoBrush);
            p.drawPath(stream);
        }
    }

    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(195,218,255,225), QColor(255,155,55,185), QColor(195,65,18,128), 1.00f, false);
    drawBHSphere(p, c, bR, QColor(175,208,255));
    drawPhotonRing(p, c, bR, sY, QColor(195,218,255,215));
    drawDiscHalf(p, c, bR, oR, sY,
                 QColor(195,218,255,240), QColor(255,155,55,205), QColor(195,65,18,148), 1.16f, true);
}

// GW150914: isolated merger remnant, no disc, no companion.
// Animated gravitational-wave ripple rings expanding from the dark sphere.
void BlackHolePreviewWidget::paintGW(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.112f;

    drawBackground(p, QColor(5,5,11), QColor(2,2,6));
    drawStarfield(p);

    // GW ripple rings
    {
        float bp = std::fmod(m_phase * 0.22f, 1.0f);
        for (int w = 0; w < 5; ++w) {
            float wp = std::fmod(bp + float(w) * 0.20f, 1.0f);
            float wr = bR * 1.55f + H * 0.37f * wp;
            int   a  = int(85.f * std::pow(1.f - wp, 1.8f));
            QPen rp(QColor(215,205,175,a), 1.2f);
            p.setPen(rp); p.setBrush(Qt::NoBrush);
            p.drawEllipse(c, wr, wr * 0.82f);
        }
    }

    QRadialGradient shadow(c, bR * 2.2f);
    shadow.setColorAt(0.0f, QColor(0,0,0,185)); shadow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(shadow); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.2f, bR * 2.2f);

    QPointF lit(c.x() - bR*0.28f, c.y() - bR*0.28f);
    QRadialGradient sphere(lit, bR * 1.5f);
    sphere.setColorAt(0.0f, QColor(12,10,20));
    sphere.setColorAt(0.6f, QColor(4, 3, 10));
    sphere.setColorAt(1.0f, Qt::black);
    p.setBrush(sphere); p.drawEllipse(c, bR, bR);

    QPen lens(QColor(195,185,155,42), bR * 0.07f);
    p.setPen(lens); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, bR * 1.50f, bR * 1.50f);
}

// Gaia BH1: dormant, quiescent. No disc, no jets.
// G-type Sun-like companion on slow wide orbit. Einstein lensing ring.
void BlackHolePreviewWidget::paintGaiaBH1(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.086f;

    drawBackground(p, QColor(8,8,16), QColor(3,3,10));
    drawStarfield(p);

    float compAngle = m_angle * 0.32f;
    float orbitR    = float(width()) * 0.285f;
    float starR     = H * 0.038f;
    drawCompanion(p, c, orbitR, starR, compAngle, QColor(255,218,135));

    // Einstein lensing ring
    for (int i = 3; i >= 1; --i) {
        float lr = bR * (1.42f + float(i) * 0.14f);
        QPen lp(QColor(228,218,198, 14*(4-i)), bR * 0.10f * float(i));
        p.setPen(lp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, lr, lr);
    }
    QPen ering(QColor(228,218,198,52), bR * 0.07f);
    p.setPen(ering); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, bR * 1.45f, bR * 1.45f);

    QRadialGradient shadow(c, bR * 2.2f);
    shadow.setColorAt(0.0f, QColor(0,0,0,162)); shadow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(shadow); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.2f, bR * 2.2f);

    QPointF lit(c.x() - bR*0.26f, c.y() - bR*0.26f);
    QRadialGradient sphere(lit, bR * 1.4f);
    sphere.setColorAt(0.0f, QColor(11,9,19));
    sphere.setColorAt(1.0f, Qt::black);
    p.setBrush(sphere); p.drawEllipse(c, bR, bR);
}

// Gaia BH2: dormant, red-giant companion on wide ~3.5yr orbit. Companion fills ~30% of frame.
void BlackHolePreviewWidget::paintGaiaBH2(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), W=float(width()), bR=H*0.086f;

    drawBackground(p, QColor(10,6,14), QColor(4,2,8));
    drawStarfield(p);

    // Red giant companion, larger than in BH1 since it's evolved
    float compAngle = m_angle * 0.18f;   // slow orbit (1277-day period)
    float orbitR    = W * 0.300f;
    float starR     = H * 0.055f;        // notably larger than BH1's companion
    drawCompanion(p, c, orbitR, starR, compAngle, QColor(255, 120, 60));

    // Faint Roche-lobe hint, ellipse stretched toward BH
    float phase = qDegreesToRadians(compAngle);
    QPointF compPos(c.x() + orbitR * std::cos(phase),
                    c.y() + orbitR * std::sin(phase) * 0.55f);
    float rl = starR * 1.28f;
    QRadialGradient rlGrad(compPos, rl * 1.6f);
    rlGrad.setColorAt(0.0f, QColor(255,100,40,14));
    rlGrad.setColorAt(1.0f, Qt::transparent);
    p.setBrush(rlGrad); p.setPen(Qt::NoPen);
    p.drawEllipse(compPos, rl * 1.6f, rl * 1.6f);

    // Lensing ring
    for (int i = 3; i >= 1; --i) {
        float lr = bR * (1.42f + float(i) * 0.14f);
        p.setPen(QPen(QColor(210,185,170, 12*(4-i)), bR * 0.10f * float(i)));
        p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, lr, lr);
    }
    p.setPen(QPen(QColor(210,185,170,46), bR * 0.07f));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, bR * 1.45f, bR * 1.45f);

    QRadialGradient shadow(c, bR * 2.2f);
    shadow.setColorAt(0.0f, QColor(0,0,0,162)); shadow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(shadow); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.2f, bR * 2.2f);

    QRadialGradient sphere(QPointF(c.x()-bR*0.26f, c.y()-bR*0.26f), bR * 1.4f);
    sphere.setColorAt(0.0f, QColor(11,9,19)); sphere.setColorAt(1.0f, Qt::black);
    p.setBrush(sphere); p.drawEllipse(c, bR, bR);
}

// Gaia BH3: most massive stellar-mass BH in Milky Way. Metal-poor subgiant companion on
// ~11.6-year orbit. Show wide slow orbit + faint Population-II star + large lensing ring.
void BlackHolePreviewWidget::paintGaiaBH3(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), W=float(width()), bR=H*0.094f;  // larger BH shadow

    drawBackground(p, QColor(6,8,14), QColor(2,3,10));
    drawStarfield(p, true);   // denser starfield for nearby location (~590 ly)

    // Metal-poor subgiant, bluish-white, dimmer
    float compAngle = m_angle * 0.08f;   // ~11.6-year period → very slow
    float orbitR    = W * 0.310f;
    float starR     = H * 0.033f;
    drawCompanion(p, c, orbitR, starR, compAngle, QColor(195, 210, 255));

    // Extra faint pop-II background haze
    QRadialGradient popii(c, W * 0.46f);
    popii.setColorAt(0.0f, QColor(120,130,200, 8));
    popii.setColorAt(1.0f, Qt::transparent);
    p.setBrush(popii); p.setPen(Qt::NoPen);
    p.drawEllipse(c, W * 0.46f, W * 0.46f);

    // Larger lensing ring (more massive BH → larger shadow radius)
    for (int i = 4; i >= 1; --i) {
        float lr = bR * (1.38f + float(i) * 0.13f);
        p.setPen(QPen(QColor(200,200,230, 10*(5-i)), bR * 0.10f * float(i)));
        p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, lr, lr);
    }
    p.setPen(QPen(QColor(200,200,230,54), bR * 0.07f));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, bR * 1.45f, bR * 1.45f);

    QRadialGradient shadow(c, bR * 2.4f);
    shadow.setColorAt(0.0f, QColor(0,0,0,172)); shadow.setColorAt(1.0f, Qt::transparent);
    p.setBrush(shadow); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.4f, bR * 2.4f);

    QRadialGradient sphere(QPointF(c.x()-bR*0.28f, c.y()-bR*0.28f), bR * 1.5f);
    sphere.setColorAt(0.0f, QColor(10,9,20)); sphere.setColorAt(1.0f, Qt::black);
    p.setBrush(sphere); p.drawEllipse(c, bR, bR);
}

// V404 Cygni: actively accreting transient X-ray binary. K-giant donor, short 6.5-day orbit.
// Show bright accretion disc, jets, and rapid flickering via m_phase.
void BlackHolePreviewWidget::paintV404Cygni(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), bR=H*0.082f;

    drawBackground(p, QColor(12,4,6), QColor(5,1,2));
    drawStarfield(p);

    // Fast-orbiting K-giant companion
    float compAngle = m_angle * 1.65f;   // 6.5-day period → fast
    float orbitR    = float(width()) * 0.240f;
    float starR     = H * 0.042f;
    drawCompanion(p, c, orbitR, starR, compAngle, QColor(255, 140, 50));

    // Hot bright accretion disc
    float sY = 0.38f;
    drawDiscHalf(p, c, bR * 1.20f, bR * 3.2f, sY,
                 QColor(255,240,180), QColor(255,210,130), QColor(200,100,40), 1.0f, true);
    drawDiscHalf(p, c, bR * 1.20f, bR * 3.2f, sY,
                 QColor(255,240,180), QColor(255,210,130), QColor(200,100,40), 0.8f, false);

    // X-ray flickering, vary brightness with phase
    float flick = 0.5f + 0.5f * std::sin(m_phase * 5.3f);
    int fa = int(40 + 60 * flick);
    QRadialGradient xray(c, bR * 2.8f);
    xray.setColorAt(0.0f, QColor(255,240,200,fa));
    xray.setColorAt(1.0f, Qt::transparent);
    p.setBrush(xray); p.setPen(Qt::NoPen);
    p.drawEllipse(c, bR * 2.8f, bR * 2.8f);

    // Relativistic jets (ejection blobs)
    drawJet(p, c, bR, m_angle, H*0.36f, bR*1.0f, QColor(200,160,255), 0.80f);

    drawBHSphere(p, c, bR, QColor(255,180,80));
    drawPhotonRing(p, c, bR, sY, QColor(255,200,100,200));
}

// Phoenix A: ultramassive ~100-billion-solar-mass BH at centre of galaxy cluster.
// Enormously large shadow, dual radio lobes, surrounding galaxy haze.
void BlackHolePreviewWidget::paintPhoenixA(QPainter &p)
{
    const QPointF c(width() * 0.50f, height() * 0.52f);
    const float H=float(height()), W=float(width()), bR=H*0.130f;  // huge shadow

    // Deep-field galaxy-cluster background
    drawBackground(p, QColor(2,2,8), QColor(0,0,4));
    // Scatter a few faint background galaxies
    {
        struct GalPos { float x, y, r; QColor col; };
        static const GalPos gals[] = {
            {0.14f, 0.18f, 0.018f, QColor(210,195,165,40)},
            {0.82f, 0.22f, 0.014f, QColor(190,175,150,32)},
            {0.71f, 0.75f, 0.020f, QColor(200,190,160,36)},
            {0.25f, 0.80f, 0.012f, QColor(215,200,170,28)},
            {0.88f, 0.55f, 0.016f, QColor(205,185,155,34)},
        };
        for (auto &g : gals) {
            QRadialGradient gl(QPointF(g.x*W, g.y*H), g.r*W*2.2f);
            gl.setColorAt(0.0f, g.col); gl.setColorAt(1.0f, Qt::transparent);
            p.setBrush(gl); p.setPen(Qt::NoPen);
            p.drawEllipse(QPointF(g.x*W, g.y*H), g.r*W*2.2f, g.r*W*2.2f);
        }
    }

    // Host galaxy haze, massive elliptical BCG
    QRadialGradient bcg(c, W * 0.42f);
    bcg.setColorAt(0.0f, QColor(220,200,160,55));
    bcg.setColorAt(0.5f, QColor(180,160,120,20));
    bcg.setColorAt(1.0f, Qt::transparent);
    p.setBrush(bcg); p.setPen(Qt::NoPen);
    p.drawEllipse(c, W * 0.42f, W * 0.30f);

    // Jet feedback radio lobes
    float lobeDist = bR * 2.6f;
    float lobeSz   = bR * 1.8f;
    float lobePulse = 0.7f + 0.3f * std::sin(m_phase * 0.9f);
    for (int side : {-1, 1}) {
        QPointF lobeC(c.x(), c.y() + side * lobeDist);
        QRadialGradient lobe(lobeC, lobeSz);
        lobe.setColorAt(0.0f, QColor(80,100,220, int(55*lobePulse)));
        lobe.setColorAt(1.0f, Qt::transparent);
        p.setBrush(lobe); p.setPen(Qt::NoPen);
        p.drawEllipse(lobeC, lobeSz, lobeSz * 0.65f);
    }

    // Jet spine
    drawJet(p, c, bR, m_angle * 0.05f, H*0.42f, bR*1.3f, QColor(100,120,255), 0.85f);
    drawJet(p, c, bR, m_angle * 0.05f + 180.f, H*0.42f, bR*1.3f, QColor(100,120,255), 0.75f);

    // Photon ring (very large)
    float discSY = 0.30f;
    drawPhotonRing(p, c, bR, discSY, QColor(255,180,100,220));

    // Massive disc, wide and thick
    drawDiscHalf(p, c, bR * 1.15f, bR * 4.5f, discSY,
                 QColor(255,230,160), QColor(255,200,120), QColor(180,80,30), 1.0f, true);
    drawDiscHalf(p, c, bR * 1.15f, bR * 4.5f, discSY,
                 QColor(255,230,160), QColor(255,200,120), QColor(180,80,30), 0.8f, false);

    drawBHSphere(p, c, bR, QColor(255,160,60));
}

// ---------------------------------------------------------------------------
// Shared helpers added for system-object and research-scenario scenes
// ---------------------------------------------------------------------------

// Draw a thin orbit ellipse (optionally dashed)
void BlackHolePreviewWidget::drawOrbitPath(QPainter &p, QPointF c,
                                            float a, float b, float angleDeg,
                                            QColor col, float alpha, bool dashed)
{
    QPen pen(QColor(col.red(), col.green(), col.blue(), int(alpha * 255.f)), 1.2f);
    if (dashed) {
        pen.setStyle(Qt::DashLine);
        pen.setDashPattern({4.0, 4.0});
    }
    p.save();
    p.setPen(pen); p.setBrush(Qt::NoBrush);
    p.translate(c);
    p.rotate(angleDeg);
    p.drawEllipse(QPointF(0,0), a, b);
    p.restore();
}

// Draw a small "context" black hole + optional thin disc for system-object scenes
void BlackHolePreviewWidget::drawContextBH(QPainter &p, QPointF c, float bhR,
                                            bool drawDisc, float discOuter,
                                            float scaleY, QColor discColor)
{
    if (drawDisc) {
        // back half
        p.save();
        p.setClipRect(0, 0, width(), int(c.y()));
        p.translate(c); p.scale(1.f, scaleY);
        for (int i = 8; i >= 1; --i) {
            float r = bhR + (discOuter - bhR) * float(i) / 8.f;
            int   a = 35 + i * 12;
            QPen rp(QColor(discColor.red(), discColor.green(), discColor.blue(), a), 2.5f);
            p.setPen(rp); p.setBrush(Qt::NoBrush);
            p.drawEllipse(QPointF(0,0), r, r);
        }
        p.restore();
        // front half
        p.save();
        p.setClipRect(0, int(c.y()), width(), height());
        p.translate(c); p.scale(1.f, scaleY);
        for (int i = 8; i >= 1; --i) {
            float r = bhR + (discOuter - bhR) * float(i) / 8.f;
            int   a = 45 + i * 14;
            QPen rp(QColor(discColor.red(), discColor.green(), discColor.blue(), a), 2.5f);
            p.setPen(rp); p.setBrush(Qt::NoBrush);
            p.drawEllipse(QPointF(0,0), r, r);
        }
        p.restore();
    }
    // sphere
    QPointF lit(c.x() - bhR * 0.3f, c.y() - bhR * 0.3f);
    QRadialGradient sg(lit, bhR * 1.5f);
    sg.setColorAt(0.0f, QColor(18,14,28)); sg.setColorAt(1.0f, Qt::black);
    p.setPen(Qt::NoPen); p.setBrush(sg);
    p.drawEllipse(c, bhR, bhR);
}

// ---------------------------------------------------------------------------
// Research Scenario paint methods
// ---------------------------------------------------------------------------

// ISCO Unstable (5M): test body slowly spiralling inward to plunge
void BlackHolePreviewWidget::paintIscoUnstable(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.10f;

    drawBackground(p, QColor(8,4,18), QColor(2,2,8));
    drawStarfield(p);
    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    // Photon ring hint
    QPen pr(QColor(255,200,80,55), 1.5f);
    p.setPen(pr); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, bR * 1.5f, bR * 1.5f);

    // ISCO radius ring (r=5M)
    QPen ir(QColor(255,80,60,80), 1.0f, Qt::DashLine);
    ir.setDashPattern({4.0, 4.0});
    p.setPen(ir); p.setBrush(Qt::NoBrush);
    float iscoR = bR * 2.5f;
    p.drawEllipse(c, iscoR, iscoR * 0.38f);

    // Spiralling path, tight inward arc
    QPainterPath spiral;
    bool started = false;
    for (int i = 0; i <= 180; ++i) {
        float t    = float(i) / 180.f;
        float rad  = iscoR * (1.f - t * 0.55f);
        float ang  = (m_angle + float(i) * 2.4f) * PI / 180.f;
        float sx   = c.x() + rad * std::cos(ang);
        float sy   = c.y() + rad * 0.38f * std::sin(ang);
        if (!started) { spiral.moveTo(sx, sy); started = true; }
        else          { spiral.lineTo(sx, sy); }
    }
    QPen sp(QColor(255,120,60,200), 1.5f);
    p.setPen(sp); p.setBrush(Qt::NoBrush);
    p.drawPath(spiral);

    // Body dot at leading end
    float ang0 = m_angle * PI / 180.f;
    float r0   = iscoR * 0.475f;
    QPointF body(c.x() + r0 * std::cos(ang0), c.y() + r0 * 0.38f * std::sin(ang0));
    QRadialGradient bg(body, 5.f); bg.setColorAt(0.f,QColor(255,160,80)); bg.setColorAt(1.f,Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(body, 5.f, 5.f);
}

// ISCO Critical (6M): marginally stable circle, slowly wobbling
void BlackHolePreviewWidget::paintIscoCritical(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.10f;

    drawBackground(p, QColor(5,10,20), QColor(2,2,10));
    drawStarfield(p);
    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    float iscoR = bR * 3.0f;
    // Glow on ISCO ring
    for (int i = 3; i >= 1; --i) {
        QPen gp(QColor(80,200,255, 20 * i), float(i) * 2.0f);
        p.setPen(gp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, iscoR, iscoR * 0.40f);
    }
    QPen ip(QColor(80,200,255,180), 1.5f);
    p.setPen(ip); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, iscoR, iscoR * 0.40f);

    // Body at ISCO, slightly perturbed radially
    float wobble = iscoR * (1.f + 0.04f * std::sin(m_phase * 2.0f));
    float ang = m_angle * PI / 180.f;
    QPointF body(c.x() + wobble * std::cos(ang), c.y() + wobble * 0.40f * std::sin(ang));
    QRadialGradient bg(body, 6.f); bg.setColorAt(0.f,QColor(120,220,255)); bg.setColorAt(1.f,Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(body, 6.f, 6.f);

    // Label annotation line
    QPen lp(QColor(80,200,255,70), 1.f, Qt::DotLine);
    p.setPen(lp);
    p.drawLine(c, QPointF(c.x() + iscoR * 0.707f, c.y() - iscoR * 0.40f * 0.707f));
}

// ISCO Stable (7M): healthy circular orbit well outside ISCO
void BlackHolePreviewWidget::paintIscoStable(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.10f;

    drawBackground(p, QColor(4,12,8), QColor(2,4,2));
    drawStarfield(p);
    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    float orbitR = bR * 3.5f;
    // Green stable orbit ring
    for (int i = 3; i >= 1; --i) {
        QPen gp(QColor(60,220,100, 18 * i), float(i) * 2.0f);
        p.setPen(gp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, orbitR, orbitR * 0.40f);
    }
    QPen ip(QColor(60,220,100,200), 1.5f);
    p.setPen(ip); p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, orbitR, orbitR * 0.40f);

    float ang = m_angle * PI / 180.f;
    QPointF body(c.x() + orbitR * std::cos(ang), c.y() + orbitR * 0.40f * std::sin(ang));
    QRadialGradient bg(body, 6.f); bg.setColorAt(0.f,QColor(100,255,140)); bg.setColorAt(1.f,Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(body, 6.f, 6.f);
}

// Photon Sphere: photons arcing around at r=3M, some captured, some escaping
void BlackHolePreviewWidget::paintPhotonSphere(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.10f;

    drawBackground(p, QColor(8,6,18), QColor(2,2,8));
    drawStarfield(p);
    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    float psR = bR * 1.5f; // photon sphere r=3M
    // Photon sphere glow
    for (int i = 4; i >= 1; --i) {
        QPen gp(QColor(255,220,80, 12 * i), float(i) * 1.8f);
        p.setPen(gp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, psR, psR);
    }

    // Six photon arcs orbiting at different phases
    const float offsets[6] = {0.f, 60.f, 120.f, 180.f, 240.f, 300.f};
    for (int ph = 0; ph < 6; ++ph) {
        bool captured = (ph < 3); // inner three captured, outer three escape
        float baseAng = (m_angle + offsets[ph]) * PI / 180.f;
        // Arc: partial circle around photon sphere
        QPainterPath arc;
        bool started = false;
        int steps = captured ? 200 : 120;
        for (int i = 0; i <= steps; ++i) {
            float t   = float(i) / float(steps);
            float ang = baseAng + t * (captured ? 2.8f : 1.4f);
            float r   = psR * (captured ? (1.f - t * 0.35f) : (1.f + t * 0.55f));
            float sx  = c.x() + r * std::cos(ang);
            float sy  = c.y() + r * std::sin(ang);
            if (!started) { arc.moveTo(sx, sy); started = true; }
            else          { arc.lineTo(sx, sy); }
        }
        QColor col = captured ? QColor(255,100,80,160) : QColor(180,255,180,140);
        p.setPen(QPen(col, 1.2f)); p.setBrush(Qt::NoBrush);
        p.drawPath(arc);
    }
}

// Radial Infall: test particle dropped from rest, coordinate velocity → 0
void BlackHolePreviewWidget::paintRadialInfall(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.10f;

    drawBackground(p, QColor(6,4,14), QColor(2,2,8));
    drawStarfield(p);
    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    // Radial trajectory line from top of widget to BH
    float dropFromY = H * 0.08f;
    QPen lp(QColor(200,200,255,80), 1.2f, Qt::DashLine);
    lp.setDashPattern({5.0, 3.0});
    p.setPen(lp);
    p.drawLine(QPointF(c.x(), dropFromY), QPointF(c.x(), c.y() - bR));

    // Coordinate distance rings (gravitational time dilation)
    float radii[4] = {bR*4.f, bR*3.0f, bR*2.0f, bR*1.3f};
    int   alphas[4]= {40, 55, 75, 100};
    for (int i = 0; i < 4; ++i) {
        QPen rp(QColor(140,140,255, alphas[i]), 1.f, Qt::DotLine);
        p.setPen(rp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, radii[i], radii[i]);
    }

    // Infalling body: slows down asymptotically near horizon
    // Position: y decreases quickly then asymptotes near BH surface
    float phase = std::fmod(m_phase * 0.55f, 1.0f);
    float factor = 1.f - std::exp(-phase * 5.f); // 0→1 asymptotic
    float yBody = dropFromY + (c.y() - bR * 1.1f - dropFromY) * factor;
    QPointF body(c.x(), yBody);
    QRadialGradient bg(body, 5.f);
    bg.setColorAt(0.f, QColor(180,180,255)); bg.setColorAt(1.f, Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(body, 5.f, 5.f);

    // Redshift colour: shifts from white→orange→red as body approaches
    float redshift = factor * factor;
    QColor bodyCol(int(220 - redshift * 100), int(220 - redshift * 160), int(255 - redshift * 220));
    QRadialGradient bc(body, 4.f);
    bc.setColorAt(0.f, bodyCol); bc.setColorAt(1.f, Qt::transparent);
    p.setBrush(bc);
    p.drawEllipse(body, 4.f, 4.f);
}

// Tidal Disruption: star on eccentric orbit stretched into debris stream
void BlackHolePreviewWidget::paintTidalDisrupt(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.45f, H * 0.52f);
    const float bR = H * 0.09f;

    drawBackground(p, QColor(12,4,8), QColor(4,2,4));
    drawStarfield(p);

    // Debris stream: elongated arc through pericentre
    QPainterPath stream;
    bool started = false;
    float aStar = bR * 4.0f;
    for (int i = 0; i <= 160; ++i) {
        float t  = float(i) / 160.f;
        float ang = (m_angle * 0.3f + t * 260.f) * PI / 180.f;
        float r  = aStar * (1.2f - 0.6f * std::cos(ang));
        float sx = c.x() + r * std::cos(ang);
        float sy = c.y() + r * 0.50f * std::sin(ang);
        if (!started) { stream.moveTo(sx, sy); started = true; }
        else          { stream.lineTo(sx, sy); }
    }
    // Outer half, fading escape arm
    QLinearGradient streamGrad(c.x(), c.y() - aStar, c.x() + aStar * 1.4f, c.y());
    streamGrad.setColorAt(0.0, QColor(255,140,80,200));
    streamGrad.setColorAt(0.6, QColor(255,80,40,100));
    streamGrad.setColorAt(1.0, Qt::transparent);
    p.setPen(QPen(QBrush(streamGrad), 2.5f)); p.setBrush(Qt::NoBrush);
    p.drawPath(stream);

    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    // Stretched star at pericentre, elongated ellipse
    float periAng = (m_angle * 0.3f + 90.f) * PI / 180.f;
    float periR   = aStar * 0.4f;
    QPointF peri(c.x() + periR * std::cos(periAng), c.y() + periR * 0.50f * std::sin(periAng));
    p.save(); p.translate(peri); p.rotate(m_angle * 0.3f);
    QRadialGradient sg(QPointF(0,0), 12.f);
    sg.setColorAt(0.f, QColor(255,200,120,220)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(QPointF(0,0), 12.f, 5.f);
    p.restore();
}

// Pulsar Inspiral: NS on decaying orbit, jets, GW chirp rings
void BlackHolePreviewWidget::paintPulsarInspiral(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.09f;

    drawBackground(p, QColor(4,6,16), QColor(2,2,8));
    drawStarfield(p);

    // GW ripple rings expanding outward
    for (int i = 0; i < 4; ++i) {
        float phase = std::fmod(m_phase + float(i) * 0.62f, 2.5f);
        float r     = bR * 1.6f + phase * bR * 3.5f;
        int   a     = int(90.f * (1.f - phase / 2.5f));
        QPen gp(QColor(100,180,255, a), 1.4f - phase * 0.4f);
        p.setPen(gp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, r, r);
    }

    drawContextBH(p, c, bR, false, 0.f, 0.f, Qt::white);

    // Decaying orbit path (shrinking spiral)
    float decay = std::fmod(m_phase * 0.2f, 1.f);
    float orbitA = bR * (3.5f - decay * 1.8f);
    drawOrbitPath(p, c, orbitA, orbitA * 0.42f, 0.f, QColor(200,220,255), 0.45f, true);

    // NS body
    float nsAng = m_angle * PI / 180.f;
    QPointF ns(c.x() + orbitA * std::cos(nsAng), c.y() + orbitA * 0.42f * std::sin(nsAng));

    // Pulsar jets (rotating with the NS)
    drawJet(p, ns, 4.f, m_angle * 2.f,        H * 0.18f, 2.5f, QColor(120,200,255), 0.7f);
    drawJet(p, ns, 4.f, m_angle * 2.f + 180.f, H * 0.18f, 2.5f, QColor(120,200,255), 0.7f);

    QRadialGradient ng(ns, 5.f);
    ng.setColorAt(0.f, QColor(180,220,255)); ng.setColorAt(1.f, Qt::transparent);
    p.setBrush(ng); p.setPen(Qt::NoPen);
    p.drawEllipse(ns, 5.f, 5.f);
}

// ---------------------------------------------------------------------------
// Sgr A* system paint methods
// ---------------------------------------------------------------------------

// Helper: draws galactic-centre background (dense star field + dusty haze)
static void sgraBg(BlackHolePreviewWidget *self, QPainter &p, QPointF c, float bR)
{
    // Dusty warm haze, draw as semi-transparent gradient wash
    QRadialGradient haze(c, bR * 7.f);
    haze.setColorAt(0.0f, QColor(60,40,20,55));
    haze.setColorAt(0.5f, QColor(30,20,8,30));
    haze.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0, 0, self->width(), self->height()), haze);
}

// S2 Analog: eccentric precessing orbit around Sgr A*
void BlackHolePreviewWidget::paintSgraS2(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.48f, H * 0.52f);
    const float bR = H * 0.078f;

    drawBackground(p, QColor(14,10,20), QColor(4,4,12));
    drawStarfield(p, true);
    sgraBg(this, p, c, bR);

    // Precessing orbit, semi-major axis oriented by m_angle * 0.06 (slow)
    float orbitA = bR * 3.8f, orbitB = orbitA * 0.30f;
    float precess = m_angle * 0.06f; // slow precession
    drawOrbitPath(p, c, orbitA, orbitB, precess, QColor(200,200,255), 0.55f, true);

    drawContextBH(p, c, bR, true, bR * 2.2f, 0.22f, QColor(100,80,60));

    // S2 star position
    float starAngle = (m_angle + precess) * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(starAngle), c.y() + orbitB * std::sin(starAngle));
    // Blue-white massive star
    QRadialGradient sg(star, 6.f);
    sg.setColorAt(0.f, QColor(200,220,255,240)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 6.f, 6.f);
    p.setBrush(QColor(220,235,255)); p.drawEllipse(star, 3.f, 3.f);
}

// S14-like close orbit: ultra-tight, stronger precession
void BlackHolePreviewWidget::paintSgraS14(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.085f;

    drawBackground(p, QColor(12,8,18), QColor(3,3,10));
    drawStarfield(p, true);
    sgraBg(this, p, c, bR);

    // Very tight orbit, small semi-major axis, faster precession
    float orbitA = bR * 2.6f, orbitB = orbitA * 0.20f;
    float precess = m_angle * 0.18f;
    // Draw several ghost orbit ellipses to show precession rosette
    for (int i = 3; i >= 1; --i) {
        drawOrbitPath(p, c, orbitA, orbitB, precess - float(i) * 22.f,
                      QColor(180,160,255), 0.18f * float(i), true);
    }
    drawOrbitPath(p, c, orbitA, orbitB, precess, QColor(220,200,255), 0.70f, false);

    drawContextBH(p, c, bR, true, bR * 1.9f, 0.22f, QColor(100,80,60));

    float ang = (m_angle + precess) * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    QRadialGradient sg(star, 5.f);
    sg.setColorAt(0.f, QColor(255,240,200,240)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 5.f, 5.f);
    p.setBrush(QColor(255,245,210)); p.drawEllipse(star, 2.5f, 2.5f);
}

// IRS 16 cluster star: wide orbit, dusty, OB star
void BlackHolePreviewWidget::paintSgraIrs16(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.070f;

    drawBackground(p, QColor(16,11,22), QColor(5,4,12));
    drawStarfield(p, true);
    sgraBg(this, p, c, bR);

    // Wide nearly circular orbit
    float orbitA = bR * 4.8f, orbitB = orbitA * 0.50f;
    drawOrbitPath(p, c, orbitA, orbitB, 15.f, QColor(180,200,230), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.22f, QColor(90,70,55));

    // OB star, bright blue-white, large
    float ang = m_angle * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    for (int i = 3; i >= 1; --i) {
        QRadialGradient gl(star, float(i) * 4.f);
        gl.setColorAt(0.f, QColor(160,200,255, 50 * i)); gl.setColorAt(1.f, Qt::transparent);
        p.setBrush(gl); p.setPen(Qt::NoPen);
        p.drawEllipse(star, float(i) * 4.f, float(i) * 4.f);
    }
    p.setBrush(QColor(210,230,255)); p.drawEllipse(star, 4.5f, 4.5f);
}

// Circumnuclear Gas Clump: tidal stretching as it orbits
void BlackHolePreviewWidget::paintSgraGasClump(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.072f;

    drawBackground(p, QColor(10,8,18), QColor(3,3,10));
    drawStarfield(p, true);
    sgraBg(this, p, c, bR);

    float orbitA = bR * 4.5f, orbitB = orbitA * 0.40f;
    drawOrbitPath(p, c, orbitA, orbitB, -12.f, QColor(120,160,200), 0.40f, true);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.22f, QColor(90,70,55));

    // Gas clump, elongated, colour based on phase (near pericentre = hotter)
    float ang = m_angle * PI / 180.f;
    float dist = std::sqrt(std::pow(orbitA * std::cos(ang),2) + std::pow(orbitB * std::sin(ang),2));
    float heating = 1.f - std::min(dist / (orbitA * 1.1f), 1.f);
    QColor clumpCol(int(80 + heating * 120), int(140 + heating * 60), int(200 - heating * 80), 200);

    QPointF clump(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    // Tidal stretching: elongated along orbit direction
    p.save(); p.translate(clump); p.rotate(m_angle);
    QRadialGradient cg(QPointF(0,0), 14.f);
    cg.setColorAt(0.f, clumpCol); cg.setColorAt(1.f, Qt::transparent);
    p.setBrush(cg); p.setPen(Qt::NoPen);
    p.drawEllipse(QPointF(0,0), 14.f, 7.f);
    p.restore();
}

// Generic S-cluster member
void BlackHolePreviewWidget::paintSgraSMember(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.078f;

    drawBackground(p, QColor(12,9,18), QColor(4,3,11));
    drawStarfield(p, true);
    sgraBg(this, p, c, bR);

    float orbitA = bR * 3.5f, orbitB = orbitA * 0.36f;
    drawOrbitPath(p, c, orbitA, orbitB, m_angle * 0.04f, QColor(200,180,230), 0.50f, true);
    drawContextBH(p, c, bR, true, bR * 2.0f, 0.22f, QColor(90,70,55));

    float ang = m_angle * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    QRadialGradient sg(star, 5.f);
    sg.setColorAt(0.f, QColor(255,255,200,230)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 5.f, 5.f);
    p.setBrush(QColor(255,255,210)); p.drawEllipse(star, 2.5f, 2.5f);
}

// ---------------------------------------------------------------------------
// TON 618 system paint methods
// ---------------------------------------------------------------------------

// Helper: TON 618 context, blazing purple-white quasar background
static void ton618Bg(QPainter &p, int w, int h)
{
    QRadialGradient qg(QPointF(w*0.5f, h*0.5f), h * 1.1f);
    qg.setColorAt(0.0f, QColor(30,18,40,80));
    qg.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,w,h), qg);
}

// Inner BLR gas clump: blazing quasar core, orbiting gas blob
void BlackHolePreviewWidget::paintTon618Blr(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(18,8,28), QColor(5,2,10));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    drawContextBH(p, c, bR, true, bR * 3.0f, 0.20f, QColor(200,160,255));

    // Photon-ionised BLR glow
    QRadialGradient blrGlow(c, bR * 3.5f);
    blrGlow.setColorAt(0.0f, QColor(160,80,255,50));
    blrGlow.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,W,H), blrGlow);

    // Orbiting gas blob
    float orbitR = bR * 2.4f;
    drawOrbitPath(p, c, orbitR, orbitR * 0.35f, 0.f, QColor(180,120,255), 0.35f, true);
    float ang = m_angle * PI / 180.f;
    QPointF blob(c.x() + orbitR * std::cos(ang), c.y() + orbitR * 0.35f * std::sin(ang));
    QRadialGradient bg(blob, 9.f);
    bg.setColorAt(0.f, QColor(220,160,255,220)); bg.setColorAt(1.f, Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(blob, 9.f, 9.f);
}

// Tidally stripped star: debris stream, stripped stellar remnant
void BlackHolePreviewWidget::paintTon618TdStar(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.44f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(16,7,24), QColor(4,2,9));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    // Debris stream arc
    QPainterPath stream;
    float a = bR * 4.5f;
    bool started = false;
    for (int i = 0; i <= 140; ++i) {
        float t   = float(i) / 140.f;
        float ang = (m_angle * 0.25f + t * 220.f) * PI / 180.f;
        float r   = a * (1.1f - 0.5f * std::cos(ang));
        float sx  = c.x() + r * std::cos(ang);
        float sy  = c.y() + r * 0.45f * std::sin(ang);
        if (!started) { stream.moveTo(sx,sy); started = true; }
        else          { stream.lineTo(sx,sy); }
    }
    QPen sp(QColor(200,140,255,160), 2.0f); p.setPen(sp); p.setBrush(Qt::NoBrush);
    p.drawPath(stream);

    drawContextBH(p, c, bR, true, bR * 2.5f, 0.20f, QColor(200,160,255));

    // Stripped remnant body
    float ang = m_angle * 0.28f * PI / 180.f;
    float r   = a * 0.52f;
    QPointF rem(c.x() + r * std::cos(ang), c.y() + r * 0.45f * std::sin(ang));
    QRadialGradient rg(rem, 7.f);
    rg.setColorAt(0.f, QColor(255,200,180,210)); rg.setColorAt(1.f, Qt::transparent);
    p.setBrush(rg); p.setPen(Qt::NoPen);
    p.drawEllipse(rem, 7.f, 7.f);
}

// Hot gas filament: TDE debris sheared by differential Keplerian rotation
void BlackHolePreviewWidget::paintTon618GasFilm(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(14,6,22), QColor(4,2,9));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    // Sheared filament: a long curved arc rotating with m_angle
    for (int layer = 0; layer < 3; ++layer) {
        QPainterPath film;
        bool started = false;
        float wid = float(layer) * 0.8f;
        for (int i = 0; i <= 200; ++i) {
            float t   = float(i) / 200.f;
            float baseR = bR * (2.0f + t * 3.5f);
            float ang  = (m_angle * 0.2f + t * 180.f + float(layer) * 8.f) * PI / 180.f;
            float sx   = c.x() + baseR * std::cos(ang);
            float sy   = c.y() + baseR * 0.32f * std::sin(ang);
            if (!started) { film.moveTo(sx,sy); started = true; }
            else          { film.lineTo(sx,sy); }
        }
        int a = 120 - layer * 35;
        p.setPen(QPen(QColor(180,120,255,a), 2.2f - float(layer)*0.6f));
        p.setBrush(Qt::NoBrush); p.drawPath(film);
    }

    drawContextBH(p, c, bR, true, bR * 2.2f, 0.20f, QColor(200,160,255));
}

// Plunging star: near-radial trajectory, direct capture
void BlackHolePreviewWidget::paintTon618Plunge(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(14,6,22), QColor(4,2,9));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    drawContextBH(p, c, bR, true, bR * 2.5f, 0.20f, QColor(200,160,255));

    // Near-radial plunge path, almost a straight line from upper-right
    float startX = c.x() + W * 0.35f, startY = c.y() - H * 0.34f;
    float phase  = std::fmod(m_phase * 0.4f, 1.0f);
    float bodyX  = startX + (c.x() - bR * 1.05f - startX) * phase;
    float bodyY  = startY + (c.y() - bR * 1.05f - startY) * phase;

    // Trajectory trail
    QPen tp(QColor(255,180,140,120), 1.5f, Qt::DashLine);
    tp.setDashPattern({4.0,3.0});
    p.setPen(tp); p.setBrush(Qt::NoBrush);
    p.drawLine(QPointF(startX, startY), QPointF(c.x(), c.y() - bR * 1.05f));

    // Body (no flare, swallowed whole)
    QPointF body(bodyX, bodyY);
    float alpha = 1.f - phase * 0.85f;
    QRadialGradient bg(body, 6.f);
    bg.setColorAt(0.f, QColor(255,220,180, int(alpha*220))); bg.setColorAt(1.f, Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(body, 6.f, 6.f);
}

// Outer stellar orbit: wide, apsidal precession
void BlackHolePreviewWidget::paintTon618OuterStar(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.058f;

    drawBackground(p, QColor(14,6,22), QColor(4,2,9));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    float orbitA = bR * 5.2f, orbitB = orbitA * 0.60f;
    float precess = m_angle * 0.03f;
    for (int i = 2; i >= 0; --i)
        drawOrbitPath(p, c, orbitA, orbitB, precess - float(i)*18.f,
                      QColor(180,140,255), 0.18f*(3-i), true);
    drawOrbitPath(p, c, orbitA, orbitB, precess, QColor(200,160,255), 0.60f, false);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.20f, QColor(200,160,255));

    float ang = (m_angle + precess) * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    QRadialGradient sg(star, 5.f);
    sg.setColorAt(0.f, QColor(255,250,210,230)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 5.f, 5.f);
    p.setBrush(QColor(255,252,220)); p.drawEllipse(star, 2.5f, 2.5f);
}

// Infalling stellar cluster: multi-point cluster, dynamical friction arc
void BlackHolePreviewWidget::paintTon618Cluster(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.058f;

    drawBackground(p, QColor(14,6,22), QColor(4,2,9));
    drawStarfield(p);
    ton618Bg(p, width(), height());

    float orbitA = bR * 5.8f, orbitB = orbitA * 0.55f;
    // Decaying orbit, slightly smaller each pass
    float decay = std::fmod(m_phase * 0.15f, 1.f);
    float curA  = orbitA * (1.f - decay * 0.18f);
    float curB  = orbitB * (1.f - decay * 0.18f);
    drawOrbitPath(p, c, curA, curB, 0.f, QColor(200,180,255), 0.50f, true);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.20f, QColor(200,160,255));

    // Cluster: scatter of stars around centroid
    float cAng = m_angle * PI / 180.f;
    QPointF clCentre(c.x() + curA * std::cos(cAng), c.y() + curB * std::sin(cAng));
    std::mt19937 rng(77u);
    std::uniform_real_distribution<float> offs(-11.f, 11.f);
    for (int i = 0; i < 12; ++i) {
        float ox = offs(rng), oy = offs(rng);
        QPointF sp(clCentre.x() + ox, clCentre.y() + oy);
        float br = 0.5f + float(rng() & 15) / 30.f;
        p.setBrush(QColor(220,200,255,int(br*200))); p.setPen(Qt::NoPen);
        p.drawEllipse(sp, 1.8f, 1.8f);
    }
}

// ---------------------------------------------------------------------------
// 3C 273 system paint methods
// ---------------------------------------------------------------------------

static void c3273Bg(QPainter &p, int w, int h, QPointF c, float bR)
{
    // Warm yellow quasar glow
    QRadialGradient qg(c, bR * 8.f);
    qg.setColorAt(0.0f, QColor(50,40,10,60));
    qg.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,w,h), qg);
}

// Inner jet-base knot: moves along angled jet
void BlackHolePreviewWidget::paint3c273JetBlob(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.52f, H * 0.54f);
    const float bR = H * 0.072f;

    drawBackground(p, QColor(18,12,4), QColor(5,4,2));
    drawStarfield(p);
    c3273Bg(p, width(), height(), c, bR);

    drawContextBH(p, c, bR, true, bR * 2.6f, 0.20f, QColor(255,200,80));

    // Angled jet direction (35 deg from vertical)
    float jetAng = -35.f;
    float jetLen = H * 0.65f;
    drawJet(p, c, bR, jetAng, jetLen, bR * 0.45f, QColor(255,220,120), 0.55f);

    // Knot position: travels outward along jet and loops
    float phase = std::fmod(m_phase * 0.55f, 1.f);
    float jRad  = jetAng * PI / 180.f;
    float kDist = bR * 1.4f + phase * jetLen * 0.75f;
    QPointF knot(c.x() + kDist * std::sin(jRad), c.y() - kDist * std::cos(jRad));
    // Knot brightness fades as it moves out
    int ka = int(220 * (1.f - phase * 0.7f));
    QRadialGradient kg(knot, 7.f);
    kg.setColorAt(0.f, QColor(255,240,120,ka)); kg.setColorAt(1.f, Qt::transparent);
    p.setBrush(kg); p.setPen(Qt::NoPen);
    p.drawEllipse(knot, 7.f, 7.f);
    p.setBrush(QColor(255,250,180,ka)); p.drawEllipse(knot, 3.f, 3.f);
}

// BLR cloud: reverberation mapping, orbiting gas
void BlackHolePreviewWidget::paint3c273Blr(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.072f;

    drawBackground(p, QColor(16,11,4), QColor(4,3,2));
    drawStarfield(p);
    c3273Bg(p, width(), height(), c, bR);

    drawContextBH(p, c, bR, true, bR * 2.6f, 0.20f, QColor(255,200,80));

    float orbitR = bR * 3.2f;
    drawOrbitPath(p, c, orbitR, orbitR * 0.38f, 20.f, QColor(255,200,100), 0.40f, true);

    float ang = m_angle * PI / 180.f;
    QPointF cloud(c.x() + orbitR * std::cos(ang), c.y() + orbitR * 0.38f * std::sin(ang));
    // Photoionised gas blob, warm yellow
    QRadialGradient cg(cloud, 10.f);
    cg.setColorAt(0.f, QColor(255,210,100,200)); cg.setColorAt(1.f, Qt::transparent);
    p.setBrush(cg); p.setPen(Qt::NoPen);
    p.drawEllipse(cloud, 10.f, 10.f);
    // Reverberation time-lag indicator line
    float lag = std::fmod(m_phase * 0.8f, 1.f);
    QPen lagPen(QColor(255,220,100, int(lag * 120)), 1.f, Qt::DotLine);
    p.setPen(lagPen);
    p.drawLine(c, cloud);
}

// Stripped dwarf remnant: compact nucleus orbiting
void BlackHolePreviewWidget::paint3c273Dwarf(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(16,11,4), QColor(4,3,2));
    drawStarfield(p);
    c3273Bg(p, width(), height(), c, bR);

    float orbitA = bR * 5.0f, orbitB = orbitA * 0.55f;
    drawOrbitPath(p, c, orbitA, orbitB, -20.f, QColor(200,180,120), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.2f, 0.20f, QColor(255,200,80));

    float ang = m_angle * PI / 180.f;
    QPointF nuc(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    // Compact stripped nucleus: tiny diffuse galaxy remnant
    QRadialGradient ng(nuc, 10.f);
    ng.setColorAt(0.0f, QColor(255,230,160,180));
    ng.setColorAt(0.4f, QColor(200,180,120,100));
    ng.setColorAt(1.0f, Qt::transparent);
    p.setBrush(ng); p.setPen(Qt::NoPen);
    p.drawEllipse(nuc, 10.f, 7.f);
    p.setBrush(QColor(255,240,180)); p.drawEllipse(nuc, 2.5f, 2.5f);
}

// Close stellar orbit: radiation pressure perturbing orbit
void BlackHolePreviewWidget::paint3c273CloseStar(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.072f;

    drawBackground(p, QColor(16,11,4), QColor(4,3,2));
    drawStarfield(p);
    c3273Bg(p, width(), height(), c, bR);

    // Radiation pressure perturbs the orbit, eccentricity oscillates
    float ecc = 0.60f + 0.10f * std::sin(m_phase * 0.7f);
    float a = bR * 4.0f, b = a * (1.f - ecc);
    drawOrbitPath(p, c, a, b, m_angle * 0.05f, QColor(255,200,100), 0.50f, true);

    drawContextBH(p, c, bR, true, bR * 2.4f, 0.20f, QColor(255,200,80));

    float ang = m_angle * PI / 180.f;
    QPointF star(c.x() + a * std::cos(ang), c.y() + b * std::sin(ang));
    QRadialGradient sg(star, 5.f);
    sg.setColorAt(0.f, QColor(255,255,200,230)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 5.f, 5.f);
    p.setBrush(QColor(255,255,210)); p.drawEllipse(star, 2.8f, 2.8f);
}

// ---------------------------------------------------------------------------
// J0529-4351 system paint methods
// ---------------------------------------------------------------------------

static void j0529Bg(QPainter &p, int w, int h, QPointF c, float bR)
{
    QRadialGradient qg(c, bR * 9.f);
    qg.setColorAt(0.0f, QColor(10,40,30,70));
    qg.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,w,h), qg);
}

// Fast accretion blob: near-Eddington, disc wind visible
void BlackHolePreviewWidget::paintJ0529FastBlob(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(4,16,12), QColor(2,5,4));
    drawStarfield(p);
    j0529Bg(p, width(), height(), c, bR);

    // Disc wind arcs expanding outward
    for (int i = 0; i < 3; ++i) {
        float ph  = std::fmod(m_phase + float(i) * 0.85f, 2.f);
        float ar  = bR * 2.2f + ph * bR * 3.0f;
        int   a   = int(80 * (1.f - ph / 2.f));
        QPen wp(QColor(80,220,160, a), 1.5f);
        p.setPen(wp); p.setBrush(Qt::NoBrush);
        p.drawEllipse(c, ar, ar * 0.30f);
    }

    drawContextBH(p, c, bR, true, bR * 2.8f, 0.20f, QColor(80,220,140));

    // Fast orbiting blob
    float orbitR = bR * 2.0f;
    drawOrbitPath(p, c, orbitR, orbitR * 0.30f, 0.f, QColor(80,220,140), 0.45f, true);
    float ang = m_angle * 1.4f * PI / 180.f; // faster orbit, near-Eddington inner disc
    QPointF blob(c.x() + orbitR * std::cos(ang), c.y() + orbitR * 0.30f * std::sin(ang));
    QRadialGradient bg(blob, 8.f);
    bg.setColorAt(0.f, QColor(100,255,180,220)); bg.setColorAt(1.f, Qt::transparent);
    p.setBrush(bg); p.setPen(Qt::NoPen);
    p.drawEllipse(blob, 8.f, 8.f);
}

// UV-bright photoionised clump
void BlackHolePreviewWidget::paintJ0529UvClump(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(4,14,10), QColor(2,4,3));
    drawStarfield(p);
    j0529Bg(p, width(), height(), c, bR);

    drawContextBH(p, c, bR, true, bR * 2.6f, 0.20f, QColor(80,220,140));

    // UV ionisation cone, faint wedge from disc
    QPainterPath cone;
    float coneHalf = 28.f * PI / 180.f;
    float cLen = H * 0.42f;
    cone.moveTo(c);
    cone.lineTo(c.x() + cLen * std::cos(-PI/2.f - coneHalf),
                c.y() + cLen * std::sin(-PI/2.f - coneHalf));
    cone.lineTo(c.x() + cLen * std::cos(-PI/2.f + coneHalf),
                c.y() + cLen * std::sin(-PI/2.f + coneHalf));
    cone.closeSubpath();
    QLinearGradient coneGrad(c, QPointF(c.x(), c.y() - cLen));
    coneGrad.setColorAt(0.f, QColor(60,220,150,60));
    coneGrad.setColorAt(1.f, Qt::transparent);
    p.fillPath(cone, coneGrad);

    // Clump orbiting in ionised region
    float orbitR = bR * 3.4f;
    drawOrbitPath(p, c, orbitR, orbitR * 0.35f, 0.f, QColor(60,220,140), 0.35f, true);
    float ang = m_angle * PI / 180.f;
    QPointF clump(c.x() + orbitR * std::cos(ang), c.y() + orbitR * 0.35f * std::sin(ang));
    QRadialGradient cg(clump, 9.f);
    cg.setColorAt(0.f, QColor(80,255,180,200)); cg.setColorAt(1.f, Qt::transparent);
    p.setBrush(cg); p.setPen(Qt::NoPen);
    p.drawEllipse(clump, 9.f, 9.f);
}

// Actively disrupting star: mass stream, partially intact body
void BlackHolePreviewWidget::paintJ0529TdeStar(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.44f, H * 0.52f);
    const float bR = H * 0.065f;

    drawBackground(p, QColor(4,14,10), QColor(2,4,3));
    drawStarfield(p);
    j0529Bg(p, width(), height(), c, bR);

    // Tidal stream
    QPainterPath stream;
    float a = bR * 4.2f;
    bool started = false;
    for (int i = 0; i <= 140; ++i) {
        float t   = float(i) / 140.f;
        float ang = (m_angle * 0.22f + t * 200.f) * PI / 180.f;
        float r   = a * (1.05f - 0.48f * std::cos(ang));
        float sx  = c.x() + r * std::cos(ang);
        float sy  = c.y() + r * 0.42f * std::sin(ang);
        if (!started) { stream.moveTo(sx,sy); started = true; }
        else          { stream.lineTo(sx,sy); }
    }
    QPen sp(QColor(80,220,150,160), 2.0f); p.setPen(sp); p.setBrush(Qt::NoBrush);
    p.drawPath(stream);

    drawContextBH(p, c, bR, true, bR * 2.4f, 0.20f, QColor(80,220,140));

    // Star body, slightly elongated, partly disrupted
    float stAng = m_angle * 0.25f * PI / 180.f;
    QPointF star(c.x() + a * 0.50f * std::cos(stAng), c.y() + a * 0.42f * 0.50f * std::sin(stAng));
    p.save(); p.translate(star); p.rotate(m_angle * 0.25f);
    QRadialGradient sg(QPointF(0,0), 11.f);
    sg.setColorAt(0.f, QColor(120,255,180,220)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(QPointF(0,0), 11.f, 6.f);
    p.restore();
}

// Outer Keplerian gas stream
void BlackHolePreviewWidget::paintJ0529GasStream(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.060f;

    drawBackground(p, QColor(4,14,10), QColor(2,4,3));
    drawStarfield(p);
    j0529Bg(p, width(), height(), c, bR);

    // Multiple spiral arm strands of in-falling gas
    for (int arm = 0; arm < 2; ++arm) {
        QPainterPath armPath;
        bool started = false;
        float offset = float(arm) * 180.f;
        for (int i = 0; i <= 180; ++i) {
            float t   = float(i) / 180.f;
            float ang = (m_angle * 0.15f + offset + t * 270.f) * PI / 180.f;
            float r   = bR * (5.5f - t * 2.5f);
            float sx  = c.x() + r * std::cos(ang);
            float sy  = c.y() + r * 0.45f * std::sin(ang);
            if (!started) { armPath.moveTo(sx,sy); started = true; }
            else          { armPath.lineTo(sx,sy); }
        }
        int a = arm == 0 ? 160 : 100;
        p.setPen(QPen(QColor(60,200,130,a), 1.8f)); p.setBrush(Qt::NoBrush);
        p.drawPath(armPath);
    }
    drawContextBH(p, c, bR, true, bR * 2.2f, 0.20f, QColor(80,220,140));
}

// Infalling stellar cluster: J0529 context
void BlackHolePreviewWidget::paintJ0529Cluster(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.058f;

    drawBackground(p, QColor(4,14,10), QColor(2,4,3));
    drawStarfield(p);
    j0529Bg(p, width(), height(), c, bR);

    float decay = std::fmod(m_phase * 0.12f, 1.f);
    float orbitA = bR * (5.6f - decay * 1.2f);
    float orbitB = orbitA * 0.52f;
    drawOrbitPath(p, c, orbitA, orbitB, 0.f, QColor(80,200,140), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.20f, QColor(80,220,140));

    float cAng = m_angle * PI / 180.f;
    QPointF clCentre(c.x() + orbitA * std::cos(cAng), c.y() + orbitB * std::sin(cAng));
    std::mt19937 rng(55u);
    std::uniform_real_distribution<float> offs(-10.f, 10.f);
    for (int i = 0; i < 12; ++i) {
        QPointF sp(clCentre.x() + offs(rng), clCentre.y() + offs(rng));
        float br = 0.5f + float(rng() & 15) / 30.f;
        p.setBrush(QColor(120,255,180,int(br*200))); p.setPen(Qt::NoPen);
        p.drawEllipse(sp, 1.8f, 1.8f);
    }
}

// ---------------------------------------------------------------------------
// M87 system paint methods
// ---------------------------------------------------------------------------

static void m87Bg(QPainter &p, int w, int h, QPointF c, float bR)
{
    QRadialGradient qg(c, bR * 9.f);
    qg.setColorAt(0.0f, QColor(20,16,8,60));
    qg.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,w,h), qg);
}

// Jet-base knot moving along M87 jet
void BlackHolePreviewWidget::paintM87JetKnot(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.52f, H * 0.54f);
    const float bR = H * 0.075f;

    drawBackground(p, QColor(12,8,4), QColor(4,3,2));
    drawStarfield(p);
    m87Bg(p, width(), height(), c, bR);

    drawContextBH(p, c, bR, true, bR * 2.5f, 0.20f, QColor(240,160,60));

    // One-sided Doppler-boosted jet toward upper-left
    float jetAng = -142.f;
    float jetLen = H * 0.70f;
    drawJet(p, c, bR, jetAng, jetLen, bR * 0.40f, QColor(255,180,60), 0.65f);

    // Counter-jet faint
    drawJet(p, c, bR, jetAng + 180.f, jetLen * 0.22f, bR * 0.25f, QColor(200,140,50), 0.20f);

    // Moving knot on jet
    float phase  = std::fmod(m_phase * 0.45f, 1.f);
    float jRad   = jetAng * PI / 180.f;
    float kDist  = bR * 1.5f + phase * jetLen * 0.80f;
    QPointF knot(c.x() + kDist * std::cos(jRad), c.y() + kDist * std::sin(jRad));
    int ka = int(230 * (1.f - phase * 0.65f));
    QRadialGradient kg(knot, 8.f);
    kg.setColorAt(0.f, QColor(255,220,100,ka)); kg.setColorAt(1.f, Qt::transparent);
    p.setBrush(kg); p.setPen(Qt::NoPen);
    p.drawEllipse(knot, 8.f, 8.f);
    p.setBrush(QColor(255,240,160,ka)); p.drawEllipse(knot, 3.5f, 3.5f);
}

// Inner stellar orbit in M87 nucleus
void BlackHolePreviewWidget::paintM87InnerStar(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.075f;

    drawBackground(p, QColor(12,8,4), QColor(4,3,2));
    drawStarfield(p);
    m87Bg(p, width(), height(), c, bR);

    float orbitA = bR * 3.6f, orbitB = orbitA * 0.45f;
    drawOrbitPath(p, c, orbitA, orbitB, 30.f, QColor(220,180,100), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.4f, 0.20f, QColor(240,160,60));

    float ang = m_angle * PI / 180.f;
    QPointF star(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    QRadialGradient sg(star, 5.f);
    sg.setColorAt(0.f, QColor(255,255,210,230)); sg.setColorAt(1.f, Qt::transparent);
    p.setBrush(sg); p.setPen(Qt::NoPen);
    p.drawEllipse(star, 5.f, 5.f);
    p.setBrush(QColor(255,255,220)); p.drawEllipse(star, 2.8f, 2.8f);
}

// Hot X-ray gas shell / AGN feedback bubble
void BlackHolePreviewWidget::paintM87HotGas(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.070f;

    drawBackground(p, QColor(10,6,2), QColor(3,2,1));
    drawStarfield(p);
    m87Bg(p, width(), height(), c, bR);

    // Radio bubble cavity: expanding shell
    for (int i = 0; i < 3; ++i) {
        float phase = std::fmod(m_phase + float(i) * 0.9f, 2.5f);
        float ar    = bR * 2.8f + phase * bR * 2.5f;
        int   a     = int(70 * (1.f - phase / 2.5f));
        for (int j = 2; j >= 1; --j) {
            QPen bp(QColor(255,120,40, a / j), float(j) * 1.5f);
            p.setPen(bp); p.setBrush(Qt::NoBrush);
            p.drawEllipse(c, ar, ar * 0.80f);
        }
    }

    // Hot ICM gas fill, X-ray orange glow
    QRadialGradient xrayGlow(c, bR * 5.5f);
    xrayGlow.setColorAt(0.0f, QColor(255,80,20,35));
    xrayGlow.setColorAt(0.5f, QColor(200,50,10,18));
    xrayGlow.setColorAt(1.0f, Qt::transparent);
    p.fillRect(QRectF(0,0,W,H), xrayGlow);

    drawContextBH(p, c, bR, true, bR * 2.2f, 0.20f, QColor(240,160,60));
}

// Inspiraling globular cluster
void BlackHolePreviewWidget::paintM87Globular(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.068f;

    drawBackground(p, QColor(10,6,2), QColor(3,2,1));
    drawStarfield(p);
    m87Bg(p, width(), height(), c, bR);

    float decay  = std::fmod(m_phase * 0.10f, 1.f);
    float orbitA = bR * (5.8f - decay * 1.0f);
    float orbitB = orbitA * 0.62f;
    drawOrbitPath(p, c, orbitA, orbitB, 0.f, QColor(200,160,80), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.0f, 0.20f, QColor(240,160,60));

    float cAng = m_angle * PI / 180.f;
    QPointF clCentre(c.x() + orbitA * std::cos(cAng), c.y() + orbitB * std::sin(cAng));
    std::mt19937 rng(33u);
    std::uniform_real_distribution<float> offs(-10.f, 10.f);
    for (int i = 0; i < 14; ++i) {
        QPointF sp(clCentre.x() + offs(rng), clCentre.y() + offs(rng));
        float br = 0.5f + float(rng() & 15) / 30.f;
        p.setBrush(QColor(255,230,160,int(br*200))); p.setPen(Qt::NoPen);
        p.drawEllipse(sp, 1.8f, 1.8f);
    }
    // Globular core
    QRadialGradient gc(clCentre, 7.f);
    gc.setColorAt(0.f, QColor(255,240,180,200)); gc.setColorAt(1.f, Qt::transparent);
    p.setBrush(gc); p.drawEllipse(clCentre, 7.f, 7.f);
}

// Infalling stripped dwarf galaxy, potential IMBH binary
void BlackHolePreviewWidget::paintM87Dwarf(QPainter &p)
{
    const float W = width(), H = height();
    const QPointF c(W * 0.50f, H * 0.52f);
    const float bR = H * 0.068f;

    drawBackground(p, QColor(10,6,2), QColor(3,2,1));
    drawStarfield(p);
    m87Bg(p, width(), height(), c, bR);

    float orbitA = bR * 5.0f, orbitB = orbitA * 0.55f;
    drawOrbitPath(p, c, orbitA, orbitB, -15.f, QColor(180,140,80), 0.45f, true);

    drawContextBH(p, c, bR, true, bR * 2.2f, 0.20f, QColor(240,160,60));

    float ang = m_angle * PI / 180.f;
    QPointF nuc(c.x() + orbitA * std::cos(ang), c.y() + orbitB * std::sin(ang));
    // Stripped dwarf nucleus, diffuse elliptical halo + compact centre
    QRadialGradient halo(nuc, 13.f);
    halo.setColorAt(0.0f, QColor(255,220,140,160));
    halo.setColorAt(0.5f, QColor(200,170,100,80));
    halo.setColorAt(1.0f, Qt::transparent);
    p.setBrush(halo); p.setPen(Qt::NoPen);
    p.drawEllipse(nuc, 13.f, 9.f);
    // Compact IMBH candidate core
    p.setBrush(QColor(255,235,165)); p.drawEllipse(nuc, 3.f, 3.f);
}
