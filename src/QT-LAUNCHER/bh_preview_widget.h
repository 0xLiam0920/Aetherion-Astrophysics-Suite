#pragma once
#include <QWidget>
#include <QTimer>
#include <QString>
#include <QVector>

/**
 * BlackHolePreviewWidget
 *
 * Animated QPainter-based preview used in the Object Library detail panel.
 * Each PreviewStyle renders a unique astrophysics-accurate animated scene.
 *
 * ── Black-hole objects ────────────────────────────────────────────────────
 *   "sgra"          Sgr A*         RIAF, dense galactic star field
 *   "m87"           M87*           one-sided relativistic jet, photon ring
 *   "ton618"        TON 618        hyperluminous disc, powerful twin jets
 *   "3c273"         3C 273         yellow-white disc, angled optical jet
 *   "j0529"         J0529-4351     super-Eddington winds / arcs
 *   "cygnusx1"      Cygnus X-1     X-ray binary, large blue OB companion
 *   "gw150914"      GW150914       merger remnant, gravitational-wave ripples
 *   "gaiabh1"       Gaia BH1       dormant, G-type companion, lensing ring
 *
 * ── Research Scenarios ───────────────────────────────────────────────────
 *   "isco_unstable" ISCO 5M        inward spiral → plunge
 *   "isco_critical" ISCO 6M        marginally stable circle
 *   "isco_stable"   ISCO 7M        stable circular orbit
 *   "photon_sphere" Photon sphere  photon arcs near r=3M
 *   "radial_infall" Radial infall  straight geodesic, time-dilation labels
 *   "tidal_disrupt" TDE            star stretched to debris stream
 *   "pulsar_inspiral" Pulsar NS    decaying spiral, GW chirp rings, jets
 *
 * ── Sgr A* system ────────────────────────────────────────────────────────
 *   "sgra_s2"       S2 analog      eccentric precessing stellar orbit
 *   "sgra_s14"      S14-like       ultra-tight orbit, strong precession
 *   "sgra_irs16"    IRS 16 star    wide orbit, OB star, dusty field
 *   "sgra_gasclump" Gas clump      gas cloud, tidal stretching
 *   "sgra_smember"  S-member       generic S-cluster star
 *
 * ── TON 618 system ───────────────────────────────────────────────────────
 *   "ton618_blr"    BLR clump      inner gas blob in quasar glow
 *   "ton618_tdstar" Tidal star     stripped star, debris stream
 *   "ton618_gasfilm" Gas filament  sheared TDE debris filament
 *   "ton618_plunge" Plunging star  near-radial capture, no TDE flare
 *   "ton618_outerstar" Outer star  wide orbit, apsidal precession
 *   "ton618_cluster" Cluster       infalling globular, dynamical friction
 *
 * ── 3C 273 system ────────────────────────────────────────────────────────
 *   "3c273_jetblob" Jet blob       knot moving along optical jet
 *   "3c273_blr"     BLR cloud      reverberation-mapping cloud
 *   "3c273_dwarf"   Dwarf remnant  tidally stripped dwarf galaxy core
 *   "3c273_closestar" Close star   radiation-pressure-perturbed orbit
 *
 * ── J0529-4351 system ────────────────────────────────────────────────────
 *   "j0529_fastblob" Fast blob     high-velocity accretion gas, disc wind
 *   "j0529_uvclump" UV clump       photoionised UV-bright cloud
 *   "j0529_tdestar" TDE star       active tidal disruption, mass stream
 *   "j0529_gasstream" Gas stream   outer Keplerian stream feeding disc
 *   "j0529_cluster" NS cluster     infalling stellar cluster
 *
 * ── M87 system ───────────────────────────────────────────────────────────
 *   "m87_jetknot"   Jet knot       knot on 6-kpc relativistic jet
 *   "m87_innerstar" Inner star     nuclear star cluster member
 *   "m87_hotgas"    Hot gas shell  X-ray ICM shell, AGN feedback bubble
 *   "m87_globular"  Globular       inspiraling globular cluster
 *   "m87_dwarf"     Dwarf galaxy   plunging stripped dwarf, IMBH candidate
 *
 *   ""/"none"  hidden
 */
class BlackHolePreviewWidget : public QWidget
{
    Q_OBJECT

public:
    enum class PreviewStyle {
        NONE,
        // ── Black-hole objects ──────────────────────────────────────────────
        SGR_A,        // Sgr A*
        M87,          // M87*
        TON618,       // TON 618
        C3_273,       // 3C 273
        J0529,        // J0529-4351
        CYGNUS_X1,    // Cygnus X-1
        GW150914,     // GW150914 remnant
        GAIA_BH1,     // Gaia BH1
        GAIA_BH2,     // Gaia BH2
        GAIA_BH3,     // Gaia BH3
        V404_CYGNI,   // V404 Cygni
        PHOENIX_A,    // Phoenix A
        // ── Research Scenarios ─────────────────────────────────────────────
        ISCO_UNSTABLE,    // unstable ISCO test (5M)
        ISCO_CRITICAL,    // critical ISCO (6M)
        ISCO_STABLE,      // stable orbit (7M)
        PHOTON_SPHERE,    // photon sphere test
        RADIAL_INFALL,    // radial geodesic infall
        TIDAL_DISRUPTION, // tidal disruption event
        PULSAR_INSPIRAL,  // inspiraling neutron star / pulsar
        // ── Sgr A* system ──────────────────────────────────────────────────
        SGRA_S2,          // S2 analog
        SGRA_S14,         // S14-like close orbit
        SGRA_IRS16,       // IRS 16 cluster star
        SGRA_GASCLUMP,    // circumnuclear gas clump
        SGRA_SMEMBER,     // generic S-cluster member
        // ── TON 618 system ─────────────────────────────────────────────────
        TON618_BLR,       // inner BLR gas clump
        TON618_TDSTAR,    // tidally stripped star
        TON618_GASFILM,   // TDE debris filament
        TON618_PLUNGE,    // plunging star (direct capture)
        TON618_OUTERSTAR, // outer stellar orbit
        TON618_CLUSTER,   // infalling stellar cluster
        // ── 3C 273 system ──────────────────────────────────────────────────
        C3273_JETBLOB,    // inner jet-base knot
        C3273_BLR,        // broad-line region cloud
        C3273_DWARF,      // tidally stripped dwarf remnant
        C3273_CLOSESTAR,  // close stellar orbit
        // ── J0529-4351 system ──────────────────────────────────────────────
        J0529_FASTBLOB,   // fast accretion blob
        J0529_UVCLUMP,    // UV-bright photoionised clump
        J0529_TDESTAR,    // actively disrupting star
        J0529_GASSTREAM,  // outer Keplerian gas stream
        J0529_CLUSTER,    // infalling stellar cluster
        // ── M87 system ─────────────────────────────────────────────────────
        M87_JETKNOT,      // relativistic jet-base knot
        M87_INNERSTAR,    // nuclear star cluster member
        M87_HOTGAS,       // hot X-ray ICM shell / feedback bubble
        M87_GLOBULAR,     // inspiraling globular cluster
        M87_DWARF,        // infalling stripped dwarf galaxy
    };

    explicit BlackHolePreviewWidget(QWidget *parent = nullptr);

    void setStyle(PreviewStyle style);
    void setStyle(const QString &key);  // "sgra","m87","ton618","3c273","j0529",
                                        // "cygnusx1","gw150914","gaiabh1", or ""/"none"

protected:
    void paintEvent(QPaintEvent *) override;

private slots:
    void tick();

private:
    PreviewStyle  m_style = PreviewStyle::NONE;
    QTimer       *m_timer = nullptr;
    float         m_angle = 0.f;   // degrees — drives disc rotation and companion orbits
    float         m_phase = 0.f;   // generic animation phase (jets, winds, ripples)
    int           m_frame = 0;

    struct Star { float x, y, r, brightness; };
    QVector<Star> m_stars;    // standard sparse field
    QVector<Star> m_gcStars;  // dense field for galactic-centre scenes

    void seedStars();

    // Black-hole object paint methods
    void paintSgrA   (QPainter &p);
    void paintM87    (QPainter &p);
    void paintTon618 (QPainter &p);
    void paint3c273  (QPainter &p);
    void paintJ0529  (QPainter &p);
    void paintCygX1  (QPainter &p);
    void paintGW     (QPainter &p);
    void paintGaiaBH1(QPainter &p);
    void paintGaiaBH2(QPainter &p);
    void paintGaiaBH3(QPainter &p);
    void paintV404Cygni(QPainter &p);
    void paintPhoenixA(QPainter &p);

    // Research scenario paint methods
    void paintIscoUnstable  (QPainter &p);
    void paintIscoCritical  (QPainter &p);
    void paintIscoStable    (QPainter &p);
    void paintPhotonSphere  (QPainter &p);
    void paintRadialInfall  (QPainter &p);
    void paintTidalDisrupt  (QPainter &p);
    void paintPulsarInspiral(QPainter &p);

    // Sgr A* system paint methods
    void paintSgraS2      (QPainter &p);
    void paintSgraS14     (QPainter &p);
    void paintSgraIrs16   (QPainter &p);
    void paintSgraGasClump(QPainter &p);
    void paintSgraSMember (QPainter &p);

    // TON 618 system paint methods
    void paintTon618Blr      (QPainter &p);
    void paintTon618TdStar   (QPainter &p);
    void paintTon618GasFilm  (QPainter &p);
    void paintTon618Plunge   (QPainter &p);
    void paintTon618OuterStar(QPainter &p);
    void paintTon618Cluster  (QPainter &p);

    // 3C 273 system paint methods
    void paint3c273JetBlob  (QPainter &p);
    void paint3c273Blr      (QPainter &p);
    void paint3c273Dwarf    (QPainter &p);
    void paint3c273CloseStar(QPainter &p);

    // J0529-4351 system paint methods
    void paintJ0529FastBlob  (QPainter &p);
    void paintJ0529UvClump   (QPainter &p);
    void paintJ0529TdeStar   (QPainter &p);
    void paintJ0529GasStream (QPainter &p);
    void paintJ0529Cluster   (QPainter &p);

    // M87 system paint methods
    void paintM87JetKnot   (QPainter &p);
    void paintM87InnerStar (QPainter &p);
    void paintM87HotGas    (QPainter &p);
    void paintM87Globular  (QPainter &p);
    void paintM87Dwarf     (QPainter &p);

    // Shared drawing utilities
    void drawBackground (QPainter &p, QColor nearColor, QColor farColor);
    void drawStarfield  (QPainter &p, bool dense = false);

    // 3D disc: draw half the disc, clipped to upper (back) or lower (front) half of widget
    void drawDiscHalf(QPainter &p, QPointF c, float bhR, float outerR, float scaleY,
                      QColor inner, QColor mid, QColor outer,
                      float alphaScale, bool front);

    // BH sphere with depth shading and disc rim-light
    void drawBHSphere  (QPainter &p, QPointF c, float r, QColor rimColor);
    void drawPhotonRing(QPainter &p, QPointF c, float r, float scaleY, QColor col);

    // Single jet cone (angleDeg = 0 points up, 180 down)
    void drawJet(QPainter &p, QPointF c, float bhR,
                 float angleDeg, float length, float baseW,
                 QColor col, float alpha);

    // Companion star; addStream draws a dashed mass-transfer arc toward bhPos
    void drawCompanion(QPainter &p, QPointF c, float orbitR, float starR,
                       float orbitAngleDeg, QColor col,
                       bool addStream = false, QPointF bhPos = {});

    // Draw an elliptical orbit path (thin dashed line)
    void drawOrbitPath(QPainter &p, QPointF c, float a, float b, float angleDeg,
                       QColor col, float alpha, bool dashed = true);

    // Small context BH sphere + optional mini disc (for system-object scenes)
    void drawContextBH(QPainter &p, QPointF c, float bhR,
                       bool drawDisc, float discOuter, float scaleY,
                       QColor discColor);
};
