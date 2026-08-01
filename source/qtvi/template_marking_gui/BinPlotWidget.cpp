// ============================================================================
// BinPlotWidget.cpp
//
// ECG and PPG are drawn at a single per-paint pixels-per-sample scale
// (BinPlotWidget::pxPerSample()), chosen so the widest trace fills the
// cell the layout gives this widget. Every bin window therefore comes out
// the same on-screen length; the full PPG tail stays visible and shares
// the ECG's time axis within the widget.
//
// ECG/PPG alignment (Patch B/C):
//   Every channel's template is sliced from the SAME real-time window --
//   [t_R_i - pad, t_R_{i+1} + pad], driven by ch1.raw -- at that channel's
//   own sample rate. Sample 0 of the ECG, PPG, and every arterial channel
//   all correspond to the same real-time instant (pad seconds before the
//   first R). So they overlay by construction: ppgStartSample() = 0, no
//   R->foot delay adjustment needed. R sits at column pad*channelRate in
//   every template.
//
//   Every channel may run at its OWN sample rate (set via setChannelRate).
//   Geometry never assumes rates match: xFromSample/sampleFromX take a
//   (startSample, ratio) pair, where ratio = rateRatio(Channel) rescales
//   that channel's own sample index into ECG-equivalent sample units --
//   the frame's reference space (see pxPerSample/totalSampleSpan, both
//   sized from ECG alone). ratio is 1.0 whenever a rate is unknown or
//   matches ECG's, so this is a no-op for the historical case where every
//   channel happened to share one rate.
//
// Std band: when a per-sample std vector is available for the trace
// (covering at least the visible samples), the widget paints a
// translucent gray polygon between mean-std and mean+std underneath
// the line. The band uses the SAME lo/hi vertical range as the trace
// so it lines up with the line at every sample. Empty std => no band,
// just the line.
// ============================================================================
#include "BinPlotWidget.hpp"
#include "feature_marks.hpp"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QStringList>
#include <algorithm>
#include <cmath>

namespace {

    // ------------------------------------------------------------------
    // Colors
    //
    // Trace colors:
    //   ECG  - dark blue
    //   PPG  - dark red
    //
    // Marker colors are grouped by trace and lightened/darkened to
    // disambiguate the individual landmarks within each group:
    //   ECG markers - shades of black / dark blue
    //   PPG markers - shades of red
    //
    // Edit these constants in one place; markerColor() does the lookup.
    // ------------------------------------------------------------------
    constexpr QColor kColorEcgTrace{ 10,  20,  90 };   // dark navy blue
    constexpr QColor kColorPpgTrace{ 130,  10,  20 };   // dark red

    // Average-trace opacity: a little transparent so overlapping traces
    // (up to ECG + PPG + 3 arterial) reveal where they cross. ~60% opaque.
    constexpr int kTraceAlpha = 150;

    // Std bands: very light -- barely noticeable at a glance, visible on a
    // close look. Low alpha does that.
    inline const QColor color_iqrband_ecg{ 90, 130, 220, 38 };   // light navy/blue
    inline const QColor color_iqrband_ppg{ 220, 120, 130, 38 };  // soft red/pink

    // Apply kTraceAlpha to a base trace color.
    inline QColor withTraceAlpha(QColor c) { c.setAlpha(kTraceAlpha); return c; }


    // ECG markers (P peak, Q begin, R peak, S end, T peak, T end).
    constexpr QColor kColorEcgPPeak{ 0,   0,   0 };
    constexpr QColor kColorEcgQBegin{ 20,  20,  60 };
    constexpr QColor kColorEcgRPeak{ 15,  15,  40 };
    constexpr QColor kColorEcgS{ 30, 35, 85 };
    constexpr QColor kColorEcgTBegin{ 50,  60, 130 };
    constexpr QColor kColorEcgTEnd{ 70,  90, 160 };

    // PPG markers (On, P50, Pk, Dc, 2, En) - shades of red, darkest to lightest.
    constexpr QColor kColorPpgOnset{ 110,   0,   0 };  // dark red
    constexpr QColor kColorPpgP50{ 150,  20,  20 };  // dark red variant
    constexpr QColor kColorPpgPeak{ 180,   0,   0 };  // red
    constexpr QColor kColorPpgDicrotic{ 220,  50,  50 };  // medium red
    constexpr QColor kColorPpgPeak2{ 235, 100, 100 };  // light red (2nd/diastolic peak)
    constexpr QColor kColorPpgT80{ 210,  30,  70 };  // pink-red (80% upslope)
    constexpr QColor kColorPpgEnd{ 200,  60,  90 };  // dark pink-red

    // Arterial markers (ABP green, ART purple, ART_PULM orange),
    // darkest-to-lightest within a group.
    constexpr QColor kColorAbp[5] = {
        {0,80,30},{0,115,45},{30,150,75},{80,185,120},{130,210,160} };
    constexpr QColor kColorArt[5] = {
        {75,20,110},{105,35,150},{140,75,185},{170,120,210},{195,160,225} };
    constexpr QColor kColorArtPulm[5] = {
        {150,70,0},{190,100,15},{215,135,45},{230,165,90},{240,195,140} };

    QColor markerColor(int m) {
        switch (m) {
        case BinPlotWidget::EcgPPeak:    return kColorEcgPPeak;
        case BinPlotWidget::EcgQBegin:   return kColorEcgQBegin;
        case BinPlotWidget::EcgRPeak:    return kColorEcgRPeak;
        case BinPlotWidget::EcgSEnd:     return kColorEcgS;
        case BinPlotWidget::EcgTBegin:    return kColorEcgTBegin;
        case BinPlotWidget::EcgTEnd:     return kColorEcgTEnd;
        case BinPlotWidget::PpgOnset:    return kColorPpgOnset;
        case BinPlotWidget::PpgP50:      return kColorPpgP50;
        case BinPlotWidget::PpgPeak:     return kColorPpgPeak;
        case BinPlotWidget::PpgDicrotic: return kColorPpgDicrotic;
        case BinPlotWidget::PpgPeak2:    return kColorPpgPeak2;
        case BinPlotWidget::PpgT80:      return kColorPpgT80;
        case BinPlotWidget::PpgEnd:      return kColorPpgEnd;
        }
        if (BinPlotWidget::markerIsAbp(m))     return kColorAbp[m - BinPlotWidget::AbpOnset];
        if (BinPlotWidget::markerIsArt(m))     return kColorArt[m - BinPlotWidget::ArtOnset];
        if (BinPlotWidget::markerIsArtPulm(m)) return kColorArtPulm[m - BinPlotWidget::ArtPulmOnset];
        return Qt::black;
    }
    const char* markerShortLabel(int m) {
        switch (m) {
        case BinPlotWidget::EcgPBegin:   return "P beg";
        case BinPlotWidget::EcgPPeak:    return "P peak";
        case BinPlotWidget::EcgQBegin:   return "Q beg";
        case BinPlotWidget::EcgRPeak:    return "R peak";
        case BinPlotWidget::EcgSEnd:     return "S end";
        case BinPlotWidget::EcgTBegin:    return "T begin";
        case BinPlotWidget::EcgTEnd:     return "T end";
        case BinPlotWidget::PpgOnset:    return "PPG On";
        case BinPlotWidget::PpgP50:      return "PPG 50%";
        case BinPlotWidget::PpgPeak:     return "PPG Peak";
        case BinPlotWidget::PpgDicrotic: return "DN";
        case BinPlotWidget::PpgPeak2:    return "PPG Peak2";
        case BinPlotWidget::PpgT80:      return "T80";
        case BinPlotWidget::PpgEnd:      return "PPG End";
        case BinPlotWidget::AbpOnset: return "aBP On";  case BinPlotWidget::AbpPeak: return "aBP Pk";
        case BinPlotWidget::AbpDicrotic: return "aBP DN"; case BinPlotWidget::AbpPeak2: return "aBP Pk2";
        case BinPlotWidget::AbpEnd: return "aBP End";
        case BinPlotWidget::ArtOnset: return "ART On";  case BinPlotWidget::ArtPeak: return "ART Pk";
        case BinPlotWidget::ArtDicrotic: return "ART DN"; case BinPlotWidget::ArtPeak2: return "ART Pk2";
        case BinPlotWidget::ArtEnd: return "ART End";
        case BinPlotWidget::ArtPulmOnset: return "APul On"; case BinPlotWidget::ArtPulmPeak: return "APul Pk";
        case BinPlotWidget::ArtPulmDicrotic: return "APul DN"; case BinPlotWidget::ArtPulmPeak2: return "APul Pk2";
        case BinPlotWidget::ArtPulmEnd: return "APul End";
        }
        return "?";
    }

    void computeVisibleRange(const std::vector<double>& v, const std::vector<double>& /*sd*/, int visN, double& lo, double& hi) {
        /* Compute the visible range of the average template plot, ignoring the std band. 5% padding is added to accout for x marks
        */
        const int n = std::min(visN, static_cast<int>(v.size()));
        lo = 0.0; hi = 1.0;
        bool have = false;
        for (int i = 0; i < n; ++i) {
            const double m = v[i];
            if (std::isnan(m)) continue;
            if (!have) { lo = m; hi = m; have = true; }
            else { lo = std::min(lo, m); hi = std::max(hi, m); }
        }
        if (!have) { lo = 0.0; hi = 1.0; return; }
        lo -= 0.1; //account for X size
        hi += 0.1;
    }

    // Draw the gray ±std band at a fixed pixels-per-sample scale using
    // the supplied (lo, hi) range. The caller must use the SAME range
    // for the trace draw so the band and line agree vertically.
    void drawIqrBand(QPainter& p, const std::vector<double>& v, const std::vector<double>& sd, double startPx, int mt, int ph,
        double pxPerSample, int visN, double lo, double hi, QColor color)
    {
        if (visN < 2 || (int)v.size() < 2) return;
        if ((int)sd.size() < visN) return;          // empty/mismatched => no band
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        // One filled polygon per contiguous non-NaN run; NaN gaps break the
        // band. Prevents QPainter NaN-arcTo warnings from padded regions.
        p.setPen(Qt::NoPen);
        p.setBrush(color);
        int i = 0;
        while (i < visN) {
            while (i < visN && (std::isnan(v[i]) || std::isnan(sd[i]))) ++i;
            const int runStart = i;
            while (i < visN && !std::isnan(v[i]) && !std::isnan(sd[i])) ++i;
            const int runEnd = i;   // exclusive
            if (runEnd - runStart < 2) continue;

            QPainterPath band;
            for (int k = runStart; k < runEnd; ++k) {
                const double x = startPx + (double)k * pxPerSample;
                const double y = mt + ph - ((v[k] + sd[k]) - lo) / r * ph;
                if (k == runStart) band.moveTo(x, y); else band.lineTo(x, y);
            }
            for (int k = runEnd - 1; k >= runStart; --k) {
                const double x = startPx + (double)k * pxPerSample;
                const double y = mt + ph - ((v[k] - sd[k]) - lo) / r * ph;
                band.lineTo(x, y);
            }
            band.closeSubpath();
            p.drawPath(band);
        }
    }

    // Draw a trace at a fixed pixels-per-sample scale using the supplied
    // (lo, hi) vertical range. Caller-supplied range so that the band
    // and the line share an axis.
    void drawTraceFixedScale(QPainter& p, const std::vector<double>& v,
        double startPx, int mt, int ph,
        double pxPerSample, const QPen& pen, int visN,
        double lo, double hi)
    {
        if (visN < 2 || (int)v.size() < 2) return;
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        p.setPen(pen);
        p.setBrush(Qt::NoBrush);

        // NaN breaks the subpath: on the next real sample, restart with
        // moveTo. Prevents lineTo-to-NaN from producing arcTo NaN warnings.
        QPainterPath path;
        bool pending_move = true;
        for (int i = 0; i < visN; ++i) {
            if (std::isnan(v[i])) { pending_move = true; continue; }
            const double x = startPx + (double)i * pxPerSample;
            const double y = mt + ph - (v[i] - lo) / r * ph;
            if (pending_move) { path.moveTo(x, y); pending_move = false; }
            else { path.lineTo(x, y); }
        }
        p.drawPath(path);
    }
}

BinPlotWidget::BinPlotWidget(int binIndex, int leadIndex,
    const QString& leadLabel, QWidget* parent)
    : QWidget(parent), m_binIndex(binIndex), m_leadIndex(leadIndex),
    m_leadLabel(leadLabel)
{
    setMinimumHeight(100);
    setMouseTracking(true);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void BinPlotWidget::setChannelRate(Channel ch, double hz) {
    const size_t i = static_cast<size_t>(ch);
    if (m_rates[i] == hz) return;
    m_rates[i] = hz;
    update();
}

int BinPlotWidget::sampleFromX(double x, double startSample, double ratio) const {
    const double pps = pxPerSample();
    const double basePx = margin_left + startSample * pps;
    return static_cast<int>(std::round((x - basePx) / (pps * ratio)));
}

double BinPlotWidget::xFromSample(int s, double startSample, double ratio) const {
    const double pps = pxPerSample();
    return margin_left + startSample * pps + s * ratio * pps;
}

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ppgIqr,
    const std::vector<double>& ecg,
    const std::vector<double>& ecgIqr,
    int pPeak, int qBegin, int rPeak, int sEnd, int tPeak, int tEnd,
    int ppgOnset, int ppgP50, int ppgPeak,
    int ppgDicrotic, int ppgPeak2, int ppgT80, int ppgEnd,
    double rPeakSample,
    int nEcgBeats,
    int nPpgBeats,
    int ppgOnsetAuto,
    int ppgPeakAuto,
    int ppgPeak2Auto,
    bool ppgPeak2FoundAuto,
    int ppgDicroticAuto,
    bool ppgDicroticFoundAuto,
    int ppgEndAuto,
    bool ppgEndFoundAuto)
{
    m_nEcgBeats = nEcgBeats;
    m_nPpgBeats = nPpgBeats;
    m_ppg = ppg;
    m_ppgIqr = ppgIqr;
    m_ecg = ecg;
    m_ecgIqr = ecgIqr;
    m_markers[EcgPPeak] = pPeak;
    m_markers[EcgQBegin] = qBegin;
    m_markers[EcgRPeak] = rPeak;
    m_markers[EcgSEnd] = sEnd;
    m_markers[EcgTBegin] = tPeak;
    m_markers[EcgTEnd] = tEnd;
    m_markers[PpgOnset] = ppgOnset;
    m_markers[PpgP50] = ppgP50;
    m_markers[PpgPeak] = ppgPeak;
    m_markers[PpgDicrotic] = ppgDicrotic;
    m_markers[PpgPeak2] = ppgPeak2;
    m_markers[PpgT80] = ppgT80;
    m_markers[PpgEnd] = ppgEnd;
    // Fall back to the current bar position if no true auto value was
    // supplied, so the glyph snapshot always has something sensible to
    // seed from.
    m_ppgOnsetAuto = (ppgOnsetAuto >= 0) ? ppgOnsetAuto : ppgOnset;
    m_ppgPeakAuto = (ppgPeakAuto >= 0) ? ppgPeakAuto : ppgPeak;
    m_ppgPeak2Auto = (ppgPeak2Auto >= 0) ? ppgPeak2Auto : ppgPeak2;
    m_ppgPeak2FoundAuto = ppgPeak2FoundAuto;
    m_ppgDicroticAuto = (ppgDicroticAuto >= 0) ? ppgDicroticAuto : ppgDicrotic;
    m_ppgDicroticFoundAuto = ppgDicroticFoundAuto;
    m_ppgEndAuto = (ppgEndAuto >= 0) ? ppgEndAuto : ppgEnd;
    m_ppgEndFoundAuto = ppgEndFoundAuto;
    m_rPeakSample = rPeakSample;
    m_hasPPG = !ppg.empty();
    m_ecgVisibleN = std::max(static_cast<int>(m_ecg.size()), 2);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));
    captureGlyphSnapshot();   // frozen PPG glyphs = the auto fields, straight through
    updateGeometry();
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setState(State s) { m_state = s; update(); }

void BinPlotWidget::setShowEcgMarkers(bool show) {
    if (m_showEcgMarkers == show) return;
    m_showEcgMarkers = show;
    update();
}

void BinPlotWidget::setShowPpgMarkers(bool show) {
    if (m_showPpgMarkers == show) return;
    m_showPpgMarkers = show;
    update();
}

void BinPlotWidget::setShowAbpMarkers(bool show) {
    if (m_showAbpMarkers == show) return;
    m_showAbpMarkers = show; update();
}
void BinPlotWidget::setShowArtMarkers(bool show) {
    if (m_showArtMarkers == show) return;
    m_showArtMarkers = show; update();
}
void BinPlotWidget::setShowArtPulmMarkers(bool show) {
    if (m_showArtPulmMarkers == show) return;
    m_showArtPulmMarkers = show; update();
}

void BinPlotWidget::setShowEcgTrace(bool show) {
    if (m_showEcgTrace == show) return;
    m_showEcgTrace = show; update();
}
void BinPlotWidget::setShowPpgTrace(bool show) {
    if (m_showPpgTrace == show) return;
    m_showPpgTrace = show; update();
}
void BinPlotWidget::setShowAbpTrace(bool show) {
    if (m_showAbpTrace == show) return;
    m_showAbpTrace = show; update();
}
void BinPlotWidget::setShowArtTrace(bool show) {
    if (m_showArtTrace == show) return;
    m_showArtTrace = show; update();
}
void BinPlotWidget::setShowArtPulmTrace(bool show) {
    if (m_showArtPulmTrace == show) return;
    m_showArtPulmTrace = show; update();
}

void BinPlotWidget::setArterialTraces(const std::vector<double>& abp,
    const std::vector<double>& art,
    const std::vector<double>& artPulm,
    const std::vector<double>& abpIqr,
    const std::vector<double>& artIqr,
    const std::vector<double>& artPulmIqr)
{
    m_abp = abp;
    m_art = art;
    m_artPulm = artPulm;
    m_abpIqr = abpIqr;
    m_artIqr = artIqr;
    m_artPulmIqr = artPulmIqr;
    update();
}

void BinPlotWidget::setMarker(Marker m, int idx) {
    m_markers[m] = idx;
    update();
}

void BinPlotWidget::setBackgroundTraces(
    const std::vector<std::pair<std::vector<double>, QColor>>& traces)
{
    m_bgTraces = traces;
    update();
}

int BinPlotWidget::visibleN(bool isEcg) const {
    return isEcg ? m_ecgVisibleN : m_ppgVisibleN;
}

double BinPlotWidget::ppgStartSample() const {
    // Under Patch B, every channel is real-time-aligned by construction:
    // PPG and all arterial channels share sample 0 with the ECG template.
    // No shift needed.
    return 0.0;
}

int BinPlotWidget::totalSampleSpan() const {
    // Frame width is the ECG's visible extent. Pulse traces (PPG, ABP,
    // ART, ART_PULM) that extend past this are visually clipped by the
    // painter's clipRect. Their markers get their own visibility bound
    // via markerTrace()->pulseClipN() so no marker can land in the
    // clipped tail.
    return std::max(m_ecgVisibleN, 2);
}

// Number of pulse-trace samples that fit before the ECG's right edge, at
// the shared scale. Pulse sample i sits at shared position ppgStartSample()+i;
// keep those with (ppgStartSample()+i) <= m_ecgVisibleN. Also serves as the
// marker-visibility bound for pulse channels: markers past this can't be
// drawn or dragged (see markerTrace).
int BinPlotWidget::pulseClipN() const {
    const int clip = m_ecgVisibleN - static_cast<int>(std::llround(ppgStartSample()));
    return std::max(0, clip);
}

double BinPlotWidget::pxPerSample() const {
    // ONE scale shared by ECG and all pulse traces, so the R->foot offset is
    // a true horizontal distance. The frame spans the widest trace; every
    // trace is drawn at this same pps. Per-channel RATE differences are
    // handled separately via rateRatio()/xFromSample's ratio argument, not
    // here -- this stays purely a sample-count-to-pixel-width scale.
    const int span = totalSampleSpan();
    const double drawW = std::max(1, width() - margin_left - margin_right);
    return (span > 1) ? (drawW / static_cast<double>(span - 1)) : 1.0;
}

// For a given marker, resolve which trace vector bounds it, whether its
// group is currently visible, the visible-sample count that bounds its
// markers, and the ECG-equivalent rate ratio for its geometry. ECG uses
// its own visible window (ratio 1.0); PPG and all arterial channels are
// foot-anchored and rescaled by their own rateRatio(). Pulse markers are
// further capped at pulseClipN() so a marker sitting in the PPG (or
// arterial) tail past the ECG's right edge is neither drawn nor draggable
// -- the tail is visually clipped, and we don't want invisible markers
// out there that the user can't see or reach. Returns false if the trace
// is empty/absent.
bool BinPlotWidget::markerTrace(int m, const std::vector<double>*& vec,
    bool& isEcg, bool& visible, int& visN, double& ratio) const
{
    if (markerIsEcg(m)) {
        vec = &m_ecg; isEcg = true; visible = m_showEcgMarkers;
        visN = m_ecgVisibleN; ratio = rateRatio(Channel::Ecg);
        return !m_ecg.empty();
    }
    isEcg = false;   // PPG and all arterial groups ride the foot-anchored geometry
    const int pulseCap = pulseClipN();
    if (markerIsPpg(m)) {
        vec = &m_ppg; visible = m_showPpgMarkers;
        visN = std::min(visiblePpgCount(static_cast<int>(m_ppg.size())), pulseCap);
        ratio = rateRatio(Channel::Ppg);
        return m_hasPPG && !m_ppg.empty();
    }
    if (markerIsAbp(m)) {
        vec = &m_abp; visible = m_showAbpMarkers;
        visN = std::min(visiblePpgCount(static_cast<int>(m_abp.size())), pulseCap);
        ratio = rateRatio(Channel::Abp);
        return !m_abp.empty();
    }
    if (markerIsArt(m)) {
        vec = &m_art; visible = m_showArtMarkers;
        visN = std::min(visiblePpgCount(static_cast<int>(m_art.size())), pulseCap);
        ratio = rateRatio(Channel::Art);
        return !m_art.empty();
    }
    if (markerIsArtPulm(m)) {
        vec = &m_artPulm; visible = m_showArtPulmMarkers;
        visN = std::min(visiblePpgCount(static_cast<int>(m_artPulm.size())), pulseCap);
        ratio = rateRatio(Channel::ArtPulm);
        return !m_artPulm.empty();
    }
    return false;
}

int BinPlotWidget::glyphAtX(double x, int& outCol, bool& outIsEcg) const {
    // Hit-test the feature GLYPHS (the X marks), as opposed to markerAtX
    // which hit-tests the draggable bars. B2 focus mode selects landmarks by
    // clicking their glyph. Returns a Marker-enum routing id (so refreshFocus
    // can label/route it) or -1 if no glyph is near x. outCol receives the
    // glyph's sample column; outIsEcg whether it rides ECG or PPG geometry.
    if (!m_glyphs.valid) return -1;

    struct GlyphHit { int col; bool isEcg; int routeMarker; };
    const double ecgRatio = rateRatio(Channel::Ecg);
    const double ppgRatio = rateRatio(Channel::Ppg);

    // Map each ECG glyph to (column, geometry, routing Marker id). The R-peak
    // glyph (which has no draggable bar) routes to EcgRPeak so the focus
    // panel can show it. T-peak glyph routes to EcgTBegin so refreshFocus
    // sends it to the JT-side panel.
    const std::vector<GlyphHit> hits = {
        { m_glyphs.ecgPBegin, true,  EcgPBegin },
        { m_glyphs.ecgPPeak,  true,  EcgPPeak  },
        { m_glyphs.ecgQ,      true,  EcgQBegin },
        { m_glyphs.ecgRPeak,  true,  EcgRPeak  },
        { m_glyphs.ecgS,      true,  EcgSEnd   },
        { m_glyphs.ecgTPeak,  true,  EcgTBegin },
        { m_glyphs.ecgTend,   true,  EcgTEnd   },
        // PPG glyphs (B2 focus mode extended to PPG). All ride PPG geometry
        // (isEcg=false). Routed to their PPG Marker ids so refreshFocus can
        // pull the PPG template + label them.
        { m_glyphs.ppgFoot,   false, PpgOnset    },
        { m_glyphs.ppgP1,     false, PpgPeak     },
        { m_glyphs.ppgDic,    false, PpgDicrotic },
        { m_glyphs.ppgP2,     false, PpgPeak2    },
        { m_glyphs.ppgEnd,    false, PpgEnd      },
        { m_glyphs.ppgP50,    false, PpgP50      },
        { m_glyphs.ppgT80,    false, PpgT80      },
    };

    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (const auto& h : hits) {
        if (h.col < 0) continue;
        // NOTE: no visibility gate here (unlike markerAtX). The feature
        // glyphs are always drawn regardless of the marker-visibility
        // toggles, so they're always selectable.
        const double startSample = h.isEcg ? 0.0 : ppgStartSample();
        const double ratio = h.isEcg ? ecgRatio : ppgRatio;
        const double gx = xFromSample(h.col, startSample, ratio);
        const double d = std::abs(x - gx);
        if (d < bestDist) {
            bestDist = d; best = h.routeMarker;
            outCol = h.col; outIsEcg = h.isEcg;
        }
    }
    return best;
}

int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        if (m == EcgRPeak) continue;   // R is auto-only: no draggable bar
        if (m == PpgPeak) continue;    // systolic peak is auto-only (shown as X)
        if (m == PpgT80) continue;     // t80 is a reactive glyph now, not draggable
        if (m == PpgP50) continue;      // p50 is a reactive glyph now, not draggable
        const std::vector<double>* vec = nullptr;
        bool isEcg = false, visible = false;
        int visN = 0;
        double ratio = 1.0;
        if (!markerTrace(m, vec, isEcg, visible, visN, ratio)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        // Every marker clamps to its trace's visible window (same rule for
        // ECG, PPG, and all arterial channels).
        if (idx >= visN) continue;
        const double mx = xFromSample(idx, isEcg ? 0.0 : ppgStartSample(), ratio);
        const double d = std::abs(x - mx);
        if (d < bestDist) { bestDist = d; best = m; }
    }
    return best;
}

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    const int h = height();
    const int w = width();
    const int ph = h - margin_top - margin_bottom;
    const double pps = pxPerSample();
    p.fillRect(rect(), Qt::white);

    //title: channel, bin number, and number of beats in the bin
    p.setPen(QColor(150, 150, 150));
    { QFont f = p.font(); f.setPointSize(8); p.setFont(f); }
    QString beatSuffix;
    if (m_nEcgBeats > 0 || m_nPpgBeats > 0) {
        QStringList parts;
        if (m_nEcgBeats > 0) parts << QString("%1 ECG beats").arg(m_nEcgBeats);
        if (m_nPpgBeats > 0) parts << QString("%1 PPG beats").arg(m_nPpgBeats);
        beatSuffix = "  " + parts.join(", ");
    }
    if (m_binIndex == 0)
        p.drawText(QRectF(margin_left, 0, w - margin_left, 2.0 * margin_top),
            Qt::AlignBottom | Qt::AlignLeft,
            QString("  Bin %1  [%2 over time in seconds]\n%3")
            .arg(m_binIndex).arg(m_leadLabel).arg(beatSuffix));
    else
        p.drawText(margin_left, 11,
            QString("  Bin %1  [%2]%3")
            .arg(m_binIndex).arg(m_leadLabel).arg(beatSuffix));

    // Y-axis rules for the normalized traces:
    //   Y-max is FIXED at 1.0 (one decimal) for every panel so bins share
    //   a common upper reference.
    //   Y-min is autoscaled from the data, but always clipped at <=0.0 so
    //   the 0.0 major tick is always inside the frame.
    // Left axis (ECG) and right axis (pulse) each follow this rule.
    double yLo = 0, yHi = 0;
    computeVisibleRange(m_ecg, m_ecgIqr, m_ecgVisibleN, yLo, yHi);
    if (yLo > 0.0) yLo = 0.0;

    // Right axis range: shared by ALL pulse traces (PPG + arterial), so they
    // sit on one common normalized scale shown on the right. Range only over
    // the RETAINED samples (clipped at the ECG's right edge).
    const int pulseClip = pulseClipN();
    double pLo = 1e300, pHi = -1e300;
    auto mergePulse = [&](const std::vector<double>& v,
        const std::vector<double>& sd, int visN) {
            const int n = std::min(visN, pulseClip);
            if (n < 1) return;
            double lo, hi; computeVisibleRange(v, sd, n, lo, hi);
            pLo = std::min(pLo, lo); pHi = std::max(pHi, hi);
        };
    if (m_hasPPG && !m_ppg.empty()) mergePulse(m_ppg, m_ppgIqr, m_ppgVisibleN);
    mergePulse(m_abp, m_abpIqr, static_cast<int>(m_abp.size()));
    mergePulse(m_art, m_artIqr, static_cast<int>(m_art.size()));
    mergePulse(m_artPulm, m_artPulmIqr, static_cast<int>(m_artPulm.size()));
    if (pLo > pHi) { pLo = 0.0; pHi = 1.0; }
    if (pLo > 0.0) pLo = 0.0;

    // Each axis keeps its own y-min (autoscaled from the DATA it displays,
    // clipped at 0). The pulse axis already only considers samples up to
    // pulseClipN() = m_ecgVisibleN via the mergePulse lambda above, so the
    // pulse range reflects only what's visible on-screen inside the ECG
    // plot's right edge. Anything past the ECG width is ignored.

    // ---- Axes: frame + ticks + dual labeled Y-axes ----
    {
        const double xAxisY = margin_top + ph;       // == h - kMB
        const double yAxisL = margin_left;            // left  (ECG)
        const double yAxisR = w - margin_right;       // right (PPG/pulse)

        p.setPen(QColor(150, 150, 150));
        p.drawLine(QPointF(yAxisL, margin_top), QPointF(yAxisL, xAxisY));   // left  y-axis
        p.drawLine(QPointF(yAxisR, margin_top), QPointF(yAxisR, xAxisY));   // right y-axis
        p.drawLine(QPointF(yAxisL, xAxisY), QPointF(yAxisR, xAxisY));       // x-axis

        QFont af = p.font(); af.setPointSize(7); p.setFont(af);

        // X ticks: time across the full frame, at the shared scale.
        const int span = totalSampleSpan();
        p.setPen(QColor(150, 150, 150));
        for (int t = 0; t <= 4; ++t) {
            const int s = static_cast<int>(std::round((double)span * t / 4));
            double x = margin_left + s * pps;
            if (x > yAxisR) x = yAxisR;
            p.drawLine(QPointF(x, xAxisY), QPointF(x, xAxisY + 3));
            QString lbl = QString::number(s / channelRate(Channel::Ecg), 'f', 2);
            p.drawText(QPointF(x - 3.0 * lbl.size(), xAxisY + 12), lbl);
            // X-axis caption, centered under the tick numbers.
            if (m_binIndex == 0) {
                p.drawText(QRectF(yAxisL, xAxisY + 13.0, yAxisR - yAxisL, 11.0),
                    Qt::AlignHCenter | Qt::AlignTop, "time (s)");
            }
        }

        // Helper: y-pixel for a given data value on a given (lo, hi) scale.
        auto yPix = [&](double val, double lo, double hi) {
            const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
            return margin_top + ph - (val - lo) / r * ph;
            };

        // LEFT y-axis: ECG (normalized). Top = yHi, bottom = yLo, plus a
        // major tick at 0.0. Fixed 1-decimal format, no scientific notation.
        // Gray to match the x-axis ticks; labels right-aligned snug to the axis
        // (right edge stops just short of the tick marks so they don't overlap).
        p.setPen(QColor(150, 150, 150));
        {
            const double yTop = yPix(yHi, yLo, yHi);
            const double yBot = yPix(yLo, yLo, yHi);
            const double yZero = yPix(0.0, yLo, yHi);
            struct Tick { double val; double y; };
            const Tick ticks[] = {
                { yHi, yTop  },
                { 0.0, yZero },
                { yLo, yBot  }
            };
            for (const Tick& t : ticks) {
                p.drawLine(QPointF(yAxisL - 3, t.y), QPointF(yAxisL, t.y));
                const QString s = QString::number(t.val, 'f', 1);
                const double tw = p.fontMetrics().horizontalAdvance(s);
                p.drawText(QPointF(std::max(1.0, yAxisL - 4.0 - tw), t.y + 3.0), s);
                const double lx = std::max(0.0, yAxisL - 4.0 - tw);
                p.drawText(QPointF(lx, t.y + 3.0), s);
            }
        }

        // RIGHT y-axis: PPG / pulse (normalized). Same rules; gray, labels
        // left-aligned snug to the right of the tick marks.
        p.setPen(QColor(150, 150, 150));
        {
            const double yTop = yPix(pHi, pLo, pHi);
            const double yBot = yPix(pLo, pLo, pHi);
            const double yZero = yPix(0.0, pLo, pHi);
            struct Tick { double val; double y; };
            const Tick ticks[] = {
                { pHi, yTop  },
                { 0.0, yZero },
                { pLo, yBot  }
            };
            for (const Tick& t : ticks) {
                p.drawLine(QPointF(yAxisR, t.y), QPointF(yAxisR + 3, t.y));
                const QString rs = QString::number(t.val, 'f', 1);
                const double rtw = p.fontMetrics().horizontalAdvance(rs);
                double rx = yAxisR + 5.0;
                if (rx + rtw > w - 1.0) rx = w - 1.0 - rtw;
                p.drawText(QPointF(rx, t.y + 3.0), rs);
            }
        }

        // Rotated axis titles naming the tracing on each side. Gray to match
        // the ticks; kept at the outer edge, just past the number labels.
        QFont tf = p.font(); tf.setPointSize(7);
        p.setFont(tf);
        p.save();
        p.setPen(QColor(150, 150, 150));
        p.translate(9, margin_top + ph / 2.0);
        p.rotate(-90);
        p.drawText(QRectF(-ph / 2.0, -9, ph, 12), Qt::AlignCenter, "ECG (norm)");
        p.restore();
        p.save();
        p.setPen(QColor(150, 150, 150));
        p.translate(w - 9, margin_top + ph / 2.0);
        p.rotate(-90);
        p.drawText(QRectF(-ph / 2.0, -9, ph, 12), Qt::AlignCenter, "PPG (norm)");
        p.restore();
    }
    p.save();
    p.setClipRect(QRectF(margin_left, margin_top,
        w - margin_left - margin_right, ph));

    // -------- Arterial traces (ABP/ART/ART_PULM) --------
    // Foot-anchored like the PPG, drawn on the SHARED right-axis range
    // (pLo,pHi) so all pulse tracings sit on one mV scale. Each channel's
    // OWN rate ratio rescales pps so a channel running at a different rate
    // than ECG still draws at the correct real-time width.
    {
        const double startPx = margin_left + ppgStartSample() * pps;
        struct ArtTrace {
            const std::vector<double>* v;
            const std::vector<double>* sd;
            QColor line;
            QColor band;
            bool show;
            double ratio;
        };
        const ArtTrace arts[] = {
            { &m_abp,     &m_abpIqr,     QColor(0, 115, 45),   QColor(80, 185, 120, 38),  m_showAbpTrace,     rateRatio(Channel::Abp) },     // green
            { &m_art,     &m_artIqr,     QColor(140, 75, 185), QColor(180, 130, 215, 38), m_showArtTrace,     rateRatio(Channel::Art) },     // purple
            { &m_artPulm, &m_artPulmIqr, QColor(215, 135, 45), QColor(235, 175, 100, 38), m_showArtPulmTrace, rateRatio(Channel::ArtPulm) }, // orange
        };
        for (const auto& a : arts) {
            if (!a.show) continue;
            const std::vector<double>& v = *a.v;
            const int visN = std::min(static_cast<int>(v.size()), pulseClip);
            if (visN < 2) continue;
            const double effPps = pps * a.ratio;
            const std::vector<double>& sd = *a.sd;
            const bool haveStd = static_cast<int>(sd.size()) >= visN;
            if (haveStd)
                drawIqrBand(p, v, sd, startPx, margin_top, ph, effPps, visN, pLo, pHi, a.band);
            drawTraceFixedScale(p, v, startPx, margin_top, ph, effPps,
                QPen(withTraceAlpha(a.line), 1.3), visN, pLo, pHi);
        }
    }

    // -------- ECG (left axis) --------
    if (m_showEcgTrace) {
        drawIqrBand(p, m_ecg, m_ecgIqr, margin_left, margin_top, ph, pps, m_ecgVisibleN, yLo, yHi, color_iqrband_ecg);
        drawTraceFixedScale(p, m_ecg, margin_left, margin_top, ph, pps,
            QPen(withTraceAlpha(kColorEcgTrace), 1.5), m_ecgVisibleN, yLo, yHi);
    }

    // -------- PPG (right/shared axis) --------
    if (m_showPpgTrace && m_hasPPG && !m_ppg.empty() && m_ppgVisibleN > 0) {
        const double startPx = margin_left + ppgStartSample() * pps;
        const double effPps = pps * rateRatio(Channel::Ppg);
        const int ppgN = std::min(m_ppgVisibleN, pulseClipN());
        std::vector<double> ppgIqrReal;
        if (static_cast<int>(m_ppgIqr.size()) >= ppgN)
            ppgIqrReal.assign(m_ppgIqr.begin(), m_ppgIqr.begin() + ppgN);

        if (ppgN >= 2) {
            drawIqrBand(p, m_ppg, ppgIqrReal, startPx, margin_top, ph, effPps,
                ppgN, pLo, pHi, color_iqrband_ppg);
            drawTraceFixedScale(p, m_ppg, startPx, margin_top, ph,
                effPps, QPen(withTraceAlpha(kColorPpgTrace), 1.5), ppgN, pLo, pHi);
        }
    }
    p.restore();


    // Clip marker bars, labels, and fiducial glyphs to the plot area so
    // nothing (including the ~4px X/O glyphs at edge samples) draws past
    // the frame boundary.
    p.save();
    p.setClipRect(QRectF(margin_left, margin_top,
        w - margin_left - margin_right, ph));

    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        if (m == EcgRPeak) continue;   // R is auto-only: no draggable bar
        if (m == PpgPeak) continue;    // systolic peak is auto-only (shown as X)
        if (m == PpgT80) continue;     // t80 is a reactive glyph now, not draggable
        if (m == PpgP50) continue;      // p50 is a reactive glyph now, not draggable
        const std::vector<double>* vec = nullptr;
        bool isEcg = false, visible = false;
        int visN = 0;
        double ratio = 1.0;
        if (!markerTrace(m, vec, isEcg, visible, visN, ratio)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        if (idx >= visN) continue;
        double mx = xFromSample(idx, isEcg ? 0.0 : ppgStartSample(), ratio);
        QPen pen(markerColor(m), 2);
        pen.setStyle(markerIsBegin(m) ? Qt::DashLine : Qt::SolidLine);
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(QPointF(mx + 2, margin_top + 8), markerShortLabel(m));
    }

    drawFeatureGlyphs(p, yLo, yHi, pLo, pHi, ph);

    p.restore();


    if (m_state == State::BadPPG) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(margin_left, margin_top, w - margin_right, h - margin_bottom);
        p.drawLine(margin_left, h - margin_bottom, w - margin_right, margin_top);
    }

    if (m_state == State::BadR) {
        p.setPen(QPen(QColor(200, 0, 0), 2));
        QFont bf = p.font(); bf.setPointSize(14); bf.setBold(true); p.setFont(bf);
        p.drawText(rect(), Qt::AlignCenter, "BAD R");
    }
}


void BinPlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        // B2 focus mode: a click on a feature GLYPH selects that landmark for
        // the focus panel (this is the intended selection target -- glyphs,
        // not bars). Checked first so glyph-only landmarks like the R-peak
        // (which has no draggable bar) are selectable.
        int gcol = -1; bool gEcg = false;
        int gm = glyphAtX(e->position().x(), gcol, gEcg);
        if (gm >= 0) {
            emit landmarkSelected(m_binIndex, m_leadIndex, gm, gcol);
            // fall through: if a draggable bar is also under the cursor, still
            // begin a drag on it (glyph selection and bar drag can coexist).
        }
        int m = markerAtX(e->position().x());
        if (m >= 0) {
            m_dragMarker = m;
            emit markerDragStarted(m_binIndex, m_leadIndex, m);
            // B2 focus mode: bar-click selects focus ONLY for arterial
            // channels (ABP/ART/ART_PULM), which have draggable bars but no X
            // glyphs so the glyph path can't reach them. ECG and PPG have
            // glyphs -- their focus selection comes from glyphAtX above, and
            // their bars are for dragging only (avoids the redundant
            // double-trigger of selecting via both bar and glyph).
            if (gm < 0 && BinPlotWidget::markerIsArterial(m))
                emit landmarkSelected(m_binIndex, m_leadIndex, m, m_markers[m]);
            return;
        }
        if (gm >= 0) return;   // glyph selected but no bar to drag
    }

    if (e->button() == Qt::RightButton) {
        switch (m_state) {
        case State::Good:
            m_state = State::BadR;
            emit badRToggled(m_binIndex, m_leadIndex, true);
            break;
        case State::BadR:
            if (m_hasPPG) {
                m_state = State::BadPPG;
                emit badRToggled(m_binIndex, m_leadIndex, false);
                emit badPPGToggled(m_binIndex, true);
            }
            else {
                m_state = State::Good;
                emit badRToggled(m_binIndex, m_leadIndex, false);
            }
            break;
        case State::BadPPG:
            m_state = State::Good;
            emit badPPGToggled(m_binIndex, false);
            break;
        }
        update();
        return;
    }
}

void BinPlotWidget::mouseMoveEvent(QMouseEvent* e) {
    if (m_dragMarker < 0) return;
    const std::vector<double>* vec = nullptr;
    bool isEcg = false, visible = false;
    int visN = 0;
    double ratio = 1.0;
    if (!markerTrace(m_dragMarker, vec, isEcg, visible, visN, ratio)) return;
    if (vec->empty()) return;
    int s = sampleFromX(e->position().x(), isEcg ? 0.0 : ppgStartSample(), ratio);
    // Every marker clamps to its trace's visible window (one rule for
    // ECG, PPG, and all arterial channels).
    s = std::clamp(s, 0, std::max(0, visN - 1));
    m_markers[m_dragMarker] = s;
    emit markerMoved(m_binIndex, m_leadIndex, m_dragMarker, s);
    // Always recompute: ECG's Q-peak/S-peak track the live ECG bars by
    // design. The frozen PPG glyphs (foot/P1/dicrotic/P2/end) are sourced
    // from the auto members, not from m_markers, so they can't drift no
    // matter which bar is dragged. T80/P50 ARE reactive to the current PPG
    // bars, so this needs to run on PPG drags too for them to track live.
    captureGlyphSnapshot();
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragMarker = -1;
}
// ============================================================================
// Feature-glyph QC layer (formerly BinPlotGlyphs.cpp)
//
// Landmark positions are captured ONCE at setData() (captureGlyphSnapshot)
// from the current markers and frozen -- they do NOT follow subsequent
// drags. drawFeatureGlyphs() just renders the stored snapshot.
//
// Marks:
//   ECG - P peak, Q begin, Q peak (reactive), R peak, S peak (reactive),
//         S end, T peak, T end. All drawn as X.
//   PPG - foot, P1, P50, dicrotic notch, P2, end. Each of P1/P50/dic/P2
//         falls back to an O glyph at a sensible midpoint if the shape
//         detection didn't find a real landmark.
// ECG glyphs use the left-axis (yLo,yHi) scale + ECG x-geometry; PPG
// glyphs use the shared right-axis (pLo,pHi) scale + foot-anchored x.
// ============================================================================

void BinPlotWidget::captureGlyphSnapshot() {
    m_glyphs = GlyphSnapshot{};   // reset to all -1 / false

    if ((int)m_ecg.size() >= 3) {
        const int N = (int)m_ecg.size();
        // Frozen (own-bar) glyphs: read from the load-time auto positions, NOT
        // the live markers, so dragging a bar leaves its X where detection put
        // it. Fall back to the bar only if no auto value was supplied (-1).
        auto frozen = [&](int autoVal, int barVal) {
            const int v = (autoVal >= 0) ? autoVal : barVal;
            return (v >= 0 && v < N) ? v : -1;
            };
        m_glyphs.ecgPBegin = frozen(m_ecgAuto.pBegin, m_markers[EcgPBegin]);
        m_glyphs.ecgPPeak = frozen(m_ecgAuto.pPeak, m_markers[EcgPPeak]);
        m_glyphs.ecgQ = frozen(m_ecgAuto.qBegin, m_markers[EcgQBegin]);
        m_glyphs.ecgS = frozen(m_ecgAuto.sEnd, m_markers[EcgSEnd]);
        m_glyphs.ecgTend = frozen(m_ecgAuto.tEnd, m_markers[EcgTEnd]);

        // R is NEVER autodetected: draw it at the passed-in R marker
        // (m_markers[EcgRPeak] = r_peak_ch = r_col, straight from peak-finding
        // through alignment). No argmax, no compute_r_wave, no window search.
        m_glyphs.ecgRPeak = m_markers[EcgRPeak];

        // Responsive (bracketed) glyph: T peak = extremum between T-begin and
        // T-end. Recomputed live from the current markers so it tracks as those
        // bars move. Q-peak/S-peak remain unused (-1), as before.
        m_glyphs.ecgTPeak = FeatureMarks::compute_t_peak(
            m_ecg, m_markers[EcgTBegin], m_markers[EcgTEnd]);
        m_glyphs.ecgQPeak = -1;
        m_glyphs.ecgSPeak = -1;
    }

    if (m_hasPPG && (int)m_ppg.size() >= 3) {
        // Frozen glyphs: read DIRECTLY from the auto fields captured in
        // setData() -- these are the single source of truth (see
        // FeatureMarks::detect_ppg_fiducials), so there's nothing to
        // recompute and nothing that can drift from the auto-seeded
        // movable bars. No freeze flag needed either: these members never
        // change after setData() loads new data.
        m_glyphs.ppgFoot = m_ppgOnsetAuto;
        m_glyphs.ppgP1 = m_ppgPeakAuto;
        m_glyphs.ppgP2 = m_ppgPeak2Auto;        m_glyphs.ppgPeak2Found = m_ppgPeak2FoundAuto;
        m_glyphs.ppgDic = m_ppgDicroticAuto;    m_glyphs.ppgNotchFound = m_ppgDicroticFoundAuto;
        m_glyphs.ppgEnd = m_ppgEndAuto;         m_glyphs.ppgEndFound = m_ppgEndFoundAuto;

        // T80 and P50 are reactive, not frozen: always recomputed from the
        // CURRENT markers (peak is auto-only/effectively fixed; onset and
        // end are still draggable), same treatment as the ECG Q-peak/
        // S-peak reactive glyphs, just applied on the PPG side. Uses the
        // same amplitude_crossing formula the auto-detection itself uses,
        // so the two never disagree on what "80% down" means.
        const int pk = m_markers[PpgPeak];
        const int en = m_markers[PpgEnd];
        const int on = m_markers[PpgOnset];
        if (pk >= 0 && en > pk) m_glyphs.ppgT80 = FeatureMarks::amplitude_crossing(m_ppg, pk, en, 0.80);
        if (on >= 0 && pk > on) m_glyphs.ppgP50 = FeatureMarks::amplitude_crossing(m_ppg, on, pk, 0.50);
    }

    m_glyphs.valid = true;
}

void BinPlotWidget::drawFeatureGlyphs(QPainter& p,
    double yLo, double yHi, double pLo, double pHi, int ph) const
{
    if (!m_glyphs.valid) return;

    auto plotY = [&](double val, double lo, double hi) {
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
        return margin_top + ph - (val - lo) / r * ph;
        };
    auto glyph = [&](double x, double y) {   // opaque black "X"
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(Qt::black, 1.25));
        const double s = 2.0;
        p.drawLine(QPointF(x - s, y - s), QPointF(x + s, y + s));
        p.drawLine(QPointF(x - s, y + s), QPointF(x + s, y - s));
        };

    if (m_showEcgTrace && (int)m_ecg.size() >= 3) {
        const std::vector<double>& v = m_ecg;
        const int N = (int)v.size();
        double baseline = 0.0;
        {
            const int bhi = std::min(10, N);
            std::vector<double> bw;
            for (int i = 0; i < bhi; ++i) if (!std::isnan(v[i])) bw.push_back(v[i]);
            if (!bw.empty()) { std::sort(bw.begin(), bw.end()); baseline = bw[bw.size() / 2]; }
        }
        auto g = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? baseline : v[idx];
            glyph(xFromSample(idx, 0.0, rateRatio(Channel::Ecg)), plotY(val, yLo, yHi));
            };
        g(m_glyphs.ecgPBegin);  // P begin (frozen)
        g(m_glyphs.ecgPPeak);   // P wave
        g(m_glyphs.ecgQ);       // Q onset
        g(m_glyphs.ecgRPeak);   // R wave
        g(m_glyphs.ecgS);       // S end
        g(m_glyphs.ecgTPeak);   // T peak (between T begin/end)
        g(m_glyphs.ecgTend);    // T end (= marker)
    }

    if (m_showPpgTrace && m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        auto g = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            glyph(xFromSample(idx, ppgStartSample(), rateRatio(Channel::Ppg)), plotY(val, pLo, pHi));
            };
        auto circ = [&](int idx) {
            if (idx < 0 || idx >= N) return;
            const double val = std::isnan(v[idx]) ? pLo : v[idx];
            const double x = xFromSample(idx, ppgStartSample(), rateRatio(Channel::Ppg));
            const double y = plotY(val, pLo, pHi);
            p.setBrush(Qt::NoBrush);
            p.setPen(QPen(Qt::black, 1.8));
            p.drawEllipse(QPointF(x, y), 4.0, 4.0);
            };
        g(m_glyphs.ppgFoot);
        g(m_glyphs.ppgP50);   // reactive: 50% onset->peak, always a computed X
        g(m_glyphs.ppgP1);
        if (m_glyphs.ppgNotchFound) g(m_glyphs.ppgDic);
        else                        circ(m_glyphs.ppgDic);
        if (m_glyphs.ppgPeak2Found) g(m_glyphs.ppgP2);
        else                        circ(m_glyphs.ppgP2);
        if (m_glyphs.ppgEndFound)   g(m_glyphs.ppgEnd);
        else                        circ(m_glyphs.ppgEnd);
        g(m_glyphs.ppgT80);   // reactive: 80% peak->end, always a computed X
    }
}