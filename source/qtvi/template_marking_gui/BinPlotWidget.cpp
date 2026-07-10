// ============================================================================
// BinPlotWidget.cpp
//
// ECG and PPG are drawn at a single per-paint pixels-per-sample scale
// (BinPlotWidget::pxPerSample()), chosen so the widest trace fills the
// cell the layout gives this widget. Every bin window therefore comes out
// the same on-screen length; the full PPG tail stays visible and shares
// the ECG's time axis within the widget.
//
// ECG/PPG alignment:
//   The ECG template is R-anchored (R sits at m_rPeakSample). The PPG
//   template is foot-anchored (its own average). To overlay them at the
//   true simultaneous instant, the PPG is shifted right so its foot lands
//   at R + m_ppgDelay, where m_ppgDelay is the measured R->foot transit
//   delay (in samples). ppgStartSample() encodes that origin; the PPG is
//   drawn from its real samples (no fabricated pad).
//
// Std band: when a per-sample std vector is available for the trace
// (covering at least the visible samples), the widget paints a
// translucent gray polygon between mean-std and mean+std underneath
// the line. The band uses the SAME lo/hi vertical range as the trace
// so it lines up with the line at every sample. Empty std => no band,
// just the line.
// ============================================================================
#include "BinPlotWidget.hpp"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
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

    //different colored STD bands - bluish for ecg and reddish for ppg
    inline const QColor color_stdband_ecg{ 90, 130, 220, 110 };   // light navy/blue
    inline const QColor color_stdband_ppg{ 220, 120, 130, 110 };  // soft red/pink


    // ECG markers (P, Q, Tb, Te) - blacks and dark blues, darkest to lightest.
    constexpr QColor kColorEcgP{ 0,   0,   0 };   // black
    constexpr QColor kColorEcgQBegin{ 20,  20,  60 };   // very dark navy
    constexpr QColor kColorEcgS{ 30, 35, 85 };   // dark navy, between Q and Tb
    constexpr QColor kColorEcgTBegin{ 40,  50, 110 };   // dark navy
    constexpr QColor kColorEcgTEnd{ 70,  90, 160 };   // medium navy

    // PPG markers (On, Pk, Dc, 50, En) - shades of red, darkest to lightest.
    constexpr QColor kColorPpgOnset{ 110,   0,   0 };  // dark red
    constexpr QColor kColorPpgPeak{ 180,   0,   0 };  // red
    constexpr QColor kColorPpgDicrotic{ 220,  50,  50 };  // medium red
    constexpr QColor kColorPpg50{ 235, 100, 100 };  // light red
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
        case BinPlotWidget::EcgP:        return kColorEcgP;
        case BinPlotWidget::EcgQBegin:   return kColorEcgQBegin;
        case BinPlotWidget::EcgSEnd:     return kColorEcgS;
        case BinPlotWidget::EcgTBegin:   return kColorEcgTBegin;
        case BinPlotWidget::EcgTEnd:     return kColorEcgTEnd;
        case BinPlotWidget::PpgOnset:    return kColorPpgOnset;
        case BinPlotWidget::PpgPeak:     return kColorPpgPeak;
        case BinPlotWidget::PpgDicrotic: return kColorPpgDicrotic;
        case BinPlotWidget::Ppg50:       return kColorPpg50;
        case BinPlotWidget::PpgEnd:      return kColorPpgEnd;
        }
        if (BinPlotWidget::markerIsAbp(m))     return kColorAbp[m - BinPlotWidget::AbpOnset];
        if (BinPlotWidget::markerIsArt(m))     return kColorArt[m - BinPlotWidget::ArtOnset];
        if (BinPlotWidget::markerIsArtPulm(m)) return kColorArtPulm[m - BinPlotWidget::ArtPulmOnset];
        return Qt::black;
    }
    const char* markerShortLabel(int m) {
        switch (m) {
        case BinPlotWidget::EcgP:        return "P peak";
        case BinPlotWidget::EcgQBegin:   return "Q beg";
        case BinPlotWidget::EcgSEnd:     return "S end";
        case BinPlotWidget::EcgTBegin:   return "T beg";
        case BinPlotWidget::EcgTEnd:     return "T end";
        case BinPlotWidget::PpgOnset:    return "PPG On";
        case BinPlotWidget::PpgPeak:     return "PPG Peak";
        case BinPlotWidget::PpgDicrotic: return "DN";
        case BinPlotWidget::Ppg50:       return "PPG 50% Upsl";
        case BinPlotWidget::PpgEnd:      return "PPG End";
        case BinPlotWidget::AbpOnset: return "aBP On";  case BinPlotWidget::AbpPeak: return "aBP Pk";
        case BinPlotWidget::AbpDicrotic: return "aBP DN"; case BinPlotWidget::Abp50: return "aBP 50%";
        case BinPlotWidget::AbpEnd: return "aBP End";
        case BinPlotWidget::ArtOnset: return "ART On";  case BinPlotWidget::ArtPeak: return "ART Pk";
        case BinPlotWidget::ArtDicrotic: return "ART DN"; case BinPlotWidget::Art50: return "ART 50%";
        case BinPlotWidget::ArtEnd: return "ART End";
        case BinPlotWidget::ArtPulmOnset: return "APul On"; case BinPlotWidget::ArtPulmPeak: return "APul Pk";
        case BinPlotWidget::ArtPulmDicrotic: return "APul DN"; case BinPlotWidget::ArtPulm50: return "APul 50%";
        case BinPlotWidget::ArtPulmEnd: return "APul End";
        }
        return "?";
    }

    // Compute the visible-range vertical bounds for a trace. If a matching
    // std vector is supplied, expand the range to include mean +- std at
    // every visible sample so the band fits inside the drawing area
    // without clipping. visN is clamped to the vector length so a caller
    // asking for more samples than exist (e.g. a stunted template from a
    // very short recording) can't seek past the end.
    void computeVisibleRange(const std::vector<double>& v,
        const std::vector<double>& sd,
        int visN,
        double& lo, double& hi)
    {
        const int n = std::min(visN, static_cast<int>(v.size()));
        if (n < 1) { lo = 0.0; hi = 1.0; return; }   // nothing to range over
        lo = *std::min_element(v.begin(), v.begin() + n);
        hi = *std::max_element(v.begin(), v.begin() + n);
        if ((int)sd.size() < n) return;
        for (int i = 0; i < n; ++i) {
            lo = std::min(lo, v[i] - sd[i]);
            hi = std::max(hi, v[i] + sd[i]);
        }
    }

    // Draw the gray ±std band at a fixed pixels-per-sample scale using
    // the supplied (lo, hi) range. The caller must use the SAME range
    // for the trace draw so the band and line agree vertically.
    void drawStdBand(QPainter& p, const std::vector<double>& v, const std::vector<double>& sd, double startPx, int mt, int ph,
        double pxPerSample, int visN, double lo, double hi, QColor color)
    {
        if (visN < 2 || (int)v.size() < 2) return;
        if ((int)sd.size() < visN) return;          // empty/mismatched => no band
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        QPainterPath band;
        // Top edge: mean + std, left to right.
        for (int i = 0; i < visN; ++i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - ((v[i] + sd[i]) - lo) / r * ph;
            if (i == 0) band.moveTo(x, y);
            else        band.lineTo(x, y);
        }
        // Bottom edge: mean - std, right to left, to close the polygon.
        for (int i = visN - 1; i >= 0; --i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - ((v[i] - sd[i]) - lo) / r * ph;
            band.lineTo(x, y);
        }
        band.closeSubpath();

        p.setPen(Qt::NoPen);
        p.setBrush(color);
        p.drawPath(band);
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

        QPainterPath path;
        for (int i = 0; i < visN; ++i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - (v[i] - lo) / r * ph;
            if (i == 0) path.moveTo(x, y);
            else        path.lineTo(x, y);
        }

        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
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

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ppgStd,
    const std::vector<double>& ecg,
    const std::vector<double>& ecgStd,
    int ecgP, int qBegin, int sEnd, int tBegin, int tEnd,
    int ppgOnset, int ppgPeak,
    int ppgDicrotic, int ppg50, int ppgEnd,
    double rPeakSample, double ppgDelay)
{
    m_ppg = ppg;
    m_ppgStd = ppgStd;
    m_ecg = ecg;
    m_ecgStd = ecgStd;
    m_markers[EcgP] = ecgP;
    m_markers[EcgQBegin] = qBegin;
    m_markers[EcgSEnd] = sEnd;
    m_markers[EcgTBegin] = tBegin;
    m_markers[EcgTEnd] = tEnd;
    m_markers[PpgOnset] = ppgOnset;
    m_markers[PpgPeak] = ppgPeak;
    m_markers[PpgDicrotic] = ppgDicrotic;
    m_markers[Ppg50] = ppg50;
    m_markers[PpgEnd] = ppgEnd;
    m_rPeakSample = rPeakSample;
    // R->foot transit delay and the foot's column in the foot-aligned PPG
    // template. Captured here (not read live from the marker) so dragging the
    // onset marker doesn't reposition the whole PPG trace.
    m_ppgDelay = ppgDelay;
    m_ppgFootIdx = (ppgOnset > 0) ? ppgOnset : 0;
    m_hasPPG = !ppg.empty();
    // Cache the per-trace visible counts so paint, hit-test, and
    // drag-clamp all see the same numbers.
    m_ecgVisibleN = std::max(static_cast<int>(m_ecg.size()), 2);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));

    // sizeHint depends on visible counts; notify layout.
    updateGeometry();
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setSampleRate(double hz) {
    if (m_sampleRate == hz) return;
    m_sampleRate = hz;
    update();
}
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
    const std::vector<double>& abpStd,
    const std::vector<double>& artStd,
    const std::vector<double>& artPulmStd)
{
    m_abp = abp;
    m_art = art;
    m_artPulm = artPulm;
    m_abpStd = abpStd;
    m_artStd = artStd;
    m_artPulmStd = artPulmStd;
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
    // PPG sample 0 sits at (R column + transit delay - foot column) so the
    // foot lands at R + delay. Clamped >= 0 to avoid drawing left of the axis.
    double s = m_rPeakSample + m_ppgDelay - static_cast<double>(m_ppgFootIdx);
    return (s > 0.0) ? s : 0.0;
}

int BinPlotWidget::totalSampleSpan() const {
    // ECG runs over samples [0, m_ecgVisibleN). The PPG is drawn shifted
    // right so its foot lands at R + delay, so it runs over
    // [0, ppgStartSample() + ppgVisibleN). The widest of the two is what
    // has to fit in the drawable width.
    const int pad = static_cast<int>(std::round(ppgStartSample()));
    const int ppgSpan = (m_ppgVisibleN > 0) ? pad + m_ppgVisibleN : 0;
    return std::max(std::max(m_ecgVisibleN, ppgSpan), 2);
}

double BinPlotWidget::pxPerSample() const {
    // Scale so the widest trace exactly fills the drawable width. ECG and
   // PPG share this one value, so they stay aligned within this widget;
   // different bins use different values, which is what makes every bin
   // window come out the same on-screen length.
    const int span = totalSampleSpan();
    const double drawW = std::max(1, width() - margin_left - margin_right);
    return (span > 1) ? (drawW / static_cast<double>(span - 1)) : 1.0;
}

int BinPlotWidget::sampleFromX(double x, bool isEcg) const {
    const double pps = pxPerSample();
    if (isEcg) {
        return static_cast<int>(std::round((x - margin_left) / pps));
    }
    else {
        const double ppgPx = margin_left + ppgStartSample() * pps;
        return static_cast<int>(std::round((x - ppgPx) / pps));
    }
}

double BinPlotWidget::xFromSample(int s, bool isEcg) const {
    const double pps = pxPerSample();
    if (isEcg) {
        return margin_left + s * pps;
    }
    else {
        const double ppgPx = margin_left + ppgStartSample() * pps;
        return ppgPx + s * pps;
    }
}

// For a given marker, resolve which trace vector bounds it, whether its
// group is currently visible, and whether it uses ECG x-geometry, plus the
// visible-sample count that bounds its markers. ECG uses its own visible
// window; PPG and all arterial channels are foot-anchored and use the same
// visiblePpgCount() rule (the whole trace), so arterial markers behave
// identically to PPG markers. Returns false if the trace is empty/absent.
bool BinPlotWidget::markerTrace(int m, const std::vector<double>*& vec,
    bool& isEcg, bool& visible, int& visN) const
{
    if (markerIsEcg(m)) {
        vec = &m_ecg; isEcg = true; visible = m_showEcgMarkers;
        visN = m_ecgVisibleN;
        return !m_ecg.empty();
    }
    isEcg = false;   // PPG and all arterial groups ride the foot-anchored geometry
    if (markerIsPpg(m)) {
        vec = &m_ppg; visible = m_showPpgMarkers;
        visN = visiblePpgCount(static_cast<int>(m_ppg.size()));
        return m_hasPPG && !m_ppg.empty();
    }
    if (markerIsAbp(m)) {
        vec = &m_abp; visible = m_showAbpMarkers;
        visN = visiblePpgCount(static_cast<int>(m_abp.size()));
        return !m_abp.empty();
    }
    if (markerIsArt(m)) {
        vec = &m_art; visible = m_showArtMarkers;
        visN = visiblePpgCount(static_cast<int>(m_art.size()));
        return !m_art.empty();
    }
    if (markerIsArtPulm(m)) {
        vec = &m_artPulm; visible = m_showArtPulmMarkers;
        visN = visiblePpgCount(static_cast<int>(m_artPulm.size()));
        return !m_artPulm.empty();
    }
    return false;
}

int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        const std::vector<double>* vec = nullptr;
        bool isEcg = false, visible = false;
        int visN = 0;
        if (!markerTrace(m, vec, isEcg, visible, visN)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        // Every marker clamps to its trace's visible window (same rule for
        // ECG, PPG, and all arterial channels).
        if (idx >= visN) continue;
        const double mx = xFromSample(idx, isEcg);
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
    const double pps = pxPerSample();   // one scale for both traces this paint

    p.fillRect(rect(), Qt::white);

    // Title
    p.setPen(Qt::black);
    { QFont f = p.font(); f.setPointSize(8); p.setFont(f); }
    p.drawText(margin_left, 11, QString("Bin %1  [%2]").arg(m_binIndex).arg(m_leadLabel));

    double yLo = 0, yHi = 0;
    computeVisibleRange(m_ecg, m_ecgStd, m_ecgVisibleN, yLo, yHi);

    // ---- Axes: frame + ticks + labels ----
    {
        const double xAxisY = margin_top + ph;       // == h - kMB
        const double yAxisX = margin_left;

        p.setPen(QColor(150, 150, 150));
        p.drawLine(QPointF(yAxisX, margin_top), QPointF(yAxisX, xAxisY));   // y-axis
        p.drawLine(QPointF(yAxisX, xAxisY), QPointF(w - margin_right, xAxisY));  // x-axis

        QFont af = p.font(); af.setPointSize(7); p.setFont(af);

        // X ticks: sample index across the full widget span.
        const int span = totalSampleSpan();
        const int nx = 4;
        for (int t = 0; t <= nx; ++t) {
            const int s = static_cast<int>(std::round((double)span * t / nx));
            double x = margin_left + s * pps;
            if (x > w - margin_right) x = w - margin_right;
            p.drawLine(QPointF(x, xAxisY), QPointF(x, xAxisY + 3));
            QString lbl;
            if (m_sampleRate > 0.0) {
                lbl = QString::number(s / m_sampleRate, 'f', 2);
                if (t == nx) lbl += " s";   // unit on the last tick only
            }
            else {
                lbl = QString::number(s);
            }
            p.drawText(QPointF(x - 3.0 * lbl.size(), xAxisY + 12), lbl);
        }

        // Y ticks: amplitude. Top of plot = hi, bottom = lo.
        for (int t = 0; t < 2; ++t) {
            const double y = (t == 0) ? margin_top : margin_top + ph;
            const double val = (t == 0) ? yHi : yLo;
            p.drawLine(QPointF(yAxisX - 3, y), QPointF(yAxisX, y));
            p.drawText(QPointF(1, y + 3), QString::number(val, 'e', 1));
        }

    }

    // -------- Arterial traces (ABP/ART/ART_PULM) --------
    // Foot-anchored like the PPG (shared ppgStartSample()). Each autoscaled
    // to its own range and gated by its own trace-visibility flag. A
    // translucent std band is drawn under the line when a matching std
    // vector is present (same treatment as ECG/PPG).
    {
        const double startPx = margin_left + ppgStartSample() * pps;
        struct ArtTrace {
            const std::vector<double>* v;
            const std::vector<double>* sd;
            QColor line;
            QColor band;
            bool show;
        };
        const ArtTrace arts[] = {
            { &m_abp,     &m_abpStd,     QColor(0, 115, 45),   QColor(80, 185, 120, 110),  m_showAbpTrace },     // green
            { &m_art,     &m_artStd,     QColor(140, 75, 185), QColor(180, 130, 215, 110), m_showArtTrace },     // purple
            { &m_artPulm, &m_artPulmStd, QColor(215, 135, 45), QColor(235, 175, 100, 110), m_showArtPulmTrace }, // orange
        };
        for (const auto& a : arts) {
            if (!a.show) continue;
            const std::vector<double>& v = *a.v;
            const int visN = static_cast<int>(v.size());
            if (visN < 2) continue;
            double lo = *std::min_element(v.begin(), v.end());
            double hi = *std::max_element(v.begin(), v.end());
            // Expand the vertical range to fit the std band (mean +/- std)
            // so it doesn't clip, mirroring computeVisibleRange for PPG.
            const std::vector<double>& sd = *a.sd;
            const bool haveStd = static_cast<int>(sd.size()) >= visN;
            if (haveStd) {
                for (int i = 0; i < visN; ++i) {
                    lo = std::min(lo, v[i] - sd[i]);
                    hi = std::max(hi, v[i] + sd[i]);
                }
            }
            if (haveStd)
                drawStdBand(p, v, sd, startPx, margin_top, ph, pps, visN, lo, hi, a.band);
            drawTraceFixedScale(p, v, startPx, margin_top, ph, pps,
                QPen(a.line, 1.3), visN, lo, hi);
        }
    }

    // -------- ECG --------
    if (m_showEcgTrace) {
        drawStdBand(p, m_ecg, m_ecgStd, margin_left, margin_top, ph, pps, m_ecgVisibleN, yLo, yHi, color_stdband_ecg);
        drawTraceFixedScale(p, m_ecg, margin_left, margin_top, ph, pps, QPen(kColorEcgTrace, 1.5), m_ecgVisibleN, yLo, yHi);
    }

    // -------- PPG --------
    if (m_showPpgTrace && m_hasPPG && !m_ppg.empty() && m_ppgVisibleN > 0) {
        // Real PPG, positioned so its foot sits at R + transit delay. No
        // fabricated lead-in: the trace begins at ppgStartSample() and
        // nothing is drawn to its left.
        const double startPx = margin_left + ppgStartSample() * pps;

        std::vector<double> ppgStdReal;
        if (static_cast<int>(m_ppgStd.size()) >= m_ppgVisibleN)
            ppgStdReal.assign(m_ppgStd.begin(),
                m_ppgStd.begin() + m_ppgVisibleN);

        double lo, hi;
        computeVisibleRange(m_ppg, ppgStdReal, m_ppgVisibleN, lo, hi);

        drawStdBand(p, m_ppg, ppgStdReal, startPx, margin_top, ph, pps,
            m_ppgVisibleN, lo, hi, color_stdband_ppg);
        drawTraceFixedScale(p, m_ppg, startPx, margin_top, ph,
            pps, QPen(kColorPpgTrace, 1.5), m_ppgVisibleN, lo, hi);
    }


    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        const std::vector<double>* vec = nullptr;
        bool isEcg = false, visible = false;
        int visN = 0;
        if (!markerTrace(m, vec, isEcg, visible, visN)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        if (idx >= visN) continue;
        double mx = xFromSample(idx, isEcg);
        QPen pen(markerColor(m), 2);
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(QPointF(mx + 2, margin_top + 8), markerShortLabel(m));
    }

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
        int m = markerAtX(e->position().x());
        if (m >= 0) {
            m_dragMarker = m;
            emit markerDragStarted(m_binIndex, m_leadIndex, m);
            return;
        }
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
    if (!markerTrace(m_dragMarker, vec, isEcg, visible, visN)) return;
    if (vec->empty()) return;
    int s = sampleFromX(e->position().x(), isEcg);
    // Every marker clamps to its trace's visible window (one rule for
    // ECG, PPG, and all arterial channels).
    s = std::clamp(s, 0, std::max(0, visN - 1));
    m_markers[m_dragMarker] = s;
    emit markerMoved(m_binIndex, m_leadIndex, m_dragMarker, s);
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragMarker = -1;
}