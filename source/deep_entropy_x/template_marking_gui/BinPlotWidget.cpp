/*
* @brief BinPlotWidget.cpp
*
* ECG and PPG are drawn at a single per-paint pixels-per-sample scale
* (BinPlotWidget::pxPerSample()), chosen so the widest trace fills the
* cell the layout gives this widget. Every bin window therefore comes out
* the same on-screen length; the full PPG tail stays visible and shares
* the ECG's time axis within the widget.
*
* ECG/PPG alignment (Patch B/C):
*   Every channel's template is sliced from the SAME real-time window --
*   [t_R_i - pad, t_R_{i+1} + pad], driven by ch1.raw -- at that channel's
*   own sample rate. Sample 0 of the ECG, PPG, and every arterial channel
*   all correspond to the same real-time instant (pad seconds before the
*   first R). So they overlay by construction: ppgStartSample() = 0, no
*   R->foot delay adjustment needed. R sits at column pad*channelRate in
*   every template.

* Every channel may run at its OWN sample rate (set via setChannelRate).
*   Geometry never assumes rates match: xFromSample/sampleFromX take a
*   (startSample, ratio) pair, where ratio = rateRatio(Channel) rescales
*   that channel's own sample index into ECG-equivalent sample units --
*   the frame's reference space (see px_per_sample/totalSampleSpan, both
*   sized from ECG alone). ratio is 1.0 whenever a rate is unknown or
*   matches ECG's, so this is a no-op for the historical case where every
*   channel happened to share one rate.
*
* ============================================================================
* Glyphs and bars
*
* A BAR is a draggable vertical line the operator positions. A GLYPH is a
* small mark the widget draws and the operator cannot touch -- markerAtX()
* never hit-tests glyphs, so a click can neither select nor move one.
*
* Glyphs come in two flavours, distinguished by when they are computed:
*
*   FROZEN   -- captured once per seeding pass by setAuto()
*               (captureGlyphSnapshot) from the bin's *_auto columns.
*               Does NOT follow subsequent drags, so dragging a bar leaves
*               its glyph where detection put it. That difference is what
*               makes the paired _autodetect / _user CSV columns meaningful.
*   REACTIVE -- never stored. reactiveGlyphs() recomputes it from the current
*               bars at every paint, via the same FeatureMarks functions the
*               CSV/bin writers use, so screen and files cannot disagree.
*
* Three glyph shapes, three meanings:
*
*   X       the detector found a real landmark
*   O       the detector fell back to a placeholder position (notch and
*           diastolic peak only -- they are the two that carry a *_found flag)
*   dash    an auto-only derivative landmark, coloured by derivative order:
*           VPG blue, APG green, JPG amber.
*
* Marks drawn:
*   ECG  - P begin, P peak, Q onset, R peak, S end, T peak (reactive), T end.
*          All X.
*   PPG  - foot, T50 (reactive), systolic peak, dicrotic notch, diastolic
*          peak, pulse end, T80 (reactive). X, except notch and diastolic
*          peak which fall back to O.
*   VPG  - u, v, w                      | dashes, behind the
*   APG  - a, b, c, d, e, f             | "Show PPG Derivative Markers"
*   JPG  - p1, p2                       | checkbox
*
* A derivative landmark reads -1 when it genuinely does not exist in the
* pulse (c, d and p2 are the common cases), and -1 simply draws nothing --
* there is no placeholder and no flag.
*
* ECG glyphs use the left-axis (yLo,yHi) scale + ECG x-geometry; PPG and
* derivative glyphs use the shared right-axis (pLo,pHi) scale + the pulse
* x-geometry.
*/

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
    // Edit these constants in one place; marker_color() does the lookup.
    // ------------------------------------------------------------------
    constexpr QColor ecg_trace_color{ 10,  20,  90 };   // dark navy blue
    constexpr QColor ppg_trace_color{ 130,  10,  20 };   // dark red

    // Average-trace opacity: a little transparent so overlapping traces
    // (up to ECG + PPG + 3 arterial) reveal where they cross. ~60% opaque.
    constexpr int trace_alpha = 150;

    // Std bands: very light -- barely noticeable at a glance, visible on a
    // close look. Low alpha does that.
    inline const QColor color_iqrband_ecg{ 90, 130, 220, 38 };   // light navy/blue
    inline const QColor color_iqrband_ppg{ 220, 120, 130, 38 };  // soft red/pink

    // Apply trace_alpha to a base trace color.
    inline QColor with_trace_alpha(QColor c) { c.setAlpha(trace_alpha); return c; }


    // ECG markers (P peak, Q begin, R peak, S end, T peak, T end).
    constexpr QColor ecg_p_peak_color{ 0,   0,   0 };
    constexpr QColor ecg_q_begin_color{ 20,  20,  60 };
    constexpr QColor ecg_r_peak_color{ 15,  15,  40 };
    constexpr QColor ecg_s_color{ 30, 35, 85 };
    constexpr QColor ecg_t_begin_color{ 50,  60, 130 };
    constexpr QColor ecg_t_end_color{ 70,  90, 160 };

    // PPG bar markers - shades of red, darkest to lightest.
    constexpr QColor ppg_onset_color{ 110,   0,   0 };  // dark red
    constexpr QColor p50_color{ 150,  20,  20 };  // dark red variant
    constexpr QColor ppg_peak_color{ 180,   0,   0 };  // red
    constexpr QColor ppg_dicrotic_color{ 220,  50,  50 };  // medium red
    constexpr QColor ppg_peak2_color{ 235, 100, 100 };  // light red (2nd/diastolic peak)
    constexpr QColor ppg_t80_color{ 210,  30,  70 };  // pink-red (80% downslope)
    constexpr QColor ppg_end_color{ 200,  60,  90 };  // dark pink-red
    // Derivative marks, one colour per derivative order. Deliberately outside
    // the PPG reds above: a dash is an auto-only landmark with no bar, and the
    // hue says so at a glance.
    constexpr QColor vpg_mark_color{ 20,  80, 170 };  // blue   (PPG')
    constexpr QColor apg_mark_color{ 20, 120,  60 };  // green  (PPG'')
    constexpr QColor jpg_mark_color{ 190, 110,   0 };  // amber  (PPG''')

    // Glyph geometry, in pixels. marker_half_size is the X's half-extent on
    // both axes; dash_half_width is the dash's half-length. The dash is longer
    // because it has one stroke to the X's two and needs the extra reach to
    // stay findable where the trace is steep.
    constexpr double marker_half_size = 2.0;
    constexpr double dash_half_width = 5.0;
    constexpr double marker_pen_size = 1.25;
    constexpr double marker_circle_radius = 4.0;
    constexpr double marker_circle_pen = 1.8;

    // Arterial markers (ABP green, ART purple, ART_PULM orange),
    // darkest-to-lightest within a group.
    constexpr QColor abp_marker_colors[5] = { {0,80,30},{0,115,45},{30,150,75},{80,185,120},{130,210,160} };
    constexpr QColor art_marker_colors[5] = { {75,20,110},{105,35,150},{140,75,185},{170,120,210},{195,160,225} };
    constexpr QColor art_pulm_marker_colors[5] = { {150,70,0},{190,100,15},{215,135,45},{230,165,90},{240,195,140} };

    QColor marker_color(int m) {
        switch (m) {
        case BinPlotWidget::EcgPPeak:    return ecg_p_peak_color;
        case BinPlotWidget::EcgQBegin:   return ecg_q_begin_color;
        case BinPlotWidget::EcgRPeak:    return ecg_r_peak_color;
        case BinPlotWidget::EcgSEnd:     return ecg_s_color;
        case BinPlotWidget::EcgTBegin:    return ecg_t_begin_color;
        case BinPlotWidget::EcgTEnd:     return ecg_t_end_color;
        case BinPlotWidget::PpgOnset:    return ppg_onset_color;
        case BinPlotWidget::PpgT50:      return p50_color;
        case BinPlotWidget::PpgPeak:     return ppg_peak_color;
        case BinPlotWidget::PpgDicrotic: return ppg_dicrotic_color;
        case BinPlotWidget::PpgPeak2:    return ppg_peak2_color;
        case BinPlotWidget::PpgT80:      return ppg_t80_color;
        case BinPlotWidget::PpgEnd:      return ppg_end_color;
        }
        if (BinPlotWidget::markerIsAbp(m))     return abp_marker_colors[m - BinPlotWidget::AbpOnset];
        if (BinPlotWidget::markerIsArt(m))     return art_marker_colors[m - BinPlotWidget::ArtOnset];
        if (BinPlotWidget::markerIsArtPulm(m)) return art_pulm_marker_colors[m - BinPlotWidget::ArtPulmOnset];
        return Qt::black;
    }
    const char* marker_short_label(int m) {
        switch (m) {
        case BinPlotWidget::EcgPBegin:   return "P beg";
        case BinPlotWidget::EcgPPeak:    return "P peak";
        case BinPlotWidget::EcgQBegin:   return "Q beg";
        case BinPlotWidget::EcgRPeak:    return "R peak";
        case BinPlotWidget::EcgSEnd:     return "S end";
        case BinPlotWidget::EcgTBegin:    return "T begin";
        case BinPlotWidget::EcgTEnd:     return "T end";
        case BinPlotWidget::PpgOnset:    return "PPG On";
        case BinPlotWidget::PpgT50:      return "PPG 50%";
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

    void compute_visible_range(const std::vector<double>& v, const std::vector<double>& /*sd*/, int visN, double& lo, double& hi) {
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
    void draw_iqr_band(QPainter& p, const std::vector<double>& v, const std::vector<double>& sd, double startPx, int mt, int ph,
        double px_per_sample, int visN, double lo, double hi, QColor color)
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
                const double x = startPx + (double)k * px_per_sample;
                const double y = mt + ph - ((v[k] + sd[k]) - lo) / r * ph;
                if (k == runStart) band.moveTo(x, y); else band.lineTo(x, y);
            }
            for (int k = runEnd - 1; k >= runStart; --k) {
                const double x = startPx + (double)k * px_per_sample;
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
    void draw_trace_fixed_scale(QPainter& p, const std::vector<double>& v,
        double startPx, int mt, int ph,
        double px_per_sample, const QPen& pen, int visN,
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
            const double x = startPx + (double)i * px_per_sample;
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
    // Every marker starts unset. Anything the seeding pass doesn't provide
    // stays -1 and is simply not drawn -- never an indeterminate column.
    std::fill(std::begin(m_markers), std::end(m_markers), -1);
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

double BinPlotWidget::xFromSample(double s, double startSample, double ratio) const {
    const double pps = pxPerSample();
    return margin_left + startSample * pps + s * ratio * pps;
}

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ppgIqr,
    const std::vector<double>& ecg,
    const std::vector<double>& ecgIqr,
    double rPeakSample,
    int nEcgBeats,
    int nPpgBeats)
{
    m_nEcgBeats = nEcgBeats;
    m_nPpgBeats = nPpgBeats;
    m_ppg = ppg;
    m_ppgIqr = ppgIqr;
    m_ecg = ecg;
    m_ecgIqr = ecgIqr;
    m_rPeakSample = rPeakSample;
    m_hasPPG = !ppg.empty();
    m_ecgVisibleN = std::max(static_cast<int>(m_ecg.size()), 2);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));
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

// Derivative marks (u/v/w, a-f, p1/p2) have their own toggle: they are
// auto-only and numerous, so they stay off while the operator works the bars.
void BinPlotWidget::setShowPpgDerivMarkers(bool show) {
    if (m_showPpgDerivMarkers == show) return;
    m_showPpgDerivMarkers = show;
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

// Reactive glyphs. No arithmetic lives here -- FeatureMarks owns the formula,
// and the CSV/bin writers call the same functions with the bar set they're
// reporting on, so the screen and the files agree by construction.
BinPlotWidget::Reactive BinPlotWidget::reactiveGlyphs() const {
    Reactive r;
    r.ecgTPeak = FeatureMarks::reactive_ecg(
        m_ecg, m_markers[EcgTBegin], m_markers[EcgTEnd]).t_peak;
    if (m_hasPPG) {
        const FeatureMarks::ReactivePpg p = FeatureMarks::reactive_ppg(
            m_ppg, m_markers[PpgOnset], m_markers[PpgPeak], m_markers[PpgEnd]);
        r.ppgT50 = p.t50;
        r.ppgT80 = p.t80;
    }
    return r;
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


int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        if (m == EcgRPeak) continue;   // R is auto-only: no draggable bar
        if (m == PpgPeak) continue;    // systolic peak is auto-only (shown as X)
        if (m == PpgT80) continue;     // t80 is a reactive glyph now, not draggable
        if (m == PpgT50) continue;      // p50 is a reactive glyph now, not draggable
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
    compute_visible_range(m_ecg, m_ecgIqr, m_ecgVisibleN, yLo, yHi);
    if (yLo > 0.0) yLo = 0.0;

    // Right axis range: shared by ALL pulse traces (PPG + arterial), so they
    // sit on one common normalized scale shown on the right. Range only over
    // the RETAINED samples (clipped at the ECG's right edge).
    const int pulseClip = pulseClipN();
    double pLo = 1e300, pHi = -1e300;
    auto merge_pulse = [&](const std::vector<double>& v,
        const std::vector<double>& sd, int visN) {
            const int n = std::min(visN, pulseClip);
            if (n < 1) return;
            double lo, hi; compute_visible_range(v, sd, n, lo, hi);
            pLo = std::min(pLo, lo); pHi = std::max(pHi, hi);
        };
    if (m_hasPPG && !m_ppg.empty()) merge_pulse(m_ppg, m_ppgIqr, m_ppgVisibleN);
    merge_pulse(m_abp, m_abpIqr, static_cast<int>(m_abp.size()));
    merge_pulse(m_art, m_artIqr, static_cast<int>(m_art.size()));
    merge_pulse(m_artPulm, m_artPulmIqr, static_cast<int>(m_artPulm.size()));
    if (pLo > pHi) { pLo = 0.0; pHi = 1.0; }
    if (pLo > 0.0) pLo = 0.0;

    // Each axis keeps its own y-min (autoscaled from the DATA it displays,
    // clipped at 0). The pulse axis already only considers samples up to
    // pulseClipN() = m_ecgVisibleN via the merge_pulse lambda above, so the
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
        auto y_pix = [&](double val, double lo, double hi) {
            const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
            return margin_top + ph - (val - lo) / r * ph;
            };

        // LEFT y-axis: ECG (normalized). Top = yHi, bottom = yLo, plus a
        // major tick at 0.0. Fixed 1-decimal format, no scientific notation.
        // Gray to match the x-axis ticks; labels right-aligned snug to the axis
        // (right edge stops just short of the tick marks so they don't overlap).
        p.setPen(QColor(150, 150, 150));
        {
            const double yTop = y_pix(yHi, yLo, yHi);
            const double yBot = y_pix(yLo, yLo, yHi);
            const double yZero = y_pix(0.0, yLo, yHi);
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
            const double yTop = y_pix(pHi, pLo, pHi);
            const double yBot = y_pix(pLo, pLo, pHi);
            const double yZero = y_pix(0.0, pLo, pHi);
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
        struct art_trace {
            const std::vector<double>* v;
            const std::vector<double>* sd;
            QColor line;
            QColor band;
            bool show;
            double ratio;
        };
        const art_trace arts[] = {
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
                draw_iqr_band(p, v, sd, startPx, margin_top, ph, effPps, visN, pLo, pHi, a.band);
            draw_trace_fixed_scale(p, v, startPx, margin_top, ph, effPps,
                QPen(with_trace_alpha(a.line), 1.3), visN, pLo, pHi);
        }
    }

    // -------- ECG (left axis) --------
    if (m_showEcgTrace) {
        draw_iqr_band(p, m_ecg, m_ecgIqr, margin_left, margin_top, ph, pps, m_ecgVisibleN, yLo, yHi, color_iqrband_ecg);
        draw_trace_fixed_scale(p, m_ecg, margin_left, margin_top, ph, pps,
            QPen(with_trace_alpha(ecg_trace_color), 1.5), m_ecgVisibleN, yLo, yHi);
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
            draw_iqr_band(p, m_ppg, ppgIqrReal, startPx, margin_top, ph, effPps,
                ppgN, pLo, pHi, color_iqrband_ppg);
            draw_trace_fixed_scale(p, m_ppg, startPx, margin_top, ph,
                effPps, QPen(with_trace_alpha(ppg_trace_color), 1.5), ppgN, pLo, pHi);
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
        if (m == PpgT50) continue;      // p50 is a reactive glyph now, not draggable
        const std::vector<double>* vec = nullptr;
        bool isEcg = false, visible = false;
        int visN = 0;
        double ratio = 1.0;
        if (!markerTrace(m, vec, isEcg, visible, visN, ratio)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        if (idx >= visN) continue;
        double mx = xFromSample(idx, isEcg ? 0.0 : ppgStartSample(), ratio);
        QPen pen(marker_color(m), 2);
        pen.setStyle(markerIsBegin(m) ? Qt::DashLine : Qt::SolidLine);
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(QPointF(mx + 2, margin_top + 8), marker_short_label(m));
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
        // B2 focus mode: focus selection is driven by the user BAR (the
        // draggable marker), NOT the automated glyph. A click on a bar selects
        // that landmark for the focus panel (using the bar's own position) and
        // begins a drag; a subsequent drag re-fires focus via markerMoved.
        int m = markerAtX(e->position().x());
        if (m >= 0) {
            m_dragMarker = m;
            emit markerDragStarted(m_binIndex, m_leadIndex, m);
            emit landmarkSelected(m_binIndex, m_leadIndex, m, m_markers[m]);
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
    double ratio = 1.0;
    if (!markerTrace(m_dragMarker, vec, isEcg, visible, visN, ratio)) return;
    if (vec->empty()) return;
    int s = sampleFromX(e->position().x(), isEcg ? 0.0 : ppgStartSample(), ratio);
    // Every marker clamps to its trace's visible window (one rule for
    // ECG, PPG, and all arterial channels).
    s = std::clamp(s, 0, std::max(0, visN - 1));
    m_markers[m_dragMarker] = s;
    emit markerMoved(m_binIndex, m_leadIndex, m_dragMarker, s);
    // No snapshot recapture: the snapshot holds only frozen autodetect
    // columns, which a drag cannot affect. The reactive glyphs (T-peak,
    // T50/T80) are recomputed inside the repaint this update() triggers, so
    // they track the bar live without anything having to be invalidated.
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragMarker = -1;
}

void BinPlotWidget::captureGlyphSnapshot(const TemplateBin& b) {
    m_glyphs = GlyphSnapshot{};
    const int c = m_leadIndex;

    if ((int)m_ecg.size() >= 3) {
        const int N = (int)m_ecg.size();
        // Glyph positions are sub-sample. The *_auto_ch fields are already
        // doubles and are now passed through unrounded -- froz() used to
        // lround them "exactly as the bar seeds do in seed_all()", but a
        // draggable BAR has to land on a sample the operator can grab, whereas
        // a drawn GLYPH does not, and rounding it put the mark up to half a
        // sample off the landmark the writers report.
        auto frozen = [&](double v) {
            return (v >= 0.0 && v <= static_cast<double>(N - 1)) ? v : -1.0;
            };
        auto froz = frozen;
        m_glyphs.ecgPBegin = froz(b.p_begin_auto_ch[c]);
        m_glyphs.ecgPPeak = froz(b.p_peak_auto_ch[c]);
        m_glyphs.ecgQ = froz(b.q_begin_auto_ch[c]);
        m_glyphs.ecgS = froz(b.s_end_auto_ch[c]);
        m_glyphs.ecgTend = froz(b.t_end_auto_ch[c]);
        m_glyphs.ecgRPeak = frozen(m_markers[EcgRPeak]);
    }

    if (m_hasPPG && (int)m_ppg.size() >= 3) {
        const int N = (int)m_ppg.size();
        // The ppg_*_auto source fields are still int in TemplateBin (a
        // versioned on-disk struct), so nothing is gained here yet -- but the
        // lambda is double so the glyph path stops being the thing that
        // narrows once those fields widen.
        auto frozen = [&](double v) {
            return (v >= 0.0 && v <= static_cast<double>(N - 1)) ? v : -1.0;
            };
        m_glyphs.ppgFoot = frozen(b.ppg_onset_auto);
        m_glyphs.ppgP1 = frozen(b.ppg_peak_auto);
        m_glyphs.ppgP2 = frozen(b.ppg_peak2_auto);    m_glyphs.ppgPeak2Found = b.ppg_peak2_found_auto;
        m_glyphs.ppgDic = frozen(b.ppg_dicrotic_auto); m_glyphs.ppgNotchFound = b.ppg_dicrotic_found_auto;
        m_glyphs.ppgEnd = frozen(b.ppg_end_auto);
        m_glyphs.vpgU = frozen(b.ppg_u_auto);
        m_glyphs.vpgV = frozen(b.ppg_v_auto);
        m_glyphs.vpgW = frozen(b.ppg_w_auto);
        m_glyphs.apgA = frozen(b.ppg_a_auto);
        m_glyphs.apgB = frozen(b.ppg_b_auto);
        m_glyphs.apgC = frozen(b.ppg_c_auto);
        m_glyphs.apgD = frozen(b.ppg_d_auto);
        m_glyphs.apgE = frozen(b.ppg_e_auto);
        m_glyphs.apgF = frozen(b.ppg_f_auto);
        m_glyphs.jpgP1 = frozen(b.ppg_p1_auto);
        m_glyphs.jpgP2 = frozen(b.ppg_p2_auto);
    }
}

void BinPlotWidget::drawFeatureGlyphs(QPainter& p,
    double yLo, double yHi, double pLo, double pHi, int ph) const
{
    // Reactive columns, recomputed every paint from the current bars.
    const Reactive rx = reactiveGlyphs();

    auto plot_y = [&](double val, double lo, double hi) {
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
        return margin_top + ph - (val - lo) / r * ph;
        };

    // The detector found a real landmark.
    auto x_glyph = [&](double x, double y) {
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(Qt::black, marker_pen_size));
        p.drawLine(QPointF(x - marker_half_size, y - marker_half_size),
            QPointF(x + marker_half_size, y + marker_half_size));
        p.drawLine(QPointF(x - marker_half_size, y + marker_half_size),
            QPointF(x + marker_half_size, y - marker_half_size));
        };
    // The detector fell back to a placeholder position. Only the notch and
    // the diastolic peak can reach this: they are the two marks that carry a
    // *_found flag.
    auto circle_glyph = [&](double x, double y) {
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(Qt::black, marker_circle_pen));
        p.drawEllipse(QPointF(x, y), marker_circle_radius, marker_circle_radius);
        };
    // A derivative (VPG/APG/JPG) landmark: horizontal dash, coloured by
    // derivative order. No draggable bar exists for these, and the colour
    // plus the shape both say so.
    auto dash_glyph = [&](double x, double y, const QColor& col) {
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(col, marker_pen_size));
        p.drawLine(QPointF(x - dash_half_width, y), QPointF(x + dash_half_width, y));
        };

    // ---- ECG --------------------------------------------------------------
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
        // One index -> plot-point rule for the whole block. False when the
        // landmark is unset (-1) or out of range, so a failed detection
        // draws nothing and no shape carries its own bounds logic.
        // idx is a sub-sample position; the amplitude is interpolated at it
        // rather than read from a rounded column, so the glyph sits on the
        // trace instead of up to half a sample beside it.
        auto point = [&](double idx, QPointF& out) {
            if (idx < 0.0 || idx > static_cast<double>(N - 1)) return false;
            const double raw = FeatureMarks::sample_at(v, idx);
            const double val = std::isnan(raw) ? baseline : raw;
            out = QPointF(xFromSample(idx, 0.0, rateRatio(Channel::Ecg)),
                plot_y(val, yLo, yHi));
            return true;
            };
        auto cross = [&](double idx) { QPointF q; if (point(idx, q)) x_glyph(q.x(), q.y()); };

        cross(m_glyphs.ecgPBegin);   // P begin
        cross(m_glyphs.ecgPPeak);    // P wave
        cross(m_glyphs.ecgQ);        // Q onset
        cross(m_glyphs.ecgRPeak);    // R wave
        cross(m_glyphs.ecgS);        // S end
        cross(rx.ecgTPeak);          // reactive: between the T-begin/T-end bars
        cross(m_glyphs.ecgTend);     // T end
    }

    // ---- PPG --------------------------------------------------------------
    if (m_showPpgTrace && m_hasPPG && (int)m_ppg.size() >= 3) {
        const std::vector<double>& v = m_ppg;
        const int N = (int)v.size();
        // Same rule as the ECG block, but PPG x-geometry (foot-anchored
        // start, PPG rate ratio) and the right-axis scale.
        auto point = [&](double idx, QPointF& out) {
            if (idx < 0.0 || idx > static_cast<double>(N - 1)) return false;
            const double raw = FeatureMarks::sample_at(v, idx);
            const double val = std::isnan(raw) ? pLo : raw;
            out = QPointF(xFromSample(idx, ppgStartSample(), rateRatio(Channel::Ppg)),
                plot_y(val, pLo, pHi));
            return true;
            };
        auto cross = [&](double idx) { QPointF q; if (point(idx, q)) x_glyph(q.x(), q.y()); };
        auto circle = [&](double idx) { QPointF q; if (point(idx, q)) circle_glyph(q.x(), q.y()); };
        auto dash = [&](double idx, const QColor& col) {
            QPointF q; if (point(idx, q)) dash_glyph(q.x(), q.y(), col);
            };
        // X when the shape detection succeeded, O when it fell back.
        auto found = [&](double idx, bool ok) { ok ? cross(idx) : circle(idx); };

        cross(m_glyphs.ppgFoot);
        cross(rx.ppgT50);            // reactive: 50% onset->peak
        cross(m_glyphs.ppgP1);       // systolic peak
        found(m_glyphs.ppgDic, m_glyphs.ppgNotchFound);
        found(m_glyphs.ppgP2, m_glyphs.ppgPeak2Found);
        cross(m_glyphs.ppgEnd);
        cross(rx.ppgT80);            // reactive: 80% peak->end
        // Derivative landmarks: auto-only, no draggable bar, so a horizontal
        // dash rather than an X, coloured by derivative order. Horizontal
        // because the PPG upstroke is near-vertical at u and a, where a
        // vertical dash lies along the trace and vanishes into it.
        //
        // Behind their own checkbox: eleven extra marks on one pulse is noise
        // while the operator is dragging bars.
        if (m_showPpgDerivMarkers) {
            dash(m_glyphs.vpgU, vpg_mark_color);
            dash(m_glyphs.vpgV, vpg_mark_color);
            dash(m_glyphs.vpgW, vpg_mark_color);
            dash(m_glyphs.apgA, apg_mark_color);
            dash(m_glyphs.apgB, apg_mark_color);
            dash(m_glyphs.apgC, apg_mark_color);
            dash(m_glyphs.apgD, apg_mark_color);
            dash(m_glyphs.apgE, apg_mark_color);
            dash(m_glyphs.apgF, apg_mark_color);
            dash(m_glyphs.jpgP1, jpg_mark_color);
            dash(m_glyphs.jpgP2, jpg_mark_color);
        }
    }
}