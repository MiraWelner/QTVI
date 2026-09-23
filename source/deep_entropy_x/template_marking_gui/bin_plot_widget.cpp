/*
* @brief BinPlotWidget.cpp
*
* Handles the drawing of a single plot on the template marking GUI.
*
*/

#include "bin_plot_widget.hpp"
#include "sample_extent.hpp"
#include "template_marking_gui\anchor_view.hpp"
#include "noise_marking_gui/annotation_types.hpp"
#include <QMenu>
#include <QAction>
#include "feature_marks.hpp"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <QStringList>
#include <algorithm>
#include <cmath>
#include <limits>


namespace {

    // The slicer's lead-in. Every pulse template is R-anchored with R1 at
    // padSeconds * channelRate (create_arterial_templates.hpp:501), and all
    // four CreatePulseTemplates call sites -- PPG in make_averaged_templates,
    // ABP/ART/ART_PULM in build_templates -- omit the argument and so take its
    // 0.3 default. Keep in step with that default.
    //
    // The ECG's R column is NOT this: alignment.hpp puts it at 0.4 * the bin's
    // longest RR, which varies per bin, so it arrives through setData instead.
    constexpr double kSlicePadSeconds = 0.4;

    // ------------------------------------------------------------------
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
    constexpr QColor ecg_p_begin_color{ 150, 80, 180 };
    constexpr QColor ecg_p_peak_color{ 180, 100, 210 };
    constexpr QColor ecg_q_onset_color{ 20,  20,  60 };
    constexpr QColor ecg_r_peak_color{ 15,  15,  40 };
    constexpr QColor ecg_s_color{ 30, 35, 85 };
    constexpr QColor ecg_t_end_color{ 70,  90, 160 };

    // PPG bar markers - shades of red, darkest to lightest.
    constexpr QColor ppg_onset_color{ 110,   0,   0 };  // dark red
    constexpr QColor t50_color{ 150,  20,  20 };  // dark red variant
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
    // ONE ALIGNMENT'S OWN BARS. Teal, the colour the R-aligned overlay used,
    // because in a forced view every bar on screen was measured on the SAME
    // waveform -- what distinguishes them is the landmark, which the label
    // says, not the alignment, which is the same for all of them. The
    // per-landmark palette above is for Automatic, where the bars genuinely
    // come from four different alignments.
    // A glyph is a small X, so it gets a small target -- tighter than
    // click_radius_around_marker, which is sized for grabbing a full-height
    // bar. Within this of the X means the glyph; outside it, the bar.
    constexpr double glyph_click_radius = 8.0;

    constexpr QColor align_bar_color{ 0, 140, 140 };

    // The landmark's letter, to follow the alignment's. Matches the old
    // overlay's "Rp"/"Rq"/"Rs"/"Rt".
    const char* align_bar_letter(int m) {
        switch (m) {
        case BinPlotWidget::EcgPBegin: return "p";
        case BinPlotWidget::EcgQBegin: return "q";
        case BinPlotWidget::EcgSEnd:   return "s";
        case BinPlotWidget::EcgTEnd:   return "t";
        }
        return "?";
    }

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
        case BinPlotWidget::EcgPBegin:   return ecg_p_begin_color;
        case BinPlotWidget::EcgPPeak:    return ecg_p_peak_color;
        case BinPlotWidget::EcgQBegin:   return ecg_q_onset_color;
        case BinPlotWidget::EcgRPeak:    return ecg_r_peak_color;
        case BinPlotWidget::EcgSEnd:     return ecg_s_color;
        case BinPlotWidget::EcgTEnd:     return ecg_t_end_color;
        case BinPlotWidget::PpgOnset:    return ppg_onset_color;
        case BinPlotWidget::PpgT50:      return t50_color;
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

    double marker_text_y_offset(int m)
    {
        switch (m) {
        case BinPlotWidget::EcgPBegin:   return 10.0;
        case BinPlotWidget::EcgPPeak:    return 20.0;
        case BinPlotWidget::EcgQBegin:   return 28.0;
        case BinPlotWidget::EcgSEnd:     return 38.0;
        case BinPlotWidget::EcgTEnd:     return 58.0;

        case BinPlotWidget::PpgOnset:    return 8.0;
        case BinPlotWidget::PpgDicrotic: return 18.0;
        case BinPlotWidget::PpgPeak2:    return 28.0;
        case BinPlotWidget::PpgEnd:      return 38.0;

        default:                         return 8.0;
        }
    }

    // Vertical range of a trace, ignoring the std band, with room for the
    // glyphs. No sample-count argument: every finite sample is inside the frame
    // now, so the range is over the whole array.
    void compute_visible_range(const std::vector<double>& v, double& lo, double& hi) {
        lo = 0.0; hi = 1.0;
        bool have = false;
        for (double m : v) {
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
        // visN is a COUNT OF SAMPLES ON THIS TRACE, and the clamp keeps it one.
        // Callers used to pass the frame width, which after the fit-to-extent
        // change could exceed the array -- and the loops below index v[i] up to
        // it. Clamping here rather than at each call site means no future caller
        // can reintroduce that.
        visN = std::min(visN, static_cast<int>(v.size()));
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
        visN = std::min(visN, static_cast<int>(v.size()));   // see draw_iqr_band
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
    // Anchors start unknown for the same reason the markers do: a channel with
    // no anchor is not drawn, rather than drawn at a guessed position.
    std::fill(std::begin(m_rAnchor), std::end(m_rAnchor), -1.0);
}

void BinPlotWidget::setChannelRate(Channel ch, double hz) {
    const size_t i = static_cast<size_t>(ch);
    if (m_rates[i] == hz) return;
    m_rates[i] = hz;
    // A pulse channel's R column is a fixed number of SECONDS into its
    // template, so it follows the rate. Derived here rather than asked of the
    // caller, because the caller would have to know kSlicePadSeconds and the
    // two would drift apart. The ECG's anchor is not derivable this way -- it is
    // 0.3 * the bin's longest RR, which only alignment knows -- so it comes in
    // through setData as rPeakSample.
    if (ch != Channel::Ecg && hz > 0.0)
        m_rAnchor[i] = kSlicePadSeconds * hz;
    recomputeFrame();
    update();
}

void BinPlotWidget::setReferenceLines(
    const std::vector<global_interval_lines::Line>& lines) {
    m_refLines = lines;
    update();
}

// ---------------------------------------------------------------------------
// The time model. Six small functions, and every position in the widget goes
// through them. See the geometry note in the header for why.
// ---------------------------------------------------------------------------

double BinPlotWidget::timeAt(Channel ch, double i) const {
    const size_t k = static_cast<size_t>(ch);
    const double rate = m_rates[k], anchor = m_rAnchor[k];
    if (!(rate > 0.0) || anchor < 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return (i - anchor) / rate;
}

double BinPlotWidget::xFromTime(double t) const {
    const double drawW = std::max(1, width() - margin_left - margin_right);
    const double span = (m_tMax - m_tMin > 1e-9) ? (m_tMax - m_tMin) : 1.0;
    return margin_left + (t - m_tMin) / span * drawW;
}

double BinPlotWidget::pxPerSecond() const {
    const double drawW = std::max(1, width() - margin_left - margin_right);
    const double span = (m_tMax - m_tMin > 1e-9) ? (m_tMax - m_tMin) : 1.0;
    return drawW / span;
}

// x(i) is affine in i, so these two numbers are all a trace draw needs -- which
// is what lets draw_iqr_band and draw_trace_fixed_scale keep their existing
// (startPx + i * px_per_sample) form.
double BinPlotWidget::channelX0(Channel ch) const {
    const double t0 = timeAt(ch, 0.0);
    return std::isnan(t0) ? static_cast<double>(margin_left) : xFromTime(t0);
}

double BinPlotWidget::channelDx(Channel ch) const {
    const double rate = m_rates[static_cast<size_t>(ch)];
    if (!(rate > 0.0)) return 0.0;   // no rate: not drawable, and the callers
    // test dx > 0 rather than duplicating this
    return pxPerSecond() / rate;     // pixels per sample on THIS channel
}

double BinPlotWidget::xFromSample(Channel ch, double i) const {
    return channelX0(ch) + i * channelDx(ch);
}

int BinPlotWidget::sampleFromX(Channel ch, double x) const {
    const double dx = channelDx(ch);
    if (!(std::abs(dx) > 1e-12)) return 0;
    return static_cast<int>(std::round((x - channelX0(ch)) / dx));
}

void BinPlotWidget::recomputeFrame() {
    // FRAME = THE UNION OF EVERY PRESENT CHANNEL'S OWN DRAWN EXTENT, in seconds
    // relative to R. Because it is a union, no channel can have a tail outside
    // it, which is what makes a per-channel clip count unnecessary rather than
    // merely inconvenient.
    auto firstFinite = [](const std::vector<double>& v) {
        for (int i = 0; i < static_cast<int>(v.size()); ++i)
            if (!std::isnan(v[i])) return i;
        return -1;
        };
    auto lastFinite = [](const std::vector<double>& v) {
        for (int i = static_cast<int>(v.size()) - 1; i >= 0; --i)
            if (!std::isnan(v[i])) return i;
        return -1;
        };

    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();

    auto add = [&](const std::vector<double>& v, Channel ch, int lastOverride) {
        if (v.size() < 2) return;
        const int f = firstFinite(v);
        const int l = (lastOverride >= 0) ? lastOverride : lastFinite(v);
        if (f < 0 || l <= f) return;
        const double t0 = timeAt(ch, f), t1 = timeAt(ch, l);
        if (std::isnan(t0) || std::isnan(t1)) return;   // no rate/anchor: not drawn
        lo = std::min(lo, t0);
        hi = std::max(hi, t1);
        };

    add(m_ecg, Channel::Ecg, -1);
    if (m_hasPPG) add(m_ppg, Channel::Ppg, -1);
    add(m_abp, Channel::Abp, -1);
    add(m_art, Channel::Art, -1);
    add(m_artPulm, Channel::ArtPulm, -1);

    // ECG span pinned to a fixed window (the union of the four alignments'
    // extents), so the axis does not rescale when the anchor changes. Folded
    // in AFTER the per-trace adds so it widens, never narrows: the drawn ECG
    // trace is one of the four whose union this is, so it always fits inside.
    if (m_ecgFrameFixed) {
        lo = std::min(lo, m_ecgFrameLo);
        hi = std::max(hi, m_ecgFrameHi);
    }

    if (!(lo < hi)) { m_tMin = 0.0; m_tMax = 1.0; return; }   // nothing drawable
    m_tMin = lo - 0.005;
    m_tMax = hi + 0.005;
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

    // The ECG's R column is the one anchor that cannot be derived from a rate:
    // alignment.hpp puts it at 0.3 * the bin's LONGEST RR, so it varies per bin
    // and only the caller knows it. Every other channel's follows its rate, in
    // setChannelRate. recomputeFrame then does the rest -- there is no view
    // width to compute here any more, because the frame is a time span over the
    // union of the channels rather than a sample count over one of them.
    m_rAnchor[static_cast<size_t>(Channel::Ecg)] = rPeakSample;

    // NEW OCCUPANT: both detections describe the previous one.
    m_detValid = false;
    m_pdetValid = false;
    recomputeFrame();
    updateGeometry();
    update();
}

// Replace ONLY the ECG trace, its band, R column and beat count, leaving the
// PPG and arterial channels and every marker exactly as they are. Used when
// Automatic alignment re-anchors the grid on a bar click: the ECG average
// changes per anchor but nothing else does, so a full setData (which would
// also re-take the unchanged pulse channels) is unnecessary -- and, crucially,
// this leaves the widget object alive, so a re-skin mid-click does not disturb
// an in-progress drag the way rebuilding the panel would.
void BinPlotWidget::setEcgData(const std::vector<double>& ecg,
    const std::vector<double>& ecgIqr,
    double rPeakSample,
    int nEcgBeats)
{
    m_nEcgBeats = nEcgBeats;
    m_ecg = ecg;
    m_ecgIqr = ecgIqr;
    m_rPeakSample = rPeakSample;
    m_rAnchor[static_cast<size_t>(Channel::Ecg)] = rPeakSample;

    m_detValid = false;   // ECG trace changed; the pulse channels did not
    recomputeFrame();
    updateGeometry();
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }

void BinPlotWidget::setEcgFrame(double tMinSec, double tMaxSec) {
    if (!(tMaxSec > tMinSec)) return;   // ignore a degenerate window
    m_ecgFrameFixed = true;
    m_ecgFrameLo = tMinSec;
    m_ecgFrameHi = tMaxSec;
    recomputeFrame();
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
    recomputeFrame();   // arterial extents are part of the frame
    update();
}

void BinPlotWidget::setMarker(Marker m, double idx) {
    // NO CLAMP. EcgPBegin used to be pinned to firstDrawnSample when it
    // arrived before the window -- "a seeded value can land in the noisy
    // pre-ECG lead-in". It landed there because the seeding detected on the
    // BIN-WIDE average while this panel draws slot 0's own per-slot average:
    // different NaN padding, different first finite sample. Both sides use the
    // per-slot array now, so a P onset cannot precede the window, and pinning
    // it would only hide a real disagreement at the edge.
    m_markers[m] = idx;
    update();
}

// Reactive glyphs. No arithmetic lives here -- FeatureMarks owns the formula,
// and the CSV/bin writers call the same functions with the bar set they're
// reporting on, so the screen and the files agree by construction.
BinPlotWidget::Reactive BinPlotWidget::reactiveGlyphs() const {
    // ---- THE REACTIVE HALF IS CACHED ON THE BARS ------------------------
    //
    // "Reactive" means it is a function of the bar positions, NOT that it is
    // cheap: compute_p_peak and compute_t_peak each run a Gaussian-weighted
    // quadratic AND cubic with a BIC choice between them. Two of those per
    // call.
    //
    // This function is called several times per paint -- drawFeatureGlyphs at
    // its top, detectedLandmarks() inside it, the hit test, the focus panel --
    // and Move-Subsequent update()s EVERY column after the dragged one on every
    // mouse-move, so the per-move cost was (columns x calls-per-paint x 2 fits)
    // and the page crawled. The bars change at most once per move, so one
    // evaluation per distinct bar set is all that is ever needed.
    //
    // The key is the four bars plus the detection identity the reactive half
    // brackets against. It is deliberately NOT time- or counter-based: a stale
    // reactive glyph is a glyph drawn somewhere the bars no longer are.
    const ReactiveKey key{
        m_markers[EcgPBegin], m_markers[EcgQBegin],
        m_markers[EcgSEnd],   m_markers[EcgTEnd],
        m_markers[PpgOnset],  m_markers[PpgPeak],
        m_markers[PpgDicrotic], m_markers[PpgEnd],
        m_bin, m_frame, m_templateIndex, m_detValid
    };
    if (m_rxValid && m_rxKey == key) return m_rx;

    Reactive r;

    // P peak between the P-onset and Q-onset bars, T peak between S-end and
    // T-end. Both track a drag of any of the four, and both come from the
    // same FeatureMarks call the CSV/bin writers use, so the screen and the
    // files cannot disagree about where a landmark is.
    // ONE DETECTION CALL, shared with the bar seeding and the focus panel.
    // ecgFiducials assembles the trace, the R column and the fit modes itself
    // from (bin, lead, slot, alignment), so this panel can no longer detect
    // from its own copies and land somewhere else: m_ecg is the bin-wide array
    // on some paths and m_rPeakSample is the DRAWING AXIS origin, R's column on
    // every alignment, which as a detection seed is off by
    // r_col(R) - r_col(anchor).
    //
    // STILL REACTIVE: the bars go in as an argument and are read fresh on every
    // repaint, so P and T peak track a drag exactly as before.
    if (m_bin) {
        // THE DETECTOR RUNS ON TRACE CHANGES, NOT ON PAINTS. This is called
        // once per repaint, and detect_template_landmarks is the most
        // expensive call in the GUI; nothing in it reacts to the bars except
        // the two bracketed peaks below.
        if (!m_detValid || m_detBin != m_bin || m_detFrame != m_frame
            || m_detSlot != m_templateIndex) {
            // ON m_ecg, THE WAVEFORM ON SCREEN. It is this slot's own average
            // (setData/setEcgData put it there via leadsForBinTemplate), it has
            // been through the notch filter, and it is amplitude-scaled -- so
            // detecting on it is the only way the notch reaches the glyphs and
            // the only way they describe the trace they are drawn on. Detecting
            // on slotView's raw array instead left the glyphs measuring an
            // un-notched signal while the operator looked at a filtered one.
            //
            // Seeded with the DISPLAYED alignment's own R column, the same one
            // the R bar uses. detect_template_landmarks refines only +-7
            // samples around the seed, so the flat R-aligned column would make
            // the R glyph miss the peak on a P/Q/J average.
            // THE SAME CALL THAT SEEDS THE BARS. ecgDetect goes through
            // slotView, so it detects on this slot's RAW stored average with
            // that alignment's own r_col -- byte for byte the inputs
            // seedSlotBars hands seed_bank_template. A bar and its glyph are
            // therefore the same number by construction.
            // NOT m_ecg: that is normalized and optionally notched, and the
            // detector is not scale-invariant (Q_MIN_DEPTH is absolute). The
            // trade is that glyphs no longer follow the notch toggle.
            m_det = ecgDetect(*m_bin, m_leadIndex, m_templateIndex, m_frame,
                m_rates[static_cast<size_t>(Channel::Ecg)],
                m_onOffsetFitMode, m_peakFitMode);

            // ecgDetect's peak fields are a SEED, not a placement.
            // placeEcgPeak converts them; the focus panel and the CSV call it
            // too, so all three show the same number.
            //
            // On m_det.tmpl -- the array the landmarks were measured on -- so
            // the fit and the seed see the same samples. NOT m_ecg, which is
            // normalized and optionally notched per the note above.
            // THE CONTEST IS RUN HERE AND NOWHERE ELSE, and the candidates it
            // produced are KEPT -- the focus panel takes them from
            // peakCandidatesFor() instead of fitting again. It used to re-fit
            // seeded from this result, which is one fit iteration more than the
            // glyph got: bestPeakExtremumFit searches only +-peak_halfwidth
            // around its seed, so two rounds land closer to the apex than one.
            // That extra round is why the focus mark sat on the peak and the
            // glyph sat short of it.
            m_peakCands = {};
            if (m_det.valid && m_det.tmpl) {
                const std::vector<double>& t = *m_det.tmpl;
                auto place = [&](EcgPeak which, double& fid) {
                    const auto pc = ecgPeakCandidates(t, which, fid, m_peakFitMode);
                    m_peakCands[static_cast<size_t>(which)] = pc;
                    // A failed fit leaves the seed standing rather than
                    // blanking a landmark the detector did find.
                    if (pc.valid && pc.placement >= 0.0) fid = pc.placement;
                    };
                place(EcgPeak::R, m_det.lm.r_peak);
                place(EcgPeak::Q, m_det.lm.q_peak);
                place(EcgPeak::P, m_det.lm.p_peak);
                place(EcgPeak::S, m_det.s_peak);
                // t_peak excluded: it is the reactive glyph, recomputed from
                // the S-end / T-end bars on every repaint.
            }

            m_detBin = m_bin;
            m_detFrame = m_frame;
            m_detSlot = m_templateIndex;
            m_detValid = true;
        }

        // The bars are still read fresh on every repaint (see above), so the
        // reactive half is two bracketed argmaxes -- which is what this call
        // now costs.
        tbank::BankMarkerSet bars;
        bars.p_begin = m_markers[EcgPBegin];
        bars.q_onset = m_markers[EcgQBegin];
        bars.s_end = m_markers[EcgSEnd];
        bars.t_end = m_markers[EcgTEnd];
        const EcgFiducials fid = ecgFiducialsFrom(m_det,
            m_rates[static_cast<size_t>(Channel::Ecg)], m_peakFitMode, bars);
        r.ecgPPeak = fid.p_peak;
        r.ecgTPeak = fid.t_peak;
    }

    if (m_hasPPG) {
        const FeatureMarks::ReactivePpg p = FeatureMarks::reactive_ppg(
            m_ppg, m_markers[PpgOnset], m_markers[PpgPeak],
            m_markers[PpgDicrotic], m_markers[PpgEnd]);
        r.ppgT50 = p.t50;
        r.ppgT80 = p.t80;
        r.ppgPeak2 = p.peak2;
    }

    // Keyed on the state read above, INCLUDING m_detValid as it is NOW: the
    // detector may have run in this very call, which changes what the reactive
    // half bracketed against.
    m_rxKey = ReactiveKey{
        m_markers[EcgPBegin], m_markers[EcgQBegin],
        m_markers[EcgSEnd],   m_markers[EcgTEnd],
        m_markers[PpgOnset],  m_markers[PpgPeak],
        m_markers[PpgDicrotic], m_markers[PpgEnd],
        m_bin, m_frame, m_templateIndex, m_detValid
    };
    m_rx = r;
    m_rxValid = true;
    return r;
}

// Resolve a marker to its channel, trace, and group visibility.
//
// No visible-sample bound and no rate ratio: both existed only because the
// frame used to be the ECG's own sample space, so a pulse marker needed
// converting into it and clipping at its right edge. The frame is now the union
// of every channel's extent in time, so a marker inside its own array is on
// screen, and its x comes from xFromSample(ch, i).
int BinPlotWidget::lastDrawnSample(Channel ch) const {
    const std::vector<double>* v = nullptr;
    switch (ch) {
    case Channel::Ecg:     v = &m_ecg;      break;
    case Channel::Ppg:     v = &m_ppg;      break;
    case Channel::Abp:     v = &m_abp;      break;
    case Channel::Art:     v = &m_art;      break;
    case Channel::ArtPulm: v = &m_artPulm;  break;
    default: return -1;
    }
    int last = sample_extent::lastFinite(*v);
    // SAME TRIM recomputeFrame applies, and it has to be the same or the wall
    // this reports is not the wall that was drawn.
    if (ch == Channel::Ecg && m_ecgIqr.size() == m_ecg.size())
        last = sample_extent::lastDrawn(m_ecg, m_ecgIqr);
    return last;
}
int BinPlotWidget::firstDrawnSample(Channel ch) const {
    const std::vector<double>* v = nullptr;
    switch (ch) {
    case Channel::Ecg:     v = &m_ecg;      break;
    case Channel::Ppg:     v = &m_ppg;      break;
    case Channel::Abp:     v = &m_abp;      break;
    case Channel::Art:     v = &m_art;      break;
    case Channel::ArtPulm: v = &m_artPulm;  break;
    default: return -1;
    }
    int first = sample_extent::firstFinite(*v);
    // SAME TRIM lastDrawnSample applies, mirrored. align_beat_matrix leaves
    // IQR == 0.0 on columns that had fewer than two beats, and recomputeFrame
    // trims them off both ends -- so a bar dropped there would sit outside the
    // drawn extent even though the sample is not NaN.
    // sample_extent::firstDrawn IS this trim -- one definition, shared with
    // compute_p_begin so the detector and the painter agree on where the
    // drawn extent starts.
    if (ch == Channel::Ecg)
        first = sample_extent::firstDrawn(m_ecg, m_ecgIqr);
    return first;
}

bool BinPlotWidget::markerTrace(int m, const std::vector<double>*& vec,
    Channel& ch, bool& visible) const
{
    if (markerIsEcg(m)) {
        vec = &m_ecg; ch = Channel::Ecg; visible = m_showEcgMarkers;
        return !m_ecg.empty();
    }
    if (markerIsPpg(m)) {
        vec = &m_ppg; ch = Channel::Ppg; visible = m_showPpgMarkers;
        return m_hasPPG && !m_ppg.empty();
    }
    if (markerIsAbp(m)) {
        vec = &m_abp; ch = Channel::Abp; visible = m_showAbpMarkers;
        return !m_abp.empty();
    }
    if (markerIsArt(m)) {
        vec = &m_art; ch = Channel::Art; visible = m_showArtMarkers;
        return !m_art.empty();
    }
    if (markerIsArtPulm(m)) {
        vec = &m_artPulm; ch = Channel::ArtPulm; visible = m_showArtPulmMarkers;
        return !m_artPulm.empty();
    }
    return false;
}


int BinPlotWidget::markerAtX(double x, double* distOut) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    if (distOut) *distOut = std::numeric_limits<double>::infinity();
    for (int m = 0; m < MarkerCount; ++m) {
        const double idx = m_markers[m];
        if (idx < 0.0) continue;
        // Auto-only marks: drawn as glyphs, never as draggable bars.
        if (m == EcgRPeak || m == EcgPPeak || m == PpgPeak || m == PpgT80
            || m == PpgT50 || m == PpgPeak2) continue;
        const std::vector<double>* vec = nullptr;
        Channel ch = Channel::Ecg;
        bool visible = false;
        if (!markerTrace(m, vec, ch, visible)) continue;
        if (!visible) continue;
        // BOUNDED BY THE PLOT WALL, not by the array. A marker in the NaN
        // tail is inside vec but outside the drawn frame, so hit-testing it
        // would hand out a bar the operator cannot see. Compared as a double:
        // casting a fractional index to int would truncate and let the wall
        // column plus a fraction through.
        const int wall = lastDrawnSample(ch);
        if (wall < 0 || idx > static_cast<double>(wall)) continue;
        const double d = std::abs(x - xFromSample(ch, idx));
        if (d < bestDist) { bestDist = d; best = m; }
    }
    if (distOut && best >= 0) *distOut = bestDist;
    return best;
}

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    const int h = height();
    const int w = width();
    const int ph = h - margin_top - margin_bottom;
    p.fillRect(rect(), Qt::white);

    // TITLE, TWO LINES: identity first, counts second. One line did not fit a
    // panel at page width -- "Bin 1  [Ch1 PQRST_A n=365]  401 ECG beats" ran
    // past the frame and clipped mid-number, which is worse than wrapping
    // because a truncated count still looks like a count.
    //
    // The counts are now BOTH labelled and both come from the caller, which
    // resolves them per template (ECG) and per bin (PPG). The old "n=" inside
    // the label was the per-template number and the trailing figure was the
    // per-bin one; nothing on screen said so.
    // TITLE IN BLACK, counts in gray. The identity line is what the operator
    // reads to know which bin and lead they are looking at, so it is not
    // secondary text; the axis labels and the beat counts around it are.
    { QFont f = p.font(); f.setPointSize(8); p.setFont(f); }

    // TWO TONES ON ONE LINE. "Bin 12" is black and the rest -- the channel,
    // the template name, how it was split, the axis hint -- is gray. Pages are
    // packed by column now, so a bin's templates can straddle a page boundary
    // and the bin number is the only thing on screen that says whether the
    // panel beside this one is the same bin. The template name is secondary to
    // that, and was competing with it at equal weight.
    const QString binPart = QString("Bin %1").arg(m_binIndex);
    QString restPart = "  " + m_leadLabel;
    // Bin 0 carries the x-axis units hint, since it is the panel whose axis is
    // labelled for the page.
    if (m_binIndex == 0) restPart += "  (time in seconds)";

    QStringList counts;
    if (m_nEcgBeats > 0) counts << QString("ECG beats %1").arg(m_nEcgBeats);
    if (m_nPpgBeats > 0) counts << QString("PPG beats %1").arg(m_nPpgBeats);

    // Baselines rather than a rect: margin_top is 20 px and two 8 pt lines are
    // ~22, so an AlignBottom rect would push the second line into the plot
    // frame. 9 and 19 keep both clear of it.
    p.setPen(Qt::black);
    p.drawText(margin_left, 9, binPart);
    // Measured with the font already set above, so the gray half starts where
    // the black half actually ended rather than at a guessed offset.
    p.setPen(QColor(150, 150, 150));
    p.drawText(margin_left + p.fontMetrics().horizontalAdvance(binPart), 9,
        restPart);
    if (!counts.isEmpty()) {
        p.setPen(QColor(150, 150, 150));
        p.drawText(margin_left, 19, counts.join("   "));
    }

    // Y-axis rules for the normalized traces:
    //   Y-max is FIXED at 1.0 (one decimal) for every panel so bins share
    //   a common upper reference.
    //   Y-min is autoscaled from the data, but always clipped at <=0.0 so
    //   the 0.0 major tick is always inside the frame.
    // Left axis (ECG) and right axis (pulse) each follow this rule.
    double yLo = 0, yHi = 0;
    compute_visible_range(m_ecg, yLo, yHi);
    if (yLo > 0.0) yLo = 0.0;

    // Right axis range: shared by ALL pulse traces (PPG + arterial), so they
    // sit on one common normalized scale shown on the right. Over the whole of
    // each array -- there is no clip point now, because the frame is the union
    // of the channels rather than the ECG's own extent, so nothing a channel
    // holds falls outside it.
    double pLo = 1e300, pHi = -1e300;
    auto merge_pulse = [&](const std::vector<double>& v) {
        if (v.size() < 2) return;
        double lo, hi; compute_visible_range(v, lo, hi);
        pLo = std::min(pLo, lo); pHi = std::max(pHi, hi);
        };
    if (m_hasPPG) merge_pulse(m_ppg);
    merge_pulse(m_abp);
    merge_pulse(m_art);
    merge_pulse(m_artPulm);
    if (pLo > pHi) { pLo = 0.0; pHi = 1.0; }
    if (pLo > 0.0) pLo = 0.0;

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

        // X ticks: SECONDS RELATIVE TO R, so 0.00 is the R peak on every
        // channel and the labels left of it are negative. The old axis was
        // samples-since-column-zero divided by the ECG rate, which meant one
        // thing on the ECG and something else on the PPG.
        p.setPen(QColor(150, 150, 150));
        for (int t = 0; t <= 4; ++t) {
            const double tsec = m_tMin + (m_tMax - m_tMin) * t / 4.0;
            double x = xFromTime(tsec);
            if (x > yAxisR) x = yAxisR;
            p.drawLine(QPointF(x, xAxisY), QPointF(x, xAxisY + 3));
            const QString lbl = QString::number(tsec, 'f', 2);
            p.drawText(QPointF(x - 3.0 * lbl.size(), xAxisY + 12), lbl);
            // X-axis caption, centered under the tick numbers.
            if (m_binIndex == 0) {
                p.drawText(QRectF(yAxisL, xAxisY + 13.0, yAxisR - yAxisL, 11.0),
                    Qt::AlignHCenter | Qt::AlignTop, "time (s, 0 = R peak)");
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
    // Drawn on the SHARED right-axis range (pLo,pHi) so every pulse tracing
    // sits on one scale. Each channel supplies its own x0/dx, which carry its
    // rate AND its R column -- so a channel running at a different rate than
    // the ECG lands at the correct real time without a ratio applied here.
    {
        struct art_trace {
            const std::vector<double>* v;
            const std::vector<double>* sd;
            QColor line;
            QColor band;
            bool show;
            Channel ch;
        };
        const art_trace arts[] = {
            { &m_abp,     &m_abpIqr,     QColor(0, 115, 45),   QColor(80, 185, 120, 38),  m_showAbpTrace,     Channel::Abp },     // green
            { &m_art,     &m_artIqr,     QColor(140, 75, 185), QColor(180, 130, 215, 38), m_showArtTrace,     Channel::Art },     // purple
            { &m_artPulm, &m_artPulmIqr, QColor(215, 135, 45), QColor(235, 175, 100, 38), m_showArtPulmTrace, Channel::ArtPulm }, // orange
        };
        for (const auto& a : arts) {
            if (!a.show) continue;
            const std::vector<double>& v = *a.v;
            const int n = static_cast<int>(v.size());
            if (n < 2) continue;
            const double x0 = channelX0(a.ch), dx = channelDx(a.ch);
            if (!(dx > 0.0)) continue;          // no rate/anchor: not drawn
            const std::vector<double>& sd = *a.sd;
            if (static_cast<int>(sd.size()) >= n)
                draw_iqr_band(p, v, sd, x0, margin_top, ph, dx, n, pLo, pHi, a.band);
            draw_trace_fixed_scale(p, v, x0, margin_top, ph, dx,
                QPen(with_trace_alpha(a.line), 1.3), n, pLo, pHi);
        }
    }

    // -------- ECG (left axis) --------
    if (m_showEcgTrace) {
        const int n = static_cast<int>(m_ecg.size());
        const double x0 = channelX0(Channel::Ecg), dx = channelDx(Channel::Ecg);
        // BAND AND LINE TAKE THE SAME x0 AND dx, so they cannot drift apart --
        // which they did when the band was reached through one set of
        // correction terms and the line through another.
        draw_iqr_band(p, m_ecg, m_ecgIqr, x0, margin_top, ph, dx, n,
            yLo, yHi, color_iqrband_ecg);

        draw_trace_fixed_scale(p, m_ecg, x0, margin_top, ph, dx,
            QPen(with_trace_alpha(ecg_trace_color), 1.5), n, yLo, yHi);
    }

    // -------- PPG (right/shared axis) --------
    if (m_showPpgTrace && m_hasPPG && m_ppg.size() >= 2) {
        const int n = static_cast<int>(m_ppg.size());
        const double x0 = channelX0(Channel::Ppg), dx = channelDx(Channel::Ppg);
        if (dx > 0.0) {
            draw_iqr_band(p, m_ppg, m_ppgIqr, x0, margin_top, ph, dx, n,
                pLo, pHi, color_iqrband_ppg);
            draw_trace_fixed_scale(p, m_ppg, x0, margin_top, ph, dx,
                QPen(with_trace_alpha(ppg_trace_color), 1.5), n, pLo, pHi);
        }
    }
    p.restore();

    global_interval_lines::paint(p, m_refLines,
        [this](double s) { return xFromSample(Channel::Ecg, s); },
        margin_top, h - margin_bottom, (int)m_ecg.size());

    // Clip marker bars, labels, and fiducial glyphs to the plot area so
    // nothing (including the ~4px X/O glyphs at edge samples) draws past
    // the frame boundary.
    p.save();
    p.setClipRect(QRectF(margin_left, margin_top,
        w - margin_left - margin_right, ph));

    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        double idx = m_markers[m];
        if (idx < 0.0) continue;
        if (m == EcgRPeak || m == EcgPPeak || m == PpgPeak || m == PpgT80
            || m == PpgT50 || m == PpgPeak2) continue;
        const std::vector<double>* vec = nullptr;
        Channel ch = Channel::Ecg;
        bool visible = false;
        if (!markerTrace(m, vec, ch, visible)) continue;
        if (!visible) continue;
        // Same wall the drag clamps to, so a bar is drawn exactly where it
        // can be grabbed.
        const int wallR = lastDrawnSample(ch);
        const int wallL = firstDrawnSample(ch);
        if (wallR < 0 || idx >(double)wallR) continue;
        const double drawIdx = (wallL >= 0 && idx < (double)wallL)
            ? (double)wallL : idx;
        const double mx = xFromSample(ch, drawIdx);
        // OVERLAY STYLE for a forced alignment's own ECG bars. Pulse and
        // arterial bars are foot-anchored and have no alignment, so they keep
        // their palette in every view.
        const bool overlay = !m_alignBadge.isEmpty() && markerIsEcg(m);
        QPen pen(overlay ? align_bar_color : marker_color(m), 2);
        pen.setStyle(overlay ? Qt::DotLine
            : (markerIsBegin(m) ? Qt::DashLine : Qt::SolidLine));
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(
            QPointF(mx + 2, margin_top + marker_text_y_offset(m)),
            overlay
            ? (m_alignBadge + QString::fromLatin1(align_bar_letter(m)))
            : QString::fromLatin1(marker_short_label(m))
        );
    }

    drawFeatureGlyphs(p, yLo, yHi, pLo, pHi, ph);

    p.restore();


    // ONE MARKING PER STATE, so the panel says which verdict it carries rather
    // than leaving it to be inferred from two overlaid marks:
    //   BadR     -> "BAD ECG"
    //   BadPPG   -> "BAD PPG"
    //   BadBoth  -> the diagonal cross
    if (m_state == State::BadBoth) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(margin_left, margin_top, w - margin_right, h - margin_bottom);
        p.drawLine(margin_left, h - margin_bottom, w - margin_right, margin_top);
    }
    else if (m_state == State::BadR || m_state == State::BadPPG) {
        p.setPen(QPen(QColor(200, 0, 0), 2));
        QFont bf = p.font(); bf.setPointSize(14); bf.setBold(true); p.setFont(bf);
        p.drawText(rect(), Qt::AlignCenter,
            (m_state == State::BadR) ? "BAD ECG" : "BAD PPG");
    }
}


void BinPlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        // ---- ONE HIT TEST OVER BARS AND GLYPHS: NEAREST WINS ------------
        //
        // Bars used to win unconditionally inside the 12 px radius, and that
        // made some glyphs unreachable rather than merely hard to hit. The P
        // peak is the worst case: it sits 30-50 ms from the P-onset BAR, and
        // 12 px at page zoom is of the order of 100 ms, so every click aimed
        // at the P peak was claimed by the onset bar -- which then also
        // recorded an operator touch on that bar (user_clicked_on_bar),
        // re-aligned the grid to P_ONSET, and armed a drag on it. When
        // compute_p_begin's `clamp(pb, fFin, pPeak)` saturates the two share a
        // column exactly and no zoom could separate them.
        //
        // Both candidates are now measured and the CLOSER one answers. A bar
        // still wins a tie (strict <), so deliberately grabbing a bar that has
        // a glyph sitting on it behaves as before.
        const double px = e->position().x();

        double barDist = std::numeric_limits<double>::infinity();
        const int mBar = markerAtX(px, &barDist);

        // GLYPHS ARE READ-ONLY, hit-tested at the positions ACTUALLY DRAWN:
        // the same detectedLandmarks() / detectedPulse() / reactiveGlyphs()
        // drawFeatureGlyphs paints from. NOT m_markers[...], which for R holds
        // r_col_raw and sits a few samples off the drawn cross -- that mismatch
        // is why clicking the R glyph used to do nothing.
        int    glyphMarker = -1;
        double glyphIdx = -1.0;
        double glyphDist = click_radius_around_marker;   // must be within radius

        // MEASURED IN 2D, because that is what separates a glyph from the bar
        // sitting on top of it. A bar is a full-height vertical line -- it has
        // no y -- while a glyph is a small X on the trace. Clicking the line
        // anywhere away from the trace means the bar; clicking the X means the
        // glyph.
        //
        // ONE TEST FOR BOTH CHANNELS. It was ECG-only, so no pulse glyph was
        // clickable at all and every PPG click resolved to a bar -- which is
        // half of why the pulse focus panel showed the bar column and the X sat
        // somewhere else.
        const double py = e->position().y();
        const bool haveY = (m_lastPh > 0);
        struct GlyphHit { int marker; double idx; };
        auto testGlyphs = [&](Channel ch, const std::vector<double>& v,
            double axLo, double axHi, const GlyphHit* g, int n) {
                const int wall = lastDrawnSample(ch);
                const double yRange = (axHi - axLo > 1e-10) ? (axHi - axLo) : 1.0;
                for (int k = 0; k < n; ++k) {
                    if (g[k].idx < 0.0) continue;
                    if (wall >= 0 && g[k].idx > static_cast<double>(wall)) continue;
                    const double dx = px - xFromSample(ch, g[k].idx);
                    double d = std::abs(dx);
                    if (haveY) {
                        // Same y the glyph was DRAWN at: the trace value at
                        // that column, on the axis the last paint used. NaN
                        // draws at the axis floor, so test it there too.
                        const double raw = FeatureMarks::sample_at(v, g[k].idx);
                        const double val = std::isnan(raw) ? axLo : raw;
                        const double gy = margin_top + m_lastPh
                            - (val - axLo) / yRange * m_lastPh;
                        const double dy = py - gy;
                        d = std::sqrt(dx * dx + dy * dy);
                    }
                    if (d < glyphDist) {
                        glyphDist = d; glyphMarker = g[k].marker; glyphIdx = g[k].idx;
                    }
                }
            };

        if (m_showEcgTrace && !m_ecg.empty()) {
            const Reactive rx = reactiveGlyphs();
            const FeatureMarks::TemplateLandmarks& lm = detectedLandmarks();
            const GlyphHit glyphs[] = {
                { EcgRPeak, lm.r_peak },
                { EcgPPeak, rx.ecgPPeak },
                { EcgQPeak, lm.q_onset_found ? lm.q_peak : -1.0 },
                { EcgTPeak, rx.ecgTPeak },
                // Transition fiducials at their DRAWN (detected) positions --
                // independent of the bars, so clickable where they're shown even
                // after the bar has been dragged elsewhere.
                { EcgPBegin, lm.p_begin },
                { EcgQBegin, lm.q_onset },
                { EcgSEnd,   lm.s_end },
                { EcgTEnd,   lm.t_end },
            };
            testGlyphs(Channel::Ecg, m_ecg, m_lastYLo, m_lastYHi,
                glyphs, (int)(sizeof(glyphs) / sizeof(glyphs[0])));
        }
        if (m_showPpgTrace && m_hasPPG && !m_ppg.empty()) {
            const Reactive rx = reactiveGlyphs();
            const FeatureMarks::PpgFiducials& pf = detectedPulse();
            const GlyphHit glyphs[] = {
                { PpgOnset,    pf.onset },
                { PpgT50,      rx.ppgT50 },
                { PpgPeak,     pf.peak },
                { PpgDicrotic, pf.dicrotic },
                { PpgPeak2,    pf.peak2 },
                { PpgT80,      rx.ppgT80 },
                { PpgEnd,      pf.end },
            };
            testGlyphs(Channel::Ppg, m_ppg, m_lastPLo, m_lastPHi,
                glyphs, (int)(sizeof(glyphs) / sizeof(glyphs[0])));
        }

        // THE GLYPH WINS ONLY ON A GENUINE 2D HIT.
        //
        // Bars and glyphs are different-shaped targets, so comparing their
        // distances directly never worked: the four transition markers exist as
        // both, and an x-only comparison let the glyph swallow every click on
        // its bar -- which silently turned "select this bar and realign the
        // grid" into "open the close-up", because a glyph click is
        // landmarkFocusOnly and user_clicked_on_bar never ran.
        //
        // In 2D the two separate on their own. The bar is a tall thin line and
        // the glyph is a small X on the trace, so a click needs to be near the
        // glyph in BOTH axes to mean the glyph. Everywhere else on the line is
        // the bar, and a glyph-only mark (R peak, P peak, Q peak, T peak) is
        // reachable because no bar competes for those pixels at all.
        const bool glyphHit = (glyphMarker >= 0)
            && (glyphDist <= glyph_click_radius || mBar < 0)
            && glyphDist < barDist;
        if (glyphHit) {
            emit landmarkFocusOnly(m_binIndex, m_leadIndex,
                m_templateIndex, glyphMarker, glyphIdx);
            return;   // focus only; no drag, no touch
        }
        if (mBar >= 0) {
            m_dragMarker = mBar;
            emit markerDragStarted(m_binIndex, m_leadIndex, mBar);
            emit landmarkSelected(m_binIndex, m_leadIndex, m_templateIndex,
                mBar, m_markers[mBar]);
            return;
        }
    }

    // Ctrl+right-click: confirm this template's CLASS (Section 4.6 bullet 3).
    //
    // Interim gesture, and it is worth saying why rather than leaving it to be
    // discovered. A class confirmation is a deliberate clinical judgment and
    // deserves a visible control -- a combo box or a toolbar -- not a modifier
    // chord. It lives here because the .ui files are generated and were not in
    // reach when this was wired, and because plain right-click is already spoken
    // for by the Good/BadR/BadPPG cycle that operators use constantly. Move it
    // to a real control at the first opportunity; the signal and everything
    // downstream of it stay exactly as they are.
    if (e->button() == Qt::RightButton
        && (e->modifiers() & Qt::ControlModifier)) {
        QMenu menu(this);
        menu.addAction(tr("Confirm class for this template"))->setEnabled(false);
        menu.addSeparator();
        for (const auto& t : annotation_types::noise_types) {
            if (t.paramEdit || t.invertEdit) continue;   // not classes
            QAction* a = menu.addAction(QString::fromUtf8(t.label));
            a->setData(t.code);
        }
        menu.addSeparator();
        QAction* normal = menu.addAction(tr("Normal / sinus (confirm as regular)"));
        normal->setData(0);   // kUnlabeled means "confirmed, no abnormal class"

        if (QAction* picked = menu.exec(e->globalPosition().toPoint())) {
            emit classConfirmRequested(m_binIndex, m_leadIndex,
                m_templateIndex, picked->data().toInt());
        }
        return;
    }

    if (e->button() == Qt::RightButton) {
        // Good -> BadR -> BadPPG (pulse only) -> BadBoth -> Good.
        switch (m_state) {
        case State::Good:
            m_state = State::BadR;
            emit badRToggled(m_binIndex, m_leadIndex, m_templateIndex, true);
            break;
        case State::BadR:
            if (m_hasPPG) {
                m_state = State::BadPPG;
                emit badRToggled(m_binIndex, m_leadIndex, m_templateIndex, false);
                emit badPPGToggled(m_binIndex, m_templateIndex, true);
            }
            else {
                m_state = State::Good;
                emit badRToggled(m_binIndex, m_leadIndex, m_templateIndex, false);
            }
            break;
        case State::BadPPG:
            m_state = State::BadBoth;
            emit badPPGToggled(m_binIndex, m_templateIndex, true);
            emit badRToggled(m_binIndex, m_leadIndex, m_templateIndex, true);
            break;
        case State::BadBoth:
            m_state = State::Good;
            emit badPPGToggled(m_binIndex, m_templateIndex, false);
            emit badRToggled(m_binIndex, m_leadIndex, m_templateIndex, false);
            break;
        }
        update();
        return;
    }
}

void BinPlotWidget::mouseMoveEvent(QMouseEvent* e) {
    if (m_dragMarker < 0) return;
    const std::vector<double>* vec = nullptr;
    Channel ch = Channel::Ecg;
    bool visible = false;
    if (!markerTrace(m_dragMarker, vec, ch, visible)) return;
    if (vec->empty()) return;
    // CLAMPED TO THE PLOT WALL, not to the array. vec->size()-1 let a drag
    // push a bar into the channel's NaN tail, past m_tMax, and the bar drew
    // outside the frame with no trace under it. lastDrawnSample is the same
    // bound recomputeFrame used to build that frame.

    const int wallL = firstDrawnSample(ch);
    const int wallR = lastDrawnSample(ch);
    if (wallL < 0 || wallR < wallL) return;
    int s = std::clamp(sampleFromX(ch, e->position().x()), wallL, wallR);
    // SAME COLUMN, NOTHING TO DO. Several pixels map to one sample at normal
    // zoom, so most drag events re-placed a bar where it already was -- each
    // one running a full propagation pass. A drag only ever produces whole
    // columns (sampleFromX), so the rounded compare loses nothing.
    if (s == static_cast<int>(std::lround(m_markers[m_dragMarker]))) return;
    m_markers[m_dragMarker] = s;
    // TEMPLATE-AWARE signal only. markerMoved carried no slot, so a drag on a
    // sub-template column was indistinguishable from one on slot 0 and wrote
    // into the bin's marker set either way -- which is why sub-templates could
    // not safely have bars at all. The receiver routes on templateIdx.
    emit markerMovedOnTemplate(m_binIndex, m_leadIndex, m_templateIndex,
        m_dragMarker, s);
    // No snapshot recapture: the snapshot holds only frozen autodetect
    // columns, which a drag cannot affect. The reactive glyphs (T-peak,
    // T50/T80) are recomputed inside the repaint this update() triggers, so
    // they track the bar live without anything having to be invalidated.
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    // NOTHING IS EMITTED HERE, deliberately. A release used to fire a
    // markerDragFinished so the owner could re-apply the whole page once per
    // gesture -- but m_dragMarker is armed by any bar CLICK, not just a drag,
    // and an automatic alignment shift IS a bar click, so every one of those
    // paid for a full-page re-detect on mouse-up. The pass had no work to do
    // either: everything that follows a bar is reactive (P and T peak through
    // reactiveGlyphs, T50/T80 likewise) or was already pushed by the
    // propagation loops' setMarker calls.
    m_dragMarker = -1;
}

// ---- THE PULSE DETECTION, ONCE PER (BIN, SLOT) -------------------------
//
// The twin of reactiveGlyphs' m_det block, and it exists for the same reason:
// one run of the detector whose answer is both PAINTED and READ. captureGlyph-
// Snapshot and overridePulseGlyphs used to sit here -- a frozen copy of the ECG
// detection plus a stored BankPulseMarkerSet pushed in from the viewer -- and
// between them they gave every landmark two positions on two invalidation
// schedules.
//
// ON THE SLOT'S OWN STORED AVERAGE, through the bank, NOT on m_ppg: m_ppg is
// the perfusion-normalised display copy, and the pulse finders are no more
// scale-invariant than the ECG ones. Same array seed_pulse_bank_template runs
// on, so the viewer's stored *_auto and this agree by construction. Columns are
// columns -- the axis is shared -- so the positions drop straight onto m_ppg.
const FeatureMarks::PpgFiducials& BinPlotWidget::detectedPulse() const {
    const double rate = m_rates[static_cast<size_t>(Channel::Ppg)];
    if (m_pdetValid && m_pdetBin == m_bin && m_pdetSlot == m_templateIndex)
        return m_pdet;

    m_pdet = FeatureMarks::PpgFiducials{};
    m_pdetBin = m_bin;
    m_pdetSlot = m_templateIndex;
    m_pdetValid = true;

    if (!m_bin || rate <= 0.0 || m_templateIndex < 0) return m_pdet;
    if (m_templateIndex >= m_bin->ppg_bank.size()) return m_pdet;
    const std::vector<double>& t =
        m_bin->ppg_bank.templates[m_templateIndex].tmpl;
    if (t.size() < 3) return m_pdet;

    m_pdet = FeatureMarks::detect_ppg_fiducials(t, static_cast<int>(t.size()),
        rate);
    return m_pdet;
}

void BinPlotWidget::drawFeatureGlyphs(QPainter& p,
    double yLo, double yHi, double pLo, double pHi, int ph) const
{
    // Reactive columns, recomputed every paint from the current bars.
    const Reactive rx = reactiveGlyphs();

    // Remembered for the hit test: see m_lastYLo. BOTH axes, because the hit
    // test covers the pulse glyphs too and they sit on the right-hand scale.
    m_lastYLo = yLo; m_lastYHi = yHi;
    m_lastPLo = pLo; m_lastPHi = pHi;
    m_lastPh = ph;

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
            out = QPointF(xFromSample(Channel::Ecg, idx), plot_y(val, yLo, yHi));
            return true;
            };
        auto cross = [&](double idx) { QPointF q; if (point(idx, q)) x_glyph(q.x(), q.y()); };
        auto circle = [&](double idx) { QPointF q; if (point(idx, q)) circle_glyph(q.x(), q.y()); };
        // X when the landmark was fitted, O when a fallback produced the position.
        // Same convention as the PPG block below, which has its own copy because
        // each block's cross/circle close over its own axis and geometry.
        auto found = [&](double idx, bool ok) { ok ? cross(idx) : circle(idx); };

        // STRAIGHT FROM THE DETECTION every reader of this panel gets.
        const FeatureMarks::TemplateLandmarks& lm = detectedLandmarks();
        cross(lm.p_begin);           // detected, like every other transition
        cross(rx.ecgPPeak);          // reactive: P-onset bar -> Q-onset bar
        found(lm.q_onset, lm.q_onset_found);
        // No Q trough means no Q peak: q_peak comes back -1 there, and the
        // circle above already says the onset was a fallback.
        if (lm.q_onset_found) cross(lm.q_peak);
        cross(lm.r_peak);            // R wave
        cross(lm.s_end);             // S end
        cross(rx.ecgTPeak);          // reactive: between the S-end/T-end bars
        cross(lm.t_end);             // T end
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
            out = QPointF(xFromSample(Channel::Ppg, idx), plot_y(val, pLo, pHi));
            return true;
            };
        auto cross = [&](double idx) { QPointF q; if (point(idx, q)) x_glyph(q.x(), q.y()); };
        auto circle = [&](double idx) { QPointF q; if (point(idx, q)) circle_glyph(q.x(), q.y()); };
        auto dash = [&](double idx, const QColor& col) {
            QPointF q; if (point(idx, q)) dash_glyph(q.x(), q.y(), col);
            };
        // X when the shape detection succeeded, O when it fell back.
        auto found = [&](double idx, bool ok) { ok ? cross(idx) : circle(idx); };

        // SAME TERMS AS THE ECG BLOCK ABOVE: this slot's own detection, the
        // one detectedPulse() hands the focus panel. These used to be
        // b.ppg_*_auto -- the BIN's marks, measured on b.ppgTemplate, which is
        // not the waveform any column draws -- pushed in by overridePulseGlyphs.
        const FeatureMarks::PpgFiducials& pf = detectedPulse();
        cross(pf.onset);
        cross(rx.ppgT50);            // reactive: 50% onset->peak
        cross(pf.peak);              // systolic peak
        found(pf.dicrotic, pf.notch_found);
        cross(pf.peak2);
        cross(pf.end);
        cross(rx.ppgT80);            // reactive: 80% peak->end
        // Derivative landmarks: auto-only, no draggable bar, so a horizontal
        // dash rather than an X, coloured by derivative order. Horizontal
        // because the PPG upstroke is near-vertical at u and a, where a
        // vertical dash lies along the trace and vanishes into it.
        //
        // Behind their own checkbox: eleven extra marks on one pulse is noise
        // while the operator is dragging bars.
        if (m_showPpgDerivMarkers) {
            dash(pf.u, vpg_mark_color);
            dash(pf.v, vpg_mark_color);
            dash(pf.w, vpg_mark_color);
            dash(pf.a, apg_mark_color);
            dash(pf.b, apg_mark_color);
            dash(pf.c, apg_mark_color);
            dash(pf.d, apg_mark_color);
            dash(pf.e, apg_mark_color);
            dash(pf.f, apg_mark_color);
            dash(pf.p1, jpg_mark_color);
            dash(pf.p2, jpg_mark_color);
        }
    }
}