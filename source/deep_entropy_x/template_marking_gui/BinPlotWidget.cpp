/*
* @brief BinPlotWidget.cpp
*
* Signals are drawn such that the widest trace fills the  cell the layout gives this widget.
* They have the same on screen length
*
* ECG/PPG alignment:
*   Every channel is drawn in SECONDS RELATIVE TO ITS OWN R, so R lands at the
*   same x on every trace by construction. Sample 0 is NOT a shared instant:
*   the ECG's R column is 0.3 * the bin's LONGEST RR (alignment.hpp) while every
*   pulse channel's is a fixed 0.3 s (create_arterial_templates.hpp), and
*   treating them as equal is what put the ECG most of a second ahead of the
*   PPG on any bin holding a pause. See the geometry note in BinPlotWidget.hpp.
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
*   O       a fallback produced the position rather than the fit -- currently
*           the Q onset alone, where the monophasic-R slope walk ran because
*           there was no Q trough to fit against. A different measurement, not
*           a worse version of the same one.
*   dash    an auto-only derivative landmark, coloured by derivative order:
*           VPG blue, APG green, JPG amber.
*
* Marks drawn:
*   ECG  - Bars (user control): P begin, Q onset, S end, T end -- one per
*          ALIGNMENT, see anchor_view.hpp. Glyphs (automatic): P peak, R peak,
*          T peak. A bar's column belongs to the alignment its close-up shows
*          and is reported only there; a glyph is measured on all four.

*   PPG  - User control: Onset, Dicrotic Notch. Automatic: 50% Rise, 80% rise, T80, Foot, Systolic Peak, Diastolic Peak
*   VPG  - u, v, w                      | dashes, behind the
*   APG  - a, b, c, d, e, f             | "Show PPG Derivative Markers"
*   JPG  - p1, p2                       | checkbox
*
*/

#include "BinPlotWidget.hpp"
#include "template_anchoring\anchor_view.hpp"
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

// anchor_view.hpp mirrors these ids as literals because it cannot include this
// header back (BinPlotWidget.hpp -> template_marking_bin_io.hpp ->
// anchor_view.hpp). Drift between the two would route a drag to the wrong
// alignment silently, so it is a build error instead.
static_assert(BinPlotWidget::EcgPBegin == anchor_view::kPBegin, "anchor_view marker id drift");
static_assert(BinPlotWidget::EcgPPeak == anchor_view::kPPeak, "anchor_view marker id drift");
static_assert(BinPlotWidget::EcgQBegin == anchor_view::kQBegin, "anchor_view marker id drift");
static_assert(BinPlotWidget::EcgRPeak == anchor_view::kRPeak, "anchor_view marker id drift");
static_assert(BinPlotWidget::EcgSEnd == anchor_view::kSEnd, "anchor_view marker id drift");
static_assert(BinPlotWidget::EcgTEnd == anchor_view::kTEnd, "anchor_view marker id drift");
// Every marker anchor_view calls a bar must be one markerAtX will actually
// hit-test, or an alignment would own a landmark the operator cannot reach.
static_assert(anchor_view::isBar(BinPlotWidget::EcgPBegin), "");
static_assert(anchor_view::isBar(BinPlotWidget::EcgQBegin), "");
static_assert(anchor_view::isBar(BinPlotWidget::EcgSEnd), "");
static_assert(anchor_view::isBar(BinPlotWidget::EcgTEnd), "");
static_assert(anchor_view::isGlyph(BinPlotWidget::EcgRPeak), "");
static_assert(anchor_view::isGlyph(BinPlotWidget::EcgPPeak), "");

namespace {

    // The slicer's lead-in. Every pulse template is R-anchored with R1 at
    // padSeconds * channelRate (create_arterial_templates.hpp:501), and all
    // four CreatePulseTemplates call sites -- PPG in make_averaged_templates,
    // ABP/ART/ART_PULM in build_templates -- omit the argument and so take its
    // 0.3 default. Keep in step with that default.
    //
    // The ECG's R column is NOT this: alignment.hpp puts it at 0.3 * the bin's
    // longest RR, which varies per bin, so it arrives through setData instead.
    constexpr double kSlicePadSeconds = 0.3;

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
    constexpr QColor ecg_q_begin_color{ 20,  20,  60 };
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
        case BinPlotWidget::EcgQBegin:   return ecg_q_begin_color;
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

void BinPlotWidget::setChannelAnchor(Channel ch, double rColumn) {
    const size_t i = static_cast<size_t>(ch);
    if (m_rAnchor[i] == rColumn) return;
    m_rAnchor[i] = rColumn;
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

    // The stored ECG array is framed on the bin's LONGEST RR, because no beat
    // may lose a sample to framing -- then the outlier filters drop the pause
    // beats that justified the width, leaving an all-NaN tail. lastFinite
    // handles that. What it does not handle is a tail supported by one or two
    // surviving beats: those columns are finite but are not a waveform.
    // align_beat_matrix only writes ecg_template_iqr when nc >= 2, so a
    // one-beat column reads exactly 0.0 -- trim those. Trims the FRAME only,
    // never m_ecg.
    int ecgLast = lastFinite(m_ecg);
    if (m_ecgIqr.size() == m_ecg.size())
        while (ecgLast > 0 && m_ecgIqr[ecgLast] == 0.0) --ecgLast;

    add(m_ecg, Channel::Ecg, ecgLast);
    if (m_hasPPG) add(m_ppg, Channel::Ppg, -1);
    add(m_abp, Channel::Abp, -1);
    add(m_art, Channel::Art, -1);
    add(m_artPulm, Channel::ArtPulm, -1);

    if (!(lo < hi)) { m_tMin = 0.0; m_tMax = 1.0; return; }   // nothing drawable
    m_tMin = lo;
    m_tMax = hi;
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

    recomputeFrame();
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
    recomputeFrame();   // arterial extents are part of the frame
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

    // P peak between the P-onset and Q-onset bars, T peak between S-end and
    // T-end. Both track a drag of any of the four, and both come from the
    // same FeatureMarks call the CSV/bin writers use, so the screen and the
    // files cannot disagree about where a landmark is.
    const FeatureMarks::ReactiveEcg e = FeatureMarks::reactive_ecg(
        m_ecg, m_markers[EcgPBegin], m_markers[EcgQBegin],
        m_markers[EcgSEnd], m_markers[EcgTEnd]);
    r.ecgPPeak = e.p_peak;
    r.ecgTPeak = e.t_peak;

    if (m_hasPPG) {
        const FeatureMarks::ReactivePpg p = FeatureMarks::reactive_ppg(
            m_ppg, m_markers[PpgOnset], m_markers[PpgPeak],
            m_markers[PpgDicrotic], m_markers[PpgEnd]);
        r.ppgT50 = p.t50;
        r.ppgT80 = p.t80;
        r.ppgPeak2 = p.peak2;
    }
    return r;
}

void BinPlotWidget::setBackgroundTraces(
    const std::vector<std::pair<std::vector<double>, QColor>>& traces)
{
    m_bgTraces = traces;
}

void BinPlotWidget::setBankTraces(
    const std::vector<std::pair<std::vector<double>, QColor>>& traces) {
    m_bankTraces = traces;
    update();   // repaint now; otherwise the overlay waits for an unrelated one
}

// Resolve a marker to its channel, trace, and group visibility.
//
// No visible-sample bound and no rate ratio: both existed only because the
// frame used to be the ECG's own sample space, so a pulse marker needed
// converting into it and clipping at its right edge. The frame is now the union
// of every channel's extent in time, so a marker inside its own array is on
// screen, and its x comes from xFromSample(ch, i).
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


int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        const int idx = m_markers[m];
        if (idx < 0) continue;
        // Auto-only marks: drawn as glyphs, never as draggable bars.
        if (m == EcgRPeak || m == EcgPPeak || m == PpgPeak || m == PpgT80
            || m == PpgT50 || m == PpgPeak2) continue;
        const std::vector<double>* vec = nullptr;
        Channel ch = Channel::Ecg;
        bool visible = false;
        if (!markerTrace(m, vec, ch, visible)) continue;
        if (!visible) continue;
        // ARRAY BOUNDS ARE THE ONLY BOUND -- see markerTrace.
        if (idx >= static_cast<int>(vec->size())) continue;
        const double d = std::abs(x - xFromSample(ch, idx));
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
    p.setPen(QColor(150, 150, 150));
    { QFont f = p.font(); f.setPointSize(8); p.setFont(f); }

    QString titleLine = QString("Bin %1  %2").arg(m_binIndex).arg(m_leadLabel);
    // Bin 0 carries the x-axis units hint, since it is the panel whose axis is
    // labelled for the page.
    if (m_binIndex == 0) titleLine += "  (time in seconds)";

    QStringList counts;
    if (m_nEcgBeats > 0) counts << QString("ECG beats %1").arg(m_nEcgBeats);
    if (m_nPpgBeats > 0) counts << QString("PPG beats %1").arg(m_nPpgBeats);

    // Baselines rather than a rect: margin_top is 20 px and two 8 pt lines are
    // ~22, so an AlignBottom rect would push the second line into the plot
    // frame. 9 and 19 keep both clear of it.
    p.drawText(margin_left, 9, titleLine);
    if (!counts.isEmpty())
        p.drawText(margin_left, 19, counts.join("   "));

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

        // ---- Section 4.6 bank overlay ----------------------------------
        // Every template in this bin's bank beyond slot 0, drawn UNDER the
        // sinus trace on the same axis and the same scale. They share slot 0's
        // geometry legitimately: alignment puts every beat's detected R at the
        // same column regardless of morphology, so x0/dx need no adjustment.
        //
        // Dashed and thinner so slot 0 still reads as the primary trace. No
        // label is drawn: an unconfirmed template must not display a class,
        // because showing a guessed one would be the display making the very
        // judgment the operator is being asked to make.
        for (const auto& bg : m_bankTraces) {
            if (static_cast<int>(bg.first.size()) < 3) continue;
            draw_trace_fixed_scale(p, bg.first, x0, margin_top, ph, dx,
                QPen(with_trace_alpha(bg.second), 1.1, Qt::DashLine),
                static_cast<int>(bg.first.size()), yLo, yHi);
        }
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
        int idx = m_markers[m];
        if (idx < 0) continue;
        if (m == EcgRPeak || m == EcgPPeak || m == PpgPeak || m == PpgT80
            || m == PpgT50 || m == PpgPeak2) continue;
        const std::vector<double>* vec = nullptr;
        Channel ch = Channel::Ecg;
        bool visible = false;
        if (!markerTrace(m, vec, ch, visible)) continue;
        if (!visible) continue;
        if (idx >= (int)vec->size()) continue;
        const double mx = xFromSample(ch, idx);
        QPen pen(marker_color(m), 2);
        pen.setStyle(markerIsBegin(m) ? Qt::DashLine : Qt::SolidLine);
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(
            QPointF(mx + 2, margin_top + marker_text_y_offset(m)),
            marker_short_label(m)
        );
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
            emit landmarkSelected(m_binIndex, m_leadIndex, m_templateIndex,
                m, m_markers[m]);
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
            m_state = State::Good;
            emit badPPGToggled(m_binIndex, m_templateIndex, false);
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
    // Clamped to the trace's own array, which is the only bound now -- see
    // markerTrace.
    int s = std::clamp(sampleFromX(ch, e->position().x()),
        0, static_cast<int>(vec->size()) - 1);
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
        // (no ecgPPeak: the P peak is REACTIVE now, bracketed by the P-onset
        //  and Q-onset bars -- see reactiveGlyphs.)
        m_glyphs.ecgQ = froz(b.q_begin_auto_ch[c]);
        m_glyphs.ecgQFound = b.q_begin_found_auto_ch[c];
        m_glyphs.ecgS = froz(b.s_end_auto_ch[c]);
        m_glyphs.ecgQPeak = froz(b.q_peak_auto_ch[c]);
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
        m_glyphs.ppgP2 = frozen(b.ppg_peak2_auto);
        m_glyphs.ppgDic = frozen(b.ppg_dicrotic_auto);
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

// Replace the BIN's pulse glyphs with this bank slot's own. captureGlyphSnapshot
// reads b.ppg_*_auto, which describes b.ppgTemplate -- the wrong waveform for
// every slot the panel actually draws. Called AFTER captureGlyphSnapshot.
void BinPlotWidget::overridePulseGlyphs(const tbank::BankPulseMarkerSet& pm) {
    if ((int)m_ppg.size() < 3) return;
    const int N = (int)m_ppg.size();
    auto froz = [&](double v) {
        return (v >= 0.0 && v <= (double)(N - 1)) ? v : -1.0;
        };
    m_glyphs.ppgFoot = froz(pm.onset_auto);
    m_glyphs.ppgP1 = froz(pm.peak_auto);
    m_glyphs.ppgDic = froz(pm.dicrotic_auto);  m_glyphs.ppgNotchFound = pm.notch_found;
    m_glyphs.ppgP2 = froz(pm.peak2_auto);
    m_glyphs.ppgEnd = froz(pm.end_auto);
    update();
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
            out = QPointF(xFromSample(Channel::Ecg, idx), plot_y(val, yLo, yHi));
            return true;
            };
        auto cross = [&](double idx) { QPointF q; if (point(idx, q)) x_glyph(q.x(), q.y()); };
        auto circle = [&](double idx) { QPointF q; if (point(idx, q)) circle_glyph(q.x(), q.y()); };
        // X when the landmark was fitted, O when a fallback produced the position.
        // Same convention as the PPG block below, which has its own copy because
        // each block's cross/circle close over its own axis and geometry.
        auto found = [&](double idx, bool ok) { ok ? cross(idx) : circle(idx); };

        cross(m_glyphs.ecgPBegin);   // P begin
        cross(rx.ecgPPeak);          // reactive: P-onset bar -> Q-onset bar
        found(m_glyphs.ecgQ, m_glyphs.ecgQFound);
        cross(m_glyphs.ecgQPeak);
        cross(m_glyphs.ecgRPeak);    // R wave
        cross(m_glyphs.ecgS);        // S end
        cross(rx.ecgTPeak);          // reactive: between the S-end/T-end bars
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

        cross(m_glyphs.ppgFoot);
        cross(rx.ppgT50);            // reactive: 50% onset->peak
        cross(m_glyphs.ppgP1);       // systolic peak
        found(m_glyphs.ppgDic, m_glyphs.ppgNotchFound);
        cross(m_glyphs.ppgP2);
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