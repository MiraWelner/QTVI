// ============================================================================
// BinPlotWidget.hpp - One ECG lead + optional PPG overlay + draggable markers
//
// Markers:
//   ECG (per channel):  P-onset, Q-begin, T-begin, T-end
//   PPG (shared/bin):   Onset, Peak, Dicrotic notch, 50% point, End
//
// Right-click cycles:  Good -> BadR -> BadPPG -> Good  (skips BadPPG if no PPG)
// Left-drag:           move whichever marker is closest to the click
//
// Drawing scale:
//   The widget fills whatever width its layout cell gives it. The
//   pixels-per-sample scale is computed at draw time (pxPerSample()) so
//   that the widest trace exactly spans the drawable width -- i.e. every
//   bin window is the SAME on-screen length regardless of how many
//   samples it holds. ECG and PPG still share that single scale within a
//   widget, so they stay temporally aligned with each other; only the
//   scale differs from one bin to the next (which is fine -- bins don't
//   need to share an x-axis, just a length).
//
// Per-sample std (gray band):
//   When std vectors matching the trace length are provided, the widget
//   paints a translucent gray polygon between mean-std and mean+std
//   underneath each line. Empty std vectors just disable the band for
//   that trace.
// ============================================================================
#pragma once
#include <QWidget>
#include <QString>
#include <QColor>
#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>
#include "template_marking_bin_io.hpp"
#include "global_interval_lines.hpp"

class QPainter;

class BinPlotWidget : public QWidget {
    Q_OBJECT
public:
    enum class State { Good, BadR, BadPPG };

    // Each enum value MUST be unique (it's used as an array index into
    // m_markers). ECG markers come first, then PPG markers, so
    // markerIsEcg / markerIsPpg can use range checks.
    enum Marker : int {
        EcgPBegin = 0,
        EcgPPeak = 1,
        EcgQBegin = 2,
        EcgRPeak = 3,
        EcgSEnd = 4,
        EcgTBegin = 5,
        EcgTEnd = 6,
        // --- PPG markers (contiguous, immediately after ECG) ---
        PpgOnset = 7,
        PpgT50 = 8,
        PpgPeak = 9,
        PpgDicrotic = 10,
        PpgPeak2 = 11,
        PpgT80 = 12,
        PpgEnd = 13,
        // --- Arterial markers ---
        AbpOnset = 14, AbpPeak = 15, AbpDicrotic = 16, AbpPeak2 = 17, AbpEnd = 18,
        ArtOnset = 19, ArtPeak = 20, ArtDicrotic = 21, ArtPeak2 = 22, ArtEnd = 23,
        ArtPulmOnset = 24, ArtPulmPeak = 25, ArtPulmDicrotic = 26,
        ArtPulmPeak2 = 27, ArtPulmEnd = 28,
        MarkerCount = 29
    };

    static bool markerIsEcg(int m) { return m >= EcgPBegin && m <= EcgTEnd; }
    static bool markerIsPpg(int m) { return m >= PpgOnset && m <= PpgEnd; }
    static bool markerIsAbp(int m) { return m >= AbpOnset && m <= AbpEnd; }
    static bool markerIsArt(int m) { return m >= ArtOnset && m <= ArtEnd; }
    void setReferenceLines(const std::vector<global_interval_lines::Line>& lines);
    static bool markerIsArtPulm(int m) { return m >= ArtPulmOnset && m <= ArtPulmEnd; }
    // Any arterial marker (ABP/ART/ART_PULM) rides the PPG x-geometry.
    static bool markerIsArterial(int m) { return m >= AbpOnset && m <= ArtPulmEnd; }

    // "Begin"/onset markers are drawn dashed; everything else (incl. "end")
    // solid.
    static bool markerIsBegin(int m) {
        return m == EcgPBegin || m == EcgQBegin || m == EcgTBegin
            || m == PpgOnset || m == AbpOnset || m == ArtOnset || m == ArtPulmOnset;
    }
    double m_rPeakSample = 0.0;   // R-peak sample index within the ECG template
    // (m_ppgDelay / m_ppgFootIdx retired in Patch C: every channel is
    // real-time-aligned by construction under Patch B slicing.)

    // Every trace this widget can draw, including ECG (the reference
    // channel that defines the frame's sample-space -- see
    // totalSampleSpan()/pxPerSample()). Add new channels here only, before
    // Count; nothing else in the geometry code needs to change.
    enum class Channel { Ecg, Ppg, Abp, Art, ArtPulm, Count };

    // ----------------------------------------------------------------------
    // Visible-range rules.
    //
    // PPG: always show the whole trace.
    // ECG: cut off before the next beat's P-wave. The cutoff is derived
    //      from the signal itself so it scales with sampling rate and
    //      heart rate:
    //
    //   1. Find this-beat R as the largest deviation from the front-half
    //      mean within samples [0, N/2).
    //   2. Find next-beat R the same way within [N/2, N).
    //   3. RR = next_R - this_R.
    //   4. Cut at next_R - kPMarginRRFrac * RR.
    //
    // At normal heart rates the P-wave starts ~15-20% of RR before QRS;
    // backing off 25% lands the cut firmly in the flat TP-segment before
    // any P-bump. Because the margin is a fraction of RR (which scales
    // linearly with sample rate), this is sample-rate independent.
    // ----------------------------------------------------------------------
    static constexpr double kPMarginRRFrac = 0.25;

    // Pixels-per-sample is no longer fixed. It's derived per-paint from
    // the widget's current width and the widest trace (see pxPerSample()),
    // so each bin window comes out the same on-screen length.

    static int visiblePpgCount(int nFull) {
        return std::max(nFull, 2);
    }

    void setTemplateIndex(int t) { m_templateIndex = t; }
    int  templateIndex() const { return m_templateIndex; }

    explicit BinPlotWidget(int binIndex, int leadIndex,
        const QString& leadLabel, QWidget* parent = nullptr);

    // Traces only. Markers and autodetect positions go in through
    // setMarker()/setAuto(), which TemplateViewerWindow::applyBinToWidget()
    // calls as one unit -- there is deliberately no second way to get a
    // position into this widget, since two paths is how the bars and the
    // glyphs drifted apart in the first place. Std vectors shorter than the
    // visible sample count are ignored at draw time (the band silently
    // disappears for that trace), so stale data is safe.
    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ppgIqr,
        const std::vector<double>& ecg,
        const std::vector<double>& ecgIqr,
        double rPeakSample,
        int nEcgBeats = 0,
        int nPpgBeats = 0);

    void setHasPPG(bool has);
    bool hasPPG() const { return m_hasPPG; }

    void setChannelRate(Channel ch, double hz);
    double channelRate(Channel ch) const { return m_rates[static_cast<size_t>(ch)]; }
    // Ratio that converts `ch`'s own sample index into ECG-equivalent
    // sample units (the frame's reference space). 1.0 for ECG itself, and
    // 1.0 whenever either rate is unknown -- the historical case, when
    // every rate happened to match, so this never changes old behavior for
    // a channel with no rate set.
    double rateRatio(Channel ch) const {
        const double r = m_rates[static_cast<size_t>(ch)];
        const double ecgR = m_rates[static_cast<size_t>(Channel::Ecg)];
        return (r > 0.0 && ecgR > 0.0) ? ecgR / r : 1.0;
    }

    void setState(State s);
    State state() const { return m_state; }

    void setMarker(Marker m, int idx);
    int  marker(Marker m) const { return m_markers[m]; }

    // Recaptures the frozen glyph snapshot and repaints. Call this LAST in a
    // seeding pass (see applyBinToWidget): m_glyphs.ecgRPeak reads the R bar,
    // which must already be set.
    void setAuto(const TemplateBin& b) { captureGlyphSnapshot(b); update(); }

    // Reactive glyphs: pure functions of the CURRENT bar positions, computed
    // on demand at paint time and never stored. T-peak therefore tracks
    // T-begin/T-end live, per drag pixel, and no caller has to remember to
    // refresh anything. The formula itself lives in FeatureMarks, shared with
    // the CSV/bin writers.
    // Sub-sample positions, straight from FeatureMarks. Were int, which
    // rounded the refined T-peak and the interpolated T50/T80 crossings on the
    // way into the paint path -- so the glyph drew up to half a sample away
    // from the value the CSV reported for the same landmark.
    struct Reactive {
        double ecgTPeak = -1.0, ppgT50 = -1.0, ppgT80 = -1.0,
            ppgPeak2 = -1.0;
    };    Reactive reactiveGlyphs() const;

    // Per-trace marker visibility. When false, that group's markers
    // are neither drawn nor hit-testable (drag-pick ignores them).
    // Both default to true.
    void setShowEcgMarkers(bool show);
    void setShowPpgMarkers(bool show);
    void setShowPpgDerivMarkers(bool show);
    void setShowAbpMarkers(bool show);
    void setShowArtMarkers(bool show);
    void setShowArtPulmMarkers(bool show);
    // Per-trace waveform visibility (independent of the markers on that
    // trace). When false, the trace line/band is not drawn, but its
    // markers may still show if their marker-visibility flag is on. All
    // default to true.
    void setShowEcgTrace(bool show);
    void setShowPpgTrace(bool show);
    void setShowAbpTrace(bool show);
    void setShowArtTrace(bool show);
    void setShowArtPulmTrace(bool show);
    // Provide the arterial trace vectors (for marker bounds/geometry). Any
    // may be empty when that channel is absent. Markers on an empty trace
    // are never drawn or hit-tested. Call before setting arterial markers.
    void setArterialTraces(const std::vector<double>& abp,
        const std::vector<double>& art,
        const std::vector<double>& artPulm,
        const std::vector<double>& abpIqr = {},
        const std::vector<double>& artIqr = {},
        const std::vector<double>& artPulmIqr = {});
    void setBackgroundTraces(const std::vector<std::pair<std::vector<double>, QColor>>& traces);
    std::vector<std::pair<std::vector<double>, QColor>> m_bgTraces;

    // Section 4.6 bank overlay: templates 1..N-1 of this (bin, channel)'s bank,
    // drawn on the ECG axis under the slot 0 trace.
    //
    // Deliberately NOT reusing m_bgTraces. That one carries the arterial
    // context traces (ABP / ART / ART_PULM), which are foot-anchored and belong
    // on a different axis -- and it is set from the plot-construction path,
    // which runs BEFORE applyBinToWidget(), so sharing the member would mean
    // whichever wrote last silently erased the other.
    void setBankTraces(const std::vector<std::pair<std::vector<double>, QColor>>& traces);
    std::vector<std::pair<std::vector<double>, QColor>> m_bankTraces;

    int  binIndex()  const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

    int  ecgVisibleN() const { return m_ecgVisibleN; }

    // Median RR in ECG samples. Call BEFORE setData: setData computes the drawn
    // width from it. -1 restores the old behaviour (axis = array length).
    void setMedianRr(int samples) { m_medianRrSamples = samples; }

    // The widget no longer dictates its width from the trace length.
    // It advertises a modest preferred width and a small minimum so the
    // grid can hand every cell an equal share of the window; the trace is
    // then scaled to fill whatever width the cell receives.
    QSize sizeHint() const override { return QSize(220, 120); }
    QSize minimumSizeHint() const override { return QSize(40, 60); }

signals:
    void markerMoved(int binIndex, int leadIndex, int marker, int newIdx);
    void markerDragStarted(int binIndex, int leadIndex, int marker);

    // B2 focus mode: emitted when the operator selects (clicks) a landmark,
    // so the owner can render that landmark's focus panel. Distinct from
    // markerDragStarted (which is about beginning a drag) -- selection fires
    // on the same click but carries the intent "show this landmark's focus
    // view," and the owner decides what to do (e.g. J-point refreshes both
    // the QRS and JT panels).
    // TEMPLATE INDEX INCLUDED, like the marker and quality signals. A panel
    // is a (bin, template) pair, and the focus view reads mean/sd/n for the
    // waveform under the landmark -- without the slot it read the BIN's
    // waveform whichever panel was clicked.
    void landmarkSelected(int binIndex, int leadIndex, int templateIdx,
        int marker, int col);

    // TEMPLATE INDEX ADDED to the marker signals. A panel is a (bin, template)
    // pair now, and without the slot the receiver cannot tell whether a drag
    // belongs to the bin's own marker set or to a bank template's -- so every
    // drag on a sub-template column would have written into slot 0's landmarks.
    void markerMovedOnTemplate(int binIndex, int leadIndex, int templateIdx,
        int marker, int newIdx);

    // TEMPLATE INDEX ADDED, for the same reason the marker signals carry it. A
    // panel is a (bin, template) pair, and without the slot the receiver could
    // only record the verdict against the BIN -- so one right-click on one
    // column crossed out every column of that bin, including morphologies the
    // operator had not looked at.
    void badRToggled(int binIndex, int leadIndex, int templateIdx, bool bad);
    void badPPGToggled(int binIndex, int templateIdx, bool bad);

    // Section 4.6 class confirmation. The operator picks ONE class from the
    // annotation_types table for the template shown in this panel; the owner
    // turns that into tbank::propagateLabel(), which attaches the label to the
    // template and from there to every beat assigned to it.
    //
    // Carries a class code only. No subtype: "the operator never types a
    // subtype index and never sees one until the bank produces it", so the
    // widget has no business knowing one exists.
    //
    // Ctrl+right-click opens the menu, because plain right-click already cycles
    // Good/BadR/BadPPG and that gesture is in operators' hands already. A
    // toolbar combo is the better home for this once the .ui file can be
    // edited -- see the note in BinPlotWidget.cpp.
    void classConfirmRequested(int binIndex, int leadIndex, int templateIdx,
        int annotationCode);
    void overridePulseGlyphs(const tbank::BankPulseMarkerSet& pm);


protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    std::vector<global_interval_lines::Line> m_refLines; //global refernce lines eg. earliest Q-onset 
    int    sampleFromX(double x, double startSample, double ratio) const;
    // s is a sub-sample position: glyph landmarks are fractional.
    double xFromSample(double s, double startSample, double ratio) const;
    // The ONLY hit-test. Bars are clickable; glyphs are display-only and are
    // deliberately not hit-tested, so a click can never select or drag an
    // automated mark.
    int    markerAtX(double x) const;
    int    visibleN(bool isEcg) const;
    // Resolve a marker's trace vector, geometry, current visibility,
    // visible-sample bound, and ECG-equivalent rate ratio. ECG, PPG, and
    // all arterial channels resolve through one path so they behave
    // identically.
    bool   markerTrace(int m, const std::vector<double>*& vec,
        bool& isEcg, bool& visible, int& visN, double& ratio) const;

    // Widest sample extent any trace needs (in samples), and the
    // pixels-per-sample that makes that extent fill the drawable width.
    int    totalSampleSpan() const;
    double pxPerSample() const;

    // View width in ECG samples, derived from the trace: the next beat's R in a
    // median template sits at the median RR, which is what the frame should be
    // sized from. The STORED array is framed on the longest RR so that no beat
    // is clipped, so without this the axis of every panel is stretched by the
    // single worst interval in the bin. Falls back to the full length when no
    // second R is identifiable.

    double ppgStartSample() const;
    // Pulse-trace samples that fit before the ECG's right edge (clip point).
    int    pulseClipN() const;

    // Feature-glyph QC marks (the black X's). Defined in BinPlotGlyphs.cpp.
    // Reads trace/marker state directly (member); takes the paint-local
    // axis ranges + plot height the glyphs need.
    void   drawFeatureGlyphs(QPainter& p,
        double yLo, double yHi, double pLo, double pHi, int ph) const;

    int m_binIndex;
    int m_leadIndex;
    // Which bank slot this panel is showing. Needed because a class
    // confirmation names a TEMPLATE, and one bin now occupies several panels.
    // Defaults to 0 so a pre-bank file behaves as it always did.
    int m_templateIndex = 0;
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ppgIqr;
    std::vector<double> m_ecg;
    std::vector<double> m_ecgIqr;

    int m_ecgVisibleN = 0;
    int m_ppgVisibleN = 0;

    // Median RR of the beats behind this template, in samples; -1 = unknown.
    // Sets the DRAWN width only -- see the note in setData. Never used to trim
    // m_ecg itself.
    int m_medianRrSamples = -1;
    int m_markers[MarkerCount];   // all -1 until seeded (filled in the ctor)

    // Hz per channel (indexed by Channel); 0 = unknown -> rateRatio()
    // falls back to 1.0 for that channel.
    std::array<double, static_cast<size_t>(Channel::Count)> m_rates{};

    // FROZEN glyph columns only -- every field here is a copy of an m_auto
    // value (bounds-checked against the trace), so nothing in this struct can
    // drift as the user drags. The reactive glyphs (ECG T-peak, PPG T50/T80)
    // are deliberately NOT here: they are recomputed on demand by
    // reactiveGlyphs(), because a stored copy of a derived value is exactly
    // what has to be manually refreshed and therefore exactly what goes
    // stale. The one bar-dependent field is ecgRPeak, and R is not draggable.
    // All sub-sample. Every one of these comes from a refined finder, so an
    // int field here re-quantised what the refinement had resolved.
    struct GlyphSnapshot {
        double ecgPBegin = -1.0, ecgPPeak = -1.0, ecgQ = -1.0, ecgRPeak = -1.0,
            ecgS = -1.0, ecgTend = -1.0;
        double ppgFoot = -1.0;    // = ppgOnset auto
        double ppgP1 = -1.0;      // = ppgPeak auto
        double ppgP2 = -1.0;      bool ppgPeak2Found = false;
        double ppgDic = -1.0;     bool ppgNotchFound = false;
        double ppgEnd = -1.0;
        double vpgU = -1.0, vpgV = -1.0, vpgW = -1.0;
        double apgA = -1.0, apgB = -1.0, apgC = -1.0, apgD = -1.0, apgE = -1.0,
            apgF = -1.0;
        double jpgP1 = -1.0, jpgP2 = -1.0;
    };

    GlyphSnapshot m_glyphs;

    // Compute the glyph snapshot from current trace + marker state.
    void captureGlyphSnapshot(const TemplateBin& b);

    // Arterial trace vectors (own sample space; drawn foot-anchored at the
    // PPG origin). Empty when the channel is absent.
    std::vector<double> m_abp;
    std::vector<double> m_art;
    std::vector<double> m_artPulm;
    // Per-sample std for each arterial trace (empty => no band drawn).
    std::vector<double> m_abpIqr;
    std::vector<double> m_artIqr;
    std::vector<double> m_artPulmIqr;

    State m_state = State::Good;
    bool  m_hasPPG = false;
    bool  m_showEcgMarkers = true;
    bool  m_showPpgMarkers = true;
    bool  m_showPpgDerivMarkers = false;
    bool  m_showAbpMarkers = true;
    bool  m_showArtMarkers = true;
    bool  m_showArtPulmMarkers = true;
    bool  m_showEcgTrace = true;
    bool  m_showPpgTrace = true;
    bool  m_showAbpTrace = true;
    bool  m_showArtTrace = true;
    bool  m_showArtPulmTrace = true;
    int   m_dragMarker = -1;
    // Slice counts (post drop-rules) fed to the median for this widget's
    // ECG channel and the PPG. Displayed in the title when non-zero.
    int   m_nEcgBeats = 0;
    int   m_nPpgBeats = 0;

    static constexpr int margin_left = 26, margin_right = 26, margin_top = 20, margin_bottom = 26;
    static constexpr int click_radius_around_marker = 12;
};