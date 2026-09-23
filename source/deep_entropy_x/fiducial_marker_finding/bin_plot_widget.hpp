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
//   The widget fills whatever width its layout cell gives it, and every
//   channel is drawn in SECONDS RELATIVE TO ITS OWN R -- so R lands at the
//   same x on every trace by construction, whatever their sample rates. The
//   frame is the union of every channel's own extent, scaled to the drawable
//   width, so each bin window comes out the same on-screen length. See the
//   geometry note in the class body for why this replaced sample-index
//   drawing with per-channel offsets and ratios.
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
#include "fiducial_marker_finding/curve_fit.hpp"   // curve_fit::FitMode / PeakFitMode

class QPainter;

class BinPlotWidget : public QWidget {
    Q_OBJECT
public:

    // Right-click cycle: Good -> BadR -> BadPPG -> BadBoth -> Good.
    // BadPPG is PULSE ONLY (ECG good) and BadBoth is both, so the two verdicts
    // are independent rather than alternatives. With no pulse channel the
    // cycle is just Good -> BadR -> Good.
    enum class State { Good, BadR, BadPPG, BadBoth };

    // Each enum value MUST be unique (it's used as an array index into
    // m_markers). ECG markers come first, then PPG markers, so
    // markerIsEcg / markerIsPpg can use range checks.
    enum Marker : int {
        EcgPBegin = 0,
        EcgPPeak = 1,
        EcgQBegin = 2,
        EcgRPeak = 3,
        EcgSEnd = 4,
        // T begin removed and the numbering CLOSED UP, not left as a hole.
        // Nothing persists a marker id -- m_touchedMarks is keyed within a
        // session and every file format stores landmarks by field, not by enum
        // value -- so the renumbering is invisible outside this process.
        // anchor_view.hpp mirrors these values and BinPlotWidget.cpp
        // static_asserts them, so a missed shift is a build error.
        EcgTEnd = 5,
        // --- PPG markers (contiguous, immediately after ECG) ---
        PpgOnset = 6,
        PpgT50 = 7,
        PpgPeak = 8,
        PpgDicrotic = 9,
        PpgPeak2 = 10,
        PpgT80 = 11,
        PpgEnd = 12,
        // --- Arterial markers ---
        AbpOnset = 13, AbpPeak = 14, AbpDicrotic = 15, AbpPeak2 = 16, AbpEnd = 17,
        ArtOnset = 18, ArtPeak = 19, ArtDicrotic = 20, ArtPeak2 = 21, ArtEnd = 22,
        ArtPulmOnset = 23, ArtPulmPeak = 24, ArtPulmDicrotic = 25,
        ArtPulmPeak2 = 26, ArtPulmEnd = 27,
        // Q-peak and T-peak glyph ids, APPENDED (not in the [0,5] ECG range) so
        // the shared 0..5 values anchor_view mirrors and static_asserts stay
        // put. They exist only so a click on those glyphs can open a read-only
        // focus view; they are never bars, never stored in a file, and
        // markerIsEcg is widened by hand to include them.
        EcgQPeak = 28, EcgTPeak = 29,
        MarkerCount = 30
    };

    static bool markerIsEcg(int m) {
        return (m >= EcgPBegin && m <= EcgTEnd)   // contiguous bars + P/R-peak glyphs
            || m == EcgQPeak || m == EcgTPeak;     // appended glyph ids (focus only)
    }
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
        return m == EcgPBegin || m == EcgQBegin
            || m == PpgOnset || m == AbpOnset || m == ArtOnset || m == ArtPulmOnset;
    }
    double m_rPeakSample = 0.0;
    // Set by setAuto; null until then. See the note there.
    const TemplateBin* m_bin = nullptr;
    AnchorType m_frame = AnchorType::R_PEAK;   // R-peak sample index within the ECG template
    // (m_ppgDelay / m_ppgFootIdx retired in Patch C: every channel is
    // real-time-aligned by construction under Patch B slicing.)

    // Every trace this widget can draw. Add new channels here only, before
    // Count; nothing else in the geometry code needs to change, because a
    // channel is fully described by its rate and its R column (see the time
    // model below).
    enum class Channel { Ecg, Ppg, Abp, Art, ArtPulm, Count };

    // Resolve a marker to its channel, trace, and group visibility. No
// visible-sample bound and no ratio: the frame is the union of every
// channel's extent, so a marker inside its own array is on screen.
    int    lastDrawnSample(Channel ch) const;
    int    firstDrawnSample(Channel ch) const;



    // ----------------------------------------------------------------------
    // GEOMETRY: EVERY CHANNEL IS DRAWN IN SECONDS RELATIVE TO ITS OWN R.
    //
    //     t(ch, i) = (i - anchor[ch]) / rate[ch]
    //     x(t)     = margin_left + (t - tMin) / (tMax - tMin) * drawW
    //
    // R lands at t = 0 on every channel, so it lands at the same x on every
    // trace by construction. Nothing to subtract, no shared origin to assert,
    // no per-channel clip.
    //
    // DO NOT ASSERT A SHARED SAMPLE 0. The ECG's R column is 0.3 * the bin's
    // LONGEST RR (alignment.hpp) while every pulse channel's is a fixed 0.3 s
    // (create_arterial_templates.hpp), so on any bin holding a pause that puts
    // the ECG most of a second ahead of the PPG.
    //
    // The frame [tMin, tMax] is the UNION of every present channel's own drawn
    // extent, which is why no clip count is needed: a channel cannot have a
    // tail outside it, so a marker inside its own array is on screen.
    //
    // A channel with no rate or no anchor is NOT DRAWN. Guessing an anchor is
    // what produced the misalignment this geometry replaces.
    // ----------------------------------------------------------------------

    // The slot selects which waveform BOTH detectors measure, so it drops both
    // caches. There is no third copy to keep in step: the glyphs are DRAWN from
    // these two detections, so nothing can be invalidated on one schedule and
    // painted on another.
    void setTemplateIndex(int t) {
        if (t == m_templateIndex) return;
        m_templateIndex = t;
        m_detValid = false;
        m_pdetValid = false;
    }
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

    // Replace ONLY the ECG trace/band/R-column/count, leaving PPG, arterial
    // channels and all markers untouched. For Automatic alignment re-anchoring
    // the grid on a bar click: the ECG average is the one thing that changes
    // per anchor, and updating in place (rather than rebuilding the panel)
    // keeps an in-progress drag alive.
    void setEcgData(const std::vector<double>& ecg,
        const std::vector<double>& ecgIqr,
        double rPeakSample,
        int nEcgBeats = 0);

    // The pulse twin of setEcgData: replace ONLY the PPG trace, its band and
    // its beat count, leaving the ECG, the arterial channels and every marker
    // where they are. For an operator re-stack on mouse-up (see
    // TemplateViewerWindow::realignPulseFromFoot): the pulse average changes
    // and nothing else does, and the bars must NOT move -- the whole point of
    // re-anchoring on the dragged foot is that the trace comes to the bar.
    void setPpgData(const std::vector<double>& ppg,
        const std::vector<double>& ppgIqr,
        int nPpgBeats = 0);
    bool hasPPG() const { return m_hasPPG; }

    // Pin the ECG channel's contribution to the x-frame to a FIXED window, in
    // seconds relative to R, instead of deriving it from the current ECG
    // trace's finite extent. The four alignments (P/Q/R/J) produce averages
    // with slightly different left/right extents, so without this the x-axis
    // rescales every time the anchor changes and the waveform appears to slide.
    // The caller passes the union of all four anchors' extents, so whichever
    // anchor is shown draws inside the same window and the axis holds still.
    // PPG and arterial channels still union in normally -- only the ECG span is
    // pinned.
    void setEcgFrame(double tMinSec, double tMaxSec);

    void setChannelRate(Channel ch, double hz);
    double channelRate(Channel ch) const { return m_rates[static_cast<size_t>(ch)]; }

    // Frame bounds in seconds relative to R; negative before it.
    double frameTMin() const { return m_tMin; }
    double frameTMax() const { return m_tMax; }

    // Pixels per SECOND -- the one number describing how zoomed a panel is,
    // shared by every channel. channelDx(ch) == pxPerSecond() / rate[ch], so a
    // per-channel pixels-per-sample is still available where one is wanted.
    double pxPerSecond() const;

    void setState(State s);
    State state() const { return m_state; }

    // SUB-SAMPLE POSITIONS. Every marker the widget holds is a double: the
    // bars come from BankMarkerSet (double) and the glyphs from the detectors
    // and reactive_* (double). A drag still lands on a whole column, because
    // sampleFromX maps a pixel to a sample -- but nothing else quantises.
    void setMarker(Marker m, double idx);
    double marker(Marker m) const { return m_markers[m]; }

    // WHICH TEMPLATE THIS PANEL IS. frame = the alignment whose waveform it
    // draws; both detectors measure on that alignment's own array, through
    // slotView, so the landmarks come back in the columns this panel plots.
    // Non-owning -- the bins outlive the panels, and setAuto runs on every
    // rebuild.
    //
    // Call order does not matter: this captures nothing. It points the panel at
    // a waveform and drops the caches; the detections are taken on demand by
    // whoever draws or reads them.
    void setAuto(const TemplateBin& b,
        AnchorType frame = AnchorType::R_PEAK) {
        // A new bin or alignment is a new waveform. Covers the in-place
        // re-skin, which changes the frame without going through setData.
        if (m_bin != &b || m_frame != frame) {
            m_detValid = false;
            m_pdetValid = false;
        }
        m_bin = &b;
        m_frame = frame;
        update();
    }

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
        // (ecgPBegin removed: it was a copy of m_det.lm.p_begin under a second
        //  name, which read as though the P onset were bar-derived like the two
        //  peaks below. It is not -- it is detected, and the snapshot holds it.)
        double ecgPPeak = -1.0;   // between the P-onset and Q-onset bars
        double ecgTPeak = -1.0;   // between the S-end and T-end bars
        double ppgT50 = -1.0, ppgT80 = -1.0, ppgPeak2 = -1.0;
    };
    Reactive reactiveGlyphs() const;

    // THE COLUMNS THE ECG X GLYPHS ARE DRAWN AT -- drawFeatureGlyphs paints
    // straight out of this, so a reader of it cannot be looking at a different
    // number from the operator.
    const FeatureMarks::TemplateLandmarks& detectedLandmarks() const {
        reactiveGlyphs();          // populates / reuses m_det
        return m_det.lm;
    }

    // THE COLUMNS THE PULSE GLYPHS ARE DRAWN AT, on the same terms: one
    // detection of this slot's own pulse average, cached on (bin, slot), drawn
    // from and read from. The pulse side used to have a setter and no getter --
    // overridePulseGlyphs pushed a stored BankPulseMarkerSet in and nothing
    // could read it back -- which is why the focus panel had no fiducial to
    // draw and fell back to the bar column.
    const FeatureMarks::PpgFiducials& detectedPulse() const;

    // The fitted curves behind one peak glyph, for a viewer that wants to DRAW
    // the contest rather than re-run it. Same cache, same call.
    const upsample_for_fit::PeakCandidates& peakCandidatesFor(EcgPeak w) const {
        reactiveGlyphs();          // populates / reuses m_det + m_peakCands
        return m_peakCands[static_cast<size_t>(w)];
    }

    // Fit models for the DETECTED glyphs (the X marks): the transition fiducials
    // are drawn at the selected model's crossing. Setting them invalidates the
    // glyph snapshot so the next paint re-detects with the new model.
    void setFitModes(curve_fit::FitMode onOffset, curve_fit::PeakFitMode peak) {
        m_onOffsetFitMode = onOffset; m_peakFitMode = peak;
        // The fit modes are detector inputs, so the detection is dropped. The
        // glyphs follow it with nothing to invalidate of their own.
        m_detValid = false; update();
    }

    // Per-trace marker visibility. When false, that group's markers
    // are neither drawn nor hit-testable (drag-pick ignores them).
    // Both default to true.
    void setShowEcgMarkers(bool show);
    void setShowPpgMarkers(bool show);

    // THIS PANEL IS SHOWING ONE ALIGNMENT'S OWN BARS. Non-empty (the
    // alignment's letter, "P"/"Q"/"R"/"J") switches the ECG bars to the
    // overlay style the R-aligned overlay used to have: dotted, all one
    // colour, each labelled with this letter plus its landmark -- "Pp", "Qq",
    // "Qs", "Jt". Empty means Automatic, where the bars are canonical copies
    // assembled from four different alignments and the per-landmark colours
    // are the cue that matters.
    void setAlignmentBadge(const char* label) {
        const QString s = label ? QString::fromLatin1(label) : QString();
        if (s == m_alignBadge) return;
        m_alignBadge = s;
        update();
    }

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

    // Section 4.6 bank overlay: templates 1..N-1 of this (bin, channel)'s bank,
    // drawn on the ECG axis under the slot 0 trace.
    //
    // Deliberately NOT reusing m_bgTraces. That one carries the arterial
    // context traces (ABP / ART / ART_PULM), which are foot-anchored and belong
    // on a different axis -- and it is set from the plot-construction path,
    // which runs BEFORE applyBinToWidget(), so sharing the member would mean
    // whichever wrote last silently erased the other.

    int  binIndex()  const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

    // The widget no longer dictates its width from the trace length.
    // It advertises a modest preferred width and a small minimum so the
    // grid can hand every cell an equal share of the window; the trace is
    // then scaled to fill whatever width the cell receives.
    QSize sizeHint() const override { return QSize(220, 120); }
    QSize minimumSizeHint() const override { return QSize(40, 60); }
    // DO NOT ADD A GLYPH SETTER HERE. Two once existed (overridePulseGlyphs,
    // overrideEcgGlyphs) writing into a frozen snapshot, and they gave every
    // landmark two positions on two invalidation schedules -- the snapshot
    // dropped on trace changes, the detections on (bin, slot, frame) -- so the
    // X on screen and the number every reader got came from two runs of the
    // same detector over different windows. Glyphs are painted from
    // detectedLandmarks() / detectedPulse(), which is where the focus panel
    // reads them, so the two are one number by construction.

signals:
    void markerDragStarted(int binIndex, int leadIndex, int marker);


    // The operator selected (clicked) a landmark: render its focus panel.
    // Fires on the same click as markerDragStarted but carries a different
    // intent -- "show this landmark's close-up" rather than "a drag is
    // beginning" -- and the owner decides what to do with it.
    //
    // TEMPLATE INDEX INCLUDED. A panel is a (bin, template) pair and the focus
    // view reads mean/sd/n for the waveform under the landmark; without the
    // slot it reads the BIN's waveform whichever panel was clicked.
    //
    // `col` IS A DOUBLE and the receiving slot must match: a Qt signal and slot
    // with different parameter types connect at runtime and then silently never
    // fire.
    void landmarkSelected(int binIndex, int leadIndex, int templateIdx,
        int marker, double col);

    // A GLYPH was clicked (read-only): open the focus panel WITHOUT recording a
    // touch or re-aligning. Glyphs (peaks and the detected transition fiducials)
    // are the detector's answer, independent of the draggable bars, so clicking
    // one must not mark the bar as operator-moved.
    void landmarkFocusOnly(int binIndex, int leadIndex, int templateIdx,
        int marker, double col);

    // TEMPLATE INDEX ADDED to the marker signals. A panel is a (bin, template)
    // pair now, and without the slot the receiver cannot tell whether a drag
    // belongs to the bin's own marker set or to a bank template's -- so every
    // drag on a sub-template column would have written into slot 0's landmarks.
    void markerMovedOnTemplate(int binIndex, int leadIndex, int templateIdx,
        int marker, int newIdx);

    // GESTURE FINISHED, and only when the bar actually moved. Emitted from
    // mouseReleaseEvent for work that is too expensive to do per mouse-move --
    // the pulse re-stack reads a beat matrix off disk and re-medians it, which
    // cannot run at drag rate.
    //
    // THE OLD markerDragFinished WAS REMOVED FOR A GOOD REASON and this is not
    // a revival of it: m_dragMarker is armed by any bar CLICK, and an automatic
    // alignment shift IS a bar click, so that signal fired a full-page
    // re-detect on every one of them. This fires only if mouseMoveEvent
    // actually changed the bar's column, so a click that moved nothing is
    // silent.
    void markerReleasedOnTemplate(int binIndex, int leadIndex, int templateIdx,
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


protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    std::vector<global_interval_lines::Line> m_refLines; //global refernce lines eg. earliest Q-onset 
    // Time of sample i on `ch`, in seconds relative to that channel's R.
    // NaN when the channel has no rate or no anchor.
    double timeAt(Channel ch, double i) const;

    // The only place the frame becomes pixels.
    double xFromTime(double t) const;

    // x(i) is AFFINE in i, so a trace draw needs only these two numbers:
    //     x(i) = channelX0(ch) + i * channelDx(ch)
    // which is what lets the two draw helpers keep their existing form.
    double channelX0(Channel ch) const;
    double channelDx(Channel ch) const;

    // Pixel x of (sub-sample) index i on `ch`, and the inverse.
    double xFromSample(Channel ch, double i) const;
    int    sampleFromX(Channel ch, double x) const;

    // The ONLY hit-test. Bars are clickable; glyphs are display-only and are
    // deliberately not hit-tested, so a click can never select or drag an
    // automated mark.
    // The nearest BAR to x, or -1. distOut receives that bar's pixel distance
    // (infinity on a miss) so the caller can weigh it against a glyph hit
    // instead of letting bars win unconditionally.
    int    markerAtX(double x, double* distOut = nullptr) const;

    bool   markerTrace(int m, const std::vector<double>*& vec,
        Channel& ch, bool& visible) const;

    // Recompute [m_tMin, m_tMax] from every present channel. Called by
    // setData, setArterialTraces, setChannelRate and, so the
    // frame can never be stale with respect to the traces.
    void   recomputeFrame();

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

    double m_markers[MarkerCount];   // all -1 until seeded (filled in the ctor)

    // Hz per channel (indexed by Channel); 0 = unknown -> channel not drawn.
    std::array<double, static_cast<size_t>(Channel::Count)> m_rates{};

    // R column per channel, in that channel's own samples; -1 = unknown.
    // Filled in the ctor with -1 and set by setChannelRate / setData.
    std::array<double, static_cast<size_t>(Channel::Count)> m_rAnchor{};

    // Frame bounds in seconds relative to R. Recomputed by recomputeFrame().
    double m_tMin = 0.0;
    double m_tMax = 1.0;

    // Fixed ECG x-extent (seconds rel R), set by setEcgFrame. When active,
    // recomputeFrame folds this in for the ECG channel instead of the current
    // trace's own finite span, so switching alignment does not move the axis.
    bool   m_ecgFrameFixed = false;
    double m_ecgFrameLo = 0.0;
    double m_ecgFrameHi = 0.0;

    // The reactive glyphs' expensive half, held while the trace is unchanged:
    // only P and T peak react to the bars (see ecgDetect). m_det.tmpl points
    // INTO the bin, so the identity fields are part of the guard -- a rebuild
    // can hand this panel another bin or slot without going through setData.
    // The peak contest for this trace, indexed by EcgPeak. Filled by the same
    // pass that places the glyphs, so the fit the panel DRAWS is the fit that
    // placed the mark -- one contest, not one per viewer.
    mutable std::array<upsample_for_fit::PeakCandidates, 5> m_peakCands{};
    mutable EcgDetection       m_det;
    mutable bool               m_detValid = false;
    mutable const TemplateBin* m_detBin = nullptr;
    mutable AnchorType         m_detFrame = AnchorType::R_PEAK;
    mutable int                m_detSlot = -1;

    // THE PULSE TWIN, same cache and same guard minus the alignment: pulse
    // channels are foot-anchored once, so there is no per-anchor pulse average
    // and no frame to key on.
    mutable FeatureMarks::PpgFiducials m_pdet;
    mutable bool               m_pdetValid = false;
    mutable const TemplateBin* m_pdetBin = nullptr;
    mutable int                m_pdetSlot = -1;

    // ---- THE REACTIVE GLYPHS' OWN CACHE --------------------------------
    //
    // reactiveGlyphs() is called several times per paint and the bracketed
    // peaks inside it are two weighted polynomial fits each, so it is not
    // free. Keyed on everything it reads: the eight bars it brackets with,
    // the template identity, and whether the detection it brackets against
    // was already valid on entry. Anything that changes one of those changes
    // the key, so there is nothing to remember to invalidate.
    struct ReactiveKey {
        double pBegin, qOnset, sEnd, tEnd;
        double ppgOnset, ppgPeak, ppgDicrotic, ppgEnd;
        const TemplateBin* bin;
        AnchorType frame;
        int slot;
        bool detValid;
        bool operator==(const ReactiveKey& o) const {
            return pBegin == o.pBegin && qOnset == o.qOnset
                && sEnd == o.sEnd && tEnd == o.tEnd
                && ppgOnset == o.ppgOnset && ppgPeak == o.ppgPeak
                && ppgDicrotic == o.ppgDicrotic && ppgEnd == o.ppgEnd
                && bin == o.bin && frame == o.frame && slot == o.slot
                && detValid == o.detValid;
        }
    };
    mutable Reactive    m_rx;
    mutable ReactiveKey m_rxKey{};
    mutable bool        m_rxValid = false;

    // Alignment letter for the overlay bar style; empty in Automatic.
    QString m_alignBadge;

    // ---- THE AXIS GEOMETRY THE LAST PAINT USED --------------------------
    //
    // A glyph is a five-pixel X at one point ON THE TRACE, so hit-testing it
    // needs its y, which needs the axis range and plot height -- and those are
    // computed in paintEvent, so they do not exist when a mouse press arrives.
    // Recorded by the paint that drew them, which also means the hit test uses
    // exactly the geometry the operator was looking at.
    //
    // m_lastPh <= 0 means nothing has been painted yet; the press handler falls
    // back to testing x alone.
    mutable double m_lastYLo = 0.0;
    mutable double m_lastYHi = 0.0;
    // The PULSE axis the same paint used. The hit test covers both channels
    // now, and a pulse glyph's y is on the right-hand scale, not the ECG's.
    mutable double m_lastPLo = 0.0;
    mutable double m_lastPHi = 0.0;
    mutable int    m_lastPh = 0;

    curve_fit::FitMode     m_onOffsetFitMode = curve_fit::FitMode::Auto;
    curve_fit::PeakFitMode m_peakFitMode = curve_fit::PeakFitMode::Auto;


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
    // DID THIS GESTURE MOVE ANYTHING. Set by mouseMoveEvent when a bar's
    // column actually changes, cleared on press and on release. It is what
    // separates "the operator dragged the foot" from "the operator, or an
    // automatic alignment shift, clicked it" -- see markerReleasedOnTemplate.
    bool  m_dragMoved = false;
    // Slice counts (post drop-rules) fed to the median for this widget's
    // ECG channel and the PPG. Displayed in the title when non-zero.
    int   m_nEcgBeats = 0;
    int   m_nPpgBeats = 0;

    static constexpr int margin_left = 26, margin_right = 26, margin_top = 20, margin_bottom = 26;
    static constexpr int click_radius_around_marker = 12;
};