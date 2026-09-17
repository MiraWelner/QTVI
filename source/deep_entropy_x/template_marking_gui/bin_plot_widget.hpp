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
#include "template_anchoring\curve_fit.hpp"   // curve_fit::FitMode / PeakFitMode

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
    // trace by construction. There is nothing to subtract, no shared origin to
    // assert, and no per-channel clip.
    //
    // WHY IT IS DONE THIS WAY. Position used to be expressed four ways at once:
    // a sample index, a per-channel start offset, a per-channel rate ratio, and
    // a clip count. Every alignment fix added another correction term, and the
    // terms disagreed -- the ECG's R column is 0.3 * the bin's LONGEST RR
    // (alignment.hpp) while every pulse channel's is a fixed 0.3 s
    // (create_arterial_templates.hpp), so asserting a shared sample 0 put the
    // ECG most of a second ahead of the PPG on any bin holding a pause. Bands
    // drifted from their own traces for the same reason: band and line reached
    // the array through different terms.
    //
    // The frame [tMin, tMax] is the UNION of every present channel's own drawn
    // extent. Because it is a union, no channel can have a tail outside it --
    // which is what makes a clip count unnecessary rather than merely
    // inconvenient, and why a marker inside its own array is on screen.
    //
    // A channel with no rate or no anchor is NOT DRAWN, deliberately: guessing
    // an anchor is what produced the misalignment this replaces.
    // ----------------------------------------------------------------------

    // The slot selects which waveform ecgDetect measures, so it invalidates
    // m_det -- but NOT m_glyphs, which is measured on m_ecg and would then
    // re-detect over the per-slot glyphs overrideEcgGlyphs pushes in.
    void setTemplateIndex(int t) {
        if (t == m_templateIndex) return;
        m_templateIndex = t;
        m_detValid = false;
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

    // Recaptures the frozen glyph snapshot and repaints. Call this LAST in a
    // seeding pass (see applyBinToWidget): m_glyphs.ecgRPeak reads the R bar,
    // which must already be set.
    // frame = the alignment whose waveform this panel is drawing. The FROZEN
    // glyphs are re-expressed in it, because each landmark is detected once
    // per alignment and the positions differ; the draggable BARS are not, and
    // stay in the R frame where they are stored and edited.
    void setAuto(const TemplateBin& b,
        AnchorType frame = AnchorType::R_PEAK) {
        // WHICH TEMPLATE THIS PANEL IS, kept so reactiveGlyphs can call the
        // shared ecgFiducials instead of assembling the detector's arguments
        // from its own copies (m_ecg, m_rPeakSample). Those copies are what
        // drifted: the bin-wide array instead of the slot's, R's column instead
        // of the alignment's. Non-owning -- the bins outlive the panels, and
        // setAuto is called on every rebuild.
        // A new bin or alignment is a new waveform. Covers the in-place
        // re-skin, which changes the frame without going through setData.
        if (m_bin != &b || m_frame != frame) m_detValid = false;
        m_bin = &b;
        m_frame = frame;
        captureGlyphSnapshot(b, frame);
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
        double ecgPBegin = -1.0;  // onset before the reactive P peak (same peak)
        double ecgPPeak = -1.0;   // between the P-onset and Q-onset bars
        double ecgTPeak = -1.0;   // between the S-end and T-end bars
        double ppgT50 = -1.0, ppgT80 = -1.0, ppgPeak2 = -1.0;
    };
    Reactive reactiveGlyphs() const;

    // Fit models for the DETECTED glyphs (the X marks): the transition fiducials
    // are drawn at the selected model's crossing. Setting them invalidates the
    // glyph snapshot so the next paint re-detects with the new model.
    void setFitModes(curve_fit::FitMode onOffset, curve_fit::PeakFitMode peak) {
        m_onOffsetFitMode = onOffset; m_peakFitMode = peak;
        // Both caches: the fit modes are detector inputs.
        m_glyphsValid = false; m_detValid = false; update();
    }

    // Per-trace marker visibility. When false, that group's markers
    // are neither drawn nor hit-testable (drag-pick ignores them).
    // Both default to true.
    void setShowEcgMarkers(bool show);
    void setShowPpgMarkers(bool show);

    // R-aligned overlay: the 4 ECG landmarks (P-onset, Q-onset, S-end, T-end)
    // as detected on the R alignment (the _R CSV columns), drawn read-only in
    // addition to the normal bars. Positions are already in the frame the panel
    // draws. setShowRMarkers toggles them; dragging one is handled via the
    // rMarkerDragStarted / rMarkerMoved signals below.
    void setShowRMarkers(bool show);
    void setRMarks(double pBegin, double qOnset, double sEnd, double tEnd);
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
    void overridePulseGlyphs(const tbank::BankPulseMarkerSet& pm);

    // Replace the BIN's ECG glyphs with columns measured on THIS panel's own
    // waveform. Counterpart to overridePulseGlyphs, and needed for exactly the
    // same reason: captureGlyphSnapshot reads the bin's *_auto_ch fields, which
    // describe the bin's anchored channel average, and a bank slot draws its
    // own. The two coincide under R alignment -- every beat in the bank has R
    // on one column by construction -- and under nothing else.
    //
    // A POD rather than FeatureMarks::TemplateLandmarks: this header is
    // included by every panel and pulling feature_marks.hpp in would drag
    // template_bank.hpp and annotation_types.hpp along with it. The caller
    // already holds the detector result and does the conversion.
    struct EcgGlyphColumns {
        double p_begin = -1.0, q_onset = -1.0, q_peak = -1.0,
            r_peak = -1.0, s_end = -1.0, t_end = -1.0;
        bool q_onset_found = false;
    };
    void overrideEcgGlyphs(const EcgGlyphColumns& g);


signals:
    // NEVER EMITTED -- superseded by markerMovedOnTemplate, which carries the
    // slot. Kept only so an existing connect() still compiles.
    void markerMoved(int binIndex, int leadIndex, int marker, int newIdx);
    void markerDragStarted(int binIndex, int leadIndex, int marker);

    // End of a bar gesture, once from mouseReleaseEvent. The owner defers its
    // one full page re-apply to this instead of paying it per mouse-move.
    void markerDragFinished(int binIndex, int leadIndex, int templateIdx,
        int marker);

    // R-aligned overlay (ecg_r_markers) is DRAGGABLE. Starting a drag on one
    // emits rMarkerDragStarted so the owner can flip the view to R while you
    // drag; each move emits rMarkerMoved with the landmark index (0..3 =
    // p_begin,q_onset,s_end,t_end) and its new column in the drawn frame.
    void rMarkerDragStarted(int binIndex, int leadIndex);
    void rMarkerMoved(int binIndex, int leadIndex, int templateIdx,
        int rIndex, double newCol);

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
    // `col` IS A DOUBLE: it is the marker's position, and the focus panel
    // measures against it. TemplateViewerWindow::onLandmarkSelected must
    // match -- a Qt signal and slot with different parameter types connect at
    // runtime and then silently never fire.
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
    int    markerAtX(double x) const;

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
        // ecgPPeak is NOT here: the P peak is REACTIVE, bracketed by the
        // P-onset and Q-onset bars, so a frozen copy would go stale the moment
        // either bar moved. It stays reactive on purpose -- see reactiveGlyphs.
        //
        // What made the reactive X land on the PR baseline was the BRACKET, not
        // the search: compute_p_begin fits the onset in a window centred on the
        // rough seed, so a seed off the P wave produced an onset past the bump
        // and the argmax between it and Q-onset had nothing but baseline to
        // pick from. detect_template_landmarks now refits the onset from the
        // re-measured peak, so the bracket contains the P wave and the reactive
        // search finds it.
        double ecgPBegin = -1.0, ecgQPeak = -1.0, ecgQ = -1.0,
            ecgRPeak = -1.0, ecgS = -1.0, ecgTend = -1.0;
        bool ecgQFound = false;
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

    // The glyph snapshot depends ONLY on the trace/bin, never on bar positions,
    // so it is recomputed only when the trace is (re)set -- on build and on an
    // alignment change (a bar click). setData/setEcgData clear this; during a
    // drag nothing re-sets the trace, so captureGlyphSnapshot returns early and
    // the expensive detect does not run per mouse-move.
    bool m_glyphsValid = false;

    // The reactive glyphs' expensive half, held while the trace is unchanged:
    // only P and T peak react to the bars (see ecgDetect). m_det.tmpl points
    // INTO the bin, so the identity fields are part of the guard -- a rebuild
    // can hand this panel another bin or slot without going through setData.
    mutable EcgDetection       m_det;
    mutable bool               m_detValid = false;
    mutable const TemplateBin* m_detBin = nullptr;
    mutable AnchorType         m_detFrame = AnchorType::R_PEAK;
    mutable int                m_detSlot = -1;

    curve_fit::FitMode     m_onOffsetFitMode = curve_fit::FitMode::Auto;
    curve_fit::PeakFitMode m_peakFitMode = curve_fit::PeakFitMode::Auto;

    // R-aligned overlay markers (ecg_r_markers). Read-only; index order is
    // P-onset, Q-onset, S-end, T-end. -1 = absent. Drawn only when shown.
    bool   m_showRMarkers = false;
    double m_rMarks[4] = { -1.0, -1.0, -1.0, -1.0 };
    int    m_dragRMark = -1;   // which R marker (0..3) is being dragged, -1 none

    // Compute the glyph snapshot from current trace + marker state.
    void captureGlyphSnapshot(const TemplateBin& b,
        AnchorType frame = AnchorType::R_PEAK);

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