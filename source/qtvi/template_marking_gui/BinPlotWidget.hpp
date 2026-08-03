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

    explicit BinPlotWidget(int binIndex, int leadIndex,
        const QString& leadLabel, QWidget* parent = nullptr);

    // Sole data setter. Pass empty std vectors to draw without a band
    // for that trace; pass -1 for any marker you don't have a position
    // for. Std vectors that are shorter than the visible sample count
    // are ignored at draw time (the band silently disappears for that
    // trace), so it's safe to call this with stale data.
    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ppgIqr,
        const std::vector<double>& ecg,
        const std::vector<double>& ecgIqr,
        int pPeak, int qBegin, int rPeak, int sEnd, int tPeak, int tEnd,
        int ppgOnset, int ppgT50, int ppgPeak,
        int ppgDicrotic, int ppgPeak2, int ppgT80, int ppgEnd,
        double rPeakSample,
        int nEcgBeats = 0,
        int nPpgBeats = 0,
        // The full set of PPG auto-detected fiducials (TemplateBin's
        // ppg_*_auto fields, all from ONE detect_ppg_fiducials() call --
        // see feature_marks.cpp). These are the single source of truth for
        // the frozen glyphs: captureGlyphSnapshot() reads them directly,
        // with no independent recompute. They're separate from
        // ppgOnset/ppgDicrotic/ppgPeak2/etc. above, which are the current
        // (possibly user-dragged) bar positions.
        int ppgOnsetAuto = -1,
        int ppgPeakAuto = -1,
        int ppgPeak2Auto = -1,
        bool ppgPeak2FoundAuto = false,
        int ppgDicroticAuto = -1,
        bool ppgDicroticFoundAuto = false,
        int ppgEndAuto = -1,
        bool ppgEndFoundAuto = false);

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

    // Load-time autodetected ECG landmark columns (the *_auto_ch fields).
    // Own-bar glyphs (P-begin, P-peak, Q-begin, S-end, T-end) are FROZEN at
    // these positions and do NOT follow their bars. -1 => fall back to the bar.
    struct EcgAuto {
        int pBegin = -1, pPeak = -1, qBegin = -1, sEnd = -1, tEnd = -1;
    };
    void setEcgAuto(const EcgAuto& a) { m_ecgAuto = a; captureGlyphSnapshot(); update(); }

    // Recompute the reactive feature glyphs from the current markers and
    // repaint. Call once after programmatically setting markers (e.g. from
    // refreshBinMarkers) so the X glyphs track the markers live.
    void refreshGlyphs() { captureGlyphSnapshot(); update(); }

    // Per-trace marker visibility. When false, that group's markers
    // are neither drawn nor hit-testable (drag-pick ignores them).
    // Both default to true.
    void setShowEcgMarkers(bool show);
    void setShowPpgMarkers(bool show);
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

    int  binIndex()  const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

    int  ecgVisibleN() const { return m_ecgVisibleN; }

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
    void landmarkSelected(int binIndex, int leadIndex, int marker, int col);

    void badRToggled(int binIndex, int leadIndex, bool bad);
    void badPPGToggled(int binIndex, bool bad);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    int    sampleFromX(double x, double startSample, double ratio) const;
    double xFromSample(int s, double startSample, double ratio) const;
    int    markerAtX(double x) const;
    // B2 focus mode: hit-test the feature glyphs (X marks). Returns a
    // Marker-enum routing id (or -1); outCol = the glyph's sample column,
    // outIsEcg = whether it uses ECG (vs PPG) geometry.
    int    glyphAtX(double x, int& outCol, bool& outIsEcg) const;
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
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ppgIqr;
    std::vector<double> m_ecg;
    std::vector<double> m_ecgIqr;

    int m_ecgVisibleN = 0;
    int m_ppgVisibleN = 0;
    int m_markers[MarkerCount];

    // Hz per channel (indexed by Channel); 0 = unknown -> rateRatio()
    // falls back to 1.0 for that channel.
    std::array<double, static_cast<size_t>(Channel::Count)> m_rates{};

    struct GlyphSnapshot {
        // All ECG glyphs read straight from the movable markers now
        // (P peak, Q onset, R peak, S end, T peak, T end). Each draws an
        // X at its marker's position; no O fallbacks (a marker is always
        // set, or the field stays -1 and the draw is skipped).
        int ecgPBegin = -1, ecgPPeak = -1, ecgQ = -1, ecgQPeak = -1, ecgRPeak = -1, ecgSPeak = -1, ecgS = -1, ecgTPeak = -1, ecgTend = -1;
        // PPG: sourced directly from the bin's single-source-of-truth
        // ppg_*_auto fields (see FeatureMarks::detect_ppg_fiducials) -- no
        // independent glyph recompute, so these can never disagree with
        // the auto-seeded movable bars (they're set once from setData()
        // and never change, since dragging a bar only touches the
        // non-auto TemplateBin fields).
        int ppgFoot = -1;    // = onset_auto
        int ppgP1 = -1;      // = peak_auto
        int ppgP2 = -1;      bool ppgPeak2Found = false;
        int ppgDic = -1;     bool ppgNotchFound = false;
        int ppgEnd = -1;     bool ppgEndFound = false;
        // T80/P50 are reactive (NOT from the auto fields): always tracked
        // live from the CURRENT onset/peak/end markers, recomputed every
        // capture -- same "reactive" treatment ECG's Q-peak/S-peak get.
        int ppgT50 = -1;
        int ppgT80 = -1;
        bool valid = false;
    };

    GlyphSnapshot m_glyphs;

    // The full set of PPG auto-detected fiducials (TemplateBin's
    // ppg_*_auto fields), captured once per setData() call. These feed the
    // frozen glyphs directly -- no recompute, no freeze flag needed, since
    // they never change after load (dragging a bar only touches the
    // separate, non-auto TemplateBin fields in m_markers[Ppg*]).
    int m_ppgOnsetAuto = -1, m_ppgPeakAuto = -1, m_ppgPeak2Auto = -1;
    bool m_ppgPeak2FoundAuto = false;
    int m_ppgDicroticAuto = -1;   bool m_ppgDicroticFoundAuto = false;
    int m_ppgEndAuto = -1;        bool m_ppgEndFoundAuto = false;

    // Frozen ECG autodetect positions (own-bar glyphs read from here).
    EcgAuto m_ecgAuto;

    // Compute the glyph snapshot from current trace + marker state.
    void captureGlyphSnapshot();

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