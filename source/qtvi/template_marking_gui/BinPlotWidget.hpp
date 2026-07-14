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
        // --- ECG markers (contiguous, starting at 0) ---
        EcgP = 0,   // P-wave onset
        EcgQBegin = 1,
        EcgSEnd = 2,   // QRS end (S-wave)
        EcgTBegin = 3,
        EcgTEnd = 4,
        // --- PPG markers (contiguous, immediately after ECG) ---
        PpgOnset = 5,
        PpgP50 = 6,   // 50% up the upslope, foot -> systolic peak
        PpgPeak = 7,
        PpgDicrotic = 8,   // dicrotic notch
        PpgPeak2 = 9,   // 2nd (diastolic) peak, after dicrotic notch
        PpgEnd = 10,   // end-of-pulse / trough after descent
        // --- Arterial markers: same 5-marker set as PPG, one group per
        //     channel (ABP, ART, ART_PULM), each contiguous. They ride the
        //     PPG x-geometry (foot-anchored at ppgStartSample()) but index
        //     into their own background trace vector. ---
        AbpOnset = 11, AbpPeak = 12, AbpDicrotic = 13, AbpPeak2 = 14, AbpEnd = 15,
        ArtOnset = 16, ArtPeak = 17, ArtDicrotic = 18, ArtPeak2 = 19, ArtEnd = 20,
        ArtPulmOnset = 21, ArtPulmPeak = 22, ArtPulmDicrotic = 23,
        ArtPulmPeak2 = 24, ArtPulmEnd = 25,
        // --- size sentinel ---
        MarkerCount = 26
    };

    // Range-based predicates. Update these bounds if you add more
    // markers to either group.
    static bool markerIsEcg(int m) { return m >= EcgP && m <= EcgTEnd; }
    static bool markerIsPpg(int m) { return m >= PpgOnset && m <= PpgEnd; }
    static bool markerIsAbp(int m) { return m >= AbpOnset && m <= AbpEnd; }
    static bool markerIsArt(int m) { return m >= ArtOnset && m <= ArtEnd; }
    static bool markerIsArtPulm(int m) { return m >= ArtPulmOnset && m <= ArtPulmEnd; }
    // Any arterial marker (ABP/ART/ART_PULM) rides the PPG x-geometry.
    static bool markerIsArterial(int m) { return m >= AbpOnset && m <= ArtPulmEnd; }
    double m_rPeakSample = 0.0;   // R-peak sample index within the ECG template
    double m_sampleRate = 0.0;   // Hz; 0 => label x-axis in samples
    // (m_ppgDelay / m_ppgFootIdx retired in Patch C: every channel is
    // real-time-aligned by construction under Patch B slicing.)

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
        const std::vector<double>& ppgStd,
        const std::vector<double>& ecg,
        const std::vector<double>& ecgStd,
        int ecgP, int qBegin, int sEnd, int tBegin, int tEnd,
        int ppgOnset, int ppgP50, int ppgPeak,
        int ppgDicrotic, int ppgPeak2, int ppgEnd,
        double rPeakSample,
        int nEcgBeats = 0,
        int nPpgBeats = 0);

    void setHasPPG(bool has);
    bool hasPPG() const { return m_hasPPG; }

    void setSampleRate(double hz);
    double sampleRate() const { return m_sampleRate; }

    void setState(State s);
    State state() const { return m_state; }

    void setMarker(Marker m, int idx);
    int  marker(Marker m) const { return m_markers[m]; }

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
        const std::vector<double>& abpStd = {},
        const std::vector<double>& artStd = {},
        const std::vector<double>& artPulmStd = {});
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

    void badRToggled(int binIndex, int leadIndex, bool bad);
    void badPPGToggled(int binIndex, bool bad);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    int    sampleFromX(double x, bool isEcg) const;
    double xFromSample(int s, bool isEcg) const;
    int    markerAtX(double x) const;
    int    visibleN(bool isEcg) const;
    // Resolve a marker's trace vector, geometry, current visibility, and
    // visible-sample bound. ECG, PPG, and all arterial channels resolve
    // through one path so they behave identically.
    bool   markerTrace(int m, const std::vector<double>*& vec,
        bool& isEcg, bool& visible, int& visN) const;

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
    std::vector<double> m_ppgStd;
    std::vector<double> m_ecg;
    std::vector<double> m_ecgStd;


    int m_ecgVisibleN = 0;
    int m_ppgVisibleN = 0;

    // Storage sized by MarkerCount so it grows automatically if you add
    // more entries to the enum. All marker slots start hidden (-1).
    int m_markers[MarkerCount] = {
        -1, -1, -1, -1, -1,       // ECG (5)
        -1, -1, -1, -1, -1, -1,   // PPG (6: onset, P50, peak, dic, peak2, end)
        -1, -1, -1, -1, -1,       // ABP
        -1, -1, -1, -1, -1,       // ART
        -1, -1, -1, -1, -1        // ART_PULM
    };

    struct GlyphSnapshot {
        int ecgP = -1, ecgR = -1, ecgT = -1, ecgQ = -1, ecgS = -1, ecgTend = -1;
        // O-fallback positions for P/R/T peaks: -1 => landmark found (draw X
        // at the ecgP/ecgR/ecgT index); >=0 => landmark absent, draw O here
        // instead. Q, S end, T end come straight from movable markers and
        // don't have fallbacks.
        int ecgPOFallback = -1, ecgROFallback = -1, ecgTOFallback = -1;
        int ppgFoot = -1, ppgP50 = -1, ppgP1 = -1, ppgDic = -1, ppgP2 = -1, ppgEnd = -1;
        // Fallback midpoints for landmarks that can be "expected but not
        // found". X drawn at the real index; O drawn at the fallback.
        int ppgP50OFallback = -1;
        int ppgP1OFallback = -1;
        int ppgP2OFallback = -1;
        bool ppgNotch = false;     // true => real notch found; false => draw 'o'
        int  ppgNoNotchO = -1;     // 'o' position (midpoint of [peak, end]) when no notch
        bool valid = false;
    };

    GlyphSnapshot m_glyphs;

    // Compute the glyph snapshot from current trace + marker state.
    void captureGlyphSnapshot();

    // Arterial trace vectors (own sample space; drawn foot-anchored at the
    // PPG origin). Empty when the channel is absent.
    std::vector<double> m_abp;
    std::vector<double> m_art;
    std::vector<double> m_artPulm;
    // Per-sample std for each arterial trace (empty => no band drawn).
    std::vector<double> m_abpStd;
    std::vector<double> m_artStd;
    std::vector<double> m_artPulmStd;

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

    static constexpr int margin_left = 40, margin_right = 40, margin_top = 16, margin_bottom = 16;
    static constexpr int click_radius_around_marker = 12;
};