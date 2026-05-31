// ============================================================================
// BinPlotWidget.h - One ECG lead + optional PPG overlay + draggable markers
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
#include <algorithm>
#include <cmath>
#include <vector>

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
        EcgTBegin = 2,
        EcgTEnd = 3,
        // --- PPG markers (contiguous, immediately after ECG) ---
        PpgOnset = 4,
        PpgPeak = 5,
        PpgDicrotic = 6,   // dicrotic notch
        Ppg50 = 7,   // 50% point
        PpgEnd = 8,   // end-of-pulse / trough after descent
        // --- size sentinel ---
        MarkerCount = 9
    };

    // Range-based predicates. Update these bounds if you add more
    // markers to either group.
    static bool markerIsEcg(int m) { return m >= EcgP && m <= EcgTEnd; }
    static bool markerIsPpg(int m) { return m >= PpgOnset && m <= PpgEnd; }

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


    static int computeEcgVisibleN(const std::vector<double>& ecg, int tEnd) {
        const int N = static_cast<int>(ecg.size());
        if (N < 4) return std::max(N, 2);
        const int half = N / 2;

        // Find the next-beat R as the largest deviation from the back-half mean.
        double mean = 0.0;
        for (int i = half; i < N; ++i) mean += ecg[i];
        mean /= (N - half);

        int nextR = half;
        double bestDev = std::abs(ecg[half] - mean);
        for (int i = half + 1; i < N; ++i) {
            const double d = std::abs(ecg[i] - mean);
            if (d > bestDev) { bestDev = d; nextR = i; }
        }

        return std::clamp(nextR, half, N);
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
        int ecgP, int qBegin, int tBegin, int tEnd,
        int ppgOnset, int ppgPeak,
        int ppgDicrotic, int ppg50, int ppgEnd,
        double rPeakSample);

    void setHasPPG(bool has);
    bool hasPPG() const { return m_hasPPG; }

    void setState(State s);
    State state() const { return m_state; }

    void setMarker(Marker m, int idx);
    int  marker(Marker m) const { return m_markers[m]; }

    // Per-trace marker visibility. When false, that group's markers
    // are neither drawn nor hit-testable (drag-pick ignores them).
    // Both default to true.
    void setShowEcgMarkers(bool show);
    void setShowPpgMarkers(bool show);

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

    // Widest sample extent any trace needs (in samples), and the
    // pixels-per-sample that makes that extent fill the drawable width.
    int    totalSampleSpan() const;
    double pxPerSample() const;

    double m_rPeakSample = 0.0;   // R-peak sample index within the ECG template

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
    int m_markers[MarkerCount] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };

    State m_state = State::Good;
    bool  m_hasPPG = false;
    bool  m_showEcgMarkers = true;
    bool  m_showPpgMarkers = true;
    int   m_dragMarker = -1;

    static constexpr int kML = 5, kMR = 5, kMT = 16, kMB = 4;
    static constexpr int kDragPx = 12;
};