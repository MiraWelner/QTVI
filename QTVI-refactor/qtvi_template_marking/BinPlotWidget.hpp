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
//   The widget uses a FIXED pixels-per-sample scale (kPxPerSample). The widget
//   reports a sizeHint() that grows to whatever width is required to draw both
//   traces in full at that scale. This is what guarantees temporal alignment
//   between the ECG and PPG -- one sample of ECG occupies the same pixel
//   width as one sample of PPG, regardless of how long each visible trace is.
// ============================================================================
#pragma once
#include <QWidget>
#include <iostream>
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

    // Fixed pixels-per-sample. Both ECG and PPG use this same value so
    // that one second of ECG occupies the same chart width as one second
    // of PPG. The widget's sizeHint grows to accommodate whichever trace
    // extends further. Tune this if bins come out too wide or too narrow.
    static constexpr double kPxPerSample = 0.4;

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

    // The old setData keeps its 5-marker positional API for compatibility.
    // The new markers default to -1 (hidden) until set explicitly via
    // setMarker(). Use setDataAll if you want to provide all of them at
    // once.
    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ecg,
        int qBegin, int tBegin, int tEnd,
        int ppgOnset, int ppgPeak, double rPeakSample);

    // Extended setter that takes every marker at once. Pass -1 for any
    // marker you don't have a position for.
    void setDataAll(const std::vector<double>& ppg,
        const std::vector<double>& ecg,
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

    int  binIndex()  const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

    int  ecgVisibleN() const { return m_ecgVisibleN; }

    // Width required to draw both traces in full at kPxPerSample.
    // Used by sizeHint and minimumSizeHint.
    int requiredWidth() const;

    QSize sizeHint() const override { return QSize(requiredWidth(), 120); }
    QSize minimumSizeHint() const override { return QSize(requiredWidth(), 60); }

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
    double m_rPeakSample = 0.0;   // ECG R-peak sample index (the PPG draw origin)


    int m_binIndex;
    int m_leadIndex;
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ecg;

    int m_ecgVisibleN = 0;
    int m_ppgVisibleN = 0;

    // Storage sized by MarkerCount so it grows automatically if you add
    // more entries to the enum. All marker slots start hidden (-1).
    int m_markers[MarkerCount] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };

    State m_state = State::Good;
    bool  m_hasPPG = false;
    int   m_dragMarker = -1;

    static constexpr int kML = 5, kMR = 5, kMT = 16, kMB = 4;
    static constexpr int kDragPx = 12;
};