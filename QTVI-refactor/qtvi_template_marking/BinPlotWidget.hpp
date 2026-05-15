// ============================================================================
// BinPlotWidget.h - One ECG lead + optional PPG overlay + 5 draggable markers
//
// Markers:
//   ECG (per channel):  Q-begin (green), T-begin (yellow), T-end (magenta)
//   PPG (shared/bin):   Onset (cyan), Peak (orange)
//
// Right-click cycles:  Good -> BadR -> BadPPG -> Good  (skips BadPPG if no PPG)
// Left-drag:           move whichever marker is closest to the click
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

    enum Marker : int {
        EcgQBegin = 0,
        EcgTBegin = 1,
        EcgTEnd = 2,
        PpgOnset = 3,
        PpgPeak = 4,
        MarkerCount = 5
    };

    static bool markerIsEcg(int m) { return m >= EcgQBegin && m <= EcgTEnd; }
    static bool markerIsPpg(int m) { return m == PpgOnset || m == PpgPeak; }

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

    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ecg,
        int qBegin, int tBegin, int tEnd,
        int ppgOnset, int ppgPeak);

    void setHasPPG(bool has);
    bool hasPPG() const { return m_hasPPG; }

    void setState(State s);
    State state() const { return m_state; }

    void setMarker(Marker m, int idx);
    int  marker(Marker m) const { return m_markers[m]; }

    int  binIndex()  const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

    int  ecgVisibleN() const { return m_ecgVisibleN; }

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

    int m_binIndex;
    int m_leadIndex;
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ecg;

    int m_ecgVisibleN = 0;
    int m_ppgVisibleN = 0;

    int m_markers[MarkerCount] = { -1, -1, -1, -1, -1 };

    State m_state = State::Good;
    bool  m_hasPPG = false;
    int   m_dragMarker = -1;

    static constexpr int kML = 5, kMR = 5, kMT = 16, kMB = 4;
    static constexpr int kDragPx = 12;
};