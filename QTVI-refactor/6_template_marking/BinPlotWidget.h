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
#include <QString>
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
    int    sampleFromX(double x, int nSamples, bool isEcg) const;
    double xFromSample(int s, int nSamples, bool isEcg) const;
    int    markerAtX(double x) const;

    int m_binIndex;
    int m_leadIndex;
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ecg;

    int m_markers[MarkerCount] = { -1, -1, -1, -1, -1 };

    State m_state = State::Good;
    bool  m_hasPPG = false;
    int   m_dragMarker = -1;

    static constexpr int kML = 5, kMR = 5, kMT = 16, kMB = 4;
    static constexpr int kDragPx = 12;
};