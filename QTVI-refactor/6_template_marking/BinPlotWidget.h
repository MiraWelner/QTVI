// ============================================================================
// BinPlotWidget.h - One ECG lead + optional PPG overlay + draggable dicrotic
//
// Right-click cycles:  Good → BadR → BadPPG → Good  (skips BadPPG if no PPG)
// Left-drag:           move dicrotic notch line
// ============================================================================
#pragma once
#include <QWidget>
#include <QString>
#include <vector>

class BinPlotWidget : public QWidget {
    Q_OBJECT
public:
    enum class State { Good, BadR, BadPPG };

    explicit BinPlotWidget(int binIndex, int leadIndex,
        const QString& leadLabel, QWidget* parent = nullptr);

    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ecg,
        int dicroticIdx);

    void setHasPPG(bool has);
    bool hasPPG() const { return m_hasPPG; }

    void setState(State s);
    State state() const { return m_state; }

    void setDicrotic(int idx);
    int  dicrotic() const { return m_dicrotic; }

    int  binIndex() const { return m_binIndex; }
    int  leadIndex() const { return m_leadIndex; }

signals:
    void dicroticMoved(int binIndex, int newIdx);
    void dicroticDragStarted(int binIndex);
    void badRToggled(int binIndex, int leadIndex, bool bad);
    void badPPGToggled(int binIndex, bool bad);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    int sampleFromX(double x, int nSamples) const;
    double xFromSample(int s, int nSamples) const;

    int m_binIndex;
    int m_leadIndex;
    QString m_leadLabel;
    std::vector<double> m_ppg;
    std::vector<double> m_ecg;
    int m_dicrotic = -1;
    State m_state = State::Good;
    bool m_hasPPG = false;
    bool m_dragging = false;

    static constexpr int kML = 5, kMR = 5, kMT = 16, kMB = 4;
};