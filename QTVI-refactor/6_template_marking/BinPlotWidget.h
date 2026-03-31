// ============================================================================
// BinPlotWidget.h — Single bin plot with PPG, ECG, and draggable dicrotic line
// ============================================================================
#pragma once
#include <QWidget>
#include <vector>

class BinPlotWidget : public QWidget {
    Q_OBJECT
public:
    explicit BinPlotWidget(int binIndex, QWidget* parent = nullptr);

    void setData(const std::vector<double>& ppg,
        const std::vector<double>& ecg,
        double alignmentPoint,
        int dicroticIdx);

    void setBadR(bool bad);
    void setBadPPG(bool bad);
    void setDicrotic(int idx);
    int  dicrotic() const { return m_dicrotic; }
    bool isBadR() const { return m_badR; }
    bool isBadPPG() const { return m_badPPG; }
    int  binIndex() const { return m_binIndex; }

signals:
    void dicroticMoved(int binIndex, int newIdx);
    void badRToggled(int binIndex, bool bad);
    void badPPGToggled(int binIndex, bool bad);
    void dicroticDragStarted(int binIndex);

protected:
    void paintEvent(QPaintEvent*) override;
    void mousePressEvent(QMouseEvent*) override;
    void mouseMoveEvent(QMouseEvent*) override;
    void mouseReleaseEvent(QMouseEvent*) override;

private:
    int sampleFromX(double x) const;
    double xFromSample(int s) const;

    int m_binIndex;
    std::vector<double> m_ppg;
    std::vector<double> m_ecg;
    double m_alignPoint = 0;
    int m_dicrotic = -1;
    bool m_badR = false;
    bool m_badPPG = false;
    bool m_dragging = false;

    // Plot margins
    static constexpr int kMarginL = 5, kMarginR = 5, kMarginT = 16, kMarginB = 4;
};