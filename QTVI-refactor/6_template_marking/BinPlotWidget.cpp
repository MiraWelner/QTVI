// ============================================================================
// BinPlotWidget.cpp
// ============================================================================
#include "BinPlotWidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <algorithm>
#include <cmath>

BinPlotWidget::BinPlotWidget(int binIndex, QWidget* parent)
    : QWidget(parent), m_binIndex(binIndex) {
    setMinimumHeight(80);
    setMouseTracking(true);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
}

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ecg,
    double alignmentPoint, int dicroticIdx) {
    m_ppg = ppg;
    m_ecg = ecg;
    m_alignPoint = alignmentPoint;
    m_dicrotic = dicroticIdx;
    update();
}

void BinPlotWidget::setBadR(bool bad) { m_badR = bad; update(); }
void BinPlotWidget::setBadPPG(bool bad) { m_badPPG = bad; update(); }
void BinPlotWidget::setDicrotic(int idx) { m_dicrotic = idx; update(); }

int BinPlotWidget::sampleFromX(double x) const {
    int plotW = width() - kMarginL - kMarginR;
    int nSamples = std::max((int)m_ppg.size(), 1);
    return static_cast<int>(std::round((x - kMarginL) * (nSamples - 1.0) / plotW));
}

double BinPlotWidget::xFromSample(int s) const {
    int plotW = width() - kMarginL - kMarginR;
    int nSamples = std::max((int)m_ppg.size(), 1);
    return kMarginL + s * plotW / (double)(nSamples - 1);
}

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    int w = width(), h = height();
    int plotW = w - kMarginL - kMarginR;
    int plotH = h - kMarginT - kMarginB;

    // Background
    p.fillRect(rect(), Qt::white);

    // Label
    p.setPen(Qt::black);
    QFont f = p.font(); f.setPointSize(8); p.setFont(f);
    p.drawText(kMarginL, 11, QString("Bin %1").arg(m_binIndex));

    if (m_ppg.empty() && m_ecg.empty()) return;

    // Draw ECG (gray, right y-axis scale)
    if (!m_ecg.empty()) {
        double eMin = *std::min_element(m_ecg.begin(), m_ecg.end());
        double eMax = *std::max_element(m_ecg.begin(), m_ecg.end());
        double eRange = (eMax - eMin > 1e-10) ? (eMax - eMin) : 1.0;

        p.setPen(QPen(QColor(180, 180, 180), 1));
        int n = (int)m_ecg.size();
        for (int i = 1; i < n; ++i) {
            double x0 = kMarginL + (i - 1.0) * plotW / (n - 1);
            double x1 = kMarginL + (double)i * plotW / (n - 1);
            double y0 = kMarginT + plotH - (m_ecg[i - 1] - eMin) / eRange * plotH;
            double y1 = kMarginT + plotH - (m_ecg[i] - eMin) / eRange * plotH;
            p.drawLine(QPointF(x0, y0), QPointF(x1, y1));
        }
    }

    // Draw PPG (red, left y-axis scale)
    if (!m_ppg.empty()) {
        double pMin = *std::min_element(m_ppg.begin(), m_ppg.end());
        double pMax = *std::max_element(m_ppg.begin(), m_ppg.end());
        double pRange = (pMax - pMin > 1e-10) ? (pMax - pMin) : 1.0;

        p.setPen(QPen(Qt::red, 1.5));
        int n = (int)m_ppg.size();
        for (int i = 1; i < n; ++i) {
            double x0 = kMarginL + (i - 1.0) * plotW / (n - 1);
            double x1 = kMarginL + (double)i * plotW / (n - 1);
            double y0 = kMarginT + plotH - (m_ppg[i - 1] - pMin) / pRange * plotH;
            double y1 = kMarginT + plotH - (m_ppg[i] - pMin) / pRange * plotH;
            p.drawLine(QPointF(x0, y0), QPointF(x1, y1));
        }

        // Dicrotic notch line (blue)
        if (m_dicrotic >= 0 && m_dicrotic < n) {
            double dx = xFromSample(m_dicrotic);
            p.setPen(QPen(Qt::blue, 2));
            p.drawLine(QPointF(dx, kMarginT), QPointF(dx, h - kMarginB));
        }
    }

    // Bad overlays
    if (m_badPPG) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(kMarginL, kMarginT, w - kMarginR, h - kMarginB);
        p.drawLine(kMarginL, h - kMarginB, w - kMarginR, kMarginT);
    }
    if (m_badR && !m_badPPG) {
        p.setPen(QPen(QColor(200, 0, 0), 2));
        QFont bf = p.font(); bf.setPointSize(14); bf.setBold(true); p.setFont(bf);
        p.drawText(rect(), Qt::AlignCenter, "BAD R");
    }
}

void BinPlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton && m_dicrotic >= 0 && !m_ppg.empty()) {
        double dx = xFromSample(m_dicrotic);
        if (std::abs(e->position().x() - dx) < 15) {
            m_dragging = true;
            emit dicroticDragStarted(m_binIndex);
            return;
        }
    }
    if (e->button() == Qt::MiddleButton) {
        m_badR = !m_badR;
        emit badRToggled(m_binIndex, m_badR);
        update();
        return;
    }
    if (e->button() == Qt::RightButton) {
        m_badPPG = !m_badPPG;
        m_badR = m_badPPG ? m_badR : m_badR; // PPG bad doesn't auto-clear R
        emit badPPGToggled(m_binIndex, m_badPPG);
        update();
        return;
    }
}

void BinPlotWidget::mouseMoveEvent(QMouseEvent* e) {
    if (m_dragging) {
        int s = sampleFromX(e->position().x());
        s = std::clamp(s, 0, (int)m_ppg.size() - 1);
        m_dicrotic = s;
        emit dicroticMoved(m_binIndex, s);
        update();
    }
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragging = false;
}