// ============================================================================
// BinPlotWidget.cpp
// ============================================================================
#include "BinPlotWidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <algorithm>
#include <cmath>

BinPlotWidget::BinPlotWidget(int binIndex, int leadIndex,
    const QString& leadLabel, QWidget* parent)
    : QWidget(parent), m_binIndex(binIndex), m_leadIndex(leadIndex),
    m_leadLabel(leadLabel)
{
    setMinimumHeight(60);
    setMouseTracking(true);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ecg, int dicroticIdx) {
    m_ppg = ppg;
    m_ecg = ecg;
    m_dicrotic = dicroticIdx;
    m_hasPPG = !ppg.empty();
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setState(State s) { m_state = s; update(); }
void BinPlotWidget::setDicrotic(int idx) { m_dicrotic = idx; update(); }

int BinPlotWidget::sampleFromX(double x, int nSamples) const {
    int pw = width() - kML - kMR;
    int n = std::max(nSamples, 2);
    return static_cast<int>(std::round((x - kML) * (n - 1.0) / pw));
}

double BinPlotWidget::xFromSample(int s, int nSamples) const {
    int pw = width() - kML - kMR;
    int n = std::max(nSamples, 2);
    return kML + s * pw / (double)(n - 1);
}

namespace {
    void drawTrace(QPainter& p, const std::vector<double>& v,
        int ml, int mt, int pw, int ph, const QPen& pen) {
        int n = (int)v.size();
        if (n < 2) return;
        double lo = *std::min_element(v.begin(), v.end());
        double hi = *std::max_element(v.begin(), v.end());
        double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;
        p.setPen(pen);
        for (int i = 1; i < n; ++i) {
            double x0 = ml + (i - 1.0) * pw / (n - 1);
            double x1 = ml + (double)i * pw / (n - 1);
            double y0 = mt + ph - (v[i - 1] - lo) / r * ph;
            double y1 = mt + ph - (v[i] - lo) / r * ph;
            p.drawLine(QPointF(x0, y0), QPointF(x1, y1));
        }
    }
}

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    int w = width(), h = height();
    int pw = w - kML - kMR;
    int ph = h - kMT - kMB;

    p.fillRect(rect(), Qt::white);

    // Label
    p.setPen(Qt::black);
    QFont f = p.font(); f.setPointSize(8); p.setFont(f);
    p.drawText(kML, 11,
        QString("Bin %1 - %2").arg(m_binIndex).arg(m_leadLabel));

    // ECG in gray
    drawTrace(p, m_ecg, kML, kMT, pw, ph, QPen(QColor(180, 180, 180), 1));

    // PPG in red
    drawTrace(p, m_ppg, kML, kMT, pw, ph, QPen(Qt::red, 1.5));

    // Dicrotic notch line (blue)
    if (!m_ppg.empty() && m_dicrotic >= 0 && m_dicrotic < (int)m_ppg.size()) {
        double dx = xFromSample(m_dicrotic, (int)m_ppg.size());
        p.setPen(QPen(Qt::blue, 2));
        p.drawLine(QPointF(dx, kMT), QPointF(dx, h - kMB));
    }

    // Bad PPG overlay (red X)
    if (m_state == State::BadPPG) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(kML, kMT, w - kMR, h - kMB);
        p.drawLine(kML, h - kMB, w - kMR, kMT);
    }

    // Bad R overlay (text)
    if (m_state == State::BadR) {
        p.setPen(QPen(QColor(200, 0, 0), 2));
        QFont bf = p.font(); bf.setPointSize(14); bf.setBold(true); p.setFont(bf);
        p.drawText(rect(), Qt::AlignCenter, "BAD R");
    }
}

void BinPlotWidget::mousePressEvent(QMouseEvent* e) {
    // Left click: start dragging dicrotic line
    if (e->button() == Qt::LeftButton && m_hasPPG
        && m_dicrotic >= 0 && !m_ppg.empty()) {
        double dx = xFromSample(m_dicrotic, (int)m_ppg.size());
        if (std::abs(e->position().x() - dx) < 15) {
            m_dragging = true;
            emit dicroticDragStarted(m_binIndex);
            return;
        }
    }

    // Right click: cycle Good -> BadR -> BadPPG -> Good
    //              (skip BadPPG if no PPG)
    if (e->button() == Qt::RightButton) {
        State prev = m_state;
        switch (m_state) {
        case State::Good:
            m_state = State::BadR;
            emit badRToggled(m_binIndex, m_leadIndex, true);
            break;
        case State::BadR:
            if (m_hasPPG) {
                m_state = State::BadPPG;
                emit badRToggled(m_binIndex, m_leadIndex, false);
                emit badPPGToggled(m_binIndex, true);
            }
            else {
                m_state = State::Good;
                emit badRToggled(m_binIndex, m_leadIndex, false);
            }
            break;
        case State::BadPPG:
            m_state = State::Good;
            emit badPPGToggled(m_binIndex, false);
            break;
        }
        update();
        return;
    }
}

void BinPlotWidget::mouseMoveEvent(QMouseEvent* e) {
    if (m_dragging && !m_ppg.empty()) {
        int s = sampleFromX(e->position().x(), (int)m_ppg.size());
        s = std::clamp(s, 0, (int)m_ppg.size() - 1);
        m_dicrotic = s;
        emit dicroticMoved(m_binIndex, s);
        update();
    }
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragging = false;
}