// ============================================================================
// BinPlotWidget.cpp
// ============================================================================
#include "BinPlotWidget.hpp"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace {
    QColor markerColor(int m) {
        switch (m) {
            // ECG markers — shades of green
        case BinPlotWidget::EcgQBegin: return QColor(20, 130, 40);
        case BinPlotWidget::EcgTBegin: return QColor(80, 180, 80);
        case BinPlotWidget::EcgTEnd:   return QColor(140, 200, 60);
            // PPG markers — shades of blue
        case BinPlotWidget::PpgOnset:  return QColor(3, 0, 255);
        case BinPlotWidget::PpgPeak:   return QColor(1, 0, 120);
        }
        return Qt::black;
    }
    const char* markerShortLabel(int m) {
        switch (m) {
        case BinPlotWidget::EcgQBegin: return "Q";
        case BinPlotWidget::EcgTBegin: return "Tb";
        case BinPlotWidget::EcgTEnd:   return "Te";
        case BinPlotWidget::PpgOnset:  return "On";
        case BinPlotWidget::PpgPeak:   return "Pk";
        }
        return "?";
    }
}

BinPlotWidget::BinPlotWidget(int binIndex, int leadIndex,
    const QString& leadLabel, QWidget* parent)
    : QWidget(parent), m_binIndex(binIndex), m_leadIndex(leadIndex),
    m_leadLabel(leadLabel)
{
    setMinimumHeight(100);
    setMouseTracking(true);
    setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void BinPlotWidget::setData(const std::vector<double>& ppg,
    const std::vector<double>& ecg,
    int qBegin, int tBegin, int tEnd,
    int ppgOnset, int ppgPeak)
{
    m_ppg = ppg;
    m_ecg = ecg;
    m_markers[EcgQBegin] = qBegin;
    m_markers[EcgTBegin] = tBegin;
    m_markers[EcgTEnd] = tEnd;
    m_markers[PpgOnset] = ppgOnset;
    m_markers[PpgPeak] = ppgPeak;
    m_hasPPG = !ppg.empty();

    // Cache the per-trace visible counts so paint, hit-test, and
    // drag-clamp all see the same numbers.
    m_ecgVisibleN = computeEcgVisibleN(m_ecg, tEnd);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));

    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setState(State s) { m_state = s; update(); }

void BinPlotWidget::setMarker(Marker m, int idx) {
    m_markers[m] = idx;
    update();
}

int BinPlotWidget::visibleN(bool isEcg) const {
    return isEcg ? m_ecgVisibleN : m_ppgVisibleN;
}

int BinPlotWidget::sampleFromX(double x, bool isEcg) const {
    const int pw = width() - kML - kMR;
    const int n = visibleN(isEcg);
    if (n < 2 || pw < 1) return 0;
    return static_cast<int>(std::round((x - kML) * (n - 1.0) / pw));
}

double BinPlotWidget::xFromSample(int s, bool isEcg) const {
    const int pw = width() - kML - kMR;
    const int n = visibleN(isEcg);
    if (n < 2) return kML;
    return kML + s * pw / static_cast<double>(n - 1);
}

int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = kDragPx + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        const bool isEcg = markerIsEcg(m);
        const auto& vec = isEcg ? m_ecg : m_ppg;
        if (vec.empty() || idx >= (int)vec.size()) continue;
        if (!isEcg && !m_hasPPG) continue;
        if (idx >= visibleN(isEcg)) continue;
        const double mx = xFromSample(idx, isEcg);
        const double d = std::abs(x - mx);
        if (d < bestDist) { bestDist = d; best = m; }
    }
    return best;
}

namespace {
    void drawTrace(QPainter& p, const std::vector<double>& v,
        int ml, int mt, int pw, int ph, const QPen& pen, int visN) {
        if (visN < 2 || (int)v.size() < 2) return;

        double lo = *std::min_element(v.begin(), v.begin() + visN);
        double hi = *std::max_element(v.begin(), v.begin() + visN);
        double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        QPainterPath path;
        for (int i = 0; i < visN; ++i) {
            double x = ml + (double)i * pw / (visN - 1);
            double y = mt + ph - (v[i] - lo) / r * ph;
            if (i == 0) path.moveTo(x, y);
            else        path.lineTo(x, y);
        }

        p.setPen(pen);
        p.drawPath(path);
    }
}

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    int w = width(), h = height();
    int pw = w - kML - kMR;
    int ph = h - kMT - kMB;

    p.fillRect(rect(), Qt::white);

    p.setPen(Qt::black);
    QFont f = p.font(); f.setPointSize(8); p.setFont(f);
    p.drawText(kML, 11,
        QString("Bin %1  [%2]").arg(m_binIndex + 1).arg(m_leadLabel));
    
    drawTrace(p, m_ecg, kML, kMT, pw, ph,
        QPen(QColor(0, 207, 34), 1.5), m_ecgVisibleN);
    drawTrace(p, m_ppg, kML, kMT, pw, ph,
        QPen(QColor(0, 0, 82), 1.5), m_ppgVisibleN);

    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        bool isEcg = markerIsEcg(m);
        const auto& vec = isEcg ? m_ecg : m_ppg;
        if (vec.empty() || idx >= (int)vec.size()) continue;
        if (!isEcg && !m_hasPPG) continue;
        if (idx >= visibleN(isEcg)) continue;
        double mx = xFromSample(idx, isEcg);
        QPen pen(markerColor(m), 2);
        p.setPen(pen);
        p.drawLine(QPointF(mx, kMT), QPointF(mx, h - kMB));
        p.drawText(QPointF(mx + 2, kMT + 8), markerShortLabel(m));
    }

    if (m_state == State::BadPPG) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(kML, kMT, w - kMR, h - kMB);
        p.drawLine(kML, h - kMB, w - kMR, kMT);
    }

    if (m_state == State::BadR) {
        p.setPen(QPen(QColor(200, 0, 0), 2));
        QFont bf = p.font(); bf.setPointSize(14); bf.setBold(true); p.setFont(bf);
        p.drawText(rect(), Qt::AlignCenter, "BAD R");
    }
}

void BinPlotWidget::mousePressEvent(QMouseEvent* e) {
    if (e->button() == Qt::LeftButton) {
        int m = markerAtX(e->position().x());
        if (m >= 0) {
            m_dragMarker = m;
            emit markerDragStarted(m_binIndex, m_leadIndex, m);
            return;
        }
    }

    if (e->button() == Qt::RightButton) {
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
    if (m_dragMarker < 0) return;
    const bool isEcg = markerIsEcg(m_dragMarker);
    const auto& vec = isEcg ? m_ecg : m_ppg;
    if (vec.empty()) return;
    int s = sampleFromX(e->position().x(), isEcg);
    s = std::clamp(s, 0, visibleN(isEcg) - 1);
    m_markers[m_dragMarker] = s;
    emit markerMoved(m_binIndex, m_leadIndex, m_dragMarker, s);
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragMarker = -1;
}