// ============================================================================
// BinPlotWidget.cpp
// ============================================================================
#include "BinPlotWidget.h"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <algorithm>
#include <cmath>

namespace {
    QColor markerColor(int m) {
        switch (m) {
            // ECG markers — shades of green
        case BinPlotWidget::EcgQBegin: return QColor(20, 130, 40);   // deep green
        case BinPlotWidget::EcgTBegin: return QColor(80, 180, 80);   // mid green
        case BinPlotWidget::EcgTEnd:   return QColor(140, 200, 60);  // yellow-green
            // PPG markers — shades of brown
        case BinPlotWidget::PpgOnset:  return QColor(120, 70, 30);   // dark brown
        case BinPlotWidget::PpgPeak:   return QColor(180, 110, 50);  // light brown
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
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setState(State s) { m_state = s; update(); }

void BinPlotWidget::setMarker(Marker m, int idx) {
    m_markers[m] = idx;
    update();
}

namespace {
    // How many samples of this trace are actually drawn. ECG is cropped to
    // the first half (the back half is the run-up to the next R-peak and
    // just confuses fiducial marking); PPG is shown in full because the
    // dicrotic notch and other features we care about live in the back half.
    inline int visibleCount(bool isEcg, int nFull) {
        if (nFull < 2) return 2;
        return isEcg ? std::max(nFull / 2, 2) : nFull;
    }
}

int BinPlotWidget::sampleFromX(double x, int nSamples, bool isEcg) const {
    int pw = width() - kML - kMR;
    int n = visibleCount(isEcg, nSamples);
    return static_cast<int>(std::round((x - kML) * (n - 1.0) / pw));
}

double BinPlotWidget::xFromSample(int s, int nSamples, bool isEcg) const {
    int pw = width() - kML - kMR;
    int n = visibleCount(isEcg, nSamples);
    return kML + s * pw / (double)(n - 1);
}

int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = kDragPx + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        bool isEcg = markerIsEcg(m);
        const auto& vec = isEcg ? m_ecg : m_ppg;
        if (vec.empty() || idx >= (int)vec.size()) continue;
        if (!isEcg && !m_hasPPG) continue;
        // Skip markers that fall outside the visible portion of their trace.
        const int visN = visibleCount(isEcg, (int)vec.size());
        if (idx >= visN) continue;
        double mx = xFromSample(idx, (int)vec.size(), isEcg);
        double d = std::abs(x - mx);
        if (d < bestDist) { bestDist = d; best = m; }
    }
    return best;
}

namespace {
    void drawTrace(QPainter& p, const std::vector<double>& v,
        int ml, int mt, int pw, int ph, const QPen& pen, bool isEcg) {
        int nFull = (int)v.size();
        if (nFull < 2) return;

        // ECG: only the first half is drawn (zoomed 2x across the plot).
        // PPG: the full trace is drawn — the back half holds features
        // (e.g. dicrotic notch) we want to be able to mark.
        int n = visibleCount(isEcg, nFull);

        double lo = *std::min_element(v.begin(), v.begin() + n);
        double hi = *std::max_element(v.begin(), v.begin() + n);
        double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        // Build the whole trace as one path so any pen pattern flows
        // continuously across the trace instead of resetting per segment.
        QPainterPath path;
        for (int i = 0; i < n; ++i) {
            double x = ml + (double)i * pw / (n - 1);
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
        QString("Bin %1 - %2").arg(m_binIndex).arg(m_leadLabel));

    drawTrace(p, m_ecg, kML, kMT, pw, ph,
        QPen(QColor(40, 40, 40), 1.5), /*isEcg=*/true);
    // PPG trace: solid brown, full length (not cropped to first half).
    drawTrace(p, m_ppg, kML, kMT, pw, ph,
        QPen(QColor(140, 80, 30), 1.5), /*isEcg=*/false);

    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        bool isEcg = markerIsEcg(m);
        const auto& vec = isEcg ? m_ecg : m_ppg;
        if (vec.empty() || idx >= (int)vec.size()) continue;
        if (!isEcg && !m_hasPPG) continue;
        // Skip markers that fall outside their trace's visible range.
        const int visN = visibleCount(isEcg, (int)vec.size());
        if (idx >= visN) continue;
        double mx = xFromSample(idx, (int)vec.size(), isEcg);
        QPen pen(markerColor(m), 2);
        // Solid lines for every marker (no dashed styles).
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
    bool isEcg = markerIsEcg(m_dragMarker);
    const auto& vec = isEcg ? m_ecg : m_ppg;
    if (vec.empty()) return;
    int s = sampleFromX(e->position().x(), (int)vec.size(), isEcg);
    // Clamp to the visible range of THIS trace — first half for ECG,
    // full length for PPG.
    const int visN = visibleCount(isEcg, (int)vec.size());
    s = std::clamp(s, 0, visN - 1);
    m_markers[m_dragMarker] = s;
    emit markerMoved(m_binIndex, m_leadIndex, m_dragMarker, s);
    update();
}

void BinPlotWidget::mouseReleaseEvent(QMouseEvent*) {
    m_dragMarker = -1;
}