// ============================================================================
// BinPlotWidget.cpp
//
// ECG and PPG are drawn at a single fixed pixels-per-sample scale
// (BinPlotWidget::kPxPerSample). The widget reports a sizeHint() that
// expands to fit whichever trace runs longer, so the full PPG tail
// (including any pre-foot / post-trough widening from the templater) is
// always visible and shares the ECG's time axis.
//
// Std band: when a per-sample std vector is available for the trace
// (covering at least the visible samples), the widget paints a
// translucent gray polygon between mean-std and mean+std underneath
// the line. The band uses the SAME lo/hi vertical range as the trace
// so it lines up with the line at every sample. Empty std => no band,
// just the line.
// ============================================================================
#include "BinPlotWidget.hpp"
#include <QPainter>
#include <QPainterPath>
#include <QMouseEvent>
#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

    // ------------------------------------------------------------------
    // Colors
    //
    // Trace colors:
    //   ECG  - dark blue
    //   PPG  - dark red
    //
    // Marker colors are grouped by trace and lightened/darkened to
    // disambiguate the individual landmarks within each group:
    //   ECG markers - shades of black / dark blue
    //   PPG markers - shades of red
    //
    // Edit these constants in one place; markerColor() does the lookup.
    // ------------------------------------------------------------------
    constexpr QColor kColorEcgTrace{ 10,  20,  90 };   // dark navy blue
    constexpr QColor kColorPpgTrace{ 130,  10,  20 };   // dark red

    // Gray std band -- semi-transparent so the trace line on top stays
    // visible at every sample. Same color for ECG and PPG bands; the
    // bands sit behind their respective traces, so the trace color
    // disambiguates which signal each band belongs to.
    constexpr QColor kColorStdBand{ 200, 200, 200, 110 };

    // ECG markers (P, Q, Tb, Te) - blacks and dark blues, darkest to lightest.
    constexpr QColor kColorEcgP{ 0,   0,   0 };   // black
    constexpr QColor kColorEcgQBegin{ 20,  20,  60 };   // very dark navy
    constexpr QColor kColorEcgTBegin{ 40,  50, 110 };   // dark navy
    constexpr QColor kColorEcgTEnd{ 70,  90, 160 };   // medium navy

    // PPG markers (On, Pk, Dc, 50, En) - shades of red, darkest to lightest.
    constexpr QColor kColorPpgOnset{ 110,   0,   0 };  // dark red
    constexpr QColor kColorPpgPeak{ 180,   0,   0 };  // red
    constexpr QColor kColorPpgDicrotic{ 220,  50,  50 };  // medium red
    constexpr QColor kColorPpg50{ 235, 100, 100 };  // light red
    constexpr QColor kColorPpgEnd{ 200,  60,  90 };  // dark pink-red

    QColor markerColor(int m) {
        switch (m) {
        case BinPlotWidget::EcgP:        return kColorEcgP;
        case BinPlotWidget::EcgQBegin:   return kColorEcgQBegin;
        case BinPlotWidget::EcgTBegin:   return kColorEcgTBegin;
        case BinPlotWidget::EcgTEnd:     return kColorEcgTEnd;
        case BinPlotWidget::PpgOnset:    return kColorPpgOnset;
        case BinPlotWidget::PpgPeak:     return kColorPpgPeak;
        case BinPlotWidget::PpgDicrotic: return kColorPpgDicrotic;
        case BinPlotWidget::Ppg50:       return kColorPpg50;
        case BinPlotWidget::PpgEnd:      return kColorPpgEnd;
        }
        return Qt::black;
    }
    const char* markerShortLabel(int m) {
        switch (m) {
        case BinPlotWidget::EcgP:        return "P";
        case BinPlotWidget::EcgQBegin:   return "Q";
        case BinPlotWidget::EcgTBegin:   return "Tb";
        case BinPlotWidget::EcgTEnd:     return "Te";
        case BinPlotWidget::PpgOnset:    return "On";
        case BinPlotWidget::PpgPeak:     return "Pk";
        case BinPlotWidget::PpgDicrotic: return "Dc";
        case BinPlotWidget::Ppg50:       return "50";
        case BinPlotWidget::PpgEnd:      return "En";
        }
        return "?";
    }

    // Compute the visible-range vertical bounds for a trace. If a matching
    // std vector is supplied, expand the range to include mean±std at
    // every visible sample so the band fits inside the drawing area
    // without clipping.
    void computeVisibleRange(const std::vector<double>& v,
        const std::vector<double>& sd,
        int visN,
        double& lo, double& hi)
    {
        lo = *std::min_element(v.begin(), v.begin() + visN);
        hi = *std::max_element(v.begin(), v.begin() + visN);
        if ((int)sd.size() < visN) return;
        for (int i = 0; i < visN; ++i) {
            lo = std::min(lo, v[i] - sd[i]);
            hi = std::max(hi, v[i] + sd[i]);
        }
    }

    // Draw the gray ±std band at a fixed pixels-per-sample scale using
    // the supplied (lo, hi) range. The caller must use the SAME range
    // for the trace draw so the band and line agree vertically.
    void drawStdBand(QPainter& p,
        const std::vector<double>& v,
        const std::vector<double>& sd,
        double startPx, int mt, int ph,
        double pxPerSample, int visN,
        double lo, double hi)
    {
        if (visN < 2 || (int)v.size() < 2) return;
        if ((int)sd.size() < visN) return;          // empty/mismatched => no band
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        QPainterPath band;
        // Top edge: mean + std, left to right.
        for (int i = 0; i < visN; ++i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - ((v[i] + sd[i]) - lo) / r * ph;
            if (i == 0) band.moveTo(x, y);
            else        band.lineTo(x, y);
        }
        // Bottom edge: mean - std, right to left, to close the polygon.
        for (int i = visN - 1; i >= 0; --i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - ((v[i] - sd[i]) - lo) / r * ph;
            band.lineTo(x, y);
        }
        band.closeSubpath();

        p.setPen(Qt::NoPen);
        p.setBrush(kColorStdBand);
        p.drawPath(band);
    }

    // Draw a trace at a fixed pixels-per-sample scale using the supplied
    // (lo, hi) vertical range. Caller-supplied range so that the band
    // and the line share an axis.
    void drawTraceFixedScale(QPainter& p, const std::vector<double>& v,
        double startPx, int mt, int ph,
        double pxPerSample, const QPen& pen, int visN,
        double lo, double hi)
    {
        if (visN < 2 || (int)v.size() < 2) return;
        const double r = (hi - lo > 1e-10) ? (hi - lo) : 1.0;

        QPainterPath path;
        for (int i = 0; i < visN; ++i) {
            double x = startPx + (double)i * pxPerSample;
            double y = mt + ph - (v[i] - lo) / r * ph;
            if (i == 0) path.moveTo(x, y);
            else        path.lineTo(x, y);
        }

        p.setPen(pen);
        p.setBrush(Qt::NoBrush);
        p.drawPath(path);
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
    const std::vector<double>& ppgStd,
    const std::vector<double>& ecg,
    const std::vector<double>& ecgStd,
    int ecgP, int qBegin, int tBegin, int tEnd,
    int ppgOnset, int ppgPeak,
    int ppgDicrotic, int ppg50, int ppgEnd,
    double rPeakSample)
{
    m_ppg = ppg;
    m_ppgStd = ppgStd;
    m_ecg = ecg;
    m_ecgStd = ecgStd;
    m_markers[EcgP] = ecgP;
    m_markers[EcgQBegin] = qBegin;
    m_markers[EcgTBegin] = tBegin;
    m_markers[EcgTEnd] = tEnd;
    m_markers[PpgOnset] = ppgOnset;
    m_markers[PpgPeak] = ppgPeak;
    m_markers[PpgDicrotic] = ppgDicrotic;
    m_markers[Ppg50] = ppg50;
    m_markers[PpgEnd] = ppgEnd;
    m_hasPPG = !ppg.empty();
    m_rPeakSample = rPeakSample;

    // Cache the per-trace visible counts so paint, hit-test, and
    // drag-clamp all see the same numbers.
    m_ecgVisibleN = computeEcgVisibleN(m_ecg, tEnd);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));

    // sizeHint depends on visible counts; notify layout.
    updateGeometry();
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

int BinPlotWidget::requiredWidth() const {
    // ECG draws from kML to kML + (visN-1)*pxPerSample.
    // PPG draws from rPeakPx to rPeakPx + (visN-1)*pxPerSample.
    // Widget must be wide enough for whichever ends further right, plus kMR.
    const double ecgRight = (m_ecgVisibleN > 1)
        ? kML + (m_ecgVisibleN - 1) * kPxPerSample
        : kML;
    const double rPeakPx = kML + m_rPeakSample * kPxPerSample;
    const double ppgRight = (m_ppgVisibleN > 1)
        ? rPeakPx + (m_ppgVisibleN - 1) * kPxPerSample
        : kML;
    return static_cast<int>(std::ceil(std::max(ecgRight, ppgRight))) + kMR;
}

int BinPlotWidget::sampleFromX(double x, bool isEcg) const {
    if (isEcg) {
        return static_cast<int>(std::round((x - kML) / kPxPerSample));
    }
    else {
        const double rPeakPx = kML + m_rPeakSample * kPxPerSample;
        return static_cast<int>(std::round((x - rPeakPx) / kPxPerSample));
    }
}

double BinPlotWidget::xFromSample(int s, bool isEcg) const {
    if (isEcg) {
        return kML + s * kPxPerSample;
    }
    else {
        const double rPeakPx = kML + m_rPeakSample * kPxPerSample;
        return rPeakPx + s * kPxPerSample;
    }
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

void BinPlotWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    const int h = height();
    const int ph = h - kMT - kMB;

    p.fillRect(rect(), Qt::white);

    p.setPen(Qt::black);
    QFont f = p.font(); f.setPointSize(8); p.setFont(f);
    p.drawText(kML, 11,
        QString("Bin %1  [%2]").arg(m_binIndex + 1).arg(m_leadLabel));

    // -------- ECG --------
    // Compute vertical range including std, draw the band first so the
    // trace line stacks on top of it.
    if (m_ecgVisibleN >= 2 && (int)m_ecg.size() >= 2) {
        double lo, hi;
        computeVisibleRange(m_ecg, m_ecgStd, m_ecgVisibleN, lo, hi);

        drawStdBand(p, m_ecg, m_ecgStd, kML, kMT, ph,
            kPxPerSample, m_ecgVisibleN, lo, hi);
        drawTraceFixedScale(p, m_ecg, kML, kMT, ph,
            kPxPerSample, QPen(kColorEcgTrace, 1.5),
            m_ecgVisibleN, lo, hi);
    }

    // -------- PPG --------
    // PPG starts at the R-peak x-position. The R-peak's sample index
    // maps to pixels via the same fixed scale, which is what gives the
    // two traces a common time axis.
    if (m_ppgVisibleN >= 2 && (int)m_ppg.size() >= 2) {
        const double rPeakPx = kML + m_rPeakSample * kPxPerSample;
        double lo, hi;
        computeVisibleRange(m_ppg, m_ppgStd, m_ppgVisibleN, lo, hi);

        drawStdBand(p, m_ppg, m_ppgStd, rPeakPx, kMT, ph,
            kPxPerSample, m_ppgVisibleN, lo, hi);
        drawTraceFixedScale(p, m_ppg, rPeakPx, kMT, ph,
            kPxPerSample, QPen(kColorPpgTrace, 1.5),
            m_ppgVisibleN, lo, hi);
    }

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
        const int w = width();
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