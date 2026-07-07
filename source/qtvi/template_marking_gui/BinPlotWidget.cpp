// ============================================================================
// BinPlotWidget.cpp
//
// ECG and PPG are drawn at a single per-paint pixels-per-sample scale
// (BinPlotWidget::pxPerSample()), chosen so the widest trace fills the
// cell the layout gives this widget. Every bin window therefore comes out
// the same on-screen length; the full PPG tail stays visible and shares
// the ECG's time axis within the widget.
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

    //different colored STD bands - yellowish for ecg and gray for ppg
    inline const QColor color_stdband_ecg{ 255, 215, 0, 110 };
    inline const QColor color_stdband_ppg{ 200, 200, 200, 110 };


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
    // std vector is supplied, expand the range to include mean +- std at
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
    void drawStdBand(QPainter& p, const std::vector<double>& v, const std::vector<double>& sd, double startPx, int mt, int ph,
        double pxPerSample, int visN, double lo, double hi, QColor color)
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
        p.setBrush(color);
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
    m_rPeakSample = rPeakSample;
    m_hasPPG = !ppg.empty();
    // Cache the per-trace visible counts so paint, hit-test, and
    // drag-clamp all see the same numbers.
    m_ecgVisibleN = computeEcgVisibleN(m_ecg, tEnd);
    m_ppgVisibleN = visiblePpgCount(static_cast<int>(m_ppg.size()));

    // sizeHint depends on visible counts; notify layout.
    updateGeometry();
    update();
}

void BinPlotWidget::setHasPPG(bool has) { m_hasPPG = has; }
void BinPlotWidget::setSampleRate(double hz) {
    if (m_sampleRate == hz) return;
    m_sampleRate = hz;
    update();
}
void BinPlotWidget::setState(State s) { m_state = s; update(); }

void BinPlotWidget::setShowEcgMarkers(bool show) {
    if (m_showEcgMarkers == show) return;
    m_showEcgMarkers = show;
    update();
}

void BinPlotWidget::setShowPpgMarkers(bool show) {
    if (m_showPpgMarkers == show) return;
    m_showPpgMarkers = show;
    update();
}

void BinPlotWidget::setMarker(Marker m, int idx) {
    m_markers[m] = idx;
    update();
}

int BinPlotWidget::visibleN(bool isEcg) const {
    return isEcg ? m_ecgVisibleN : m_ppgVisibleN;
}

int BinPlotWidget::totalSampleSpan() const {
    // ECG runs over samples [0, m_ecgVisibleN). The PPG is drawn shifted
    // right by the R-peak offset, so it runs over [0, rPeak + ppgVisibleN).
    // The widest of the two is what has to fit in the drawable width.
    const int pad = static_cast<int>(std::round(m_rPeakSample));
    const int ppgSpan = (m_ppgVisibleN > 0) ? pad + m_ppgVisibleN : 0;
    return std::max(std::max(m_ecgVisibleN, ppgSpan), 2);
}

double BinPlotWidget::pxPerSample() const {
     // Scale so the widest trace exactly fills the drawable width. ECG and
    // PPG share this one value, so they stay aligned within this widget;
    // different bins use different values, which is what makes every bin
    // window come out the same on-screen length.
    const int span = totalSampleSpan();
    const double drawW = std::max(1, width() - margin_left - margin_right);
    return (span > 1) ? (drawW / static_cast<double>(span - 1)) : 1.0;
}

int BinPlotWidget::sampleFromX(double x, bool isEcg) const {
    const double pps = pxPerSample();
    if (isEcg) {
        return static_cast<int>(std::round((x - margin_left) / pps));
    }
    else {
        const double rPeakPx = margin_left + m_rPeakSample * pps;
        return static_cast<int>(std::round((x - rPeakPx) / pps));
    }
}

double BinPlotWidget::xFromSample(int s, bool isEcg) const {
    const double pps = pxPerSample();
    if (isEcg) {
        return margin_left + s * pps;
    }
    else {
        const double rPeakPx = margin_left + m_rPeakSample * pps;
        return rPeakPx + s * pps;
    }
}

int BinPlotWidget::markerAtX(double x) const {
    int best = -1;
    double bestDist = click_radius_around_marker + 1.0;
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        const bool isEcg = markerIsEcg(m);
        if (isEcg ? !m_showEcgMarkers : !m_showPpgMarkers) continue;
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
    const int w = width();
    const int ph = h - margin_top - margin_bottom;
    const double pps = pxPerSample();   // one scale for both traces this paint

    p.fillRect(rect(), Qt::white);

    // Title
    p.setPen(Qt::black);
    { QFont f = p.font(); f.setPointSize(8); p.setFont(f); }
    p.drawText(margin_left, 11, QString("Bin %1  [%2]").arg(m_binIndex).arg(m_leadLabel));

    double yLo = 0, yHi = 0;
    computeVisibleRange(m_ecg, m_ecgStd, m_ecgVisibleN, yLo, yHi);

    // ---- Axes: frame + ticks + labels ----
    {
        const double xAxisY = margin_top + ph;       // == h - kMB
        const double yAxisX = margin_left;

        p.setPen(QColor(150, 150, 150));
        p.drawLine(QPointF(yAxisX, margin_top), QPointF(yAxisX, xAxisY));   // y-axis
        p.drawLine(QPointF(yAxisX, xAxisY), QPointF(w - margin_right, xAxisY));  // x-axis

        QFont af = p.font(); af.setPointSize(7); p.setFont(af);

        // X ticks: sample index across the full widget span.
        const int span = totalSampleSpan();
        const int nx = 4;
        for (int t = 0; t <= nx; ++t) {
            const int s = static_cast<int>(std::round((double)span * t / nx));
            double x = margin_left + s * pps;
            if (x > w - margin_right) x = w - margin_right;
            p.drawLine(QPointF(x, xAxisY), QPointF(x, xAxisY + 3));
            QString lbl;
            if (m_sampleRate > 0.0) {
                lbl = QString::number(s / m_sampleRate, 'f', 2);
                if (t == nx) lbl += " s";   // unit on the last tick only
            }
            else {
                lbl = QString::number(s);
            }
            p.drawText(QPointF(x - 3.0 * lbl.size(), xAxisY + 12), lbl);
        }

        // Y ticks: amplitude. Top of plot = hi, bottom = lo.
        for (int t = 0; t < 2; ++t) {
            const double y = (t == 0) ? margin_top : margin_top + ph;
            const double val = (t == 0) ? yHi : yLo;
            p.drawLine(QPointF(yAxisX - 3, y), QPointF(yAxisX, y));
            p.drawText(QPointF(1, y + 3), QString::number(val, 'e', 1));
        }
        
    }

    // -------- ECG --------
    drawStdBand(p, m_ecg, m_ecgStd, margin_left, margin_top, ph, pps, m_ecgVisibleN, yLo, yHi, color_stdband_ecg);
    drawTraceFixedScale(p, m_ecg, margin_left, margin_top, ph, pps, QPen(kColorEcgTrace, 1.5), m_ecgVisibleN, yLo, yHi);

    // -------- PPG --------
    if (m_hasPPG && !m_ppg.empty() && m_ppgVisibleN > 0) {
        const int pad = static_cast<int>(std::round(m_rPeakSample));
        const int visN = pad + m_ppgVisibleN;

        std::vector<double> ppgPadded;
        ppgPadded.reserve(visN);
        ppgPadded.insert(ppgPadded.end(), pad, m_ppg[0]);
        ppgPadded.insert(ppgPadded.end(), m_ppg.begin(),
            m_ppg.begin() + m_ppgVisibleN);

        std::vector<double> stdPadded;
        if ((int)m_ppgStd.size() >= m_ppgVisibleN) {
            stdPadded.reserve(visN);
            stdPadded.insert(stdPadded.end(), pad, 0.0);
            stdPadded.insert(stdPadded.end(), m_ppgStd.begin(),
                m_ppgStd.begin() + m_ppgVisibleN);
        }

        double lo, hi;
        computeVisibleRange(ppgPadded, stdPadded, visN, lo, hi);

        drawStdBand(p, ppgPadded, stdPadded, margin_left, margin_top, ph, pps, visN, lo, hi, color_stdband_ppg);
        drawTraceFixedScale(p, ppgPadded, margin_left, margin_top, ph,
            pps, QPen(kColorPpgTrace, 1.5), visN, lo, hi);
    }
    

    QFont smallF = p.font(); smallF.setPointSize(7); p.setFont(smallF);
    for (int m = 0; m < MarkerCount; ++m) {
        int idx = m_markers[m];
        if (idx < 0) continue;
        bool isEcg = markerIsEcg(m);
        if (isEcg ? !m_showEcgMarkers : !m_showPpgMarkers) continue;
        const auto& vec = isEcg ? m_ecg : m_ppg;
        if (vec.empty() || idx >= (int)vec.size()) continue;
        if (!isEcg && !m_hasPPG) continue;
        if (idx >= visibleN(isEcg)) continue;
        double mx = xFromSample(idx, isEcg);
        QPen pen(markerColor(m), 2);
        p.setPen(pen);
        p.drawLine(QPointF(mx, margin_top), QPointF(mx, h - margin_bottom));
        p.drawText(QPointF(mx + 2, margin_top + 8), markerShortLabel(m));
    }

    if (m_state == State::BadPPG) {
        p.setPen(QPen(Qt::red, 4));
        p.drawLine(margin_left, margin_top, w - margin_right, h - margin_bottom);
        p.drawLine(margin_left, h - margin_bottom, w - margin_right, margin_top);
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