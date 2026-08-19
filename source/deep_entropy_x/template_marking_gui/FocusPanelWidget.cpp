#include "FocusPanelWidget.hpp"
#include "anchor_fit.hpp"

#include <QPainter>
#include <QPainterPath>
#include <QPaintEvent>
#include <algorithm>
#include <cmath>
#include <limits>

FocusPanelWidget::FocusPanelWidget(QWidget* parent) : QWidget(parent) {
    setMinimumHeight(180);
    setMinimumWidth(240);
}

void FocusPanelWidget::setFocus(const std::vector<double>& mean,
    const std::vector<double>& sd,
    int nBeats,
    int landmarkCol,
    const QString& label,
    int halfWindowSamples,
    int framingBias)
{
    m_mean = mean;
    m_sd = sd;
    m_nBeats = nBeats;
    m_landmarkCol = landmarkCol;
    m_label = label;
    m_half = std::max(4, halfWindowSamples);
    m_framingBias = framingBias;
    m_active = (landmarkCol >= 0 && !mean.empty());
    update();
}

void FocusPanelWidget::clearFocus() {
    m_active = false;
    m_mean.clear();
    m_sd.clear();
    m_nBeats = 0;
    m_landmarkCol = -1;
    update();
}

// Fitted curve over [lo, hi] via anchor_fit's BIC model selection, sampled
// per column. NaN outside [lo, hi].
std::vector<double> FocusPanelWidget::fittedCurve(int lo, int hi) const {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> out(m_mean.size(), NaN);
    if (lo < 0 || hi >= (int)m_mean.size() || hi - lo < 3) return out;
    // Fit against the mean trace within the zoom window (the anchored
    // average is what the operator sees and edits against).
    const auto fit = anchor_fit::selectAnchorModel(m_mean, lo, hi);
    if (!fit.f) return out;
    for (int i = lo; i <= hi; ++i) out[i] = fit.f(static_cast<double>(i));
    return out;
}

void FocusPanelWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.fillRect(rect(), QColor(250, 250, 250));

    const int mt = 24;                 // top margin (header)
    const int mb = 12, ml = 8, mr = 8;
    const int ph = height() - mt - mb;
    const int pw = width() - ml - mr;

    // Header label.
    p.setPen(QColor(60, 60, 60));
    p.drawText(QRect(ml, 4, width() - ml - mr, 18),
        Qt::AlignLeft | Qt::AlignVCenter,
        m_active ? m_label : QStringLiteral("Focus: (no landmark selected)"));

    if (!m_active || ph < 20 || pw < 20) return;

    // Zoom window [lo, hi] around the landmark, clamped to the trace.
    const int N = static_cast<int>(m_mean.size());
    // Framing bias shifts the window so the landmark sits toward an edge:
    // bias -1 -> landmark near the RIGHT edge (it ends this segment), so the
    // window extends mostly to the LEFT; bias +1 -> landmark near the LEFT
    // edge (starts this segment), window extends mostly RIGHT; 0 = centered.
    const int biasShift = (m_framingBias == 0) ? 0
        : (m_framingBias < 0 ? +(m_half * 3 / 4) : -(m_half * 3 / 4));
    const int center = m_landmarkCol + biasShift;
    const int lo = std::max(0, center - m_half);
    const int hi = std::min(N - 1, center + m_half);
    const int visN = hi - lo + 1;
    if (visN < 2) return;

    // 95% CI half-width per column: 1.96 * se, se = sd / sqrt(nBeats).
    // (Spec B2.) Guard nBeats >= 1; if a column's sd is NaN or n < 1, the
    // band is suppressed there (NaN), matching drawIqrBand's NaN handling.
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double sqrtN = (m_nBeats > 0) ? std::sqrt((double)m_nBeats) : 0.0;
    std::vector<double> ci(N, NaN);
    const bool haveSd = ((int)m_sd.size() == N) && sqrtN > 0.0;
    if (haveSd) {
        for (int i = lo; i <= hi; ++i) {
            if (std::isnan(m_sd[i])) continue;
            ci[i] = 1.96 * (m_sd[i] / sqrtN);
        }
    }

    // Vertical range: cover mean +/- CI across the window (fall back to the
    // mean's own range if the band is absent), with a little padding.
    double vlo = std::numeric_limits<double>::infinity();
    double vhi = -std::numeric_limits<double>::infinity();
    for (int i = lo; i <= hi; ++i) {
        if (std::isnan(m_mean[i])) continue;
        const double band = std::isnan(ci[i]) ? 0.0 : ci[i];
        vlo = std::min(vlo, m_mean[i] - band);
        vhi = std::max(vhi, m_mean[i] + band);
    }
    if (!(vhi > vlo)) return;
    const double padv = 0.08 * (vhi - vlo);
    vlo -= padv; vhi += padv;
    const double vr = (vhi - vlo > 1e-12) ? (vhi - vlo) : 1.0;

    const double pxPerSample = (double)pw / (double)(visN - 1);
    const double startPx = ml - (double)lo * pxPerSample;   // so column `lo` maps to x=ml

    auto xOf = [&](int col) { return startPx + (double)col * pxPerSample; };
    auto yOf = [&](double val) { return mt + ph - (val - vlo) / vr * ph; };

    // ---- 95% confidence band (filled, translucent) ----
    if (haveSd) {
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(70, 130, 180, 70));   // steel blue, translucent
        // One polygon per contiguous non-NaN run (same NaN handling as
        // BinPlotWidget::drawIqrBand).
        int i = lo;
        while (i <= hi) {
            while (i <= hi && (std::isnan(m_mean[i]) || std::isnan(ci[i]))) ++i;
            const int runStart = i;
            while (i <= hi && !std::isnan(m_mean[i]) && !std::isnan(ci[i])) ++i;
            const int runEnd = i;   // exclusive
            if (runEnd - runStart < 2) continue;
            QPainterPath band;
            for (int k = runStart; k < runEnd; ++k) {
                if (k == runStart) band.moveTo(xOf(k), yOf(m_mean[k] + ci[k]));
                else               band.lineTo(xOf(k), yOf(m_mean[k] + ci[k]));
            }
            for (int k = runEnd - 1; k >= runStart; --k)
                band.lineTo(xOf(k), yOf(m_mean[k] - ci[k]));
            band.closeSubpath();
            p.drawPath(band);
        }
    }

    // ---- mean trace (center line) ----
    {
        QPen pen(QColor(40, 40, 40)); pen.setWidthF(1.6);
        p.setPen(pen); p.setBrush(Qt::NoBrush);
        QPainterPath path; bool pend = true;
        for (int k = lo; k <= hi; ++k) {
            if (std::isnan(m_mean[k])) { pend = true; continue; }
            if (pend) { path.moveTo(xOf(k), yOf(m_mean[k])); pend = false; }
            else      path.lineTo(xOf(k), yOf(m_mean[k]));
        }
        p.drawPath(path);
    }

    // ---- fitted curve (from anchor_fit) ----
    {
        const std::vector<double> fit = fittedCurve(lo, hi);
        QPen pen(QColor(200, 60, 60)); pen.setWidthF(1.4); pen.setStyle(Qt::DashLine);
        p.setPen(pen); p.setBrush(Qt::NoBrush);
        QPainterPath path; bool pend = true;
        for (int k = lo; k <= hi; ++k) {
            if (std::isnan(fit[k])) { pend = true; continue; }
            if (pend) { path.moveTo(xOf(k), yOf(fit[k])); pend = false; }
            else      path.lineTo(xOf(k), yOf(fit[k]));
        }
        p.drawPath(path);
    }

    // ---- landmark marker (vertical line at landmarkCol) ----
    if (m_landmarkCol >= lo && m_landmarkCol <= hi) {
        QPen pen(QColor(20, 20, 20)); pen.setWidthF(1.0); pen.setStyle(Qt::DotLine);
        p.setPen(pen);
        const double x = xOf(m_landmarkCol);
        p.drawLine(QPointF(x, mt), QPointF(x, mt + ph));
    }

    // ---- footer: nBeats (the n in the CI) ----
    p.setPen(QColor(120, 120, 120));
    p.drawText(QRect(ml, mt + ph - 14, pw, 12),
        Qt::AlignRight | Qt::AlignVCenter,
        QStringLiteral("n=%1  (95% CI: mean +/- 1.96 sd/sqrt(n))").arg(m_nBeats));
}