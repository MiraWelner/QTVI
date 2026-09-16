#include "focus_panel_widget.hpp"
#include "subsample_refine.hpp" 
#include "template_anchoring\curve_fit.hpp"

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

void FocusPanelWidget::setSdMs(const std::vector<double>& sdMs,
    const std::vector<uint8_t>& floorMask,
    const std::vector<double>& deriv,
    double slopeFloor)
{
    m_sdMs = sdMs;
    m_floorMask = floorMask;
    m_deriv = deriv;
    m_slopeFloor = slopeFloor;
    update();
}

void FocusPanelWidget::clearFocus() {
    m_active = false;
    m_mean.clear();
    m_sd.clear();
    m_sdMs.clear();
    m_deriv.clear();
    m_floorMask.clear();
    m_nBeats = 0;
    m_landmarkCol = -1;
    m_fitKind = FitKind::Transition;
    m_transCands = subsample_refine::TransitionCandidates{};
    update();
}

// Fitted curve over [lo, hi] via anchor_fit's BIC model selection, sampled
// per column. NaN outside [lo, hi].
std::vector<double> FocusPanelWidget::peakCurve(int lo, int hi, bool cubic) const {
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> out(m_mean.size(), NaN);
    if (m_landmarkCol < 0 || m_landmarkCol >= (int)m_mean.size()) return out;
    if (lo < 0 || hi >= (int)m_mean.size() || hi - lo < 3) return out;

    // Fit the peak LOCALLY (tight window, where the parabola/cubic is a valid
    // model and the fitter's residual guard passes), then EXTRAPOLATE that
    // polynomial across the whole visible window via eval(). This always yields
    // a real curve -- the fit can't be rejected for not matching the wave's
    // non-parabolic shoulders, because it never sees them -- and draws it wide.
    const subsample_refine::ExtremumFit fit = cubic
        ? subsample_refine::asymmetricExtremumFit(m_mean, m_landmarkCol, 12.0)
        : subsample_refine::symmetricExtremumFit(m_mean, m_landmarkCol, 4.0);
    if (fit.order < 2) return out;   // degenerate solve: nothing to draw
    for (int i = lo; i <= hi; ++i) out[i] = fit.eval(static_cast<double>(i));
    return out;
}

std::vector<FocusPanelWidget::Candidate>
FocusPanelWidget::candidateCurves(int lo, int hi) const {
    std::vector<Candidate> out;
    if (lo < 0 || hi >= (int)m_mean.size() || hi - lo < 3) return out;

    // Peak fit/draw window: at least kWindowHalfWidth, but widened to the
    // landmark's own sigma so a BROAD peak (P, sigma 12) is seen over enough
    // samples to have real curvature. Over the fixed +-7 a broad wave fits with
    // a~=0 and the residual guard rejects it -> no polynomial -> nothing drawn.
    const int peakHw = std::max(subsample_refine::kWindowHalfWidth,
        static_cast<int>(std::lround(m_peakSigma)));

    auto evalExtremum = [&](const subsample_refine::ExtremumFit& f) {
        std::vector<double> c(m_mean.size(),
            std::numeric_limits<double>::quiet_NaN());
        if (f.order >= 2) {
            const int a = std::max(lo, f.seed - peakHw);
            const int b = std::min(hi, f.seed + peakHw);
            for (int i = a; i <= b; ++i) c[i] = f.eval(static_cast<double>(i));
        }
        return c;
        };
    auto evalModel = [&](const curve_fit::FitResult& f) {
        std::vector<double> c(m_mean.size(),
            std::numeric_limits<double>::quiet_NaN());
        if (f.f) for (int i = lo; i <= hi; ++i) c[i] = f.f(static_cast<double>(i));
        return c;
        };

    if (m_fitKind == FitKind::PeakQuadratic || m_fitKind == FitKind::PeakCubic) {
        if (m_landmarkCol < 0 || m_landmarkCol >= (int)m_mean.size()) return out;
        const auto q = subsample_refine::symmetricExtremumFit(m_mean, m_landmarkCol, m_peakSigma, peakHw);
        const auto c = subsample_refine::asymmetricExtremumFit(m_mean, m_landmarkCol, m_peakSigma, peakHw);
        const auto win = subsample_refine::bestPeakExtremumFit(m_mean, m_landmarkCol, m_peakSigma, peakHw, m_panelPeakMode);
        const bool winDegenerate = (win.order < 2);   // quad AND cubic both failed
        const bool fpIsWinner = winDegenerate
            || win.type == subsample_refine::CurveType::FIVE_POINT;   // forced 5-point
        if (q.order >= 2) out.push_back({ evalExtremum(q), !fpIsWinner && q.type == win.type, QStringLiteral("Quadratic") });
        if (c.order >= 2) out.push_back({ evalExtremum(c), !fpIsWinner && c.type == win.type, QStringLiteral("Cubic") });
        // ALWAYS draw the 5-point parabola so a broad/flat peak (e.g. P) still
        // has a visible curve. Green when it IS the placement (both models
        // degenerated) or the operator forced it, red otherwise.
        const auto fp5 = subsample_refine::fivePointParabolaFit(m_mean, m_landmarkCol);
        if (fp5.order >= 2) out.push_back({ evalExtremum(fp5), fpIsWinner, QStringLiteral("5-pt parabola") });
    }
    else {
        // Prefer the EXACT candidates the detector fit (supplied via
        // setTransitionCandidates): sample-indexed closures over the upsampled
        // one-sided window the panel can't reconstruct itself. Fall back to a
        // live re-fit over the visible window only when none were supplied.
        if (m_transCands.valid) {
            static const char* kNames[4] = { "Piecewise", "Sigmoid", "Fractional", "Cubic Spline" };
            for (int k = 0; k < 4; ++k) {
                const auto& fn = m_transCands.curve[k];
                std::vector<double> c(m_mean.size(),
                    std::numeric_limits<double>::quiet_NaN());
                if (fn) for (int i = lo; i <= hi; ++i) c[i] = fn(static_cast<double>(i));
                out.push_back({ c, k == m_transCands.winner, QString::fromUtf8(kNames[k]) });
            }
            return out;
        }
        // Same three candidates and same winner as selectAnchorModel, over the
        // visible window (fallback; see header note on the transition window).
        const auto pw = curve_fit::fitPiecewiseLinear(m_mean, lo, hi);
        const auto sig = curve_fit::fitSigmoid(m_mean, lo, hi, pw);
        const auto frac = curve_fit::fitFractionalPolynomial(m_mean, lo, hi);
        const auto win = curve_fit::selectAnchorModel(m_mean, lo, hi);
        out.push_back({ evalModel(pw),  pw.type == win.type, QStringLiteral("Piecewise") });
        out.push_back({ evalModel(sig), sig.type == win.type, QStringLiteral("Sigmoid") });
        out.push_back({ evalModel(frac), frac.type == win.type, QStringLiteral("Fractional") });
    }
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

    // Band half-width per column: THE SD ITSELF, +/- 1 sd.
    //
    // It used to be a 95% CI on the mean: 1.96 * sd / sqrt(nBeats). With
    // n=857 that divides the sd by ~29, so the band drew at ~7% of the actual
    // spread and looked like a hairline even where the sd was large. A CI on
    // the mean answers "how well do we know the average beat", which is not
    // the question this panel is for -- the spread of the beats is.
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> ci(N, NaN);
    const bool haveSd = ((int)m_sd.size() == N);
    if (haveSd) {
        for (int i = lo; i <= hi; ++i) {
            if (std::isnan(m_sd[i])) continue;
            ci[i] = m_sd[i];
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

    // ---- slope-floor shading -------------------------------------------
    // Columns where the local |dV/dt| was clamped at the floor: flat regions
    // and peak tops, where SD/slope would otherwise divide by ~0. The msec
    // value there is a lower bound, so the columns are washed out rather than
    // left looking like the rest.
    if ((int)m_floorMask.size() == N) {
        p.setPen(Qt::NoPen);
        p.setBrush(QColor(200, 90, 60, 28));
        int i = lo;
        while (i <= hi) {
            while (i <= hi && !m_floorMask[i]) ++i;
            const int runStart = i;
            while (i <= hi && m_floorMask[i]) ++i;
            if (i > runStart) {
                const double x0 = xOf(runStart);
                const double x1 = xOf(i - 1) + pxPerSample;
                p.drawRect(QRectF(x0, mt, std::max(1.0, x1 - x0), ph));
            }
        }
    }

    // ---- +/- 1 sd band (filled, translucent) ----
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

    // ---- fitted curve(s): every tested model; winner green, others red ----
    {
        auto drawCurve = [&](const std::vector<double>& fit, QColor col) {
            QPen pen(col); pen.setWidthF(1.4); pen.setStyle(Qt::DashLine);
            p.setPen(pen); p.setBrush(Qt::NoBrush);
            QPainterPath path; bool pend = true;
            for (int k = lo; k <= hi; ++k) {
                if (std::isnan(fit[k])) { pend = true; continue; }
                if (pend) { path.moveTo(xOf(k), yOf(fit[k])); pend = false; }
                else      path.lineTo(xOf(k), yOf(fit[k]));
            }
            p.drawPath(path);
            };
        // Draw losers first, winner last so the green sits on top.
        const std::vector<Candidate> cands = candidateCurves(lo, hi);
        for (const Candidate& c : cands)
            if (!c.selected) drawCurve(c.curve, QColor(200, 60, 60));   // red
        for (const Candidate& c : cands)
            if (c.selected)  drawCurve(c.curve, QColor(0, 150, 0));     // green

        // Name the winning model, in the winner's green, top-right.
        for (const Candidate& c : cands)
            if (c.selected && !c.label.isEmpty()) {
                p.setPen(QColor(0, 130, 0));
                p.drawText(QRect(ml, 4, width() - ml - mr, 18),
                    Qt::AlignRight | Qt::AlignVCenter, c.label);
                break;
            }
    }

    // ---- fiducial marker (dotted vertical line) ----
    // For a transition, this is WHERE THE SELECTED MODEL places the mark -- its
    // own crossing (cross[winner]), so the dotted line tracks the green curve as
    // you switch models. Falls back to the stored landmark column otherwise.
    double fidCol = static_cast<double>(m_landmarkCol);
    if (m_transCands.valid && m_transCands.winner >= 0 && m_transCands.winner < 4
        && m_transCands.cross[m_transCands.winner] >= 0.0)
        fidCol = m_transCands.cross[m_transCands.winner];
    if (fidCol >= lo && fidCol <= hi) {
        QPen pen(QColor(20, 20, 20)); pen.setWidthF(1.0); pen.setStyle(Qt::DotLine);
        p.setPen(pen);
        const double x = xOf(fidCol);
        p.drawLine(QPointF(x, mt), QPointF(x, mt + ph));
    }

    // ---- footer: the two inputs to the msec SD, so the equation is visible ----
    // sd_ms = raw_sd / slope x (1000/fs). The header prints the resulting sd_ms;
    // here we show raw_sd and the DENOMINATOR ACTUALLY USED at the bar. Where the
    // column's own slope is below the floor (see the orange/pink shading), the
    // divide used the floor, not the near-zero slope -- so print the floor and
    // mark it, or raw_sd/slope would look like a divide-by-nothing. Slope is
    // printed at 3 decimals so a small-but-nonzero value is not shown as 0.0.
    p.setPen(QColor(120, 120, 120));
    QString foot = QStringLiteral("(band: mean +/- 1 sd)");
    if (m_landmarkCol >= 0) {
        const double rawSd = (m_landmarkCol < (int)m_sd.size())
            ? m_sd[m_landmarkCol] : std::numeric_limits<double>::quiet_NaN();
        const double slope = (m_landmarkCol < (int)m_deriv.size())
            ? m_deriv[m_landmarkCol] : std::numeric_limits<double>::quiet_NaN();
        const bool floored = (m_landmarkCol < (int)m_floorMask.size())
            && m_floorMask[m_landmarkCol];
        const double denom = floored ? m_slopeFloor : slope;
        const QString sdStr = std::isfinite(rawSd)
            ? QString::number(rawSd, 'f', 4) : QStringLiteral("--");
        const QString dStr = std::isfinite(denom)
            ? QString::number(denom, 'f', 3) : QStringLiteral("--");
        foot = floored
            ? QStringLiteral("raw sd=%1  slope=%2 /sample (floor)").arg(sdStr, dStr)
            : QStringLiteral("raw sd=%1  slope=%2 /sample").arg(sdStr, dStr);
    }
    p.drawText(QRect(ml, mt + ph - 14, pw, 12),
        Qt::AlignRight | Qt::AlignVCenter, foot);
}