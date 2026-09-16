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

    // drawHw defaults to the peak window. Each model is drawn only over the
    // span it was FIT on: the 5-point parabola over +-4 rather than +-sigma,
    // where a 5-sample parabola extrapolates straight off the top of the panel.
    auto evalExtremum = [&](const subsample_refine::peak_fit& f, int drawHw = -1) {
        std::vector<double> c(m_mean.size(),
            std::numeric_limits<double>::quiet_NaN());
        if (f.order >= 2) {
            const int hw = (drawHw > 0) ? drawHw : peakHw;
            const int a = std::max(lo, f.seed - hw);
            const int b = std::min(hi, f.seed + hw);
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
        // DRAW-ONLY fits (applyGuard=false). The guarded versions collapse a
        // residual-rejected quadratic/cubic to FIVE_POINT with order 0 and no
        // coefficients, and the order>=2 test below then dropped them -- which
        // is why, on a broad peak where BOTH were rejected, the 5-point
        // parabola was the only curve on screen. The guard decides what may
        // PLACE the mark; it should not decide what is VISIBLE, since the
        // rejected curve is exactly what shows why the fallback was taken.
        const auto qD = subsample_refine::quadratic_fit(
            m_mean, m_landmarkCol, m_peakSigma, peakHw, /*applyGuard=*/false);
        const auto cD = subsample_refine::cubic_fit(
            m_mean, m_landmarkCol, m_peakSigma, peakHw, /*applyGuard=*/false);
        const auto fp5 = subsample_refine::fivePointParabolaFit(m_mean, m_landmarkCol);

        // The GUARDED contest, honouring the Fit-Peaks radio: this is the model
        // that actually places the mark, and win.position is where. The winner
        // is read off the returned TYPE rather than re-derived from the radio,
        // so a forced model that degenerated shows its fallback as green
        // instead of colouring a curve that placed nothing.
        const auto win = subsample_refine::bestPeakExtremumFit(
            m_mean, m_landmarkCol, m_peakSigma, peakHw, m_panelPeakMode);
        int winIdx = 2;   // 0=quadratic, 1=cubic, 2=five-point
        switch (win.type) {
        case subsample_refine::CurveType::QUADRATIC: winIdx = 0; break;
        case subsample_refine::CurveType::CUBIC:     winIdx = 1; break;
        default:                                     winIdx = 2; break;
        }

        auto push = [&](const subsample_refine::peak_fit& f, int idx,
            const QString& name, int drawHw = -1) {
                if (f.order < 2) return;   // genuinely no polynomial: nothing to draw
                Candidate cd;
                cd.curve = evalExtremum(f, drawHw);
                cd.selected = (idx == winIdx);
                cd.label = name;
                cd.position = (idx == winIdx) ? win.position : f.position;
                out.push_back(std::move(cd));
            };
        push(qD, 0, QStringLiteral("Quadratic"));
        push(cD, 1, QStringLiteral("Cubic"));
        // ALWAYS draw the 5-point parabola so a broad/flat peak (e.g. P) still
        // has a visible curve. Green when it IS the placement (both models
        // degenerated) or the operator forced it, red otherwise.
        push(fp5, 2, QStringLiteral("5-pt parabola"), 4);
    }
    else {
        // Prefer the EXACT candidates the detector fit (supplied via
        // setTransitionCandidates): sample-indexed closures over the upsampled
        // one-sided window the panel can't reconstruct itself. Fall back to a
        // live re-fit over the visible window only when none were supplied.
        if (m_transCands.valid) {
            static const char* kNames[5] = { "Piecewise", "Sigmoid", "Fractional", "Cubic Spline", "Cubic" };
            for (int k = 0; k < 5; ++k) {
                const auto& fn = m_transCands.curve[k];
                std::vector<double> c(m_mean.size(),
                    std::numeric_limits<double>::quiet_NaN());
                if (fn) for (int i = lo; i <= hi; ++i) c[i] = fn(static_cast<double>(i));
                Candidate cd;
                cd.curve = std::move(c);
                cd.selected = (k == m_transCands.winner);
                cd.label = QString::fromUtf8(kNames[k]);
                // cross[k] is THIS model's own fiducial crossing, already in
                // sample coordinates (transitionAnchor fills it, and its own
                // comment says the focus dotted line should use cross[winner]).
                // It was computed and then never read: the panel drew its line
                // at the integer bar column, so switching the on/offset model
                // recoloured the curves and moved nothing.
                cd.position = (m_transCands.cross[k] >= 0.0)
                    ? m_transCands.cross[k]
                    : std::numeric_limits<double>::quiet_NaN();
                out.push_back(std::move(cd));
            }
            return out;
        }
        // Same three candidates and same winner as selectBestFit, over the
        // visible window (fallback; see header note on the transition window).
        const auto pw = curve_fit::fitPiecewiseLinear(m_mean, lo, hi);
        const auto sig = curve_fit::fitSigmoid(m_mean, lo, hi, pw);
        const auto frac = curve_fit::fitFractionalPolynomial(m_mean, lo, hi);
        const auto win = curve_fit::selectBestFit(m_mean, lo, hi);
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

    const int focuspanel_top_margin = 40;
    const int focuspanel_bottom_margin = 12;
    const int ml = 8, mr = 8;
    const int ph = height() - focuspanel_top_margin - focuspanel_bottom_margin;
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

    // SUB-SAMPLE. Was int, so every call site silently truncated -- including
    // the fiducial, which sits at a fractional column.
    auto xOf = [&](double col) { return startPx + col * pxPerSample; };
    auto yOf = [&](double val) { return focuspanel_top_margin + ph - (val - vlo) / vr * ph; };

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
                p.drawRect(QRectF(x0, focuspanel_top_margin, std::max(1.0, x1 - x0), ph));
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
        auto drawCurve = [&](const std::vector<double>& fit, QColor col, Qt::PenStyle style) {
            QPen pen(col); pen.setWidthF(1.4); pen.setStyle(style);
            p.setPen(pen); p.setBrush(Qt::NoBrush);
            QPainterPath path; bool pend = true;
            for (int k = lo; k <= hi; ++k) {
                if (std::isnan(fit[k])) { pend = true; continue; }
                if (pend) { path.moveTo(xOf(k), yOf(fit[k])); pend = false; }
                else      path.lineTo(xOf(k), yOf(fit[k]));
            }
            p.drawPath(path);
            };
        // Losers first, winner (green) on top.
        const std::vector<Candidate> cands = candidateCurves(lo, hi);
        for (const Candidate& c : cands)
            if (!c.selected) drawCurve(c.curve, QColor(200, 60, 60), Qt::DashLine);
        for (const Candidate& c : cands)
            if (c.selected)  drawCurve(c.curve, QColor(0, 150, 0), Qt::DashLine);

        // Name the winning model, in the winner's green, top-right.
        for (const Candidate& c : cands)
            if (c.selected && !c.label.isEmpty()) {
                p.setPen(QColor(120, 120, 120));
                p.drawText(QRect(ml, 22, width() - ml - mr, 16),
                    Qt::AlignLeft | Qt::AlignVCenter, c.label);
                break;
            }

        // ---- fiducial marker: gray vertical dotted line, at the column the
        // SELECTED model places the landmark on --------------------------
        //
        // Same single gray line as before; the only change is WHERE. It used to
        // be drawn at m_landmarkCol -- the integer bar column, which does not
        // depend on the fit at all -- so changing either fit-model radio
        // recoloured the curves and left the line sitting still.
        double winPos = std::numeric_limits<double>::quiet_NaN();
        for (const Candidate& c : cands)
            if (c.selected) { winPos = c.position; break; }
        // The bar column stays the fallback: a selected model with no placement
        // to report (a rejected fit with no vertex, or the transition re-fit
        // path, which has no crossing to hand back) leaves the line exactly
        // where it has always been rather than dropping it off the panel.
        const double fidCol = std::isfinite(winPos)
            ? winPos : static_cast<double>(m_landmarkCol);
        m_lastFidCol = fidCol;
        if (fidCol >= lo && fidCol <= hi) {
            QPen pen(QColor(130, 130, 130)); pen.setWidthF(1.2); pen.setStyle(Qt::DotLine);
            p.setPen(pen);
            const double x = xOf(fidCol);
            p.drawLine(QPointF(x, focuspanel_top_margin), QPointF(x, focuspanel_top_margin + ph));
        }
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
        // The fitted fiducial's own column, so a sub-pixel move is still
        // readable as a change.
        if (m_lastFidCol >= 0.0)
            foot = QStringLiteral("fid=%1  ").arg(m_lastFidCol, 0, 'f', 2) + foot;
    }
    p.drawText(QRect(ml, focuspanel_top_margin + ph - 14, pw, 12),
        Qt::AlignRight | Qt::AlignVCenter, foot);
}