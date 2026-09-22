#include "focus_panel_widget.hpp"
#include "subsample_refine.hpp"

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

    // A NEW LANDMARK INVALIDATES EVERY DERIVED ARRAY, the supplied fits
    // included. They describe the PREVIOUS landmark on the PREVIOUS waveform,
    // and every caller that has replacements supplies them immediately after
    // this call -- so clearing here costs those callers nothing.
    //
    // THE FITS WERE NOT CLEARED HERE, and only clearFocus did it. A caller
    // that supplied none therefore inherited the last caller's: select an ECG
    // Q-onset, then a PPG landmark, and the panel drew the ECG's five
    // transition closures over the pulse, named an ECG model in the header, and
    // put the dotted fiducial at the ECG's detector column. Clearing is the
    // fix, not a cleanup -- a panel holding another landmark's fit has no
    // correct way to draw itself.
    m_sdMs.clear();
    m_floorMask.clear();
    m_deriv.clear();
    m_slopeFloor = 0.0;
    m_lastFidCol = -1.0;
    m_detectorFid = -1.0;
    m_transCands = subsample_refine::TransitionCandidates{};
    m_peakCands = subsample_refine::PeakCandidates{};

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

void FocusPanelWidget::setDetectorFiducial(double col) {
    m_detectorFid = col;
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
    m_peakCands = subsample_refine::PeakCandidates{};
    m_detectorFid = -1.0;   // stale for the new landmark until re-supplied
    update();
}

std::vector<FocusPanelWidget::Candidate>
FocusPanelWidget::candidateCurves(int lo, int hi) const {
    std::vector<Candidate> out;
    if (m_fitKind == FitKind::None) return out;   // nothing fitted it
    if (lo < 0 || hi >= (int)m_mean.size() || hi - lo < 3) return out;

    // THE DETECTOR'S WINDOW, EXACTLY -- which is now per landmark, handed in
    // with the sigma by setFitKind rather than read off a global constant. Any
    // other value makes this panel fit a different span from the one that
    // placed the mark, and then the dotted fiducial, the curves and the fid=
    // readout all describe a position nothing else in the system uses.
    //
    // THIS IS NOT THE max(<global constant>, round(sigma)) THAT WAS REVERTED.
    // That widened the panel's fit ALONE, leaving the detector at +-7, so the
    // drawn curve and the placed mark came from different fits. The widths are
    // now in subsample_refine::peak_halfwidth and the detector reads the same
    // table, so the two cannot diverge -- the sameness is what mattered, not
    // the number 7.
    // -1 when no peak is in focus; nothing below reads it on that path, since
    // evalExtremum is only called from the peak branch. Floored at 3 anyway so
    // a bad value cannot produce a degenerate window.
    const int peakHw = std::max(3, m_peakHalfWidth);

    // How far the curves are DRAWN. Still allowed to exceed the fit so a broad
    // wave shows its shape, but the floor is now the landmark's own fit window
    // rather than the global constant -- otherwise a P fitted over +-24 would
    // be drawn over +-12 and look NARROWER than it was fitted.
    const int drawSpan = std::max(peakHw,
        static_cast<int>(std::lround(m_peakSigma)));

    // drawHw defaults to the peak window. Each model is drawn only over the
    // span it was FIT on: the 5-point parabola over +-4 rather than +-sigma,
    // where a 5-sample parabola extrapolates straight off the top of the panel.
    auto evalExtremum = [&](const subsample_refine::peak_fit& f, int drawHw = -1) {
        std::vector<double> c(m_mean.size(),
            std::numeric_limits<double>::quiet_NaN());
        if (f.order >= 2) {
            const int hw = (drawHw > 0) ? drawHw : drawSpan;
            const int a = std::max(lo, f.seed - hw);
            const int b = std::min(hi, f.seed + hw);
            for (int i = a; i <= b; ++i) c[i] = f.eval(static_cast<double>(i));
        }
        return c;
        };
    if (m_fitKind == FitKind::PeakQuadratic || m_fitKind == FitKind::PeakCubic) {
        // NOTHING IS FITTED HERE. The three curves, the winner and the
        // placement all arrive in m_peakCands, run once by
        // subsample_refine::peakCandidates. The panel used to re-run the whole
        // contest on a seed it chose itself, so the drawn curve and the marked
        // line came from two different fits of the same data.
        if (!m_peakCands.valid || m_peakCands.winner < 0) return out;

        auto push = [&](int idx, const QString& name, int drawHw = -1) {
            const subsample_refine::peak_fit& f = m_peakCands.draw[idx];
            if (f.order < 2) return;   // genuinely no polynomial: nothing to draw
            Candidate cd;
            cd.curve = evalExtremum(f, drawHw);
            cd.selected = (idx == m_peakCands.winner);
            cd.label = f.guardFailed
                ? name + QStringLiteral(" (rejected)") : name;
            // The winner reports the placement every consumer uses; the
            // losers report their own vertex.
            cd.position = (idx == m_peakCands.winner)
                ? m_peakCands.placement : f.position;
            out.push_back(std::move(cd));
            };
        push(0, QStringLiteral("Quadratic"));
        push(1, QStringLiteral("Cubic"));
        // ALWAYS draw the 5-point parabola so a broad/flat peak (e.g. P) still
        // has a visible curve. Green when it IS the placement (both models
        // degenerated) or the operator forced it, red otherwise.
        push(2, QStringLiteral("5-pt parabola"), 4);
    }

    else {
        // THE EXACT CANDIDATES THE DETECTOR FIT, or nothing. Supplied via
        // setTransitionCandidates: sample-indexed closures over the upsampled
        // one-sided window the panel cannot reconstruct itself.
        //
        // A live re-fit used to stand behind this -- fitPiecewiseLinear /
        // fitSigmoid / fitFractionalPolynomial over the visible window -- and
        // it was worse than drawing nothing. It answered a different question
        // from the detector's (a different span, a different weighting) and it
        // was reached by every landmark with no contest behind it, including
        // every PULSE landmark, where it fitted onset/offset models across a
        // whole pulse and drew three curves nothing had placed. The panel
        // renders supplied answers; it does not produce them.
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
        }
    }
    return out;
}

void FocusPanelWidget::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.fillRect(rect(), QColor(250, 250, 250));

    // TWO HEADER LINES, TWO FOOTER LINES. The header carries the landmark and
    // its alignment on line 1 and the selected fit model on line 2; the sd
    // breakdown sits below the plot rather than inside it, so nothing overlaps
    // the trace and nothing shares a line with the landmark name.
    const int kHeadLine = 18;          // landmark + alignment
    const int kSubLine = 15;           // selected model name
    const int kFootLine = 14;          // one sd line
    const int mt = 4 + kHeadLine + kSubLine;          // top margin (both header lines)
    const int mb = 6 + 2 * kFootLine;                 // bottom margin (both sd lines)
    const int ml = 8, mr = 8;
    const int ph = height() - mt - mb;
    const int pw = width() - ml - mr;

    // Header line 1: the landmark and its alignment.
    p.setPen(QColor(60, 60, 60));
    p.drawText(QRect(ml, 4, width() - ml - mr, kHeadLine),
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

        // HEADER LINE 2: the selected model's name, in GRAY, on its own line
        // under the landmark. It was green and right-aligned on the landmark's
        // own line, which read as part of the landmark's name and competed with
        // it. The winning CURVE keeps its green -- that is what the colour is
        // for, and it is unambiguous next to the red losers.
        if (m_fitKind == FitKind::None) {
            // The line that would name the winning model says why there is
            // none, rather than sitting blank and reading as "still loading".
            // GENERIC, because two things reach it: an ECG Q with no trough
            // (compute_q_onset's R-upstroke fallback) and every PULSE landmark,
            // none of which is placed by a model contest. Both are "the
            // detector's position, with no curve behind it".
            p.setPen(QColor(120, 120, 120));
            p.drawText(QRect(ml, 4 + kHeadLine, width() - ml - mr, kSubLine),
                Qt::AlignLeft | Qt::AlignVCenter,
                QStringLiteral("no fit - detector position"));
        }
        else for (const Candidate& c : cands)
            if (c.selected && !c.label.isEmpty()) {
                p.setPen(QColor(120, 120, 120));
                p.drawText(QRect(ml, 4 + kHeadLine, width() - ml - mr, kSubLine),
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
        // THE DETECTOR'S POSITION, not a re-fit's. See m_detectorFid: the
        // panel cannot reproduce the detector's vertex because it does not know
        // the integer seed the detector fitted around. The transition path is
        // exempt -- its cross[winner] comes from transitionAnchor itself, so it
        // already IS the detector's answer -- and the selected candidate's own
        // position is the last resort when nothing was supplied.
        double winPos = std::numeric_limits<double>::quiet_NaN();
        if (m_detectorFid >= 0.0) {
            winPos = m_detectorFid;
        }
        else {
            for (const Candidate& c : cands)
                if (c.selected) { winPos = c.position; break; }
        }
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
            p.drawLine(QPointF(x, mt), QPointF(x, mt + ph));
        }
    }

    // ---- THE SD, TWO LINES, BELOW THE PLOT -----------------------------
    //
    // Line 1 is the number: sd in milliseconds at the bar's own column.
    // Line 2 is the equation behind it -- sd_ms = raw_sd / slope x (1000/fs) --
    // so the result and its two inputs are both visible without reading the
    // header. Where the column's own slope is below the floor (the orange/pink
    // shading) the divide used the floor, not the near-zero slope, so the floor
    // is printed and flagged; otherwise raw_sd/slope would look like a
    // divide-by-nothing. Slope is at 3 decimals so a small-but-nonzero value is
    // not shown as 0.0.
    //
    // BELOW the plot, not inside it: this used to draw at mt + ph - 14, i.e.
    // over the bottom of the trace and the band.
    {
        const int fy = mt + ph + 4;
        p.setPen(QColor(120, 120, 120));

        QString l1, l2;
        if (m_landmarkCol >= 0) {
            const double NaNv = std::numeric_limits<double>::quiet_NaN();
            const double sdMs = (m_landmarkCol < (int)m_sdMs.size())
                ? m_sdMs[m_landmarkCol] : NaNv;
            const double rawSd = (m_landmarkCol < (int)m_sd.size())
                ? m_sd[m_landmarkCol] : NaNv;
            const double slope = (m_landmarkCol < (int)m_deriv.size())
                ? m_deriv[m_landmarkCol] : NaNv;
            const bool floored = (m_landmarkCol < (int)m_floorMask.size())
                && m_floorMask[m_landmarkCol];
            const double denom = floored ? m_slopeFloor : slope;

            l1 = std::isfinite(sdMs)
                ? QStringLiteral("sd = %1 ms").arg(sdMs, 0, 'f', 1)
                : QStringLiteral("sd = --");

            const QString sdStr = std::isfinite(rawSd)
                ? QString::number(rawSd, 'f', 4) : QStringLiteral("--");
            const QString dStr = std::isfinite(denom)
                ? QString::number(denom, 'f', 3) : QStringLiteral("--");
            l2 = floored
                ? QStringLiteral("raw sd = %1   slope = %2 /sample (floor)").arg(sdStr, dStr)
                : QStringLiteral("raw sd = %1   slope = %2 /sample").arg(sdStr, dStr);
        }

        if (!l1.isEmpty())
            p.drawText(QRect(ml, fy, pw, kFootLine),
                Qt::AlignLeft | Qt::AlignVCenter, l1);
        if (!l2.isEmpty())
            p.drawText(QRect(ml, fy + kFootLine, pw, kFootLine),
                Qt::AlignLeft | Qt::AlignVCenter, l2);
    }
}