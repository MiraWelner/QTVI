/**
 * @file   grid_overlay.cpp
 * @brief  See grid_overlay.hpp.
 */
#include "grid_overlay.hpp"
#include "gui_handler.h"

#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <cmath>
#include <limits>

namespace {
    // Vertical spacings (in seconds along the x-axis). The horizontal
    // spacings are derived from these on every refresh so cells render
    // as visual squares regardless of the y-axis range or units.
    constexpr double kMinorSpacingSec = 0.04;
    constexpr double kMajorSpacingSec = 0.2;
}

pulse_overlay::pulse_overlay(noise_marking_gui* gui, const QStringList& labels)
    : QObject(gui)
    , m_gui(gui)
    , m_labels(labels)
    // QColor(r, g, b, alpha): alpha 0..255 (0 = invisible, 255 = solid).
    // Pen width is in pixels and can be fractional.
    // Major lines (every 0.2 s on x, ~equivalent in y) are more opaque
    // AND thicker than minor lines (every 0.04 s) so the eye can easily
    // pick them out as a coarse ruler against the finer subdivisions.
    , m_minorPen(QColor(255, 0, 0, 80), 0.6)
    , m_majorPen(QColor(255, 0, 0, 130), 1.0)
{
    for (const QString& label : m_labels) {
        ChannelSeries cs;
        auto makeSeries = [this](QLineSeries*& s, const QPen& pen) {
            s = new QLineSeries();
            s->setPen(pen);
            // OpenGL OFF: the GL backend ignores pen alpha and clamps
            // width to >= 1, so translucency / sub-pixel width would be
            // invisible. The overlay has only a few hundred points per
            // series, so software rendering is plenty fast.
            s->setUseOpenGL(false);
            };
        makeSeries(cs.vMinor, m_minorPen);
        makeSeries(cs.vMajor, m_majorPen);
        makeSeries(cs.hMinor, m_minorPen);
        makeSeries(cs.hMajor, m_majorPen);

        m_series.insert(label, cs);
    }
}

pulse_overlay::~pulse_overlay() {
    for (const ChannelSeries& cs : m_series) {
        // The host chart may still hold references; detach first so
        // ~QChart doesn't try to delete them again.
        for (QLineSeries* s : { cs.vMinor, cs.vMajor, cs.hMinor, cs.hMajor }) {
            if (s && s->chart()) s->chart()->removeSeries(s);
            delete s;
        }
    }
    m_series.clear();
}

QList<QLineSeries*> pulse_overlay::seriesForLabel(const QString& label) const {
    QList<QLineSeries*> out;
    auto it = m_series.constFind(label);
    if (it == m_series.constEnd()) return out;
    if (it->vMinor) out.append(it->vMinor);
    if (it->vMajor) out.append(it->vMajor);
    if (it->hMinor) out.append(it->hMinor);
    if (it->hMajor) out.append(it->hMajor);
    return out;
}

void pulse_overlay::setEnabled(bool on) {
    if (on == m_enabled) return;
    m_enabled = on;

    if (!on) {
        // Detach every series from its host chart so the red lines vanish.
        // Don't delete; we want to be able to re-attach in refresh() if
        // the user toggles back on.
        for (const ChannelSeries& cs : m_series) {
            for (QLineSeries* s : { cs.vMinor, cs.vMajor, cs.hMinor, cs.hMajor }) {
                if (s && s->chart()) s->chart()->removeSeries(s);
            }
        }
    }
    // On re-enable, refresh() (called by the GUI right after this) re-
    // attaches the series and rebuilds the geometry.
}

void pulse_overlay::hookChartIfNeeded(const QString& label, QChartView* cv,
    QValueAxis* yAxis)
{
    // Guard against double-subscription: refresh() can run many times
    // per second and walks every label, so we'd otherwise stack up
    // duplicate connections every redraw. QPointer-based membership
    // also drops dangling entries silently when a chart is destroyed.
    auto containsLive = [](const QList<QPointer<QObject>>& v, QObject* p) {
        for (const auto& q : v) if (q.data() == p) return true;
        return false;
        };

    QChart* chart = cv->chart();
    if (chart && !containsLive(m_hookedCharts, chart)) {
        // plotAreaChanged fires on resize, on margin change, and on
        // axis-range change after layout settles. That's exactly when
        // pixels-per-x-second or pixels-per-y-unit could shift, so
        // it's the right hook for re-squaring.
        connect(chart, &QChart::plotAreaChanged,
            this, [this](const QRectF&) { refresh(); });
        m_hookedCharts.append(QPointer<QObject>(chart));
    }

    if (yAxis && !containsLive(m_hookedAxes, yAxis)) {
        // The y-axis range often changes WITHOUT a plotAreaChanged
        // fire (e.g. the channel autoscales but the chart's pixel
        // rectangle didn't move). Hook rangeChanged so y-spacing
        // tracks autoscale.
        connect(yAxis, &QValueAxis::rangeChanged,
            this, [this](qreal, qreal) { refresh(); });
        m_hookedAxes.append(QPointer<QObject>(yAxis));
    }

    // Sweep dead pointers so the membership lists don't grow forever
    // as the host wipes and recreates charts on data reload.
    auto sweep = [](QList<QPointer<QObject>>& v) {
        for (int i = v.size() - 1; i >= 0; --i)
            if (!v[i]) v.removeAt(i);
        };
    sweep(m_hookedCharts);
    sweep(m_hookedAxes);

    (void)label;  // reserved for future per-channel diagnostics
}

void pulse_overlay::refresh() {
    if (!m_enabled) return;
    if (m_inRefresh) return;     // see m_inRefresh comment in header
    m_inRefresh = true;
    // RAII clear so any early return / exception still releases the guard.
    struct Guard {
        bool* p;
        ~Guard() { *p = false; }
    } guard{ &m_inRefresh };

    // --- Geometry builders --------------------------------------------------
    // Both directions use the same trick: pack N parallel line segments
    // into one polyline with NaN breakpoints so Qt Charts renders them
    // as disconnected strokes within a single series. Cheaper than one
    // QLineSeries per gridline.
    //
    // Minor and major are computed from independent loops with their own
    // spacings rather than "every 5th minor is a major" because that
    // approach is fragile when the starting coordinate isn't aligned to
    // a major. Doing them independently keeps each grid snapped to its
    // own spacing.
    const double nan = std::numeric_limits<double>::quiet_NaN();

    auto buildVertical = [&](double tStart, double tEnd, double spacing,
        double yMin, double yMax) -> QList<QPointF> {
            QList<QPointF> pts;
            if (tEnd <= tStart || spacing <= 0.0) return pts;
            const double firstX = std::ceil(tStart / spacing) * spacing;
            const int approxN = static_cast<int>((tEnd - tStart) / spacing) + 2;
            pts.reserve(approxN * 3);
            for (double x = firstX; x <= tEnd; x += spacing) {
                pts.append({ x, yMin });
                pts.append({ x, yMax });
                pts.append({ nan, nan });
            }
            return pts;
        };

    auto buildHorizontal = [&](double yMin, double yMax, double spacing,
        double tStart, double tEnd) -> QList<QPointF> {
            QList<QPointF> pts;
            if (yMax <= yMin || spacing <= 0.0) return pts;
            const double firstY = std::ceil(yMin / spacing) * spacing;
            const int approxN = static_cast<int>((yMax - yMin) / spacing) + 2;
            pts.reserve(approxN * 3);
            for (double y = firstY; y <= yMax; y += spacing) {
                pts.append({ tStart, y });
                pts.append({ tEnd,   y });
                pts.append({ nan,    nan });
            }
            return pts;
        };

    for (const QString& label : m_labels) {
        auto it = m_series.find(label);
        if (it == m_series.end()) continue;
        QLineSeries* vMinor = it->vMinor;
        QLineSeries* vMajor = it->vMajor;
        QLineSeries* hMinor = it->hMinor;
        QLineSeries* hMajor = it->hMajor;
        if (!vMinor || !vMajor || !hMinor || !hMajor) continue;

        // The host chart wipes and re-creates its axes on every redraw,
        // so any axis attachments we had are stale.
        for (QLineSeries* s : { vMinor, vMajor, hMinor, hMajor })
            if (s->chart()) s->chart()->removeSeries(s);

        if (!m_gui->isChannelActive(label)) continue;
        QChartView* cv = m_gui->chartViewForSignalLabel(label);
        if (!cv || !cv->chart()) continue;

        auto hAxes = cv->chart()->axes(Qt::Horizontal);
        auto vAxes = cv->chart()->axes(Qt::Vertical);
        if (hAxes.isEmpty() || vAxes.isEmpty()) continue;
        auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first());
        if (!yAxis) continue;

        // Make sure we'll re-refresh when this chart resizes or its
        // y-axis autoscales. Idempotent.
        hookChartIfNeeded(label, cv, yAxis);

        const double yMin = yAxis->min();
        const double yMax = yAxis->max();
        const double tStart = m_gui->currentStartTime();
        const double tEnd = tStart + m_gui->windowDuration();

        // --- Square-cell y-spacings -------------------------------------
        // To make a grid cell render as a visual square, the spacing in
        // y-units must equal the spacing in x-seconds scaled by the
        // axes' pixels-per-unit ratio:
        //
        //   pixels_per_second = plotArea.width()  / (tEnd - tStart)
        //   pixels_per_yunit  = plotArea.height() / (yMax - yMin)
        //
        //   yUnits_per_xSecond = pixels_per_second / pixels_per_yunit
        //                      = (plotArea.width  * (yMax - yMin))
        //                      / (plotArea.height * (tEnd - tStart))
        //
        //   ySpacing = xSpacing * yUnits_per_xSecond
        //
        // If the plot area or axis range is degenerate, we skip the
        // horizontal grid this frame -- next refresh() will catch it
        // once Qt finishes laying the chart out.
        const QRectF plot = cv->chart()->plotArea();
        const double dx = tEnd - tStart;
        const double dy = yMax - yMin;
        double minorSpacingY = 0.0;
        double majorSpacingY = 0.0;
        if (plot.width() > 0.0 && plot.height() > 0.0 && dx > 0.0 && dy > 0.0) {
            const double yUnitsPerXSec =
                (plot.width() * dy) / (plot.height() * dx);
            minorSpacingY = kMinorSpacingSec * yUnitsPerXSec;
            majorSpacingY = kMajorSpacingSec * yUnitsPerXSec;
        }

        // --- Bake geometry --------------------------------------------------
        vMinor->replace(buildVertical(tStart, tEnd, kMinorSpacingSec, yMin, yMax));
        vMajor->replace(buildVertical(tStart, tEnd, kMajorSpacingSec, yMin, yMax));
        hMinor->replace(buildHorizontal(yMin, yMax, minorSpacingY, tStart, tEnd));
        hMajor->replace(buildHorizontal(yMin, yMax, majorSpacingY, tStart, tEnd));

        // Attach order matters for z-order: later additions paint on
        // top. We want major (more opaque) over minor, and a consistent
        // vertical-vs-horizontal stacking. Horizontal-then-vertical so
        // the vertical lines (the primary time ruler) win at crossings.
        auto attach = [&](QLineSeries* s) {
            cv->chart()->addSeries(s);
            s->attachAxis(hAxes.first());
            s->attachAxis(yAxis);
            };
        attach(hMinor);
        attach(hMajor);
        attach(vMinor);
        attach(vMajor);
    }

}