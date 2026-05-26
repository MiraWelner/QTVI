/**
 * @file   chart_overlay.cpp
 * @brief  See chart_overlay.hpp.
 */
#include "chart_overlay.hpp"
#include "gui_handler.h"

#include <QtCharts/QChartView>
#include <QtCharts/QChart>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <cmath>
#include <limits>

namespace {
    constexpr double kMinorSpacingSec = 0.04;
    constexpr double kMajorSpacingSec = 0.2;
}

pulse_overlay::pulse_overlay(noise_marking_gui* gui, const QStringList& labels)
    : QObject(gui)
    , m_gui(gui)
    , m_labels(labels)
    // QColor(r, g, b, alpha): alpha 0..255 (0 = invisible, 255 = solid).
    // Pen width is in pixels and can be fractional.
    // Major lines (every 0.5 s) are more opaque AND thicker than minor
    // lines (every 0.1 s) so the eye can easily pick them out as a
    // coarse ruler against the finer 0.1 s subdivisions.
    , m_minorPen(QColor(255, 0, 0, 90), 1.0)
    , m_majorPen(QColor(255, 0, 0, 160), 1.6)
{
    for (const QString& label : m_labels) {
        ChannelSeries cs;
        cs.minor = new QLineSeries();
        cs.minor->setPen(m_minorPen);
        // OpenGL OFF: the GL backend ignores pen alpha and clamps width
        // to >= 1, so translucency / sub-pixel width would be invisible.
        // The overlay has only a few hundred points, so software rendering
        // is plenty fast.
        cs.minor->setUseOpenGL(false);

        cs.major = new QLineSeries();
        cs.major->setPen(m_majorPen);
        cs.major->setUseOpenGL(false);

        m_series.insert(label, cs);
    }
}

pulse_overlay::~pulse_overlay() {
    for (const ChannelSeries& cs : m_series) {
        // The host chart may still hold references; detach first so
        // ~QChart doesn't try to delete them again.
        if (cs.minor && cs.minor->chart()) cs.minor->chart()->removeSeries(cs.minor);
        if (cs.major && cs.major->chart()) cs.major->chart()->removeSeries(cs.major);
        delete cs.minor;
        delete cs.major;
    }
    m_series.clear();
}

QList<QLineSeries*> pulse_overlay::seriesForLabel(const QString& label) const {
    QList<QLineSeries*> out;
    auto it = m_series.constFind(label);
    if (it == m_series.constEnd()) return out;
    if (it->minor) out.append(it->minor);
    if (it->major) out.append(it->major);
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
            if (cs.minor && cs.minor->chart()) cs.minor->chart()->removeSeries(cs.minor);
            if (cs.major && cs.major->chart()) cs.major->chart()->removeSeries(cs.major);
        }
    }
    // On re-enable, refresh() (called by the GUI right after this) re-
    // attaches the series and rebuilds the geometry.
}

void pulse_overlay::refresh() {
    if (!m_enabled) return;

    // Build the (x-coords) for one channel's worth of minor and major lines
    // over [tStart, tEnd], then bake them into the two polylines (vertical
    // segments separated by NaN breakpoints, so Qt Charts renders each
    // segment as its own line within a single series).
    //
    // Minor and major are computed from independent loops with their own
    // spacings rather than "every 5th minor is a major" because that
    // approach is fragile when tStart isn't aligned to a major. Doing
    // them independently keeps each grid snapped to its own spacing.
    auto buildPoints = [](double tStart, double tEnd, double spacing,
        double yMin, double yMax) -> QList<QPointF> {
            QList<QPointF> pts;
            if (tEnd <= tStart || spacing <= 0.0) return pts;

            // Snap the first line to a multiple of `spacing` so the grid
            // doesn't shimmer as the window scrolls.
            const double firstX = std::ceil(tStart / spacing) * spacing;

            const int approxN = static_cast<int>((tEnd - tStart) / spacing) + 2;
            pts.reserve(approxN * 3);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (double x = firstX; x <= tEnd; x += spacing) {
                pts.append({ x, yMin });
                pts.append({ x, yMax });
                pts.append({ nan, nan });
            }
            return pts;
        };

    for (const QString& label : m_labels) {
        auto it = m_series.find(label);
        if (it == m_series.end()) continue;
        QLineSeries* minor = it->minor;
        QLineSeries* major = it->major;
        if (!minor || !major) continue;

        // The host chart wipes and re-creates its axes on every redraw,
        // so any axis attachments we had are stale.
        if (minor->chart()) minor->chart()->removeSeries(minor);
        if (major->chart()) major->chart()->removeSeries(major);

        if (!m_gui->isChannelActive(label)) continue;
        QChartView* cv = m_gui->chartViewForSignalLabel(label);
        if (!cv || !cv->chart()) continue;

        auto hAxes = cv->chart()->axes(Qt::Horizontal);
        auto vAxes = cv->chart()->axes(Qt::Vertical);
        if (hAxes.isEmpty() || vAxes.isEmpty()) continue;
        auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first());
        if (!yAxis) continue;

        const double yMin = yAxis->min();
        const double yMax = yAxis->max();
        const double tStart = m_gui->currentStartTime();
        const double tEnd = tStart + m_gui->windowDuration();

        minor->replace(buildPoints(tStart, tEnd, kMinorSpacingSec, yMin, yMax));
        major->replace(buildPoints(tStart, tEnd, kMajorSpacingSec, yMin, yMax));

        // Add minor first, then major, so the major lines render on top.
        cv->chart()->addSeries(minor);
        minor->attachAxis(hAxes.first());
        minor->attachAxis(yAxis);

        cv->chart()->addSeries(major);
        major->attachAxis(hAxes.first());
        major->attachAxis(yAxis);
    }
}