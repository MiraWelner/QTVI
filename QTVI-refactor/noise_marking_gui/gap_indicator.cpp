/**
 * @file   gap_indicator.cpp
 * @brief  See gap_indicator.hpp.
 */
#include "gap_indicator.hpp"
#include "gui_handler.h"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QValueAxis>
#include <QColor>
#include <QFont>
#include <QPen>

#include <algorithm>

namespace {

    // ------------------------------------------------------------------
    // Tunables. Adjust freely.
    // ------------------------------------------------------------------

    // Any jump between consecutive raw timestamps longer than this
    // counts as a gap. 60 s is the practical floor for "device dropped
    // out for at least a minute" in CHAOS recordings.
    constexpr double kGapThresholdSec = 60.0;

    // Bracket-line style.
    const QColor  COLOR_GAP_LINE = QColor(120, 120, 120, 200);
    constexpr double kBracketWidthPx = 1.5;

    // Label style.
    const QColor  COLOR_GAP_LABEL = QColor(60, 60, 60, 255);
    const QFont   GAP_LABEL_FONT = QFont("Arial", 9, QFont::Bold);

    // True iff a raw vector contains real data (not the (-1, -1) sentinel
    // file_to_bin writes for missing channels).
    bool isRawUsable(const QVector<QPointF>& v) {
        return v.size() >= 2 && !(v.size() == 1 && v[0].x() == -1.0);
    }

}  // namespace

gap_indicator::gap_indicator(noise_marking_gui* gui) : m_gui(gui) {}

gap_indicator::~gap_indicator() {
    clearSeries();
}

void gap_indicator::rescan() {
    m_gaps.clear();
    if (!m_gui) return;

    // The recording's gap structure is shared across channels (a real
    // dropout affects all signals at once), so we only need to scan one
    // raw vector. Prefer ECG1, fall back to other markable channels.
    const QVector<QPointF>* src = nullptr;
    if (isRawUsable(m_gui->m_ecg1Raw)) src = &m_gui->m_ecg1Raw;
    else if (isRawUsable(m_gui->m_ecg2Raw)) src = &m_gui->m_ecg2Raw;
    else if (isRawUsable(m_gui->m_ecg3Raw)) src = &m_gui->m_ecg3Raw;
    else if (isRawUsable(m_gui->m_ppgRaw))  src = &m_gui->m_ppgRaw;
    else if (isRawUsable(m_gui->m_abpRaw))  src = &m_gui->m_abpRaw;
    if (!src) return;

    for (int i = 1; i < src->size(); ++i) {
        const double dt = (*src)[i].x() - (*src)[i - 1].x();
        if (dt > kGapThresholdSec) {
            m_gaps.append({ (*src)[i - 1].x(), (*src)[i].x() });
        }
    }
}

void gap_indicator::clearSeries() {
    for (auto* s : m_series) {
        if (!s) continue;
        if (s->chart()) s->chart()->removeSeries(s);
        delete s;
    }
    m_series.clear();
}

void gap_indicator::refresh() {
    if (!m_gui || m_gaps.isEmpty()) return;

    const double viewStart = m_gui->currentStartTime();
    const double viewEnd = viewStart + m_gui->windowDuration();

    for (const QString& label : noise_marking_gui::markableChannelLabels()) {
        if (!m_gui->isChannelActive(label)) continue;
        QChartView* cv = m_gui->chartViewForSignalLabel(label);
        if (!cv || !cv->chart()) continue;

        QChart* chart = cv->chart();
        auto hAxes = chart->axes(Qt::Horizontal);
        auto vAxes = chart->axes(Qt::Vertical);
        if (hAxes.isEmpty() || vAxes.isEmpty()) continue;
        auto* yAx = qobject_cast<QValueAxis*>(vAxes.first());
        if (!yAx) continue;

        const double yLo = yAx->min();
        const double yHi = yAx->max();
        const double yMid = 0.5 * (yLo + yHi);

        for (const auto& gap : m_gaps) {
            // Visible-window intersection. Outside-the-window gaps are
            // skipped entirely.
            if (gap.second < viewStart) continue;
            if (gap.first > viewEnd)   continue;

            auto addBracket = [&](double x) {
                if (x < viewStart || x > viewEnd) return;
                auto* line = new QLineSeries();
                line->setPen(QPen(COLOR_GAP_LINE, kBracketWidthPx, Qt::DashLine));
                // OpenGL OFF: GL backend ignores dash patterns and pen
                // alpha. The bracket has 2 points; software rendering
                // is trivially cheap.
                line->setUseOpenGL(false);
                line->append(x, yLo);
                line->append(x, yHi);
                chart->addSeries(line);
                line->attachAxis(hAxes.first());
                line->attachAxis(yAx);
                m_series.append(line);
                };
            addBracket(gap.first);
            addBracket(gap.second);

            // Label at the visible midpoint so it's always shown when
            // the gap touches the window, even if the gap is wider than
            // the window (in which case both brackets are off-screen).
            const double visStart = std::max(gap.first, viewStart);
            const double visEnd = std::min(gap.second, viewEnd);
            const double midX = 0.5 * (visStart + visEnd);

            // QScatterSeries with an invisible marker is the cheapest
            // way to place freely-positioned text in chart coordinates;
            // Qt Charts has no first-class annotation primitive.
            auto* labelSeries = new QScatterSeries();
            labelSeries->setMarkerSize(0.1);
            labelSeries->setColor(Qt::transparent);
            labelSeries->setBorderColor(Qt::transparent);
            labelSeries->setPointLabelsVisible(true);
            labelSeries->setPointLabelsColor(COLOR_GAP_LABEL);
            labelSeries->setPointLabelsFont(GAP_LABEL_FONT);
            labelSeries->setPointLabelsClipping(false);
            const double dur = gap.second - gap.first;
            labelSeries->setPointLabelsFormat(
                QString("%1 s gap").arg(dur, 0, 'f', 1));
            labelSeries->append(midX, yMid);
            chart->addSeries(labelSeries);
            labelSeries->attachAxis(hAxes.first());
            labelSeries->attachAxis(yAx);
            m_series.append(labelSeries);
        }
    }
}