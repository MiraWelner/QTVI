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
#include "chart_utils.hpp"

namespace {
    const int min_length_for_gap = 30; //timestamp distance between the raw samples constitutes a 'gap' to be labled in the gui
    const QColor  gap_marker_color = QColor(255, 20, 147, 255);
    QString gap_lable_fmt(double sec) {
        if (sec < 60.0)        return QString("%1 sec gap").arg(sec, 0, 'f', 1);
        if (sec < 3600.0)      return QString("%1 min gap").arg(sec / 60.0, 0, 'f', 1);
        return                       QString("%1 hour gap").arg(sec / 3600.0, 0, 'f', 2);
    }

}  // namespace

gap_indicator::gap_indicator(noise_marking_gui* gui) : m_gui(gui) {}

gap_indicator::~gap_indicator() {
    clearSeries();
}

void gap_indicator::rescan() {
    m_gaps.clear();
    if (!m_gui) return;

    // Pick the source vector to scan. All channels share the same gap
    // structure (they come from the same .dat), so one is enough. Prefer
    // ECG1; fall back through the other markable channels.
    // We ALSO need that channel's native rate to convert sample indices
    // into chart-x positions. After the GUI's raw-timestamp rewrite, the
    // post-gap sample sits at chart x = idx / native_rate, regardless of
    // its real wall-clock stamp. That same chart x is correct for every
    // other channel too, since they're all rate-aligned in real time.
    int srcCh = -1;
    const QVector<QPointF>* src = nullptr;
    if (!is_missing_signal(m_gui->m_ecg1)) { src = &m_gui->m_ecg1Raw; srcCh = noise_marking_gui::CH_ECG1; }
    else if (!is_missing_signal(m_gui->m_ecg2)) { src = &m_gui->m_ecg2Raw; srcCh = noise_marking_gui::CH_ECG2; }
    else if (!is_missing_signal(m_gui->m_ecg3)) { src = &m_gui->m_ecg3Raw; srcCh = noise_marking_gui::CH_ECG3; }
    else if (!is_missing_signal(m_gui->m_ppg)) { src = &m_gui->m_ppgRaw;  srcCh = noise_marking_gui::CH_PPG; }
    else if (!is_missing_signal(m_gui->m_abp)) { src = &m_gui->m_abpRaw;  srcCh = noise_marking_gui::CH_ABP; }
    if (!src || srcCh < 0) return;

    const float nativeHz = m_gui->channel_native_rates[srcCh];
    if (nativeHz <= 0.0f) return;
    const double dt = 1.0 / nativeHz;

    // Scan the REAL timestamps for jumps. Caller must have arranged that
    // this runs BEFORE the timestamp rewrite (gui_handler.cpp does so).
    for (int i = 1; i < src->size(); ++i) {
        const double realDt = (*src)[i].x() - (*src)[i - 1].x();
        if (realDt > min_length_for_gap) {
            // Bracket sits at the index-time position of the post-gap
            // sample. Sub-pixel-wide visually; effectively one vertical
            // line between the last pre-gap sample and the first post-
            // gap sample.
            GapEntry e;
            e.xPos = i * dt;
            e.durationSec = realDt;
            m_gaps.append(e);
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
            if (gap.xPos < viewStart || gap.xPos > viewEnd) continue;

            // One dashed vertical line at the gap position.
            auto* line = new QLineSeries();
            line->setPen(QPen(gap_marker_color, 1, Qt::DotLine));
            line->setUseOpenGL(true);
            line->append(gap.xPos, yLo);
            line->append(gap.xPos, yHi);
            chart->addSeries(line);
            line->attachAxis(hAxes.first());
            line->attachAxis(yAx);
            m_series.append(line);

            // Duration label. QScatterSeries with an invisible marker is
            // the cheapest way to place freely-positioned text in chart
            // coordinates; Qt Charts has no first-class annotation
            // primitive. Place at the gap x, at vertical mid-axis.
            auto* labelSeries = new QScatterSeries();
            labelSeries->setMarkerSize(0.1);
            labelSeries->setColor(Qt::transparent);
            labelSeries->setBorderColor(Qt::transparent);
            labelSeries->setPointLabelsVisible(true);
            labelSeries->setPointLabelsColor(gap_marker_color);
            labelSeries->setPointLabelsFont(cv->font());
            labelSeries->setPointLabelsClipping(false);
            labelSeries->setPointLabelsFormat(gap_lable_fmt(gap.durationSec));
            labelSeries->append(gap.xPos, yMid);
            chart->addSeries(labelSeries);
            labelSeries->attachAxis(hAxes.first());
            labelSeries->attachAxis(yAx);
            m_series.append(labelSeries);
        }
    }
}