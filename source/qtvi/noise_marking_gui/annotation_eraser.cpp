/**
 * @file   annotation_eraser.cpp
 * @brief  See annotation_eraser.h.
 */
#include "annotation_eraser.h"
#include "gui_handler.h"
#include "user_annotation_handler.h"
#include "logging/user_mark_log.hpp"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QPoint>
#include <memory>

annotation_eraser::annotation_eraser(noise_marking_gui* gui) : m_gui(gui) {}

bool annotation_eraser::handleRightClick(QChartView* cv, const QPoint& viewportPos) {
    if (!m_gui || !cv || !cv->chart()) return false;

    const QString label = m_gui->signalLabelForChartView(cv);
    if (label.isEmpty() || !m_gui->isChannelActive(label)) return false;

    const double localX = cv->chart()->mapToValue(viewportPos).x();
    const double globalOffset =
        m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once;
    const double gt = localX + globalOffset;

    // 1) Annotations first (the primary highlight).
    GenExcStruct& exc = m_gui->m_genExc;
    int hit = -1;
    for (int i = exc.noiseExc.size() - 1; i >= 0; --i) {
        if (exc.data_type[i] != label) continue;
        if (gt >= exc.noiseExc[i].first && gt <= exc.noiseExc[i].second) { hit = i; break; }
    }
    if (hit >= 0) {
        const double remStart = exc.noiseExc[hit].first;
        const double remEnd = exc.noiseExc[hit].second;

        exc.noiseExc.removeAt(hit);
        exc.data_type.removeAt(hit);
        exc.marking_type.removeAt(hit);

        m_gui->m_noiseManager = std::make_unique<annotation_handler>();
        for (int i = 0; i < exc.noiseExc.size(); ++i) {
            const double sr = m_gui->sampleRateForSignal(exc.data_type[i]);
            m_gui->m_noiseManager->addSegment(
                static_cast<int>(exc.noiseExc[i].first * sr),
                static_cast<int>(exc.noiseExc[i].second * sr),
                exc.data_type[i].toStdString(),
                exc.marking_type[i].toStdString(),
                sr);
        }
        if (m_gui->m_beatLog) m_gui->m_beatLog->removeInRange(beat_log::channelForLabel(label), remStart, remEnd);
        m_gui->handle_data_plot();
        return true;
    }

    // 2) No annotation under the cursor -> try the param-override rectangles.
    //    Removing one reverts detection in that span to the config defaults, so
    //    drop the stale logged beats there; the redraw re-detects/re-logs.
    auto eraseOverride = [&](QVector<noise_marking_gui::ParamOverride>& v) -> bool {
        for (int i = v.size() - 1; i >= 0; --i) {
            if (v[i].channel != label) continue;
            if (gt >= v[i].start && gt <= v[i].end) {
                const double s = v[i].start, e = v[i].end;
                v.removeAt(i);
                if (m_gui->m_beatLog) m_gui->m_beatLog->removeInRange(beat_log::channelForLabel(label), s, e);
                return true;
            }
        }
        return false;
        };
    if (eraseOverride(m_gui->m_thresholdOverrides) || eraseOverride(m_gui->m_blankingOverrides)|| eraseOverride(m_gui->m_invertOverrides)) {
        m_gui->handle_data_plot();
        return true;
    }
    return false;
}