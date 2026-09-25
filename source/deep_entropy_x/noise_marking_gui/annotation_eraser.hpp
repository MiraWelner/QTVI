/**
 * @file   annotation_eraser.hpp
 * @brief  Right-click-to-delete for annotated regions on the markable charts.
 *
 *         A right-click inside a highlighted annotation removes that annotation
 *         (the most-recently-added one when several overlap, so repeated clicks
 *         peel a stack off one at a time, topmost first). The marking list
 *         (m_genExc) is the source of truth -- and now the only store --
 *         exactly as undo/clear treat it: the matching entry is dropped, the
 *         now-stale logged beats in the erased span are removed (the redraw
 *         re-detects/re-logs whatever is actually there now), and the dialog is
 *         redrawn.
 */
#pragma once

#include "gui_handler.h"
#include "logging/user_mark_log.hpp"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QPoint>
#include <QString>
#include <QVector>
#include <string>

class annotation_eraser {
public:
    explicit annotation_eraser(noise_marking_gui* gui) : m_gui(gui) {}

    annotation_eraser(const annotation_eraser&) = delete;
    annotation_eraser& operator=(const annotation_eraser&) = delete;

    bool handleRightClick(QChartView* cv, const QPoint& viewportPos) {
        //delete region on right click
        if (!m_gui || !cv || !cv->chart()) return false;

        const QString label = m_gui->signalLabelForChartView(cv);
        if (label.isEmpty() || !m_gui->isChannelActive(label)) return false;
        if (!noise_marking_gui::isMarkableChannel(label)) return false;

        const double localX = cv->chart()->mapToValue(viewportPos).x();
        const double globalOffset =
            m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once;
        const double gt = localX + globalOffset;

        AllFileMarkings& exc = m_gui->m_genExc;
        const std::string labelStd = label.toStdString();
        int hit = -1;
        for (int i = exc.marks.size() - 1; i >= 0; --i) {
            if (exc.marks[i].channel != labelStd) continue;
            if (gt >= exc.marks[i].start && gt <= exc.marks[i].end) { hit = i; break; }
        }
        if (hit >= 0) {
            const double remStart = exc.marks[hit].start;
            const double remEnd = exc.marks[hit].end;
            exc.marks.removeAt(hit);
            if (m_gui->m_beatLog) {
                m_gui->m_beatLog->removeInRange(beat_log::channelForLabel(label), remStart, remEnd);
            }
            m_gui->rebuildParamIndex();
            m_gui->handle_data_plot();
            return true;
        }

        return false;
    }

private:
    noise_marking_gui* m_gui = nullptr;
};