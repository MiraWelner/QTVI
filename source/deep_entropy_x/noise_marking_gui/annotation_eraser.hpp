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
 *
 *         Friend of noise_marking_gui so it can edit m_genExc and trigger a
 *         redraw, mirroring gap_indicator's access.
 *
 *         HEADER-ONLY, and it includes gui_handler.h rather than
 *         forward-declaring: the whole body dereferences m_gui, so it needs the
 *         complete type anyway. No cycle, because gui_handler.h forward-declares
 *         annotation_eraser and holds it behind a unique_ptr; it never includes
 *         this file. Both translation units that use the eraser
 *         (gui_handler.cpp, user_marking_handler.cpp) already include
 *         gui_handler.h, so this costs them nothing.
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

    /// Delete the annotation under a right-click. @p viewportPos is the click
    /// position in the chart view's viewport coordinates (as delivered to
    /// eventFilter). Returns true if an annotation was erased (the caller
    /// should then treat the mouse event as handled).
    bool handleRightClick(QChartView* cv, const QPoint& viewportPos) {
        if (!m_gui || !cv || !cv->chart()) return false;

        const QString label = m_gui->signalLabelForChartView(cv);
        if (label.isEmpty() || !m_gui->isChannelActive(label)) return false;

        const double localX = cv->chart()->mapToValue(viewportPos).x();
        const double globalOffset =
            m_gui->current_chunk_index * noise_marking_gui::seconds_in_memory_at_once;
        const double gt = localX + globalOffset;

        // 1) Annotations first (the primary highlight).
        GenExcStruct& exc = m_gui->m_genExc;
        const std::string labelStd = label.toStdString();
        int hit = -1;
        for (int i = exc.marks.size() - 1; i >= 0; --i) {
            if (exc.marks[i].channel != labelStd) continue;
            if (gt >= exc.marks[i].start && gt <= exc.marks[i].end) { hit = i; break; }
        }
        if (hit >= 0) {
            const double remStart = exc.marks[hit].start;
            const double remEnd = exc.marks[hit].end;

            // One removeAt, and nothing to rebuild. This used to remove from
            // three of five parallel vectors -- leaving threshold and blanking
            // one entry long, which made consistent() false and silently
            // stopped every override being applied thereafter -- and then
            // reconstruct the whole sample-indexed second store from scratch.
            // Both problems were the second store.
            exc.marks.removeAt(hit);
            if (m_gui->m_beatLog)
                m_gui->m_beatLog->removeInRange(
                    beat_log::channelForLabel(label), remStart, remEnd);
            m_gui->handle_data_plot();
            return true;
        }

        // 2) No annotation under the cursor -> try the param-override rectangles.
        //    Removing one reverts detection in that span to the config defaults,
        //    so drop the stale logged beats there; the redraw re-detects/re-logs.
        auto eraseOverride = [&](QVector<noise_marking_gui::ParamOverride>& v) -> bool {
            for (int i = v.size() - 1; i >= 0; --i) {
                if (v[i].channel != label) continue;
                if (gt >= v[i].start && gt <= v[i].end) {
                    const double s = v[i].start, e = v[i].end;
                    v.removeAt(i);
                    if (m_gui->m_beatLog)
                        m_gui->m_beatLog->removeInRange(
                            beat_log::channelForLabel(label), s, e);
                    return true;
                }
            }
            return false;
            };
        if (eraseOverride(m_gui->m_thresholdOverrides)
            || eraseOverride(m_gui->m_blankingOverrides)
            || eraseOverride(m_gui->m_invertOverrides)) {
            m_gui->handle_data_plot();
            return true;
        }
        return false;
    }

private:
    noise_marking_gui* m_gui = nullptr;
};