/**
 * @file   annotation_eraser.cpp
 * @brief  See annotation_eraser.h.
 */
#include "annotation_eraser.h"
#include "gui_handler.h"
#include "user_annotation_handler.h"
#include "beat_log.hpp"

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QPoint>
#include <memory>

annotation_eraser::annotation_eraser(noise_marking_gui* gui) : m_gui(gui) {}

bool annotation_eraser::handleRightClick(QChartView* cv, const QPoint& viewportPos) {
    if (!m_gui || !cv || !cv->chart()) return false;

    const QString label = m_gui->signalLabelForChartView(cv);
    if (label.isEmpty() || !m_gui->isChannelActive(label)) return false;

    // Annotations are stored in GLOBAL seconds, so convert the click the same
    // way the marking code does: chunk-local x from the chart + chunk offset.
    const double localX = cv->chart()->mapToValue(viewportPos).x();
    const double globalOffset =
        m_gui->m_currentChunkIndex * noise_marking_gui::seconds_in_memory_at_once;
    const double gt = localX + globalOffset;

    GenExcStruct& exc = m_gui->m_genExc;

    // Most-recently-added annotation on this channel whose [start,end] contains
    // the click. Back-to-front so overlapping marks peel off topmost first.
    int hit = -1;
    for (int i = exc.noiseExc.size() - 1; i >= 0; --i) {
        if (exc.data_type[i] != label) continue;
        if (gt >= exc.noiseExc[i].first && gt <= exc.noiseExc[i].second) { hit = i; break; }
    }
    if (hit < 0) return false;   // not inside any annotation -> let the event pass

    const double remStart = exc.noiseExc[hit].first;
    const double remEnd = exc.noiseExc[hit].second;

    exc.noiseExc.removeAt(hit);
    exc.data_type.removeAt(hit);
    exc.marking_type.removeAt(hit);

    // Rebuild the annotation manager from the survivors -- the same
    // "exception list is source of truth" reconstruction undo/clear use.
    m_gui->m_noiseManager = std::make_unique<annotation_handler>(m_gui->m_ecgSR);
    for (int i = 0; i < exc.noiseExc.size(); ++i) {
        const double sr = m_gui->sampleRateForSignal(exc.data_type[i]);
        m_gui->m_noiseManager->addSegment(
            static_cast<int>(exc.noiseExc[i].first * sr),
            static_cast<int>(exc.noiseExc[i].second * sr),
            exc.data_type[i].toStdString(),
            exc.marking_type[i].toStdString());
    }

    // Drop stale logged beats in the erased span; the redraw re-detects and
    // re-logs (a removed noise span starts detecting again; a removed
    // arrhythmia tag re-logs its beats with markType 0).
    if (m_gui->m_beatLog) {
        auto beatCh = [](const QString& l) -> beat_log::ChannelIdx {
            if (l == "ECG1") return beat_log::ECG1;
            if (l == "ECG2") return beat_log::ECG2;
            if (l == "ECG3") return beat_log::ECG3;
            if (l == "PPG")  return beat_log::PPG;
            return beat_log::ABP;
            };
        m_gui->m_beatLog->removeInRange(beatCh(label), remStart, remEnd);
    }

    m_gui->handle_data_plot();
    return true;
}