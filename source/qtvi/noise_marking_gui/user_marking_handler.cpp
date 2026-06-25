/**
 * @file   user_marking_handler.cpp
 * @brief  Marking state machine, event filter, and mouse handling.
 *         Parallels user_annotation_handler (data) and
 *         user_control_handler (buttons) for the marking UI layer.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "beat_log.hpp"
#include "annotation_types.hpp"
#include "annotation_eraser.h"

#include <QtCharts/QAreaSeries>
#include <QtCharts/QLineSeries>
#include <QMouseEvent>
#include <algorithm>

 // ============================================================================
 // Marking operations
 // ============================================================================

void noise_marking_gui::finalizeMarking(QChartView* /*cv*/, double endX, const QString& signalLabel)
{
    clearDragPreview();
    ChannelMarkingState& state = markStateFor(signalLabel);
    const double sr = sampleRateForSignal(signalLabel);
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    const double globalEnd = endX + globalOffset;
    const double globalStart = state.globalStartTime;

    if (annotation_types::isParamEdit(m_currentMarkingType)) {
        state.phase = MarkPhase::Idle;
        clearStartMarker(state);
        updateEcgMarkButtonStyle();
        if (m_markAllMode != MarkAllMode::None) {
            bool allIdle = true;
            for (const QString& lbl : markableChannelLabels())
                if (isChannelActive(lbl) && markStateFor(lbl).phase != MarkPhase::Idle) {
                    allIdle = false;
                    break;
                }
            if (allIdle) { m_markAllMode = MarkAllMode::None; updateAllChannelButtonStates(); }
        }
        finalizeParamEdit(signalLabel, globalStart, globalEnd);
        return;
    }

    const double snappedS = std::round(std::min(globalStart, globalEnd) * sr) / sr;
    const double snappedE = std::round(std::max(globalStart, globalEnd) * sr) / sr;

    m_noiseManager->addSegment(
        static_cast<int>(snappedS * sr),
        static_cast<int>(snappedE * sr),
        signalLabel.toStdString(),
        m_currentMarkingType.toStdString());

    m_genExc.noiseExc.append({ snappedS, snappedE });
    m_genExc.data_type.append(signalLabel);
    m_genExc.marking_type.append(m_currentMarkingType);

    // Any mark may change this channel's beats in the span: noise suppresses
    // detection, and any annotation makes the beats re-detect against an
    // annotation-free reference and re-tag with the new marking_type. The log
    // is append-only, so drop the stale rows here; the handle_data_plot() below
    // re-logs whatever still detects (nothing in a noise span) with the right
    // tag. Without this, beats logged before the mark linger as "None".
    if (m_beatLog) {
        auto beatCh = [](const QString& l) -> beat_log::ChannelIdx {
            if (l == "ECG1") return beat_log::ECG1;
            if (l == "ECG2") return beat_log::ECG2;
            if (l == "ECG3") return beat_log::ECG3;
            if (l == "PPG_ACCEL")  return beat_log::PPG;
            return beat_log::ABP;
            };
        m_beatLog->removeInRange(beatCh(signalLabel), snappedS, snappedE);
    }

    state.phase = MarkPhase::Idle;
    clearStartMarker(state);
    updateEcgMarkButtonStyle();


    if (m_markAllMode != MarkAllMode::None) {
        bool allIdle = true;
        for (const QString& lbl : markableChannelLabels())
            if (isChannelActive(lbl) && markStateFor(lbl).phase != MarkPhase::Idle) {
                allIdle = false; 
                break;
            }
        if (allIdle) { m_markAllMode = MarkAllMode::None; updateAllChannelButtonStates(); }
    }
    handle_data_plot();
}

void noise_marking_gui::cancelMarking(const QString& signalLabel) {
    clearDragPreview();
    ChannelMarkingState& state = markStateFor(signalLabel);
    state.phase = MarkPhase::Idle;
    clearStartMarker(state);
    updateEcgMarkButtonStyle();
}

void noise_marking_gui::beginMarking(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    if (state.phase != MarkPhase::Idle) { cancelMarking(signalLabel); return; }
    m_currentMarkingType = ui->marking_type->currentText();
    state.phase = MarkPhase::WaitingForStart;
    updateEcgMarkButtonStyle();
}

void noise_marking_gui::beginStopPhase(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    if (state.phase != MarkPhase::WaitingForEnd) return;
    state.phase = MarkPhase::WaitingForStop;
    updateEcgMarkButtonStyle();
}

void noise_marking_gui::beginMarkingAll() {
    if (m_markAllMode == MarkAllMode::All) { 
        exitAllMarkModes(); 
        updateAllChannelButtonStates(); 
        return; 
    }
    exitAllMarkModes();
    m_markAllMode = MarkAllMode::All;
    m_currentMarkingType = ui->marking_type->currentText();
    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        markStateFor(label).phase = MarkPhase::WaitingForStart;
    }
    updateAllChannelButtonStates();
}


void noise_marking_gui::beginStopPhaseAll() {
    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::WaitingForEnd) continue;
        state.phase = MarkPhase::WaitingForStop;
    }
    updateAllChannelButtonStates();
}

void noise_marking_gui::beginMarkingEcgAll() {
    static const QStringList kEcg{ "ECG1", "ECG2", "ECG3" };
    if (m_markAllMode == MarkAllMode::Ecg) { 
        exitAllMarkModes(); 
        updateAllChannelButtonStates();
        return;
    }
    exitAllMarkModes();
    m_markAllMode = MarkAllMode::Ecg;
    m_currentMarkingType = ui->marking_type->currentText();
    for (const QString& label : kEcg) {
        if (!isChannelActive(label)) continue;
        markStateFor(label).phase = MarkPhase::WaitingForStart;
    }
    updateAllChannelButtonStates();
}

void noise_marking_gui::beginStopPhaseEcgAll() {
    static const QStringList kEcg{ "ECG1", "ECG2", "ECG3" };
    for (const QString& label : kEcg) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::WaitingForEnd) continue;
        state.phase = MarkPhase::WaitingForStop;
    }
    updateAllChannelButtonStates();
}

void noise_marking_gui::toggleMark(const QString& signalLabel) {
    switch (markStateFor(signalLabel).phase) {
    case MarkPhase::Idle:            beginMarking(signalLabel);   break;
    case MarkPhase::WaitingForStart: cancelMarking(signalLabel);  break;  // un-arm
    case MarkPhase::WaitingForEnd:   beginStopPhase(signalLabel); break;  // arm the end click
    case MarkPhase::WaitingForStop:  cancelMarking(signalLabel);  break;  // un-arm
    }
}

void noise_marking_gui::toggleMarkAll() {
    if (m_markAllMode != MarkAllMode::All) { beginMarkingAll(); return; }
    bool anyWaitingEnd = false;
    for (const QString& lbl : markableChannelLabels())
        if (isChannelActive(lbl) && markStateFor(lbl).phase == MarkPhase::WaitingForEnd) { anyWaitingEnd = true; break; }
    if (anyWaitingEnd) beginStopPhaseAll();
    else               beginMarkingAll();
}


void noise_marking_gui::toggleMarkEcgAll() {
    static const QStringList kEcg{ "ECG1", "ECG2", "ECG3" };
    if (m_markAllMode != MarkAllMode::Ecg) { beginMarkingEcgAll(); return; }
    bool anyWaitingEnd = false;
    for (const QString& lbl : kEcg)
        if (isChannelActive(lbl) && markStateFor(lbl).phase == MarkPhase::WaitingForEnd) { anyWaitingEnd = true; break; }
    if (anyWaitingEnd) beginStopPhaseEcgAll();
    else               beginMarkingEcgAll();
}

// ============================================================================
// Marker line helpers
// ============================================================================

void noise_marking_gui::clearStartMarker(ChannelMarkingState& state) {
    if (state.startMarkerLine && state.startMarkerLine->chart()) {
        state.startMarkerLine->chart()->removeSeries(state.startMarkerLine);
        delete state.startMarkerLine;
        state.startMarkerLine = nullptr;
    }
}

void noise_marking_gui::showStartMarker(QChartView* cv, double xValue,
    ChannelMarkingState& state, const QColor& color, QPushButton* stopBtn)
{
    clearStartMarker(state);
    state.startMarkerLine = new QLineSeries();
    state.startMarkerLine->setPen(QPen(color, 2, Qt::DashLine));

    auto* yAxis = qobject_cast<QValueAxis*>(cv->chart()->axes(Qt::Vertical).first());
    state.startMarkerLine->append(xValue, yAxis->min());
    state.startMarkerLine->append(xValue, yAxis->max());

    cv->chart()->addSeries(state.startMarkerLine);
    state.startMarkerLine->attachAxis(cv->chart()->axes(Qt::Horizontal).first());
    state.startMarkerLine->attachAxis(yAxis);
    if (stopBtn) stopBtn->setEnabled(true);
}

void noise_marking_gui::updateDragPreview(QChartView* cv, double x0, double x1,
    const QColor& color)
{
    if (!cv || !cv->chart()) return;
    auto vAxes = cv->chart()->axes(Qt::Vertical);
    auto hAxes = cv->chart()->axes(Qt::Horizontal);
    if (vAxes.isEmpty() || hAxes.isEmpty()) return;
    auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first());
    if (!yAxis) return;

    const double lo = std::min(x0, x1);
    const double hi = std::max(x0, x1);
    const double yTop = yAxis->max();
    const double yBot = yAxis->min();

    // Lazily create the rectangle on the first move of this drag; afterwards
    // just move its two edges. Translucent fill so it reads as a preview.
    if (!m_dragPreview) {
        auto* upper = new QLineSeries();
        auto* lower = new QLineSeries();
        upper->append(lo, yTop); upper->append(hi, yTop);
        lower->append(lo, yBot); lower->append(hi, yBot);
        m_dragPreview = new QAreaSeries(upper, lower);
        QColor fill = color;
        m_dragPreview->setBrush(fill);
        m_dragPreview->setPen(Qt::NoPen);
        cv->chart()->addSeries(m_dragPreview);
        m_dragPreview->attachAxis(hAxes.first());
        m_dragPreview->attachAxis(yAxis);
    }
    else {
        m_dragPreview->upperSeries()->replace({ {lo, yTop}, {hi, yTop} });
        m_dragPreview->lowerSeries()->replace({ {lo, yBot}, {hi, yBot} });
    }
}

void noise_marking_gui::clearDragPreview() {
    if (!m_dragPreview) return;
    if (m_dragPreview->chart()) m_dragPreview->chart()->removeSeries(m_dragPreview);
    // QAreaSeries owns its upper/lower line series, so deleting it frees them.
    delete m_dragPreview;
    m_dragPreview = nullptr;
    m_dragPreviewLabel.clear();
}

void noise_marking_gui::restoreMarkingMarkers() {
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    const double chunkEnd = globalOffset + seconds_in_memory_at_once;

    auto restore = [&](const QString& label, ChannelMarkingState& state) {
        if (state.phase != MarkPhase::WaitingForEnd
            && state.phase != MarkPhase::WaitingForStop) return;
        state.startMarkerLine = nullptr;
        QChartView* cv = chartViewForSignalLabel(label);
        if (!cv || !isChannelActive(label)) return;
        if (state.globalStartTime >= globalOffset
            && state.globalStartTime <= chunkEnd) {
            showStartMarker(cv, state.globalStartTime - globalOffset, state,
                colorForSignal(label), stopButtonForSignal(label));
        }
        };
    restore("ECG1", m_markState_ecg1); restore("ECG2", m_markState_ecg2);
    restore("ECG3", m_markState_ecg3); restore("PPG_ACCEL", m_markState_ppg);
    restore("ABP", m_markState_abp);
}

// ============================================================================
// Event filter
// ============================================================================

bool noise_marking_gui::eventFilter(QObject* watched, QEvent* event) {
    auto* viewport = qobject_cast<QWidget*>(watched);
    if (!viewport) return QDialog::eventFilter(watched, event);
    auto* cv = qobject_cast<QChartView*>(viewport->parent());
    if (!cv || !cv->chart()) return QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseButtonPress) {
        auto* me = static_cast<QMouseEvent*>(event);
        // Right-click deletes the annotation under the cursor; if the click
        // wasn't inside one, fall through to normal handling.
        if (me->button() == Qt::RightButton && m_annotationEraser
            && m_annotationEraser->handleRightClick(cv, me->pos()))
            return true;
        return handleMousePress(cv, viewport, me)
            ? true : QDialog::eventFilter(watched, event);
    }

    if (event->type() == QEvent::MouseMove && m_isDragging) {
        auto* me = static_cast<QMouseEvent*>(event);
        const QString moveLabel = signalLabelForChartView(cv);
        if (!moveLabel.isEmpty() && moveLabel == m_dragSignalLabel) {
            double curX = cv->chart()->mapToValue(me->position().toPoint()).x();
            curX = std::clamp(curX, current_start_time,
                current_start_time + visible_window_size);
            const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
            const double startLocal =
                markStateFor(m_dragSignalLabel).globalStartTime - globalOffset;
            m_dragPreviewLabel = m_dragSignalLabel;
            updateDragPreview(cv, startLocal, curX,
                annotation_types::colorFor(m_currentMarkingType));
        }
        return true;
    }

    if (event->type() == QEvent::MouseButtonRelease) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() == Qt::LeftButton && m_isDragging) {
            if (m_draggedViewport) {
                m_draggedViewport->releaseMouse();
                m_draggedViewport = nullptr;
            }
            m_isDragging = false;

            const QString label = signalLabelForChartView(cv);
            if (!label.isEmpty() && label == m_dragSignalLabel) {
                double endX = cv->chart()->mapToValue(me->pos()).x();
                endX = std::clamp(endX, current_start_time, current_start_time + visible_window_size);
                const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

                ChannelMarkingState& state = markStateFor(label);
                const double localStart = state.globalStartTime - globalOffset;
                if (std::abs(endX - localStart) > 0.1)
                    finalizeMarking(cv, endX, label);   // clears preview itself (line 25)
                else
                    clearDragPreview();   // too short to commit -> drop preview
            }
            else {
                clearDragPreview();           // released elsewhere -> drop preview
            }
        }
    }
    return QDialog::eventFilter(watched, event);
}

bool noise_marking_gui::handleMousePress(QChartView* cv, QWidget* viewport, QMouseEvent* me)
{
    if (me->button() != Qt::LeftButton) return false;

    double clickedX = cv->chart()->mapToValue(me->pos()).x();
    const double chunkDur = totalChunkDuration();
    clickedX = std::clamp(clickedX, 0.0, chunkDur);

    const QString label = signalLabelForChartView(cv);
    if (!label.isEmpty()) {
        if (!isChannelActive(label)) return false;
        const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

        if (single_ecg_marker_clicked) {
            ChannelMarkingState& st = markStateFor(label);
            st.globalStartTime = clickedX + globalOffset;
            showStartMarker(cv, clickedX, st, colorForSignal(label), stopButtonForSignal(label));
            st.phase = MarkPhase::WaitingForStop;
            single_ecg_marker_clicked = false;
            updateEcgMarkButtonStyle();
            return true;
        }

        if(m_markAllMode != MarkAllMode::None){
            bool anyStart = false, anyEnd = false, anyStop = false;
            for (const QString& lbl : markableChannelLabels()) {
                if (!isChannelActive(lbl)) continue;
                MarkPhase p = markStateFor(lbl).phase;
                if (p == MarkPhase::WaitingForStart) anyStart = true;
                if (p == MarkPhase::WaitingForEnd)   anyEnd = true;
                if (p == MarkPhase::WaitingForStop)  anyStop = true;
            }
            if (anyStop) {
                for (const QString& lbl : markableChannelLabels()) {
                    if (!isChannelActive(lbl)) continue;
                    if (markStateFor(lbl).phase != MarkPhase::WaitingForStop) continue;
                    finalizeMarking(chartViewForSignalLabel(lbl), clickedX, lbl);
                }
                return true;
            }
            if (anyStart || anyEnd) {
                for (const QString& lbl : markableChannelLabels()) {
                    if (!isChannelActive(lbl)) continue;
                    ChannelMarkingState& st = markStateFor(lbl);
                    if (st.phase != MarkPhase::WaitingForStart
                        && st.phase != MarkPhase::WaitingForEnd) continue;
                    st.globalStartTime = clickedX + globalOffset;
                    showStartMarker(chartViewForSignalLabel(lbl), clickedX, st,
                        colorForSignal(lbl), stopButtonForSignal(lbl));
                    st.phase = MarkPhase::WaitingForStop;
                }
                updateAllChannelButtonStates();
                return true;
            }
        }

        ChannelMarkingState& state = markStateFor(label);
        switch (state.phase) {
        case MarkPhase::WaitingForStart:
        case MarkPhase::WaitingForEnd:
            state.globalStartTime = clickedX + globalOffset;
            showStartMarker(cv, clickedX, state, colorForSignal(label),
                stopButtonForSignal(label));
            state.phase = MarkPhase::WaitingForStop;
            updateEcgMarkButtonStyle();
            return true;
        case MarkPhase::WaitingForStop:
            finalizeMarking(cv, clickedX, label);
            return true;
        case MarkPhase::Idle:
            m_isDragging = true;
            m_dragStartPos = me->pos();
            m_dragSignalLabel = label;
            state.globalStartTime = clickedX + globalOffset;
            if (!m_draggedViewport) {
                m_draggedViewport = viewport;
                m_draggedViewport->grabMouse();
            }
            return true;
        }
    }

    const bool sleepPresent = sleep_data_present(m_sleepStages);
    const bool isNavChart =
        (cv == ui->ecg_ampogram_axis)
        || (cv == ui->ppg_ampogram_axis)
        || (cv == ui->hyp_resp_axis && sleepPresent);

    if (isNavChart) {
        const double globalClickX = cv->chart()->mapToValue(me->pos()).x();
        const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
        const double localTarget = globalClickX - globalOffset;
        const double maxStart = std::max(0.0, chunkDur - visible_window_size);
        current_start_time = std::clamp(
            localTarget - visible_window_size / 2.0, 0.0, maxStart);
        handle_data_plot();
        updateAmpogramCursor();
        return true;
    }
    return false;
}