/**
 * @file   user_marking_handler.cpp
 * @brief  Marking state machine, event filter, and mouse handling.
 *         Parallels user_annotation_handler (data) and
 *         user_control_handler (buttons) for the marking UI layer.
 */

#include "gui_handler.h"
#include "chart_utils.hpp"

#include <QtCharts/QValueAxis>
#include <QtCharts/QLineSeries>
#include <QMouseEvent>
#include <algorithm>

 // ============================================================================
 // Marking operations
 // ============================================================================

void noise_marking_gui::finalizeMarking(QChartView* /*cv*/, double endX,
    const QString& signalLabel)
{
    ChannelMarkingState& state = markStateFor(signalLabel);
    const double sr = sampleRateForSignal(signalLabel);
    const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
    const double globalEnd = endX + globalOffset;
    const double globalStart = state.globalStartTime;

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

    state.phase = MarkPhase::Idle;
    clearStartMarker(state);

    if (QPushButton* b = startButtonForSignal(signalLabel)) b->setStyleSheet("");
    if (QPushButton* b = stopButtonForSignal(signalLabel)) {
        b->setStyleSheet(""); b->setEnabled(false);
    }

    if (m_markAllActive) {
        bool allIdle = true;
        for (const QString& lbl : markableChannelLabels())
            if (isChannelActive(lbl) && markStateFor(lbl).phase != MarkPhase::Idle)
            {
                allIdle = false; break;
            }
        if (allIdle) {
            m_markAllActive = false;
            ui->start_all_mark->setStyleSheet("");
            ui->stop_all_mark->setStyleSheet("");
            ui->stop_all_mark->setEnabled(false);
        }
    }
    updateNoiseHighlights();
}

void noise_marking_gui::cancelMarking(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    state.phase = MarkPhase::Idle;
    clearStartMarker(state);
    if (QPushButton* b = startButtonForSignal(signalLabel)) b->setStyleSheet("");
    if (QPushButton* b = stopButtonForSignal(signalLabel)) {
        b->setStyleSheet(""); b->setEnabled(false);
    }
}

void noise_marking_gui::beginMarking(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    if (state.phase != MarkPhase::Idle) { cancelMarking(signalLabel); return; }
    m_currentMarkingType = ui->marking_type->currentText();
    state.phase = MarkPhase::WaitingForStart;
    if (QPushButton* b = startButtonForSignal(signalLabel))
        b->setStyleSheet("background-color: #f39c12; color: white;");
}

void noise_marking_gui::beginMarkingAll() {
    if (m_markAllActive) {
        for (const QString& label : markableChannelLabels())
            if (isChannelActive(label)) cancelMarking(label);
        m_markAllActive = false;
        ui->start_all_mark->setStyleSheet("");
        ui->stop_all_mark->setStyleSheet("");
        ui->stop_all_mark->setEnabled(false);
        return;
    }
    m_markAllActive = true;
    m_currentMarkingType = ui->marking_type->currentText();
    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::Idle) cancelMarking(label);
        state.phase = MarkPhase::WaitingForStart;
        if (QPushButton* b = startButtonForSignal(label))
            b->setStyleSheet("background-color: #f39c12; color: white;");
    }
    ui->start_all_mark->setStyleSheet("background-color: #f39c12; color: white;");
}

void noise_marking_gui::beginStopPhaseAll() {
    for (const QString& label : markableChannelLabels()) {
        if (!isChannelActive(label)) continue;
        ChannelMarkingState& state = markStateFor(label);
        if (state.phase != MarkPhase::WaitingForEnd) continue;
        state.phase = MarkPhase::WaitingForStop;
        if (QPushButton* b = stopButtonForSignal(label))
            b->setStyleSheet("background-color: #e74c3c; color: white;");
        if (QPushButton* b = startButtonForSignal(label)) b->setStyleSheet("");
    }
    ui->start_all_mark->setStyleSheet("");
    ui->stop_all_mark->setStyleSheet("background-color: #e74c3c; color: white;");
}

void noise_marking_gui::beginStopPhase(const QString& signalLabel) {
    ChannelMarkingState& state = markStateFor(signalLabel);
    if (state.phase != MarkPhase::WaitingForEnd) return;
    state.phase = MarkPhase::WaitingForStop;
    if (QPushButton* b = stopButtonForSignal(signalLabel))
        b->setStyleSheet("background-color: #e74c3c; color: white;");
    if (QPushButton* b = startButtonForSignal(signalLabel)) b->setStyleSheet("");
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

void noise_marking_gui::restoreMarkingMarkers() {
    const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
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
    restore("ECG3", m_markState_ecg3); restore("PPG", m_markState_ppg);
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

    if (event->type() == QEvent::MouseButtonPress)
        return handleMousePress(cv, viewport, static_cast<QMouseEvent*>(event))
        ? true : QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseMove && m_isDragging) return true;

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
                endX = std::clamp(endX, m_currentStartTime,
                    m_currentStartTime + m_windowDuration);
                ChannelMarkingState& state = markStateFor(label);
                const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
                const double localStart = state.globalStartTime - globalOffset;
                if (std::abs(endX - localStart) > 0.1)
                    finalizeMarking(cv, endX, label);
            }
            m_dragSignalLabel.clear();
            return true;
        }
    }
    return QDialog::eventFilter(watched, event);
}

bool noise_marking_gui::handleMousePress(QChartView* cv, QWidget* viewport,
    QMouseEvent* me)
{
    if (me->button() != Qt::LeftButton) return false;

    double clickedX = cv->chart()->mapToValue(me->pos()).x();
    const double chunkDur = totalChunkDuration();
    clickedX = std::clamp(clickedX, 0.0, chunkDur);

    const QString label = signalLabelForChartView(cv);
    if (!label.isEmpty()) {
        if (!isChannelActive(label)) return false;
        const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;

        if (m_markAllActive) {
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
                    st.phase = MarkPhase::WaitingForEnd;
                }
                ui->stop_all_mark->setEnabled(true);
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
            state.phase = MarkPhase::WaitingForEnd;
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

    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    const bool isNavChart =
        (cv == ui->ecg_ampogram_axis)
        || (cv == ui->ppg_ampogram_axis)
        || (cv == ui->hyp_accel_resp_axis && sleepPresent);

    if (isNavChart) {
        const double globalClickX = cv->chart()->mapToValue(me->pos()).x();
        const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
        const double localTarget = globalClickX - globalOffset;
        const double maxStart = std::max(0.0, chunkDur - m_windowDuration);
        m_currentStartTime = std::clamp(
            localTarget - m_windowDuration / 2.0, 0.0, maxStart);
        handle_data_plot();
        updateAmpogramCursor();
        return true;
    }
    return false;
}