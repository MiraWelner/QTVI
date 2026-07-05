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

    const double snappedS = std::round(std::min(globalStart, globalEnd) * sr) / sr;
    const double snappedE = std::round(std::max(globalStart, globalEnd) * sr) / sr;

    m_noiseManager->addSegment(
        static_cast<int>(snappedS * sr), static_cast<int>(snappedE * sr),
        signalLabel.toStdString(), m_currentMarkingType.toStdString(), sr);
    m_genExc.noiseExc.append({ snappedS, snappedE });
    m_genExc.data_type.append(signalLabel);
    m_genExc.marking_type.append(m_currentMarkingType);

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
}
void noise_marking_gui::commitMarkingSpan(const QString& clickedLabel,
    double globalStart, double globalEnd)
{
    const QStringList chans = scopeChannels(clickedLabel);
    if (chans.isEmpty()) { clearDragPreview(); return; }

    if (annotation_types::isParamEdit(m_currentMarkingType)) {
        for (const QString& ch : chans) {            // clear any click-click start markers
            ChannelMarkingState& st = markStateFor(ch);
            st.phase = MarkPhase::Idle;
            clearStartMarker(st);
        }
        finalizeParamEdit(chans, globalStart, globalEnd);   // one dialog, all channels, redraws
        return;
    }

    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    for (const QString& ch : chans) {
        markStateFor(ch).globalStartTime = globalStart;     // same span on every scope channel
        finalizeMarking(chartViewForSignalLabel(ch), globalEnd - globalOffset, ch);
    }
    handle_data_plot();
}
void noise_marking_gui::cancelMarking(const QString& signalLabel) {
    clearDragPreview();
    ChannelMarkingState& state = markStateFor(signalLabel);
    state.phase = MarkPhase::Idle;
    clearStartMarker(state);
    updateMarkingButtons();
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

    const double lo = std::min(x0, x1), hi = std::max(x0, x1);
    const double yTop = yAxis->max(), yBot = yAxis->min();

    QAreaSeries*& prev = m_dragPreviews[cv];   // inserts nullptr if this chart has none yet
    if (!prev) {
        auto* upper = new QLineSeries(); auto* lower = new QLineSeries();
        upper->append(lo, yTop); upper->append(hi, yTop);
        lower->append(lo, yBot); lower->append(hi, yBot);
        prev = new QAreaSeries(upper, lower);
        prev->setBrush(color); prev->setPen(Qt::NoPen);
        cv->chart()->addSeries(prev);
        prev->attachAxis(hAxes.first()); prev->attachAxis(yAxis);
    }
    else {
        prev->upperSeries()->replace({ {lo, yTop}, {hi, yTop} });
        prev->lowerSeries()->replace({ {lo, yBot}, {hi, yBot} });
    }
}

void noise_marking_gui::clearDragPreview() {
    for (auto* p : m_dragPreviews) {
        if (!p) continue;
        if (p->chart()) p->chart()->removeSeries(p);
        delete p;   // QAreaSeries owns its upper/lower line series
    }
    m_dragPreviews.clear();
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
    restore("ABP", m_markState_abp); restore("ART", m_markState_art);
    restore("ART_PULM", m_markState_art_pulm);
}

// ============================================================================
// Event filter
// ============================================================================

bool noise_marking_gui::eventFilter(QObject* watched, QEvent* event) {
    auto* viewport = qobject_cast<QWidget*>(watched);
    if (!viewport) return QDialog::eventFilter(watched, event);
    auto* cv = qobject_cast<QChartView*>(viewport->parent());
    if (!cv || !cv->chart()) return QDialog::eventFilter(watched, event);

    if (event->type() == QEvent::MouseButtonDblClick) {
        auto* me = static_cast<QMouseEvent*>(event);
        if (me->button() == Qt::LeftButton && editParamOverrideAt(cv, me->pos()))
            return true;
        return QDialog::eventFilter(watched, event);
    }

    if (event->type() == QEvent::MouseButtonPress) {
        auto* me = static_cast<QMouseEvent*>(event);
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
            const QColor col = annotation_types::colorFor(m_currentMarkingType);
            for (const QString& ch : scopeChannels(m_dragSignalLabel))
                if (QChartView* ccv = chartViewForSignalLabel(ch))
                    updateDragPreview(ccv, startLocal, curX, col);
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
                const double startLocal =
                    markStateFor(m_dragSignalLabel).globalStartTime - globalOffset;

                if (endX != startLocal) {              // moved -> drag commit
                    commitMarkingSpan(label,
                        markStateFor(label).globalStartTime, endX + globalOffset);
                }
                else if (m_markArmed) {                               // click-click: first click
                    clearDragPreview();
                    for (const QString& ch : scopeChannels(label)) {
                        ChannelMarkingState& st = markStateFor(ch);
                        st.phase = MarkPhase::WaitingForStop;
                        if (QChartView* ccv = chartViewForSignalLabel(ch))
                            showStartMarker(ccv, st.globalStartTime - globalOffset, st,
                                colorForSignal(ch), nullptr);
                    }
                }
                else {
                    clearDragPreview();
                }
            }
            else {
                clearDragPreview();           // released elsewhere -> drop preview
            }
        }
    }
    return QDialog::eventFilter(watched, event);
}

bool noise_marking_gui::handleMousePress(QChartView* cv, QWidget* viewport, QMouseEvent* me)
{;
    if (me->button() != Qt::LeftButton) return false;

    double clickedX = cv->chart()->mapToValue(me->pos()).x();
    const double chunkDur = totalChunkDuration();
    clickedX = std::clamp(clickedX, 0.0, chunkDur);

    const QString label = signalLabelForChartView(cv);
    if (!label.isEmpty()) {
        if (!isChannelActive(label)) return false;       // drag needs no arming
        const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

        // Second click of an armed click-click: commit start..here across scope.
        if (m_markArmed && markStateFor(label).phase == MarkPhase::WaitingForStop) {
            commitMarkingSpan(label,
                markStateFor(label).globalStartTime, clickedX + globalOffset);
            return true;
        }

        // Otherwise begin a potential drag / first click: record the start on
        // every scope channel; the release decides drag-commit vs first click.
        for (const QString& ch : scopeChannels(label))
            markStateFor(ch).globalStartTime = clickedX + globalOffset;
        m_isDragging = true;
        m_dragStartPos = me->pos();
        m_dragSignalLabel = label;
        if (!m_draggedViewport) { m_draggedViewport = viewport; m_draggedViewport->grabMouse(); }
        return true;
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