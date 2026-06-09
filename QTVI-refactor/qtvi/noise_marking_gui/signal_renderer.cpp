/**
 * @file   signal_renderer.cpp
 * @brief  Signal chart rendering for the noise-marking GUI:
 *           - setupHypnogram        (sleep-stage overview)
 *           - handle_ampogram_plot  (amplitude variability overview)
 *           - renderWindowedChart   (core windowed signal renderer)
 *           - handle_data_plot      (main per-window redraw)
 *           - updateAmpogramCursor
 *           - updateNoiseHighlights
 *           - peaksForWindow / peaksForBpmWindow
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "simple_peak_finder.hpp"
#include "theme/theme.h"

#include <QtCharts/QAreaSeries>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QCategoryAxis>
#include <QFont>
#include <QSet>
#include <algorithm>
#include <cmath>

struct markable_data_series {
    //A series such as ECG, PPG, that can be annotated
    const QVector<double>* data;
    QColor color;
    const QVector<QPointF>* rawData = nullptr;
};

namespace {

    QCategoryAxis* make_time_labled_xaxis(double startLocal, double duration,
        double globalOffset, bool labelsVisible)
    {
        auto* xAxis = new QCategoryAxis();
        xAxis->setRange(startLocal, startLocal + duration);
		xAxis->append(QString::fromStdString("(HH:MM:SS)"), 0);
        const double fracs[] = { 1.0 / 6.0, 0.5, 5.0 / 6.0 };
        for (double f : fracs) {
            double localVal = startLocal + f * duration;
            xAxis->append(formatHMS(globalOffset + localVal), localVal);
        }
        xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        xAxis->setTruncateLabels(false);
        xAxis->setGridLineVisible(false);
        xAxis->setLabelsVisible(labelsVisible);
        if (!labelsVisible) xAxis->setVisible(false);
        return xAxis;
    }

    std::pair<double, double> renderWindowedChart(
        QChartView* view,
        const QList<markable_data_series>& serieses,
        QList<QLineSeries*>& persistentLines,
        double currentStartTime, double windowDuration,
        double globalOffset, double ecgSR,
        bool labelsVisible,
        bool useScatterMode,
        bool forceLineForUpsampled,
        double yScale = 1.0,
        bool compactTime = false)
    {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();
        chart->legend()->hide();
        chart->setMargins(QMargins(0, 0, 0, 4));

        QSet<QAbstractSeries*> persistentSet;
        for (auto* ln : persistentLines) persistentSet.insert(ln);

        for (auto* s : chart->series()) {
            chart->removeSeries(s);
            if (!persistentSet.contains(s)) delete s;
        }
        for (auto* a : chart->axes()) { chart->removeAxis(a); delete a; }

        QCategoryAxis* xAxis = make_time_labled_xaxis(
            currentStartTime, windowDuration, globalOffset, labelsVisible);
        chart->addAxis(xAxis, Qt::AlignBottom);

        //Create a Y axis for every plot
        auto* yAxis = new QValueAxis();
        yAxis->setLineVisible(true);
        yAxis->setGridLineVisible(false);
        yAxis->setVisible(true);
        yAxis->setLabelsVisible(true);
        yAxis->setLabelFormat("%.1e"); 
        yAxis->setTickCount(2);
        chart->addAxis(yAxis, Qt::AlignLeft);
        const QColor axisGray = yAxis->linePenColor();
        yAxis->setLabelsColor(axisGray);
        xAxis->setLabelsColor(axisGray);

        double gMin = 1e9, gMax = -1e9;

        struct PendingRaw { const QVector<QPointF>* rawData; QColor color; double center; };
        QList<PendingRaw> rawsToAdd;

        for (int slot = 0; slot < serieses.size(); ++slot) {
            const auto& d = serieses[slot];
            if (!d.data || isMissingSignal(*d.data)) continue;
            const bool hasRaw = d.rawData && isRawUsable(*d.rawData);

            auto ensurePersistentLine = [&]() -> QLineSeries* {
                while (persistentLines.size() <= slot)
                    persistentLines.append(nullptr);
                if (!persistentLines[slot]) {
                    auto* ln = new QLineSeries();
                    ln->setUseOpenGL(true);
                    persistentLines[slot] = ln;
                }
                persistentLines[slot]->setPen(QPen(d.color, 1));
                return persistentLines[slot];
                };

            QXYSeries* plotSeries = nullptr;
            if (forceLineForUpsampled || !useScatterMode) {
                QLineSeries* ln = ensurePersistentLine();
                ln->setVisible(true); chart->addSeries(ln); plotSeries = ln;
            }
            else if (!hasRaw) {
                auto* sc = new QScatterSeries();
                sc->setColor(d.color); sc->setBorderColor(Qt::transparent);
                sc->setMarkerSize(2.0); sc->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                sc->setUseOpenGL(true); chart->addSeries(sc); plotSeries = sc;
            }
            else {
                if (slot < persistentLines.size() && persistentLines[slot])
                    persistentLines[slot]->setVisible(false);
            }

            int startIdx = std::clamp(static_cast<int>(currentStartTime * ecgSR),
                0, static_cast<int>(d.data->size() - 1));
            int endIdx = std::clamp(static_cast<int>((currentStartTime + windowDuration) * ecgSR),
                0, static_cast<int>(d.data->size()));

            //for scaling, use median as the center, so the scaling is robust to noise
            std::vector<double> winVals;
            winVals.reserve(endIdx - startIdx);
            for (int i = startIdx; i < endIdx; ++i)
                winVals.push_back((*d.data)[i]);
            for (const QPointF& p : *d.rawData) {
                if (p.x() < currentStartTime) continue;
                if (p.x() > currentStartTime + windowDuration) break;
                winVals.push_back(p.y());
                if (p.y() < gMin) gMin = p.y();
                if (p.y() > gMax) gMax = p.y();
            }
            double center = 0.0;
            if (!winVals.empty()) {
                const auto mid = winVals.begin() + winVals.size() / 2;
                std::nth_element(winVals.begin(), mid, winVals.end());
                center = *mid;
            }

            QList<QPointF> pts;
            if (plotSeries) pts.reserve(endIdx - startIdx);
            for (int i = startIdx; i < endIdx; ++i) {
                const double raw = (*d.data)[i];
                if (raw < gMin) gMin = raw;
                if (raw > gMax) gMax = raw;
                if (plotSeries) {
                    const double scaled = (raw - center) * yScale + center;
                    pts.append({ static_cast<double>(i) / ecgSR, scaled });
                }
            }
            if (plotSeries) {
                plotSeries->replace(pts);
                plotSeries->attachAxis(xAxis);
                plotSeries->attachAxis(yAxis);
            }
            if (hasRaw) rawsToAdd.append({ d.rawData, d.color, center });
        }

        QList<QScatterSeries*> scattersForTopReorder;
        for (const auto& r : rawsToAdd) {
            auto* rawScatter = new QScatterSeries();
            rawScatter->setColor(Qt::black);
            rawScatter->setBorderColor(Qt::transparent);
            rawScatter->setMarkerSize(3.0);
            rawScatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
            rawScatter->setUseOpenGL(true);
            chart->addSeries(rawScatter);
            scattersForTopReorder.append(rawScatter);

            QList<QPointF> rawPts;
            rawPts.reserve(std::min<int>(r.rawData->size(), 4096));
            for (const QPointF& p : *r.rawData) {
                if (p.x() < currentStartTime) continue;
                if (p.x() > currentStartTime + windowDuration) break;
                rawPts.append({ p.x(), (p.y() - r.center) * yScale + r.center });
            }
            rawScatter->replace(rawPts);
            rawScatter->attachAxis(xAxis);
            rawScatter->attachAxis(yAxis);
        }

        for (auto* scatter : scattersForTopReorder) {
            chart->removeSeries(scatter);
            chart->addSeries(scatter);
            scatter->attachAxis(xAxis);
            scatter->attachAxis(yAxis);
        }

        setPaddedYRange(yAxis, gMin, gMax);
        return { gMin, gMax };
    }

} // namespace

QVector<QPointF> noise_marking_gui::display_peaks_in_window(const QString& label) const {
    /*
        If the user has clicked the peak finding check box, this function shows the peaks
    */
    if (!m_showPeaks) return {};
    data_channel_features r = channelRefs(label);
    auto params = simple_peak_finder::paramsFor(label);
    params.threshold = r.threshold_box->value();
    params.blanking_period = r.blanking_period_box->value();
    if (label == "PPG" || label == "ABP")
        return simple_peak_finder::findPeaksDerivative(
            *r.dataRaw, m_currentStartTime,
            m_currentStartTime + m_windowDuration, params);
    return simple_peak_finder::findPeaks(
        *r.dataRaw, m_currentStartTime,
        m_currentStartTime + m_windowDuration, params);
}

QVector<QPointF> noise_marking_gui::get_bpm(const QString& label, double& outDuration) const
{
    /*
        If the user has clicked the peak finding check box, disply the beats per minute.
        Different from display_peaks_in_window because it never takes data from a a shorter 
        window than 10s even if the window being displayed is shorter than 10s
    */
    outDuration = 0.0;
    if (!m_showPeaks) return {};
    data_channel_features r = channelRefs(label);
    constexpr double kMinBpmWindowSec = 10.0;
    const double tVisEnd = m_currentStartTime + m_windowDuration;
    double tStart = m_currentStartTime;
    if (m_windowDuration < kMinBpmWindowSec)
        tStart = std::max(0.0, tVisEnd - kMinBpmWindowSec);
    outDuration = tVisEnd - tStart;
    auto params = simple_peak_finder::paramsFor(label);
    params.threshold = m_cfg.height_threshold_percent;
    params.blanking_period = m_cfg.blanking_period;
    if (label == "PPG" || label == "ABP")
        return simple_peak_finder::findPeaksDerivative(*r.dataRaw, tStart, tVisEnd, params);
    return simple_peak_finder::findPeaks(*r.dataRaw, tStart, tVisEnd, params);
}

void noise_marking_gui::setupHypnogram() {
    if (m_sleepSR <= 0.0 || !sleepDataPresent(m_sleepStages)) return;

    auto* chart = ui->hyp_accel_resp_axis->chart();
    if (m_cvpCursorBar && m_cvpCursorBar->chart() == chart)
        chart->removeSeries(m_cvpCursorBar);
    for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
    m_hypnoStageSeries.clear();

    struct Stage { int value; QColor color; };
    const QList<Stage> stages = {
        {0, Qt::black}, {1, Qt::darkGreen}, {2, Qt::blue}, {3, Qt::cyan}, {4, Qt::red}
    };

    const double dt = 1.0 / m_sleepSR;
    const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;

    for (const auto& st : stages) {
        auto* s = new QScatterSeries();
        s->setColor(st.color); s->setMarkerSize(3.0);
        s->setPen(Qt::NoPen);
        s->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
        for (int i = 0; i < m_sleepStages.size(); ++i) {
            if (static_cast<int>(m_sleepStages[i]) == st.value)
                s->append(globalOffset + i * dt + dt / 2.0, st.value);
        }
        chart->addSeries(s);
        m_hypnoStageSeries.append(s);
    }

    for (QAbstractAxis* axis : chart->axes()) {
        chart->removeAxis(axis);
    }

    chart->setMargins(QMargins(0, 0, 0, 5));
    chart->setTitle("Sleep stages");
    chart->setTitleFont(Theme::chartTitleFont());
    chart->setTitleBrush(Qt::black);

    auto* xAxis = new QCategoryAxis();
    xAxis->setRange(globalOffset, globalOffset + seconds_in_memory_at_once);
    const int startHour = static_cast<int>(globalOffset / 3600.0);
    for (int h = 0; h <= 8; h += 2)
        xAxis->append(QString::number(startHour + h) + 'h', globalOffset + h * 3600.0);
    xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    xAxis->setTruncateLabels(false);
    xAxis->setGridLineVisible(false);
    xAxis->setLabelsFont(Theme::chartAxisFont());
    xAxis->setLabelsVisible(true);

    auto* yAxis = new QCategoryAxis();
    for (const auto& st : stages) yAxis->append("", st.value + 0.4);
    yAxis->setRange(-0.5, 4.5); yAxis->setReverse(true);
    yAxis->setVisible(false); yAxis->setGridLineVisible(false);

    chart->addAxis(xAxis, Qt::AlignBottom);
    chart->addAxis(yAxis, Qt::AlignLeft);
    for (auto* s : chart->series()) { s->attachAxis(xAxis); s->attachAxis(yAxis); }

    if (m_hypnoCursorBar) {
        chart->removeSeries(m_hypnoCursorBar);
        chart->addSeries(m_hypnoCursorBar);
        m_hypnoCursorBar->attachAxis(xAxis);
        m_hypnoCursorBar->attachAxis(yAxis);
    }
}

void noise_marking_gui::ampogram(double range) {
    /*
        creates the ampograms on the top right, which show the difference between min and max in a given range across the current
        8 hour period
    */
    const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;

    auto calculate_amplitude = [range, globalOffset](
        const QVector<double>& data, double sr)
        {
            QList<QPointF> pts;
            if (data.isEmpty() || sr <= 0.0) return pts;
            double duration = data.size() / sr;
            for (double t = 0; t <= duration - range; t += range) {
                int s = static_cast<int>(t * sr);
                int e = static_cast<int>((t + range) * sr);
                auto [mi, ma] = std::minmax_element(data.begin() + s, data.begin() + e);
                pts.append({ globalOffset + t, *ma - *mi });
            }
            return pts;
        };

    auto create_plot = [globalOffset](
        QChartView* view, QLineSeries* series, const QList<QPointF>& pts,
        QLineSeries* cursor, const QColor& color,
        const QString& title, bool showLabels)
        {
            series->replace(pts);
            series->setPen(QPen(color, 1));

            auto* chart = view->chart();
            for (QAbstractAxis* axis : chart->axes()) {
                chart->removeAxis(axis);
            }
            chart->legend()->hide();
            chart->setMargins(QMargins(0, 0, 10, 5));
            if (!title.isEmpty()) {
                chart->setTitle(title);
                chart->setTitleFont(Theme::chartTitleFont());
            }
            else {
                chart->setTitle(QString());
            }

            //add the x axis with hour lables
            auto* x_axis = new QCategoryAxis();
            x_axis->setRange(globalOffset, globalOffset + seconds_in_memory_at_once);
			x_axis->append(QString::fromStdString("(h)"), globalOffset+1800);
            const int startHour = static_cast<int>(globalOffset / 3600);
            for (int h = 1; h <= seconds_in_memory_at_once/ 3600; ++h) {
                x_axis->append(QString::number(startHour + h), globalOffset + h * 3600);
            }
            x_axis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
            x_axis->setTruncateLabels(false);
            x_axis->setGridLineVisible(false);
            x_axis->setLabelsVisible(showLabels);
            chart->addAxis(x_axis, Qt::AlignBottom);
            series->attachAxis(x_axis);
            cursor->attachAxis(x_axis);
            
			//add the y axis with two ticks and a bit of padding
            double yMin = 0, yMax = 1.0;
            if (!pts.isEmpty()) {
                auto [mi, ma] = std::minmax_element(pts.begin(), pts.end(),
                    [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
                double pad = std::max(0.5, (ma->y() - mi->y()) * 0.05);
                yMin = mi->y() - pad; yMax = ma->y() + pad;
            }
            auto* yAxis = new QValueAxis();
            yAxis->setRange(yMin, yMax);
            yAxis->setVisible(true);
            yAxis->setLabelsVisible(true);
            yAxis->setLineVisible(true);
            yAxis->setGridLineVisible(false);
            yAxis->setLabelFormat("%.1e");
            yAxis->setTickCount(2);
            chart->addAxis(yAxis, Qt::AlignLeft);
            const QColor axisGray = yAxis->linePenColor();
            yAxis->setLabelsColor(axisGray);
            x_axis->setLabelsColor(axisGray);
            series->attachAxis(yAxis);
            if (cursor) cursor->attachAxis(yAxis);
        };

    auto ecg1Pts = calculate_amplitude(m_ecg1, m_ecgSR);
    auto ecg2Pts = calculate_amplitude(m_ecg2, m_ecgSR);
    auto ecg3Pts = calculate_amplitude(m_ecg3, m_ecgSR);

    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    const bool ppgAmpHasLabels = !sleepPresent
        && ui->ppg_ampogram_axis && !ui->ppg_ampogram_axis->isHidden();
    const bool ecgAmpHasLabels = !sleepPresent && !ppgAmpHasLabels;

    auto* chart = ui->ecg_ampogram_axis->chart();
    if (ecg2_ampogram_series->chart() == chart) chart->removeSeries(ecg2_ampogram_series);
    if (ecg3_ampogram_series->chart() == chart) chart->removeSeries(ecg3_ampogram_series);

    create_plot(ui->ecg_ampogram_axis, ecg1_ampogram_series,
        ecg1Pts, m_ecgCursorBar, COLOR_ECG1, "ECG Amp-O-Gram", ecgAmpHasLabels);

    auto* xAxis = chart->axes(Qt::Horizontal).first();
    auto* yAxis = qobject_cast<QValueAxis*>(chart->axes(Qt::Vertical).first());

    auto attachExtra = [&](QLineSeries* series, const QList<QPointF>& pts, const QColor& color) {
        if (pts.isEmpty()) return;
        series->replace(pts); series->setPen(QPen(color, 1));
        chart->addSeries(series);
        series->attachAxis(xAxis); series->attachAxis(yAxis);
        };
    attachExtra(ecg2_ampogram_series, ecg2Pts, COLOR_ECG2);
    attachExtra(ecg3_ampogram_series, ecg3Pts, COLOR_ECG3);

    auto allPts = ecg1Pts + ecg2Pts + ecg3Pts;
    if (yAxis && !allPts.isEmpty()) {
        auto [mi, ma] = std::minmax_element(allPts.begin(), allPts.end(),
            [](const QPointF& a, const QPointF& b) { return a.y() < b.y(); });
        setPaddedYRange(yAxis, mi->y(), ma->y());
    }

    create_plot(ui->ppg_ampogram_axis, ppg_ampogram_series,
        calculate_amplitude(m_ppg, m_ppgSR),
        m_ppgCursorBar, COLOR_PPG, "PPG Amp-O-Gram", ppgAmpHasLabels);
}

// ============================================================================
// Main data plot
// ============================================================================

void noise_marking_gui::handle_data_plot() {
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area;
    }
    m_highlights.clear();
    if (m_gapIndicator) m_gapIndicator->clearSeries();

    m_markState_ecg1.startMarkerLine = nullptr;
    m_markState_ecg2.startMarkerLine = nullptr;
    m_markState_ecg3.startMarkerLine = nullptr;
    m_markState_ppg.startMarkerLine = nullptr;
    m_markState_abp.startMarkerLine = nullptr;

    QChartView* xLabelOwnerRight = nullptr;
    for (auto* cv : { ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
                      ui->ppg_axis, ui->accel_or_abp_axis }) {
        if (cv && !cv->isHidden()) xLabelOwnerRight = cv;
    }

    auto keepFor = [&](QChartView* cv) -> QList<QAbstractSeries*> {
        QList<QAbstractSeries*> keep;
        if (!cv) return keep;
        auto it = m_persistentLines.constFind(cv);
        if (it != m_persistentLines.constEnd())
            for (auto* ln : it.value()) if (ln) keep.append(ln);
        if (m_pulseOverlay) {
            const QString label = signalLabelForChartView(cv);
            if (!label.isEmpty())
                for (QLineSeries* s : m_pulseOverlay->seriesForLabel(label))
                    if (s) keep.append(s);
        }
        return keep;
        };

    wipeChartContent(ui->ecg_axis_1->chart(), keepFor(ui->ecg_axis_1));
    wipeChartContent(ui->ecg_axis_2->chart(), keepFor(ui->ecg_axis_2));
    wipeChartContent(ui->ecg_axis_3->chart(), keepFor(ui->ecg_axis_3));
    wipeChartContent(ui->ppg_axis->chart(), keepFor(ui->ppg_axis));
    if (ui->accel_or_abp_axis && !isMissingSignal(m_abp))
        wipeChartContent(ui->accel_or_abp_axis->chart(),
            keepFor(ui->accel_or_abp_axis));

    const bool sleepPresent = sleepDataPresent(m_sleepStages);
    if (!sleepPresent && ui->hyp_accel_resp_axis)
        wipeChartContent(ui->hyp_accel_resp_axis->chart(),
            keepFor(ui->hyp_accel_resp_axis));
    if (ui->cvp_axis)
        wipeChartContent(ui->cvp_axis->chart(), keepFor(ui->cvp_axis));

    auto plotMarkable = [&](const QString& label) {
        if (!isChannelActive(label)) return;
        const data_channel_features r = channelRefs(label);
        if (!r.chartView || !r.upsampled_data || !r.state) return;

        const double sr = r.sampleRate ? *r.sampleRate : m_ecgSR;
        const QVector<QPointF> emptyRaw;
        const QVector<QPointF>& rawData = r.dataRaw ? *r.dataRaw : emptyRaw;

        double nativeHz = sr;
        const bool hasRealRaw = rawData.size() >= 2
            && rawData.last().x() > rawData.first().x()
            && !(rawData.size() == 1 && rawData[0].x() == -1.0);
        if (hasRealRaw)
            nativeHz = (rawData.size() - 1)
            / (rawData.last().x() - rawData.first().x());

        const double pxPerSec = (m_windowDuration > 0.0)
            ? r.chartView->chart()->plotArea().width() / m_windowDuration : 0.0;
        const double pxPerSample = (nativeHz > 0.0) ? pxPerSec / nativeHz : 0.0;
        r.chartView->setProperty("signalName", label);
        r.chartView->setProperty("nativeHz", nativeHz);
        r.chartView->chart()->setTitleFont(Theme::chartTitleFont());

        const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
        renderWindowedChart(
            r.chartView, { { r.upsampled_data, r.color, &rawData } },
            m_persistentLines[r.chartView],
            m_currentStartTime, m_windowDuration, globalOffset, sr,
            r.chartView == xLabelOwnerRight,
            m_plotMode == PlotMode::Scatter,
            m_plotMode == PlotMode::Line,
            yScaleForSignal(label), true);

        // --- pin axis to chunk-wide extremes when drift filtering is OFF ---
        if (!m_filterBaselineDrift) {
            double glo = 1e9, ghi = -1e9;
            for (double v : *r.upsampled_data) { if (v < glo) glo = v; if (v > ghi) ghi = v; }
            for (const QPointF& p : rawData) {
                if (p.x() == -1.0) continue;            // skip the (-1,-1) sentinel
                if (p.y() < glo) glo = p.y();
                if (p.y() > ghi) ghi = p.y();
            }
            auto vAxes = r.chartView->chart()->axes(Qt::Vertical);
            if (!vAxes.isEmpty()) {
                if (auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first()))
                    setPaddedYRange(yAxis, glo, ghi);   // reuse the same 5% padding helper
            }
        }
        double bpm = -1.0;
        if (m_showPeaks) {
            double dur = 0.0;
            const QVector<QPointF> bpmPeaks = get_bpm(label, dur);
            if (dur > 0.0) bpm = bpmPeaks.size() * 60.0 / dur;
        }
        r.chartView->setProperty("bpm", bpm);
        r.chartView->chart()->setTitle(formatChartTitle(label, nativeHz, pxPerSample, bpm));

        const QVector<QPointF> peaks = display_peaks_in_window(label);
        if (peaks.isEmpty()) return;

        const double yScale = yScaleForSignal(label);
        QList<QPointF> scaledPeaks;
        scaledPeaks.reserve(peaks.size());

        if (std::abs(yScale - 1.0) > 1e-9) {
            std::vector<double> vals;
            for (const QPointF& p : rawData) {
                if (p.x() < m_currentStartTime) continue;
                if (p.x() > m_currentStartTime + m_windowDuration) break;
                vals.push_back(p.y());
            }
            double center = 0.0;
            if (!vals.empty()) {
                const auto mid = vals.begin() + vals.size() / 2;
                std::nth_element(vals.begin(), mid, vals.end());
                center = *mid;
            }
            for (const QPointF& p : peaks)
                scaledPeaks.append({ p.x(), (p.y() - center) * yScale + center });
        }
        else {
            for (const QPointF& p : peaks) scaledPeaks.append(p);
        }

        auto* peakSeries = new QScatterSeries();
        peakSeries->setColor(Qt::red); peakSeries->setMarkerSize(8.0);
        peakSeries->setMarkerShape(QScatterSeries::MarkerShapeTriangle);
        peakSeries->setUseOpenGL(true); peakSeries->replace(scaledPeaks);
        r.chartView->chart()->addSeries(peakSeries);
        peakSeries->attachAxis(r.chartView->chart()->axes(Qt::Horizontal).first());
        peakSeries->attachAxis(r.chartView->chart()->axes(Qt::Vertical).first());
        };

    plotMarkable("ECG1"); plotMarkable("ECG2"); plotMarkable("ECG3");
    plotMarkable("PPG");
    if (ui->accel_or_abp_axis && !isMissingSignal(m_abp)) plotMarkable("ABP");

    auto plotDisplay = [&](QChartView* view, const QString& title,
        const QList<markable_data_series>& serieses)
        {
            if (!view || !view->chart()) return;
            view->chart()->setTitle(title);
            view->chart()->setTitleFont(Theme::chartTitleFont());
            const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;

            renderWindowedChart(view, serieses, m_persistentLines[view],
                m_currentStartTime, m_windowDuration, globalOffset, m_ecgSR,
                true, m_plotMode == PlotMode::Scatter, false,
                1.0, true);
        };

    if (ui->cvp_axis && !isMissingSignal(m_cvp)) {
        plotDisplay(ui->cvp_axis, "CVP", { { &m_cvp, COLOR_CVP, &m_cvpRaw } });
    }

    if (ui->hyp_accel_resp_axis) {
        const bool anyAccel = !isMissingSignal(m_accelX)
            || !isMissingSignal(m_accelY) || !isMissingSignal(m_accelZ);
        if (anyAccel) {
            plotDisplay(ui->hyp_accel_resp_axis, "ACCEL", {
                { &m_accelX, COLOR_ACCEL_X, &m_accelXRaw },
                { &m_accelY, COLOR_ACCEL_Y, &m_accelYRaw },
                { &m_accelZ, COLOR_ACCEL_Z, &m_accelZRaw },
                });
        }
        else if (!sleepPresent && !isMissingSignal(m_resp)) {
            plotDisplay(ui->hyp_accel_resp_axis, "RESP",
                { { &m_resp, COLOR_RESP, &m_respRaw } });
        }
    }

    updateNoiseHighlights();
    if (m_pulseOverlay) m_pulseOverlay->refresh();
    if (m_gapIndicator)  m_gapIndicator->refresh();
}

// ============================================================================
// Cursor & highlight updates
// ============================================================================

void noise_marking_gui::updateAmpogramCursor() {
    auto draw = [this](QChartView* view, QLineSeries* cursor) {
        if (!view || !cursor) return;
        auto axes = view->chart()->axes(Qt::Vertical);
        if (axes.isEmpty()) return;
        auto* yAxis = qobject_cast<QValueAxis*>(axes.first());
        if (!yAxis) return;
        double x = m_currentChunkIndex * seconds_in_memory_at_once
            + m_currentStartTime + m_windowDuration / 2.0;
        cursor->replace({ {x, yAxis->min()}, {x, yAxis->max()} });
        };
    draw(ui->ecg_ampogram_axis, m_ecgCursorBar);
    draw(ui->ppg_ampogram_axis, m_ppgCursorBar);
    if (sleepDataPresent(m_sleepStages))
        draw(ui->hyp_accel_resp_axis, m_hypnoCursorBar);
}

void noise_marking_gui::updateNoiseHighlights() {
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area;
    }
    m_highlights.clear();

    struct ChartAxes {
        QChart* chart = nullptr; QAbstractAxis* xAxis = nullptr; QValueAxis* yAxis = nullptr;
    };
    QMap<QString, ChartAxes> axesMap;
    for (const QString& lbl : markableChannelLabels()) {
        if (!isChannelActive(lbl)) continue;
        auto* cv = chartViewForSignalLabel(lbl);
        if (!cv) continue;
        ChartAxes ca;
        ca.chart = cv->chart();
        auto hAxes = ca.chart->axes(Qt::Horizontal);
        auto vAxes = ca.chart->axes(Qt::Vertical);
        ca.xAxis = hAxes.isEmpty() ? nullptr : hAxes.first();
        ca.yAxis = vAxes.isEmpty() ? nullptr
            : qobject_cast<QValueAxis*>(vAxes.first());
        axesMap[lbl] = ca;
    }

    const double globalOffset = m_currentChunkIndex * seconds_in_memory_at_once;
    const double viewStart = m_currentStartTime;
    const double viewEnd = viewStart + m_windowDuration;

    for (const auto& seg : m_noiseManager->getSegments()) {
        QString segLabel = QString::fromStdString(seg.label);
        if (!axesMap.contains(segLabel)) continue;
        const double sr = sampleRateForSignal(segLabel);
        const double segStart = seg.startSample / sr - globalOffset;
        const double segEnd = seg.endSample / sr - globalOffset;
        if (segEnd < viewStart || segStart > viewEnd) continue;

        const double ds = std::max(segStart, viewStart);
        const double de = std::min(segEnd, viewEnd);
        const QColor color = MARKING_COLORS.value(
            QString::fromStdString(seg.marking_type), QColor(0, 0, 0, 100));

        const ChartAxes& ca = axesMap[segLabel];
        if (!ca.chart || !ca.xAxis || !ca.yAxis) continue;

        auto* upper = new QLineSeries();
        auto* lower = new QLineSeries();
        upper->append({ {ds, ca.yAxis->max()}, {de, ca.yAxis->max()} });
        lower->append({ {ds, ca.yAxis->min()}, {de, ca.yAxis->min()} });

        auto* area = new QAreaSeries(upper, lower);
        area->setBrush(color); area->setPen(Qt::NoPen);
        ca.chart->addSeries(area);
        area->attachAxis(ca.xAxis); area->attachAxis(ca.yAxis);
        m_highlights.append(area);
    }
}