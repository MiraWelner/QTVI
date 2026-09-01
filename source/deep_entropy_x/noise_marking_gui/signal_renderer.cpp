/**
 * @file   signal_renderer.cpp
 * @brief  Signal chart rendering for the noise-marking GUI:
 *           - setupHypnogram          (sleep-stage overview)
 *           - ampogram                (amplitude-variability overview)
 *           - renderWindowedChart     (core windowed signal renderer)
 *           - detectPeaks / display_peaks_in_window / get_bpm
 *           - handle_data_plot        (main per-window redraw)
 *           - updateAmpogramCursor
 *           - updateNoiseHighlights
 */

#include "gui_handler.h"
#include "chart_utils.hpp"
#include "gui_peak_finder.hpp"
#include "logging/user_mark_log.hpp"
#include "theme/theme.h"
#include "annotation_types.hpp"


#include <QtCharts/QAreaSeries>
#include <QOpenGLWidget>
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
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <QtCharts/QLegendMarker>
#include <iostream>

namespace {

    int firstRawAtOrAfter(const QVector<QPointF>& raw, double target) {
        /*
            First index whose x >= target. rawData x is monotonic (index-time after the
            loader's rewrite), so we can binary-search the window start instead of
            linear-skipping every sample before it.
        */
        int lo = 0, hi = raw.size();
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (raw[mid].x() < target) lo = mid + 1; else hi = mid;
        }
        return lo;
    }

    // One eligible annotation reduced to what post-tagging needs: where it ends
    // (chunk-local s) and the tag its following beat should carry.
    struct PostSpan { double end; int tag; };

    // Tag each peak that is the FIRST detected beat after an eligible
    // annotation's end. `peaks` are this window's detected peaks, ascending in
    // x. A peak p[k] is a post beat for span s when s.end lies strictly between
    // the previous peak and p[k] (so no beat sits between the annotation end and
    // p[k]) and s.end is inside the visible window (so the gap is fully on
    // screen -- a beat whose annotation ended off the left edge isn't guessed
    // at). When several eligible annotations qualify, the nearest one wins.
    std::vector<int> tagPostBeats(const QVector<QPointF>& peaks,
        const std::vector<PostSpan>& spans, double detStart, double detEnd) {
        std::vector<int> tags(peaks.size(), 0);
        for (int k = 0; k < peaks.size(); ++k) {
            const double prevX = (k > 0) ? peaks[k - 1].x() : -1e300;
            const double px = peaks[k].x();
            double bestEnd = -1e300;
            int bestTag = 0;
            for (const PostSpan& s : spans) {
                if (s.end <= prevX || s.end >= px) continue;        // a beat lies between
                if (s.end < detStart || s.end > detEnd) continue;   // annotation end off screen
                if (s.end > bestEnd) { bestEnd = s.end; bestTag = s.tag; }
            }
            if (bestTag != 0) tags[k] = bestTag;
        }
        return tags;
    }

    // Per-channel annotation spans reduced to (start, end, markCode), in GLOBAL
    // seconds, sorted by start, with a prefix-max of ends. codeAt(gt) answers
    // "which annotation covers the beat at global time gt" in O(log n): binary-
    // search the last span starting at/<= gt, then walk back only while a span
    // could still reach gt (the prefix-max short-circuits the common no-mark
    // case). Built once per channel per redraw.
    //
    // Keeping this in one struct means the normal and accel logging paths
    // share the exact same span lookup, so they can't drift apart.
    struct MarkSpanIndex {
        struct Span { double s, e; int code; };
        std::vector<Span>   spans;
        std::vector<double> prefMaxEnd;

        int codeAt(double gt) const {
            int lo = 0, hi = static_cast<int>(spans.size());
            while (lo < hi) { const int mid = (lo + hi) >> 1; if (spans[mid].s <= gt) lo = mid + 1; else hi = mid; }
            if (lo == 0 || prefMaxEnd[lo - 1] < gt) return 0;   // nothing starting <= gt reaches gt
            for (int j = lo - 1; j >= 0; --j) {                 // only runs when the beat is inside a mark
                if (spans[j].e >= gt) return spans[j].code;
                if (prefMaxEnd[j] < gt) break;
            }
            return 0;
        }
    };

    MarkSpanIndex buildMarkSpanIndex(const annotation_handler* mgr,
        const std::string& label, double sr) {
        MarkSpanIndex idx;
        if (!mgr) return idx;
        for (const auto& seg : mgr->getSegments()) {
            if (seg.label != label) continue;
            idx.spans.push_back({ seg.startSample / sr, seg.endSample / sr,
                                  annotation_types::markCode(seg.marking_type) });
        }
        std::sort(idx.spans.begin(), idx.spans.end(),
            [](const MarkSpanIndex::Span& a, const MarkSpanIndex::Span& b) { return a.s < b.s; });
        idx.prefMaxEnd.resize(idx.spans.size());
        double running = -1e300;
        for (size_t i = 0; i < idx.spans.size(); ++i) {
            running = std::max(running, idx.spans[i].e);
            idx.prefMaxEnd[i] = running;
        }
        return idx;
    }

    QCategoryAxis* make_time_labled_xaxis(double startLocal, double duration,
        double globalOffset, bool labelsVisible)
    {
        auto* xAxis = new QCategoryAxis();
        xAxis->setRange(startLocal, startLocal + duration);
        xAxis->append(QString::fromStdString("(HH:MM:SS)"), startLocal);

        // Ticks at 1/4, 1/2, 3/4 of the window. get_timestamp shows tenths,
        // so labels stay distinct at every width (including the 1 s view).
        const double fracs[] = { 0.25, 0.5, 0.75 };
        for (double f : fracs) {
            const double t = startLocal + f * duration;
            xAxis->append(get_timestamp(globalOffset + t), t);
        }

        xAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
        xAxis->setTruncateLabels(false);
        xAxis->setGridLineVisible(false);
        xAxis->setLabelsVisible(labelsVisible);
        return xAxis;
    }

    std::pair<double, double> renderWindowedChart(
        QChartView* view,
        const QList<markable_data_series>& serieses,
        QList<QLineSeries*>& persistentLines,
        QList<QScatterSeries*>& persistentRawScatter,
        double currentStartTime, double windowDuration,
        double globalOffset, double ecgSR,
        bool labelsVisible,
        bool useScatterMode,
        bool forceLineForUpsampled,
        double yScale = 1.0)
    {
        if (!view || !view->chart()) return { 1e9, -1e9 };
        QChart* chart = view->chart();
        chart->legend()->hide();
        chart->setMargins(QMargins(0, 0, 0, 4));

        QSet<QAbstractSeries*> persistentSet;
        for (auto* ln : persistentLines) persistentSet.insert(ln);
        for (auto* sc : persistentRawScatter) if (sc) persistentSet.insert(sc);
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
            if (!d.data || is_missing_signal(*d.data)) continue;
            const bool hasRaw = d.rawData;

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
            int endIdx = std::clamp(
                static_cast<int>((currentStartTime + windowDuration) * ecgSR) + 1,
                0, static_cast<int>(d.data->size()));

            //for scaling, use median as the center, so the scaling is robust to noise
            std::vector<double> winVals;
            winVals.reserve(endIdx - startIdx);
            for (int i = startIdx; i < endIdx; ++i) {
                double v = (*d.data)[i];
                if (std::isnan(v)) continue;
                winVals.push_back(v);
            }
            const double winEnd = currentStartTime + windowDuration;
            for (int i = firstRawAtOrAfter(*d.rawData, currentStartTime);
                i < d.rawData->size(); ++i) {
                const QPointF& p = (*d.rawData)[i];
                if (p.x() > winEnd) break;
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
                if (std::isnan(raw)) continue;   // gap: don't feed NaN to the OpenGL line
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

        for (int ri = 0; ri < rawsToAdd.size(); ++ri) {
            const auto& r = rawsToAdd[ri];
            while (persistentRawScatter.size() <= ri) persistentRawScatter.append(nullptr);
            QScatterSeries*& rawScatter = persistentRawScatter[ri];
            if (!rawScatter) {
                rawScatter = new QScatterSeries();
                rawScatter->setBorderColor(Qt::transparent);
                rawScatter->setMarkerSize(3.0);
                rawScatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
                rawScatter->setUseOpenGL(true);
            }
            rawScatter->setColor(Qt::black);

            QList<QPointF> rawPts;
            const int firstIdx = firstRawAtOrAfter(*r.rawData, currentStartTime);
            const double winEnd = currentStartTime + windowDuration;

            // Count what's in the window, then stride so we draw at most kMaxDots
            // markers. Scatter markers are the dominant paint cost; past a few
            // thousand they can't be visually resolved and just stall the GPU.
            int lastIdx = firstIdx;
            while (lastIdx < r.rawData->size() && (*r.rawData)[lastIdx].x() <= winEnd)
                ++lastIdx;
            const int inWindow = lastIdx - firstIdx;
            constexpr int kMaxDots = 3000;
            const int stride = std::max(1, inWindow / kMaxDots);

            rawPts.reserve(std::min(inWindow, kMaxDots) + 1);
            for (int i = firstIdx; i < lastIdx; i += stride) {
                const QPointF& p = (*r.rawData)[i];
                rawPts.append({ p.x(), (p.y() - r.center) * yScale + r.center });
            }

            rawScatter->replace(rawPts);
            rawScatter->setVisible(true);

            if (rawScatter->chart()) chart->removeSeries(rawScatter);
            chart->addSeries(rawScatter);
            rawScatter->attachAxis(xAxis);
            rawScatter->attachAxis(yAxis);

            // Raw-sample dots are a data-authenticity overlay, not a distinct series
            // the reader should think of by name -- keep them out of the legend so
            // only the trace lines show up there.
            for (auto* marker : chart->legend()->markers(rawScatter))
                marker->setVisible(false);

        }
        // A frame with fewer overlays than a previous one: detach the leftovers
        // (kept alive for next time, just not shown).
        for (int ri = rawsToAdd.size(); ri < persistentRawScatter.size(); ++ri)
            if (persistentRawScatter[ri] && persistentRawScatter[ri]->chart())
                chart->removeSeries(persistentRawScatter[ri]);

        set_padded_y_range(yAxis, gMin, gMax);
        return { gMin, gMax };
    }

    // Walk back from `end`, counting only UNMARKED time, until `need` seconds of
    // it are collected (or we reach 0). Returns the resulting reference start.
    // `spans` are this channel's marked intervals (chunk-local s); they may
    // overlap and need not be sorted. This is how a mark gets skipped so the
    // reference always reaches a full 10 s of clean data.
    double reachBackUnmarked(double end, double need,
        std::vector<std::pair<double, double>> spans) {
        std::vector<std::pair<double, double>> m;
        {
            std::vector<std::pair<double, double>> c;
            for (auto& e : spans) {
                double a = std::max(0.0, e.first), b = std::min(end, e.second);
                if (a < b) c.push_back({ a, b });
            }
            std::sort(c.begin(), c.end());
            for (auto& s : c) {
                if (!m.empty() && s.first <= m.back().second)
                    m.back().second = std::max(m.back().second, s.second);
                else m.push_back(s);
            }
        }
        double pos = end;
        int j = static_cast<int>(m.size()) - 1;
        while (need > 0.0 && pos > 0.0) {
            while (j >= 0 && m[j].first >= pos) --j;
            if (j >= 0 && m[j].second >= pos) { pos = m[j].first; --j; continue; }
            const double lo = (j >= 0) ? m[j].second : 0.0;
            const double cleanLen = pos - lo;
            if (cleanLen >= need) { pos -= need; need = 0.0; }
            else { need -= cleanLen; pos = lo; }
        }
        return pos;
    }

} // namespace

// Reference window (chunk-local seconds) for peak-finder statistics: the 60 s
// preceding the visible window, or the 60 s following it when the visible
// window starts within the first 60 s of the chunk.
std::pair<double, double>
noise_marking_gui::statsWindow(double detStart, double detEnd) const {
    constexpr double kStatsWindowSec = gui_peak_finder::previous_seconds_to_train_on;
    if (detStart >= kStatsWindowSec)
        return { detStart - kStatsWindowSec, detStart };          // preceding 10 s
    double end = detEnd + kStatsWindowSec;                        // following 10 s
    const double chunkDur = totalChunkDuration();
    if (chunkDur > 0.0) end = std::min(end, chunkDur);
    return { detEnd, end };
}

// Core peak detection over an explicit [detStart, detEnd] window (chunk-local
// seconds). All annotation-aware behaviour lives here -- exclusions, per-region
// overrides, within-annotation reference, and the post-beat two-pass -- so any
// caller (the on-screen peaks and the BPM readout) detects identically and
// differs only in the window it asks for.
QVector<QPointF> noise_marking_gui::detectPeaks(const QString& label,
    double detStart, double detEnd, std::vector<int>* outPostTags) const {
    data_channel_features r = channelRefs(label);
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

    // Marked regions on THIS channel (chunk-local seconds). Two exclusion
    // lists are built here and passed to the peak finder as its do_not_learn_from_region
    // and no_peaks_in_region parameters:
    //   do_not_learn_from_region:  R peaks in this region are not used to set threshold
    //   no_peaks_in_region:        This has been marked with R peak noise and no R peaks detected here.
    //   withinSpans: non-noise annotations longer than the reference window;
    //                beats inside them take stats from WITHIN the annotation
    //                rather than reaching back to clean pre-annotation data.
    //   postSpans:   the five post-eligible annotations (AF/SVT/VT/PVC/PAC),
    //                reduced to (end, tag) for marking the beat that follows.
    constexpr double kRefSec = gui_peak_finder::previous_seconds_to_train_on;
    const double sr = sampleRateForSignal(label);
    std::vector<std::pair<double, double>> do_not_learn_from_region, no_peaks_in_region, withinSpans;
    std::vector<PostSpan> postSpans;
    const std::string labelStd = label.toStdString();   // convert once, not per segment
    if (m_noiseManager) {
        for (const auto& seg : m_noiseManager->getSegments()) {
            if (seg.label != labelStd) continue;        // std::string compare, no per-seg QString alloc
            const double s = seg.startSample / sr - globalOffset;
            const double e = seg.endSample / sr - globalOffset;
            // Most markings are excluded from the reference-window stats
            // (threshold gate + mean R-R). Types flagged includeInThreshold
            // (e.g. Minor Noise) are kept in those stats, so they have no
            // effect on detection or on where the following R peaks land.
            const bool inclThr = annotation_types::includeInThreshold(seg.marking_type);
            if (!inclThr)
                do_not_learn_from_region.push_back({ s, e });
            if (annotation_types::suppressesDetection(seg.marking_type))
            {
                no_peaks_in_region.push_back({ s, e });
            }
            else if (!inclThr && e - s > kRefSec)
            {
                withinSpans.push_back({ s, e });
            }
            const int tag = annotation_types::postCode(seg.marking_type);
            if (tag != 0) postSpans.push_back({ e, tag });
        }
    }

    // Reference window = the previous 10 seconds of UNMARKED data: start at the
    // window's left edge and walk back, skipping any annotated spans, until a
    // full 10 s of clean data has been gathered. Near the chunk start (no 10 s
    // behind us) fall back to statsWindow's following window.
    // (refStart/refEnd are now advisory only -- the finders compute a per-beat
    // reference internally -- but are still passed for signature compatibility.)
    double refStart, refEnd; bool refPreceding;
    {
        if (detStart >= kRefSec) {
            refEnd = detStart; refPreceding = true;
            refStart = reachBackUnmarked(detStart, kRefSec, do_not_learn_from_region);
        }
        else {
            const auto [s, e] = statsWindow(detStart, detEnd);
            refStart = s; refEnd = e; refPreceding = (s < detStart);
        }
    }

    // Per-peak parameter accessors, backed by an O(log n) piecewise-constant
    // index built once from THIS channel's overrides (in chunk-local seconds).
    // A window with no nearby overrides pays nothing for overrides elsewhere
    // in the chunk.
    gui_peak_finder::ParamIndex thrIdx, blkIdx;
    {
        std::vector<std::tuple<double, double, double>> tSeg, bSeg;
        for (const ParamOverride& o : m_thresholdOverrides)
            if (o.channel == label)
                tSeg.push_back({ o.start - globalOffset, o.end - globalOffset, o.value });
        for (const ParamOverride& o : m_blankingOverrides)
            if (o.channel == label)
                bSeg.push_back({ o.start - globalOffset, o.end - globalOffset, o.value });
        thrIdx = gui_peak_finder::ParamIndex::build(m_cfg.threshold, std::move(tSeg));
        blkIdx = gui_peak_finder::ParamIndex::build(m_cfg.blanking_period, std::move(bSeg));
    }
    auto thrFn = [thrIdx](double localT) { return thrIdx.at(localT); };
    auto blkFn = [blkIdx](double localT) { return blkIdx.at(localT); };

    const bool baseInverted = invertedForSignal(label);
    std::vector<std::pair<double, double>> invSpans;                // chunk-local
    for (const ParamOverride& o : m_invertOverrides)
        if (o.channel == label)
            invSpans.push_back({ o.start - globalOffset, o.end - globalOffset });
    auto sgnFn = [baseInverted, invSpans](double localT) -> double {
        bool inv = baseInverted;
        for (const auto& s : invSpans)
            if (localT >= s.first && localT <= s.second) { inv = !inv; break; }
        return inv ? -1.0 : 1.0;
        };
    auto runFinder = [&](const std::vector<std::pair<double, double>>& refEx) {
        // Pressure/pulse waveforms (PPG, ABP, ART, ART_PULM) are systolic
        // upstroke + rounded apex -> derivative detector. Only the three ECG
        // leads use the local-max detector (and thus the per-sample sign fn).
        if (label == "PPG" || label == "ABP" || label == "ART" || label == "ART_PULM")
            return gui_peak_finder::findPeaksDerivative(
                *r.dataRaw, detStart, detEnd, refStart, refEnd, thrFn, blkFn,
                refPreceding, refEx, no_peaks_in_region, withinSpans);
        return gui_peak_finder::findPeaks(
            *r.dataRaw, detStart, detEnd, refStart, refEnd, thrFn, blkFn,
            refPreceding, sgnFn, refEx, no_peaks_in_region, withinSpans);   // sgnFn, not scalar
        };


    // Pass 1: detect with annotation spans excluded from the reference.
    QVector<QPointF> peaks = runFinder(do_not_learn_from_region);

    // A post beat (first beat after AF/SVT/VT/PVC/PAC) must not contribute to
    // any other beat's reference. Tag the post beats ONCE here; use the tags
    // both to exclude those beats (a hair-wide window around each, then a
    // re-detect so their compensatory pause / altered amplitude can't skew
    // neighbouring gates) and to hand back to the caller for colouring and
    // logging. The re-run only fires when a post beat is in view; the common
    // no-arrhythmia case pays nothing.
    std::vector<int> tags;
    if (!postSpans.empty()) {
        tags = tagPostBeats(peaks, postSpans, detStart, detEnd);
        constexpr double kPostEps = 0.06;   // ~one QRS/systole half-width
        std::vector<std::pair<double, double>> refEx2 = do_not_learn_from_region;
        bool any = false;
        for (int k = 0; k < peaks.size(); ++k)
            if (tags[k] != 0) {
                refEx2.push_back({ peaks[k].x() - kPostEps, peaks[k].x() + kPostEps });
                any = true;
            }
        if (any) {
            QVector<QPointF> peaks2 = runFinder(refEx2);
            // Re-tag only if excluding the post beats actually changed the peak
            // set (rare); otherwise the pass-1 tags still line up one-to-one.
            if (peaks2 != peaks) tags = tagPostBeats(peaks2, postSpans, detStart, detEnd);
            peaks = std::move(peaks2);
        }
    }

    if (outPostTags) {
        if (tags.size() != static_cast<size_t>(peaks.size()))
            tags.assign(peaks.size(), 0);   // no eligible annotations in view
        *outPostTags = std::move(tags);
    }
    return peaks;
}

// Padding (seconds) added on each side of the requested window before
// detection, then cropped away. The reference/gate is already window-
// independent (it reads fixed chunk-anchored frames), so the only window-
// sensitive parts are short-range -- the local-max neighbour, the blanking
// run-up, and the post-beat tag distance, all under ~2 s. Detecting this far
// past each edge makes every *visible* beat interior to the detection range,
// so its detection and red/blue classification are independent of where the
// window edge falls. 5 s is comfortably above all those reaches.
static constexpr double kDetectMargin = 5.0;

// On-screen peaks: detect over a padded interval, then crop to the visible
// window. Cropping (not a narrower detection) is what makes the result
// reproducible regardless of scroll position.
QVector<QPointF> noise_marking_gui::display_peaks_in_window(const QString& label, std::vector<int>* outPostTags) const {
    const double winStart = current_start_time;
    const double winEnd = current_start_time + visible_window_size;
    const double chunkDur = totalChunkDuration();
    const double detStart = std::max(0.0, winStart - kDetectMargin);
    const double detEnd = (chunkDur > 0.0) ? std::min(chunkDur, winEnd + kDetectMargin)
        : winEnd + kDetectMargin;

    std::vector<int> allTags;
    const QVector<QPointF> all = detectPeaks(label, detStart, detEnd,
        outPostTags ? &allTags : nullptr);

    QVector<QPointF> out;
    out.reserve(all.size());
    std::vector<int> tags;
    if (outPostTags) tags.reserve(all.size());
    for (int i = 0; i < all.size(); ++i) {
        if (all[i].x() < winStart || all[i].x() > winEnd) continue;   // drop the padding
        out.append(all[i]);
        if (outPostTags) tags.push_back(i < static_cast<int>(allTags.size()) ? allTags[i] : 0);
    }
    if (outPostTags) *outPostTags = std::move(tags);
    return out;
}

QVector<QPointF> noise_marking_gui::get_bpm(const QString& label, double& outDuration) const
{
    /*
        Beats-per-minute readout. Uses the SAME detection as the on-screen peaks
        (detectPeaks -- same exclusions, per-region overrides, within-annotation
        reference, and post-beat handling), so the rate always agrees with the
        beats that are drawn. The only difference is the window: BPM never uses a
        span shorter than kMinBpmWindowSec, so a rate is still stable when the
        visible window is only a second or two wide.
    */
    constexpr double kMinBpmWindowSec = 10.0;
    const double tVisEnd = current_start_time + visible_window_size;
    double tStart = current_start_time;
    if (visible_window_size < kMinBpmWindowSec)
        tStart = std::max(0.0, tVisEnd - kMinBpmWindowSec);
    outDuration = tVisEnd - tStart;

    // Same padded-detect-then-crop as display_peaks_in_window, so the rate
    // agrees with the drawn beats and doesn't shift with scroll position.
    const double chunkDur = totalChunkDuration();
    const double detStart = std::max(0.0, tStart - kDetectMargin);
    const double detEnd = (chunkDur > 0.0) ? std::min(chunkDur, tVisEnd + kDetectMargin)
        : tVisEnd + kDetectMargin;

    const QVector<QPointF> all = detectPeaks(label, detStart, detEnd);
    QVector<QPointF> out;
    out.reserve(all.size());
    for (const QPointF& p : all)
        if (p.x() >= tStart && p.x() <= tVisEnd) out.append(p);
    return out;
}

void noise_marking_gui::setupHypnogram() {
    if (m_sleepSR <= 0.0 || !sleep_data_present(m_sleepStages)) return;

    auto* chart = ui->hyp_resp_axis->chart();
    if (m_cvpCursorBar && m_cvpCursorBar->chart() == chart)
        chart->removeSeries(m_cvpCursorBar);
    for (auto* s : m_hypnoStageSeries) { chart->removeSeries(s); delete s; }
    m_hypnoStageSeries.clear();

    struct Stage { int value; QColor color; const char* name; };
    const QList<Stage> stages = {
        {0, Qt::black,     "Wake" },
        {1, Qt::darkGreen, "NREM1"},
        {2, Qt::blue,      "NREM2"},
        {3, Qt::cyan,      "NREM3"},
        {4, Qt::red,       "REM"  }
    };

    const double dt = 1.0 / m_sleepSR;
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

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

    // The five rows above are the AASM collapse (N3 and N4 merged, REM last),
    // which is what MESA bins already carry. Compumedics staging -- what the
    // SHHS annotation XML actually contains -- numbers the stages 0=Wake,
    // 1=N1, 2=N2, 3=N3, 4=N4, 5=REM, 9=unscored. A file whose staging was NOT
    // collapsed on the way into the .bin therefore loses every REM epoch here
    // and draws N4 on the row LABELLED REM, which looks plausible rather than
    // broken. The collapse belongs in file_to_bin next to the XML parse, so
    // this reports the mismatch instead of guessing a remap at draw time.
    // (Codes below 0 are placeholder padding and are not counted.)
    {
        int unplotted = 0;
        int firstCode = 0;
        for (double v : m_sleepStages) {
            const int code = static_cast<int>(v);
            if (code > 4) {
                if (unplotted++ == 0) firstCode = code;
            }
        }
        if (unplotted > 0)
            std::cerr << "[hypnogram] " << m_cfg.dataset_type << ": " << unplotted
            << " of " << m_sleepStages.size() << " epochs in this chunk carry a "
            "stage code above 4 (first: " << firstCode << ") and are not drawn -- "
            "staging may not be collapsed to Wake/N1/N2/N3/REM\n";
    }

    for (QAbstractAxis* axis : chart->axes()) {
        chart->removeAxis(axis);
    }

    chart->setMargins(QMargins(0, 0, 12, 5));
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
    yAxis->setLabelsPosition(QCategoryAxis::AxisLabelsPositionOnValue);
    for (const auto& st : stages)
        yAxis->append(st.name, st.value);     // label drawn exactly at the row value
    yAxis->setRange(-0.5, 4.5);
    yAxis->setStartValue(-0.5);
    yAxis->setReverse(true);
    yAxis->setTruncateLabels(false);
    yAxis->setVisible(true);
    yAxis->setGridLineVisible(false);
    yAxis->setLabelsFont(Theme::chartAxisFont());

    chart->addAxis(xAxis, Qt::AlignBottom);
    chart->addAxis(yAxis, Qt::AlignLeft);
    for (auto* s : chart->series()) { s->attachAxis(xAxis); s->attachAxis(yAxis); }
    const QColor axisGray = yAxis->linePenColor();
    xAxis->setLabelsColor(axisGray);
    yAxis->setLabelsColor(axisGray);
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
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;

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
            x_axis->append(QString::fromStdString("(h)"), globalOffset + 1800);
            const int startHour = static_cast<int>(globalOffset / 3600);
            for (int h = 1; h <= seconds_in_memory_at_once / 3600; ++h) {
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

    auto ecg1Pts = calculate_amplitude(m_ecg1, channel_upsampled_rates[CH_ECG1]);
    auto ecg2Pts = calculate_amplitude(m_ecg2, channel_upsampled_rates[CH_ECG1]);
    auto ecg3Pts = calculate_amplitude(m_ecg3, channel_upsampled_rates[CH_ECG1]);

    const bool sleepPresent = sleep_data_present(m_sleepStages);
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
        set_padded_y_range(yAxis, mi->y(), ma->y());
    }

    create_plot(ui->ppg_ampogram_axis, ppg_ampogram_series,
        calculate_amplitude(m_ppg, channel_upsampled_rates[CH_PPG]),
        m_ppgCursorBar, COLOR_PPG, "PPG Amp-O-Gram", ppgAmpHasLabels);
}

void noise_marking_gui::plot_nonmarkable(QChartView* view, const QString& title, const QList<markable_data_series>& serieses, double sampling_rate)
{
    if (!view || !view->chart()) return;
    view->chart()->setTitle(title);
    view->chart()->setTitleFont(Theme::chartTitleFont());
    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    renderWindowedChart(view, serieses, m_persistentLines[view],
        m_persistentRawScatter[view],
        current_start_time, visible_window_size, globalOffset, sampling_rate,
        true, m_plotMode == PlotMode::Scatter,
        m_plotMode == PlotMode::Line, 1.0);
}

void noise_marking_gui::determine_which_nonmarkable_charts_to_plot() {
    /**
    * @brief Fill the three shared non-markable slots (cvp_eeg_axis, hyp_resp_axis,
    *        pacemaker_axis), choosing what goes in each by dataset:
    *        BITTIUM -> temperature / marker / pacemaker events,
    *        CHAOS   -> CVP / RESP,
    *        SHHS1|SHHS2 -> AIRFLOW+THOR+ABDO / (hypnogram) / -- ,
    *        MESA    -> (hypnogram) only.
    *        SHHS's SaO2 is NOT one of these three: it has its own sao2_axis at
    *        the bottom of the main plot column and is plotted here too, at the
    *        end of the SHHS branch.
    *        Visibility for the same slots is set in loadChunkFromFile; the two
    *        must agree or a chart is shown empty (or hidden with data in it).
    */
    if (m_cfg.dataset_type == "BITTIUM") {
        if (ui->cvp_eeg_axis && !is_missing_signal(m_temp))
            plot_nonmarkable(ui->cvp_eeg_axis, "Temperature",
                { { &m_temp, COLOR_TEMP, &m_tempRaw } }, channel_upsampled_rates[CH_TEMP]);
        if (ui->hyp_resp_axis && !is_missing_signal(m_marker))
            plot_nonmarkable(ui->hyp_resp_axis, "Marker",
                { { &m_marker, COLOR_MARKER, &m_markerRaw } }, channel_upsampled_rates[CH_MARKER]);
        if (ui->pacemaker_axis && !is_missing_signal(m_pacemaker))
            plot_nonmarkable(ui->pacemaker_axis, "Pacemaker Events",
                { { &m_pacemaker, COLOR_MARKER, &m_pacemakerRaw } }, channel_upsampled_rates[CH_PACEMAKER_EVENT]);
    }
    if (m_cfg.dataset_type == "CHAOS") {
        if (ui->cvp_eeg_axis && m_cvpRaw.size() >= 2)
            plot_nonmarkable(ui->cvp_eeg_axis, "CVP",
                { { &m_cvp, COLOR_CVP, &m_cvpRaw } }, channel_upsampled_rates[CH_CVP]);
        if (ui->hyp_resp_axis && !is_missing_signal(m_resp))
            plot_nonmarkable(ui->hyp_resp_axis, "RESP",
                { { &m_resp, COLOR_RESP, &m_respRaw } }, channel_upsampled_rates[CH_RESP]);
    }
    if (m_cfg.dataset_type == "SHHS1" || m_cfg.dataset_type == "SHHS2") {
        // AIRFLOW, THOR RES and ABDO RES share cvp_eeg_axis. They are all in
        // arbitrary units, so one y axis is not a unit clash -- the shape is
        // what a reviewer reads off this chart -- but they must share a TIME
        // base, because plot_nonmarkable takes a single rate for the whole
        // chart. Config gives the three the same upsample rate in every SHHS
        // row so far; a channel that disagrees is left out rather than drawn
        // at the wrong times.
        const float respRate = (channel_upsampled_rates[CH_THOR] > 0.0f)
            ? channel_upsampled_rates[CH_THOR]
            : channel_upsampled_rates[CH_FLOW];
        QList<markable_data_series> resp;
        auto addResp = [&](const QVector<double>& sig, const QVector<QPointF>& raw,
            int ch, const QColor& color) {
                if (is_missing_signal(sig)) return;
                if (channel_upsampled_rates[ch] != respRate) {
                    std::cerr << "[shhs] channel " << ch << " upsampled to "
                        << channel_upsampled_rates[ch] << " Hz but the respiratory "
                        "chart is on " << respRate << " Hz; not plotted\n";
                    return;
                }
                resp.append({ &sig, color, &raw });
            };
        addResp(m_flow, m_flowRaw, CH_FLOW, COLOR_FLOW);
        addResp(m_thor, m_thorRaw, CH_THOR, COLOR_THOR);
        addResp(m_abdo, m_abdoRaw, CH_ABDO, COLOR_ABDO);
        if (ui->cvp_eeg_axis && !resp.isEmpty() && respRate > 0.0f)
            plot_nonmarkable(ui->cvp_eeg_axis, "AIRFLOW / THOR / ABDO", resp, respRate);

        if (ui->sao2_axis && !is_missing_signal(m_spo2))
            plot_nonmarkable(ui->sao2_axis, "SaO2",
                { { &m_spo2, COLOR_SPO2, &m_spo2Raw } }, channel_upsampled_rates[CH_SPO2]);
    }
}

void noise_marking_gui::handle_data_plot() {
    clearDragPreview();
    for (auto* area : m_highlights) {
        if (area->chart()) area->chart()->removeSeries(area);
        delete area;
    }
    m_highlights.clear();
    for (const QString& lbl : markableChannelLabels()) {
        markStateFor(lbl).startMarkerLine = nullptr;
    }

    // The timestamp x-axis belongs to the bottom-most ACTIVE markable chart.
    // Walk markableChannelLabels() (top-to-bottom order) and keep the last
    // active one, so the labeled axis rides down to ART / ART_PULM when they
    // are present and falls back to PPG/ABP otherwise.
    QChartView* xLabelOwnerRight = nullptr;
    for (const QString& lbl : markableChannelLabels()) {
        if (!isChannelActive(lbl)) continue;
        if (QChartView* cv = chartViewForSignalLabel(lbl))
            xLabelOwnerRight = cv;   // last active wins
    }

    auto keepFor = [&](QChartView* cv) -> QList<QAbstractSeries*> {
        QList<QAbstractSeries*> keep;
        if (!cv) return keep;
        auto it = m_persistentLines.constFind(cv);
        if (it != m_persistentLines.constEnd())
            for (auto* ln : it.value()) if (ln) keep.append(ln);
        auto sit = m_persistentRawScatter.constFind(cv);
        if (sit != m_persistentRawScatter.constEnd())
            for (auto* sc : sit.value()) if (sc) keep.append(sc);
        if (m_pulseOverlay) {
            const QString label = signalLabelForChartView(cv);
            if (!label.isEmpty())
                for (QLineSeries* s : m_pulseOverlay->seriesForLabel(label))
                    if (s) keep.append(s);
        }
        return keep;
        };

    wipe_chart(ui->ecg_axis_1->chart(), keepFor(ui->ecg_axis_1));
    wipe_chart(ui->ecg_axis_2->chart(), keepFor(ui->ecg_axis_2));
    wipe_chart(ui->ecg_axis_3->chart(), keepFor(ui->ecg_axis_3));
    wipe_chart(ui->ppg_axis->chart(), keepFor(ui->ppg_axis));
    if (ui->accel_axis && !is_missing_signal(m_accelX))
        wipe_chart(ui->accel_axis->chart(), keepFor(ui->accel_axis));
    if (ui->abp_axis && !is_missing_signal(m_abp))
        wipe_chart(ui->abp_axis->chart(),
            keepFor(ui->abp_axis));
    if (ui->art_axis && !is_missing_signal(m_art))
        wipe_chart(ui->art_axis->chart(), keepFor(ui->art_axis));
    if (ui->art_pulm_axis && !is_missing_signal(m_artPulm))
        wipe_chart(ui->art_pulm_axis->chart(), keepFor(ui->art_pulm_axis));
    if (ui->kors_matrix && !is_missing_signal(m_vcg))
        wipe_chart(ui->kors_matrix->chart(), keepFor(ui->kors_matrix));

    const bool sleepPresent = sleep_data_present(m_sleepStages);
    if (!sleepPresent && ui->hyp_resp_axis)
        wipe_chart(ui->hyp_resp_axis->chart(),
            keepFor(ui->hyp_resp_axis));
    if (ui->cvp_eeg_axis)
        wipe_chart(ui->cvp_eeg_axis->chart(), keepFor(ui->cvp_eeg_axis));
    if (ui->pacemaker_axis)
        wipe_chart(ui->pacemaker_axis->chart(), keepFor(ui->pacemaker_axis));
    // Every chart redrawn per window has to be wiped here or its series pile up
    // across redraws and the trace stops tracking the scroll position.
    if (ui->sao2_axis)
        wipe_chart(ui->sao2_axis->chart(), keepFor(ui->sao2_axis));

    auto plotMarkable = [&](const QString& label) {
        if (!isChannelActive(label)) return;
        data_channel_features r = channelRefs(label);
        if (!r.chartView || !r.upsampled_data || !r.state) return;

        const QVector<QPointF> emptyRaw;
        const QVector<QPointF>& rawData = r.dataRaw ? *r.dataRaw : emptyRaw;
        QList<markable_data_series> serieses;
        if (label == "ACCEL") {
            const QVector<QPointF>& yRaw = m_accelYRaw.isEmpty() ? emptyRaw : m_accelYRaw;
            const QVector<QPointF>& zRaw = m_accelZRaw.isEmpty() ? emptyRaw : m_accelZRaw;
            if (!is_missing_signal(m_accelX))
                serieses.append({ &m_accelX, COLOR_ACCEL_X, &rawData });   // X's raw via channelRefs
            if (!is_missing_signal(m_accelY))
                serieses.append({ &m_accelY, COLOR_ACCEL_Y, &yRaw });
            if (!is_missing_signal(m_accelZ))
                serieses.append({ &m_accelZ, COLOR_ACCEL_Z, &zRaw });
            if (serieses.isEmpty()) return;
        }
        else {
            serieses.append({ r.upsampled_data, r.color, &rawData });
        }

        // "Original Frequency" comes straight from config (channel_native_rates,
        // carried on r.nativeRate) -- not inferred from the raw span, which is
        // wrong for gappy/mostly-absent channels (e.g. ART).
        double nativeHz = (r.nativeRate > 0.0) ? r.nativeRate : r.sampleRate;
        const double upHz = r.sampleRate;

        const double pxPerSec = (visible_window_size > 0.0)
            ? r.chartView->chart()->plotArea().width() / visible_window_size : 0.0;
        const double pxPerSample = pxPerSec / nativeHz;
        r.chartView->setProperty("signalName", label);
        r.chartView->setProperty("nativeHz", nativeHz);
        r.chartView->setProperty("upHz", r.sampleRate);

        const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
        renderWindowedChart(
            r.chartView, serieses,
            m_persistentLines[r.chartView],
            m_persistentRawScatter[r.chartView],
            current_start_time, visible_window_size, globalOffset, r.sampleRate,
            r.chartView == xLabelOwnerRight,
            m_plotMode == PlotMode::Scatter,
            m_plotMode == PlotMode::Line,
            yScaleForSignal(label));


        if (label == "PPG") {
            const QVector<double>& v = *r.upsampled_data;
            double lo = 1e18, hi = -1e18; int nan = 0, n = 0;
            for (double x : v) { if (std::isnan(x)) { ++nan; continue; } lo = std::min(lo, x); hi = std::max(hi, x); ++n; }
            auto vAxes = r.chartView->chart()->axes(Qt::Vertical);
            double amin = 0, amax = 0;
            if (!vAxes.isEmpty()) if (auto* ya = qobject_cast<QValueAxis*>(vAxes.first())) { amin = ya->min(); amax = ya->max(); }
        }

        // --- y-range: frozen (Fix Scale) wins; else chunk-wide pin when
        //     drift filtering is off; else the per-window autoscale that
        //     renderWindowedChart already applied. ---
        auto fixedIt = m_fixedYRange.constFind(label);
        if (fixedIt != m_fixedYRange.constEnd()) {
            auto vAxes = r.chartView->chart()->axes(Qt::Vertical);
            if (!vAxes.isEmpty())
                if (auto* yAxis = qobject_cast<QValueAxis*>(vAxes.first()))
                    yAxis->setRange(fixedIt->first, fixedIt->second);   // no re-pad; use captured range as-is
        }
        

        // ACCEL has no R-peaks: render the trace (done above), then log the accel
        // value at each ECG1 beat time into its own beat_log column, markType
        // taken from ACCEL's own annotations. No detection, no red/blue markers,
        // no bpm.
        if (label == "ACCEL") {

            r.chartView->setProperty("bpm", QVariant());
            r.chartView->chart()->setTitle(get_chart_title(label, nativeHz, pxPerSample, -1.0, r.sampleRate));

            {
                const char* names[] = { "X", "Y", "Z" };
                int nameIdx = 0;
                for (auto* s : r.chartView->chart()->series()) {
                    if (auto* ln = qobject_cast<QLineSeries*>(s)) {
                        if (nameIdx < 3) ln->setName(names[nameIdx++]);
                    }
                }
                r.chartView->chart()->legend()->setVisible(true);
                r.chartView->chart()->legend()->setAlignment(Qt::AlignRight);
            }

            if (m_beatLog) {
                const MarkSpanIndex markIdx =
                    buildMarkSpanIndex(m_noiseManager.get(), label.toStdString(), r.sampleRate);
                // Anchor rows to ECG1 beat times (accel produces no beats of its own).
                const QVector<QPointF> ecgPeaks = display_peaks_in_window("ECG1");
                for (const QPointF& p : ecgPeaks) {
                    const double gt = p.x() + globalOffset;
                    double accelY = 0.0;
                    const int idx = static_cast<int>(std::round(p.x() * r.sampleRate));
                    if (idx >= 0 && idx < m_accelX.size()) accelY = m_accelX[idx];
                    m_beatLog->logPeak(beat_log::ACCEL, gt, accelY,
                        blankingAt(label, gt), thresholdAt(label, gt),
                        markIdx.codeAt(gt), 0, 0);
                }
            }
            return;   // accel chart drawn above; no detection/markers
        }

        std::vector<int> postTags;
        const QVector<QPointF> peaks = display_peaks_in_window(label, &postTags);
        double bpm = 0.0;
        if (visible_window_size >= 10.0) {            // keep in sync with kMinBpmWindowSec in get_bpm
            if (visible_window_size > 0.0) bpm = peaks.size() * 60.0 / visible_window_size;
        }
        else {
            double dur = 0.0;
            const QVector<QPointF> bpmPeaks = get_bpm(label, dur);
            if (dur > 0.0) bpm = bpmPeaks.size() * 60.0 / dur;
        }
        r.chartView->setProperty("bpm", bpm);
        r.chartView->chart()->setTitle(get_chart_title(label, nativeHz, pxPerSample, bpm, r.sampleRate));

        // Log every detected beat: value, the blanking/threshold in effect,
        // which annotation (if any) covers it, its post-arrhythmia tag, and
        // whether the channel is inverted there.
        if (m_beatLog) {
            const beat_log::ChannelIdx ch = beat_log::channelForLabel(label);
            const MarkSpanIndex markIdx =
                buildMarkSpanIndex(m_noiseManager.get(), label.toStdString(), r.sampleRate);
            for (int k = 0; k < peaks.size(); ++k) {
                const double gt = peaks[k].x() + globalOffset;
                m_beatLog->logPeak(ch, gt, peaks[k].y(),
                    blankingAt(label, gt), thresholdAt(label, gt),
                    markIdx.codeAt(gt), postTags[k],
                    invertedAt(label, gt) ? 1 : 0);
            }
        }

        const double yScale = yScaleForSignal(label);
        QList<QPointF> scaledPeaks;
        scaledPeaks.reserve(peaks.size());

        if (std::abs(yScale - 1.0) > 1e-9) {
            std::vector<double> vals;
            for (const QPointF& p : rawData) {
                if (p.x() < current_start_time) continue;
                if (p.x() > current_start_time + visible_window_size) break;
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

        // Normal beats render red; post-arrhythmia beats render blue. Two
        // series so the colors don't bleed (a QScatterSeries is one color).
        QList<QPointF> redPts, bluePts;
        for (int k = 0; k < scaledPeaks.size(); ++k)
            ((k < (int)postTags.size() && postTags[k] != 0) ? bluePts : redPts)
            .append(scaledPeaks[k]);

        auto addPeakSeries = [&](const QList<QPointF>& pts, const QColor& color) {
            if (pts.isEmpty()) return;
            auto* s = new QScatterSeries();
            s->setColor(color); s->setMarkerSize(8.0);
            s->setMarkerShape(QScatterSeries::MarkerShapeTriangle);
            s->setUseOpenGL(true); s->replace(pts);
            r.chartView->chart()->addSeries(s);
            s->attachAxis(r.chartView->chart()->axes(Qt::Horizontal).first());
            s->attachAxis(r.chartView->chart()->axes(Qt::Vertical).first());
            };
        addPeakSeries(redPts, Qt::red);
        addPeakSeries(bluePts, Qt::blue);
        };
        plotMarkable("ECG1"); plotMarkable("ECG2"); plotMarkable("ECG3"); plotMarkable("VCG"); plotMarkable("PPG");  plotMarkable("ACCEL");
    if (ui->abp_axis && !is_missing_signal(m_abp)) plotMarkable("ABP");
    if (ui->art_axis && !is_missing_signal(m_art)) plotMarkable("ART");
    if (ui->art_pulm_axis && !is_missing_signal(m_artPulm)) plotMarkable("ART_PULM");

    determine_which_nonmarkable_charts_to_plot();
    updateNoiseHighlights();
    if (m_pulseOverlay) m_pulseOverlay->refresh();
    syncChunkScrollBar();

    // GL-accelerated series sit in a QOpenGLWidget overlaying the plot area,
    // which otherwise swallows clicks before they reach the viewport's event
    // filter. Let presses fall through so marking/erase land on the first click.
    for (QChartView* cv : { ui->ecg_axis_1, ui->ecg_axis_2, ui->ecg_axis_3,
                            ui->ppg_axis, ui->accel_axis, ui->abp_axis }) {
        if (!cv) continue;
        if (auto* gl = cv->findChild<QOpenGLWidget*>())
            gl->setAttribute(Qt::WA_TransparentForMouseEvents, true);
    }
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
        double x = current_chunk_index * seconds_in_memory_at_once
            + current_start_time + visible_window_size / 2.0;
        cursor->replace({ {x, yAxis->min()}, {x, yAxis->max()} });
        };
    draw(ui->ecg_ampogram_axis, m_ecgCursorBar);
    draw(ui->ppg_ampogram_axis, m_ppgCursorBar);
    if (sleep_data_present(m_sleepStages))
        draw(ui->hyp_resp_axis, m_hypnoCursorBar);
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

    const double globalOffset = current_chunk_index * seconds_in_memory_at_once;
    const double viewStart = current_start_time;
    const double viewEnd = viewStart + visible_window_size;

    // Per-channel sample rate keyed by std::string so off-screen segments
    // (the vast majority at high marking counts) are filtered with no
    // per-segment QString allocation; the QString is built only for segments
    // that actually fall in the visible window.
    std::unordered_map<std::string, double> srByLabel;
    for (const QString& lbl : markableChannelLabels())
        if (isChannelActive(lbl))
            srByLabel.emplace(lbl.toStdString(), sampleRateForSignal(lbl));

    for (const auto& seg : m_noiseManager->getSegments()) {
        auto srIt = srByLabel.find(seg.label);
        if (srIt == srByLabel.end()) continue;          // not an active channel; no alloc
        const double sr = srIt->second;
        const double segStart = seg.startSample / sr - globalOffset;
        const double segEnd = seg.endSample / sr - globalOffset;
        if (segEnd < viewStart || segStart > viewEnd) continue;   // off-screen; no alloc
        const QString segLabel = QString::fromStdString(seg.label);   // visible segments only

        const double ds = std::max(segStart, viewStart);
        const double de = std::min(segEnd, viewEnd);
        const QColor color = annotation_types::colorFor(QString::fromStdString(seg.marking_type));

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

    // Translucent rectangle behind any region where a per-channel threshold or
    // blanking override is in effect, so the user gets confirmation the edit
    // registered. Gray = threshold override, blue = blanking override.
    auto drawOverride = [&](const ParamOverride& o, const QColor& fill) {
        if (!axesMap.contains(o.channel)) return;
        const double s = o.start - globalOffset;
        const double e = o.end - globalOffset;
        if (e < viewStart || s > viewEnd) return;
        const double ds = std::max(s, viewStart);
        const double de = std::min(e, viewEnd);
        const ChartAxes& ca = axesMap[o.channel];
        if (!ca.chart || !ca.xAxis || !ca.yAxis) return;

        auto* upper = new QLineSeries();
        auto* lower = new QLineSeries();
        upper->append({ {ds, ca.yAxis->max()}, {de, ca.yAxis->max()} });
        lower->append({ {ds, ca.yAxis->min()}, {de, ca.yAxis->min()} });
        auto* area = new QAreaSeries(upper, lower);
        area->setBrush(fill);
        area->setPen(Qt::NoPen);
        ca.chart->addSeries(area);
        area->attachAxis(ca.xAxis); area->attachAxis(ca.yAxis);
        m_highlights.append(area);
        };
    for (const ParamOverride& o : m_thresholdOverrides)
        drawOverride(o, QColor(200, 200, 200, 70));   // gray = threshold override
    for (const ParamOverride& o : m_invertOverrides)
        drawOverride(o, QColor(180, 100, 0, 70));   // orange = inversion region
}