/**
 * @file   gui_peak_finder.hpp
 * @brief  Header-only minimal peak finder for visual overlay on the
 *         noise-marking GUI. Operates on raw native-rate (t, v) pairs.
 *
 *         Reference statistics (the amplitude level for the gate, and the mean
 *         inter-peak / R-R interval for blanking) come from a SEPARATE
 *         reference window [refStart, refEnd] -- typically the 10 s preceding
 *         the visible window -- while detection runs over the detection window
 *         [detStart, detEnd]. A wide reference window keeps the gate and
 *         refractory steady even when the visible window is a second or two
 *         wide; if the reference window is empty/too short the detection window
 *         is used as the reference instead.
 *
 *         Gate upper level: when usePeakMedian is true (the preceding-window
 *         case) the gate's top is the MEDIAN amplitude of the R peaks detected
 *         in the reference window -- robust to a lone artifact spike that would
 *         otherwise inflate the absolute max and raise the gate past the real
 *         beats. When usePeakMedian is false (the subsequent-window fallback,
 *         used when there is no preceding window) the absolute max is used.
 *
 *         The threshold and blanking are supplied as PER-PEAK accessors --
 *         functions of (chunk-local) time -- so a regional override applies to
 *         exactly the beats inside its span: each candidate is gated with the
 *         threshold at its own time, and each peak is blanked with the blanking
 *         value at its own time. Callers wanting a constant just pass a
 *         function returning a constant.
 *
 *           - threshold(t):  height gate at threshold(t) of (gateTop - vMin).
 *           - blanking(t):   multiplier on the reference mean R-R; peaks closer
 *                            than blanking(t) * meanRR to the previous kept one
 *                            are dropped.
 *
 *         findPeaks            (ECG): local-maximum detector.
 *         findPeaksDerivative  (PPG/ABP): squared-derivative systolic detector.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <QVector>
#include <QPointF>
#include <QString>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

namespace simple_peak_finder {

    // True if t lies within any [start, end] span.
    inline bool inExcludedSpan(const std::vector<std::pair<double, double>>& spans, double t) {
        for (const auto& s : spans)
            if (t >= s.first && t <= s.second) return true;
        return false;
    }


    // Reference amplitude [min, max], excluding samples inside marked spans.
    // Falls back to the full window when exclusion leaves < 2 samples or a
    // degenerate range. Callers still check vMax > vMin.
    inline std::pair<double, double> referenceRange(
        const QVector<QPointF>& rawPairs, int refFirst, int refLast,
        const std::vector<std::pair<double, double>>& refExcluded) {
        double vMin = 1e300, vMax = -1e300;
        int kept = 0;
        for (int i = refFirst; i <= refLast; ++i) {
            if (inExcludedSpan(refExcluded, rawPairs[i].x())) continue;
            const double v = rawPairs[i].y();
            vMin = std::min(vMin, v); vMax = std::max(vMax, v); ++kept;
        }
        if (kept < 2 || vMax <= vMin) {
            vMin = 1e300; vMax = -1e300;
            for (int i = refFirst; i <= refLast; ++i) {
                const double v = rawPairs[i].y();
                vMin = std::min(vMin, v); vMax = std::max(vMax, v);
            }
        }
        return { vMin, vMax };
    }

    // Drop reference peaks inside marked spans, unless that empties the set.
    inline QVector<QPointF> cleanReferencePeaks(QVector<QPointF> refPeaks,
        const std::vector<std::pair<double, double>>& refExcluded) {
        if (refExcluded.empty()) return refPeaks;
        QVector<QPointF> kept;
        kept.reserve(refPeaks.size());
        for (const QPointF& p : refPeaks)
            if (!inExcludedSpan(refExcluded, p.x())) kept.append(p);
        return kept.isEmpty() ? refPeaks : kept;
    }

    // Mean spacing between consecutive peaks, skipping any interval that overlaps an
    // excluded span (so a removed beat's gap isn't miscounted as one R-R). With no
    // spans this telescopes to the plain meanInterval().
    inline double meanIntervalExcludingSpans(const QVector<QPointF>& peaks,
        const std::vector<std::pair<double, double>>& spans) {
        double sum = 0.0; int cnt = 0;
        for (int i = 1; i < peaks.size(); ++i) {
            const double a = peaks[i - 1].x(), b = peaks[i].x();
            bool straddles = false;
            for (const auto& s : spans)
                if (s.first <= b && a <= s.second) { straddles = true; break; }
            if (straddles) continue;
            sum += b - a; ++cnt;
        }
        return cnt > 0 ? sum / cnt : 0.0;
    }

    // Per-peak parameter accessor: maps a (chunk-local) time to a parameter
    // value (threshold fraction or blanking fraction).
    using ParamFn = std::function<double(double)>;

    // [first, last] indices of rawPairs whose x falls inside [tStart, tEnd].
    // Returns {-1, -1} when nothing falls in range. rawPairs is monotone in x.
    inline std::pair<int, int> indexRange(const QVector<QPointF>& rawPairs,
        double tStart, double tEnd) {
        int firstIdx = -1, lastIdx = -1;
        for (int i = 0; i < rawPairs.size(); ++i) {
            const double t = rawPairs[i].x();
            if (t < tStart) continue;
            if (t > tEnd)   break;
            if (firstIdx < 0) firstIdx = i;
            lastIdx = i;
        }
        return { firstIdx, lastIdx };
    }

    // Mean spacing of a (time-sorted) peak set. 0 if fewer than 2 peaks.
    inline double meanInterval(const QVector<QPointF>& peaks) {
        if (peaks.size() < 2) return 0.0;
        return (peaks.last().x() - peaks.first().x()) / (peaks.size() - 1);
    }

    // Median amplitude (y) of a peak set. 0 if empty.
    inline double medianY(const QVector<QPointF>& peaks) {
        if (peaks.isEmpty()) return 0.0;
        std::vector<double> amps;
        amps.reserve(peaks.size());
        for (const QPointF& p : peaks) amps.push_back(p.y());
        auto m = amps.begin() + amps.size() / 2;
        std::nth_element(amps.begin(), m, amps.end());
        return *m;
    }

    // Drop peaks closer than blanking(t) * refMeanInterval to the one kept
    // before them. refMeanInterval is the reference window's mean R-R, NOT
    // derived from `peaks`. A non-positive interval = no blanking.
    inline QVector<QPointF> applyBlanking(const QVector<QPointF>& peaks,
        const ParamFn& blanking, double refMeanInterval) {
        if (peaks.size() < 2 || refMeanInterval <= 0.0) return peaks;
        QVector<QPointF> out;
        double lastT = -1e300;
        for (const QPointF& c : peaks) {
            const double blank = blanking(c.x()) * refMeanInterval;
            if (blank > 0.0 && c.x() - lastT < blank) continue;
            out.append(c);
            lastT = c.x();
        }
        return out;
    }

    // First index whose x >= val (rawPairs is sorted ascending in x).
    inline int lowerIdx(const QVector<QPointF>& rp, double val) {
        int lo = 0, hi = rp.size();
        while (lo < hi) { const int mid = (lo + hi) >> 1; if (rp[mid].x() < val) lo = mid + 1; else hi = mid; }
        return lo;
    }

    // Walk back from anchor `t`, counting only UNMARKED time, until `need`
    // seconds have been gathered (or we reach 0). Returns the reference start
    // time. `spans` are this channel's marked intervals; may overlap / be unsorted.
    inline double reachBackTime(double t, double need,
        const std::vector<std::pair<double, double>>& spans) {
        std::vector<std::pair<double, double>> m;
        {
            std::vector<std::pair<double, double>> c;
            for (const auto& e : spans) {
                const double a = std::max(0.0, e.first), b = std::min(t, e.second);
                if (a < b) c.push_back({ a, b });
            }
            std::sort(c.begin(), c.end());
            for (const auto& s : c) {
                if (!m.empty() && s.first <= m.back().second) m.back().second = std::max(m.back().second, s.second);
                else m.push_back(s);
            }
        }
        double pos = t;
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

    // Reference stats for the local-max (ECG) detector, taken from the `need`
    // clean seconds immediately BEFORE time `t` (anchored to the beat, never the
    // viewport): the amplitude floor vMin, the gate's upper level (median
    // reference-peak amplitude, or the window max if none), and the mean R-R.
    struct RefStats { double vMin = 0.0, gateTop = 0.0, meanRR = 0.0; bool ok = false; };
    inline RefStats refStatsLocalMax(const QVector<QPointF>& rp, double t, double need,
        const std::vector<std::pair<double, double>>& refExcluded, const ParamFn& threshold) {
        RefStats rs;
        const double rStart = reachBackTime(t, need, refExcluded);
        const int a = lowerIdx(rp, rStart);
        const int b = lowerIdx(rp, t);          // samples strictly before the beat
        if (b - a < 3) return rs;

        double vMin = 1e300, vMax = -1e300; int kept = 0;
        for (int i = a; i < b; ++i) {
            if (inExcludedSpan(refExcluded, rp[i].x())) continue;
            const double v = rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); ++kept;
        }
        if (kept < 2 || vMax <= vMin) {     // exclusion left too little -> use full span
            vMin = 1e300; vMax = -1e300;
            for (int i = a; i < b; ++i) { const double v = rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); }
        }
        if (vMax <= vMin) return rs;

        QVector<QPointF> refPeaks;
        const double spn = vMax - vMin;
        for (int i = a + 1; i < b - 1; ++i) {
            const double tt = rp[i].x();
            if (inExcludedSpan(refExcluded, tt)) continue;
            const double gate = vMin + threshold(tt) * spn;
            const double v = rp[i].y();
            if (v >= gate && v >= rp[i - 1].y() && v > rp[i + 1].y()) refPeaks.append({ tt, v });
        }
        rs.vMin = vMin;
        rs.meanRR = meanIntervalExcludingSpans(refPeaks, refExcluded);
        rs.gateTop = refPeaks.isEmpty() ? vMax : medianY(refPeaks);
        rs.ok = true;
        return rs;
    }

    inline QVector<QPointF> findPeaks(const QVector<QPointF>& rawPairs,
        double detStart, double detEnd,
        double refStart, double refEnd,
        const ParamFn& threshold,
        const ParamFn& blanking,
        bool usePeakMedian,
        const std::vector<std::pair<double, double>>& refExcluded = {},
        const std::vector<std::pair<double, double>>& detExcluded = {})
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 3 || detStart >= detEnd) return out;
        (void)refStart; (void)refEnd; (void)usePeakMedian;   // reference is now per-beat

        // Each candidate is gated by stats from the 10 clean seconds preceding
        // IT (refStatsLocalMax), never anything tied to the viewport -- so a
        // given beat detects identically no matter where the window is panned.
        // The sweep starts a little before the visible window so a beat at the
        // left edge is blanked against its true predecessor, not treated as the
        // first beat in view.
        constexpr double kRefSec = 10.0;
        constexpr double kBlankMargin = 2.0;   // s; >= any plausible blank interval
        const int i0 = std::max(1, lowerIdx(rawPairs, detStart - kBlankMargin));
        const int i1 = std::min(static_cast<int>(rawPairs.size()) - 1,
            lowerIdx(rawPairs, detEnd) + 1);

        QVector<QPointF> prelim;
        std::vector<double> prelimRR;
        for (int i = i0; i < i1; ++i) {
            const double v = rawPairs[i].y();
            if (!(v >= rawPairs[i - 1].y() && v > rawPairs[i + 1].y())) continue;   // local max
            const double t = rawPairs[i].x();
            if (inExcludedSpan(detExcluded, t)) continue;        // no detection in noise spans
            const RefStats rs = refStatsLocalMax(rawPairs, t, kRefSec, refExcluded, threshold);
            if (!rs.ok) continue;
            const double gate = rs.vMin + threshold(t) * (rs.gateTop - rs.vMin);
            if (v >= gate) { prelim.append({ t, v }); prelimRR.push_back(rs.meanRR); }
        }

        // Per-beat blanking: drop a beat within blanking(t) * meanRR(t) of the
        // last kept one. meanRR is the beat's own reference value.
        QVector<QPointF> kept;
        double lastT = -1e300;
        for (int k = 0; k < prelim.size(); ++k) {
            const double blank = blanking(prelim[k].x()) * prelimRR[k];
            if (blank > 0.0 && prelim[k].x() - lastT < blank) continue;
            kept.append(prelim[k]); lastT = prelim[k].x();
        }
        for (const QPointF& p : kept)
            if (p.x() >= detStart && p.x() <= detEnd) out.append(p);
        return out;
    }

    /*
     PPG / ABP peak finder via first-derivative method.

     PPG and ABP have rounded systolic peaks that the threshold-based
     findPeaks() handles poorly. The derivative method finds the steep
     systolic upstroke (a sharp peak in the squared derivative) and walks
     forward to the rounded apex. The amplitude level, the adaptive
     squared-derivative upstroke gate, and the blanking inter-R interval are
     taken from the reference window; the height gate and blanking are
     evaluated per peak. Like findPeaks, the gate's upper level is the median
     detected-peak amplitude when usePeakMedian is set, else the abs max.
     */
    inline QVector<QPointF> findPeaksDerivative(const QVector<QPointF>& rawPairs,
        double detStart, double detEnd,
        double refStart, double refEnd,
        const ParamFn& threshold,
        const ParamFn& blanking,
        bool usePeakMedian,
        const std::vector<std::pair<double, double>>& refExcluded = {},
        const std::vector<std::pair<double, double>>& detExcluded = {})
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 4 || detStart >= detEnd) return out;

        // --- Reference window indices (fall back to det window) ---
        int refFirst, refLast;
        {
            auto rr = indexRange(rawPairs, refStart, refEnd);
            refFirst = rr.first; refLast = rr.second;
            if (refFirst < 0 || refLast - refFirst < 3) {
                auto dd = indexRange(rawPairs, detStart, detEnd);
                refFirst = dd.first; refLast = dd.second;
                if (refFirst < 0 || refLast - refFirst < 3) return out;
            }
        }

        const std::pair<double, double> vr =
            referenceRange(rawPairs, refFirst, refLast, refExcluded);
        const double vMin = vr.first, vMax = vr.second;
        if (vMax <= vMin) return out;

        // Adaptive upstroke gate (90th-pct squared derivative) from reference window.
        std::vector<double> d2ref;
        d2ref.reserve(refLast - refFirst);
        for (int i = refFirst; i < refLast; ++i) {
            if (inExcludedSpan(refExcluded, rawPairs[i].x())
                || inExcludedSpan(refExcluded, rawPairs[i + 1].x())) continue;
            const double dv = rawPairs[i + 1].y() - rawPairs[i].y();
            d2ref.push_back(dv * dv);
        }
        if (d2ref.empty()) {            // exclusion removed everything; use full window
            for (int i = refFirst; i < refLast; ++i) {
                const double dv = rawPairs[i + 1].y() - rawPairs[i].y();
                d2ref.push_back(dv * dv);
            }
        }
        if (d2ref.empty()) return out;
        std::sort(d2ref.begin(), d2ref.end());
        const double d90 = d2ref[(int)(0.9 * d2ref.size())];
        if (d90 <= 0.0) return out;                  // flat reference window
        const double upstrokeGate = 0.6 * d90;

        // Systolic-peak collector. `topLevel` sets the height-gate span.
        auto collect = [&](int a, int b, double topLevel) -> QVector<QPointF> {
            QVector<QPointF> peaks;
            const double sp = topLevel - vMin;
            if (sp <= 0.0) return peaks;
            std::vector<double> d2;
            d2.reserve(b - a);
            for (int i = a; i < b; ++i) {
                const double dv = rawPairs[i + 1].y() - rawPairs[i].y();
                d2.push_back(dv * dv);
            }
            const int N = (int)d2.size();
            for (int k = 1; k < N - 1; ++k) {
                if (d2[k] < upstrokeGate) continue;
                if (!(d2[k] >= d2[k - 1] && d2[k] > d2[k + 1])) continue;

                int rawK = a + k;
                const double tUpstroke = rawPairs[rawK].x();
                int peakIdx = rawK;
                double peakVal = rawPairs[rawK].y();
                for (int j = rawK + 1; j <= b; ++j) {
                    if (rawPairs[j].x() - tUpstroke > 0.2) break;  // 200 ms cap
                    const double v = rawPairs[j].y();
                    if (v >= peakVal) { peakVal = v; peakIdx = j; }
                    else break;
                }
                const double t = rawPairs[peakIdx].x();
                const double gate = vMin + threshold(t) * sp;
                if (peakVal >= gate)
                    peaks.append({ t, peakVal });
            }
            return peaks;
            };

        const QVector<QPointF> refPeaks =
            cleanReferencePeaks(collect(refFirst, refLast, vMax), refExcluded);
        const double refMeanRR = meanIntervalExcludingSpans(refPeaks, refExcluded);
        const double gateTop = (usePeakMedian && !refPeaks.isEmpty())
            ? medianY(refPeaks) : vMax;

        // --- Detection in the visible window, gated against gateTop ---
        auto det = indexRange(rawPairs, detStart, detEnd);
        if (det.first < 0 || det.second - det.first < 3) return out;

        // Don't detect R peaks inside noise/artifact regions (detExcluded);
        // other annotation types detect normally. Filter before blanking so the
        // refractory bridges the gap using the real beats on either side.
        QVector<QPointF> detPeaks = collect(det.first, det.second, gateTop);
        if (!detExcluded.empty())
            detPeaks.erase(std::remove_if(detPeaks.begin(), detPeaks.end(),
                [&](const QPointF& p) { return inExcludedSpan(detExcluded, p.x()); }),
                detPeaks.end());
        return applyBlanking(detPeaks, blanking, refMeanRR);
    }

}  // namespace simple_peak_finder