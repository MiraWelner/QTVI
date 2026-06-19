/**
 * @file   gui_peak_finder.hpp
 * @brief  A peak finder for running within the noise marking GUI. Two constants:
 *         blanking (minimum time between peaks) and threshold (minimum height),
 *         both per-peak functions of time. Reference stats (gate level + mean
 *         R-R) come from a fixed-length frame: the signal is split into frames
 *         of `previous_seconds_to_train_on` seconds anchored at the chunk start,
 *         and a beat in frame m takes its reference from the nearest EARLIER
 *         frame that contains no marked span (frames with noise in them are
 *         skipped over whole; frame 0 is the floor). Inversion is supplied per
 *         call as `sgn` (+1 upright, -1 inverted) -- the ECG charts set it from
 *         their per-channel "Lead Reversed" checkbox.
 *
 *         findPeaks            (ECG): local-extremum detector (max, or min when inverted).
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

namespace gui_peak_finder {

    inline constexpr double previous_seconds_to_train_on = 5.0;

    // True if t lies within any [start, end] span.
    inline bool inExcludedSpan(const std::vector<std::pair<double, double>>& spans, double t) {
        for (const auto& s : spans)
            if (t >= s.first && t <= s.second) return true;
        return false;
    }

    // Reference amplitude [min, max], excluding samples inside marked spans.
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

    // Mean spacing between consecutive peaks, skipping intervals overlapping a span.
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

    using ParamFn = std::function<double(double)>;

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

    inline double meanInterval(const QVector<QPointF>& peaks) {
        if (peaks.size() < 2) return 0.0;
        return (peaks.last().x() - peaks.first().x()) / (peaks.size() - 1);
    }

    inline double medianY(const QVector<QPointF>& peaks) {
        if (peaks.isEmpty()) return 0.0;
        std::vector<double> amps;
        amps.reserve(peaks.size());
        for (const QPointF& p : peaks) amps.push_back(p.y());
        auto m = amps.begin() + amps.size() / 2;
        std::nth_element(amps.begin(), m, amps.end());
        return *m;
    }

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

    inline int lowerIdx(const QVector<QPointF>& rp, double val) {
        int lo = 0, hi = rp.size();
        while (lo < hi) { const int mid = (lo + hi) >> 1; if (rp[mid].x() < val) lo = mid + 1; else hi = mid; }
        return lo;
    }
    // Merge + sort excluded spans once, clipped to [0, tMax], dropping empties.
    // The result is what reachBackWalk consumes; computing it once per detection
    // call instead of once per candidate beat is the optimization.
    inline std::vector<std::pair<double, double>> mergeExcluded(
        double tMax, const std::vector<std::pair<double, double>>& spans) {
        std::vector<std::pair<double, double>> c;
        c.reserve(spans.size());
        for (const auto& e : spans) {
            const double a = std::max(0.0, e.first), b = std::min(tMax, e.second);
            if (a < b) c.push_back({ a, b });
        }
        std::sort(c.begin(), c.end());
        std::vector<std::pair<double, double>> m;
        for (const auto& s : c) {
            if (!m.empty() && s.first <= m.back().second)
                m.back().second = std::max(m.back().second, s.second);
            else m.push_back(s);
        }
        return m;
    }

    inline bool spanContaining(const std::vector<std::pair<double, double>>& spans,
        double t, double& s, double& e) {
        for (const auto& sp : spans)
            if (t >= sp.first && t <= sp.second) { s = sp.first; e = sp.second; return true; }
        return false;
    }

    // Reference stats for the local-max (ECG) detector, from the reference
    // FRAME for t (see file header). `sgn` (+1 upright, -1 inverted) is supplied
    // by the caller (per-channel checkbox); all amplitude reads run in that frame.
    struct RefStats { double vMin = 0.0, gateTop = 0.0, meanRR = 0.0; bool ok = false; };
    inline RefStats refStatsLocalMax(const QVector<QPointF>& rp, double t, double need,
        const std::vector<std::pair<double, double>>& refExcluded,
        const std::vector<std::pair<double, double>>& withinSpans,
        const ParamFn& threshold,
        double sgn = 1.0) {
        RefStats rs;
        double aS, aE; const bool within = spanContaining(withinSpans, t, aS, aE);
        int a, b;
        if (within) { a = lowerIdx(rp, std::max(aS, t - need)); b = lowerIdx(rp, t); }
        else {
            // Frame-based reference. The chunk is split into fixed `need`-second
            // frames anchored at t = 0. A beat in frame m = floor(t / need)
            // references the nearest EARLIER frame with no marked span in it;
            // any frame containing noise is skipped over whole (no partial use,
            // no reaching back for clean seconds). Frame 0 is the floor -- if
            // nothing clean precedes the beat, frame 0 is used as-is.
            auto frameHasNoise = [&](double lo, double hi) {
                for (const auto& s : refExcluded)
                    if (s.first < hi && lo < s.second) return true;
                return false;
                };
            const long m = static_cast<long>(t / need);   // t >= 0 => floor(t/need)
            long k = m - 1;
            while (k > 0 && frameHasNoise(k * need, (k + 1) * need)) --k;
            if (k < 0) k = 0;
            a = lowerIdx(rp, k * need);
            b = lowerIdx(rp, (k + 1) * need);
        }
        if (b - a < 3) return rs;
        auto excluded = [&](double tt) { return within ? false : inExcludedSpan(refExcluded, tt); };

        double vMin = 1e300, vMax = -1e300; int kept = 0;
        for (int i = a; i < b; ++i) {
            if (excluded(rp[i].x())) continue;
            const double v = sgn * rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); ++kept;
        }
        if (kept < 2 || vMax <= vMin) {     // exclusion left too little -> use full span
            vMin = 1e300; vMax = -1e300;
            for (int i = a; i < b; ++i) { const double v = sgn * rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); }
        }
        if (vMax <= vMin) return rs;

        // The threshold is piecewise-constant in time, changing only at override
        // edges. Probe it at both ends of the reference window: if they match,
        // no edge cuts the window, so reuse that one value for every sample (one
        // lookup instead of one per sample). If they differ, an override starts
        // or ends inside the window -- fall back to the exact per-sample lookup.
        const double thrLo = threshold(rp[a].x());
        const double thrHi = threshold(rp[b - 1].x());
        const bool thrUniform = (thrLo == thrHi);

        QVector<QPointF> refPeaks;
        const double spn = vMax - vMin;
        for (int i = a + 1; i < b - 1; ++i) {
            const double tt = rp[i].x();
            if (excluded(tt)) continue;
            const double thr = thrUniform ? thrLo : threshold(tt);
            const double gate = vMin + thr * spn;
            const double v = sgn * rp[i].y();
            if (v >= gate && v >= sgn * rp[i - 1].y() && v > sgn * rp[i + 1].y()) refPeaks.append({ tt, v });
        }
        rs.vMin = vMin;
        rs.meanRR = within ? meanInterval(refPeaks)
            : meanIntervalExcludingSpans(refPeaks, refExcluded);
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
        double sgn = 1.0,
        const std::vector<std::pair<double, double>>& refExcluded = {},
        const std::vector<std::pair<double, double>>& detExcluded = {},
        const std::vector<std::pair<double, double>>& withinSpans = {})
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 3 || detStart >= detEnd) return out;
        (void)refStart; (void)refEnd; (void)usePeakMedian;   // reference is per-beat
        const std::vector<std::pair<double, double>> refMerged = mergeExcluded(1e300, refExcluded);


        constexpr double kBlankMargin = 2.0;   // s; >= any plausible blank interval
        const int i0 = std::max(1, lowerIdx(rawPairs, detStart - kBlankMargin));
        const int i1 = std::min(static_cast<int>(rawPairs.size()) - 1,
            lowerIdx(rawPairs, detEnd) + 1);

        QVector<QPointF> prelim;
        std::vector<double> prelimRR;
        for (int i = i0; i < i1; ++i) {
            const double yv = sgn * rawPairs[i].y();   // upright-frame amplitude
            if (!(yv >= sgn * rawPairs[i - 1].y() && yv > sgn * rawPairs[i + 1].y())) continue;  // local max upright (= min when inverted)
            const double t = rawPairs[i].x();
            if (inExcludedSpan(detExcluded, t)) continue;
            const RefStats rs = refStatsLocalMax(rawPairs, t, previous_seconds_to_train_on, refMerged, withinSpans, threshold, sgn);
            if (!rs.ok) continue;
            const double gate = rs.vMin + threshold(t) * (rs.gateTop - rs.vMin);
            if (yv >= gate) { prelim.append({ t, rawPairs[i].y() }); prelimRR.push_back(rs.meanRR); }
        }

        // Per-beat blanking: drop a beat within blanking(t) * meanRR(t) of the last kept one.
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

    // Systolic-peak collector over rp[a, b): squared-derivative maxima above
    // upstrokeGate, walked forward <=200 ms to the rounded apex.
    inline QVector<QPointF> collectSystolic(const QVector<QPointF>& rp, int a, int b,
        double upstrokeGate, double vMin, double topLevel, const ParamFn& threshold) {
        QVector<QPointF> peaks;
        const double sp = topLevel - vMin;
        const int n = static_cast<int>(rp.size());
        if (sp <= 0.0 || b > n - 1) b = std::min(b, n - 1);
        if (sp <= 0.0 || b - a < 3) return peaks;
        // Same per-frame threshold hoist as refStatsLocalMax: probe both ends of
        // the window; reuse the value unless an override edge cuts the window.
        const double thrLo = threshold(rp[a].x());
        const double thrHi = threshold(rp[b - 1].x());
        const bool thrUniform = (thrLo == thrHi);
        for (int k = a + 1; k < b - 1; ++k) {
            const double dkm = rp[k].y() - rp[k - 1].y();
            const double dk = rp[k + 1].y() - rp[k].y();
            const double dkp = rp[k + 2].y() - rp[k + 1].y();
            const double d2km = dkm * dkm, d2k = dk * dk, d2kp = dkp * dkp;
            if (d2k < upstrokeGate) continue;
            if (!(d2k >= d2km && d2k > d2kp)) continue;
            const double tUp = rp[k].x();
            int apexIdx = k; double apexVal = rp[k].y();
            for (int j = k + 1; j <= b && j < n; ++j) {
                if (rp[j].x() - tUp > 0.2) break;
                const double v = rp[j].y();
                if (v >= apexVal) { apexVal = v; apexIdx = j; }
                else break;
            }
            const double t = rp[apexIdx].x();
            const double thr = thrUniform ? thrLo : threshold(t);
            if (apexVal >= vMin + thr * sp) peaks.append({ t, apexVal });
        }
        return peaks;
    }

    struct RefStatsD { double vMin = 0.0, gateTop = 0.0, meanRR = 0.0, upstrokeGate = 0.0; bool ok = false; };
    inline RefStatsD refStatsDeriv(const QVector<QPointF>& rp, double t, double need,
        const std::vector<std::pair<double, double>>& refExcluded,
        const std::vector<std::pair<double, double>>& withinSpans,
        const ParamFn& threshold) {
        RefStatsD rs;
        double aS, aE; const bool within = spanContaining(withinSpans, t, aS, aE);
        int a, b;
        if (within) { a = lowerIdx(rp, std::max(aS, t - need)); b = lowerIdx(rp, t); }
        else {
            // Frame-based reference (see refStatsLocalMax for the full rationale):
            // the nearest earlier frame with no marked span in it; noisy frames
            // are skipped whole; frame 0 is the floor.
            auto frameHasNoise = [&](double lo, double hi) {
                for (const auto& s : refExcluded)
                    if (s.first < hi && lo < s.second) return true;
                return false;
                };
            const long m = static_cast<long>(t / need);   // t >= 0 => floor(t/need)
            long k = m - 1;
            while (k > 0 && frameHasNoise(k * need, (k + 1) * need)) --k;
            if (k < 0) k = 0;
            a = lowerIdx(rp, k * need);
            b = lowerIdx(rp, (k + 1) * need);
        }
        if (b - a < 4) return rs;
        auto excluded = [&](double tt) { return within ? false : inExcludedSpan(refExcluded, tt); };

        double vMin = 1e300, vMax = -1e300; int kept = 0;
        for (int i = a; i < b; ++i) {
            if (excluded(rp[i].x())) continue;
            const double v = rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); ++kept;
        }
        if (kept < 2 || vMax <= vMin) {
            vMin = 1e300; vMax = -1e300;
            for (int i = a; i < b; ++i) { const double v = rp[i].y(); vMin = std::min(vMin, v); vMax = std::max(vMax, v); }
        }
        if (vMax <= vMin) return rs;

        std::vector<double> d2;
        d2.reserve(b - a);
        for (int i = a; i < b - 1; ++i) {
            if (excluded(rp[i].x()) || excluded(rp[i + 1].x())) continue;
            const double dv = rp[i + 1].y() - rp[i].y(); d2.push_back(dv * dv);
        }
        if (d2.empty()) {
            for (int i = a; i < b - 1; ++i) { const double dv = rp[i + 1].y() - rp[i].y(); d2.push_back(dv * dv); }
        }
        if (d2.empty()) return rs;
        std::sort(d2.begin(), d2.end());
        const double d90 = d2[(int)(0.9 * d2.size())];
        if (d90 <= 0.0) return rs;
        rs.upstrokeGate = 0.6 * d90;

        const QVector<QPointF> refRaw = collectSystolic(rp, a, b, rs.upstrokeGate, vMin, vMax, threshold);
        const QVector<QPointF> refPeaks = within ? refRaw : cleanReferencePeaks(refRaw, refExcluded);
        rs.vMin = vMin;
        rs.meanRR = within ? meanInterval(refPeaks)
            : meanIntervalExcludingSpans(refPeaks, refExcluded);
        rs.gateTop = refPeaks.isEmpty() ? vMax : medianY(refPeaks);
        rs.ok = true;
        return rs;
    }

    inline QVector<QPointF> findPeaksDerivative(const QVector<QPointF>& rawPairs,
        double detStart, double detEnd,
        double refStart, double refEnd,
        const ParamFn& threshold,
        const ParamFn& blanking,
        bool usePeakMedian,
        const std::vector<std::pair<double, double>>& refExcluded = {},
        const std::vector<std::pair<double, double>>& detExcluded = {},
        const std::vector<std::pair<double, double>>& withinSpans = {})
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 4 || detStart >= detEnd) return out;
        (void)refStart; (void)refEnd; (void)usePeakMedian;

        constexpr double kBlankMargin = 2.0;
        const int n = static_cast<int>(rawPairs.size());
        const int i0 = std::max(1, lowerIdx(rawPairs, detStart - kBlankMargin));
        const int i1 = std::min(n - 2, lowerIdx(rawPairs, detEnd) + 1);

        QVector<QPointF> prelim;
        std::vector<double> prelimRR;
        for (int k = i0; k < i1; ++k) {
            const double dkm = rawPairs[k].y() - rawPairs[k - 1].y();
            const double dk = rawPairs[k + 1].y() - rawPairs[k].y();
            const double dkp = rawPairs[k + 2].y() - rawPairs[k + 1].y();
            const double d2km = dkm * dkm, d2k = dk * dk, d2kp = dkp * dkp;
            if (!(d2k >= d2km && d2k > d2kp)) continue;

            const double tUp = rawPairs[k].x();
            int apexIdx = k; double apexVal = rawPairs[k].y();
            for (int j = k + 1; j < n; ++j) {
                if (rawPairs[j].x() - tUp > 0.2) break;
                const double v = rawPairs[j].y();
                if (v >= apexVal) { apexVal = v; apexIdx = j; }
                else break;
            }
            const double t = rawPairs[apexIdx].x();
            if (inExcludedSpan(detExcluded, t)) continue;
            const RefStatsD rs = refStatsDeriv(rawPairs, t, previous_seconds_to_train_on, refExcluded, withinSpans, threshold);
            if (!rs.ok) continue;
            if (d2k < rs.upstrokeGate) continue;
            const double gate = rs.vMin + threshold(t) * (rs.gateTop - rs.vMin);
            if (apexVal >= gate) { prelim.append({ t, apexVal }); prelimRR.push_back(rs.meanRR); }
        }

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

}  // namespace gui_peak_finder