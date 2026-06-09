/**
 * @file   simple_peak_finder.hpp
 * @brief  Header-only minimal peak finder for visual overlay on the
 *         noise-marking GUI. Operates on raw native-rate (t, v) pairs.
 *
 *         Both finders share the same two post-filters:
 *           - threshold:       height gate at `threshold` of the vMin..vMax span.
 *           - blanking_period: multiplier on the mean inter-peak interval;
 *                              peaks closer than that to the previous kept one
 *                              are dropped.
 *
 *         findPeaks            (ECG): local-extremum detector.
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
#include <vector>

namespace simple_peak_finder {

    struct peakfinding_features {
        double threshold; //the span between min and max is calculated, multiplied by threshold, and R peaks must be above threshold
        double blanking_period; //this is multiplied by the mean r-r distance, and a new R peak cannot occur the blanking period distance after previous one
    };

    inline peakfinding_features paramsFor(const QString& /*label*/) { return {}; }

    // Drop peaks closer than blanking_period * (mean inter-peak interval) to the
    // one kept before them. Shared by both finders.
    inline QVector<QPointF> applyBlanking(const QVector<QPointF>& peaks, double blankingFrac) {
        QVector<QPointF> out;
        if (peaks.size() < 2) return peaks;
        const double blank = blankingFrac
            * (peaks.last().x() - peaks.first().x()) / (peaks.size() - 1);
        double lastT = -1e300;
        for (const QPointF& c : peaks) {
            if (c.x() - lastT < blank) continue;
            out.append(c);
            lastT = c.x();
        }
        return out;
    }

    inline QVector<QPointF> findPeaks(const QVector<QPointF>& rawPairs,
        double tStart, double tEnd,
        const peakfinding_features& p)
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 3 || tStart >= tEnd) return out;

        // Locate the index range inside [tStart, tEnd].
        int firstIdx = -1, lastIdx = -1;
        for (int i = 0; i < rawPairs.size(); ++i) {
            const double t = rawPairs[i].x();
            if (t < tStart) continue;
            if (t > tEnd)   break;
            if (firstIdx < 0) firstIdx = i;
            lastIdx = i;
        }
        if (firstIdx < 0 || lastIdx - firstIdx < 2) return out;

        // min / max / median (median only decides polarity).
        std::vector<double> vals;
        vals.reserve(lastIdx - firstIdx + 1);
        double vMin = 1e300, vMax = -1e300;
        for (int i = firstIdx; i <= lastIdx; ++i) {
            const double v = rawPairs[i].y();
            vals.push_back(v);
            vMin = std::min(vMin, v);
            vMax = std::max(vMax, v);
        }
        if (vMax <= vMin) return out;
        auto mid = vals.begin() + vals.size() / 2;
        std::nth_element(vals.begin(), mid, vals.end());
        const bool wantMax = (*mid - vMin) <= (vMax - *mid);   // upright if median near min

        // Gate at `threshold` of the full span; collect local extrema past it.
        const double span = vMax - vMin;
        const double gate = wantMax ? vMin + p.threshold * span
            : vMax - p.threshold * span;
        QVector<QPointF> peaks;
        for (int i = firstIdx + 1; i < lastIdx; ++i) {
            const double v = rawPairs[i].y();
            const bool hit = wantMax
                ? (v >= gate && v >= rawPairs[i - 1].y() && v > rawPairs[i + 1].y())
                : (v <= gate && v <= rawPairs[i - 1].y() && v < rawPairs[i + 1].y());
            if (hit) peaks.append({ rawPairs[i].x(), v });
        }

        return applyBlanking(peaks, p.blanking_period);
    }

    /*
     PPG / ABP peak finder via first-derivative method.

     PPG and ABP have rounded systolic peaks that the threshold-based
     findPeaks() handles poorly. The derivative method finds the steep
     systolic upstroke (a sharp peak in the squared derivative) and walks
     forward to the rounded apex. The same threshold (height gate) and
     blanking (mean-interval) post-filters as findPeaks are then applied.
     */
    inline QVector<QPointF> findPeaksDerivative(const QVector<QPointF>& rawPairs,
        double tStart, double tEnd,
        const peakfinding_features& p)
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 4 || tStart >= tEnd) return out;

        // Locate the index range inside [tStart, tEnd].
        int firstIdx = -1, lastIdx = -1;
        for (int i = 0; i < rawPairs.size(); ++i) {
            const double t = rawPairs[i].x();
            if (t < tStart) continue;
            if (t > tEnd)   break;
            if (firstIdx < 0) firstIdx = i;
            lastIdx = i;
        }
        if (firstIdx < 0 || lastIdx - firstIdx < 3) return out;

        // min / max for the height threshold.
        double vMin = 1e300, vMax = -1e300;
        for (int i = firstIdx; i <= lastIdx; ++i) {
            const double v = rawPairs[i].y();
            vMin = std::min(vMin, v);
            vMax = std::max(vMax, v);
        }
        if (vMax <= vMin) return out;
        const double gate = vMin + p.threshold * (vMax - vMin);

        // Squared first derivative over the window.
        std::vector<double> d2;
        d2.reserve(lastIdx - firstIdx);
        for (int i = firstIdx; i < lastIdx; ++i) {
            const double dv = rawPairs[i + 1].y() - rawPairs[i].y();
            d2.push_back(dv * dv);
        }
        if (d2.empty()) return out;

        // Adaptive gate on the squared derivative (90th pct).
        std::vector<double> d2sorted = d2;
        std::sort(d2sorted.begin(), d2sorted.end());
        const double d90 = d2sorted[(int)(0.9 * d2sorted.size())];
        if (d90 <= 0.0) return out;          // flat window
        const double upstrokeGate = 0.6 * d90;

        // Walk the derivative; collect systolic peaks past the height gate.
        QVector<QPointF> peaks;
        const int N = (int)d2.size();
        for (int k = 1; k < N - 1; ++k) {
            if (d2[k] < upstrokeGate) continue;
            if (!(d2[k] >= d2[k - 1] && d2[k] > d2[k + 1])) continue;

            int rawK = firstIdx + k;
            const double tUpstroke = rawPairs[rawK].x();
            int peakIdx = rawK;
            double peakVal = rawPairs[rawK].y();
            for (int j = rawK + 1; j <= lastIdx; ++j) {
                if (rawPairs[j].x() - tUpstroke > 0.2) break;  // 200 ms cap
                const double v = rawPairs[j].y();
                if (v >= peakVal) { peakVal = v; peakIdx = j; }
                else break;
            }
            if (peakVal >= gate)
                peaks.append({ rawPairs[peakIdx].x(), peakVal });
        }

        return applyBlanking(peaks, p.blanking_period);
    }

}  // namespace simple_peak_finder