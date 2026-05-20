/**
 * @file   simple_peak_finder.hpp
 * @brief  Header-only minimal peak finder for visual overlay on the
 *         noise-marking GUI. Operates on raw native-rate (t, v) pairs.
 *
 *         Algorithm:
 *           1. Check polarity (is the signal flipped?)
 *           2. If upright: find every local max above 2/3 of (min..max)
 *           3. If inverted: find every local min below 1/3 of (min..max)
 *
 *         A small refractory suppresses plateau jitter where a notched
 *         peak produces two adjacent local extrema.
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

    struct Params {
        double thresholdFrac = 2.0 / 3.0;  //< gate at this fraction of min->max span
        double refractorySec = 0.05;       //< suppress duplicates within this window
    };

    inline Params paramsFor(const QString& /*label*/) { return {}; }

    inline QVector<QPointF> findPeaks(const QVector<QPointF>& rawPairs,
        double tStart, double tEnd,
        const Params& p)
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 3 || tStart >= tEnd) return out;

        // ---- 1. Locate the index range inside [tStart, tEnd]. ----
        int firstIdx = -1, lastIdx = -1;
        for (int i = 0; i < rawPairs.size(); ++i) {
            const double t = rawPairs[i].x();
            if (t < tStart) continue;
            if (t > tEnd)   break;
            if (firstIdx < 0) firstIdx = i;
            lastIdx = i;
        }
        if (firstIdx < 0 || lastIdx - firstIdx < 2) return out;

        // ---- 2. min, max, median for polarity + gate. ----
        std::vector<double> vals;
        vals.reserve(lastIdx - firstIdx + 1);
        double vMin = 1e300, vMax = -1e300;
        for (int i = firstIdx; i <= lastIdx; ++i) {
            const double v = rawPairs[i].y();
            vals.push_back(v);
            if (v < vMin) vMin = v;
            if (v > vMax) vMax = v;
        }
        if (vMax <= vMin) return out;
        auto mid = vals.begin() + vals.size() / 2;
        std::nth_element(vals.begin(), mid, vals.end());
        const double vMedian = *mid;

        // ---- 3. Is the signal inverted? ----
        // Inverted means the dominant deflections point downward, which
        // shows up as the median being closer to the max than the min.
        const bool inverted = (vMedian - vMin) > (vMax - vMedian);

        // ---- 4. One pass. Branch per polarity. ----
        const double span = vMax - vMin;
        double lastPeakTime = -1e300;

        if (!inverted) {
            // Upright: local max with v >= min + 2/3*span.
            // `>=` on the left of the strict-greater test allows
            // flat-topped digital peaks (consecutive samples at the
            // same ADC value at the apex) to register.
            const double gate = vMin + p.thresholdFrac * span;
            for (int i = firstIdx + 1; i < lastIdx; ++i) {
                const double v = rawPairs[i].y();
                if (v < gate) continue;
                if (!(v >= rawPairs[i - 1].y() && v > rawPairs[i + 1].y())) continue;

                const double t = rawPairs[i].x();
                if (t - lastPeakTime < p.refractorySec) {
                    if (!out.isEmpty() && v > out.last().y()) {
                        out.last() = QPointF(t, v);
                        lastPeakTime = t;
                    }
                    continue;
                }
                out.append(QPointF(t, v));
                lastPeakTime = t;
            }
        }
        else {
            // Inverted: local min with v <= min + 1/3*span (the mirror).
            const double gate = vMin + (1.0 - p.thresholdFrac) * span;
            for (int i = firstIdx + 1; i < lastIdx; ++i) {
                const double v = rawPairs[i].y();
                if (v > gate) continue;
                if (!(v <= rawPairs[i - 1].y() && v < rawPairs[i + 1].y())) continue;

                const double t = rawPairs[i].x();
                if (t - lastPeakTime < p.refractorySec) {
                    if (!out.isEmpty() && v < out.last().y()) {
                        out.last() = QPointF(t, v);
                        lastPeakTime = t;
                    }
                    continue;
                }
                out.append(QPointF(t, v));
                lastPeakTime = t;
            }
        }

        return out;
    }

    /*
     PPG / ABP peak finder via first-derivative method.
    
     PPG and ABP have rounded systolic peaks that the threshold-based
     findPeaks() handles poorly (the dicrotic notch often clears the
     gate, baseline drift moves the reference around, the flat apex
     makes local-max testing unreliable). The derivative method
     sidesteps all of this:
    
       1. Compute the squared first derivative. The systolic upstroke
          -- the steepest part of the rise -- becomes a sharp,
          well-isolated peak in the derivative even though the apex
          itself is rounded.
       2. Find local maxima of the squared derivative above an
          adaptive gate.
       3. For each upstroke, walk forward in the raw signal until
          the signal starts descending. That's the systolic peak.
    
     Robust to baseline drift (derivative kills DC), insensitive to
     the dicrotic notch (notch's upstroke is much weaker than systolic
     upstroke), and naturally produces one detection per beat.
     */
    inline QVector<QPointF> findPeaksDerivative(const QVector<QPointF>& rawPairs,
        double tStart, double tEnd,
        const Params& p)
    {
        QVector<QPointF> out;
        if (rawPairs.size() < 4 || tStart >= tEnd) return out;

        // ---- 1. Locate the index range inside [tStart, tEnd]. ----
        int firstIdx = -1, lastIdx = -1;
        for (int i = 0; i < rawPairs.size(); ++i) {
            const double t = rawPairs[i].x();
            if (t < tStart) continue;
            if (t > tEnd)   break;
            if (firstIdx < 0) firstIdx = i;
            lastIdx = i;
        }
        if (firstIdx < 0 || lastIdx - firstIdx < 3) return out;

        // ---- 2. Compute squared first derivative over the window. ----
        // d2[k] corresponds to the gap between rawPairs[firstIdx+k] and
        // rawPairs[firstIdx+k+1]. Squaring rectifies and emphasizes
        // steep transitions.
        std::vector<double> d2;
        d2.reserve(lastIdx - firstIdx);
        for (int i = firstIdx; i < lastIdx; ++i) {
            const double dv = rawPairs[i + 1].y() - rawPairs[i].y();
            d2.push_back(dv * dv);
        }
        if (d2.empty()) return out;

        // ---- 3. Adaptive gate on the squared derivative. ----
        // Use the 90th percentile of d2 as the reference for "what a
        // strong upstroke looks like in this window." This is robust to
        // a single artifact spike (which would dominate `max`) and to
        // long quiet stretches (which would deflate `mean`).
        std::vector<double> d2sorted = d2;
        std::sort(d2sorted.begin(), d2sorted.end());
        const double d90 = d2sorted[(int)(0.9 * d2sorted.size())];
        if (d90 <= 0.0) return out;          // flat window
        const double upstrokeGate = 0.6 * d90;

        // ---- 4. Walk the derivative; find upstroke peaks. ----
        // For each local max in d2 above the gate, scan forward in the
        // raw signal until it stops climbing -- that's the systolic peak.
        double lastPeakTime = -1e300;
        const int N = (int)d2.size();

        for (int k = 1; k < N - 1; ++k) {
            if (d2[k] < upstrokeGate) continue;
            if (!(d2[k] >= d2[k - 1] && d2[k] > d2[k + 1])) continue;

            // d2[k] is the upstroke. Raw index for the START of that
            // upstroke is firstIdx+k; the peak is shortly AFTER it.
            // Walk forward in the raw signal until the value stops
            // increasing. Cap the search at ~200 ms of samples so we
            // don't run away on a noisy descent.
            int rawK = firstIdx + k;
            const double tUpstroke = rawPairs[rawK].x();

            // Bound the forward scan by both index count and time.
            int peakIdx = rawK;
            double peakVal = rawPairs[rawK].y();
            for (int j = rawK + 1; j <= lastIdx; ++j) {
                if (rawPairs[j].x() - tUpstroke > 0.2) break;  // 200 ms cap
                const double v = rawPairs[j].y();
                if (v >= peakVal) {
                    peakVal = v;
                    peakIdx = j;
                }
                else {
                    break;  // signal turned around; we're past the peak
                }
            }

            const double tPeak = rawPairs[peakIdx].x();
            const double vPeak = rawPairs[peakIdx].y();

            // Refractory: skip if we just emitted one. PPG never beats
            // faster than ~3 Hz, so 0.3 s is a safe floor.
            const double ppgRefractory = std::max(p.refractorySec, 0.3);
            if (tPeak - lastPeakTime < ppgRefractory) {
                // Within refractory: keep the higher of the two.
                if (!out.isEmpty() && vPeak > out.last().y()) {
                    out.last() = QPointF(tPeak, vPeak);
                    lastPeakTime = tPeak;
                }
                continue;
            }

            out.append(QPointF(tPeak, vPeak));
            lastPeakTime = tPeak;
        }

        return out;
    }

}  // namespace simple_peak_finder