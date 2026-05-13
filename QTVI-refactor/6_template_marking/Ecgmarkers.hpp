#pragma once
//
// Simple ECG fiducial detectors for the template-marking viewer.
//
// All three operate on a single beat-aligned ECG template, with the R-peak
// position already known (from ChannelMethodTemplate::alignment_point). Each
// returns -1 if the template is too short or the search fails.
//
// These are deliberately simple and deterministic. They're meant to seed the
// UI with a "close enough" starting point that the user drags into final
// position. Robustness over cleverness.
//
// All indices are sample indices into the template vector.
//

#include <vector>
#include <cmath>
#include <algorithm>

namespace ecg_markers {

    // ---- Helpers --------------------------------------------------------

    namespace detail {
        // Centered 3-point slope. Edges clipped.
        inline double slope_at(const std::vector<double>& v, int i) {
            const int N = static_cast<int>(v.size());
            if (i <= 0)        return v[1] - v[0];
            if (i >= N - 1)    return v[N - 1] - v[N - 2];
            return 0.5 * (v[i + 1] - v[i - 1]);
        }

        // Smoothed slope over a small window — less twitchy on real templates.
        inline double slope_window(const std::vector<double>& v, int i, int w = 3) {
            const int N = static_cast<int>(v.size());
            int lo = std::max(0, i - w);
            int hi = std::min(N - 1, i + w);
            if (hi <= lo) return 0.0;
            return (v[hi] - v[lo]) / static_cast<double>(hi - lo);
        }

        // Refine an R-peak hint to the actual argmax in +/-30 samples.
        // alignment_point_raw is reliable but not pixel-exact on the
        // visual R-spike; this snaps to the local maximum so downstream
        // logic ("Q is before R", "T-wave starts >80ms after R") uses the
        // visually-correct R.
        inline int refine_r(const std::vector<double>& beat, int hint) {
            const int N = static_cast<int>(beat.size());
            if (hint < 0 || hint >= N) return hint;
            const int lo = std::max(0, hint - 30);
            const int hi = std::min(N - 1, hint + 30);
            int rIdx = lo;
            double rVal = beat[lo];
            for (int i = lo; i <= hi; ++i) {
                if (beat[i] > rVal) { rVal = beat[i]; rIdx = i; }
            }
            return rIdx;
        }
    }

    // ---- Q-wave begin ---------------------------------------------------
    //
    // The Q-wave is a small negative dip immediately before R. Q-begin is
    // the bottom of that dip -- the local minimum just left of R.
    //
    // Approach: find the lowest sample in a short window before R. Simple
    // and robust against small wiggles in the QRS upstroke.
    // ---- Q-wave begin ---------------------------------------------------
    //
    // alignment_point_raw consistently lands ~10 samples before the visual
    // R-peak in our templates -- it points at the foot of the QRS upstroke,
    // which is exactly what we want for Q-begin. Use it directly.
    inline int detect_q_begin(const std::vector<double>& beat, int rPeakHint) {
        const int N = static_cast<int>(beat.size());
        if (rPeakHint < 0 || rPeakHint >= N) return -1;
        return rPeakHint;
    }

    // ---- T-wave begin ---------------------------------------------------
    //
    // Forward from R: skip the QRS downstroke through the S trough, then
    // find where the trace starts climbing again toward the T hump.
    //
    // Approach: locate the S trough (first local minimum after R), then
    // walk forward until slope turns clearly positive.
    inline int detect_t_begin(const std::vector<double>& beat, int rPeakHint) {
        const int N = static_cast<int>(beat.size());
        if (rPeakHint < 0 || rPeakHint >= N - 4) return -1;

        const int rPeak = detail::refine_r(beat, rPeakHint);

        // T-wave doesn't physiologically start within ~80ms of R at normal
        // heart rates. At 1 kHz that's 80 samples.
        const int minDelay = std::min(80, (N - rPeak) / 3);
        const int afterQRS = rPeak + minDelay;
        if (afterQRS >= N - 4) return -1;

        // Locate S trough inside [rPeak, afterQRS] — the QRS-region minimum.
        int sIdx = rPeak;
        double sVal = beat[rPeak];
        for (int i = rPeak + 1; i <= afterQRS; ++i) {
            if (beat[i] < sVal) { sVal = beat[i]; sIdx = i; }
        }

        const double qrsAmp = std::max(1e-9, beat[rPeak] - sVal);
        const double thresh = 0.02 * qrsAmp;

        for (int i = afterQRS; i < N - 2; ++i) {
            if (detail::slope_window(beat, i, 3) > thresh) return i;
        }
        return -1;
    }

    // ---- T-wave end (tangent method, simplified) ------------------------
    //
    // Find the T-peak, then on its descending limb find the inflection
    // (most negative slope). The tangent at that inflection extrapolated
    // to baseline gives T-end. We approximate baseline as the median of
    // the last 10% of the beat, which is usually the isoelectric tail.
    //
    // Falls back to "last sample where slope flattens" if the tangent
    // intersection lies beyond the beat or before the inflection.
    inline int detect_t_end(const std::vector<double>& beat, int rPeakHint) {
        const int N = static_cast<int>(beat.size());
        if (rPeakHint < 0 || rPeakHint >= N - 10) return -1;

        const int rPeak = detail::refine_r(beat, rPeakHint);

        const int minDelay = std::min(80, (N - rPeak) / 3);
        const int searchStart = rPeak + minDelay;
        if (searchStart >= N - 4) return -1;

        // T peak: highest sample in the post-ST region.
        // (Assumes upright T-wave; inverted-T leads would need abs().)
        int tPeak = searchStart;
        double tVal = beat[searchStart];
        for (int i = searchStart; i < N; ++i) {
            if (beat[i] > tVal) { tVal = beat[i]; tPeak = i; }
        }
        if (tPeak >= N - 4) return -1;

        // Baseline: median of the iso-electric ST window between R and the
        // T-peak. That window sits AFTER the S-trough and BEFORE the T-peak
        // starts climbing -- it's the flattest region.
        const int isoStart = rPeak + (minDelay / 2);
        const int isoEnd = std::max(isoStart + 3, tPeak - (tPeak - rPeak) / 4);
        std::vector<double> iso(beat.begin() + isoStart,
            beat.begin() + std::min(isoEnd, N));
        double baseline = tVal;  // fallback
        if (iso.size() >= 3) {
            std::nth_element(iso.begin(),
                iso.begin() + iso.size() / 2, iso.end());
            baseline = iso[iso.size() / 2];
        }

        // Walk forward from T peak until the trace returns to within 10%
        // of (T-peak - baseline) above baseline -- that's T-end.
        const double tHeight = std::max(1e-9, tVal - baseline);
        const double endThresh = baseline + 0.10 * tHeight;
        for (int i = tPeak + 1; i < N - 1; ++i) {
            if (beat[i] <= endThresh) return i;
        }
        return N - 1;
    }

} // namespace ecg_markers