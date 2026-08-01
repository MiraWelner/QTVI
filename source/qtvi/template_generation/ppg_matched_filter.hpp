#pragma once
//
// ppg_matched_filter.hpp
//
// Standalone matched-filter PPG pulse detector, per spec: a COMPLEMENTARY
// quality check that runs alongside the primary PPG template pipeline, not
// fused into it. It locates pulses entirely from the PPG signal's own first
// derivative -- no borrowed ECG R-peaks, no dependency on how the primary
// template was built -- so it can flag morphologically inconsistent pulses
// independently of whatever anchoring the main pipeline uses.
//
// Three stages (per spec):
//   (1) Initial pulse locations = maxima of the first derivative (the steep
//       upstroke of each pulse).
//   (2) Construct a template by averaging (column-wise median) the beats
//       sliced at those initial detections.
//   (3) Accept/reject each candidate by NORMALIZED ERROR against the
//       template: ||beat - template|| / ||template|| <= threshold
//       (default 0.05 = 5%). Accepted beats define the final pulse set.
//
// The derivative-peak finder mirrors fp_findpeaks_simple in
// find_foot_pulseox.hpp (kept local here so this header stands alone).
//
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

namespace ppg_matched_filter {

    struct Result {
        std::vector<int>    initialLocs;   // stage 1: derivative-max locations
        std::vector<int>    acceptedLocs;  // stage 3: locations passing the error test
        std::vector<double> templ;         // stage 2: the constructed template
        std::vector<double> normError;     // per-initial-candidate normalized error
    };

    // Local-maxima on a 1-D array (matches fp_findpeaks_simple semantics:
    // strict rise, non-strict fall, plateau-centered).
    inline void findPeaksSimple(const std::vector<double>& d, std::vector<int>& locs) {
        locs.clear();
        if (d.size() < 3) return;
        for (size_t i = 1; i + 1 < d.size(); ++i) {
            if (std::isnan(d[i])) continue;
            if (d[i] > d[i - 1] && d[i] >= d[i + 1]) {
                size_t j = i;
                while (j + 1 < d.size() && d[j] == d[j + 1]) ++j;
                if (j + 1 < d.size() && d[j] > d[j + 1]) {
                    locs.push_back(static_cast<int>(i + (j - i) / 2));
                    i = j;
                }
            }
        }
    }

    // Stage 1: first-derivative maxima. minSep suppresses multiple detections
    // within one pulse (keep the strongest derivative peak per refractory
    // window); pass minSep<=0 to disable.
    inline std::vector<int> derivativePulseLocations(const std::vector<double>& sig,
        int minSep) {
        const int n = static_cast<int>(sig.size());
        std::vector<int> locs;
        if (n < 3) return locs;

        // First difference (central where possible). NaN-safe: a diff touching
        // a NaN is left as NaN so findPeaksSimple skips it.
        std::vector<double> d(n, std::numeric_limits<double>::quiet_NaN());
        for (int i = 1; i < n - 1; ++i) {
            if (std::isnan(sig[i - 1]) || std::isnan(sig[i + 1])) continue;
            d[i] = 0.5 * (sig[i + 1] - sig[i - 1]);
        }

        std::vector<int> raw;
        findPeaksSimple(d, raw);
        if (minSep <= 0 || raw.empty()) return raw;

        // Enforce a refractory minimum separation: within minSep samples keep
        // only the peak with the largest derivative value.
        std::vector<int> kept;
        for (int L : raw) {
            if (kept.empty() || L - kept.back() >= minSep) {
                kept.push_back(L);
            }
            else if (d[L] > d[kept.back()]) {
                kept.back() = L;   // stronger upstroke wins the window
            }
        }
        return kept;
    }

    // Slice a fixed-width beat centered on loc (loc-pre .. loc+post-1).
    // Out-of-range samples -> NaN. Length = pre+post.
    inline std::vector<double> sliceBeat(const std::vector<double>& sig,
        int loc, int pre, int post) {
        const int n = static_cast<int>(sig.size());
        std::vector<double> beat(pre + post, std::numeric_limits<double>::quiet_NaN());
        for (int k = 0; k < pre + post; ++k) {
            const int idx = loc - pre + k;
            if (idx >= 0 && idx < n) beat[k] = sig[idx];
        }
        return beat;
    }

    // Column-wise NaN-skipping median across beats -> template.
    inline std::vector<double> buildTemplate(const std::vector<std::vector<double>>& beats,
        int width) {
        std::vector<double> tmpl(width, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> col;
        for (int c = 0; c < width; ++c) {
            col.clear();
            for (const auto& b : beats)
                if (c < (int)b.size() && !std::isnan(b[c])) col.push_back(b[c]);
            if (col.empty()) continue;
            const size_t mid = col.size() / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            if (col.size() % 2) tmpl[c] = hi;
            else tmpl[c] = 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi);
        }
        return tmpl;
    }

    // Normalized error of one beat vs the template, over columns where both
    // are non-NaN:  ||beat - templ|| / ||templ||.  Returns +inf if the
    // template has zero norm or no overlap.
    inline double normalizedError(const std::vector<double>& beat,
        const std::vector<double>& templ) {
        double num = 0.0, den = 0.0;
        int overlap = 0;
        const int w = std::min(beat.size(), templ.size());
        for (int c = 0; c < w; ++c) {
            if (std::isnan(beat[c]) || std::isnan(templ[c])) continue;
            const double e = beat[c] - templ[c];
            num += e * e;
            den += templ[c] * templ[c];
            ++overlap;
        }
        if (overlap == 0 || den <= 0.0) return std::numeric_limits<double>::infinity();
        return std::sqrt(num / den);
    }

    // Full three-stage matched filter. Runs independently of the primary
    // template pipeline -- callers use this as a complementary QC signal
    // (e.g. flag/count rejected pulses), not as the source of the template
    // itself.
    //   sig       : the PPG signal (one contiguous record).
    //   pre, post : beat window around each detection (samples before/after).
    //   minSep    : refractory min separation between detections (samples).
    //   threshold : normalized-error accept cutoff (default 0.05 = 5%).
    inline Result detect(const std::vector<double>& sig,
        int pre, int post, int minSep,
        double threshold = 0.05) {
        Result r;

        // Stage 1: initial locations from first-derivative maxima.
        r.initialLocs = derivativePulseLocations(sig, minSep);
        if (r.initialLocs.empty()) return r;

        // Stage 2: template from the initial detections.
        const int width = pre + post;
        std::vector<std::vector<double>> beats;
        beats.reserve(r.initialLocs.size());
        for (int L : r.initialLocs) beats.push_back(sliceBeat(sig, L, pre, post));
        r.templ = buildTemplate(beats, width);

        // Stage 3: accept/reject each candidate by normalized error.
        r.normError.reserve(beats.size());
        for (size_t i = 0; i < beats.size(); ++i) {
            const double err = normalizedError(beats[i], r.templ);
            r.normError.push_back(err);
            if (err <= threshold) r.acceptedLocs.push_back(r.initialLocs[i]);
        }
        return r;
    }

} // namespace ppg_matched_filter