#pragma once
//
// pulse_matched_filter.hpp
//
// Matched-filter beat QC for PPG and arterial (ABP/ART/ART_PULM) signals.
// Provides:
//   (1) derivativePulseLocations -- first-derivative-maximum pulse census
//       (steepest-upstroke detector) used by the arterial foot-anchored
//       builder to locate pulses from the signal's own derivative.
//   (2) buildTemplate -- column-wise median template across given beats.
//   (3) normalizedError -- per-beat ||beat - template|| / ||template||.
//       Callers apply their own accept/reject threshold (see
//       create_arterial_templates.hpp).
//
// (Consolidated from the former ppg_matched_filter.hpp. The unused 3-stage
// detect() entry point, its Result struct, and sliceBeat were dropped as
// dead code -- the production pipeline never called them.)
//
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

namespace pulse_matched_filter {

    // Local-maxima on a 1-D array (strict rise, non-strict fall, plateau-
    // centered).
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

    // Stage 1: first-derivative maxima (steepest-upstroke pulse census).
    // minSep suppresses multiple detections within one pulse (keep the
    // strongest derivative peak per refractory window); pass minSep<=0 to
    // disable.
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

} // namespace pulse_matched_filter