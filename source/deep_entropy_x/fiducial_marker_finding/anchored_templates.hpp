#pragma once
#include "template_marking_gui\feature_marks.hpp"   // AnchorType, AnchorLocator, make_anchor_locator
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

struct AnchoredTemplate {
    AnchorType anchor;
    std::vector<double> mean;   // column-wise NaN-skipping median
    std::vector<double> iqr;    // column-wise IQR (q3-q1); robust, not std
    // B2 focus-mode stats (spec's confidence band uses true mean +/- 1.96*se,
    // se = sd/sqrt(n)). Kept SEPARATE from the median/iqr above so the main
    // view's median center line is unaffected -- these are additive, only the
    // focus panel reads them.
    std::vector<double> arithMean;   // column-wise NaN-skipping arithmetic mean
    std::vector<double> sd;          // column-wise NaN-skipping sample sd (ddof=1)
    std::vector<int>    colCount;    // per-column count of non-NaN beats (the n in se)
    int refColumn = -1;
    int nBeats = 0;
};

// Re-anchor retained clean beats onto `type`'s landmark and average.
// width = beats' shared width; refColumn = column the landmark lands on
// (pass locate(ref_beat_of_median_length) so the reference beat doesn't move).
inline AnchoredTemplate makeAnchoredAverage(
    const std::vector<std::vector<double>>& beats,
    const AnchorLocator& locate,
    AnchorType type, int width, int refColumn)
{
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> aligned;
    aligned.reserve(beats.size());
    for (const auto& beat : beats) {
        // Sub-sample locator, whole-sample shift: rounds, as the previous
        // int-returning locator did inside the finder.
        const int lm = static_cast<int>(std::lround(locate(beat)));
        if (lm < 0) continue;
        std::vector<double> row(width, NaN);
        const int shift = refColumn - lm;
        for (int i = 0; i < (int)beat.size(); ++i) {
            const int j = i + shift;
            if (j >= 0 && j < width) row[j] = beat[i];
        }
        aligned.push_back(std::move(row));
    }

    AnchoredTemplate out;
    out.anchor = type;
    out.refColumn = refColumn;
    out.nBeats = static_cast<int>(aligned.size());
    out.mean.assign(width, NaN);
    out.iqr.assign(width, 0.0);
    out.arithMean.assign(width, NaN);
    out.sd.assign(width, 0.0);
    out.colCount.assign(width, 0);

    std::vector<double> col;
    col.reserve(aligned.size());
    for (int c = 0; c < width; ++c) {
        col.clear();
        for (const auto& r : aligned)
            if (!std::isnan(r[c])) col.push_back(r[c]);
        const size_t nc = col.size();
        if (nc == 0) continue;

        // Arithmetic mean + sample sd (ddof=1) for the focus-mode band, over
        // the same non-NaN column values used for the median below.
        out.colCount[c] = static_cast<int>(nc);
        double sum = 0.0;
        for (double v : col) sum += v;
        const double m = sum / static_cast<double>(nc);
        out.arithMean[c] = m;
        if (nc >= 2) {
            double ss = 0.0;
            for (double v : col) { const double d = v - m; ss += d * d; }
            out.sd[c] = std::sqrt(ss / static_cast<double>(nc - 1));
        }

        const size_t mid = nc / 2;
        std::nth_element(col.begin(), col.begin() + mid, col.end());
        const double hi = col[mid];
        out.mean[c] = (nc % 2 == 0)
            ? 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi)
            : hi;
        if (nc >= 2) {
            const size_t q1i = nc / 4;
            const size_t q3i = (3 * nc) / 4;
            std::nth_element(col.begin(), col.begin() + q1i, col.end());
            const double q1 = col[q1i];
            std::nth_element(col.begin() + q1i, col.begin() + q3i, col.end());
            const double q3 = col[q3i];
            out.iqr[c] = q3 - q1;
        }
    }
    return out;
}

// Per-bin lazy cache (CHAOS spec Step 4). Point cleanBeats at this bin's
// retained [beat][sample] matrix; ref_beat_of_median_length is the median-length
// reference beat (same object the Q pass uses). width = cleanBeats.front().size().
struct AnchorCache {
    const std::vector<std::vector<double>>* cleanBeats = nullptr;
    const std::vector<double>* ref_beat_of_median_length = nullptr;
    int r_col = -1, width = 0;
    double fs = 0.0;
    std::map<AnchorType, AnchoredTemplate> cache;

    const AnchoredTemplate& get(AnchorType a) {
        auto it = cache.find(a);
        if (it != cache.end()) return it->second;
        auto locate = make_anchor_locator(a, r_col, fs);
        const int refColumn = static_cast<int>(std::lround(
            locate(*ref_beat_of_median_length)));
        AnchoredTemplate t = makeAnchoredAverage(*cleanBeats, locate, a, width, refColumn);
        return cache.emplace(a, std::move(t)).first->second;
    }
    void invalidate(AnchorType a) { cache.erase(a); }   // on landmark edit
};