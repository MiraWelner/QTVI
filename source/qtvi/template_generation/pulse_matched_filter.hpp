#pragma once
//
// pulse_matched_filter.hpp
//
// Matched-filter beat QC for PPG and arterial (ABP/ART/ART_PULM) signals.
// Beats are located elsewhere (R-anchored, via alignment::extract_ppg_
// beats_and_align) and handed in already sliced; this file only provides
// stages 2-3 of the original 3-stage spec:
//   (2) Construct a template by averaging (column-wise median) the given
//       beats.
//   (3) Per-beat NORMALIZED ERROR against that template:
//       ||beat - template|| / ||template||. Callers apply their own
//       accept/reject threshold (see create_arterial_templates.hpp).
//
// Stage 1 (self-detecting pulses from the signal's own first derivative)
// and the combined detect() entry point were removed as dead code -- the
// production pipeline never calls them; beats always arrive pre-located.
//
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

namespace pulse_matched_filter {

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