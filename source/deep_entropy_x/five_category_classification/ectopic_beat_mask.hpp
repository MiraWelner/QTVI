#pragma once
/**
 * @file   ectopic_beat_mask.hpp
 * @brief  Which beats must stay out of the template median, so the trace the
 *         operator marks is not an average of sinus and ectopy.
 *
 *         THE RHYTHM VERDICT IS NOT COMPUTED HERE. It is assigned in
 *         alignment.hpp, after the slice and before any pruning, and compacted
 *         through every apply_mask() alongside `beats` -- so by the time this
 *         file sees the surviving beats, each one already carries its verdict
 *         in the same index space. This file receives that as `rhythmFlagged`
 *         and adds only the morphology route.
 *
 *         That ordering matters because alignment's FIRST Tukey pass rejects
 *         on RR LENGTH, and a premature beat is short by definition. Ask "was
 *         this beat premature" after that pass and you are asking about the
 *         beats that survived not being premature. Ask before it and every
 *         sliced beat has an answer that outlives the pruning.
 *
 *         TWO ROUTES OUT, OR-ED, both from Section 4.6:
 *
 *           RHYTHM. Premature -- RR(t) < 0.80 * median of the trailing ten --
 *           OR voted premature by the 5-of-8 rule. The vote widens the filter
 *           rather than narrowing it: inside a run the trailing median has
 *           itself gone short, so an intra-run beat stops reading as premature
 *           on its own and only its neighbours give it away. An isolated PVC
 *           needs no vote; the filter already has it.
 *
 *           MORPHOLOGY. r < 0.85 against the pass-1 median. Catches ectopy
 *           that is not premature at all: an interpolated PVC, a fusion beat,
 *           bigeminy settled into a regular rhythm. Needs a template, so it
 *           necessarily runs after the first median -- unlike the rhythm
 *           route, which needs nothing but the peaks.
 *
 *         PACs GO TOO. 4.6 excludes the whole abnormal-rhythm category. A PAC
 *         is premature and supraventricular, so it would pass the correlation
 *         floor and only the rhythm route takes it.
 *
 *         FAIL-SAFES DIFFER BY ROUTE. The morphology route is capped at
 *         maxDropFrac: a third of beats scoring below 0.85 against their own
 *         median means a bad median or failed alignment, and masking on that
 *         empties the template instead of cleaning it. The rhythm route is not
 *         capped that way -- in atrial fibrillation a large share of beats
 *         really are premature, and that is the reading, not a runaway. Both
 *         floor at minKeepBeats.
 *
 * @date   2026-08-25
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace ectopic_mask {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    struct Params {
        double corrFloor = 0.85;   // 4.6, verbatim (ECG)
        double maxDropFrac = 0.30;   // morphology route only
        int    minKeepBeats = 8;
        int    minBeatsToRun = 12;
    };

    namespace detail {
        inline double medianOf(std::vector<double> v) {
            v.erase(std::remove_if(v.begin(), v.end(),
                [](double x) { return !std::isfinite(x); }), v.end());
            if (v.empty()) return kNaN;
            std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
            return v[v.size() / 2];
        }

        inline double pearson(const std::vector<double>& a,
            const std::vector<double>& b)
        {
            const std::size_t n = std::min(a.size(), b.size());
            double sa = 0.0, sb = 0.0;
            std::size_t m = 0;
            for (std::size_t i = 0; i < n; ++i) {
                if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
                sa += a[i]; sb += b[i]; ++m;
            }
            if (m < 3) return kNaN;
            const double ma = sa / m, mb = sb / m;
            double num = 0.0, va = 0.0, vb = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                if (!std::isfinite(a[i]) || !std::isfinite(b[i])) continue;
                const double da = a[i] - ma, db = b[i] - mb;
                num += da * db; va += da * da; vb += db * db;
            }
            if (va <= 0.0 || vb <= 0.0) return kNaN;
            return num / std::sqrt(va * vb);
        }
    } // namespace detail

    // ---------------------------------------------------------------------
    // Stage 2: the mask over the surviving beats
    // ---------------------------------------------------------------------
    struct Result {
        std::vector<char>   ectopic;     // 1 = keep out of the median
        std::vector<double> corr;        // per beat, against the pass-1 median
        int    nEctopic = 0, nRhythm = 0, nMorph = 0;
        bool   applied = false;
        bool   morphRouteUsed = false;
        const char* note = "not_run";
    };

    // `beats`        the surviving, aligned, NaN-padded beat matrix.
    // `refTemplate`  the pass-1 median over those beats.
    // `rhythmFlagged` per SURVIVING beat, from alignment's compacted
    //                 premature/voted vectors. Empty disables the rhythm
    //                 route, leaving morphology alone.
    inline Result build(const std::vector<std::vector<double>>& beats,
        const std::vector<double>& refTemplate,
        const std::vector<char>& rhythmFlagged,
        const Params& p = {})
    {
        Result r;
        const std::size_t n = beats.size();
        r.ectopic.assign(n, 0);
        r.corr.assign(n, kNaN);

        if (n == 0 || refTemplate.empty()) { r.note = "no_beats"; return r; }
        if (static_cast<int>(n) < p.minBeatsToRun) {
            r.note = "too_few_beats_to_gate";
            return r;
        }

        // ---- route 1: the rhythm flags, computed pre-pruning -------------
        const bool haveRhythm = (rhythmFlagged.size() == n);
        if (haveRhythm) {
            for (std::size_t i = 0; i < n; ++i)
                if (rhythmFlagged[i]) { r.ectopic[i] = 1; ++r.nRhythm; }
        }

        // ---- route 2: the 4.6 morphology floor --------------------------
        int morphOnly = 0;
        std::vector<char> morphFlag(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            r.corr[i] = detail::pearson(beats[i], refTemplate);
            if (!std::isfinite(r.corr[i])) continue;   // unscorable stays in
            if (r.corr[i] < p.corrFloor) {
                morphFlag[i] = 1;
                if (!r.ectopic[i]) ++morphOnly;
            }
        }
        // The cap is on what the morphology route ADDS, so a record already
        // largely flagged on rhythm (AF) does not switch it off.
        if (morphOnly > 0
            && static_cast<double>(morphOnly) <= p.maxDropFrac * static_cast<double>(n))
        {
            r.morphRouteUsed = true;
            for (std::size_t i = 0; i < n; ++i)
                if (morphFlag[i] && !r.ectopic[i]) { r.ectopic[i] = 1; ++r.nMorph; }
        }

        for (std::size_t i = 0; i < n; ++i) r.nEctopic += r.ectopic[i];

        if (r.nEctopic == 0) {
            r.note = haveRhythm ? "no_ectopic" : "no_ectopic_no_rhythm_input";
            return r;
        }
        if (static_cast<int>(n) - r.nEctopic < p.minKeepBeats) {
            std::fill(r.ectopic.begin(), r.ectopic.end(), 0);
            r.nEctopic = r.nRhythm = r.nMorph = 0;
            r.morphRouteUsed = false;
            r.note = "would_leave_too_few";
            return r;
        }

        r.applied = true;
        r.note = haveRhythm ? "applied" : "applied_morphology_only";
        return r;
    }

} // namespace ectopic_mask