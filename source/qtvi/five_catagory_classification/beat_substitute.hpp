#pragma once
//
// beat_substitute.hpp
//
// Spec 4.6, "beat substitution, EWMA with alpha = 1/8".
//
// STATUS: stubs, except ewma_update() which is the spec's own arithmetic and
// is short enough to be complete.
//
// ---------------------------------------------------------------------------
// The spec conflates two different operations under one name, so they are
// separate functions here:
//
//   ewma_update()  -- a RUNNING TEMPLATE. avg = (1-a)*avg + a*beat, applied
//                     per sample index, so every column of the template is an
//                     independent EWMA across beat index. This is what the
//                     spec's snippet actually computes.
//
//   fill_gap()     -- a REPLACEMENT BEAT for a rejected one, so downstream
//                     code that assumes an unbroken beat sequence still works.
//
// They must not be the same call. The spec's snippet takes `current` and mixes
// 12.5% of it into the average -- if `current` is the beat that was just
// rejected, the rejection has been undone: at alpha = 1/8 that beat still
// retains 34% of its original weight eight updates later, and 12% after
// sixteen. ewma_update() is therefore only ever called on ACCEPTED beats;
// a rejected beat's slot is filled from its neighbours by fill_gap().
//
// What alpha = 1/8 buys, in beats:
//   newest beat weight      12.5%      weight k beats back = a*(1-a)^k
//   last 8 beats            65.6%
//   last 16 beats           88.2%
//   half-life               5.2 beats
//   mean age of the data    7.0 beats
//   ~equivalent to a        15-point simple moving average
// The 1/8 is chosen because (x + (y-x)/8) is a shift, not a divide.
//
// NOTE this is a MEAN, streaming and recency-weighted, whereas
// CreateEcgTemplates builds its templates from a per-sample MEDIAN over all
// kept beats -- batch and outlier-resistant. They are two different estimators
// with different behaviour, not variations of one. Which of them feeds the
// dynamic master template selection of spec 9.5 is still unresolved.
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace beat_substitute {

    constexpr double kDefaultAlpha = 0.125;   // 1/8, per spec 4.6

    // Diagnostics for one EWMA state, so a caller can tell a warmed-up
    // template from one that is still converging out of its seed.
    struct RunningTemplate {
        std::vector<double> avg;          // the template itself
        int    n_updates = 0;             // accepted beats folded in so far
        double alpha = kDefaultAlpha;
        bool   seeded = false;            // false => avg is not yet meaningful

        // Warm-up guard: until roughly 1/alpha beats have been folded in, the
        // average still carries a visible share of whatever seeded it.
        bool warm() const { return seeded && n_updates >= static_cast<int>(std::lround(1.0 / alpha)); }
    };

    // Seed. The spec gives no initialisation, and it needs one: starting from a
    // zero vector makes the template climb toward the real morphology over
    // ~8 beats, so everything before that is wrong. Seed with the first
    // accepted beat instead of zeros.
    inline void seed(RunningTemplate& rt, const std::vector<double>& first_accepted,
        double alpha = kDefaultAlpha)
    {
        rt.avg = first_accepted;
        rt.alpha = alpha;
        rt.n_updates = first_accepted.empty() ? 0 : 1;
        rt.seeded = !first_accepted.empty();
    }

    // The spec's arithmetic: out[j] = (1-a)*avgOld[j] + a*current[j].
    // Call ONLY with an accepted (category-1) beat -- see the header note.
    // Length mismatch is not an error here (beats can differ by a sample after
    // alignment): the overlap updates, any tail of avgOld is left alone, and a
    // longer `current` is ignored past the end. NaN samples in `current` are
    // skipped so one bad sample can't poison a column.
    inline void ewma_update(RunningTemplate& rt, const std::vector<double>& current)
    {
        if (current.empty()) return;
        if (!rt.seeded) { seed(rt, current, rt.alpha); return; }
        const size_t n = std::min(rt.avg.size(), current.size());
        const double a = rt.alpha;
        for (size_t j = 0; j < n; ++j) {
            if (std::isnan(current[j])) continue;
            if (std::isnan(rt.avg[j])) { rt.avg[j] = current[j]; continue; }
            rt.avg[j] += a * (current[j] - rt.avg[j]);   // == (1-a)*avg + a*current
        }
        ++rt.n_updates;
    }

    // Free-function form, matching the spec's signature, for callers that keep
    // their own state. Same caveat: accepted beats only.
    inline std::vector<double> ewma_update(const std::vector<double>& avg_old,
        const std::vector<double>& current,
        double alpha = kDefaultAlpha)
    {
        RunningTemplate rt; rt.avg = avg_old; rt.alpha = alpha; rt.seeded = !avg_old.empty();
        ewma_update(rt, current);
        return rt.avg;
    }

    // ---- Gap filling ------------------------------------------------------
    // A rejected beat leaves a hole. This produces a stand-in so the sequence
    // stays evenly populated -- "a smooth blend, not a copy", per the
    // acceptance test -- without letting the rejected samples back in.
    //
    // TODO: blend the nearest accepted beat on each side (weighted by distance
    //       in beat index), falling back to the running template when only one
    //       side exists, and to the running template alone at the record
    //       edges. Must NOT read the rejected beat.
    inline std::vector<double> fill_gap(const std::vector<std::vector<double>>& /*beats*/,
        const std::vector<char>& /*accepted*/,
        int /*t*/,
        const RunningTemplate& rt)
    {
        return rt.avg;   // placeholder: the running template
    }

    // Whole-record pass: walk the beats in order, fold accepted ones into the
    // running template, fill the rejected slots. Returns the substituted beat
    // series; `rt` is left holding the final template.
    // TODO
    inline std::vector<std::vector<double>> substitute_all(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& /*accepted*/,
        RunningTemplate& /*rt*/,
        double /*alpha*/ = kDefaultAlpha)
    {
        return beats;
    }

}   // namespace beat_substitute