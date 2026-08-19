#pragma once
//
// beat_substitute.hpp
//
// Spec 4.6, "beat substitution, EWMA with alpha = 1/8".
//
// ---------------------------------------------------------------------------
// The spec conflates two different operations under one name, so they are
// separate here:
//
//   ewma_update()  -- a RUNNING TEMPLATE. avg = (1-a)*avg + a*beat, per sample
//                     index, so every column is an independent EWMA across beat
//                     index. This is what the spec's snippet actually computes.
//
//   fill_gap()     -- a REPLACEMENT BEAT for a rejected one, so downstream code
//                     that assumes an unbroken beat sequence still works.
//
// They must not be the same call. The spec's snippet takes `current` and mixes
// 12.5% of it into the average -- if `current` is the beat that was just
// rejected, the rejection has been undone: at alpha = 1/8 that beat still holds
// 34% of its original weight eight updates later and 12% after sixteen. So
// ewma_update() is only ever called on ACCEPTED beats, and a rejected slot is
// filled from its neighbours by fill_gap(), which never reads the bad beat.
//
// What alpha = 1/8 buys, in beats:
//   newest beat weight      12.5%      weight k beats back = a*(1-a)^k
//   last 8 beats            65.6%
//   last 16 beats           88.2%
//   half-life               5.2 beats
//   mean age of the data    7.0 beats
//   ~equivalent to a        15-point simple moving average
// The 1/8 is chosen because x + (y-x)/8 is a shift, not a divide.
//
// NOTE this is a MEAN, streaming and recency-weighted, whereas
// CreateEcgTemplates builds its templates from a per-sample MEDIAN over all kept
// beats -- batch and outlier-resistant. Two different estimators, not variations
// of one. Which feeds the dynamic master template selection of spec 9.5 is still
// unresolved.
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace beat_substitute {

    constexpr double kDefaultAlpha = 0.125;   // 1/8, per spec 4.6

    struct RunningTemplate {
        std::vector<double> avg;
        int    n_updates = 0;          // accepted beats folded in so far
        double alpha = kDefaultAlpha;
        bool   seeded = false;

        // Until ~1/alpha beats are in, the average still carries a visible share
        // of whatever seeded it.
        bool warm() const {
            return seeded && n_updates >= static_cast<int>(std::lround(1.0 / alpha));
        }
    };

    // Seed. The spec gives no initialisation and it needs one: starting from
    // zeros makes the template climb toward the real morphology over ~8 beats,
    // so everything before that is wrong. Seed with the first accepted beat.
    inline void seed(RunningTemplate& rt, const std::vector<double>& first_accepted,
        double alpha = kDefaultAlpha)
    {
        rt.avg = first_accepted;
        rt.alpha = alpha;
        rt.n_updates = first_accepted.empty() ? 0 : 1;
        rt.seeded = !first_accepted.empty();
    }

    // The spec's arithmetic. Call ONLY with an accepted beat.
    //
    // Length mismatch is tolerated (beats can differ by a sample after
    // alignment): the overlap updates, any tail of avg is left alone, a longer
    // `current` is ignored past the end. NaN samples in `current` are skipped so
    // one bad sample can't poison a column.
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

    // Free-function form matching the spec's signature, for callers keeping their
    // own state. Same caveat: accepted beats only.
    inline std::vector<double> ewma_update(const std::vector<double>& avg_old,
        const std::vector<double>& current,
        double alpha = kDefaultAlpha)
    {
        RunningTemplate rt;
        rt.avg = avg_old; rt.alpha = alpha; rt.seeded = !avg_old.empty();
        ewma_update(rt, current);
        return rt.avg;
    }

    // ---- Gap filling ------------------------------------------------------
    // A rejected beat leaves a hole. This produces a stand-in so the sequence
    // stays evenly populated -- "a smooth blend, not a copy", per the acceptance
    // test -- without letting the rejected samples back in.
    //
    // Blends the nearest accepted beat on each side, weighted by distance in
    // beat index: a hole one beat from a good beat on the left and five from one
    // on the right leans 5:1 toward the left. With only one side available it
    // blends that neighbour with the running template; with neither it returns
    // the running template alone.
    //
    // Never reads beats[t].
    inline std::vector<double> fill_gap(const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        int t,
        const RunningTemplate& rt)
    {
        const int n = static_cast<int>(beats.size());
        if (t < 0 || t >= n) return rt.avg;

        int L = -1, R = -1;
        for (int i = t - 1; i >= 0; --i)
            if (i < static_cast<int>(accepted.size()) && accepted[i] && !beats[i].empty()) { L = i; break; }
        for (int i = t + 1; i < n; ++i)
            if (i < static_cast<int>(accepted.size()) && accepted[i] && !beats[i].empty()) { R = i; break; }

        const std::vector<double>* a = nullptr;
        const std::vector<double>* b = nullptr;
        double wa = 0.0, wb = 0.0;

        if (L >= 0 && R >= 0) {
            const double dl = static_cast<double>(t - L), dr = static_cast<double>(R - t);
            a = &beats[L]; b = &beats[R];
            wa = dr / (dl + dr); wb = dl / (dl + dr);       // nearer neighbour weighs more
        }
        else if (L >= 0 && rt.seeded) { a = &beats[L]; b = &rt.avg; wa = 0.5; wb = 0.5; }
        else if (R >= 0 && rt.seeded) { a = &beats[R]; b = &rt.avg; wa = 0.5; wb = 0.5; }
        else if (L >= 0) return beats[L];
        else if (R >= 0) return beats[R];
        else return rt.avg;

        const size_t len = std::max(a->size(), b->size());
        std::vector<double> out(len, std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < len; ++j) {
            const bool ha = j < a->size() && !std::isnan((*a)[j]);
            const bool hb = j < b->size() && !std::isnan((*b)[j]);
            if (ha && hb) out[j] = wa * (*a)[j] + wb * (*b)[j];
            else if (ha)  out[j] = (*a)[j];
            else if (hb)  out[j] = (*b)[j];
        }
        return out;
    }

    // Whole-record pass, in beat order: fold accepted beats into the running
    // template, fill the rejected slots from their neighbours. Returns the
    // substituted series; `rt` is left holding the final template.
    //
    // Two-pass, deliberately: the running template is built from accepted beats
    // FIRST, so a gap near the start is filled against a warm template rather
    // than a cold one.
    inline std::vector<std::vector<double>> substitute_all(
        const std::vector<std::vector<double>>& beats,
        const std::vector<char>& accepted,
        RunningTemplate& rt,
        double alpha = kDefaultAlpha)
    {
        rt = RunningTemplate{};
        rt.alpha = alpha;

        for (size_t i = 0; i < beats.size(); ++i)
            if (i < accepted.size() && accepted[i]) ewma_update(rt, beats[i]);

        std::vector<std::vector<double>> out = beats;
        for (size_t i = 0; i < beats.size(); ++i) {
            const bool ok = (i < accepted.size() && accepted[i]);
            if (!ok) out[i] = fill_gap(beats, accepted, static_cast<int>(i), rt);
        }
        return out;
    }

    // How far a substituted beat sits from the original -- the acceptance test's
    // "a smooth blend, not a copy". 0 means it IS a copy.
    inline double substitution_distance(const std::vector<double>& original,
        const std::vector<double>& substituted)
    {
        const size_t n = std::min(original.size(), substituted.size());
        if (n == 0) return std::numeric_limits<double>::quiet_NaN();
        double s = 0.0; size_t m = 0;
        for (size_t j = 0; j < n; ++j) {
            if (std::isnan(original[j]) || std::isnan(substituted[j])) continue;
            const double d = original[j] - substituted[j];
            s += d * d; ++m;
        }
        return (m > 0) ? std::sqrt(s / m) : std::numeric_limits<double>::quiet_NaN();
    }

}   // namespace beat_substitute