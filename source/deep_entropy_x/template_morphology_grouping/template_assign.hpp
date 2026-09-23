#pragma once
/**
 * @file   template_assign.hpp
 * @brief  Assignment, spawning, merge/cap-raise, and pass-2 refinement for the
 *         Section 4.6 template bank. State lives in template_bank.hpp.
 *
 *         TWO PASSES, AND WHY. Pass 1 is causal: a beat is scored against the
 *         templates that exist AT THAT MOMENT, because that is the only thing
 *         "open a new template" and "subtype index in order of first
 *         appearance" can mean. That makes the result depend on beat order --
 *         the first exemplar of a morphology becomes its seed, so a noisy
 *         first exemplar is a poor centroid and later clean beats of the same
 *         morphology may fail the floor against it and spawn a redundant
 *         template. Pass 2 removes the order dependence by rescoring every
 *         beat against the FINAL templates and reassigning. Nothing physically
 *         prevented this: bins are processed retrospectively, so all the
 *         beats were always in hand.
 *
 *         WHAT PASS 2 MUST NOT DO. It must not renumber. spawn_seq is assigned
 *         in pass 1 and carried through unchanged, and subtype indices are
 *         issued once and immutable, so reassignment cannot cause PVC-2 to
 *         become PVC-1 between runs. It also does not spawn: pass 1 decided
 *         how many morphologies the bin contains, and letting pass 2 add more
 *         would make the passes disagree about that count with no rule for
 *         which wins.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <vector>

#include "template_bank.hpp"

namespace tbank {

    // ---------------------------------------------------------------------
    // Correlation
    // ---------------------------------------------------------------------

    struct CorrResult {
        double r = std::numeric_limits<double>::quiet_NaN();
        int    n_overlap = 0;
        bool   scorable() const { return n_overlap >= kMinOverlapColumns && !std::isnan(r); }
    };

    // Pearson r over the columns where both vectors are non-NaN.
    //
    // Returns NaN rather than 0.0 when unscorable. alignment.hpp's pearson()
    // returns 0.0 for too-few-samples and for zero variance, which collides
    // with a genuine correlation of zero -- harmless there, because 0.0 sits
    // below its corr_min either way, but wrong here: 0.0 is below every match
    // floor, so an unscorable beat would silently SPAWN A TEMPLATE. NaN keeps
    // "no match" and "cannot be scored" distinct.
    //
    // Note that r is invariant to scale and offset, so two morphologies
    // differing mainly in amplitude score near 1.0 and land on the same
    // template. Every beat here has been PQ-leveled by alignment Pass 3, so
    // amplitude differences are real signal rather than baseline drift, and
    // alignment itself does not trust r alone -- it pairs corr_min = 0.30 with
    // an RMS-deviation-in-SDs criterion for exactly this reason. Adding a
    // second criterion here would be a second assignment threshold, which the
    // spec forbids, so this is left as r alone. If amplitude-only variants
    // turn out to merge in practice, that is a spec question, not a local fix.
    inline CorrResult correlate(const std::vector<double>& a,
        const std::vector<double>& b, int lo = 0, int hi = std::numeric_limits<int>::max())
    {
        CorrResult out;
        double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
        int n = 0;
        const size_t w = std::min(a.size(), b.size());
        const size_t kLo = (lo > 0) ? static_cast<size_t>(lo) : 0;
        const size_t kHi = (hi < 0) ? 0
            : std::min(w, static_cast<size_t>(hi) + 1);
        for (size_t k = kLo; k < kHi; ++k) {
            if (std::isnan(a[k]) || std::isnan(b[k])) continue;
            sa += a[k]; sb += b[k];
            saa += a[k] * a[k]; sbb += b[k] * b[k]; sab += a[k] * b[k];
            ++n;
        }
        out.n_overlap = n;
        if (n < kMinOverlapColumns) return out;

        const double ma = sa / n, mb = sb / n;
        const double cov = sab / n - ma * mb;
        const double va = saa / n - ma * ma, vb = sbb / n - mb * mb;
        if (va <= 0.0 || vb <= 0.0) return out;   // flat vector: r undefined
        out.r = cov / std::sqrt(va * vb);
        return out;
    }

    // ---------------------------------------------------------------------
    // The band-match score (Section 4.6's metric)
    // ---------------------------------------------------------------------

    struct BandResult {
        // THE CORRELATION, in [-1, 1]. This is what every floor is compared
        // against; see bandMatch below for why it stopped being a fraction.
        double score = std::numeric_limits<double>::quiet_NaN();
        int    n_overlap = 0;
        double r = std::numeric_limits<double>::quiet_NaN();      // == score
        // Fraction of comparable columns inside the corridor: what `score`
        // used to be. Kept as a DIAGNOSTIC and filled only for the winning
        // template, by assignBeat. Nothing routes on it.
        double frac_in_band = std::numeric_limits<double>::quiet_NaN();
        bool   scorable() const {
            return n_overlap >= kMinOverlapColumns && !std::isnan(score);
        }
    };

    // The old score, kept as a report. A beat well inside its group's corridor
    // and a beat that correlates well but sits at twice the amplitude are
    // different situations, and after the change below only this number can
    // tell them apart -- so it is still computed, once per beat instead of once
    // per beat per template.
    inline double fracInBand(const std::vector<double>& beat,
        const BankTemplate& t)
    {
        if (t.band_lo.empty() || t.band_hi.empty())
            return std::numeric_limits<double>::quiet_NaN();
        const size_t w = std::min(beat.size(),
            std::min(t.band_lo.size(), t.band_hi.size()));
        int n = 0, inside = 0;
        for (size_t k = 0; k < w; ++k) {
            const double lo = t.band_lo[k], hi = t.band_hi[k];
            if (std::isnan(beat[k]) || std::isnan(lo) || std::isnan(hi)) continue;
            ++n;
            if (beat[k] >= lo && beat[k] <= hi) ++inside;
        }
        if (n < kMinOverlapColumns)
            return std::numeric_limits<double>::quiet_NaN();
        return static_cast<double>(inside) / static_cast<double>(n);
    }

    // Fraction of the beat's comparable samples lying inside the template's
    // per-column corridor. This is the quantity Section 4.6 thresholds at 0.85
    // (ECG) and 0.80 (PPG), read as a fraction rather than a percentage --
    // design note 1 in template_bank.hpp.
    //
    // Columns where either side is NaN are not comparable and are excluded from
    // both numerator and denominator. Excluding them from the denominator too is
    // the part that matters: counting an incomparable column as a miss would
    // score a short beat down for being short, which is a length judgment
    // wearing a morphology judgment's clothes.
    //
    // A template with no corridor at all (no members, nothing inherited) is
    // unscorable rather than scored 0. kUnscorable and "matched nothing" have to
    // stay distinct or an incomparable beat spawns a template, which is the same
    // trap correlate() documents.
    // ---- THE SCORE IS A CORRELATION AGAINST THE TEMPLATE -----------------
    //
    // "Morphology correlation thresholds for exclusion: r < 0.85 (ECG), r <
    // 0.80 (PPG)" -- the spec's own words, and now the metric.
    //
    // IT USED TO BE THE FRACTION OF COLUMNS INSIDE THE CORRIDOR, and that
    // fraction cannot distinguish two morphologies. A beat frame is ~1.8 RR
    // wide -- on a 1000 Hz record, over seven thousand columns -- and a QRS
    // occupies eighty of them. Everything else is baseline, T wave and
    // diastole, where a sinus beat and a ventricular beat sit in the same place
    // and both fall inside the corridor. An ectopic beat missing the corridor
    // on EVERY QRS column still scores 7080/7160 = 0.989, so it could not fail
    // an 0.85 floor -- and raising the floor did not help: it still could not
    // fail, while ordinary sinus beats began to. On a real record this showed
    // up as ectopic beats assigned to the sinus template and then removed
    // afterwards for prematurity: separated by the RHYTHM, which the bank is
    // supposed to be blind to, instead of by their shape.
    //
    // Correlation measures deviations about the mean, so a long flat stretch
    // contributes almost no variance and the QRS excursion dominates on its
    // own -- with no threshold on which columns count. Measured on a synthetic
    // record holding two morphologies: r = 0.14 between them, 0.996 and 0.998
    // within each.
    //
    // WHAT THIS GIVES UP. r is invariant to scale and offset, so a beat of the
    // same shape at twice the amplitude scores 1.0 and joins, where the
    // corridor would have rejected it. A real loss, and a deliberate one: the
    // alternative is a second assignment criterion, which the spec forbids in
    // as many words -- "Do not introduce a second, looser assignment
    // threshold". The amplitude Tukey fence runs after the partition and
    // catches a grossly mis-scaled beat as an outlier inside its own group.
    //
    // THE CORRIDOR STILL MATTERS, just not here. It is what the seeded slot 0
    // carries over from Phase 1, what young groups inherit, what is drawn as
    // the band under a template, and what fracInBand reports. It no longer
    // decides membership.
    //
    // AND IT NO LONGER GATES SCORABILITY. The old version returned unscorable
    // when a template had no corridor -- the state of every one-member group --
    // so a freshly spawned group could not accept its own second beat and
    // spawned again. Correlation needs only the waveform.
    inline BandResult bandMatch(const std::vector<double>& beat,
        const BankTemplate& t)
    {
        BandResult out;
        if (t.tmpl.empty()) return out;

        // Through the single correlate() every other consumer already uses --
        // cleanGroups, NSVT's cross-bin identity, the substitution band. A
        // second implementation here would be a second definition of "the same
        // shape", and they would disagree first on the beats that matter.
        // MORPHOLOGY SPLIT WINDOW. When the template carries a corr_halfwin
        // (samples), score only over [r_col - halfwin, r_col + halfwin] -- i.e.
        // +-0.5 s around the R peak (ECG) or the systolic peak (PPG). Only this
        // split metric is windowed; correlate's other callers (cleanGroups,
        // NSVT, substitution) pass no window and stay full-width. halfwin <= 0
        // or r_col < 0 means no restriction (reloaded banks, older callers).
        CorrResult cr;
        if (t.corr_halfwin > 0 && t.r_col >= 0) {
            cr = correlate(beat, t.tmpl,
                t.r_col - t.corr_halfwin, t.r_col + t.corr_halfwin);
        }
        else {
            cr = correlate(beat, t.tmpl);
        }
        out.n_overlap = cr.n_overlap;
        out.score = cr.r;
        out.r = cr.r;
        // frac_in_band stays NaN: it costs another pass over the overlap for
        // every beat against every template, and only the winner's value is
        // ever read. assignBeat fills it there.
        return out;
    }

    // ---------------------------------------------------------------------
    // Template recomputation
    // ---------------------------------------------------------------------

    // Column-wise NaN-skipping median over the template's members, matching
    // create_ecg_templates.hpp's medianOver(). Deliberately NOT the alpha =
    // 1/8 EWMA from the beat substitution section: an EWMA is a recursion
    // whose value depends on arrival order, and pass 2 exists to remove order
    // dependence. Median-over-members also keeps the archive property that a
    // template is exactly reconstructible from the per-beat flags.
    // `floor_corridor` supplies slot 0's per-column half-width, inherited when
    // this template has too few members for its own 2.5/97.5 corridor to be an
    // estimate of anything (design note 1). Pass nullptr for slot 0 itself, or
    // for a template known to have enough members.
    inline void recomputeTemplate(BankTemplate& t,
        const std::vector<std::vector<double>>& beats,
        int width,
        const std::vector<double>* floor_corridor = nullptr)
    {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        t.tmpl.assign(width, NaN);
        t.tmpl_iqr.assign(width, NaN);
        t.band_lo.assign(width, NaN);
        t.band_hi.assign(width, NaN);
        t.corridor_inherited = false;
        if (t.members.empty() || width <= 0) return;

        // Own corridor only when there are enough members for percentiles to
        // mean something; otherwise widths come from the floor below.
        const bool own_corridor = t.memberCount() >= kMinMembersForCorridor;
        t.corridor_inherited = !own_corridor;

        std::vector<double> col;
        col.reserve(t.members.size());
        for (int c = 0; c < width; ++c) {
            col.clear();
            for (uint32_t m : t.members) {
                if (m >= beats.size()) continue;
                const auto& b = beats[m];
                if (c < static_cast<int>(b.size()) && !std::isnan(b[c]))
                    col.push_back(b[c]);
            }
            if (col.empty()) continue;
            const size_t n = col.size();

            // nth_element, not sort. Only three order statistics are needed per
            // column (median, Q1, Q3), and this runs once per column per
            // template per recompute -- on a 991-member slot 0 across a 200
            // column axis that is the single hottest loop in the pass. Partial
            // selection is O(n) against sort's O(n log n) and measured roughly
            // 3x faster here at these sizes.
            //
            // The nth_element calls are ordered low-to-high so each one only
            // has to partition the range the previous one left, rather than the
            // whole column again.
            const size_t iq1 = n / 4;
            const size_t imid = n / 2;
            const size_t iq3 = std::min(n - 1, (3 * n) / 4);

            std::nth_element(col.begin(), col.begin() + iq1, col.end());
            const double q1 = col[iq1];
            std::nth_element(col.begin() + iq1, col.begin() + imid, col.end());
            const double hi_mid = col[imid];
            std::nth_element(col.begin() + imid, col.begin() + iq3, col.end());
            const double q3 = col[iq3];

            if (n % 2) {
                t.tmpl[c] = hi_mid;
            }
            else {
                // Even n: the median averages the two central values, and the
                // lower one is the max of everything below imid -- already
                // partitioned there by the nth_element calls above, so no
                // further selection is needed.
                const double lo_mid =
                    *std::max_element(col.begin(), col.begin() + imid);
                t.tmpl[c] = 0.5 * (lo_mid + hi_mid);
            }
            // IQR as the spread measure, consistent with the *_iqr fields
            // already carried alongside every template in TemplateBin.
            t.tmpl_iqr[c] = q3 - q1;

            // --- the 2.5/97.5 corridor -------------------------------------
            if (own_corridor) {
                // A COLUMN TOO THIN FOR ITS OWN PERCENTILES FALLS BACK; IT
                // DOES NOT GO SILENT.
                //
                // own_corridor is a MEMBER count. n is this COLUMN's finite
                // count, and those are different numbers: beats are NaN-padded,
                // so a 900-member template still has edge columns reached by two
                // or three of them.
                //
                // The four order statistics below run in ascending order so each
                // partitions only what the previous left, which requires
                //     i_lo <= i_lo1 <= i_hi <= i_hi1
                // and at n == 2 that inverts -- i_lo1 is 1 while i_hi is 0, so
                // the third nth_element gets first > nth. Debug asserts "vector
                // iterator range transposed"; Release is undefined and the
                // introselect loop need not terminate. It presented as a hang
                // partway through the per-bin split.
                //
                // SKIPPING THE COLUMN IS NOT THE FIX, and was tried: band_lo and
                // band_hi are initialised to NaN above, and the inheritance block
                // at the bottom of this function only runs for a young TEMPLATE.
                // So a skipped column inside a mature template keeps NaN,
                // bandMatch cannot score it, and beats stop matching on the
                // padded edges -- they spawn instead of joining, and the cohorts
                // fragment into slots too thin to draw. That is what emptied the
                // pulse panels.
                //
                // So the column takes the same inherited half-width a young
                // template takes: slot 0's spread where there is one, a fraction
                // of this template's own peak-to-peak otherwise.
                // corridorInflation gets the COLUMN's n, because it is this
                // column's centre that was estimated from n values.
                if (n < static_cast<size_t>(kMinMembersForCorridor)) {
                    double half = std::numeric_limits<double>::quiet_NaN();
                    if (floor_corridor && c < static_cast<int>(floor_corridor->size()))
                        half = (*floor_corridor)[c];
                    if (std::isnan(half) || half <= 0.0) {
                        double plo = std::numeric_limits<double>::infinity();
                        double phi = -std::numeric_limits<double>::infinity();
                        for (int k = 0; k < width; ++k) {
                            if (std::isnan(t.tmpl[k])) continue;
                            plo = std::min(plo, t.tmpl[k]);
                            phi = std::max(phi, t.tmpl[k]);
                        }
                        const double ptp = (phi > plo) ? (phi - plo) : 0.0;
                        half = kFallbackCorridorFrac * ptp;
                    }
                    if (half > 0.0) {
                        half *= corridorInflation(static_cast<int>(n));
                        t.band_lo[c] = t.tmpl[c] - half;
                        t.band_hi[c] = t.tmpl[c] + half;
                    }
                    continue;
                }

                // Percentiles by position with linear interpolation, over the
                // column's values. col is only partially ordered by the
                // nth_element calls above, so this sorts -- it runs only for
                // templates that have earned their own corridor, and slot 0
                // dominates the cost either way.
                // nth_element, NOT sort -- for the reason the median block
                // above states and which the first version of this block
                // ignored. A std::sort here runs per column, per template, per
                // recompute, inside the hottest loop in the pass. Four extra
                // order statistics, ordered ascending so each partitions only
                // what the previous left.
                const double x_lo = 0.025 * (static_cast<double>(n) - 1.0);
                const double x_hi = 0.975 * (static_cast<double>(n) - 1.0);
                const size_t i_lo = static_cast<size_t>(std::floor(x_lo));
                const size_t i_lo1 = std::min(n - 1, i_lo + 1);
                const size_t i_hi = static_cast<size_t>(std::floor(x_hi));
                const size_t i_hi1 = std::min(n - 1, i_hi + 1);
                const double f_lo = x_lo - static_cast<double>(i_lo);
                const double f_hi = x_hi - static_cast<double>(i_hi);
                std::nth_element(col.begin(), col.begin() + i_lo, col.end());
                const double p2a = col[i_lo];
                std::nth_element(col.begin() + i_lo, col.begin() + i_lo1, col.end());
                const double p2b = col[i_lo1];
                std::nth_element(col.begin() + i_lo1, col.begin() + i_hi, col.end());
                const double p97a = col[i_hi];
                std::nth_element(col.begin() + i_hi, col.begin() + i_hi1, col.end());
                const double p97b = col[i_hi1];
                const double p2 = p2a * (1.0 - f_lo) + p2b * f_lo;
                const double p97 = p97a * (1.0 - f_hi) + p97b * f_hi;
                const double mid = 0.5 * (p2 + p97);
                const double half = 0.5 * (p97 - p2)
                    * corridorInflation(t.memberCount());
                t.band_lo[c] = mid - half;
                t.band_hi[c] = mid + half;
            }
        }

        // --- inherited corridor ---------------------------------------------
        // Centred on THIS template's own median, widened by slot 0's spread.
        // Centring on the inheritor and not the donor is the whole point: the
        // shape being scored against is this template's, only the tolerance is
        // borrowed. Centring on slot 0 would make every young template a
        // restatement of sinus and no beat would ever fail to match it.
        if (!own_corridor) {
            for (int c = 0; c < width; ++c) {
                if (std::isnan(t.tmpl[c])) continue;
                double half = std::numeric_limits<double>::quiet_NaN();
                if (floor_corridor && c < static_cast<int>(floor_corridor->size()))
                    half = (*floor_corridor)[c];
                if (std::isnan(half) || half <= 0.0) {
                    // Neither this template nor slot 0 has a spread here. Fall
                    // back to a fraction of the template's own amplitude, which
                    // is unitless and survives normalization.
                    double lo = std::numeric_limits<double>::infinity();
                    double hi = -std::numeric_limits<double>::infinity();
                    for (int k = 0; k < width; ++k) {
                        if (std::isnan(t.tmpl[k])) continue;
                        lo = std::min(lo, t.tmpl[k]);
                        hi = std::max(hi, t.tmpl[k]);
                    }
                    const double ptp = (hi > lo) ? (hi - lo) : 0.0;
                    half = kFallbackCorridorFrac * ptp;
                    if (half <= 0.0) continue;   // flat template: leave NaN
                }
                half *= corridorInflation(t.memberCount());
                t.band_lo[c] = t.tmpl[c] - half;
                t.band_hi[c] = t.tmpl[c] + half;
            }
        }
    }

    // ---------------------------------------------------------------------
    // (THE PER-CHANNEL ASSIGNMENT ENGINE WAS HERE AND IS GONE:
    //  seedCorridorHalfWidth, mergeTemplates, AssignOutcome, assignBeat,
    //  recomputeAll(TemplateBank) and refinePass -- pass 1 and pass 2 over one
    //  channel's beats.
    //
    //  Its driver was bin_pipeline::runChannel, which partitioned each channel
    //  independently; that died when the partition became joint, and these had
    //  had no caller since. They are not the same code as jbank's assignSlice /
    //  mergeGroups / recomputeAll: those key on a SLICE, so a beat is one beat
    //  across all four channels, which is the whole point and is exactly what
    //  this file could not express.
    //
    //  WHAT SURVIVES ABOVE IS THE SHARED MEASUREMENT: correlate, fracInBand,
    //  bandMatch and recomputeTemplate operate on one channel's array and know
    //  nothing about how beats were grouped, so jbank scores through them and
    //  there is one definition of "how well does this beat match" rather than
    //  two that could drift.)
    // ---------------------------------------------------------------------

}  // namespace tbank