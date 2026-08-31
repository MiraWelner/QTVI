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
        const std::vector<double>& b)
    {
        CorrResult out;
        double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
        int n = 0;
        const size_t w = std::min(a.size(), b.size());
        for (size_t k = 0; k < w; ++k) {
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
        double score = std::numeric_limits<double>::quiet_NaN();  // fraction in [0,1]
        int    n_overlap = 0;
        double r = std::numeric_limits<double>::quiet_NaN();      // reported only
        bool   scorable() const {
            return n_overlap >= kMinOverlapColumns && !std::isnan(score);
        }
    };

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
    inline BandResult bandMatch(const std::vector<double>& beat,
        const BankTemplate& t)
    {
        BandResult out;
        if (t.band_lo.empty() || t.band_hi.empty()) return out;

        const size_t w = std::min(beat.size(),
            std::min(t.band_lo.size(), t.band_hi.size()));
        int n = 0, inside = 0;
        for (size_t k = 0; k < w; ++k) {
            const double lo = t.band_lo[k], hi = t.band_hi[k];
            if (std::isnan(beat[k]) || std::isnan(lo) || std::isnan(hi)) continue;
            ++n;
            if (beat[k] >= lo && beat[k] <= hi) ++inside;
        }
        out.n_overlap = n;
        if (n < kMinOverlapColumns) return out;
        out.score = static_cast<double>(inside) / static_cast<double>(n);
        out.r = correlate(beat, t.tmpl).r;   // reported alongside, never routed on
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
                // Percentiles by position with linear interpolation, over the
                // column's values. col is only partially ordered by the
                // nth_element calls above, so this sorts -- it runs only for
                // templates that have earned their own corridor, and slot 0
                // dominates the cost either way.
                std::sort(col.begin(), col.end());
                auto pct = [&](double p) {
                    const double x = p * (static_cast<double>(n) - 1.0) / 100.0;
                    const size_t i0 = static_cast<size_t>(std::floor(x));
                    const size_t i1 = std::min(n - 1, i0 + 1);
                    const double f = x - static_cast<double>(i0);
                    return col[i0] * (1.0 - f) + col[i1] * f;
                    };
                const double mid = 0.5 * (pct(2.5) + pct(97.5));
                const double half = 0.5 * (pct(97.5) - pct(2.5))
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

    // Slot 0's per-column corridor half-width, for young templates to inherit.
    // Half the 2.5/97.5 width, so an inheritor gets a corridor of the same total
    // span as sinus has, centred on its own median.
    inline std::vector<double> seedCorridorHalfWidth(const TemplateBank& bank) {
        std::vector<double> out;
        if (bank.size() == 0) return out;
        const BankTemplate& s0 = bank.templates[0];
        out.assign(s0.band_lo.size(), std::numeric_limits<double>::quiet_NaN());
        for (size_t c = 0; c < s0.band_lo.size() && c < s0.band_hi.size(); ++c) {
            if (std::isnan(s0.band_lo[c]) || std::isnan(s0.band_hi[c])) continue;
            out[c] = 0.5 * (s0.band_hi[c] - s0.band_lo[c]);
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // Merge
    // ---------------------------------------------------------------------

    // Absorbs b into a and erases b. Callers must have established that a and
    // b are not both confirmed -- findMergePair() reports that case with
    // both_confirmed set instead of offering it as a merge.
    //
    // Label inheritance: unlabeled-into-confirmed keeps the confirmed label
    // and its already-issued subtype. This is not an inference from
    // morphological similarity -- the operator confirmed a beat that is now a
    // member of the merged template, so the confirmation still holds for the
    // beats it was made about. The lower spawn_seq survives so that "order of
    // first appearance" is unaffected by merge history.
    // `defer_recompute` matters more here than anywhere else. Absorbing a
    // 152-member template into a 996-member one means sorting 1148 values per
    // column across the whole axis, and a saturated bin logged 29 to 67 merges
    // -- so recomputing on every merge dominated the entire pass. Deferring is
    // safe for the same reason it is safe in assignBeat: recomputeAll() runs
    // between the passes, and pass 2 rescores everything against the result.
    //
    // The exception is a merge whose product must be scoreable immediately,
    // which does not arise during pass 1: the merged template is only ever
    // scored against on a LATER beat, by which point nothing has read its
    // stale samples. Callers outside a deferred pass should pass false.
    inline void mergeTemplates(TemplateBank& bank, int a, int b,
        const std::vector<std::vector<double>>& beats,
        int width,
        bool defer_recompute = true)
    {
        if (a == b || a < 0 || b < 0) return;
        if (a > b) std::swap(a, b);

        BankTemplate& ta = bank.templates[a];
        BankTemplate& tb = bank.templates[b];

        // NO LABEL INHERITANCE. This used to copy a confirmed label and its
        // subtype onto the merged product when one side was confirmed. Section
        // 4.6 forbids it: a label reaches every beat in its template, so the
        // absorbed template's beats would become PVC on the evidence of a
        // correlation, and "do not infer a label from morphological similarity
        // to a labeled template" is exactly that operation. The old defense here
        // was that the operator's confirmation still holds for the beats it was
        // made about -- true, and silent about the other half of the merged
        // membership, which was never confirmed at all.
        //
        // findMergePair() now blocks any pair containing a confirmed template,
        // so this branch is unreachable rather than merely unused. It is checked
        // anyway, because a caller reaching mergeTemplates() by another route
        // must not quietly relabel beats -- a wrong PVC label in the archive is
        // indistinguishable from a right one afterwards.
        if (ta.confirmed() || tb.confirmed()) return;

        ta.spawn_seq = std::min(ta.spawn_seq, tb.spawn_seq);

        ta.members.insert(ta.members.end(), tb.members.begin(), tb.members.end());
        std::sort(ta.members.begin(), ta.members.end());
        ta.members.erase(std::unique(ta.members.begin(), ta.members.end()),
            ta.members.end());

        bank.templates.erase(bank.templates.begin() + b);
        if (!defer_recompute)
            recomputeTemplate(bank.templates[a], beats, width);
    }

    // ---------------------------------------------------------------------
    // Pass 1: causal assignment
    // ---------------------------------------------------------------------

    struct AssignOutcome {
        int32_t template_id = kNoMatch;
        // The band-match score this assignment was decided on, and Pearson r
        // alongside it for reporting only. Both are kept because they disagree
        // informatively: r is blind to amplitude and the corridor is not, so
        // high r with low score is an amplitude outlier of a known shape.
        double  score = std::numeric_limits<double>::quiet_NaN();
        double  r = std::numeric_limits<double>::quiet_NaN();
        bool    spawned = false;
        bool    merged = false;
        bool    cap_raised = false;
        bool    unscorable = false;
    };

    // Scores `beat` against every template currently in the bank and assigns
    // it to the best match at or above `floor`. Below the floor against all of
    // them, the beat belongs to nothing present and seeds a new template.
    //
    // ONE THRESHOLD does both jobs. A looser assignment threshold would file
    // beats into templates they do not match, which is the variance inflation
    // this section exists to remove, hidden inside a cluster instead of
    // visible as a bad median.
    //
    // Best match means argmax, not first-above-floor: with several templates
    // in the bank a beat can clear the floor against more than one.
    // `defer_recompute` skips rebuilding the template after appending a member,
    // which is the difference between a fast pass and an unusable one.
    // recomputeTemplate() sorts every column over every member, so calling it
    // per assignment costs O(beats * members * width * log members) -- with ~850
    // members and a ~200 column axis that is on the order of 10^8 sorted
    // elements per channel per bin, times three channels, times every bin.
    //
    // Deferring costs almost nothing in accuracy because pass 2 rescores every
    // beat against the final templates anyway; incremental recomputation was
    // only ever refining the reference that pass 2 replaces. It also has a small
    // benefit: slot 0 stays the trimmed sinus median for the whole of pass 1
    // rather than drifting as beats accumulate, so every spawn decision in the
    // pass is judged against the same reference.
    //
    // A SPAWNED template is still built immediately -- it has exactly one
    // member, so the cost is trivial, and subsequent beats cannot be scored
    // against it otherwise.
    inline AssignOutcome assignBeat(TemplateBank& bank,
        uint32_t beat_idx,
        const std::vector<std::vector<double>>& beats,
        int width,
        bool is_ppg,
        uint64_t bin_index,
        int channel,
        std::vector<CapRaiseEvent>* events,
        BinCounts* counts,
        bool defer_recompute = true)
    {
        AssignOutcome out;
        if (beat_idx >= beats.size()) return out;
        const std::vector<double>& beat = beats[beat_idx];
        const double floor = bank.matchFloorFor(is_ppg);

        int    best = -1;
        double best_score = -std::numeric_limits<double>::infinity();
        double best_r = std::numeric_limits<double>::quiet_NaN();
        bool   any_scorable = false;

        for (int i = 0; i < bank.size(); ++i) {
            const BandResult br = bandMatch(beat, bank.templates[i]);
            if (!br.scorable()) continue;
            any_scorable = true;
            if (br.score > best_score) {
                best_score = br.score; best_r = br.r; best = i;
            }
        }

        // Unscorable against everything: too little axis overlap for r to mean
        // anything. Do NOT spawn -- a beat that cannot be compared is not
        // evidence of a new morphology.
        if (!any_scorable && bank.size() > 0) {
            out.unscorable = true;
            out.template_id = kUnscorable;
            if (counts) ++counts->n_unscorable;
            return out;
        }

        if (best >= 0 && best_score >= floor) {
            bank.templates[best].members.push_back(beat_idx);
            if (!defer_recompute) {
                const std::vector<double> fc = seedCorridorHalfWidth(bank);
                recomputeTemplate(bank.templates[best], beats, width,
                    best == 0 ? nullptr : &fc);
            }
            ++bank.assigned_beats;
            out.template_id = best;
            out.score = best_score;
            out.r = best_r;
            return out;
        }

        // --- spawn -------------------------------------------------------
        // At cap, free a slot by merging the two closest templates. When those
        // two are both confirmed, raise the cap for this bin and log instead:
        // merging them would not lose a beat, it would change a clinical
        // finding, silently, by collapsing polymorphic ectopy into
        // monomorphic. Since the bin verdict counts only confirmed templates,
        // that merge would destroy the sole evidence for polymorphy.
        if (bank.atCap()) {
            // Closeness in the SAME metric the section thresholds, so "the two
            // closest templates" and "this beat matches that template" are the
            // same kind of statement. Symmetrized by max because band-match is
            // directional: x's samples against y's corridor is not y's against
            // x's, and a wide template swallows a narrow one asymmetrically.
            const MergeCandidate mc = findMergePair(bank,
                [](const BankTemplate& x, const BankTemplate& y) {
                    const double a = bandMatch(x.tmpl, y).score;
                    const double b = bandMatch(y.tmpl, x).score;
                    if (std::isnan(a)) return b;
                    if (std::isnan(b)) return a;
                    return std::max(a, b);
                });

            if (mc.valid() && !mc.both_confirmed) {
                mergeTemplates(bank, mc.a, mc.b, beats, width);
                out.merged = true;
                if (counts) {
                    ++counts->n_merges;
                    if (mc.garbage_pair) ++counts->n_merges_garbage;
                    else                 ++counts->n_merges_real;
                }
            }
            else {
                const int32_t old_cap = bank.effective_cap;
                ++bank.effective_cap;
                out.cap_raised = true;
                if (counts) ++counts->n_cap_raises;
                if (events && mc.valid()) {
                    CapRaiseEvent e;
                    e.bin_index = bin_index;
                    e.channel = channel;
                    e.old_cap = old_cap;
                    e.new_cap = bank.effective_cap;
                    e.template_a = mc.a;
                    e.template_b = mc.b;
                    e.closeness = mc.closeness;
                    e.label_a = bank.templates[mc.a].label_code;
                    e.label_b = bank.templates[mc.b].label_code;
                    events->push_back(e);
                }
            }
        }

        BankTemplate t;
        t.spawn_seq = bank.next_spawn_seq++;
        t.members.push_back(beat_idx);
        // A one-member template has no corridor of its own, so it inherits slot
        // 0's spread. Without this the spec's metric cannot grow a template it
        // has just opened -- see design note 1.
        {
            const std::vector<double> fc = seedCorridorHalfWidth(bank);
            recomputeTemplate(t, beats, width, &fc);
        }
        bank.templates.push_back(std::move(t));
        ++bank.assigned_beats;

        out.template_id = bank.size() - 1;
        out.spawned = true;
        out.score = (best >= 0) ? best_score : std::numeric_limits<double>::quiet_NaN();
        out.r = best_r;
        if (counts) ++counts->n_spawns;
        return out;
    }

    // Rebuilds every template from its members. Call once after a deferred
    // pass 1, before pass 2 scores anything against them.
    // Slot 0 first, then everything else against slot 0's fresh corridor -- the
    // order matters, because an inheritor must not borrow a stale spread.
    inline void recomputeAll(TemplateBank& bank,
        const std::vector<std::vector<double>>& beats,
        int width)
    {
        if (bank.size() == 0) return;
        recomputeTemplate(bank.templates[0], beats, width);
        const std::vector<double> fc = seedCorridorHalfWidth(bank);
        for (int i = 1; i < bank.size(); ++i)
            recomputeTemplate(bank.templates[i], beats, width, &fc);
    }

    // ---------------------------------------------------------------------
    // Pass 2: refinement
    // ---------------------------------------------------------------------

    // Rescores every beat against the FINAL templates and reassigns to the
    // best match, removing pass 1's order dependence.
    //
    // Does not spawn: pass 1 decided how many morphologies the bin contains.
    // Does not renumber: spawn_seq and subtype are carried through untouched.
    // A beat that clears the floor against nothing keeps its pass-1 template
    // -- it already is the best available account of that beat, and dropping
    // it would silently shrink the archive's beat count.
    //
    // Returns the number of beats that changed template. Empty templates are
    // erased afterwards, which shifts indices, so `assignment` is rewritten
    // through a remap rather than left holding stale ids.
    inline uint32_t refinePass(TemplateBank& bank,
        const std::vector<std::vector<double>>& beats,
        std::vector<int32_t>& assignment,
        int width,
        bool is_ppg,
        BinCounts* counts)
    {
        if (bank.size() == 0) return 0;

        // With one template there is nowhere to reassign to: every beat's best
        // match is the template it is already in, so the whole scoring sweep and
        // the recompute that follows it are provably no-ops. On a clean record
        // that is most bins -- 8 of the first 20 on the first real record tested
        // -- so this early-out is not a micro-optimisation, it removes pass 2
        // entirely for the common case.
        if (bank.size() == 1) return 0;

        const double floor = bank.matchFloorFor(is_ppg);

        std::vector<std::vector<uint32_t>> new_members(bank.size());
        uint32_t moved = 0;

        for (uint32_t bi = 0; bi < assignment.size(); ++bi) {
            const int32_t cur = assignment[bi];
            if (cur == kUnscorable || cur < 0) continue;
            if (bi >= beats.size()) continue;

            int    best = -1;
            double best_score = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < bank.size(); ++i) {
                const BandResult br = bandMatch(beats[bi], bank.templates[i]);
                if (!br.scorable()) continue;
                if (br.score > best_score) { best_score = br.score; best = i; }
            }

            int32_t dest = cur;
            if (best >= 0 && best_score >= floor) dest = best;
            if (dest >= static_cast<int32_t>(new_members.size())) dest = cur;

            if (dest != cur) { ++moved; }
            new_members[dest].push_back(bi);
        }

        for (int i = 0; i < bank.size(); ++i)
            bank.templates[i].members = std::move(new_members[i]);
        recomputeAll(bank, beats, width);

        // Drop templates that lost every member. A confirmed template cannot
        // be dropped: the operator's judgment about those beats stands even if
        // refinement moved them, and erasing it would retire a subtype index
        // that the archive already refers to.
        std::vector<int32_t> remap(bank.size(), -1);
        std::vector<BankTemplate> kept;
        kept.reserve(bank.templates.size());
        for (int i = 0; i < bank.size(); ++i) {
            if (bank.templates[i].members.empty() && !bank.templates[i].confirmed())
                continue;
            remap[i] = static_cast<int32_t>(kept.size());
            kept.push_back(std::move(bank.templates[i]));
        }
        bank.templates = std::move(kept);

        uint32_t total = 0;
        for (const auto& t : bank.templates) total += t.memberCount();
        bank.assigned_beats = total;

        for (auto& a : assignment) {
            if (a < 0) continue;
            a = (a < static_cast<int32_t>(remap.size())) ? remap[a] : -1;
        }
        // Reassignment above moved beats between templates; rebuild the
        // per-beat ids from the surviving membership so the two cannot drift.
        for (int i = 0; i < bank.size(); ++i)
            for (uint32_t m : bank.templates[i].members)
                if (m < assignment.size()) assignment[m] = i;

        if (counts) counts->n_reassigned_pass2 += moved;
        return moved;
    }

}  // namespace tbank
