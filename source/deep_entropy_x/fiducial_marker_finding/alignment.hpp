#pragma once
//
// alignment.hpp
//
// Per-bin beat alignment for ECG and PPG. Runs BEFORE any normalization
// or template averaging.
//
//   1) Slice one beat per RR window as 0.5*RR before R to 1.3*RR after R
//      (integer-truncated) -- percent_interval_preceeding_rpeak /
//      percent_interval_following_rpeak below, which are the only definition.
//      This line said 0.4 and bin_plot_widget's copy said 0.3; neither was
//      ever the value, and the pulse time axis was built on the guess.
//   2) Tukey outlier rejection: ECG is based on RR length and distance from max to min, PPG is based on 50% upslope location and distance from max to min.
//   3) Horizontal align the ECG such that R peaks are aligned, if it is the Q peak screening then after align by Q peak
//   4) Vertical DC shift:
//        - ECG: match each beat's PR-baseline mean to the reference beat's.
//        - PPG: match each beat's foot-baseline mean (a small window around
//               its own foot column) to the reference beat's.
//

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <functional>
#include "fiducial_marker_finding/feature_marks.hpp"


namespace alignment {

    // HOW FAR A BEAT'S VERTICAL LEVELING GOT, not which segment was chosen.
    //
    // Both passes are applied when both are usable. Pass 2 levels on TP (this
    // beat's T-end to the NEXT beat's P-onset) as a coarse move; pass 3 then
    // re-measures PQ (P-end to Q-onset) on the already-shifted beat and
    // finalizes the zero there. PQ is the authoritative reference and runs
    // last, so a beat that reaches it sits on PQ.
    //
    //   PQ   -- the PQ finalize ran. TP may have run first; this does NOT mean
    //           PQ was used instead of TP.
    //   TP   -- only TP ran; PQ's landmarks were not usable.
    //   NONE -- neither was usable. The beat is left un-shifted, and
    //           create_ecg_templates.hpp excludes it from the per-sample
    //           amplitude aggregation (median/std) rather than let an
    //           unreliable, unadjusted amplitude contribute.
    enum class BaselineSource { TP, PQ, NONE };
    constexpr double percent_interval_preceeding_rpeak = 0.5; //how far before the R peak the snip goes, in terms of percent of the RR interval length
    constexpr double percent_interval_following_rpeak = 1.4;   //how far after the R peak the snip goes, in terms of percent of the RR interval length

    // Sample counts for a given RR (integer-truncated).
    inline int64_t rr_before_samples(int64_t rr) {
        return static_cast<int64_t>(percent_interval_preceeding_rpeak * rr);
    }
    inline int64_t rr_after_samples(int64_t rr) {
        return static_cast<int64_t>(percent_interval_following_rpeak * rr);
    }

    struct ecg_beat_set {
        //a set of individually-aligned beats that will become a template.
        std::vector<std::vector<double>> beats;
        std::vector<size_t> r_indices;
        std::vector<int>    rr_lens;
        int    median_length = -1;
        size_t total_beats = 0;
        int r_aligned_col = -1;
        int q_aligned_col = -1;
        int ref_beat_index = -1;
        // Parallel to `beats`: which segment supplied each beat's DC
        // baseline (TP/PQ/NONE).
        std::vector<BaselineSource> baseline_source;
        // Parallel to `beats`: the Section 4.6 rhythm verdict, assigned after
        // slicing and before anything else, so it exists for every sliced beat.
        // The RR-length fence in jbank::cleanGroups rejects on the LONG side
        // only for this reason: a premature beat is short by definition, so a
        // short-side rejection would throw ectopics out as length outliers for
        // outlier reasons, with no flag and no record.
        //
        // premature: RR(i) < 0.80 * median of the trailing ten.
        // voted:     the 5-of-8 rule over the raw premature flags -- the
        //            middle of a run, where the trailing median has itself
        //            gone short so the beat stops reading as premature alone.
        std::vector<char> premature;
        std::vector<char> voted;

        // ---- ROW POSITION AT SLICE TIME, parallel to `beats` --------------
        //
        // NOT AN R-PAIR ORDINAL, which is what slice_index below is for. This
        // is seeded as the identity over `beats` AFTER the slicing loop has
        // already skipped pairs, so it counts survivors, not pairs. It was
        // called original_index, which collided with PpgBeatSet's field of that
        // name -- and that one does hold ordinals.
        std::vector<size_t> row_at_slice_time;

        // ---- THE R-PAIR ORDINAL, parallel to `beats` ----------------------
        //
        // slice_index[k] is the value of the slicing loop counter that produced
        // beat k: the ordinal of the pair (rPeaks[i], rPeaks[i+1]). It is the
        // ONLY quantity in this struct that means the same thing on every
        // channel, because every channel's slicer is driven by the same
        // ch1.raw R-peak vector.
        //
        // NOT THE SAME AS row_at_slice_time, and the difference is silent. The
        // slicing loop SKIPS pairs -- rr <= 3 samples, and rr > 4 s, which is a
        // dropout gap rather than a beat -- with `continue`, before anything is
        // pushed. So `beats` is already compacted against the R-pair list by
        // the time row_at_slice_time is seeded as the identity over it. On a bin
        // where nothing was skipped the two coincide; on a bin with one dropout
        // gap every later beat's ordinal is short by one, and the error grows
        // with each skip.
        //
        // That is exactly the quantity jbank::ChannelSet keys on, and using
        // row_at_slice_time in its place pairs an ECG complex with an unrelated
        // pulse, progressively further through the bin. The pulse slicer below
        // has always recorded this (PpgBeatSet::original_index, from
        // raw_slice); the ECG path did not, which is why the two could not
        // actually be joined on a bin containing a gap.
        std::vector<uint32_t> slice_index;
        // Counts as they stood BEFORE pruning, so the flag survivors can be
        // compared against what the record actually contained.
        int n_premature_presliced = 0;
        int n_voted_presliced = 0;

        // Parallel to `beats`: PQ_level - TP_level whenever BOTH pass
        // estimates succeeded for that beat (NaN otherwise). QC metric per
        // spec I-1; not used to drive any decision, just recorded.
        std::vector<double> tp_pq_delta;
        // Parallel to `beats`: the per-pass vertical DC shifts applied by the
        // two-pass leveling. tp_shift = Pass 2 (TP) amount; pq_shift =
        // Pass 3 (PQ finalize) amount. NaN where that pass didn't apply.
        // Returned for the per-beat move log.
        std::vector<double> tp_shift;
        std::vector<double> pq_shift;
    };

    // ---- A SPLICED BIN, AND WHY A LONG RR IS DROPPED ---------------------
    //
    // An annealed bin is not one continuous stretch of recording. The annealer
    // excises the noise-marked regions and redistributes the surviving good
    // fragments between neighbouring bins, then builds each bin's signal by
    // CONCATENATING those fragments end to end (see anneal_handler's writer and
    // the ecg_bin_indexs it emits). So two samples adjacent in this array can
    // be seconds apart in the recording.
    //
    // WHICH MAKES AN RR ACROSS A JOIN MEANINGLESS. The loop below computes
    // rr = rPeaks[i+1] - rPeaks[i] and treats it as a cardiac interval. For an
    // R-pair straddling a splice that number is the distance across the join,
    // and because fragments are redistributed between bins it can be large --
    // so the slicer would cut a 0.5*rr before / 1.3*rr after window from it and
    // produce one "beat" spanning several cardiac cycles. That is the multi-QRS
    // template with the over-long axis, and the column median out past the real
    // RR is that slice's LATER complexes.
    //
    // THE RR CAP IS WHAT HANDLES IT, and it is the only thing that does. An
    // earlier attempt derived the splice positions from ecg_bin_indexs and
    // dropped pairs crossing one (FragmentSeams / pair_crosses_seam, both
    // removed): switching it on emptied every bin, because the positions it
    // derived did not describe this signal and nearly every pair looked
    // spliced. The cap below covers the same case from the other direction --
    // a join shows up as a long RR, and a long RR is dropped whatever caused
    // it -- without needing to know where the joins are.
    //
    inline ecg_beat_set extract_beats_and_align(const std::vector<double>& signal, const std::vector<size_t>& rPeaks, double fs) {
        ecg_beat_set out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) {
            fprintf(stderr, "[align] no beats: signal=%lld rPeaks=%zu\n",
                static_cast<long long>(N), rPeaks.size());
            return out;
        }

        // ==================================================================
        // ALIGNMENT DOES NOT PRUNE, AT ALL. It slices, flags rhythm, aligns
        // horizontally and levels vertically. Every R-pair that survived the
        // slicing guards above is in `beats` when this function returns, so
        // row_at_slice_time is the identity over it.
        //
        // THE TUKEY FENCES LIVE IN jbank::cleanGroups, deliberately. The RR,
        // amplitude, R-location and wave-score fences all reject beats that do
        // not resemble the population, and before partitioning the population
        // is every morphology mixed together -- a PVC is an outlier by
        // construction and would be deleted before the bank could ever see it.
        // That ordering error accounted for the two-member templates, the
        // over-segmentation into _B/_C/_D, and the ECG and PPG survivor counts
        // disagreeing by a factor of ninety on the same R-peaks.
        //
        // cleanGroups runs the same 1.5*IQR convention (keep_within_tukey, in
        // stats_utils.hpp) per GROUP, after the partition is final, so
        // each fence asks whether a beat is an outlier among beats of its OWN
        // shape. It measures amplitude, R-location and template correlation
        // itself and owns the per-beat verdict that reaches the morphology CSV.
        // Alignment used to compute its own copies of those three fences and
        // record verdicts nothing read; they are gone.
        // ==================================================================

        // MAX RR = 2.5 s, a hard drop. An R-pair longer than this is not a
        // beat: it is a missed detection, or a gap left where the annealer
        // excised noise. The slice cut from one spans
        // 2.5 * (percent_interval_preceeding_rpeak +
        // percent_interval_following_rpeak) = 2.5 * 1.9 = 4.75 s, which still
        // covers several cardiac cycles -- the cap bounds the damage rather
        // than guaranteeing one complex per row.
        //
        // IT ALSO SETS THE BIN'S FRAME WIDTH. shared_w is sized on
        // max_rr_len, so ONE oversized survivor stretches every template in
        // the bin; a sparse bin with a handful of real detections could
        // produce a many-hundred-second template. That is why this is a cap
        // and not left to the Tukey fences, which measure but no longer prune.
        //
        // THIS IS ALSO THE ONLY THING THAT HANDLES A SPLICED R-PAIR. See the
        // block above: a join shows up as a long RR, and a long RR is dropped
        // whatever caused it.
        const int64_t kMaxBeatSamplesEcg = (fs > 0.0) ? static_cast<int64_t>(2.5 * fs) : 0;
        // ---- slice every beat ------------------------------------------
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t rr = static_cast<int64_t>(rPeaks[i + 1]) - r0;
            if (rr <= 3) continue;
            if (kMaxBeatSamplesEcg > 0 && rr > kMaxBeatSamplesEcg) continue;

            const int64_t before = rr_before_samples(rr);
            const int64_t after = rr_after_samples(rr);
            const int64_t len = before + after;
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;

            std::vector<double> beat(static_cast<size_t>(len),
                std::numeric_limits<double>::quiet_NaN());
            const int64_t cs = std::max<int64_t>(0, start);
            const int64_t ce = std::min<int64_t>(N, end);
            for (int64_t k = cs; k < ce; ++k)
                beat[static_cast<size_t>(k - start)] = signal[static_cast<size_t>(k)];

            out.beats.push_back(std::move(beat));
            out.r_indices.push_back(rPeaks[i]);
            out.rr_lens.push_back(static_cast<int>(rr));
            // Pushed in the SAME statement group as the beat and after every
            // `continue` above -- the pulse slicer's raw_slice carries the same
            // note for the same reason: a beat and its ordinal appended
            // separately desynchronise at the first skipped pair, and the
            // result still looks like a valid parallel pair.
            out.slice_index.push_back(static_cast<uint32_t>(i));
        }
        if (out.beats.empty()) return out;
        out.baseline_source.assign(out.beats.size(), BaselineSource::NONE);

        // ---- Section 4.6 rhythm verdict, parallel to `beats` -------------
        // out.rr_lens[i] is beat i's own interval, so the intervals this reads
        // are exact: one per beat, none missing. Assigned here, on the full
        // sliced set, so every beat has a verdict whatever cleanGroups later
        // removes.
        {
            const size_t nb = out.beats.size();
            out.premature.assign(nb, 0);
            // Seeded alongside the rhythm flags, so every sliced beat has an
            // entry.
            out.row_at_slice_time.resize(out.beats.size());
            for (size_t i = 0; i < out.row_at_slice_time.size(); ++i) out.row_at_slice_time[i] = i;
            out.voted.assign(nb, 0);
            if (nb >= 12) {
                // Beat t is premature when the interval BEFORE it is short.
                // rr_lens[i] holds the interval AFTER beat i (R[i+1] - R[i]),
                // because that is the span the slice covers -- so beat t's
                // PRECEDING interval is rr_lens[t-1]. Testing rr_lens[t]
                // directly flags the beat that precedes each PVC instead of
                // the PVC, which is an off-by-one that looks exactly like a
                // detector that "nearly works".
                for (size_t t = 11; t < nb; ++t) {
                    std::vector<double> w(out.rr_lens.begin() + (t - 11),
                        out.rr_lens.begin() + (t - 1));   // the ten before it
                    std::sort(w.begin(), w.end());
                    const double med = w[w.size() / 2];
                    if (med > 0.0 && out.rr_lens[t - 1] < 0.80 * med)
                        out.premature[t] = 1;
                }
                // The vote reads the RAW flags, never its own output: feeding
                // it back would let one run grow along the whole record.
                for (size_t t = 0; t < nb; ++t) {
                    const size_t lo = (t > 4) ? t - 4 : 0;
                    const size_t hi = std::min(nb, t + 4);
                    int c = 0;
                    for (size_t i = lo; i < hi; ++i) c += out.premature[i];
                    if (c >= 5) out.voted[t] = 1;
                }
                for (size_t i = 0; i < nb; ++i) {
                    out.n_premature_presliced += out.premature[i];
                    if (out.voted[i] && !out.premature[i]) ++out.n_voted_presliced;
                }
            }
        }

        out.total_beats = out.beats.size();

        // the Sangala document says use mode, in this case we use median length
        {
            std::vector<int> lens = out.rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
        }

        // ---- reference beat = first median-length beat -----------------
        // Used only as the PR-baseline DC-alignment reference. No RMS shape
        // clustering: the final template is the column-wise median.
        for (size_t i = 0; i < out.beats.size(); ++i) {
            if (out.rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
            }
        }
        // ---- Pass 1: R-align on shared axis ----------------------------
        // FRAMED ON THE MAXIMUM RR. This is a STORAGE width, and it is sized so
        // that NO BEAT LOSES A SAMPLE: the longest beat in the bin has to fit,
        // or its tail is gone from the matrix and from every template built
        // from it.
        //
        int max_rr_len = 0;
        for (int L : out.rr_lens)
            if (L > max_rr_len) max_rr_len = L;
        if (max_rr_len <= 0) return out;

        const int R_anchor = static_cast<int>(rr_before_samples(max_rr_len));
        const int shared_w = R_anchor + static_cast<int>(rr_after_samples(max_rr_len));
        if (shared_w <= 0) return out;

        const double NaND = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> aligned;
        aligned.reserve(out.beats.size());
        for (size_t i = 0; i < out.beats.size(); ++i) {
            const auto& b = out.beats[i];
            const int L = static_cast<int>(b.size());
            const int r_in_beat = static_cast<int>(rr_before_samples(out.rr_lens[i]));
            const int prepend = R_anchor - r_in_beat;
            std::vector<double> a(shared_w, NaND);
            for (int k = 0; k < L; ++k) {
                const int dst = prepend + k;
                if (dst >= 0 && dst < shared_w) a[dst] = b[k];
            }
            aligned.push_back(std::move(a));
        }
        out.beats = std::move(aligned);
        out.r_aligned_col = R_anchor;


        // ---- Passes 2 and 3: two-pass TP-then-PQ vertical DC alignment --
        // Pass 2 (TP): this beat's own T-end -> the NEXT beat's P-onset.
        // Pass 3 (PQ): this beat's own P-end -> this beat's own Q-onset.
        // PQ is the higher-priority reference (spec): when BOTH estimates
        // are available, the beat is shifted using PQ's level, not TP's --
        // TP is used only when PQ's landmarks aren't usable. Both are always
        // attempted whenever possible (not short-circuited on TP success)
        // so the TP-to-PQ delta can be recorded as a per-beat QC metric even
        // on beats where PQ ends up winning.
        //
        // Both estimators report level as the MEDIAN over the flattest
        // (lowest-variance) min_w-wide sub-window of their respective
        // landmark-bounded segment -- flatness search still uses variance
        // (cheapest way to find "flat"), only the reported level is median
        // rather than mean.
        //
        // The next beat's own R, in this beat's own SHARED aligned column
        // space, sits at R_anchor + rr_lens[i] (this beat's own R was moved
        // to R_anchor by Pass 1; the next R was the same distance further in
        // the original slice, so the same shift lands it there) -- pure
        // arithmetic, no extra detection needed to locate it.
        //
        // Degrades gracefully: baseline_source[i] records which segment (if
        // any) supplied this beat's baseline; tp_pq_delta[i] records
        // PQ_level - TP_level whenever both succeeded (NaN otherwise). NONE
        // means neither was usable; that beat is left un-shifted, and
        // CreateEcgTemplates.hpp excludes it from the per-sample amplitude
        // aggregation (median/std) rather than silently contribute an
        // unreliable, unadjusted amplitude.
        if (!out.beats.empty())
        {
            const int min_w = std::max(1, static_cast<int>(std::lround(0.010 * fs)));  // 10 ms window
            // RR-fraction bounds for the TP baseline window, and ms-before-R
            // bounds for the PQ baseline window. Conservative defaults; tune
            // against real recordings if the leveling looks off.
            const double kTpLoFrac = 0.55;   // start of TP flat, fraction of RR after R
            const double kTpHiFrac = 0.85;   // end of TP flat, before next P
            const double kPqPreRMs = 80.0;   // PQ window opens this many ms before R
            const double kPqGuardMs = 20.0;  // and closes this many ms before R (guard vs Q)

            // Median over the lowest-variance min_w-wide sub-window of
            // [lo, hi] (inclusive, in this beat's own aligned column space).
            auto flattest_median = [&](const std::vector<double>& beat, int lo, int hi)
                -> std::pair<double, bool>
                {
                    const int hi_bound = std::min(hi, static_cast<int>(beat.size()) - min_w);
                    if (hi_bound < lo)
                        return { std::numeric_limits<double>::quiet_NaN(), false };

                    double best_var = std::numeric_limits<double>::infinity();
                    int best_start = -1;
                    for (int s = lo; s <= hi_bound; ++s) {
                        double sum = 0.0, sumsq = 0.0; int n = 0;
                        for (int k = s; k < s + min_w; ++k) {
                            const double v = beat[k];
                            if (std::isnan(v)) continue;
                            sum += v; sumsq += v * v; ++n;
                        }
                        if (n < static_cast<int>(min_w * 0.7)) continue;
                        const double mean = sum / n;
                        const double var = sumsq / n - mean * mean;
                        if (var < best_var) { best_var = var; best_start = s; }
                    }
                    if (best_start < 0)
                        return { std::numeric_limits<double>::quiet_NaN(), false };

                    std::vector<double> win;
                    win.reserve(min_w);
                    for (int k = best_start; k < best_start + min_w; ++k)
                        if (!std::isnan(beat[k])) win.push_back(beat[k]);
                    if (win.empty())
                        return { std::numeric_limits<double>::quiet_NaN(), false };
                    // Called three times per beat (both_levels twice, then the
                    // PQ re-measure after pass 2's shift), so per bin per method
                    // per channel this is ~3x the beat count. One order
                    // statistic, so nth_element.
                    const size_t m = win.size() / 2;
                    std::nth_element(win.begin(), win.begin() + m, win.end());
                    const double hi_mid = win[m];
                    const double med = (win.size() % 2 == 0)
                        ? 0.5 * (*std::max_element(win.begin(), win.begin() + m)
                            + hi_mid)
                        : hi_mid;
                    return { med, true };
                };

            /* Pass 2: TP Window. The baseline segment is the TP window ESTIMATED
            * by range:
            * lo = R_anchor + 0.55 · RR
            * hi = R_anchor + 0.85 · RR
            * And then the flattest 10ms wide segment is located via a minimum variance scan
            */
            auto tp_window = [&](size_t i) -> std::pair<int, int> {
                const int N = static_cast<int>(out.beats[i].size());
                const int rr = out.rr_lens[i];
                if (rr <= 0) return { -1, -1 };
                const int lo = R_anchor + static_cast<int>(std::lround(kTpLoFrac * rr));
                const int hi = R_anchor + static_cast<int>(std::lround(kTpHiFrac * rr));
                const int clo = std::max(0, lo);
                const int chi = std::min(N, hi);
                if (chi - clo < 3) return { -1, -1 };
                return { clo, chi };
                };

            // Pass 3, PQ window: Segment is always 80 ms before R to 20 ms before R
            auto pq_window = [&](size_t i) -> std::pair<int, int> {
                const int N = static_cast<int>(out.beats[i].size());
                const int lo = R_anchor - static_cast<int>(std::lround(kPqPreRMs * 0.001 * fs));
                const int hi = R_anchor - static_cast<int>(std::lround(kPqGuardMs * 0.001 * fs));
                const int clo = std::max(0, lo);
                const int chi = std::min(N, hi);
                if (chi - clo < 3) return { -1, -1 };
                return { clo, chi };
                };

            // Compute both pass estimates for beat i (never short-circuited,
            // so the delta can be recorded even when PQ will be used).
            auto both_levels = [&](size_t i) -> std::tuple<double, bool, double, bool> {
                double tpLvl = std::numeric_limits<double>::quiet_NaN(), pqLvl = tpLvl;
                bool tpOk = false, pqOk = false;
                const auto [tlo, thi] = tp_window(i);
                if (tlo >= 0) { const auto [lvl, ok] = flattest_median(out.beats[i], tlo, thi); tpLvl = lvl; tpOk = ok; }
                const auto [plo, phi] = pq_window(i);
                if (plo >= 0) { const auto [lvl, ok] = flattest_median(out.beats[i], plo, phi); pqLvl = lvl; pqOk = ok; }
                return { tpLvl, tpOk, pqLvl, pqOk };
                };

            /*Two-pass vertical alignment. - Pass 2: level each beat on the TP segment (end of T to just before onset of
            next P). Estimate level with median over the flattest sub-window bounded by fitted landmarks. - Pass
            3: finalize the zero on the PQ segment (end of P to immediately before Q onset). PQ is the higher-priority reference.*/

            // ---- EVERY BEAT'S LEVELS, MEASURED ONCE -------------------
            //
            // both_levels is two flattest_median scans and each of those is a
            // variance sweep over its window, so it is the expensive call in
            // this function. Held here rather than recomputed: the median
            // below needs every beat's levels, and so does the apply loop
            // further down. Two calls per beat total, which is what the
            // previous arrangement already cost -- one for the reference, one
            // per beat in the apply loop, plus a borrow pass when the
            // reference came up short.
            const size_t nb = out.beats.size();
            std::vector<double> tpLvlOf(nb, std::numeric_limits<double>::quiet_NaN());
            std::vector<double> pqLvlOf(nb, std::numeric_limits<double>::quiet_NaN());
            std::vector<char>   tpOkOf(nb, 0), pqOkOf(nb, 0);
            for (size_t i = 0; i < nb; ++i) {
                const auto [tpLvl, tpOk, pqLvl, pqOk] = both_levels(i);
                tpLvlOf[i] = tpLvl;  tpOkOf[i] = tpOk ? 1 : 0;
                pqLvlOf[i] = pqLvl;  pqOkOf[i] = pqOk ? 1 : 0;
            }

            // ---- THE TARGET IS A MEDIAN OVER THE BIN, NOT ONE BEAT ----
            //
            // WHAT THIS REPLACED, AND WHY. The target used to be the levels of
            // ONE beat -- ref_beat_index, the first beat whose RR happened to
            // equal the bin's median RR. That beat was selected for framing
            // (median_length is what sizes the matrix so no beat is clipped)
            // and then reused as the level every other beat is matched to,
            // which it was never chosen for: it could be the noisiest beat in
            // the bin, or an ectopic one with an ordinary interval, and nothing
            // downstream could tell. A borrow pass existed underneath it for
            // when that one beat had no usable segment, which is a fallback for
            // a problem the median does not have.
            //
            // The choice of target is arbitrary in one sense -- shift every
            // beat by the same extra amount and they are still mutually
            // aligned, only the bin's absolute level moves -- so the median is
            // taken to keep that level where a single-beat target put it,
            // rather than moving the whole amplitude scale. Levelling to zero
            // would work equally well for the alignment and is arguably the
            // better convention (the PQ segment IS the isoelectric zero), but
            // it shifts every raw template amplitude and anything reading one.
            auto median_of = [](std::vector<double> v) -> double {
                if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
                const size_t m = v.size() / 2;
                std::nth_element(v.begin(), v.begin() + m, v.end());
                const double hi = v[m];
                if (v.size() % 2) return hi;
                return 0.5 * (*std::max_element(v.begin(), v.begin() + m) + hi);
                };
            std::vector<double> tpPool, pqPool;
            tpPool.reserve(nb);  pqPool.reserve(nb);
            for (size_t i = 0; i < nb; ++i) {
                if (tpOkOf[i]) tpPool.push_back(tpLvlOf[i]);
                if (pqOkOf[i]) pqPool.push_back(pqLvlOf[i]);
            }
            const double refTpTarget = median_of(tpPool);
            const double refPqTarget = median_of(pqPool);
            const bool haveAnyRef = !std::isnan(refTpTarget) || !std::isnan(refPqTarget);

            // Which segment the passes below level against, for the log line.
            const BaselineSource refSrc =
                !std::isnan(refPqTarget) ? BaselineSource::PQ
                : !std::isnan(refTpTarget) ? BaselineSource::TP
                : BaselineSource::NONE;

            out.tp_pq_delta.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());

            // Diagnostic accumulators (spec-neutral: recorded only, doesn't
            // change any beat).
            int diag_pq = 0, diag_tp = 0, diag_none = 0;
            std::vector<double> diag_shifts;
            diag_shifts.reserve(out.beats.size());

            out.tp_shift.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());
            out.pq_shift.assign(out.beats.size(), std::numeric_limits<double>::quiet_NaN());
            if (haveAnyRef) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    // From the pass above, not re-measured.
                    const double tpLvl = tpLvlOf[i];
                    const double pqLvl = pqLvlOf[i];
                    const bool   tpOk = tpOkOf[i] != 0;
                    const bool   pqOk = pqOkOf[i] != 0;
                    if (tpOk && pqOk) out.tp_pq_delta[i] = pqLvl - tpLvl;

                    // Pass 2 (TP): coarse level, applied only if BOTH this
                    // beat and the reference have a usable TP.
                    double applied = 0.0;
                    bool leveled = false;
                    if (tpOk && !std::isnan(refTpTarget)) {
                        const double d = tpLvl - refTpTarget;
                        for (double& v : out.beats[i]) if (!std::isnan(v)) v -= d;
                        applied += d; leveled = true;
                        out.tp_shift[i] = d;
                    }
                    // Pass 3 (PQ finalize): authoritative zero. Recompute PQ
                    // level AFTER pass 2's shift, match it to the reference PQ.
                    if (!std::isnan(refPqTarget)) {
                        const auto [plo, phi] = pq_window(i);
                        if (plo >= 0) {
                            const auto [lvl2, ok2] = flattest_median(out.beats[i], plo, phi);
                            if (ok2) {
                                const double d = lvl2 - refPqTarget;
                                for (double& v : out.beats[i]) if (!std::isnan(v)) v -= d;
                                applied += d; leveled = true;
                                out.pq_shift[i] = d;
                                out.baseline_source[i] = BaselineSource::PQ;
                            }
                        }
                    }
                    if (!leveled) { out.baseline_source[i] = BaselineSource::NONE; continue; }
                    if (out.baseline_source[i] != BaselineSource::PQ)
                        out.baseline_source[i] = BaselineSource::TP;   // TP-only
                    diag_shifts.push_back(applied);
                    if (out.baseline_source[i] == BaselineSource::PQ) ++diag_pq;
                    else ++diag_tp;
                }
                diag_none = (int)out.beats.size() - diag_pq - diag_tp;
            }
            // One line per aligned matrix. Shows whether shifts are sane
            // (small std, tight range) or wild (large std, big range).
            {
                const char* refTag = (refSrc == BaselineSource::PQ) ? "PQ"
                    : (refSrc == BaselineSource::TP) ? "TP" : "NONE";
                double smin = std::numeric_limits<double>::infinity();
                double smax = -std::numeric_limits<double>::infinity();
                double ssum = 0.0, ssumsq = 0.0;
                for (double s : diag_shifts) {
                    if (s < smin) smin = s;
                    if (s > smax) smax = s;
                    ssum += s; ssumsq += s * s;
                }
                const size_t nS = diag_shifts.size();
                const double smean = (nS > 0) ? ssum / nS : 0.0;
                const double sstd = (nS > 1)
                    ? std::sqrt(std::max(0.0, ssumsq / nS - smean * smean)) : 0.0;
                double smed = 0.0;
                if (!diag_shifts.empty()) {
                    auto tmp = diag_shifts;
                    std::sort(tmp.begin(), tmp.end());
                    smed = tmp[tmp.size() / 2];
                }
            }
        }
        out.q_aligned_col = out.r_aligned_col;
        return out;
    }

    // =========================================================================
    // Q-align a bin's cached (R-aligned) snippets and re-median.
    //
    // Mirrors R-alignment's structure: R-align picks a reference "median
    // snippet" (the first beat whose finite length equals the median) and
    // aligns to it. Here we do the same, but on Q. The alignment reference is
    // the location of the R peak on the beat of median length.  
    // find_q only looks in the R-upstroke region. We find Q on that median, then
    // shift every snippet so its own Q lands on that marker, and re-median the
    // shifted snippets into the Q template. The snippets stay full length
    // (never cropped); the window only bounds where the marker is found.
    //
    // Because the snippets are re-aligned on Q, each beat's R moves by its own
    // shift, so the Q template's R spike sits at R_anchor + median(shift). That
    // column is returned in r_col -- computed from the applied shifts, not a
    // window search (which would risk landing on Q or S).
    //
    // R_anchor is the R-pass R column; fs is the ECG rate; refMedian is the
    // R-aligned template used to locate the Q marker. compute_iqr is retained
    // for call-site compatibility but no longer gates anything -- the gray-band
    // IQR (a robust spread, not a standard deviation) is always filled.
    // q_aligned_col is the marker Q landed on; r_col is the R fiducial (R_anchor).
    struct aligned_beats {
        std::vector<double> tmpl;
        std::vector<double> iqr;
        std::vector<std::vector<double>> beats;
        int q_aligned_col = -1;
        int r_col = -1;
    };

    // exclude_from_median: one entry per row of beatsIn, non-zero meaning the
    // row is SHIFTED but not averaged. Null or short means average everything,
    // which is the behaviour this had before the parameter existed.
    //
    // WHY IT HAS TO BE PASSED IN. This function re-shifts an already-sliced
    // matrix and takes a column median; it has no verdicts of its own and never
    // will, because it does not slice, level, or look for a landmark on a beat.
    // extract_beats_and_align computes all of that -- baseline_source, the four
    // Tukey fences, prematurity -- and deliberately prunes NOTHING, because 4.6
    // requires flagged beats retained in the record.
    //
    // So the exclusion used to be applied at the point of averaging, in
    // create_ecg_templates' `usable` gate. When the per-anchor averages moved
    // here, the gate did not move with them: this median drew on every beat,
    // including the ones every CSV already reported as excluded. A dropped R
    // detection is the visible case -- its rr spans several cardiac cycles so
    // its slice is 1.8x that, it is the ONLY row with samples in the far tail,
    // and the column median out there is its later QRS complexes at full
    // amplitude. Bin 5 slot A drew four R peaks over a 2.84 s axis on 986
    // members while _templates.csv reported 28 of them excluded.
    //
    // SHIFTED, NOT DROPPED. The row stays in the matrix and keeps its place, so
    // nothing downstream that indexes by row is disturbed and the beat is still
    // in the record. Only the median skips it.
    inline aligned_beats align_beat_matrix(
        const std::vector<std::vector<double>>& beatsIn,
        int R_anchor, double fs, bool compute_iqr,
        const std::vector<double>& ref_beat_of_median_length,
        const std::function<double(const std::vector<double>&)>& locate,
        // The per-beat correlation floor, tbank::morphThresholdEcg() from
        // config.csv. PASSED IN, not read here, so this header does not depend
        // on template_bank.hpp. No default -- a silent 0.0 would disable the
        // guard and let every beat through however badly it correlated.
        double corrFloor,
        const std::vector<char>* exclude_from_median = nullptr)
    {
        aligned_beats res;
        if (beatsIn.empty()) return res;

        std::vector<std::vector<double>> beats = beatsIn;
        const int Wsh = static_cast<int>(beats.front().size());

        // Anchor = the landmark on the column-wise median of all beats (for
        // R_PEAK that is the R template's own R). Every snippet is then shifted
        // so its own landmark lands there.
        //
        // ---- SUB-SAMPLE, AND THAT IS THE WHOLE POINT ---------------------
        //
        // This used to lround BOTH ends and shift by whole samples. The
        // per-beat variability of P onset, Q onset or the J point RELATIVE TO R
        // is a few milliseconds -- under one sample at typical rates -- so
        // every shift rounded to zero, no beat moved, and all four anchor
        // alignments came out bit-identical to the R one. The locators return
        // genuine fractional positions (291.31, 302.916, ...); lround threw
        // away exactly the precision that distinguishes the alignments.
        //
        // So the shift is a double and the resample is linear between the two
        // straddling samples. Linear rather than cubic deliberately: the shifts
        // are a fraction of a sample, where linear and cubic agree to within
        // far less than the beat-to-beat spread being measured, and linear
        // cannot overshoot -- it will not invent a taller R peak than the beat
        // actually had.
        const double markerD = locate(ref_beat_of_median_length);
        const int marker = (markerD >= 0.0)
            ? static_cast<int>(std::lround(markerD)) : -1;
        res.q_aligned_col = (marker >= 0) ? marker : R_anchor;

        std::vector<int> shifts;   // per-beat shift, rounded, for R's new column
        std::vector<double> locatedPos;   // DIAG: raw per-beat located landmark (sub-sample)
        if (markerD >= 0.0) {
            const double NaNv = std::numeric_limits<double>::quiet_NaN();
            shifts.reserve(beats.size());
            // ---- PER-BEAT SHIFT BY CROSS-CORRELATION --------------------
            //
            // This ran `locate(beats[i])` -- the FULL landmark finder, with its
            // Q-peak search, sigma-4 refine, 4x cubic upsample and four-model
            // BIC selection -- once per beat, and used the result for exactly
            // one thing: markerD - mi. Cost scaled with window width rather
            // than beat count (~20.9 ms/beat at 40 samples, ~268.9 ms at 200),
            // which is why the P pass ran ~6.9x the Q pass.
            //
            // The landmark is measured ONCE, on the reference beat, above
            // (markerD). A beat only owes the alignment its offset from that,
            // and xcorrShift returns it at ~0.011 ms/beat, agreeing with the
            // per-beat fit to 0.009 ms. Per-beat model switching goes with it:
            // BIC picking different models on different beats put the ~5 ms
            // inter-model bias difference into the beat-to-beat spread.
            //
            // WINDOW: centred on the reference landmark, +-halfWin. It has to
            // contain the transition being aligned and nothing else -- a window
            // holding two transitions locks onto the larger one.
            //
            // FLOOR: tbank::morphThresholdEcg(), the config.csv split threshold
            // the bank already uses to decide whether a beat belongs to a
            // morphology. Below it there is NO estimate and the beat is
            // skipped, exactly as a locator returning -1 was skipped.
            const int halfWin = std::max(5,
                static_cast<int>(std::lround(0.030 * fs)));   // +-30 ms
            const int xlo = static_cast<int>(std::lround(markerD)) - halfWin;
            const int xhi = static_cast<int>(std::lround(markerD)) + halfWin;
            size_t nNoCorr = 0;
            for (size_t i = 0; i < beats.size(); ++i) {
                const double shiftD = upsample_for_fit::xcorrShift(
                    beats[i], ref_beat_of_median_length, xlo, xhi, corrFloor);
                if (!std::isfinite(shiftD)) { ++nNoCorr; continue; }
                locatedPos.push_back(markerD - shiftD);   // DIAG: implied landmark
                shifts.push_back(static_cast<int>(std::lround(shiftD)));

                // A shift under a thousandth of a sample is not worth a
                // resample -- it would only add interpolation error.
                if (std::abs(shiftD) < 1e-3) continue;

                // dst = k + shiftD, so the value landing on integer column j
                // comes from source position j - shiftD. Reading BACKWARDS from
                // the destination is what makes this a resample rather than a
                // scatter: every output column gets exactly one value, with no
                // holes where a forward scatter would skip a column and no
                // collisions where two sources would round to the same one.
                std::vector<double> a(Wsh, NaNv);
                for (int j = 0; j < Wsh; ++j) {
                    const double src = static_cast<double>(j) - shiftD;
                    const int s0 = static_cast<int>(std::floor(src));
                    const int s1 = s0 + 1;
                    if (s0 < 0 || s1 >= Wsh) continue;   // off the window: stays NaN
                    const double v0 = beats[i][s0];
                    const double v1 = beats[i][s1];
                    // NaN-safe: a beat is NaN-padded at its edges, and
                    // interpolating across that boundary would smear the pad
                    // inward. Either neighbour missing leaves the column NaN,
                    // which the column-wise median below already skips.
                    if (std::isnan(v0) || std::isnan(v1)) continue;
                    const double f = src - static_cast<double>(s0);
                    a[j] = v0 + f * (v1 - v0);
                }
                beats[i] = std::move(a);
            }
        }

        // Column-wise median (+ optional IQR spread) over the Q-aligned beats.
        // One gather per column, reused scratch, nth_element for the median.
        res.tmpl.assign(Wsh, std::numeric_limits<double>::quiet_NaN());
        res.iqr.assign(Wsh, 0.0);

        // Hoisted out of the column loop: the same rows are skipped at every
        // column, so the test is per row and not per (row, column).
        auto skipRow = [&](size_t r) {
            return exclude_from_median && r < exclude_from_median->size()
                && (*exclude_from_median)[r] != 0;
            };

        std::vector<double> col;
        col.reserve(beats.size());
        for (int c = 0; c < Wsh; ++c) {
            col.clear();
            for (size_t r = 0; r < beats.size(); ++r) {
                if (skipRow(r)) continue;
                const double v = beats[r][c];
                if (!std::isnan(v)) col.push_back(v);
            }
            const size_t nc = col.size();
            if (nc == 0) continue;

            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            res.tmpl[c] = (nc % 2 == 0)
                ? 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi)
                : hi;

            if (nc >= 2) {
                // Per-column standard deviation (matches the QC noise metric,
                // which also uses std). Computed over the same non-NaN values
                // used for the median. Field is still named 'iqr' downstream
                // (serializer/viewer) but now carries std.
                double mean = 0.0;
                for (double v : col) mean += v;
                mean /= static_cast<double>(nc);
                double ss = 0.0;
                for (double v : col) { const double d = v - mean; ss += d * d; }
                res.iqr[c] = std::sqrt(ss / static_cast<double>(nc - 1));  // sample std
            }
        }

        // Passed-in R, tracked through the Q-align shift (median of the applied shifts).
        int med_shift = 0;
        if (!shifts.empty()) {
            std::sort(shifts.begin(), shifts.end());
            med_shift = shifts[shifts.size() / 2];
        }
        int rc = std::clamp(R_anchor + med_shift, 0, std::max(0, Wsh - 1));

        // Snap to the local peak within +/-5 ms of that passed-in position.
        // The window is far too tight to reach Q or S, so it only cleans up
        // sub-window drift of the R spike.
        const int width_to_curve_fit_r_peak = std::max(1, static_cast<int>(std::lround(0.040 * fs)));
        const int rlo = std::max(0, rc - width_to_curve_fit_r_peak);
        const int rhi = std::min(Wsh - 1, rc + width_to_curve_fit_r_peak);
        double rbest = -std::numeric_limits<double>::infinity();
        for (int i = rlo; i <= rhi; ++i) {
            if (!std::isnan(res.tmpl[i]) && res.tmpl[i] > rbest) { rbest = res.tmpl[i]; rc = i; }
        }
        res.r_col = rc;

        res.beats = std::move(beats);
        return res;
    }

    // =========================================================================
    struct PpgBeatSet {
        std::vector<std::vector<double>> beats;   // NaN-padded, 50%-upslope-aligned
        std::vector<int> peak_cols;               // per-beat systolic peak column (varies)
        std::vector<int> foot_cols;               // per-beat foot column (varies)

        // ---- WHERE THIS BEAT'S OWN R PEAK LANDED IN THE SHARED FRAME -----
        //
        // THE ONLY WAY TO PUT A PULSE ON AN R-RELATIVE TIME AXIS, and it was
        // being thrown away. Every consumer that needs "how long after the
        // QRS did this happen" -- the viewer's shared ECG/pulse time axis, the
        // construction-time peak/foot brackets in
        // create_arterial_templates -- was instead ASSUMING a fixed lead-in of
        // padSeconds * rate, and there is no fixed lead-in:
        //
        //   * a beat is sliced at r0 - rr_before_samples(rr), i.e.
        //     0.5 * THAT BEAT'S OWN RR, so its R column varies with its RR;
        //   * every beat is then shifted so its up50 lands on up50_anchor, so
        //     the R columns move again, by a different amount each.
        //
        // The consequence of the assumption was a pulse drawn at a time offset
        // of roughly (0.5 * RR_median - padSeconds) seconds -- zero near
        // 75 bpm, +200 ms at 50 bpm, -100 ms at 100 bpm. It looked like a
        // constant that needed tuning and was a quantity that had to be
        // measured.
        //
        // Parallel to `beats`, like peak_cols and foot_cols, and subject to
        // the same caveat: a heavily shifted beat can have its R clipped off
        // the left edge, so an entry may be < 0. Take the MEDIAN over the
        // beats that actually built a template, not the mean.
        std::vector<int> r_cols;

        // R-PAIR ORDINAL of each surviving beat, parallel to `beats`.
        //
        // WHY IT HAS TO BE CARRIED. This set is pruned twice -- beats with no
        // detectable peak/foot/up50 are skipped at slice time, and peak-column
        // outliers are removed in pass 1 -- so row k of `beats` is NOT R-pair k.
        // Every consumer that needs to say "this pulse and that QRS are the same
        // heartbeat" needs the ordinal, and it was being discarded here, which
        // is why a partition shared between ECG and PPG could not be built at
        // all: there was no key to join them on.
        //
        // The ordinal is the index i of the R-pair (rPeaks[i], rPeaks[i+1]) this
        // beat was sliced from, so it indexes the SAME R-peak vector the ECG
        // slicer is driven by. That is what makes it a shared key rather than
        // just another local index.
        std::vector<uint32_t> original_index;
        int    median_length = -1;
        size_t total_beats = 0;

        // ---- WHY THE SURVIVORS ARE FEWER THAN THE R-PAIRS -----------------
        //
        // This set is pruned before the pulse QC filter downstream ever sees
        // it, and the prunes were uncounted -- so a channel retaining 7.5% of
        // its pulses gave no way to tell whether the loss happened HERE, at a
        // fiducial that could not be found, or later at the fit-error
        // threshold. Two different fixes, and no evidence for choosing.
        //
        // Each counter names the exact `continue` that dropped the beat:
        //   n_dropped_rr    the R-pair was unusable -- too short, or over the
        //                   4 s hard cap, which is a dropout gap not a beat
        //   n_dropped_peak  no systolic upstroke peak in the search window
        //   n_dropped_foot  no trough between the window start and that peak
        //   n_dropped_up50  no clean 50% crossing on [foot, peak]
        //
        // n_slices is the denominator they are all against: R-pairs offered,
        // which is what "how much of the channel survived" has to be measured
        // over. total_beats above counts what made it through, not what was
        // tried.
        size_t n_slices = 0;
        size_t n_dropped_rr = 0;
        size_t n_dropped_peak = 0;
        size_t n_dropped_foot = 0;
        size_t n_dropped_up50 = 0;
        //   n_dropped_peak_col  the peak-column rejection below: this pulse's
        //                   systolic peak sits 5% of its own RR or further from
        //                   the median peak column of the first 100 beats.
        //
        // THIS ONE IS USUALLY THE LARGEST AND WAS THE ONE NOT COUNTED. A
        // report showing rr=0 peak=0 foot=0 up50=0 next to "kept 614 of 1223"
        // says every drop was accounted for and none of them were, which is
        // worse than no report: it sent the investigation to the fiducial
        // detectors, which turned out to be finding everything.
        size_t n_dropped_peak_col = 0;

        // The fence the peak-column rejection actually used, so the width is
        // reportable rather than inferred from how many it rejected. A wide
        // fence on a bin that still lost most of its pulses means the peak
        // timing is genuinely scattered; a narrow one means it is not, and the
        // losses are outliers.
        double peak_col_fence_lo = std::numeric_limits<double>::quiet_NaN();
        double peak_col_fence_hi = std::numeric_limits<double>::quiet_NaN();
        int    up50_aligned_col = -1;             // shared column all half-height points land on
        int    foot_aligned_col = -1;             // feet scatter (per-beat); not a shared column
        int    peak_aligned_col = -1;             // peaks scatter (per-beat); not a shared column
        int    ref_beat_index = -1;
    };

    inline PpgBeatSet extract_ppg_beats_and_align(const std::vector<double>& signal, const std::vector<size_t>& rPeaks, double fs)
    {
        PpgBeatSet out;
        const int64_t N = static_cast<int64_t>(signal.size());
        if (N == 0 || rPeaks.size() < 2) return out;

        // `r` is the beat's own R column inside `data` -- rr_before_samples of
        // its RR. Carried per beat because it IS per beat; apply_mask moves
        // the whole struct, so it cannot desynchronise from the waveform.
        struct Raw { std::vector<double> data; int peak; int foot; int up50; int r; };
        std::vector<Raw> raw;
        std::vector<int> rr_lens;
        // R-pair ordinals, parallel to `raw`. See PpgBeatSet::original_index.
        std::vector<uint32_t> raw_slice;
        raw.reserve(rPeaks.size());

        // Compact the parallel (raw, rr_lens) vectors, keeping only entries
        // where keep[i] is true.
        auto apply_mask = [&](const std::vector<bool>& keep) {
            std::vector<Raw> filt;
            std::vector<int> filt_lens;
            filt.reserve(raw.size());
            filt_lens.reserve(raw.size());
            for (size_t i = 0; i < raw.size(); ++i)
                if (keep[i]) {
                    filt.push_back(std::move(raw[i]));
                    filt_lens.push_back(rr_lens[i]);
                }
            raw.swap(filt);
            rr_lens.swap(filt_lens);
            };

        // MAX RR = 1.5 s, a hard drop. An R-pair longer than this is a
        // dropout or artifact gap between R-peaks, not a beat, and it also
        // sizes the shared window for the WHOLE bin (Pass 1: shared_w =
        // up50_anchor + max_tail) -- one such survivor NaN-pads every beat out
        // to tens of thousands of columns. Counted into n_dropped_rr with the
        // other RR rejections.
        //
        // TIGHTER THAN THE ECG SLICER, which caps at 2.5 s. The two paths are
        // joined on the R-pair ordinal but do not drop the same pairs: a 2.0 s
        // pair yields an ECG beat and no pulse beat. Deliberate or not, it is
        // the current behaviour.
        const int64_t kMaxBeatSamples = (fs > 0.0) ? static_cast<int64_t>(1.5 * fs) : 0;

        // ---- slice + per-beat peak/foot --------------------------------
        // n_slices is the R-pair count this loop was OFFERED, recorded before
        // any drop, so every counter below has a denominator.
        out.n_slices = (rPeaks.size() > 1) ? rPeaks.size() - 1 : 0;
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            const int64_t rr = static_cast<int64_t>(rPeaks[i + 1]) - r0;
            if (rr <= 3) { ++out.n_dropped_rr; continue; }
            if (kMaxBeatSamples > 0 && rr > kMaxBeatSamples) {
                ++out.n_dropped_rr; continue;
            }
            // Same splice drop as the ECG slicer, counted with the other RR
            // drops so the denominator stays honest.


            const int64_t before = rr_before_samples(rr);
            const int64_t after = rr_after_samples(rr);
            const int64_t len = before + after;
            const int64_t start = r0 - before;
            const int64_t end = r0 + after;

            std::vector<double> beat(static_cast<size_t>(len),
                std::numeric_limits<double>::quiet_NaN());
            const int64_t cs = std::max<int64_t>(0, start);
            const int64_t ce = std::min<int64_t>(N, end);
            for (int64_t k = cs; k < ce; ++k)
                beat[static_cast<size_t>(k - start)] = signal[static_cast<size_t>(k)];

            const int r_col = static_cast<int>(before);
            // Peak search is bounded to THIS beat's own R-R window: from
            // r_col+1 to the next R (which sits at r_col + rr). The full
            // beat slice extends past next R (beat length = 1.8*rr, next R
            // at 1.3*rr), so argmax-over-whole-beat would easily land on
            // the NEXT beat's peak whenever it's taller. That mislocated
            // "peak" then anchors up50 detection on the next beat's
            // upstroke, and up50-alignment then shifts every beat such
            // that individual beats' data effectively ends at their own
            // peak in the shared frame -- producing the peak-cutoff
            // plummet in the displayed template.
            const int peakSearchEnd = std::min(
                static_cast<int>(r_col + rr),
                static_cast<int>(beat.size()));
            // Upstroke-located FIRST peak, not the tallest sample in the window.
            // This is the site that matters most: the peak found here brackets
            // the foot and the up50 half-height crossing below, and a beat whose
            // up50 cannot be found is DISCARDED. With argmax, a pulse whose
            // reflected wave exceeds systole had its up50 searched on the
            // notch-to-P2 rise, which frequently has no clean single crossing --
            // so those beats were dropped, and the surviving count collapsed.
            const int peak = FeatureMarks::detect_ppg_upstroke_peak(beat, r_col + 1,
                peakSearchEnd);
            if (peak < 0) { ++out.n_dropped_peak; continue; }

            // Foot and up50 through the SAME primitives the display fiducials
            // use (FeatureMarks::trough_in / amplitude_crossing), so the
            // alignment axis and the markers drawn on it are defined
            // identically. Both were hand-rolled here, which is why a fix to
            // one never reached the other.
            const int foot = FeatureMarks::trough_in(beat, 0, peak);
            if (foot < 0) { ++out.n_dropped_foot; continue; }

            // 50%-upslope (half-height) crossing on [foot, peak]. This is the
            // horizontal alignment fiducial: it sits on the steep upstroke, so
            // its column is well-localized (unlike the flat apex or the shallow
            // foot). Interpolated first-upward-crossing, and it can FAIL --
            // both properties matter here and neither is provided by
            // amplitude_crossing, which is for display markers.
            // first_crossing returns a sub-sample crossing now; the up50 axis
            // is integer-column, so it rounds here.
            const double up50D = FeatureMarks::first_crossing(beat, foot, peak, 0.50);
            const int up50 = (up50D >= 0.0)
                ? static_cast<int>(std::lround(up50D)) : -1;
            if (up50 < 0) { ++out.n_dropped_up50; continue; }   // no upslope

            // NOTE: peak/foot are stored raw (no subsample_refine call) --
            // they only feed the peak-position rejection check below and
            // the (unread by any caller) peak_cols/foot_cols output,
            // neither of which needs sub-sample precision. The expensive
            // 4x-upsample fit-and-select refinement was pure overhead here;
            // the real fiducials used downstream (outPeakCol/outFootCol)
            // are recomputed independently from the final median template.
            raw.push_back({ std::move(beat), peak, foot, up50, r_col });
            rr_lens.push_back(static_cast<int>(rr));
            // Pushed HERE, in the same statement group as the beat itself, and
            // after every `continue` above. A beat and its ordinal have to be
            // appended together or the two vectors silently desynchronise at
            // the first dropped beat -- and the result still looks like a valid
            // parallel pair.
            raw_slice.push_back(static_cast<uint32_t>(i));
        }
        if (raw.empty()) return out;

        // ---- representative (median) length ----------------------------
        // Middle element of the sorted lengths (a length some beat has).
        {
            std::vector<int> lens = rr_lens;
            std::sort(lens.begin(), lens.end());
            out.median_length = lens[lens.size() / 2];
        }
        for (size_t i = 0; i < raw.size(); ++i) {
            if (rr_lens[i] == out.median_length) {
                out.ref_beat_index = static_cast<int>(i);
                break;
            }
        }

        // ---- Pass 1: 50%-upslope align on a shared axis (shift + NaN) ---
        // Anchor every beat's half-height point to the MEDIAN up50 column
        // among survivors. Beats with up50 > anchor get their leading samples
        // clipped (dst < 0 dropped by the guard below); this is accepted.
        std::vector<int> up50s;
        up50s.reserve(raw.size());
        for (const auto& r : raw) up50s.push_back(r.up50);
        std::sort(up50s.begin(), up50s.end());
        const int up50_anchor = up50s[up50s.size() / 2];   // median

        int max_tail = 0;         // max (beat_len - up50) over survivors
        for (const auto& r : raw) {
            const int tail = static_cast<int>(r.data.size()) - r.up50;
            if (tail > max_tail) max_tail = tail;
        }
        const int shared_w = up50_anchor + max_tail;
        if (shared_w <= 0) return out;
        const double NaND = std::numeric_limits<double>::quiet_NaN();

        out.beats.reserve(raw.size());
        out.peak_cols.reserve(raw.size());
        out.foot_cols.reserve(raw.size());
        out.r_cols.reserve(raw.size());
        out.original_index.reserve(raw.size());
        for (size_t ri = 0; ri < raw.size(); ++ri) {
            const auto& b = raw[ri];
            if (ri < raw_slice.size()) out.original_index.push_back(raw_slice[ri]);
            const int prepend = up50_anchor - b.up50;   // may be < 0 now
            std::vector<double> a(shared_w, NaND);
            for (int k = 0; k < (int)b.data.size(); ++k) {
                const int dst = prepend + k;
                if (dst >= 0 && dst < shared_w) a[dst] = b.data[k];   // clips left overflow
            }
            out.beats.push_back(std::move(a));
            out.peak_cols.push_back(prepend + b.peak);   // may be < 0 for clipped beats
            out.foot_cols.push_back(prepend + b.foot);   // may be < 0 for clipped beats
            // The same shift the waveform got, applied to the R column: this
            // is the one number that ties the pulse frame back to the QRS.
            out.r_cols.push_back(prepend + b.r);
        }
        out.up50_aligned_col = up50_anchor;


        out.total_beats = out.beats.size();
        if (out.beats.empty()) return out;

        // ---- Pass 2: min-baseline vertical DC match ---------------------
        // Match each beat's baseline (windowed mean around its OWN minimum
        // sample) to the reference beat's. We use the per-beat argmin rather
        // than the stored foot column: after median-anchoring, a beat's foot
        // can be clipped off the left edge (foot_cols[i] < 0), but the argmin
        // over surviving samples is always a valid in-range column. For
        // unclipped beats the min IS the foot; for clipped beats it's the
        // lowest surviving point, a good-enough baseline proxy. Windowed mean
        // (not the single argmin sample) avoids order-statistic bias. A
        // constant vertical shift can't disturb horizontal alignment.
        if (out.ref_beat_index >= 0
            && out.ref_beat_index < static_cast<int>(out.beats.size()))
        {
            const int fb_w = std::max(1, out.median_length / 50);
            auto min_baseline = [&](size_t i) -> double {
                const auto& beat = out.beats[i];
                int mc = -1;
                double mv = std::numeric_limits<double>::infinity();
                for (int k = 0; k < (int)beat.size(); ++k) {
                    const double v = beat[k];
                    if (!std::isnan(v) && v < mv) { mv = v; mc = k; }
                }
                if (mc < 0) return std::numeric_limits<double>::quiet_NaN();
                const int lo = std::max(0, mc - fb_w);
                const int hi = std::min(shared_w, mc + fb_w + 1);
                double sum = 0.0; int n = 0;
                for (int k = lo; k < hi; ++k) {
                    const double v = beat[k];
                    if (!std::isnan(v)) { sum += v; ++n; }
                }
                return n >= 1 ? sum / n
                    : std::numeric_limits<double>::quiet_NaN();
                };

            const double target = min_baseline(static_cast<size_t>(out.ref_beat_index));
            if (!std::isnan(target)) {
                for (size_t i = 0; i < out.beats.size(); ++i) {
                    const double b_base = min_baseline(i);
                    if (std::isnan(b_base)) continue;
                    const double d = b_base - target;
                    if (d == 0.0) continue;
                    for (double& v : out.beats[i])
                        if (!std::isnan(v)) v -= d;
                }
            }
        }

        return out;
    }

}   // namespace alignment