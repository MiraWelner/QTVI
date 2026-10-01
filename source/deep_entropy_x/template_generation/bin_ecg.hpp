/**
 * @file   bin_ecg.hpp
 * @brief  Create the temporally defined bins from the median of the ecg beats over a temporal range.
 *         The bins are not split based on morpohlogy, that comes later in the templating
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-27
 */
#pragma once

#include "fiducial_marker_finding/alignment.hpp"
#include "template_generation/seed_pool.hpp"
#include "template_generation/template_io.hpp"
#include <chrono>
#include <cstdio>
#include <atomic>
#include <deque>
#include <fstream>
#include <string>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <array>
#include <vector>
#include <cstdint>


struct EcgChannelResult {
    std::vector<std::vector<double>> ecgTemplates_raw;
    std::vector<std::vector<double>> ecgTemplates_raw_std;   // parallel to ecgTemplates_raw
    std::vector<std::vector<double>> ecgTemplates_squared;
    std::vector<std::vector<double>> ecgTemplates_absval;
    std::vector<std::vector<double>> ecgTemplates_unfiltered;
    std::vector<int> ref_index_raw;

    std::vector<int> r_col_raw;
    std::vector<int> r_col_squared;
    std::vector<int> r_col_absval;
    std::vector<int> r_col_unfiltered;

    std::vector<size_t> n_beats_raw;//the viewer displays the number of beats contributing to template for each channel

    std::vector<std::vector<std::vector<double>>> kept_beats_raw;
    std::vector<std::vector<uint8_t>> kept_rhythm_raw;
    std::vector<uint8_t> seed_basis_raw;
    std::vector<std::vector<double>> tp_shift_raw;
    std::vector<std::vector<double>> pq_shift_raw;
    std::vector<std::vector<int>> kept_of_aligned_raw;   // for the move log
};

struct EcgTemplateResult {
    EcgChannelResult ch1;
    EcgChannelResult ch2;
    EcgChannelResult ch3;
    std::array<std::vector<std::vector<size_t>>, 3> kept_index;
};
// ---- THE PER-BEAT MOVE LOG: <stem>_beat_moves.csv ---------------------------
//
// One row per ECG beat the raw method sliced, per channel per bin, with every
// move the build made to it:
//
//   tp_mv_shift, pq_mv_shift   VERTICAL. The two-stage leveling in
//                              extract_beats_and_align: the amount SUBTRACTED
//                              from the beat by the TP pass and by the PQ pass.
//                              NaN where that pass did not apply.
//   kept_row                   The beat's row in per_channel_beats[ch][bin] --
//                              the local row space bank members index. -1 for
//                              a beat the raw method did not keep (baseline
//                              NONE), which therefore took part in no anchor
//                              alignment.
//   <A>_ms_shift               HORIZONTAL, one column per anchor alignment
//                              (r, p, q, j): the sub-sample cross-correlation
//                              shift align_beat_matrix applied to this beat to
//                              put its own landmark on the template's, in ms.
//                              POSITIVE MOVES THE BEAT LATER. NaN where the
//                              beat got no estimate (below the correlation
//                              floor, or the landmark was not found on the
//                              template so nothing in that bin moved).
//
// TWO STAGES, ONE WRITE. The vertical shifts exist at the end of the fast
// build (CreateEcgTemplatesFast) and the horizontal ones only after
// analysis_job::prepare has run every anchor (alignTemplatesFromCache), so the
// build STASHES its half here and prepare calls write() once both are in.
// set() clears both halves, so a record can never be written with the previous
// record's numbers.
//
// Destination set once per record from prepare before the build. Empty
// dir/stem => no log. stash_vertical and horizontal_bins are called
// single-threaded; the pointer horizontal_bins returns is written from the
// parallel bin loop, one (bin, channel) element per iteration, which is
// race-free because it is sized before the loop starts.
namespace ecg_move_log {
    inline std::string g_dir;
    inline std::string g_stem;

    // [bin] -> this channel's per-aligned-beat vertical shifts and the kept
    // row each aligned beat became.
    struct VerticalBin {
        std::vector<double> tp;
        std::vector<double> pq;
        std::vector<int>    kept_row;
    };
    inline std::array<std::vector<VerticalBin>, 3> g_vertical;

    // One anchor's shifts, [bin][channel][kept_row], in samples. A deque so a
    // pointer handed out for one anchor survives the next anchor being added.
    struct HorizontalAnchor {
        std::string label;
        std::vector<std::array<std::vector<double>, 3>> bins;
    };
    inline std::deque<HorizontalAnchor> g_horizontal;

    inline void reset() {
        for (auto& v : g_vertical) v.clear();
        g_horizontal.clear();
    }
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem; reset();
    }

    // The fast build's half. Per bin: tp/pq parallel to the aligned beats,
    // and kept_row mapping each aligned beat to its kept row (or -1).
    inline void stash_vertical(int chIdx,
        std::vector<std::vector<double>> tp,
        std::vector<std::vector<double>> pq,
        std::vector<std::vector<int>> kept_row)
    {
        if (chIdx < 0 || chIdx >= 3) return;
        const size_t n = tp.size();
        std::vector<VerticalBin>& dst = g_vertical[chIdx];
        dst.assign(n, {});
        for (size_t b = 0; b < n; ++b) {
            dst[b].tp = std::move(tp[b]);
            if (b < pq.size())       dst[b].pq = std::move(pq[b]);
            if (b < kept_row.size()) dst[b].kept_row = std::move(kept_row[b]);
        }
    }

    // One anchor's store, sized to nBins. Call once per anchor BEFORE the
    // parallel bin loop; write (*out)[bin][chIdx] = shifts inside it.
    inline std::vector<std::array<std::vector<double>, 3>>*
        horizontal_bins(const std::string& label, size_t nBins)
    {
        for (HorizontalAnchor& h : g_horizontal)
            if (h.label == label) {
                h.bins.assign(nBins, {});
                return &h.bins;
            }
        g_horizontal.push_back(HorizontalAnchor{ label, {} });
        g_horizontal.back().bins.assign(nBins, {});
        return &g_horizontal.back().bins;
    }

    // Both halves, one file. fs is the ECG rate the horizontal shifts are
    // converted to ms with. Clears the store afterwards.
    inline void write(double fs) {
        if (g_dir.empty() || g_stem.empty()) { reset(); return; }
        std::ofstream f(g_dir + "/" + g_stem + "_ecg_alignment_shifts.csv", std::ios::trunc);
        if (!f) { reset(); return; }
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        const double msPerSample = (fs > 0.0) ? 1000.0 / fs : kNaN;

        f << "stem,channel,bin,beat,tp_mv_shift,pq_mv_shift,kept_row";
        for (const HorizontalAnchor& h : g_horizontal) f << ',' << h.label << "_ms_shift";
        f << '\n';

        static const char* kChan[3] = { "CH1", "CH2", "CH3" };
        for (int c = 0; c < 3; ++c) {
            const std::vector<VerticalBin>& bins = g_vertical[c];
            for (size_t b = 0; b < bins.size(); ++b) {
                const VerticalBin& vb = bins[b];
                for (size_t k = 0; k < vb.tp.size(); ++k) {
                    const double pqv = (k < vb.pq.size()) ? vb.pq[k] : kNaN;
                    const int kr = (k < vb.kept_row.size()) ? vb.kept_row[k] : -1;
                    f << g_stem << ',' << kChan[c] << ',' << b << ',' << k
                        << ',' << vb.tp[k] << ',' << pqv << ',' << kr;
                    for (const HorizontalAnchor& h : g_horizontal) {
                        double v = kNaN;
                        if (kr >= 0 && b < h.bins.size()) {
                            const std::vector<double>& row = h.bins[b][c];
                            if (static_cast<size_t>(kr) < row.size())
                                v = row[static_cast<size_t>(kr)] * msPerSample;
                        }
                        f << ',' << v;
                    }
                    f << '\n';
                }
            }
        }
        reset();
    }
}

struct SingleMethodResult {
    vector<double> ecgTemplate;
    vector<double> ecg_template_std;   // empty for methods that don't compute std
    int r_col = -1;   // true R column (alignment's r_aligned_col)
    int median_rr_samples = -1;
    // Verdict per beat handed downstream, parallel to out_kept_beats:
    //   0 NORMAL     1 PVC (premature)     2 VOTED_PVC (5-of-8 vote)
    std::vector<uint8_t> kept_rhythm;

    // How the Phase 1 reference pool was chosen, and what it cost: the counts
    // excluded for prematurity and for the vote, and the ectopic fraction of the
    // candidates BEFORE the fallback ladder. A bin whose basis is not
    // SINUS_ONLY has a reference that is not purely sinus, and that has to be
    // visible rather than inferred later.
    seed_pool::SeedSelection seed;
    int ref_beat_index = -1;
    size_t n_beats = 0;
    // Per-beat per-stage vertical DC shifts (TP stage, PQ stage) from the
    // two-stage leveling, surfaced for the move log.
    vector<double> tp_shift;
    vector<double> pq_shift;
    // Parallel to tp_shift / pq_shift (the ALIGNED beats): the kept row each
    // aligned beat became, or -1 for one excluded from `usable`. The move log
    // joins the anchor alignments' per-row shifts through it.
    std::vector<int> kept_of_aligned;



    // THE JOIN KEY: kept_idx[k] is the R-PAIR ORDINAL that produced the beat at
    // slot k of out_kept_beats. It is alignment::ecg_beat_set::slice_index,
    // carried through this channel's own pruning.
    //
    // WHY AN ORDINAL AND NOT A ROW. Two separate reasons, and conflating them
    // is how this comment used to mislead:
    //
    //   (1) A row is not an ordinal on any channel. The slicing loop skips
    //       R-pairs (rr <= 3 samples, and rr > 2.5 s, which is a dropout or a
    //       splice rather than a beat) with `continue` before anything is
    //       pushed, so row k is the k-th survivor, not pair k. Anything that
    //       looks a beat up in R-peak space -- RR intervals, mark codes, the
    //       pulse channel's slice list -- needs the ordinal.
    //
    //   (2) Rows are not comparable ACROSS channels. Every channel is sliced
    //       from the same ch1.raw peaks and so skips the same pairs, but each
    //       then prunes independently (shape QC, baseline_source == NONE, and
    //       the pulse path's own outlier rules at its own rate). So slot k of
    //       CH1 and slot k of PPG are different heartbeats. The ordinal is a
    //       shared identity because the R-peak vector is shared; the row is
    //       local bookkeeping.
    //
    // The row has no consumer left: it existed for the morphology writers,
    // whose columns are now slices, and they resolve a waveform through jbank's
    // slice -> row map instead. One map, one meaning.
    std::vector<size_t> kept_idx;
};

static inline SingleMethodResult build_ecg_template_for_method(const vector<double>& ecgSignal, const vector<size_t>& rpeaks,
    double ecgRate, vector<vector<double>>* out_kept_beats = nullptr, bool compute_std = false) {
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ecg_template_std = {};
    res.r_col = -1;

    if (rpeaks.size() < 2 || ecgSignal.empty() || ecgRate <= 0.0) return res;

    const alignment::ecg_beat_set aligned =
        alignment::extract_beats_and_align(ecgSignal, rpeaks, ecgRate);
    if (aligned.beats.empty() || aligned.median_length <= 0) return res;

    res.tp_shift = aligned.tp_shift;   // surface for the move log
    res.pq_shift = aligned.pq_shift;

    // Beats with baseline_source == NONE had neither a usable TP nor PQ
    // isoelectric reference, so their DC level is untrustworthy -- exclude
    // them from every amplitude-dependent aggregate below (median template,
    // std band, and the surviving-beats QC capture), same as the Tukey/
    // wave-score rejections that already ran upstream in extract_beats_and_
    // align(). NOTE: this changes prior behavior -- previously an
    // unavailable baseline meant "use the beat un-shifted anyway"; now it
    // means "exclude it entirely" per spec. If baseline_source is empty or
    // mismatched in length (e.g. ref_beat_index was invalid so Pass 3 never
    // ran), fall back to using every beat unfiltered rather than silently
    // producing an empty template.
    const bool haveSrc = aligned.baseline_source.size() == aligned.beats.size();
    std::vector<const std::vector<double>*> usable;
    std::vector<size_t> usableIdx;          // parallel to `usable`, into aligned.*
    usable.reserve(aligned.beats.size());
    usableIdx.reserve(aligned.beats.size());
    for (size_t i = 0; i < aligned.beats.size(); ++i) {
        if (haveSrc && aligned.baseline_source[i] == alignment::BaselineSource::NONE) continue;
        usable.push_back(&aligned.beats[i]);
        usableIdx.push_back(i);
    }
    if (usable.empty()) {   // every beat's baseline was NONE -- fail safe, don't zero the template
        usable.reserve(aligned.beats.size());
        usableIdx.reserve(aligned.beats.size());
        for (size_t i = 0; i < aligned.beats.size(); ++i) {
            usable.push_back(&aligned.beats[i]);
            usableIdx.push_back(i);
        }
    }

    const size_t maxLen = usable.front()->size();   // shared-axis width

    // R column: the detected-R fiducial the template was built around, straight
    // from alignment (every beat's detected R lands at r_aligned_col). Passed
    // through as-is -- no re-detection (a window search would grab Q or S).
    res.r_col = aligned.r_aligned_col;
    res.median_rr_samples = aligned.median_length;   // display width; see struct
    res.ref_beat_index = aligned.ref_beat_index;

    // Column-wise NaN-skipping median over the aligned beats => template.
    auto medianOver = [&](const std::vector<const std::vector<double>*>& set) {
        std::vector<double> tmpl(maxLen, NaN);
        for (size_t c = 0; c < maxLen; ++c) {
            std::vector<double> col;
            col.reserve(set.size());
            for (const auto* sl : set) {
                const double v = (*sl)[c];
                if (!std::isnan(v)) col.push_back(v);
            }
            if (col.empty()) continue;
            // nth_element, NOT sort. This is the hottest loop in the template
            // build: once per column, per method, per channel, per bin -- for a
            // 1.8*RR axis that is ~1800 columns times two fast methods times
            // three channels times every bin. A full sort is O(n log n) to
            // extract ONE order statistic; partial selection gets it in O(n).
            // Measured 17.4 ms -> 6.0 ms per pass at 1800 columns x 300 beats,
            // bit-identical output.
            //
            // The even case also needs the element below the midpoint, and
            // nth_element has already partitioned everything below imid to its
            // left -- so max_element over that prefix finds it with no second
            // selection.
            const size_t nc = col.size();
            const size_t imid = nc / 2;
            std::nth_element(col.begin(), col.begin() + imid, col.end());
            const double hi_mid = col[imid];
            tmpl[c] = (nc % 2 == 0)
                ? 0.5 * (*std::max_element(col.begin(), col.begin() + imid) + hi_mid)
                : hi_mid;
        }
        return tmpl;
        };

    // ---- THE ECTOPIC MASK, WHICH HAD NO CALLER UNTIL NOW -----------------
    //
    // seed_pool.hpp exists to select the beats allowed to form the reference,
    // and nothing called it. The median above was over `usable`, filtered on
    // exactly one condition -- baseline_source != NONE -- with no rhythm test,
    // while alignment.hpp EXEMPTS flagged beats from its pruning specifically so
    // that "the ectopic mask (create_ecg_templates.hpp)" could exclude them
    // here. The net effect was the opposite of the intent: the flags rescued
    // ectopic beats from RR-length pruning and then nothing kept them out of the
    // reference.
    //
    // "Exclude PVCs and artifact from reference calculations while retaining
    // them with flags" is the 4.5 clause, and this is the reference. The beats
    // are all still captured, still partitioned, still written.
    //
    // THIS IS ALSO WHERE THE ORDER RULE PUTS IT. The rhythm flags may gate what
    // the Phase 1 REFERENCE is built from; they may not gate the partition. The
    // bank is seeded with this waveform and then scores every beat against it
    // rhythm-blind, so a PVC still gets compared, still fails 0.85, and still
    // opens its own template -- which it cannot do if the thing it is compared
    // against is half PVC.
    //
    // ASYMMETRIC COSTS, so the gate is strict. A misclassified beat is one wrong
    // number. A contaminated reference damages the template, the corridor built
    // from the same pool, and every feature in the bin, and it compounds: a
    // wider corridor admits the next ectopic beat more easily.
    std::vector<uint8_t> rhythm_of_slot;
    {
        const bool haveFlags = aligned.premature.size() == aligned.beats.size()
            && aligned.voted.size() == aligned.beats.size();
        rhythm_of_slot.assign(usableIdx.size(), 0);
        if (haveFlags)
            for (size_t k = 0; k < usableIdx.size(); ++k) {
                const size_t ai = usableIdx[k];
                rhythm_of_slot[k] = aligned.premature[ai] ? 1u
                    : (aligned.voted[ai] ? 2u : 0u);
            }
    }
    std::vector<uint32_t> slotIdx(usableIdx.size());
    for (uint32_t k = 0; k < slotIdx.size(); ++k) slotIdx[k] = k;

    // No operator marks exist at build time, so `category` is left empty and
    // every beat reads REGULAR: the selection rests on the rhythm flags alone.
    // That is the design working -- morphology does the sorting, marks only
    // supply labels later.
    const seed_pool::SeedSelection sel =
        seed_pool::selectSeedPool(slotIdx, rhythm_of_slot, {});
    res.seed = sel;

    std::vector<const std::vector<double>*> reference;
    reference.reserve(sel.members.size());
    for (const uint32_t k : sel.members)
        if (k < usable.size()) reference.push_back(usable[k]);

    // selectSeedPool never returns an empty pool when candidates exist, and its
    // fallback ladder LABELS a contaminated pool rather than hiding it -- so a
    // bin where ectopy is the majority still gets a reference, and
    // SeedSelection::basis says it is not a clean one.
    res.ecgTemplate = medianOver(reference.empty() ? usable : reference);
    if (out_kept_beats) {
        out_kept_beats->clear();
        out_kept_beats->reserve(usable.size());
        for (const auto* sl : usable) out_kept_beats->push_back(*sl);
    }
    res.kept_of_aligned.assign(aligned.beats.size(), -1);
    for (size_t k = 0; k < usableIdx.size(); ++k)
        if (usableIdx[k] < res.kept_of_aligned.size())
            res.kept_of_aligned[usableIdx[k]] = static_cast<int>(k);
    res.kept_idx.resize(usableIdx.size());
    for (size_t k = 0; k < usableIdx.size(); ++k) {
        const size_t ai = usableIdx[k];
        res.kept_idx[k] = (ai < aligned.slice_index.size())
            ? static_cast<size_t>(aligned.slice_index[ai]) : ai;
    }
    {
        const bool haveFlags = aligned.premature.size() == aligned.beats.size()
            && aligned.voted.size() == aligned.beats.size();
        res.kept_rhythm.assign(usableIdx.size(), 0);
        if (haveFlags) {
            for (size_t k = 0; k < usableIdx.size(); ++k) {
                const size_t ai = usableIdx[k];
                // Premature wins over voted: direct evidence over inferred.
                res.kept_rhythm[k] = aligned.premature[ai] ? 1u
                    : (aligned.voted[ai] ? 2u : 0u);
            }
        }
    }

    res.n_beats = usable.size();

    if (compute_std) {//compute std is only defined for the raw templates which are displayed
        const std::vector<const std::vector<double>*>& spreadSet = reference.empty() ? usable : reference;
        res.ecg_template_std.assign(maxLen, 0.0);
        std::vector<double> col;
        col.reserve(spreadSet.size());
        for (size_t c = 0; c < maxLen; ++c) {
            col.clear();
            for (const auto* sl : spreadSet)
                if (!std::isnan((*sl)[c])) col.push_back((*sl)[c]);
            const size_t nc = col.size();
            if (nc < 2) continue;
            double mean = 0.0;
            for (double v : col) mean += v;
            mean /= static_cast<double>(nc);
            double sumsq = 0.0;
            for (double v : col) sumsq += (v - mean) * (v - mean);
            res.ecg_template_std[c] = std::sqrt(sumsq / static_cast<double>(nc - 1));   // ddof = 1
        }
    }
    return res;
}

static inline void init_channel_result(EcgChannelResult& cr, size_t n) {
    cr.ecgTemplates_raw.resize(n);
    cr.ecgTemplates_raw_std.resize(n);
    cr.ecgTemplates_squared.resize(n);
    cr.ecgTemplates_absval.resize(n);
    cr.ecgTemplates_unfiltered.resize(n);
    cr.ref_index_raw.resize(n, -1);

    cr.r_col_raw.resize(n, -1);
    cr.r_col_squared.resize(n, -1);
    cr.r_col_absval.resize(n, -1);
    cr.r_col_unfiltered.resize(n, -1);

    cr.n_beats_raw.resize(n, 0);

    cr.kept_beats_raw.resize(n);
    cr.kept_rhythm_raw.resize(n);
    cr.seed_basis_raw.assign(n, static_cast<uint8_t>(seed_pool::SeedBasis::EMPTY));
    cr.tp_shift_raw.resize(n);
    cr.pq_shift_raw.resize(n);
    cr.kept_of_aligned_raw.resize(n);
}

/**
 * @brief  Process all 4 methods for one channel in one bin.
 *
 *         Only the raw method computes std (the other three are never
 *         displayed in the viewer). Only ch1 captures the surviving raw
 *         beats for QC output.
 *
 * @param cr             Channel result accumulator
 * @param bins           All bins (for pairs access)
 * @param i              Current bin index
 * @param ecgSignal      The signal used for raw/sq/abs methods (may be preprocessed)
 * @param origSignal     The original unfiltered ECG signal for this channel
 * @param ch             Channel R-peaks struct
 * @param capture_raw_beats  If true, capture the surviving aligned beats
 *                           from the "raw" method into cr.kept_beats_raw[i].
 */
 // FAST methods: raw (the displayed one, with std) + unfiltered. These are
// everything the viewer renders. Captures the ch1 raw beats for QC.
//
// Patch A change: ecgRate is threaded in for pair-window slicing, and
// `masterPeaks` (== bin.ch1.raw) drives the slicing for every channel so
// every channel's template covers the same real-time window.
static inline void process_channel_fast(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const vector<double>& origSignal,
    const vector<size_t>& masterPeaks,
    double ecgRate,
    bool capture_raw_beats = false,
    int channel_index = 0,
    // ALIGNED -> CAPTURED SLOT for this bin, or null when the caller does not
    // need it. An OUT-PARAM rather than a field on EcgChannelResult: the map is
    // consumed by the morphology writers in this function's caller and nowhere
    // else, so a member would widen a type in another header for one local use.
    // CAPTURED SLOT -> R-PAIR SLICE for this bin, or null when the caller does
    // not need it.
    std::vector<size_t>* out_kept_idx = nullptr)
{
    const auto& bin = bins[i];

    // Method 1: raw (detection signal + master R-peaks). Only method with std.
    vector<vector<double>>* capture =
        (capture_raw_beats && i < cr.kept_beats_raw.size())
        ? &cr.kept_beats_raw[i] : nullptr;
    auto raw_res = build_ecg_template_for_method(
        ecgSignal, masterPeaks, ecgRate,
        capture, /*compute_iqr=*/true);
    if (out_kept_idx) *out_kept_idx = std::move(raw_res.kept_idx);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ecgTemplates_raw_std[i] = raw_res.ecg_template_std;
    cr.r_col_raw[i] = raw_res.r_col;
    cr.n_beats_raw[i] = raw_res.n_beats;
    if (i < cr.kept_rhythm_raw.size())
        cr.kept_rhythm_raw[i] = std::move(raw_res.kept_rhythm);
    if (i < cr.seed_basis_raw.size())
        cr.seed_basis_raw[i] = static_cast<uint8_t>(raw_res.seed.basis);
    cr.ref_index_raw[i] = raw_res.ref_beat_index;
    if (i < cr.tp_shift_raw.size()) {
        cr.tp_shift_raw[i] = std::move(raw_res.tp_shift);   // distinct i -> race-free
        cr.pq_shift_raw[i] = std::move(raw_res.pq_shift);
        cr.kept_of_aligned_raw[i] = std::move(raw_res.kept_of_aligned);
    }

    // Method 4: unfiltered (original ECG signal + master R-peaks). No std.
    auto unfilt_res = build_ecg_template_for_method(
        origSignal, masterPeaks, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.r_col_unfiltered[i] = unfilt_res.r_col;
}

// SLOW methods: squared + absval. Not displayed by the viewer; safe to
// compute off the critical path. Writes only the squared/absval fields of
// cr (which process_channel_fast leaves untouched).
//
// Patch A change: same as fast pass, master ch1 peaks drive slicing. The
// per-channel preprocessed signals are still used (squared/absval variants
// have different amplitudes), but they're indexed at the master R-peak
// positions since R sample indices map 1:1 across ECG preprocessing.
static inline void process_channel_slow(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const ChannelRPeaks& ch,
    const vector<size_t>& masterPeaks,
    double ecgRate)
{
    const auto& bin = bins[i];


    // Method 2: squared (squared signal + master R-peaks). No std.
    const auto& sq_sig = ch.squared_signal.empty() ? ecgSignal : ch.squared_signal;
    auto sq_res = build_ecg_template_for_method(
        sq_sig, masterPeaks, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.r_col_squared[i] = sq_res.r_col;

    // Method 3: absval (abs-value signal + master R-peaks). No std.
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(abs_sig, masterPeaks, ecgRate, nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.r_col_absval[i] = abs_res.r_col;
}


// FAST pass: raw + unfiltered templates (everything the viewer needs).
// Leaves the squared/absval vectors sized-but-empty for CreateEcgTemplatesSlow.
inline EcgTemplateResult CreateEcgTemplatesFast(
    const vector<output_binfile_data>& bins,
    double ecgRate)
{
    size_t n = bins.size();
    EcgTemplateResult res;
    init_channel_result(res.ch1, n);
    init_channel_result(res.ch2, n);
    init_channel_result(res.ch3, n);

    // ALIGNED -> CAPTURED SLOT, [channel][bin][slot]. Local, because the only
    // consumer is the morphology write in this function's post-loop slot.
    // Pre-sized so the parallel loop only ever writes its own element.
    // Per channel per bin: capturedSlot -> R-pair slice. The join key the
    // Section 4.6 partition is built on.
    std::array<std::vector<std::vector<size_t>>, 3> keptIdx;
    for (auto& k : keptIdx) k.resize(n);

    // ---- PER-BIN PROGRESS, deliberately NOT one-shot --------------------
    // alignment.hpp's [ALIGN-4S-CAP-ACTIVE] used a process-lifetime
    // `static bool printedOnce`: it fires on the first file of a run and never
    // again, so on file 2..N its ABSENCE says nothing about how far that file
    // got. Diagnosing from that absence is how the previous attempt concluded
    // the wrong stage. Counted, atomic and flushed, so a stall names its bin.
    std::atomic<int> _done{ 0 };

    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];
        const auto& master = bin.ch1.raw;   // ch1.raw drives every channel's slicing
        process_channel_fast(res.ch1, bins, i, bin.ecgSignal, bin.ecgSignal,
            master, ecgRate, /*capture_raw_beats=*/true, /*channel_index=*/0,
            &keptIdx[0][i]);
        // Only build ch2/ch3 templates when the channel is REAL: both the
        // signal is present AND R-peak detection actually found something.
        // file_to_bin fills absent channels with placeholder vectors (see
        // "0.0 = channel absent" in file_to_bin.hpp), so `.empty()` alone
        // isn't a reliable "channel exists" signal -- but a truly absent
        // channel will never produce R-peaks (zero-variance placeholder
        // trips run_rpeak_detection's std_dev==0 noisy flag). Requiring
        // bin.ch2.raw non-empty catches those.
        if (!bin.ecgSignal2.empty() && !bin.ch2.raw.empty())
            process_channel_fast(res.ch2, bins, i, bin.ecgSignal2, bin.ecgSignal2,
                master, ecgRate, /*capture_raw_beats=*/true, /*channel_index=*/1,
                &keptIdx[1][i]);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_fast(res.ch3, bins, i, bin.ecgSignal3, bin.ecgSignal3,
                master, ecgRate, /*capture_raw_beats=*/true, /*channel_index=*/2,
                &keptIdx[2][i]);
    }

    // STASHED, NOT WRITTEN: the horizontal half of the move log does not exist
    // until prepare has run every anchor. See ecg_move_log. Copied, because
    // res is returned with these fields and nothing says no one reads them.
    ecg_move_log::stash_vertical(0, res.ch1.tp_shift_raw, res.ch1.pq_shift_raw, res.ch1.kept_of_aligned_raw);
    ecg_move_log::stash_vertical(1, res.ch2.tp_shift_raw, res.ch2.pq_shift_raw, res.ch2.kept_of_aligned_raw);
    ecg_move_log::stash_vertical(2, res.ch3.tp_shift_raw, res.ch3.pq_shift_raw, res.ch3.kept_of_aligned_raw);

    // The join key, surfaced so the partition and the archive read the SAME map
    // rather than two copies that can drift. Moved, not copied: keptIdx dies
    // with this function otherwise.
    for (int c = 0; c < 3; ++c) res.kept_index[c] = keptIdx[c];
    return res;
}

// SLOW pass: fills the squared/absval templates onto an EcgTemplateResult
// that has already been sized (e.g. by CreateEcgTemplatesFast, or by
// init_channel_result). Touches only squared/absval fields.
inline void CreateEcgTemplatesSlow(
    const vector<output_binfile_data>& bins,
    double ecgRate,
    EcgTemplateResult& res)
{
    size_t n = bins.size();
    int max_threads = std::min(8, (int)n);
#pragma omp parallel for schedule(dynamic) num_threads(max_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& bin = bins[i];
        const auto& master = bin.ch1.raw;
        process_channel_slow(res.ch1, bins, i, bin.ecgSignal, bin.ch1, master,
            ecgRate);
        // Same "real channel" gate as the fast pass (see comment there).
        if (!bin.ecgSignal2.empty() && !bin.ch2.raw.empty())
            process_channel_slow(res.ch2, bins, i, bin.ecgSignal2, bin.ch2, master,
                ecgRate);
        if (!bin.ecgSignal3.empty() && !bin.ch3.raw.empty())
            process_channel_slow(res.ch3, bins, i, bin.ecgSignal3, bin.ch3, master,
                ecgRate);
    }
}