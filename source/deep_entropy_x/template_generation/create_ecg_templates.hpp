/**
 * @file   CreateEcgTemplates.hpp
 * @brief  Create ECG templates for each bin using EnsembleTemplate.
 *         Builds templates from 3 channels x 4 preprocessing methods
 *         (raw, squared, absval, unfiltered).
 *
 *         The "unfiltered" method uses the original ECG signal (ecgSignal,
 *         ecgSignal2, ecgSignal3) with the raw R-peaks to build a template
 *         from the signal before any preprocessing or filtering.
 *
 *         For the "raw" method we also capture two extra outputs:
 *           - The surviving aligned beats that contributed to the template
 *             (ch1 only -- downstream code writes these out for QC).
 *           - The per-sample std across those beats (all channels), which
 *             the viewer draws as a gray band around the displayed
 *             template. The other three methods (squared/absval/unfiltered)
 *             don't get std computed since the viewer never displays them.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-27
 */
#pragma once

#include "template_structs.hpp"
#include "template_marking_gui/alignment.hpp"
#include "template_morphology_grouping/bin_pipeline.hpp"
// selectSeedPool: the ectopic mask on the Phase 1 reference.
#include "template_morphology_grouping/seed_pool.hpp"
#include "template_morphology_grouping/morphology_csv.hpp"
#include <chrono>
#include <cstdio>
#include <atomic>
#include <fstream>
#include <string>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

 // Beat-move log destination. Set once from main/post_process before the build
 // (read-only afterwards, and the writer runs single-threaded post-loop, so no
 // race). Empty dir/stem => no log.
namespace ecg_move_log {
    inline std::string g_dir;
    inline std::string g_stem;
    inline void set(const std::string& dir, const std::string& stem) { g_dir = dir; g_stem = stem; }

    // Write one channel's per-bin/per-beat two-stage vertical shifts. Called
    // single-threaded (after the parallel build loop). `first` truncates +
    // writes the header; later channels append.
    inline void write_channel(const char* channel,
        const std::vector<std::vector<double>>& tp,
        const std::vector<std::vector<double>>& pq, bool first) {
        if (g_dir.empty() || g_stem.empty()) return;
        std::ofstream f(g_dir + "/" + g_stem + "_beat_moves.csv",
            first ? std::ios::trunc : std::ios::app);
        if (!f) return;
        if (first) f << "stem,channel,bin,beat,tp_mv_shift,pq_mv_shift\n";
        for (size_t b = 0; b < tp.size(); ++b)
            for (size_t k = 0; k < tp[b].size(); ++k) {
                const double pqv = (b < pq.size() && k < pq[b].size())
                    ? pq[b][k] : std::numeric_limits<double>::quiet_NaN();
                f << g_stem << ',' << channel << ',' << b << ',' << k
                    << ',' << tp[b][k] << ',' << pqv << '\n';
            }
    }
}

struct SingleMethodResult {
    vector<double> ecgTemplate;
    vector<double> ecg_template_iqr;   // empty for methods that don't compute std
    double ppg_alignment_point;
    int r_col = -1;   // true R column (alignment's r_aligned_col)

    // MEDIAN RR of the beats this template was built from, in samples.
    //
    // A DISPLAY WIDTH, NOT A STORAGE WIDTH. The beat matrix is framed on the
    // bin's LONGEST RR so that no beat ever loses a sample -- that part must
    // not change. But one 2.8 s pause then makes the array 5.6 s wide, and the
    // viewer sized its x-axis from the array length, so a normal 0.9 s beat was
    // drawn into a sixth of the panel. This is the number the axis should use
    // instead. It is carried in memory only and is NOT serialized: the viewer
    // receives TemplateBin directly, so no file format changes.
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



    // THE JOIN KEY: kept_idx[k] is the R-PAIR SLICE that produced the beat at
    // slot k of out_kept_beats.
    //
    // IT USED TO BE THE ALIGNED ROW (exactly usableIdx), and the difference is
    // silent. Each channel prunes independently, so slot k of CH1 and slot k of
    // PPG are different heartbeats and joining them needs an index they share;
    // the aligned row is not one, because the slicer SKIPS R-pairs (rr <= 3
    // samples, and rr > 4 s, which is a dropout gap rather than a beat) before
    // anything is pushed. On a bin where nothing was skipped the row and the
    // ordinal coincide; on a bin with one gap every later beat's ordinal is
    // short by one, and the error grows with each skip -- so a partition keyed
    // on it pairs each QRS with a later heartbeat's pulse, further off the
    // deeper into the bin you go.
    //
    // The aligned row has no consumer left: it existed for the morphology
    // writers, whose columns are now slices, and they resolve a waveform
    // through jbank's slice -> row map instead. One map, one meaning.
    std::vector<size_t> kept_idx;
};

static inline SingleMethodResult build_ecg_template_for_method(const vector<double>& ecgSignal, const vector<size_t>& rpeaks, const vector<vector<double>>& pairs,
    double ecgRate, vector<vector<double>>* out_kept_beats = nullptr, bool compute_iqr = false,
    // build_bank / bin_index / channel went with the per-channel bank. Left
    // unnamed rather than deleted so the existing call sites keep compiling,
    // and so passing `true` cannot quietly rebuild it.
    bool = false, uint64_t = 0, int = 0) {
    SingleMethodResult res;
    res.ecgTemplate = {};
    res.ecg_template_iqr = {};
    res.ppg_alignment_point = NaN;
    res.r_col = -1;

    if (rpeaks.size() < 2 || ecgSignal.empty() || ecgRate <= 0.0) return res;

    const alignment::ecg_beat_set aligned = alignment::extract_beats_and_align(ecgSignal, rpeaks, ecgRate);
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

    // ---- NO PER-CHANNEL BANK IS BUILT HERE ----------------------------
    //
    // build_ecg_template_for_method used to run bin_pipeline::runChannel over
    // aligned.beats and hand the result out as SingleMethodResult::bank_out,
    // which became EcgChannelResult::bank_out_raw. That was a partition of ONE
    // channel's beats, and Section 4.6 has exactly one partition, shared by the
    // three ECG leads and the pulse: jbank::buildBinBank, driven per bin from
    // make_averaged_templates.hpp.
    //
    // Its last consumer was the morphology archive, and that archive now reads
    // the joint projection. Keeping it would leave a second grouping of the same
    // beats in memory with nothing marking it as the stale one -- which is how
    // the screen and the files came to describe different partitions.
    //
    // IT WAS ALSO THE LAST PLACE THE ORDER RAN BACKWARDS. Its slot 0 came from
    // seed_pool::selectSeedPool filtered by the Tukey verdict, so prematurity
    // and Tukey steered the partition. The required order is partition first,
    // then remove premature, then Tukey on what is left, and with this gone
    // there is no code left that does it the other way.

    // ---- capture the beats handed downstream, with their rhythm verdicts --
    if (out_kept_beats) {
        out_kept_beats->clear();
        out_kept_beats->reserve(usable.size());
        for (const auto* sl : usable) out_kept_beats->push_back(*sl);
    }
    // Composed HERE, the one place aligned.slice_index and usableIdx are both
    // in scope: slot -> aligned row -> R-pair slice. usableIdx stays local,
    // because the aligned row is only needed to index aligned.* inside this
    // function (kept_rhythm below does exactly that).
    //
    // Falls back to the aligned row when alignment supplied no slice map at all
    // -- an input predating the map -- which keeps the old behaviour rather than
    // emitting zeros that would read as "every beat is slice 0".
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

    // Per-sample robust spread over the same aligned-beat matrix: per-sample
    // STD (ddof=1) -- changed from IQR per spec step 7. Used to draw the
    // gray band under the raw template.
    // NOTE: field/param names say "iqr" but this now computes the
    // per-sample STD (ddof=1), not the interquartile range -- changed
    // per spec step 7. Renaming ecg_template_iqr/compute_iqr throughout
    // the codebase (TemplateTypes.hpp, BinPlotWidget, TemplateBinIO,
    // template_io, ...) is a separate, larger follow-up; left as-is here
    // to keep this change to the computation only.
    if (compute_iqr) {
        // OVER THE MASKED REFERENCE POOL, the same set the median above used.
        // Computed over `usable` it described the spread of sinus AND ectopy
        // while the median described sinus alone, so the band drawn under the
        // template was the wrong width for the line it was drawn under -- and in
        // bigeminy it was roughly the distance between the two morphologies.
        const std::vector<const std::vector<double>*>& spreadSet =
            reference.empty() ? usable : reference;
        res.ecg_template_iqr.assign(maxLen, 0.0);
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
            res.ecg_template_iqr[c] = std::sqrt(sumsq / static_cast<double>(nc - 1));   // ddof = 1
        }
    }

    // PPG transit delay (median foot - R across paired beats), unchanged.
    if (!pairs.empty()) {
        std::vector<double> diffs;
        for (const auto& p : pairs) {
            if (p.size() >= 2 && p[0] >= 0 && p[1] >= 0 && p[0] != p[1])
                diffs.push_back(p[0] - p[1]);
        }
        if (!diffs.empty()) res.ppg_alignment_point = median(diffs);
    }
    return res;
}

static inline void init_channel_result(EcgChannelResult& cr, size_t n) {
    cr.ecgTemplates_raw.resize(n);
    cr.ecgTemplates_raw_iqr.resize(n);
    cr.ecgTemplates_squared.resize(n);
    cr.ecgTemplates_absval.resize(n);
    cr.ecgTemplates_unfiltered.resize(n);
    cr.ref_index_raw.resize(n, -1);

    cr.ppg_alignment_point_raw.resize(n, NaN);
    cr.ppg_alignment_point_squared.resize(n, NaN);
    cr.ppg_alignment_point_absval.resize(n, NaN);
    cr.ppg_alignment_point_unfiltered.resize(n, NaN);

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
    // Stage marker inside a single bin, so a stall that never reaches the
    // per-bin line below is still localized to a bin, a channel and a method.
    // Reports the R-peak count too, since that is the input size driving cost.
    std::fprintf(stderr, "  [ecgfast]   bin %llu ch%d: align+template, %zu rpeaks\n",
        (unsigned long long)i, channel_index, masterPeaks.size());
    std::fflush(stderr);

    auto raw_res = build_ecg_template_for_method(
        ecgSignal, masterPeaks, bin.pairs, ecgRate,
        capture, /*compute_iqr=*/true,
        /*unused*/false, static_cast<uint64_t>(i), channel_index);
    if (out_kept_idx) *out_kept_idx = std::move(raw_res.kept_idx);
    cr.ecgTemplates_raw[i] = raw_res.ecgTemplate;
    cr.ecgTemplates_raw_iqr[i] = raw_res.ecg_template_iqr;
    cr.ppg_alignment_point_raw[i] = raw_res.ppg_alignment_point;
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
    }

    // Method 4: unfiltered (original ECG signal + master R-peaks). No std.
    auto unfilt_res = build_ecg_template_for_method(
        origSignal, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_unfiltered[i] = unfilt_res.ecgTemplate;
    cr.ppg_alignment_point_unfiltered[i] = unfilt_res.ppg_alignment_point;
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
        sq_sig, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_squared[i] = sq_res.ecgTemplate;
    cr.ppg_alignment_point_squared[i] = sq_res.ppg_alignment_point;
    cr.r_col_squared[i] = sq_res.r_col;

    // Method 3: absval (abs-value signal + master R-peaks). No std.
    const auto& abs_sig = ch.absval_signal.empty() ? ecgSignal : ch.absval_signal;
    auto abs_res = build_ecg_template_for_method(
        abs_sig, masterPeaks, bin.pairs, ecgRate,
        nullptr, /*compute_iqr=*/false);
    cr.ecgTemplates_absval[i] = abs_res.ecgTemplate;
    cr.ppg_alignment_point_absval[i] = abs_res.ppg_alignment_point;
    cr.r_col_absval[i] = abs_res.r_col;
}

// Original all-4-methods entry point, preserved by composition.
static inline void process_channel(
    EcgChannelResult& cr,
    const vector<output_binfile_data>& bins,
    size_t i,
    const vector<double>& ecgSignal,
    const vector<double>& origSignal,
    const ChannelRPeaks& ch,
    const vector<size_t>& masterPeaks,
    double ecgRate,
    bool capture_raw_beats = false,
    int channel_index = 0)
{
    process_channel_fast(cr, bins, i, ecgSignal, origSignal, masterPeaks,
        ecgRate, capture_raw_beats, channel_index);
    process_channel_slow(cr, bins, i, ecgSignal, ch, masterPeaks,
        ecgRate);
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

        const int _d = ++_done;
        std::fprintf(stderr, "  [ecgfast] bin %d/%zu done (rpeaks=%zu)\n",
            _d, n, master.size());
        std::fflush(stderr);
    }

    // Single-threaded, post-loop: write the per-beat vertical move log
    // (CH1 truncates+headers, CH2/CH3 append).
    // ---- EVERY POST-LOOP STEP TIMED AND NAMED --------------------------
    // The parallel loop's last output is "[ecgfast] bin N/N done", and
    // everything after it is single-threaded and silent: seven separate writes
    // with no output between them. A stall anywhere in here is indistinguishable
    // from a stall in the loop's final bin, which is where two diagnoses in a
    // row went wrong. Each step announces itself before starting and reports
    // elapsed time after, flushed.
    auto _t = std::chrono::steady_clock::now();
    auto _begin = [&](const char* what) {
        std::fprintf(stderr, "  [postloop] starting %s ...\n", what);
        std::fflush(stderr);
        };
    auto _step = [&](const char* what) {
        const auto now = std::chrono::steady_clock::now();
        std::fprintf(stderr, "  [postloop] %-22s %9.1f ms\n", what,
            std::chrono::duration<double, std::milli>(now - _t).count());
        std::fflush(stderr);
        _t = now;
        };

    _begin("move log CH1-3");
    ecg_move_log::write_channel("CH1", res.ch1.tp_shift_raw, res.ch1.pq_shift_raw, /*first=*/true);
    ecg_move_log::write_channel("CH2", res.ch2.tp_shift_raw, res.ch2.pq_shift_raw, /*first=*/false);
    ecg_move_log::write_channel("CH3", res.ch3.tp_shift_raw, res.ch3.pq_shift_raw, /*first=*/false);
    _step("move log CH1-3");

    // The join key, surfaced so the partition and the archive read the SAME map
    // rather than two copies that can drift. Moved, not copied: keptIdx dies
    // with this function otherwise.
    for (int c = 0; c < 3; ++c) res.kept_index[c] = keptIdx[c];

    // ---- THE MORPHOLOGY WRITERS ARE NOT CALLED HERE ANY MORE -------------
    //
    // writeTemplates / writeBeatsBin / writeTemplatesBin ran in this post-loop
    // slot, from blocks built on EcgChannelResult::bank_out_raw. They cannot
    // stay: the archive has to describe the JOINT partition, and the joint bank
    // does not exist yet at this point in the run -- it is built per bin in
    // make_averaged_templates.hpp, after this function returns. Writing here
    // meant the archive described a per-channel partition while the screen and
    // the serialized templates described the joint one, with nothing in either
    // file saying they were different partitions of the same beats.
    //
    // They now run at the end of GenerateTemplatesFast, from the projection,
    // and they include a PPG block -- which this call site could not produce at
    // all, because the pulse channel is not an EcgChannelResult.

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

inline EcgTemplateResult CreateEcgTemplates(
    const vector<output_binfile_data>& bins,
    double ecgRate)
{
    EcgTemplateResult res = CreateEcgTemplatesFast(bins, ecgRate);
    CreateEcgTemplatesSlow(bins, ecgRate, res);
    return res;
}