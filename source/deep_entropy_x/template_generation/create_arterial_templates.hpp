/**
 * @file   create_arterial_templates.hpp
 * @brief  ONE pipeline: R-anchored pulse templates for every pulse channel.
 *
 *         PPG, ABP, ART and ART_PULM are all sliced on the same real-time
 *         windows the ECG uses -- [t_R_i - pad, t_R_{i+1} + pad], driven by
 *         ch1.raw R-peaks (ECG-frame samples), converted to this channel's own
 *         samples via the rate ratio (channelRate / ecgRate). See
 *         CreatePulseTemplates / build_pulse_template_pair_windowed.
 *
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>
#include "fiducial_marker_finding/alignment.hpp"
#include "template_structs.hpp"
#include "template_generation/normalize_template_amplitude.hpp"

struct PPGTemplatesResult {
    vector<vector<double>> templates;   // [bin][sample]
    vector<vector<double>> iqrs;        // [bin][sample], same shape as templates
    vector<vector<vector<double>>> kept; // [bin][beat][sample] retained snips
    vector<int> peakCol;                // [bin] systolic peak column (R1..R2)
    vector<int> footCol;                // [bin] foot column (R1..peak)

    // R-PAIR ORDINAL of each retained snip: keptSlices[bin][beat] is the index
    // of the R-pair that snip was sliced from, parallel to kept[bin].
    //
    // This is the join key between a pulse and a QRS. `kept` is pruned twice in
    // the aligner and once more by the fit-error filter here, so kept[bin][k] is
    // NOT R-pair k, and without the ordinal there is no way to say a PPG beat
    // and an ECG beat are the same heartbeat. Any consumer treating the two
    // channels as views of one beat needs it.
    vector<vector<uint32_t>> keptSlices;
};

/**
 * @brief  Shared "slice a channel by ch1.raw R-pairs" pulse-averager.
 *
 *         For every consecutive R-pair (R_i, R_{i+1}) in `masterPeaksEcg`
 *         (ECG-frame sample indices from bin.ch1.raw), compute the real-time
 *         window [t_R_i - pad, t_R_{i+1} + pad] and pull it out of `signal`
 *         (which lives in this channel's own sample space at `channelRate`).
 *         The R sample index maps into channel space as
 *         `round(r_ecg * channelRate / ecgRate)`.
 *
 *         Slices are aligned so that column 0 is `pad` seconds before the
 *         first R -- so R_first sits at column `padSamplesCh` in every
 *         slice. Variable-length RR => variable slice length; short slices
 *         contribute NaN past their real end. Column-wise NaN-skipping
 *         median => the template.
 *
 *         Shared by CreatePulseTemplates for every channel (PPG, ABP, ART,
 *         ART_PULM) via a member-pointer for the signal.
 */
 // ==========================================================================
 // PULSE QC: ONE ERROR THRESHOLD, FROM config.csv
 // ==========================================================================
 //
 // A candidate pulse is kept when its normalized foot-to-foot fit error against
 // the bin's median reference is below this fraction:
 //
 //     err = || beat - reference || / || reference ||   over the f2f window
 //
 // So 0.10 means "within 10% of the reference by RMS". Set from config.csv as a
 // PERCENT (ppg_fit_error_pct), converted once here.
 //
 // WHY IT IS A RUNTIME VALUE. It is the single number that decides how much of
 // the pulse channel survives, and it needs to differ by dataset: an arterial
 // line is far more repeatable than a sleep-study pulse-ox, and a threshold
 // tuned on one throws away most of the other. On a MESA record 10% retained
 // 7.5% of the channel.
 //
 // WORTH KNOWING WHAT THE METRIC IS BLIND TO. This error is SCALE-SENSITIVE: a
 // pulse of identical shape with 15% more amplitude scores 0.15 and is rejected
 // at a 10% threshold. Pulse amplitude modulates with respiration and vasomotion
 // as a matter of course, so part of what this threshold controls is tolerance
 // to normal amplitude variation rather than to shape. The ECG side judges shape
 // by correlation, which is immune to exactly that. Raising the percentage is a
 // workaround for the metric, not a fix to it.
namespace pulse_qc {

    inline constexpr double kDefaultFitErrorFraction = 0.10;   // 10%

    namespace detail { inline double g_fit_error = kDefaultFitErrorFraction; }

    inline double fitErrorFraction() { return detail::g_fit_error; }

    // `pct` is a PERCENT: 10 means 10%. Returns false and changes nothing when
    // it is outside (0, 100].
    //
    // A blank config cell parses to 0.0 through the loader, and a threshold of
    // 0 admits no pulse at all -- every bin degenerate, no pulse template
    // anywhere, and the only symptom a channel that quietly vanished. So an
    // unusable value leaves the default in place.
    inline bool setFitErrorPct(double pct) {
        if (!(pct > 0.0 && pct <= 100.0)) return false;
        detail::g_fit_error = pct / 100.0;
        return true;
    }

}  // namespace pulse_qc

// One bin's worth of output. WAS SEVEN OUT-PARAMETERS -- three vector refs,
// two int refs and an optional pointer -- written into six parallel arrays of
// PPGTemplatesResult by the single caller, which then had to clear all six by
// hand in its catch block. Add a seventh output to that arrangement and the
// catch has to learn about it too, or a bin survives half-written: a template
// with stale keptSlices, or a peakCol from an attempt that threw.
//
// Returned by value instead, so the empty state IS the default state and the
// failure path is one assignment.
struct PulseTemplateBin {
    std::vector<double> tmpl;                 // column-wise NaN-skipping median
    std::vector<double> iqr;                  // local-ratio IQR about footCol
    std::vector<std::vector<double>> kept;    // [beat][sample] retained snips
    std::vector<uint32_t> keptSlices;         // R-pair ordinal per retained snip
    int peakCol = -1;                         // systolic peak column
    int footCol = -1;                         // foot column
};

static inline PulseTemplateBin build_pulse_template_pair_windowed(
    const std::vector<double>& signal,
    double channelRate,
    const std::vector<size_t>& masterPeaksEcg,
    double ecgRate,
    double padSeconds,
    // For the [pulseqc] line only. Passed rather than inferred because this
    // function has no other way to name the bin it is working on, and a
    // retention report that cannot say WHICH bin is nearly useless.
    size_t bin_index = 0)
{
    PulseTemplateBin out;

    if (signal.empty() || masterPeaksEcg.size() < 2 ||
        channelRate <= 0.0 || ecgRate <= 0.0) return out;

    const double scale = channelRate / ecgRate;

    // Convert ch1.raw R-peaks from ECG samples to this channel's samples.
    std::vector<size_t> peaksCh;
    peaksCh.reserve(masterPeaksEcg.size());
    for (size_t r : masterPeaksEcg)
        peaksCh.push_back(static_cast<size_t>(std::llround(
            static_cast<double>(r) * scale)));

    // Per-bin peak-aligned + foot-vertical-aligned beat matrix.
    const auto aligned = alignment::extract_ppg_beats_and_align(signal, peaksCh, channelRate);
    if (aligned.beats.empty()) return out;

    // ---- Matched-filter QC, two-pass, per spec:
    //   (a) build a REFERENCE template as the column-wise NaN-skipping
    //       median across ALL candidate beats. The median is robust to
    //       outliers without needing to pick a fixed "seed" count.
    //   (b) score every candidate against the reference by normalized
    //       error ||beat - ref|| / ||ref||; keep beats whose error is
    //       below 5%.
    // The final template below is then rebuilt from the survivors,
    // giving the two-pass: median-of-all -> reject high-error ->
    // re-median. Same wave-score pruning logic as ECG, adapted to PPG's
    // normalized-error metric. Falls back to keeping everything if the
    // filter would otherwise reject the whole set (degenerate reference).
    std::vector<std::vector<double>> filteredBeats;
    // The survivor ROW INDICES, kept alongside the waveforms. filteredBeats
    // loses them -- it is a copy of the rows, not a view of them -- and the
    // fiducial block below needs them to read aligned.peak_cols/foot_cols for
    // the same beats the template was medianed from. Same loop fills both, so
    // they cannot fall out of step.
    std::vector<size_t> survivorsForMarks;
    // Diagnostic accumulators, populated inside the QC block.
    int diag_ref_col_early = 0, diag_ref_col_mid = 0, diag_ref_col_late = 0;
    int diag_ref_defined_cols = 0;
    int diag_input_beats = 0, diag_survivors = 0;
    double diag_err_min = std::numeric_limits<double>::infinity();
    double diag_err_max = -std::numeric_limits<double>::infinity();
    std::vector<double> diag_all_errs;
    {
        const int w = static_cast<int>(aligned.beats.front().size());
        diag_input_beats = static_cast<int>(aligned.beats.size());

        // (a) reference = column-wise median across ALL candidates.
        std::vector<double> reference(w, NaN);
        std::vector<int> col_counts(w, 0);   // for diagnostic
        for (int c = 0; c < w; ++c) {
            std::vector<double> col;
            col.reserve(aligned.beats.size());
            for (const auto& sl : aligned.beats)
                if (c < (int)sl.size() && !std::isnan(sl[c])) col.push_back(sl[c]);
            col_counts[c] = static_cast<int>(col.size());
            if (col.empty()) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            reference[c] = (nc % 2)
                ? col[mid]
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + col[mid]);
            ++diag_ref_defined_cols;
        }
        // Sample the column-count profile at three positions.
        diag_ref_col_early = col_counts[w / 8];
        diag_ref_col_mid = col_counts[w / 2];
        diag_ref_col_late = col_counts[(7 * w) / 8];

        // ---- Foot-to-foot region of the reference pulse (spec steps 2-3:
        // the rejection math runs on the single pulse, foot to foot -- NOT
        // on the full ECG-length window). All beats are up50-aligned to
        // shared columns, so the reference's landmarks are every beat's
        // landmarks: systolic peak = argmax; first foot = reference min on
        // [0, peak]; second foot = reference min on [peak, w). The output
        // window/template stays the full ECG length -- only the accept/
        // reject error below is restricted to [firstFoot, secondFoot]. The
        // ECG-length window carries padding (previous pulse's tail before
        // foot 1, next pulse's onset past foot 2) that varies beat-to-beat;
        // scoring the whole window would let that padding dominate the error
        // and reject clean pulses.
        // Systolic peak of the reference pulse from the upstroke, not argmax:
        // the two feet below are found relative to it, so a peak on the
        // reflected wave puts foot 1 at the dicrotic notch and narrows the
        // foot-to-foot error window to the back half of the pulse.
        int refPeak = FeatureMarks::detect_ppg_upstroke_peak(reference, 0, w);

        int f2fLo = 0, f2fHi = w;   // safe default: whole window
        if (refPeak > 0 && refPeak < w - 1) {
            // Both feet through the shared trough primitive.
            const int fl = FeatureMarks::trough_in(reference, 0, refPeak);
            const int fr = FeatureMarks::trough_in(reference, refPeak, w - 1);
            if (fl >= 0 && fr > fl) { f2fLo = fl; f2fHi = fr + 1; }   // incl. 2nd foot
        }

        // Normalized error restricted to [f2fLo, f2fHi): ||beat - ref|| /
        // ||ref|| over that column band, non-NaN overlap only. Same formula
        // as ppg_deriv's deleted normalizedError was, just windowed to the
        // foot-to-foot span.
        auto footToFootError = [&](const std::vector<double>& bt) -> double {
            double num = 0.0, den = 0.0; int overlap = 0;
            const int hi = std::min<int>(f2fHi,
                std::min<int>(static_cast<int>(bt.size()),
                    static_cast<int>(reference.size())));
            for (int c = std::max(0, f2fLo); c < hi; ++c) {
                if (std::isnan(bt[c]) || std::isnan(reference[c])) continue;
                const double e = bt[c] - reference[c];
                num += e * e; den += reference[c] * reference[c]; ++overlap;
            }
            if (overlap == 0 || den <= 0.0) return std::numeric_limits<double>::infinity();
            return std::sqrt(num / den);
            };

        // (b) per-pulse accept/reject on the foot-to-foot error.
        // SURVIVORS ARE TRACKED BY ROW INDEX, not by copying waveforms. The
        // row index is what carries the R-pair ordinal
        // (aligned.original_index), so a filter that only accumulates
        // waveforms discards the join key -- which is what this loop used to
        // do, and the reason a partition shared with the ECG channels could
        // not be built at all.
        std::vector<size_t> survivorRows;
        survivorRows.reserve(aligned.beats.size());
        diag_all_errs.reserve(aligned.beats.size());
        for (size_t k = 0; k < aligned.beats.size(); ++k) {
            const double err = footToFootError(aligned.beats[k]);
            diag_all_errs.push_back(err);
            if (std::isfinite(err)) {
                if (err < diag_err_min) diag_err_min = err;
                if (err > diag_err_max) diag_err_max = err;
            }
            if (err < pulse_qc::fitErrorFraction()) survivorRows.push_back(k);
        }
        // THE ESCALATION LADDER IS GONE. It ran 10% -> 20% -> 50%, taking the
        // first tier that reached a survivor floor, and it was the wrong shape
        // of fix twice over.
        //
        // It hid the problem. The tiers only fired when survivors fell below a
        // COUNT, so a bin keeping 64 pulses of 820 cleared the floor and
        // reported success -- and on a real MESA record this filter was
        // retaining 7.5% of the channel overall, 2030 pulses of 26970, with the
        // ladder never firing anywhere.
        //
        // And when it did fire it made the result worse quietly. A bin relaxed
        // to 50% builds its reference from pulses that failed the 10% test, so
        // two bins in the same record could have templates fitted to
        // populations selected by different standards, with nothing written down
        // to say which.
        //
        // One threshold now, from config.csv (ppg_fit_error_pct), applied
        // uniformly to every bin. A bin that keeps very few pulses keeps very
        // few, and says so on the line below rather than being rescued into
        // looking fine.

        // DEGENERATE REFERENCE -> KEEP ALL, the same fallback the arterial
        // twin has at the bottom of this file. All three tiers can come back
        // empty when footToFootError() is non-finite for every candidate, which
        // happens when the reference template is all-NaN or the pulse signal is
        // flat or dead. The original accumulation loop happened never to leave
        // filteredBeats empty, so the unguarded beatsForTemplate.front() below
        // was latently wrong; tracking survivor ROWS made it reachable and it
        // faulted. Keeping all candidates preserves the old outcome -- a
        // template built from everything, which the QC diagnostics then report
        // as low quality -- rather than dropping the bin entirely.
        // ---- ZERO SURVIVORS MEANS ZERO, NOT EVERYTHING -------------------
        //
        // This used to refill survivorRows with every candidate when the
        // threshold rejected them all. It INVERTED THE FILTER: the worse the
        // bin, the more pulses it kept. On a real record the bins reporting
        // "100.0% kept" had error medians of 0.27, 0.36, 0.52, even 1.12 --
        // every one of them a bin where nothing passed and the guard then
        // admitted the lot, garbage included, while a healthy bin next to it
        // kept 12% and reported an error median of 0.105.
        //
        // A bin that fails the filter now produces no pulse template. The
        // out-params were cleared at entry and the columns are -1, which is
        // exactly the state the callers already read as "no pulse for this
        // bin" -- the same state the catch(...) in CreatePulseTemplates
        // produces. The [pulseqc] line above says how many were offered and
        // how many passed, so the bin is accounted for rather than silently
        // absent.
        if (survivorRows.empty()) {
            // Median computed here rather than reused: the shared
            // diag_err_median is declared after this block, and returning
            // early means it is never reached.
            double medHere = std::numeric_limits<double>::quiet_NaN();
            {
                std::vector<double> fin;
                fin.reserve(diag_all_errs.size());
                for (const double e : diag_all_errs)
                    if (std::isfinite(e)) fin.push_back(e);
                if (!fin.empty()) {
                    std::sort(fin.begin(), fin.end());
                    medHere = fin[fin.size() / 2];
                }
            }
            return out;
        }

        // Waveform and ordinal appended in the SAME loop, so they cannot fall
        // out of step.
        filteredBeats.reserve(survivorRows.size());
        out.keptSlices.reserve(survivorRows.size());
        const bool ordinalsUsable =
            aligned.original_index.size() == aligned.beats.size();
        survivorsForMarks.reserve(survivorRows.size());
        for (const size_t k : survivorRows) {
            filteredBeats.push_back(aligned.beats[k]);
            survivorsForMarks.push_back(k);
            out.keptSlices.push_back(ordinalsUsable
                ? aligned.original_index[k] : static_cast<uint32_t>(k));
        }
        // Said out loud rather than papered over. Falling back to the row index
        // produces a mapping of the right SHAPE and the wrong CONTENT, and every
        // consumer downstream would treat it as a valid join key.
        if (!ordinalsUsable)
        {
            diag_survivors = static_cast<int>(filteredBeats.size());
        }
    }
    const auto& beatsForTemplate = filteredBeats;

    // Peak-column distribution of surviving beats, argmax of each (skip NaN).
    // If they cluster tightly, the template's median peak has a lot of support;
    // if the count of contributing beats drops fast past that column, that's
    // the cutoff signature.
    int diag_peak_min = std::numeric_limits<int>::max();
    int diag_peak_max = std::numeric_limits<int>::min();
    std::vector<int> diag_peak_cols;
    diag_peak_cols.reserve(beatsForTemplate.size());
    for (const auto& b : beatsForTemplate) {
        int pk = -1; double pv = -std::numeric_limits<double>::infinity();
        for (int c = 0; c < (int)b.size(); ++c)
            if (!std::isnan(b[c]) && b[c] > pv) { pv = b[c]; pk = c; }   // see note below
        if (pk >= 0) {
            diag_peak_cols.push_back(pk);
            if (pk < diag_peak_min) diag_peak_min = pk;
            if (pk > diag_peak_max) diag_peak_max = pk;
        }
    }
    int diag_peak_median = -1;
    if (!diag_peak_cols.empty()) {
        auto tmp = diag_peak_cols;
        std::sort(tmp.begin(), tmp.end());
        diag_peak_median = tmp[tmp.size() / 2];
    }
    // Median normalized error (rank quality metric).
    double diag_err_median = std::nan("");
    if (!diag_all_errs.empty()) {
        std::vector<double> tmp;
        tmp.reserve(diag_all_errs.size());
        for (double e : diag_all_errs) if (std::isfinite(e)) tmp.push_back(e);
        if (!tmp.empty()) {
            std::sort(tmp.begin(), tmp.end());
            diag_err_median = tmp[tmp.size() / 2];
        }
    }

    const size_t maxLen = beatsForTemplate.front().size();

    // Column-wise NaN-skipping median => template.
    out.tmpl.assign(maxLen, NaN);
    for (size_t c = 0; c < maxLen; ++c) {
        std::vector<double> col;
        col.reserve(beatsForTemplate.size());
        for (const auto& sl : beatsForTemplate) {
            const double v = sl[c];
            if (!std::isnan(v)) col.push_back(v);
        }
        if (col.empty()) continue;
        // Median via nth_element (O(k)) rather than a full sort (O(k log k));
        // this column loop runs maxLen times over up to beats.size() values,
        // and in messy bins the escalated QC can keep the whole bin, so the
        // full sort here was a hot path.
        const size_t nc = col.size();
        const size_t mid = nc / 2;
        std::nth_element(col.begin(), col.begin() + mid, col.end());
        const double hi = col[mid];
        out.tmpl[c] = (nc % 2)
            ? hi
            : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi);
    }

    // ---- PPG fiducials FROM THE ALIGNER'S OWN PER-BEAT MARKS -------------
    //
    // TWO WRONG VERSIONS PRECEDED THIS ONE, and both failed the same way: they
    // decided WHERE to look without anything entitling them to.
    //
    //   (1) "R1 lands at column padSeconds*channelRate by construction", then
    //       peak = argmax over [R1, R1+RR], foot = argmin over [R1, peak].
    //       True of the SLICER's frame, false of this one:
    //       extract_ppg_beats_and_align shifts every beat so its 50%-upslope
    //       crossing lands on the median up50 column, so the shared axis has
    //       its origin at up50_aligned_col and R1 is at no fixed column at
    //       all. When padSeconds*channelRate landed past the true foot the
    //       window [r1, peak] collapsed onto the apex and argmin returned a
    //       sample beside the peak -- the peak was marked as the foot.
    //
    //   (2) detect_ppg_fiducials over the WHOLE template. Right detector,
    //       no bracket: the window spans [R_i - pad, R_i+1 + pad], so it
    //       opens inside the PREVIOUS pulse's diastolic tail. The first
    //       upstroke peak in it can be that tail's rebound, and the walk-back
    //       from there lands the foot near the left edge of the window.
    //
    // The aligner already measured both marks on every beat, in THIS frame,
    // with the same primitives (detect_ppg_upstroke_peak bounded to the beat's
    // own R-R window, then trough_in back from it), and shifted them onto the
    // shared axis alongside the samples: peak_cols[k] = prepend + peak,
    // foot_cols[k] = prepend + foot. Those were computed and then read by
    // nobody. The template is the column-wise median of the same beats, so the
    // median of their marks is the mark of the median -- no search, no window,
    // no constant, and nothing that can drift from the alignment axis because
    // it IS the alignment axis.
    //
    // NEGATIVES EXCLUDED, not clamped. A beat whose up50 sat right of the
    // anchor had its leading samples clipped, so its foot column can be < 0
    // (see the `prepend` note in extract_ppg_beats_and_align). Clamping those
    // to 0 would drag the median toward the left edge; they are simply not
    // evidence about where the foot is.
    //
    // ONLY THE SURVIVORS vote. survivorRows indexes aligned.beats, and the
    // template was built from exactly those rows, so the marks and the
    // waveform describe one population.
    // Computed BEFORE the spread below, since the spread needs out.footCol.
    {
        const int N = static_cast<int>(out.tmpl.size());
        auto medianOf = [](std::vector<int> v) -> int {
            if (v.empty()) return -1;
            std::sort(v.begin(), v.end());
            return v[v.size() / 2];
            };
        std::vector<int> pks, fts;
        pks.reserve(survivorsForMarks.size());
        fts.reserve(survivorsForMarks.size());
        for (const size_t k : survivorsForMarks) {
            if (k < aligned.peak_cols.size() && aligned.peak_cols[k] >= 0
                && aligned.peak_cols[k] < N) pks.push_back(aligned.peak_cols[k]);
            if (k < aligned.foot_cols.size() && aligned.foot_cols[k] >= 0
                && aligned.foot_cols[k] < N) fts.push_back(aligned.foot_cols[k]);
        }
        const int pk = medianOf(std::move(pks));
        const int ft = medianOf(std::move(fts));
        // FOOT BEFORE PEAK or neither is trusted. The two medians are taken
        // independently, so a bin whose beats disagree badly enough to invert
        // them has no usable foot -- and -1 is the state every caller already
        // reads as "no pulse fiducial for this bin", rather than a pair that
        // would make local_ratio_iqr integrate about a maximum.
        if (pk >= 0 && ft >= 0 && ft < pk) {
            out.peakCol = pk;
            out.footCol = ft;
        }
    }

    // Per-sample robust spread, in the SAME units the mean trace will
    // eventually be displayed/exported in (minus only the final /ref
    // division, which normalize_features::scale_array_by_ref applies at
    // display time -- never here). Each beat is first converted to its own
    // local perfusion-index ratio using ITS OWN foot (out.footCol), per the
    // documented algorithm -- never a median/global foot -- then the
    // cross-beat IQR is taken of those local-ratio values. This is NOT the
    // same as taking the IQR of raw amplitudes, because the per-sample
    // transform's slope differs beat-to-beat (each beat has its own foot).
    out.iqr = (out.footCol >= 0)
        ? normalize_features::local_ratio_iqr(beatsForTemplate, out.footCol)
        : std::vector<double>(maxLen, 0.0);

    // Retain aligned per-beat slices for downstream (snips CSV, etc).
    out.kept.reserve(beatsForTemplate.size());
    for (const auto& sl : beatsForTemplate) out.kept.push_back(sl);

    // padSeconds is genuinely unused HERE now. It sizes the slice upstream and
    // nothing in this function depends on where that puts R1 any more, which is
    // the point of the change above -- so the 0.4 in this function's signature
    // disagreeing with the 0.25 in CreatePulseTemplates' doc block can no
    // longer move a fiducial. Worth reconciling anyway, since the slice width
    // still comes from it.
    (void)padSeconds; (void)channelRate; (void)ecgRate;

    // THE SUCCESS PATH HAD NO RETURN. Every early exit above returns
    // `out` explicitly, so the one path that actually builds a template
    // fell off the end of a non-void function -- undefined behaviour, and
    // MSVC reports it as C4715 (a WARNING, which is why a build that looked
    // clean shipped it). What the caller received was a return object that
    // kept the two ints and lost all four vectors, so every bin came back
    // with a plausible peakCol/footCol and an EMPTY waveform: template_good
    // stayed false for the whole record, ppg_template_good with it, and the
    // pulse channel went missing from the joint bank, the panels, the
    // templates file and the viewer -- while this function did all the work
    // and charged the 500 ms for it.
    return out;
}

/**
 * @brief  PPG templates for every bin. Slicing is driven by bin.ch1.raw
 *         (ECG-frame R-peaks), consistent with the ECG templater.
 *
 * @param bins        Input bins.
 * @param ecgRate     ECG sample rate (for R-peak time base).
 * @param ppgRate     PPG sample rate.
 * @param padSeconds  0.25 (matches CreateEcgTemplates).
 */
inline PPGTemplatesResult CreatePulseTemplates(
    const vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double ecgRate,
    double channelRate,
    double padSeconds = 0.4)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.iqrs.assign(n, {});
    out.kept.assign(n, {});
    out.keptSlices.assign(n, {});
    out.peakCol.assign(n, -1);
    out.footCol.assign(n, -1);

    if (channelRate <= 0.0) return out;   // channel absent from this dataset

    int ppg_threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(ppg_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& b = bins[i];
        if (b.bad_segment || (b.*sigMember).empty() || b.ch1.raw.size() < 2)
            continue;
        try {
            PulseTemplateBin r = build_pulse_template_pair_windowed(
                b.*sigMember, channelRate, b.ch1.raw, ecgRate, padSeconds, i);
            out.templates[i] = std::move(r.tmpl);
            out.iqrs[i] = std::move(r.iqr);
            out.kept[i] = std::move(r.kept);
            out.keptSlices[i] = std::move(r.keptSlices);
            out.peakCol[i] = r.peakCol;
            out.footCol[i] = r.footCol;
        }
        catch (...) {
            // NOTHING TO CLEAR. The six assignments above are the last thing
            // the try does, so a throw leaves this bin at the empty state the
            // arrays were initialised to -- and that state cannot go stale as
            // fields are added, which is what the by-hand version could not
            // promise (a bin whose waveforms were discarded but whose
            // keptSlices survived would present a join key pointing at beats
            // that are no longer there).
        }
    }

    return out;
}