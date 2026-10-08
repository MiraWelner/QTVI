/**
 * @file   bin_pulse.hpp
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
#include "../stats_utils.hpp"   // pearson, the shared primitive
#include "fiducial_marker_finding/alignment.hpp"
#include "template_generation/normalize_template_amplitude.hpp"
#include "logging/ptt_logging.hpp"

struct PPGTemplatesResult {
    std::vector<std::vector<double>> templates;   // [bin][sample]
    std::vector<std::vector<double>> iqrs;        // [bin][sample], same shape as templates
    std::vector<std::vector<std::vector<double>>> kept; // [bin][beat][sample] retained snips
    std::vector<int> peakCol;                // [bin] systolic peak column (R1..R2)
    std::vector<int> footCol;                // [bin] foot column (R1..peak)
    std::vector<double> notch_from_peak;     // [bin] median (notch - peak), samples
    std::vector<int> notch_n;                // [bin] member beats that had one
    std::vector<double> notchCol;            // [bin] peakCol + notch_from_peak, -1 = none
    // [bin] the template's MEASURED R column, -1 where unmeasurable. The
    // pulse has no R-relative time axis without it; see
    // alignment::PpgBeatSet::r_cols for why it is not a constant.
    std::vector<int> rCol;

    // R-PAIR ORDINAL of each retained snip: keptSlices[bin][beat] is the index
    // of the R-pair that snip was sliced from, parallel to kept[bin].
    //
    // This is the join key between a pulse and a QRS. `kept` is pruned twice in
    // the aligner and once more by the fit-error filter here, so kept[bin][k] is
    // NOT R-pair k, and without the ordinal there is no way to say a PPG beat
    // and an ECG beat are the same heartbeat. Any consumer treating the two
    // channels as views of one beat needs it.
    std::vector<std::vector<uint32_t>> keptSlices;

    // [bin] every measured pulse beat's R / peak / 50% times, for
    // <stem>_ptt.csv. Includes beats the QC filter rejected (kept = false).
    std::vector<std::vector<ptt_log::Beat>> ptt;
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
 // PULSE QC: ONE CORRELATION FLOOR, FROM config.csv
 // ==========================================================================
 //
 // A candidate pulse is kept when its SHAPE correlates with the bin's median
 // reference over the foot-to-foot window:
 //
 //     r = pearson(beat, reference)   over [f2fLo, f2fHi)
 //
 // kept iff r >= corrFloor(). Set from config.csv as a CORRELATION
 // (pulse_qc_corr_floor), not a percent.
 //
 // WAS A NORMALIZED RMS ERROR, ||beat - reference|| / ||reference||, kept
 // below 10%. That metric was SCALE-SENSITIVE: a pulse of identical shape with
 // 15% more amplitude scored 0.15 and was rejected at a 10% threshold -- and
 // pulse amplitude modulates with respiration and vasomotion as a matter of
 // course, so most of what the threshold controlled was tolerance to normal
 // amplitude variation rather than to shape. Raising the percentage was a
 // workaround for the metric, not a fix to it. Correlation judges shape alone,
 // which is what the ECG side has always done -- and through the same
 // stats_utils::pearson the bank's correlate() uses, so the two sides now
 // measure similarity with one implementation.
 //
 // THE NEW BLIND SPOT IS THE MIRROR OF THE OLD ONE, and it is not small.
 // Pearson is invariant to BOTH scale and offset, so a pulse of identical
 // shape at twice the amplitude, or riding a large baseline offset, now
 // passes. This gate no longer constrains amplitude at all. If amplitude needs
 // constraining it needs its own bound; a lower r will not do it, because r is
 // not measuring amplitude to begin with.
 //
 // EXPECT RETENTION TO MOVE, NOT DRIFT. r >= 0.95 is not the complement of
 // err < 0.05 -- they are different criteria. At 10% error a MESA record
 // retained 7.5% of the channel (2030 pulses of 26970). Correlation admits
 // everything the error metric rejected on amplitude alone, so the figure
 // should rise substantially, and wants re-checking per dataset rather than
 // assuming.
 //
 // STILL ONE THRESHOLD FOR FOUR CHANNELS. CreatePulseTemplates is called per
 // channel through a member pointer (PPG, ABP, ART, ART_PULM), and an arterial
 // line is far more repeatable than a sleep-study pulse-ox. If the four
 // diverge under correlation this wants to become per-channel, the way the
 // morphology thresholds are.
namespace pulse_qc {

    inline constexpr double kDefaultCorrFloor = 0.95;

    namespace detail { inline double g_corr = kDefaultCorrFloor; }

    inline double corrFloor() { return detail::g_corr; }

    // `r` is a CORRELATION: 0.95 means 0.95. Returns false and changes nothing
    // when it is outside (0, 1] -- note the range, which is NOT the (0, 100]
    // that setFitErrorPct took.
    //
    // A blank config cell parses to 0.0 through the loader, and a floor of 0
    // admits every pulse whatever its shape -- the mirror of the old hazard,
    // where a threshold of 0 admitted none. So an unusable value leaves the
    // default in place.
    inline bool setCorrFloor(double r) {
        if (!(r > 0.0 && r <= 1.0)) return false;
        detail::g_corr = r;
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
    std::vector<double> std;                  // local-ratio IQR about footCol
    std::vector<std::vector<double>> kept;    // [beat][sample] retained snips
    std::vector<uint32_t> keptSlices;         // R-pair ordinal per retained snip
    std::vector<ptt_log::Beat> ptt;           // every measured beat; see ptt_log
    int peakCol = -1;                         // systolic peak column
    int footCol = -1;                         // foot column
    // THE DICROTIC NOTCH, in this template's own columns. NOT detected on the
    // template: E-5 has no per-template procedure, and the one it does have
    // runs on 4-second windows of the recording. This is the MEDIAN of the
    // member beats' own notch columns (PpgBeatSet::notch_cols), which is the
    // same reduction the waveform itself gets -- the template is a column
    // median over these beats, so its notch is the median of theirs.
    //
    // -1 when no member beat had a notch. notch_n is how many did, i.e. this
    // template's detectability, so a notch resting on few beats is visible as
    // such rather than looking as firm as one resting on all of them.
    // Median (notch - systolic peak) over the member beats, in samples.
    // NaN when no member beat had a notch. An OFFSET, not a column: see the
    // note where it is computed.
    double notch_from_peak = std::numeric_limits<double>::quiet_NaN();
    int    notch_n = 0;
    // The same notch as a COLUMN of this template, = peakCol + notch_from_peak.
    // -1 when there is no notch or it falls outside the template.
    double notchCol = -1.0;
    // THE TEMPLATE'S R COLUMN, median of the surviving beats' own R columns
    // in the shared frame (alignment::PpgBeatSet::r_cols). -1 = not
    // measurable, which every consumer must read as "this pulse has no
    // R-relative time axis" rather than substituting a constant. See the note
    // on r_cols for what the constant used to cost.
    int rCol = -1;
};

static inline PulseTemplateBin build_pulse_template_pair_windowed(
    const std::vector<double>& signal,
    double channelRate,
    const std::vector<size_t>& masterPeaksEcg,
    double ecgRate,
    // padSeconds IS GONE, and its absence is the fix. It was never passed to
    // extract_ppg_beats_and_align -- the only thing that decides where a beat
    // starts -- so it could not have been the lead-in it was documented as.
    // It was used in exactly one place, the [R1, R2] bracket below, where it
    // stood in for the R column that out.rCol now measures.
    // For the [pulseqc] line only. Passed rather than inferred because this
    // function has no other way to name the bin it is working on, and a
    // retention report that cannot say WHICH bin is nearly useless.
    size_t bin_index = 0,
    // This channel's confident ECG-to-channel lag in ms (channel_offset), or
    // 0. Moves the pulse window each beat's landmarks are searched in; the
    // signal itself is never shifted. See alignment::pulseLandmarks.
    double lagMs = 0.0)
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
    // The same R detections, unrounded, so the frame's R-to-50% distance is
    // the transit time on the ECG clock (see extract_ppg_beats_and_align).
    std::vector<double> peaksChExact;
    peaksChExact.reserve(masterPeaksEcg.size());
    for (size_t r : masterPeaksEcg)
        peaksChExact.push_back(static_cast<double>(r) * scale);

    const int lagSamples = (lagMs > 0.0)
        ? static_cast<int>(std::lround(lagMs * 0.001 * channelRate)) : 0;
    const auto aligned = alignment::extract_ppg_beats_and_align(
        signal, peaksCh, channelRate, &peaksChExact, lagSamples);
    if (aligned.beats.empty()) return out;

    // ---- Matched-filter QC, two-pass, per spec:
    //   (a) build a REFERENCE template as the column-wise NaN-skipping
    //       median across ALL candidate beats. The median is robust to
    //       outliers without needing to pick a fixed "seed" count.
    //   (b) score every candidate against the reference by Pearson r over
    //       the foot-to-foot window; keep beats with r >= corrFloor().
    // The final template below is then rebuilt from the survivors,
    // giving the two-pass: median-of-all -> reject low-correlation ->
    // re-median. Same shape as the ECG side, and now the same metric --
    // correlation, not the normalized error this used to use.
    // A bin where nothing passes produces NO pulse template; see the
    // zero-survivors note below for why it must not fall back to keeping
    // everything.
    std::vector<std::vector<double>> filteredBeats;
    // Diagnostic accumulators, populated inside the QC block.
    int diag_ref_col_early = 0, diag_ref_col_mid = 0, diag_ref_col_late = 0;
    int diag_ref_defined_cols = 0;
    int diag_input_beats = 0, diag_survivors = 0;
    // CORRELATIONS NOW, so the WORST pulse is the MINIMUM, not the maximum.
    // The initialisers are unchanged (+inf / -inf) because they are still
    // seeded to lose their first comparison; only the reading of them flips.
    double diag_corr_min = std::numeric_limits<double>::infinity();
    double diag_corr_max = -std::numeric_limits<double>::infinity();
    std::vector<double> diag_all_corrs;
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

        // Pearson r restricted to [f2fLo, f2fHi), non-NaN overlap only.
        // stats_utils::pearson carries NO policy -- no minimum overlap, no
        // substitute value for failure -- so both decisions are made here.
        //
        // NaN, NOT 0.0, FOR UNDEFINED: fewer than kMinCorrOverlap finite pairs,
        // or a flat window on either side. Zero is a legitimate correlation (an
        // uncorrelated pulse), so returning it for "could not measure" merges
        // two different rejections into one number and makes the diagnostic
        // median meaningless. Both still fail the gate, since NaN >= floor is
        // false. This is the opposite convention from pearsonSQI, which returns
        // 0.0 on purpose because its score is summed and thresholded.
        constexpr int kMinCorrOverlap = 8;
        auto footToFootCorr = [&](const std::vector<double>& bt) -> double {
            const int hi = std::min<int>(f2fHi,
                std::min<int>(static_cast<int>(bt.size()),
                    static_cast<int>(reference.size())));
            const PearsonResult pr =
                pearson(bt, reference, std::max(0, f2fLo), hi);
            if (pr.n_overlap < kMinCorrOverlap) return NaN;
            return pr.r;   // already NaN when undefined
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
        diag_all_corrs.reserve(aligned.beats.size());
        for (size_t k = 0; k < aligned.beats.size(); ++k) {
            const double r = footToFootCorr(aligned.beats[k]);
            diag_all_corrs.push_back(r);
            if (std::isfinite(r)) {
                if (r < diag_corr_min) diag_corr_min = r;
                if (r > diag_corr_max) diag_corr_max = r;
            }
            // >= NOT >: the floor is inclusive, and NaN fails either way.
            if (r >= pulse_qc::corrFloor()) survivorRows.push_back(k);
        }

        // ---- TRANSIT TIMES, BEFORE ANY EARLY RETURN ---------------------
        //
        // Every beat the slicer measured, survivors or not, so a bin whose
        // filter rejects everything still reports its transit times. R is the
        // ECG detection this beat was sliced from, on the ECG clock; the pulse
        // landmarks are on the pulse clock. See ptt_logging.hpp.
        {
            const size_t nb = aligned.beats.size();
            const bool haveAbs = aligned.peak_abs.size() == nb
                && aligned.up50_abs.size() == nb
                && aligned.foot_abs.size() == nb
                && aligned.maxup_abs.size() == nb
                && aligned.notch_abs.size() == nb
                && aligned.original_index.size() == nb;
            if (haveAbs) {
                std::vector<char> isKept(nb, 0);
                for (const size_t k : survivorRows) if (k < nb) isKept[k] = 1;
                out.ptt.reserve(nb);
                for (size_t k = 0; k < nb; ++k) {
                    const uint32_t sl = aligned.original_index[k];
                    if (sl >= masterPeaksEcg.size()) continue;
                    ptt_log::Beat bt;
                    bt.slice = sl;
                    bt.kept = isKept[k] != 0;
                    bt.r_peak_distance_from_binstart_in_s =
                        static_cast<double>(masterPeaksEcg[sl]) / ecgRate;
                    bt.peak_s = aligned.peak_abs[k] / channelRate;
                    bt.t50_s = aligned.up50_abs[k] / channelRate;
                    bt.foot_s = aligned.foot_abs[k] / channelRate;
                    bt.maxup_s = aligned.maxup_abs[k] / channelRate;   // NaN stays NaN
                    bt.notch_s = aligned.notch_abs[k] / channelRate;
                    out.ptt.push_back(bt);
                }
            }
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
        // One threshold now, from config.csv (pulse_qc_corr_floor), applied
        // uniformly to every bin. A bin that keeps very few pulses keeps very
        // few, and says so on the line below rather than being rescued into
        // looking fine.

        // DEGENERATE REFERENCE -> KEEP ALL, the same fallback the arterial
        // twin has at the bottom of this file. All three tiers can come back
        // empty when footToFootCorr() is below the floor or NaN for every candidate, which
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
            // diag_corr_median is declared after this block, and returning
            // early means it is never reached.
            double medHere = std::numeric_limits<double>::quiet_NaN();
            {
                std::vector<double> fin;
                fin.reserve(diag_all_corrs.size());
                for (const double e : diag_all_corrs)
                    if (std::isfinite(e)) fin.push_back(e);
                if (!fin.empty()) {
                    std::sort(fin.begin(), fin.end());
                    medHere = fin[fin.size() / 2];
                }
            }
            return out;
        }

        // ---- THE TEMPLATE'S R COLUMN: FIXED BY THE FRAME ------------------
        //
        // The aligner builds the frame from R outward (R at the median-RR
        // beat's own R column, the 50% alignment one median ptt_t50 after it),
        // so R's column is a definition, not something to measure off the
        // survivors. It stays put; the pulse is placed relative to it.
        //
        // Fallback, for a set built without it: the median of the survivors'
        // shifted R columns, clipped beats (r_cols < 0) excluded rather than
        // clamped -- a clamp to 0 would drag the median toward the pad.
        if (aligned.r_aligned_col >= 0) {
            out.rCol = aligned.r_aligned_col;
        }
        else {
            std::vector<int> rc;
            rc.reserve(survivorRows.size());
            if (aligned.r_cols.size() == aligned.beats.size())
                for (const size_t k : survivorRows)
                    if (aligned.r_cols[k] >= 0) rc.push_back(aligned.r_cols[k]);
            if (!rc.empty()) {
                std::sort(rc.begin(), rc.end());
                out.rCol = rc[rc.size() / 2];
            }
        }

        // ---- THE NOTCH, AS AN OFFSET FROM THE SYSTOLIC PEAK --------------
        //
        // NOT detected on the template: E-5 has no per-template procedure, and
        // the one it has runs on 4-second windows of the RECORDING. The
        // windowed pass already measured a notch for each of these beats.
        //
        // STORED AS THE MEDIAN (notch - peak) OFFSET, NOT AS A COLUMN, and
        // that distinction is the whole point. A column median would be in the
        // beats' own frame, while out.peakCol below is detected on the
        // TEMPLATE -- two different estimators of "the peak". Measured on MESA
        // 3014843 they differ by about 50 samples, so a notch stored as a
        // column landed 66 ms after the drawn peak when the beats themselves
        // said 155 ms. The offset is a physiological interval and survives the
        // change of reference; the column does not.
        //
        // The caller adds it to whichever peak it is drawing. notch_n is the
        // contributor count, i.e. this template's detectability: a notch
        // resting on 1 of 30 beats is as "found" as one on 30 of 30, and this
        // is the only field that says which.
        {
            std::vector<double> offs;
            offs.reserve(survivorRows.size());
            const bool usable = aligned.notch_cols.size() == aligned.beats.size()
                && aligned.peak_cols.size() == aligned.beats.size();
            if (usable)
                for (const size_t k : survivorRows)
                    if (std::isfinite(aligned.notch_cols[k]) && aligned.peak_cols[k] >= 0)
                        offs.push_back(aligned.notch_cols[k]
                            - static_cast<double>(aligned.peak_cols[k]));
            out.notch_n = static_cast<int>(offs.size());
            if (!offs.empty()) {
                std::sort(offs.begin(), offs.end());
                const std::size_t m = offs.size() / 2;
                out.notch_from_peak = (offs.size() % 2) ? offs[m]
                    : 0.5 * (offs[m - 1] + offs[m]);
            }
        }

        // Waveform and ordinal appended in the SAME loop, so they cannot fall
        // out of step.
        filteredBeats.reserve(survivorRows.size());
        out.keptSlices.reserve(survivorRows.size());
        const bool ordinalsUsable =
            aligned.original_index.size() == aligned.beats.size();
        for (const size_t k : survivorRows) {
            filteredBeats.push_back(aligned.beats[k]);
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
    // Median correlation (rank quality metric). HIGHER IS BETTER NOW -- this
    // was a median normalized error, where lower was. Anything ranking bins on
    // it has to flip its comparison.
    double diag_corr_median = std::nan("");
    if (!diag_all_corrs.empty()) {
        std::vector<double> tmp;
        tmp.reserve(diag_all_corrs.size());
        for (double e : diag_all_corrs) if (std::isfinite(e)) tmp.push_back(e);
        if (!tmp.empty()) {
            std::sort(tmp.begin(), tmp.end());
            diag_corr_median = tmp[tmp.size() / 2];
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

    // ---- Deterministic PPG fiducials from the real R-peaks -------------
    // R1 IS out.rCol, MEASURED, not padSeconds*channelRate. R2 = R1 + one RR
    // interval. We compute the systolic peak as the max in [R1, R2] (exactly
    // one pulse -> no risk of grabbing a later pulse) and the foot as the min
    // in [R1, peak].
    //
    // THE OLD BRACKET WAS OFF BY (0.5*RR - padSeconds) SECONDS, the same error
    // the viewer's time axis had, from the same false premise -- so on a slow
    // bin the window opened ~200 ms early and the "max in [R1, R2]" could
    // still be climbing the previous pulse's decay.
    //
    // NO FALLBACK. rCol < 0 means no surviving beat had a locatable R in the
    // shared frame; peakCol and footCol then stay -1, which is the state every
    // caller already reads as "no pulse fiducials for this bin". Substituting
    // a column here is what the pad was.
    if (out.rCol >= 0) {
        const int N = static_cast<int>(out.tmpl.size());
        const int r1 = std::clamp(out.rCol, 0, std::max(0, N - 1));
        // Median RR in ECG samples -> channel samples.
        std::vector<double> gaps;
        gaps.reserve(masterPeaksEcg.size());
        for (size_t k = 1; k < masterPeaksEcg.size(); ++k)
            gaps.push_back(static_cast<double>(masterPeaksEcg[k] - masterPeaksEcg[k - 1]));
        int rrCh = 0;
        if (!gaps.empty()) {
            std::sort(gaps.begin(), gaps.end());
            const double medGapEcg = gaps[gaps.size() / 2];
            rrCh = static_cast<int>(std::llround(medGapEcg * scale));
        }
        const int r2 = (rrCh > 0) ? std::min(N - 1, r1 + rrCh) : (N - 1);

        if (N > 0 && r2 > r1) {
            int pk = r1; double pmax = -std::numeric_limits<double>::infinity();
            for (int i = r1; i <= r2; ++i)
                if (!std::isnan(out.tmpl[i]) && out.tmpl[i] > pmax) {
                    pmax = out.tmpl[i]; pk = i;
                }
            int ft = r1; double fmin = std::numeric_limits<double>::infinity();
            for (int i = r1; i <= pk; ++i)
                if (!std::isnan(out.tmpl[i]) && out.tmpl[i] < fmin) {
                    fmin = out.tmpl[i]; ft = i;
                }
            out.peakCol = pk;
            out.footCol = ft;
        }
    }

    // ---- THE NOTCH'S COLUMN, RELATIVE TO THE TEMPLATE'S OWN PEAK ---------
    //
    // Now that out.peakCol is known, the median offset becomes a column in the
    // frame the template is drawn in. This is what keeps the notch bar the
    // measured interval away from the peak the operator can see.
    if (out.peakCol >= 0 && std::isfinite(out.notch_from_peak)) {
        const double c = out.peakCol + out.notch_from_peak;
        const int N2 = static_cast<int>(out.tmpl.size());
        out.notchCol = (c >= 0.0 && c <= N2 - 1) ? c : -1.0;
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
    out.std = (out.footCol >= 0)
        ? normalize_features::local_ratio_std(beatsForTemplate, out.footCol)
        : std::vector<double>(maxLen, 0.0);

    // Retain aligned per-beat slices for downstream (snips CSV, etc).
    out.kept.reserve(beatsForTemplate.size());
    for (const auto& sl : beatsForTemplate) out.kept.push_back(sl);

    (void)channelRate; (void)ecgRate;

    // THE SUCCESS PATH HAD NO RETURN. Every early exit above returns `out`
    // explicitly; the one path that builds a template fell off the end of a
    // non-void function. MSVC reports that as C4715 -- a WARNING -- so it
    // shipped. The caller got a return object that kept the two ints and lost
    // all four vectors, which is why every bin reported a plausible
    // peakCol/footCol with an EMPTY waveform and the pulse channel vanished.
    //
    // This function returned void and wrote through out-parameters before it
    // was refactored to return by value; that shape had no return to lose.
    return out;
}

/**
 * @brief  PPG templates for every bin. Slicing is driven by bin.ch1.raw
 *         (ECG-frame R-peaks), consistent with the ECG templater.
 *
 * @param bins        Input bins.
 * @param ecgRate     ECG sample rate (for R-peak time base).
 * @param ppgRate     PPG sample rate.
 *
 * NO padSeconds PARAMETER. It documented a fixed lead-in that the slicer
 * never had (see build_pulse_template_pair_windowed), and the doc value
 * (0.25), the default (0.4) and bin_plot_widget's copy of it (0.3, then 0.4)
 * had all drifted apart -- three numbers for a quantity that is measured per
 * bin and reported as PulseTemplateBin::rCol. Every call site already omitted
 * the argument, so removing it is source-compatible.
 */
inline PPGTemplatesResult CreatePulseTemplates(
    const std::vector<output_binfile_data>& bins,
    std::vector<double> output_binfile_data::* sigMember,
    double ecgRate,
    double channelRate,
    // This channel's pulse-window lag; see build_pulse_template_pair_windowed.
    double lagMs = 0.0)
{
    size_t n = bins.size();
    PPGTemplatesResult out;
    out.templates.assign(n, {});
    out.iqrs.assign(n, {});
    out.kept.assign(n, {});
    out.keptSlices.assign(n, {});
    out.ptt.assign(n, {});
    out.peakCol.assign(n, -1);
    out.notch_from_peak.assign(n, std::numeric_limits<double>::quiet_NaN());
    out.notchCol.assign(n, -1.0);
    out.notch_n.assign(n, 0);
    out.footCol.assign(n, -1);
    out.rCol.assign(n, -1);

    if (channelRate <= 0.0) return out;   // channel absent from this dataset

    int ppg_threads = std::min(8, static_cast<int>(n > 0 ? n : 1));
#pragma omp parallel for schedule(dynamic) num_threads(ppg_threads)
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const auto& b = bins[i];
        if ((b.*sigMember).empty() || b.ch1.raw.size() < 2)
            continue;
        try {
            PulseTemplateBin r = build_pulse_template_pair_windowed(
                b.*sigMember, channelRate, b.ch1.raw, ecgRate, i, lagMs);
            out.templates[i] = std::move(r.tmpl);
            out.iqrs[i] = std::move(r.std);
            out.kept[i] = std::move(r.kept);
            out.keptSlices[i] = std::move(r.keptSlices);
            out.ptt[i] = std::move(r.ptt);
            out.peakCol[i] = r.peakCol;
            out.notch_from_peak[i] = r.notch_from_peak;
            out.notchCol[i] = r.notchCol;
            out.notch_n[i] = r.notch_n;
            out.footCol[i] = r.footCol;
            out.rCol[i] = r.rCol;
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