#pragma once
/*
The ECG Normalization algorithm is as follows:

    1. Calculate the total QRS vector magnitude for each beat to control for cardiac axis rotation:
        RS_peak(t) = abs(R_peak(t)) + abs(S_peak(t))
    2. Find the global reference by taking the median across all bins for each individual:
        Global_Ref_person = median(RS_peak(t))
    3. Normalize any amplitude feature like P, R, or T waves using the equation:
        Feature_peak_norm_abs = Feature_peak(t) / Global_Ref_person


The PPG Normalization algorithm is as follows:
    1. First, calculate the local PI for each beat:
        PI(t) = ((systolic_peak(t) - diastolic_trough(t)) / abs(diastolic_trough(t))) * 100.
    2. Find the global reference by taking the median PI across the entire recording for that individual:
        Global_Ref_person = median(PI(t)).
    3. Normalize your amplitude feature by converting it to its local baseline ratio first and then
       dividing by the global reference:
        Feature_peak_norm_abs = Feature_Local_Ratio(t) / Global_Ref_person.

*/

#include "fiducial_marker_finding\template_marking_bin_io.hpp"
#include "fiducial_marker_finding\global_intervals.hpp"
#include "fiducial_marker_finding\vcg_signal_average.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include <utility>
#include <string>
#include <fstream>

namespace normalize_features {

    inline double median_finite(std::vector<double> v) {
        // Remove NaN and Inf values, then compute the median of the remaining finite values.
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::nan("");
        const size_t mid = v.size() / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        const double a = v[mid];
        if (v.size() % 2 == 1) return a;
        // Even count: pair with max of the lower half.
        double b = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < mid; ++i) if (v[i] > b) b = v[i];
        return 0.5 * (a + b);
    }

    inline double sample_y(const std::vector<double>& v, double idx) {
        // Amplitude at a SUB-SAMPLE marker position, NaN if unavailable.
        // INTERPOLATED, not subscripted: every marker reaching this is a double
        // now, and rounding here would undo the widening. Callers passing a
        // whole column (a loop index, an integer footIdx) convert implicitly
        // and get that exact sample back.
        if (idx < 0.0 || idx > static_cast<double>(v.size()) - 1.0) return std::nan("");
        const double y = FeatureMarks::sample_at(v, idx);
        return std::isnan(y) ? std::nan("") : y;
    }

    inline double ecg_norm(double raw, double ref) {
        //ref = median over bins of |R| + |S|
        if (!std::isfinite(ref) || ref == 0.0 || std::isnan(raw)) return raw;
        return raw / ref;
    }

    inline double calculate_perfusion_index(double y, double foot_y) {
        //PI = 100 * (y - foot) / |foot|
        if (std::isnan(y) || std::isnan(foot_y) || std::abs(foot_y) < 1e-12)
            return std::nan("");
        return 100.0 * (y - foot_y) / std::abs(foot_y);
    }

    inline double pulse_norm(double y, double foot_y, double ref) {
        //divide PPG by the global reference (PI)
        if (!std::isfinite(ref) || ref == 0.0) return std::nan("");
        const double lr = calculate_perfusion_index(y, foot_y);
        return std::isnan(lr) ? lr : lr / std::abs(ref);
    }


    inline const tbank::template_of_all_signals* ecg_slot(const time_bin& b, int ch, int slot) {
        //given a bin, lead, and template slot number, get pointer to template
        if (ch < 0 || ch >= 3 || slot < 0) return nullptr;
        const auto& t = b.ecg_bank[ch].templates;
        return (static_cast<size_t>(slot) < t.size()) ? &t[slot] : nullptr;
    }

    inline double qrs_height_for_template(const time_bin& b, int ch, int slot, double sampleRateHz)
    {
        //given a bin, lead, and template slot number, get qrs complex height (|R|+|S|)
        if (b.bad_segment || b.bad_r_ch[ch]) return std::nan("");
        const tbank::template_of_all_signals* tp = ecg_slot(b, ch, slot);
        if (!tp || tp->tmpl.empty()) return std::nan("");
        const std::vector<double>& ecg = tp->tmpl;

        const tbank::BankMarkerSet& rmk = b.slotMarks(ch, slot, AnchorType::R_PEAK);
        const double rIdx = (rmk.r_peak_auto >= 0.0)
            ? rmk.r_peak_auto : static_cast<double>(tp->r_col);
        if (rIdx < 0.0) return std::nan("");
        // (reactive_ecg removed: it existed to supply computeEcgFeatures' old
        //  p_peak parameter, which the body never read. |R|+|S| is all this
        //  function wants from f.)
        EcgFeatures f = computeEcgFeatures(ecg, rmk.q_onset, rIdx, rmk.s_end, rmk.t_end, sampleRateHz, b.polarity.sign(ch));
        const double ry = sample_y(ecg, f.r_idx);
        const double sy = sample_y(ecg, f.s_idx);
        if (std::isnan(ry) || std::isnan(sy)) return std::nan("");
        return std::abs(ry) + std::abs(sy);
    }

    inline double compute_ecg_global_ref(const std::vector<time_bin>& bins, int ch, int slot, double sampleRateHz) {
        //median of all QRS heights for a given lead - called global ref in scientific literature
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) vals.push_back(qrs_height_for_template(b, ch, slot, sampleRateHz));
        return median_finite(std::move(vals));
    }

    // ------------------------------------------------------------------
    // Pulse channel accessor (0=PPG, 1=ABP, 2=ART, 3=ART_PULM).
    // ------------------------------------------------------------------
    struct PulseChannel {
        const std::vector<double>* trace;
        // SUB-SAMPLE POSITIONS, matching TemplateBin's widened pulse fields.
        // pulseChan() brace-initialises these, and brace init refuses to
        // narrow -- which is why one int here produced forty-odd errors.
        double foot_idx;
        double peak_idx;
        uint8_t issue;   // 0 = ok, 1 = user-bad, 2 = absent
        // Added for the area reference below, which needs the far bracket of
        // the wave and not just its peak. dicrotic_idx is the systolic/
        // diastolic divide, end_idx the end of the wave; both are -1 on
        // channels or bins where the notch was not found, which the area
        // reference treats as "fall back to end" and then "skip this bin".
        double dicrotic_idx;
        double peak2_idx;
        double end_idx;
    };

    // ---- THE PPG SLOT THIS BIN IS REPRESENTED BY -------------------------
    //
    // The pulse bars are per slot, so a bin-level perfusion index needs a slot
    // chosen. THE DOMINANT MORPHOLOGY: most surviving members, skipping
    // empties and thin cohorts. PI is a physiological amplitude feeding a
    // subject-wide median, so the representative pulse is wanted -- and slot
    // indices are not ordered by population, so taking slot 0 could let a rare
    // ectopic column set the reference.
    //
    // -1 when the bin has no usable pulse slot, which compute_pulse_global_ref
    // skips exactly as it skipped an empty trace.
    inline int dominantPpgSlot(const time_bin& b) {
        int best = -1; int bestN = 0;
        for (int t = 0; t < b.ppg_bank.size(); ++t) {
            const tbank::template_of_all_signals& ps = b.ppg_bank.templates[t];
            if (ps.tmpl.empty()) continue;
            if (ps.tooFewBeats(/*is_ppg=*/true)) continue;
            const int n = ps.cleanCount();
            if (n > bestN) { bestN = n; best = t; }
        }
        return best;
    }

    inline PulseChannel pulseChan(const time_bin& b, int which) {
        switch (which) {
        case 0: {
            // THE SLOT'S OWN WAVEFORM WITH THE SLOT'S OWN BARS. Pairing a
            // per-slot bar with b.ppgTemplate -- the bin-wide average -- would
            // measure a foot found on one waveform against a different one,
            // which is the error the bars moved per slot to prevent.
            //
            // peak2 is DERIVED: BankPulseMarkerSet carries the three bars and
            // the five detector columns, and peak2 comes back from
            // reactive_ppg bracketed by them.
            const int t = dominantPpgSlot(b);
            if (t < 0)
                return { &b.ppgTemplate, -1.0, -1.0, uint8_t(2), -1.0, -1.0, -1.0 };
            const tbank::template_of_all_signals& ps = b.ppg_bank.templates[t];
            const tbank::BankPulseMarkerSet& pm = ps.pulse_marks;
            const FeatureMarks::ReactivePpg rp = FeatureMarks::update_ppg_markings(
                ps.tmpl, pm.onset, pm.peak_auto, pm.dicrotic, pm.end);
            const uint8_t issue = ps.badPulseMarked() ? uint8_t(1) : b.bad_ppg;
            return { &ps.tmpl, pm.onset, pm.peak_auto, issue,
                     pm.dicrotic, rp.peak2, pm.end };
        }
        case 1: return { &b.abpTemplate,     b.abp_onset,    b.abp_peak,    b.abp_issue,
                         b.abp_dicrotic,     b.abp_peak2,    b.abp_end };
        case 2: return { &b.artTemplate,     b.art_onset,    b.art_peak,    b.art_issue,
                         b.art_dicrotic,     b.art_peak2,    b.art_end };
        default: return { &b.artPulmTemplate, b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_issue,
                         b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };
        }
    }

    inline constexpr int kNumPulseCh = 4;

    inline double compute_pulse_global_ref(const std::vector<time_bin>& bins, int which)
    {
        //PI(bin) = 100 * (peak_y - foot_y) / |foot_y|
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            const PulseChannel pc = pulseChan(b, which);
            if (pc.issue != 0) continue;
            if (pc.trace->empty()) continue;
            const double foot_y = sample_y(*pc.trace, pc.foot_idx);
            const double peak_y = sample_y(*pc.trace, pc.peak_idx);
            vals.push_back(calculate_perfusion_index(peak_y, foot_y));
        }
        return median_finite(std::move(vals));
    }

    // ------------------------------------------------------------------
    // Whole-trace normalization. These are the ONLY place a raw ECG/pulse
    // trace (mean template, individual beat, or a precomputed spread like
    // IQR) should be converted to normalized units -- callers (viewer,
    // CSV export, anywhere else) must call these rather than reimplement
    // the divide/ratio math locally.
    // ------------------------------------------------------------------

    // Divide every sample of `raw` by `ref`, with the same guards as
    // ecg_norm. This is also the correct final step for pulse channels:
    // once a trace is already in "local ratio" units (see
    // local_ratio_iqr below), dividing by Global_Ref_person is a plain
    // scalar divide, identical in form to the ECG step.
    inline std::vector<double> scale_array_by_ref(const std::vector<double>& raw, double ref) {
        std::vector<double> out(raw.size());
        for (size_t i = 0; i < raw.size(); ++i) out[i] = ecg_norm(raw[i], ref);
        return out;
    }
    // Spread (IQR, sd) in RAW pulse units -> normalized units.


    // ==================================================================
    // PULSE TRACE NORMALIZATION: FULL TRANSFORM OR NOT SHOWN
    // ==================================================================
    //
    // A pulse trace is drawn only in normalized units:
    //     100 * (y - foot) / |foot| / |ref|
    // There is NO display fallback. When the foot cannot be found (or sits at
    // zero, so it cannot be divided by) or the reference is unusable, the trace
    // is NOT SHOWN -- these functions return an EMPTY vector.
    //
    // EMPTY, NEVER A VECTOR OF NaN. A NaN-filled trace is what caused the old
    // invisible failure: the draw skipped every sample, the axis fell back to
    // 0..1, and drawFeatureGlyphs put the pulse fiducials on the axis FLOOR, so
    // an empty panel looked like a real detection. An empty result gives the
    // caller something to test: it must treat it exactly like hasPPG == false
    // -- no trace, no band, and NO GLYPHS.
    //
    // NOTE THE REFERENCE IS ONE OF THE CONDITIONS. ref comes from
    // compute_pulse_global_ref, which reads the BIN-level ppgTemplate rather
    // than the ppg_bank slot being drawn; when that comes back NaN, every pulse
    // panel in the record is hidden, not just the ones missing a foot.
    //
    // The band uses the same test and the same scale (without the shift, since
    // a spread is a difference), so a band is shown exactly when its trace is.
    inline bool pulse_trace_normalizable(double foot_y, double ref) {
        return !std::isnan(foot_y) && std::abs(foot_y) >= 1e-12
            && std::isfinite(ref) && ref != 0.0;
    }

    inline std::vector<double> scale_pulse_spread_by_ref(const std::vector<double>& raw, double foot_y, double ref) {
        // The spread is in raw units; scale it by the trace's factor,
        // 100 / |foot| / |ref|. Empty when the trace itself is not shown.
        if (!pulse_trace_normalizable(foot_y, ref)) return {};
        const double k = 100.0 / std::abs(foot_y) / std::abs(ref);
        std::vector<double> out(raw.size());
        for (size_t i = 0; i < raw.size(); ++i)
            out[i] = std::isnan(raw[i]) ? raw[i] : raw[i] * k;
        return out;
    }
    inline std::vector<double> normalize_ecg_trace(const std::vector<double>& raw, double ref) {
        return scale_array_by_ref(raw, ref);
    }

    // Pulse: local ratio (per-sample, using THIS trace's own foot) then
    // divide by ref. Works for the mean template or any individual beat --
    // never uses a median/global foot value, per the documented algorithm.
    // EMPTY when the foot or the reference is unusable: not shown (see above).
    //
    // footIdx IS A DOUBLE. It was an int, and every pulse foot reaching it is
    // sub-sample -- BankPulseMarkerSet::onset is a double, and
    // normalize_ppg_or_similar's own parameter is a double -- so the int
    // narrowed the foot to a whole column twice on the way in (once at the
    // viewer's `int ppgFootIdx`, once here) and then handed it to sample_y,
    // which interpolates. Two truncations to reach a function that did not
    // need either.
    inline std::vector<double> normalize_pulse_trace(const std::vector<double>& raw, double footIdx, double ref) {
        const double foot_y = sample_y(raw, footIdx);
        if (!pulse_trace_normalizable(foot_y, ref)) return {};
        std::vector<double> out(raw.size());
        for (size_t i = 0; i < raw.size(); ++i) out[i] = pulse_norm(raw[i], foot_y, ref);
        return out;
    }


    // ------------------------------------------------------------------
    // Cross-beat spread helpers, computed once at template-build time from
    // the raw aligned beats -- NOT from individual beats retained downstream
    // (that overlay-beat machinery has been removed; these summary
    // statistics are all that's kept).
    //
    // EVERY SPREAD IN THIS PIPELINE IS A PER-SAMPLE STD (ddof = 1). There is no
    // interquartile range left anywhere except keep_within_tukey's own fence,
    // which is a different thing entirely. The template_bank and ppg_realign
    // producers that used to write q3 - q1 into the same fields now write an SD
    // too, so a field's contents no longer depend on which code path filled it.
    // ------------------------------------------------------------------

    // ECG: pulse_norm-equivalent step is a plain scalar divide, so taking
    // the spread of raw amplitudes and dividing by ref later
    // (scale_array_by_ref) is exact -- no restructuring needed relative to
    // the raw computation.
    // (Kept here only as a named entry point so build-time code doesn't
    // need to hand-roll the loop.)
    inline std::vector<double> raw_amplitude_iqr(const std::vector<std::vector<double>>& rawBeats) {
        if (rawBeats.empty()) return {};
        size_t maxLen = 0;
        for (const auto& bt : rawBeats) maxLen = std::max(maxLen, bt.size());
        std::vector<double> sd(maxLen, 0.0);
        std::vector<double> col;
        col.reserve(rawBeats.size());
        for (size_t c = 0; c < maxLen; ++c) {
            col.clear();
            for (const auto& bt : rawBeats)
                if (c < bt.size() && !std::isnan(bt[c])) col.push_back(bt[c]);
            const size_t n = col.size();
            if (n < 2) continue;
            double mean = 0.0;
            for (double v : col) mean += v;
            mean /= static_cast<double>(n);
            double sumsq = 0.0;
            for (double v : col) sumsq += (v - mean) * (v - mean);
            sd[c] = std::sqrt(sumsq / static_cast<double>(n - 1));   // ddof = 1
        }
        return sd;
    }

    // Pulse: unlike ECG, the per-sample transform's slope varies beat-to-
    // beat (each beat has its own foot_y), so taking the spread of raw
    // values and dividing by a single factor afterward is NOT equivalent to
    // the documented algorithm. Convert each beat to its own local-ratio
    // trace FIRST (own foot, no global/median foot), take the cross-beat
    // spread of that, and defer only the final /Global_Ref_person to
    // scale_array_by_ref() at display/export time -- exactly mirroring how
    // the ECG spread defers its /ref step.
    inline std::vector<double> local_ratio_iqr(const std::vector<std::vector<double>>& rawBeats, int footIdx) {
        if (rawBeats.empty()) return {};
        std::vector<std::vector<double>> ratioBeats;
        ratioBeats.reserve(rawBeats.size());
        for (const auto& bt : rawBeats) {
            const double footY = sample_y(bt, footIdx);
            std::vector<double> r(bt.size());
            for (size_t i = 0; i < bt.size(); ++i) r[i] = calculate_perfusion_index(bt[i], footY);
            ratioBeats.push_back(std::move(r));
        }
        return raw_amplitude_iqr(ratioBeats);   // same cross-beat STD mechanics, different input units
    }

    // ==================================================================
    // Length / area / volume triad
    // ==================================================================
    // Three summary statistics over one feature's sample window [lo, hi],
    // all in raw sample-index units (no fs or amplitude-scale conversion
    // applied here -- same "defer the unit conversion to the caller"
    // convention as scale_array_by_ref / the *_iqr fields above). Multiply
    // by 1/fs and/or an amplitude scale afterward if physical units are
    // needed.
    //
    // Placed here, ahead of Section 5.2, because Option B/C below call
    // segment_area/segment_volume directly -- C++ has no forward
    // declaration for free functions used before their definition in the
    // same translation unit, so these must come first textually.

    // Curve length: cumulative Euclidean distance between consecutive
    // samples, one sample-index unit of run per step. A NaN sample breaks
    // the run at that step (a gap contributes nothing, rather than a
    // phantom straight line jumping across it).
    // SUB-SAMPLE BOUNDS. Every landmark that reaches these integrators is a
    // fractional double (BankMarkerSet, FeatureMarks, TemplateBin's pulse
    // fields), so the window they describe has fractional ends and rounding
    // them here would undo the widening -- a half-sample error on each end of
    // a 40-sample QRS is a 2.5% error in the area that feeds
    // Global_Ref_person. The sample interval containing each end contributes
    // only its overlapped fraction. Whole-number bounds reduce to the old
    // whole-sample sum exactly, so integer callers are unaffected.
    //
    // NOT VIA sample_y, deliberately. Each endpoint is interpolated from the
    // two samples of the interval being integrated, which this loop has
    // already checked are both non-NaN. sample_y would read the endpoint as a
    // position in the ARRAY: at a whole-number e it brackets [e, e+1], so a
    // NaN one sample PAST the window returns NaN for a sample that is
    // perfectly good, and that NaN then propagates through the whole sum
    // instead of being skipped as a gap.
    inline double segment_length(const std::vector<double>& v, double lo, double hi) {
        lo = std::max(0.0, lo);
        hi = std::min(hi, static_cast<double>(v.size()) - 1.0);
        if (!(hi > lo)) return std::nan("");
        double len = 0.0;
        bool any = false;
        const int first = static_cast<int>(std::floor(lo));
        const int last = static_cast<int>(std::ceil(hi)) - 1;
        for (int i = first; i <= last; ++i) {
            const double a = v[i], b = v[i + 1];
            if (std::isnan(a) || std::isnan(b)) continue;   // gap: no phantom line across it
            const double s = std::max(lo, static_cast<double>(i));
            const double e = std::min(hi, static_cast<double>(i) + 1.0);
            if (!(e > s)) continue;
            const double slope = b - a;
            const double dy = (e - s) * slope;             // (a + (e-i)*slope) - (a + (s-i)*slope)
            len += std::sqrt((e - s) * (e - s) + dy * dy);
            any = true;
        }
        return any ? len : std::nan("");
    }

    // Trapezoidal area under v over [lo, hi]. `absolute` rectifies before
    // integrating -- the conventional way to report QRS/T-wave area, since a
    // biphasic complex would otherwise partially cancel itself in a signed
    // integral. NaN samples are skipped (that trapezoid contributes nothing,
    // rather than propagating NaN across the whole sum). Sub-sample bounds,
    // as segment_length above -- see the note there.
    //
    // `baseline` is subtracted from each endpoint INSIDE the loop, which is
    // what lets a caller integrate a foot-zeroed pulse without materializing
    // a zeroed copy of the trace. It cannot be done by subtracting
    // baseline*(hi-lo) from the result afterward: NaN trapezoids are skipped,
    // so the width actually integrated is not hi-lo whenever the window
    // contains a gap, and the correction would over-subtract by exactly the
    // skipped width.
    inline double segment_area(const std::vector<double>& v, double lo, double hi,
        bool absolute = true, double baseline = 0.0) {
        lo = std::max(0.0, lo);
        hi = std::min(hi, static_cast<double>(v.size()) - 1.0);
        if (!(hi > lo)) return std::nan("");
        double area = 0.0;
        bool any = false;
        const int first = static_cast<int>(std::floor(lo));
        const int last = static_cast<int>(std::ceil(hi)) - 1;
        for (int i = first; i <= last; ++i) {
            const double v0 = v[i], v1 = v[i + 1];
            if (std::isnan(v0) || std::isnan(v1)) continue;
            const double s = std::max(lo, static_cast<double>(i));
            const double e = std::min(hi, static_cast<double>(i) + 1.0);
            if (!(e > s)) continue;
            const double slope = v1 - v0;
            double a = v0 + (s - static_cast<double>(i)) * slope - baseline;
            double b = v0 + (e - static_cast<double>(i)) * slope - baseline;
            // Rectified at the ENDPOINTS, as before: a trapezoid straddling
            // zero is still approximated by |f| at its two ends rather than
            // split at the crossing. Unchanged convention, stated because
            // fractional ends make it easier to mistake for exact.
            if (absolute) { a = std::abs(a); b = std::abs(b); }
            area += 0.5 * (a + b) * (e - s);   // trapezoid, width e-s
            any = true;
        }
        return any ? area : std::nan("");
    }

    // Spatial "volume": trapezoidal integral, over [lo, hi], of the 3-lead
    // vector magnitude sqrt(ch1^2+ch2^2+ch3^2) -- the 3-D analogue of
    // segment_area, one level up from Option C's peak-magnitude reference.
    // ch1/ch2/ch3 MUST already be on the shared R-relative axis (same
    // length, same offset origin) before calling this -- exactly the axis
    // vcg_signal_average.hpp's loopFromTemplates/perBeatLoops already
    // produce. This function does not align them; it only integrates.
    inline double segment_volume(const std::vector<double>& ch1, const std::vector<double>& ch2,
        const std::vector<double>& ch3, double lo, double hi) {
        const double n = static_cast<double>(std::min({ ch1.size(), ch2.size(), ch3.size() }));
        lo = std::max(0.0, lo);
        hi = std::min(hi, n - 1.0);
        if (!(hi > lo)) return std::nan("");
        // MAGNITUDE OF THE INTERPOLATED SAMPLES, not an interpolation of the
        // magnitude: each lead is interpolated at the same sub-sample position
        // within the interval and the vector is formed there, which is the
        // same order of operations the whole-sample version used. Interval-
        // local, for the reason given on segment_length.
        double vol = 0.0;
        bool any = false;
        const int first = static_cast<int>(std::floor(lo));
        const int last = static_cast<int>(std::ceil(hi)) - 1;
        for (int i = first; i <= last; ++i) {
            const double x0 = ch1[i], x1 = ch1[i + 1];
            const double y0 = ch2[i], y1 = ch2[i + 1];
            const double z0 = ch3[i], z1 = ch3[i + 1];
            if (std::isnan(x0) || std::isnan(x1) || std::isnan(y0)
                || std::isnan(y1) || std::isnan(z0) || std::isnan(z1)) continue;
            const double s = std::max(lo, static_cast<double>(i));
            const double e = std::min(hi, static_cast<double>(i) + 1.0);
            if (!(e > s)) continue;
            auto mag = [&](double t) -> double {
                const double f = t - static_cast<double>(i);
                const double x = x0 + f * (x1 - x0);
                const double y = y0 + f * (y1 - y0);
                const double z = z0 + f * (z1 - z0);
                return std::sqrt(x * x + y * y + z * z);
                };
            vol += 0.5 * (mag(s) + mag(e)) * (e - s);
            any = true;
        }
        return any ? vol : std::nan("");
    }

    // ==================================================================
    // Section 5.2 -- Global reference, Options A/B/C, and the CV check
    // ==================================================================
    //
    // Three ways to reduce beats to a single reference scalar. A and B are
    // PER SLOT -- each bank slot is its own reference, measured on its own
    // template -- and median across bins for that slot. C still fuses the
    // bin-level templates and is not yet per slot:
    //
    //   A (existing, above) : median(|R_peak| + |S_peak|)      -- two samples
    //   B (below)           : median(QRS area, Q-onset..J-point) -- integrates
    //                         the whole complex, so a wide-but-modest QRS and
    //                         a narrow-but-tall one sharing |R|+|S| are no
    //                         longer equivalent.
    //   C (below)           : median(peak spatial vector magnitude) -- fuses
    //                         all three ECG leads into one 3-D vector first
    //                         (reusing the SAME cross-channel R-relative
    //                         alignment vcg_signal_average.hpp already solves
    //                         for the VCG loop -- see global_intervals.hpp's
    //                         "ALIGNMENT" note on why raw column indices from
    //                         different channels cannot be combined directly),
    //                         so cardiac-axis rotation is controlled for by
    //                         geometry rather than by a |R|+|S| proxy.
    //
    // ASSUMPTION (flagged, not silently decided): the spec names Options B
    // and C without defining their measurement window. B integrates
    // Q-onset..J-point (the QRS complex) because that is the same window
    // Option A samples from (R and S both fall inside it). C's window is a
    // fixed pre/post sample margin around R (preSamples/postSamples,
    // defaulted below), because the spatial loop needs a window before any
    // per-bin QRS onset/offset can be measured FROM it (global_intervals.hpp
    // reduces per-lead onsets that are not yet known when the vector is being
    // built for that measurement's own reference). Widen the defaults if a
    // program's QRS is unusually broad.

    // Option B: area-based reference for ONE channel and ONE slot. Same rule as
    // Option A -- each slot is its own reference, measured on its own template
    // with its own R-pass markers -- swapping the |R|+|S| reduction for the
    // QRS's rectified area (Q-onset -> J-point / s_end).
    inline double slot_qrs_area(const time_bin& b, int ch, int slot)
    {
        if (b.bad_segment || b.bad_r_ch[ch]) return std::nan("");
        const tbank::template_of_all_signals* tp = ecg_slot(b, ch, slot);
        if (!tp || tp->tmpl.empty()) return std::nan("");
        const tbank::BankMarkerSet& rmk = b.slotMarks(ch, slot, AnchorType::R_PEAK);
        // STRAIGHT THROUGH AS DOUBLES. These are BankMarkerSet's
        // sub-sample landmarks and segment_area integrates over
        // fractional bounds, so there is nothing to round; -1 means
        // absent, which the qBegin < 0 test below rejects.
        const double qBegin = rmk.q_onset;
        const double jPoint = rmk.s_end;   // S_END == J_POINT (AnchorType comment)
        if (qBegin < 0.0 || jPoint <= qBegin) return std::nan("");
        return segment_area(tp->tmpl, qBegin, jPoint, /*absolute=*/true);
    }

    // Median over bins of slot `slot`'s own QRS area. Same caveat and same
    // required-slot signature as compute_ecg_global_ref.
    inline double compute_ecg_global_ref_area(const std::vector<time_bin>& bins, int ch, int slot,
        double sampleRateHz)
    {
        (void)sampleRateHz;   // kept for signature parity with Option A
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) vals.push_back(slot_qrs_area(b, ch, slot));
        return median_finite(std::move(vals));
    }

    // Option C: spatial vector-magnitude Global_Ref_person, fusing all three
    // ECG channels. Builds the R-relative 3-lead loop the SAME way
    // vcg_signal_average.hpp's save-time path does (each channel read at ITS
    // OWN r_col + offset -- see vcg_signal_average.hpp's "AXIS" note), so this
    // does not re-derive cross-channel alignment; it reuses the one already
    // proven for the VCG loop. Global_Ref_person is the median, across bins,
    // of each bin's peak spatial magnitude sqrt(x^2+y^2+z^2).
    inline double compute_ecg_global_ref_spatial(const std::vector<time_bin>& bins,
        int preSamples = 40, int postSamples = 60)
    {
        std::vector<double> peaks;
        peaks.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            const vcg_avg::Loop loop = vcg_avg::loopFromTemplates(b, preSamples, postSamples);
            if (loop.pts.empty()) continue;
            double peak = 0.0;
            bool any = false;
            for (const auto& p : loop.pts) {
                if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) continue;
                const double mag = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
                if (mag > peak) peak = mag;
                any = true;
            }
            if (any) peaks.push_back(peak);
        }
        return median_finite(std::move(peaks));
    }


    // Median absolute deviation, NaN-skipping, matching median_finite's
    // convention (used only by cvFlag below, so kept local to this file
    // rather than promoted to stats_utils.hpp).
    inline double mad_of(const std::vector<double>& x) {
        const double m = median_finite(x);
        if (std::isnan(m)) return std::nan("");
        std::vector<double> absdev;
        absdev.reserve(x.size());
        for (double v : x) if (!std::isnan(v)) absdev.push_back(std::abs(v - m));
        return median_finite(std::move(absdev));
    }

    // CV check (5.2): flags a subject/bin whose QRS-reference values are too
    // dispersed relative to their own Global_Ref_person to trust the ratio
    // normalization below -- CV = MAD / Global_Ref_person, flagged above
    // 0.15. false (not flagged) when gref is unusable, since there is then
    // nothing to compare the dispersion against.
    inline bool cv_flag(const std::vector<double>& qrsRef, double gref) {
        if (!(gref > 0.0) || std::isnan(gref)) return false;
        const double m = mad_of(qrsRef);
        if (std::isnan(m)) return false;
        return (m / gref) > 0.15;
    }

    // ==================================================================
    // Section 5.3 -- Ratio normalization
    // ==================================================================
    // Feature_peak_norm_abs = Feature_peak / Global_Ref_person. This is
    // exactly ecg_norm's single-value form (same guards: an unusable ref or
    // a NaN feature passes the raw value through unchanged rather than
    // dividing by something meaningless) -- named separately here because
    // Section 5.3 refers to it as its own step, applied to whichever
    // Global_Ref_person Option A/B/C above produced.
    inline double ratio_norm(double featurePeak, double gref) { return ecg_norm(featurePeak, gref); }



    inline double pct_scale(double ratio, double p2, double p98) {
        //Places ratio between p2 and p98, turns that position into a score from 0 to 100 (p2 gives 0, p98 gives 100), and caps anything outside that range at 0 or 100.
        const double range = p98 - p2;
        if (std::isnan(ratio) || !(range > 0.0)) return std::nan("");
        return std::clamp((ratio - p2) / range * 100.0, 0.0, 100.0);
    }


    // ==================================================================
    // Heart-rate-proportional beat segmentation, PQ-zeroed
    // ==================================================================
    // Distinct from alignment.hpp's extract_beats_and_align: that function
    // slices at FIXED proportions (0.3 RR before / 1.5 RR after, see its
    // percent_interval_preceeding_rpeak / percent_interval_following_rpeak
    // constants) baked in for template building, and prefers the TP segment
    // over PQ for its two-stage DC leveling. This slicer is a separate,
    // purpose-built segmenter for the length/area/volume feature work above:
    // 0.25 RR before R / 0.75 RR after (per spec), PQ ONLY as the vertical
    // zero (never TP), and an explicit minimum-yield gate the template
    // slicer does not have.
    //
    // PQ baseline reuses FeatureMarks' own P/Q detectors directly (P-end via
    // seed_p_peak + detect_p_end, Q-onset via compute_q_onset) rather than
    // re-deriving isoelectric detection -- the same P-end -> Q-onset window
    // alignment.hpp's Stage-2 PQ leveling comment describes.
    struct ProportionalBeat {
        std::vector<double> samples;   // R at column rCol; PQ-zeroed when pqBaseline is not NaN
        int    rCol = -1;
        int    rrLen = -1;             // this beat's own RR, in samples
        // THIS BEAT'S CHANNEL POLARITY, from LeadPolarity::sign(lead), stamped
        // by the slicer. Carried on the beat rather than threaded through
        // qrs_window_of / build_feature_time_series / _3ch, because those three
        // are handed a beat with no channel index -- and a ProportionalBeat
        // already carries its own rCol and rrLen, so its own sign belongs here
        // too. Defaults to upright so an un-stamped beat behaves as before.
        double sgn = 1.0;
        double pqBaseline = std::numeric_limits<double>::quiet_NaN();   // subtracted DC level; NaN if PQ unavailable
    };

    struct ProportionalBeatSet {
        std::vector<ProportionalBeat> beats;
        int  nExpected = 0;     // record duration / the record's own median RR, +1
        bool sufficient = false;   // beats.size() >= 0.5 * nExpected (the "at least 50%" gate)
    };

    inline ProportionalBeatSet segment_beats_proportional(
        const std::vector<double>& ecg, const std::vector<size_t>& rPeaks, double fs,
        double sgn,
        double beforeFrac = 0.25, double afterFrac = 0.75)
    {
        ProportionalBeatSet out;
        const int64_t N = static_cast<int64_t>(ecg.size());
        if (N == 0 || rPeaks.size() < 2 || !(fs > 0.0)) return out;

        // Expected beat count from the record's own median RR -- the same
        // "one number per record" role median_finite plays everywhere else
        // in this file, just over RR instead of an amplitude/ratio.
        std::vector<double> rrAll;
        rrAll.reserve(rPeaks.size() - 1);
        for (size_t i = 0; i + 1 < rPeaks.size(); ++i)
            rrAll.push_back(static_cast<double>(rPeaks[i + 1] - rPeaks[i]));
        const double medRR = median_finite(rrAll);
        if (!(medRR > 0.0)) return out;
        const double durationSamples =
            static_cast<double>(rPeaks.back() - rPeaks.front());
        out.nExpected = static_cast<int>(std::lround(durationSamples / medRR)) + 1;

        out.beats.reserve(rPeaks.size());
        for (size_t i = 0; i < rPeaks.size(); ++i) {
            const int64_t r0 = static_cast<int64_t>(rPeaks[i]);
            // This beat's OWN RR: to the next R, or (last beat only) reused
            // from the previous interval, since there is no "next" for it.
            const int64_t rr = (i + 1 < rPeaks.size())
                ? static_cast<int64_t>(rPeaks[i + 1]) - r0
                : (i > 0 ? r0 - static_cast<int64_t>(rPeaks[i - 1]) : -1);
            if (rr <= 3) continue;

            const int64_t before = static_cast<int64_t>(beforeFrac * rr);
            const int64_t after = static_cast<int64_t>(afterFrac * rr);
            const int64_t len = before + after;
            const int64_t start = r0 - before, end = r0 + after;
            if (len <= 0) continue;

            ProportionalBeat pb;
            pb.samples.assign(static_cast<size_t>(len), std::numeric_limits<double>::quiet_NaN());
            const int64_t cs = std::max<int64_t>(0, start);
            const int64_t ce = std::min<int64_t>(N, end);
            for (int64_t k = cs; k < ce; ++k)
                pb.samples[static_cast<size_t>(k - start)] = ecg[static_cast<size_t>(k)];
            pb.rCol = static_cast<int>(before);
            pb.rrLen = static_cast<int>(rr);
            pb.sgn = sgn;

            // PQ isoelectric zero, in this beat's own local (sliced)
            // coordinates: seed_p_peak / detect_p_end / compute_q_onset all take
            // an r_idx relative to the array they are handed, which pb.rCol
            // already is.
            //
            // seed_p_peak, NOT the landmark. The reported P peak is
            // compute_p_peak, bracketed by the P-onset and Q-onset bars -- but
            // there are no bars here, this is a per-beat slice with no operator
            // marks, and all that is wanted is the rough position that opens
            // detect_p_end's search. The seed is exactly that and nothing else
            // reads it.
            const double qOnD = FeatureMarks::find_q_onset(pb.samples, fs, pb.rCol, sgn);
            const double pPeakD = FeatureMarks::find_p_peak(pb.samples, 0.0, qOnD, fs);
            const int pEnd = FeatureMarks::find_p_end(pb.samples, pb.rCol, fs, 1.0, pPeakD);
            // compute_q_onset's monophasic-R path can return r_idx itself, which
            // would run the PQ window into the R upstroke. Require a real gap.
            const int qGuard = pb.rCol - static_cast<int>(std::lround(0.020 * fs));
            const int qBegin = (qOnD >= 0.0)
                ? static_cast<int>(std::lround(qOnD)) : -1;
            if (pEnd >= 0 && qBegin > pEnd && qBegin <= qGuard) {
                std::vector<double> pq(pb.samples.begin() + pEnd, pb.samples.begin() + qBegin);
                const double base = median_finite(pq);
                if (!std::isnan(base)) {
                    for (double& s : pb.samples) if (!std::isnan(s)) s -= base;
                    pb.pqBaseline = base;
                }
            }
            out.beats.push_back(std::move(pb));
        }

        out.sufficient = out.nExpected > 0
            && static_cast<double>(out.beats.size()) >= 0.5 * out.nExpected;
        return out;
    }

    // ==================================================================
    // Length / area / volume TIME SERIES over segmented beats (5.5)
    // ==================================================================
    // The final assembly the spec's "Build the length, area, and volume
    // time series for each feature" clause asks for: given beats already
    // segmented by segment_beats_proportional (0.25/0.75 RR, PQ-zeroed),
    // compute the three per-segment measures on EACH beat's feature window,
    // producing one value per beat in time order -- i.e. how that feature's
    // morphology evolves across the record, rather than collapsed to a
    // single template.
    //
    // Feature window = the QRS complex (Q-onset -> J-point), auto-detected
    // per beat with the same FeatureMarks detectors segment_beats_
    // proportional already uses for its PQ zero. length and area are
    // PER CHANNEL. volume is inherently 3-lead (segment_volume integrates
    // the vector magnitude), so it needs all three channels segmented from
    // the SAME R-peaks -- pass the three ProportionalBeatSets and it uses
    // the beats at matching indices, sampled at matching R-relative offsets.
    //
    // NaN entries mark beats where the QRS window couldn't be located (or,
    // for volume, where the three beats' windows didn't overlap) -- the
    // series stays index-aligned with the input beats rather than silently
    // shrinking, so a caller can still line each value up with its beat.
    struct FeatureTimeSeries {
        std::vector<double> length;   // per beat, over that beat's QRS window
        std::vector<double> area;     // per beat
        std::vector<double> volume;   // per beat, 3-lead (empty if <3 channels given)
        bool sufficient = false;      // carried through from the segmentation gate
    };

    // Locate a beat's QRS window [q_onset, j_point] in its own local
    // coordinates (rCol is R). Returns {-1,-1} if either landmark is
    // unavailable, which the callers treat as "skip this beat" (NaN).
    // SUB-SAMPLE, not rounded: the detectors return fractional positions and
    // the three integrators take fractional bounds, so the window is carried
    // at the precision it was measured at. -1 is still the absent sentinel.
    inline std::pair<double, double> qrs_window_of(const ProportionalBeat& pb, double fs) {
        if (pb.rCol < 0 || pb.samples.empty() || !(fs > 0.0)) return { -1.0, -1.0 };
        const double qOnset = FeatureMarks::find_q_onset(pb.samples, fs, pb.rCol, pb.sgn);
        const double jPoint = FeatureMarks::find_j_point(pb.samples, fs, pb.rCol, pb.sgn);
        if (std::isnan(qOnset) || std::isnan(jPoint)) return { -1.0, -1.0 };
        if (jPoint <= qOnset) return { -1.0, -1.0 };
        return { qOnset, jPoint };
    }

    // Single-channel: length + area series (volume left empty). Use when
    // only one lead is available or wanted.
    inline FeatureTimeSeries build_feature_time_series(
        const ProportionalBeatSet& beats, double fs)
    {
        FeatureTimeSeries out;
        out.sufficient = beats.sufficient;
        out.length.reserve(beats.beats.size());
        out.area.reserve(beats.beats.size());
        for (const ProportionalBeat& pb : beats.beats) {
            const auto [lo, hi] = qrs_window_of(pb, fs);
            if (lo < 0.0) {
                out.length.push_back(std::nan(""));
                out.area.push_back(std::nan(""));
                continue;
            }
            out.length.push_back(segment_length(pb.samples, lo, hi));
            out.area.push_back(segment_area(pb.samples, lo, hi, /*absolute=*/true));
        }
        return out;
    }

    // Three-channel: length + area (from ch1, the reference lead) AND the
    // 3-lead volume series. The three sets MUST be segmented from the same
    // R-peaks (so beats at index i correspond and share an rCol); volume
    // integrates the vector magnitude over ch1's QRS window, sampled at the
    // same R-relative offset in each channel's beat.
    inline FeatureTimeSeries build_feature_time_series_3ch(
        const ProportionalBeatSet& ch1, const ProportionalBeatSet& ch2,
        const ProportionalBeatSet& ch3, double fs)
    {
        FeatureTimeSeries out;
        out.sufficient = ch1.sufficient;
        const size_t n = std::min({ ch1.beats.size(), ch2.beats.size(), ch3.beats.size() });
        out.length.reserve(n);
        out.area.reserve(n);
        out.volume.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            const ProportionalBeat& b1 = ch1.beats[i];
            const auto [lo, hi] = qrs_window_of(b1, fs);
            if (lo < 0.0) {
                out.length.push_back(std::nan(""));
                out.area.push_back(std::nan(""));
                out.volume.push_back(std::nan(""));
                continue;
            }
            out.length.push_back(segment_length(b1.samples, lo, hi));
            out.area.push_back(segment_area(b1.samples, lo, hi, /*absolute=*/true));

            // Re-register ch2/ch3 onto ch1's R column so the window's sample
            // offsets line up across leads: shift each so its rCol sits at
            // b1.rCol, then integrate the magnitude over [lo, hi]. A beat
            // whose window falls outside a channel's samples contributes NaN
            // there and segment_volume skips it.
            auto shifted = [&](const ProportionalBeatSet& set) {
                std::vector<double> v(b1.samples.size(), std::nan(""));
                if (i >= set.beats.size()) return v;
                const ProportionalBeat& b = set.beats[i];
                const int delta = b1.rCol - b.rCol;   // map b's rCol -> b1's rCol
                for (int k = 0; k < static_cast<int>(b.samples.size()); ++k) {
                    const int dst = k + delta;
                    if (dst >= 0 && dst < static_cast<int>(v.size())) v[dst] = b.samples[k];
                }
                return v;
                };
            const std::vector<double> v2 = shifted(ch2);
            const std::vector<double> v3 = shifted(ch3);
            out.volume.push_back(segment_volume(b1.samples, v2, v3, lo, hi));
        }
        return out;
    }

    // ==================================================================
    // Pulse beat segmentation: foot-anchored, foot-zeroed, ALIGN AND
    // NORMALIZE BEFORE ANY BEAT IS REMOVED
    // ==================================================================
    // The pulse counterpart to segment_beats_proportional, and it is a
    // separate function rather than a parameterization of it for three
    // reasons, none cosmetic:
    //
    //   1. THE ANCHOR IS THE FOOT, not R. Pulse transit time means the wave
    //      for beat n arrives well after that beat's R -- 100-300 ms at the
    //      finger, and it varies with vascular tone within one record. An
    //      R-anchored window would put the systolic peak at a different
    //      column in every beat, which is precisely the smearing the anchored
    //      alignment work exists to remove.
    //   2. THE ZERO IS THE FOOT, not a PQ isoelectric segment. There is no
    //      isoelectric interval in a pulse wave; the foot IS the baseline,
    //      which is why pulse_norm subtracts it per sample.
    //   3. The window runs foot -> foot (one full pulse interval), so it needs
    //      no before/after split of the beat interval.
    //
    // ---- ORDER OF OPERATIONS, WHICH IS THE POINT OF THIS FUNCTION ----
    //
    // Align, then normalize, THEN reject. Not any other order:
    //
    //   * Rejecting before ALIGNING compares samples that are not the same
    //      phase of the wave. Two identical beats offset by 40 ms of transit
    //      time look maximally different sample-for-sample, so a
    //      shape-based outlier test rejects the beats whose transit time
    //      moved, i.e. exactly the physiology being measured.
    //   * Rejecting before NORMALIZING compares raw amplitudes across a
    //      drifting DC baseline. A PPG foot wanders with respiration and with
    //      any change in LED gain, so an amplitude test on un-zeroed traces
    //      rejects on baseline position rather than on pulse size -- and it
    //      does so periodically, at the respiratory rate, which looks like a
    //      real signal in whatever survives.
    //   * Both orderings also bias the SURVIVORS: the reference and the
    //      average are then built from a subset chosen by drift, so the
    //      reference moves with the artifact it was supposed to be immune to.
    //
    // So every beat is sliced, foot-aligned, foot-zeroed and (if a reference
    // is supplied) PI-scaled first. Only then are the exclusion statistics
    // computed, and they are computed on the NORMALIZED samples.
    //
    // Excluded beats are FLAGGED, NOT DROPPED. They stay in .beats, in order,
    // with excluded=true and a reason, so a caller can audit what went and
    // recompute with a different threshold without re-slicing. Anything
    // consuming this for an average must skip excluded beats -- see
    // sqi_weighted_average, which takes a weight per beat and is the intended
    // consumer (weight 0 is the graceful way to express an exclusion).
    struct PulseBeat {
        std::vector<double> samples;   // foot at footCol; foot-zeroed, PI-scaled if ref supplied
        int    footCol = -1;           // column of this beat's foot within samples
        int    peakCol = -1;           // systolic peak, in the same local coordinates
        int    ppLen = -1;             // this beat's own foot-to-foot interval, in samples
        double footBaseline = std::numeric_limits<double>::quiet_NaN();  // the subtracted DC level
        bool   excluded = false;
        const char* exclusionReason = nullptr;   // static string, or nullptr when kept
    };

    struct PulseBeatSet {
        std::vector<PulseBeat> beats;
        int    nExpected = 0;     // record duration / the record's own median foot-to-foot, +1
        int    nKept = 0;         // beats.size() minus the excluded ones
        bool   sufficient = false;  // nKept >= 0.5 * nExpected -- the "at least 50%" gate,
        // evaluated on KEPT beats after normalization, since a
        // gate counting beats that normalization later discards
        // would pass bins that have no usable data.
        double refUsed = std::numeric_limits<double>::quiet_NaN();
    };


    // Length + area series over the pulse wave, per beat.
    //
    // VOLUME IS DELIBERATELY ABSENT, and not because it was skipped. On the ECG
    // side volume is a SPATIAL quantity: ch1/ch2/ch3 are three roughly
    // orthogonal projections of one cardiac dipole, so the triple integral of
    // (x,y,z) has a physical meaning -- that is the same construct
    // compute_ecg_global_ref_spatial and vcg_avg::Loop rest on. PPG, ABP, ART
    // and ART_PULM are four different arteries measured by different
    // transducers, not three axes of one vector, so a product of them is
    // dimensionally a number with no referent. If a cross-channel measure is
    // wanted here it should be named for what it is (e.g. a transit-time or
    // augmentation relationship between two named sites), not called a volume.
    //
    // Window: foot -> end of wave, which for a foot-anchored beat of uniform
    // length is the whole slice. Excluded beats yield NaN so the series stays
    // index-aligned with .beats, matching FeatureTimeSeries' convention.
    struct PulseFeatureTimeSeries {
        std::vector<double> length;   // arc length of the wave
        std::vector<double> area;     // signed area above the foot
        std::vector<double> amplitude;  // peak above foot; the PI-equivalent per beat
        bool sufficient = false;
    };


    // ==================================================================
    // Length / area / volume time-series CSV (5.5), from raw per-bin data
    // ==================================================================
    // Runs the full 0.25/0.75-RR proportional segmentation + PQ-zero +
    // length/area/volume triad on the ACTUAL raw per-channel ECG and its
    // detected R-peaks (output_binfile_data), one row per (bin, channel,
    // beat). This is the concrete, testable output for the spec clause
    // "Build the length, area, and volume time series for each feature".
    //
    // Volume is 3-lead, so it is written on the ch1 row of each beat (the
    // channels are segmented from their own R-peaks and co-registered on R
    // inside build_feature_time_series_3ch); ch2/ch3 rows leave volume
    // blank. `sufficient` (the >=50%-expected-beats gate) is written per
    // (bin, channel) so a reader can drop under-sampled bins.
    //
    // `Bins` is any range of output_binfile_data (e.g. job.peakResults):
    // needs .ecgSignal/.ecgSignal2/.ecgSignal3 and .ch1/.ch2/.ch3.raw.
    template <class Bins>
    inline bool writeFeatureTimeSeriesCsv(const std::string& path,
        const std::string& subjectId, const Bins& bins, double ecgFs,
        const LeadPolarity& pol)
    {
        std::ofstream f(path, std::ios::trunc);
        if (!f) return false;
        f << "subject_id,bin_index,channel,beat_index,rr_len,pq_baseline,"
            "sufficient,qrs_length,qrs_area,qrs_volume\n";
        f.setf(std::ios::fixed);
        f.precision(6);

        auto num = [&](double v) { return std::isnan(v) ? std::string("") : std::to_string(v); };

        int bi = 0;
        for (const auto& b : bins) {
            const std::vector<double>* sig[3] =
            { &b.ecgSignal, &b.ecgSignal2, &b.ecgSignal3 };
            const std::vector<std::size_t>* rp[3] =
            { &b.ch1.raw, &b.ch2.raw, &b.ch3.raw };

            normalize_features::ProportionalBeatSet segs[3];
            for (int c = 0; c < 3; ++c)
                segs[c] = normalize_features::segment_beats_proportional(*sig[c], *rp[c], ecgFs,
                    pol.sign(c));

            // 3-lead series (length/area from ch1 + cross-lead volume), plus
            // per-channel length/area for ch2/ch3 from their own segments.
            const normalize_features::FeatureTimeSeries ts3 =
                normalize_features::build_feature_time_series_3ch(segs[0], segs[1], segs[2], ecgFs);
            normalize_features::FeatureTimeSeries perCh[3];
            for (int c = 0; c < 3; ++c)
                perCh[c] = normalize_features::build_feature_time_series(segs[c], ecgFs);

            for (int c = 0; c < 3; ++c) {
                const auto& seg = segs[c];
                for (size_t k = 0; k < seg.beats.size(); ++k) {
                    const double len = (k < perCh[c].length.size()) ? perCh[c].length[k] : std::nan("");
                    const double area = (k < perCh[c].area.size()) ? perCh[c].area[k] : std::nan("");
                    const double vol = (c == 0 && k < ts3.volume.size()) ? ts3.volume[k] : std::nan("");
                    f << subjectId << ',' << bi << ",CH" << (c + 1) << ',' << k << ','
                        << seg.beats[k].rrLen << ',' << num(seg.beats[k].pqBaseline) << ','
                        << (seg.sufficient ? 1 : 0) << ','
                        << num(len) << ',' << num(area) << ',' << num(vol) << '\n';
                }
            }
            ++bi;
        }
        return static_cast<bool>(f);
    }

}   // namespace normalize_features