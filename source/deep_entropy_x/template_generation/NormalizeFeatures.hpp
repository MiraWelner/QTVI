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

#include "template_marking_gui\template_marking_bin_io.hpp"
#include "template_marking_gui\global_intervals.hpp"
#include "template_marking_gui\vcg_signal_average.hpp"

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

    inline double sample_y(const std::vector<double>& v, int idx) {
        // Look up sample y-value at marker index, NaN if unavailable.
        if (idx < 0 || idx >= static_cast<int>(v.size())) return std::nan("");
        const double y = v[idx];
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

    // ------------------------------------------------------------------
    // ECG per-channel Global_Ref = median over bins of |R_peak| + |S_peak|.
    // ------------------------------------------------------------------
    inline double compute_ecg_global_ref(const std::vector<TemplateBin>& bins, int ch, double sampleRateHz)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            if (b.bad_r_ch[ch]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;

            // Normalization reference is a stable per-subject quantity, so it
            // always reads the R-pass markers, never the current anchor's.
            // Slot 0's marks, through the one accessor -- which selects the
            // lead, so the old per-lead subscripts are gone. The reference is a
            // per-subject quantity measured on the sinus seed, so slot 0 is the
            // right slot as well as the only one this ever read.
            const tbank::BankMarkerSet& rmk =
                b.slotMarks(ch, 0, AnchorType::R_PEAK);
            EcgFeatures f = computeEcgFeatures(ecg,
                rmk.p_peak, rmk.q_begin, b.r_peak_ch[ch],
                rmk.s_end, rmk.t_end, sampleRateHz);
            const double ry = sample_y(ecg, f.r_idx);
            const double sy = sample_y(ecg, f.s_idx);
            if (std::isnan(ry) || std::isnan(sy)) continue;
            vals.push_back(std::abs(ry) + std::abs(sy));
        }
        return median_finite(std::move(vals));
    }

    // ------------------------------------------------------------------
    // Pulse channel accessor (0=PPG, 1=ABP, 2=ART, 3=ART_PULM).
    // ------------------------------------------------------------------
    struct PulseChannel {
        const std::vector<double>* trace;
        int   foot_idx;
        int   peak_idx;
        uint8_t issue;   // 0 = ok, 1 = user-bad, 2 = absent
        // Added for the area reference below, which needs the far bracket of
        // the wave and not just its peak. dicrotic_idx is the systolic/
        // diastolic divide, end_idx the end of the wave; both are -1 on
        // channels or bins where the notch was not found, which the area
        // reference treats as "fall back to end" and then "skip this bin".
        int   dicrotic_idx;
        int   peak2_idx;
        int   end_idx;
    };

    inline PulseChannel pulseChan(const TemplateBin& b, int which) {
        switch (which) {
        case 0: return { &b.ppgTemplate,     b.ppg_onset,    b.ppg_peak,    b.bad_ppg,
                         b.ppg_dicrotic,     b.ppg_peak2,    b.ppg_end };
        case 1: return { &b.abpTemplate,     b.abp_onset,    b.abp_peak,    b.abp_issue,
                         b.abp_dicrotic,     b.abp_peak2,    b.abp_end };
        case 2: return { &b.artTemplate,     b.art_onset,    b.art_peak,    b.art_issue,
                         b.art_dicrotic,     b.art_peak2,    b.art_end };
        default: return { &b.artPulmTemplate, b.art_pulm_onset, b.art_pulm_peak, b.art_pulm_issue,
                         b.art_pulm_dicrotic, b.art_pulm_peak2, b.art_pulm_end };
        }
    }

    // Column-name stem for the four pulse channels, so the archive and the CSV
    // writers agree on spelling without each keeping its own switch.
    inline const char* pulseChanName(int which) {
        switch (which) {
        case 0:  return "ppg";
        case 1:  return "abp";
        case 2:  return "art";
        default: return "art_pulm";
        }
    }
    inline constexpr int kNumPulseCh = 4;

    inline double compute_pulse_global_ref(const std::vector<TemplateBin>& bins, int which)
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


    inline std::vector<double> scale_pulse_spread_by_ref(const std::vector<double>& raw, double foot_y, double ref) {
        //this scales things so the std band in the ppg is visible in the viewer.
        //The spread is computed in raw units, then scaled to normalized units by dividing by the global reference (median PI) and the local foot value.
        std::vector<double> out(raw.size());
        const bool ok = std::isfinite(ref) && ref != 0.0
            && !std::isnan(foot_y) && std::abs(foot_y) >= 1e-12;
        for (size_t i = 0; i < raw.size(); ++i)
            out[i] = (ok && !std::isnan(raw[i]))
            ? (100.0 * raw[i] / std::abs(foot_y)) / std::abs(ref)
            : std::numeric_limits<double>::quiet_NaN();
        return out;
    }
    inline std::vector<double> normalize_ecg_trace(const std::vector<double>& raw, double ref) {
        return scale_array_by_ref(raw, ref);
    }

    // Pulse: local ratio (per-sample, using THIS trace's own foot) then
    // divide by ref. Works for the mean template or any individual beat --
    // never uses a median/global foot value, per the documented algorithm.
    inline std::vector<double> normalize_pulse_trace(const std::vector<double>& raw, int footIdx, double ref) {
        const double footY = sample_y(raw, footIdx);
        std::vector<double> out(raw.size());
        for (size_t i = 0; i < raw.size(); ++i) out[i] = pulse_norm(raw[i], footY, ref);
        return out;
    }

    // ------------------------------------------------------------------
    // Cross-beat spread helpers, computed once at template-build time from
    // the raw aligned beats -- NOT from individual beats retained downstream
    // (that overlay-beat machinery has been removed; these summary
    // statistics are all that's kept).
    //
    // NOTE: despite the name/callers still saying "iqr" (raw_amplitude_iqr,
    // local_ratio_iqr, ecg_template_iqr, *_iqr columns in the CSV/bin
    // export), this now computes per-sample STD (ddof=1), not a true
    // interquartile range -- made consistent with the ECG side's step-7
    // change in create_ecg_templates.hpp, so every channel's "*_iqr" column
    // holds the same statistic. Renaming these identifiers throughout the
    // codebase (TemplateTypes.hpp, BinPlotWidget, TemplateBinIO,
    // template_io, TemplateViewerWindow's CSV header, ...) is a separate,
    // larger follow-up; left as-is here to keep this change to the
    // computation only.
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
    inline double segment_length(const std::vector<double>& v, int lo, int hi) {
        lo = std::max(0, lo);
        hi = std::min(hi, static_cast<int>(v.size()) - 1);
        if (hi <= lo) return std::nan("");
        double len = 0.0;
        bool any = false;
        for (int i = lo; i < hi; ++i) {
            if (std::isnan(v[i]) || std::isnan(v[i + 1])) continue;
            const double dy = v[i + 1] - v[i];
            len += std::sqrt(1.0 + dy * dy);
            any = true;
        }
        return any ? len : std::nan("");
    }

    // Trapezoidal area under v over [lo, hi]. `absolute` rectifies before
    // integrating -- the conventional way to report QRS/T-wave area, since a
    // biphasic complex would otherwise partially cancel itself in a signed
    // integral. NaN samples are skipped (that trapezoid contributes nothing,
    // rather than propagating NaN across the whole sum).
    inline double segment_area(const std::vector<double>& v, int lo, int hi, bool absolute = true) {
        lo = std::max(0, lo);
        hi = std::min(hi, static_cast<int>(v.size()) - 1);
        if (hi <= lo) return std::nan("");
        double area = 0.0;
        bool any = false;
        for (int i = lo; i < hi; ++i) {
            double a = v[i], b = v[i + 1];
            if (std::isnan(a) || std::isnan(b)) continue;
            if (absolute) { a = std::abs(a); b = std::abs(b); }
            area += 0.5 * (a + b);   // trapezoid, unit width
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
        const std::vector<double>& ch3, int lo, int hi) {
        const int n = static_cast<int>(std::min({ ch1.size(), ch2.size(), ch3.size() }));
        lo = std::max(0, lo);
        hi = std::min(hi, n - 1);
        if (hi <= lo) return std::nan("");
        auto mag = [&](int i) -> double {
            const double x = ch1[i], y = ch2[i], z = ch3[i];
            if (std::isnan(x) || std::isnan(y) || std::isnan(z)) return std::nan("");
            return std::sqrt(x * x + y * y + z * z);
            };
        double vol = 0.0;
        bool any = false;
        for (int i = lo; i < hi; ++i) {
            const double a = mag(i), b = mag(i + 1);
            if (std::isnan(a) || std::isnan(b)) continue;
            vol += 0.5 * (a + b);
            any = true;
        }
        return any ? vol : std::nan("");
    }

    // ==================================================================
    // Section 5.2 -- Global reference, Options A/B/C, and the CV check
    // ==================================================================
    //
    // Three ways to reduce one subject's beats to a single scalar
    // Global_Ref_person, all median-across-bins the same way Option A
    // (compute_ecg_global_ref, above) always has:
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

    // Option B: area-based Global_Ref_person for ONE channel. Mirrors
    // compute_ecg_global_ref's per-channel shape and bad_r_ch/bad_segment
    // gating exactly, swapping the |R|+|S| reduction for the QRS's
    // rectified area (Q-onset -> J-point / s_end).
    inline double compute_ecg_global_ref_area(const std::vector<TemplateBin>& bins, int ch, double sampleRateHz)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            if (b.bad_r_ch[ch]) continue;
            const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
            const auto& ecg = chs[ch]->ecgTemplate_raw;
            if (ecg.empty()) continue;

            // Same rule as Option A: the reference is a stable per-subject
            // quantity, so it always reads the R-pass markers.
            const tbank::BankMarkerSet& rmk =
                b.slotMarks(ch, 0, AnchorType::R_PEAK);
            const int qBegin = rmk.q_begin;
            const int jPoint = rmk.s_end;   // S_END == J_POINT (AnchorType comment)
            if (qBegin < 0 || jPoint <= qBegin) continue;
            const double area = segment_area(ecg, qBegin, jPoint, /*absolute=*/true);
            if (!std::isnan(area)) vals.push_back(area);
        }
        return median_finite(std::move(vals));
    }

    // Option C: spatial vector-magnitude Global_Ref_person, fusing all three
    // ECG channels. Builds the R-relative 3-lead loop the SAME way
    // vcg_signal_average.hpp's save-time path does (each channel read at ITS
    // OWN r_col + offset -- see vcg_signal_average.hpp's "AXIS" note), so this
    // does not re-derive cross-channel alignment; it reuses the one already
    // proven for the VCG loop. Global_Ref_person is the median, across bins,
    // of each bin's peak spatial magnitude sqrt(x^2+y^2+z^2).
    inline double compute_ecg_global_ref_spatial(const std::vector<TemplateBin>& bins,
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

    // ------------------------------------------------------------------
    // Option B for the pulse channels: area-based Global_Ref_person
    // ------------------------------------------------------------------
    // The pulse analogue of compute_ecg_global_ref_area, and it differs from it
    // in the one way that matters: the ECG version integrates a RECTIFIED
    // signal, because the QRS is biphasic and a signed integral would let the
    // Q and S deflections cancel the R. A pulse wave is monophasic above its
    // own foot, so the foot is subtracted and the integral is SIGNED -- taking
    // |y| here instead would fold diastolic undershoot back into the area as if
    // it were more pulse.
    //
    // That foot subtraction is also what makes this a pulse reference rather
    // than a scaled ECG one: pulse amplitude is only meaningful relative to the
    // beat's own baseline (see pulse_norm / calculate_perfusion_index), and the
    // absolute DC level of a PPG trace is an arbitrary function of LED gain.
    //
    // Brackets: foot -> end of wave. Falls back to the dicrotic notch when the
    // wave end was not found, and skips the bin when neither is available --
    // integrating to a guessed endpoint would make the reference a function of
    // how far the guess ran.
    inline double compute_pulse_global_ref_area(const std::vector<TemplateBin>& bins, int which)
    {
        std::vector<double> vals;
        vals.reserve(bins.size());
        for (const auto& b : bins) {
            if (b.bad_segment) continue;
            const PulseChannel pc = pulseChan(b, which);
            if (pc.issue != 0) continue;
            if (pc.trace->empty()) continue;
            const int lo = pc.foot_idx;
            int hi = (pc.end_idx > lo) ? pc.end_idx : pc.dicrotic_idx;
            if (lo < 0 || hi <= lo) continue;
            if (hi >= static_cast<int>(pc.trace->size())) hi = static_cast<int>(pc.trace->size()) - 1;
            if (hi <= lo) continue;

            // Foot-zeroed, then integrated. segment_area's absolute=true is
            // deliberately NOT used -- see the note above.
            const double footY = sample_y(*pc.trace, lo);
            if (std::isnan(footY)) continue;
            std::vector<double> zeroed(static_cast<size_t>(hi - lo + 1));
            for (int k = lo; k <= hi; ++k)
                zeroed[static_cast<size_t>(k - lo)] = sample_y(*pc.trace, k) - footY;
            const double a = segment_area(zeroed, 0, static_cast<int>(zeroed.size()) - 1,
                /*absolute=*/false);
            if (!std::isnan(a)) vals.push_back(a);
        }
        return median_finite(std::move(vals));
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

    // THE PULSE FORM IS NOT THE SAME FUNCTION, and calling ratio_norm on a
    // pulse feature is a silent error rather than a missing feature: it would
    // return a plausible number computed without subtracting the foot.
    //
    // ECG is measured from an isoelectric baseline that is already ~0 after PQ
    // leveling, so ratio normalization is a scalar divide. A pulse wave sits on
    // an arbitrary DC offset -- LED gain for PPG, transducer zero for the
    // arterial lines -- so its amplitude only means anything relative to that
    // beat's own foot. This is the same asymmetry that makes pulse_norm take
    // three arguments where ecg_norm takes two, and that makes
    // local_ratio_iqr exist separately from raw_amplitude_iqr.
    //
    // footY is THIS beat's (or this trace's) own foot, not a per-subject
    // constant: the foot drifts with respiration and vasomotion, and using a
    // single subject-wide baseline would push that drift into the feature.
    inline double pulse_ratio_norm(double featurePeak, double footY, double ref) {
        return pulse_norm(featurePeak, footY, ref);
    }

    // ==================================================================
    // Section 5.4 -- Percentile scaling
    // ==================================================================
    // Maps a ratio-normalized value onto 0-100 using the subject's OWN 2nd
    // and 98th percentile ratio values (p2/p98), clamped at both ends so an
    // outlier beyond the calibration range saturates rather than escaping
    // the scale. NaN when the calibration range itself is degenerate
    // (p98 <= p2) or the input ratio is NaN -- there is no meaningful
    // position on a zero-width or undefined scale.
    inline double pct_scale(double ratio, double p2, double p98) {
        const double range = p98 - p2;
        if (std::isnan(ratio) || !(range > 0.0)) return std::nan("");
        return std::clamp((ratio - p2) / range * 100.0, 0.0, 100.0);
    }

    // ==================================================================
    // Section 5.5 -- SQI-weighted signal average
    // ==================================================================
    // Per-sample weighted mean of `beats`, each weighted by its own SQI
    // score, so a handful of low-quality beats can no longer pull the
    // average as hard as a high-quality one. This is the WEIGHTED
    // alternative to the plain per-sample median create_ecg_templates.hpp
    // uses for the displayed template (that one is unweighted by design --
    // a robust order statistic, not a quality-aware mean); the two serve
    // different purposes and neither replaces the other.
    //
    // Hardened relative to the literal 5.5 draft: beats are allowed to be
    // ragged (indexed only up to their own length, not a hardcoded W) and
    // both NaN samples and non-positive/NaN weights are skipped rather than
    // propagating into every column of the average.
    inline std::vector<double> sqi_weighted_average(
        const std::vector<std::vector<double>>& beats, const std::vector<double>& sqi)
    {
        const size_t nBeats = beats.size();
        if (nBeats == 0 || sqi.size() != nBeats) return {};
        size_t W = 0;
        for (const auto& bt : beats) W = std::max(W, bt.size());
        if (W == 0) return {};

        std::vector<double> num(W, 0.0), den(W, 0.0);
        for (size_t t = 0; t < nBeats; ++t) {
            const double w = sqi[t];
            if (!(w > 0.0) || std::isnan(w)) continue;
            const auto& bt = beats[t];
            for (size_t j = 0; j < bt.size(); ++j) {
                if (std::isnan(bt[j])) continue;
                num[j] += w * bt[j];
                den[j] += w;
            }
        }
        std::vector<double> out(W, std::nan(""));
        for (size_t j = 0; j < W; ++j) if (den[j] > 0.0) out[j] = num[j] / den[j];
        return out;
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
        double pqBaseline = std::numeric_limits<double>::quiet_NaN();   // subtracted DC level; NaN if PQ unavailable
    };

    struct ProportionalBeatSet {
        std::vector<ProportionalBeat> beats;
        int  nExpected = 0;     // record duration / the record's own median RR, +1
        bool sufficient = false;   // beats.size() >= 0.5 * nExpected (the "at least 50%" gate)
    };

    inline ProportionalBeatSet segment_beats_proportional(
        const std::vector<double>& ecg, const std::vector<size_t>& rPeaks, double fs,
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
            const double pPeakD = FeatureMarks::seed_p_peak(pb.samples, pb.rCol, fs);
            const int pEnd = FeatureMarks::detect_p_end(pb.samples, pb.rCol, fs, pPeakD);
            const double qOnD = FeatureMarks::compute_q_onset(pb.samples, fs, pb.rCol);
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
    inline std::pair<int, int> qrs_window_of(const ProportionalBeat& pb, double fs) {
        if (pb.rCol < 0 || pb.samples.empty() || !(fs > 0.0)) return { -1, -1 };
        const double qOnsetD = FeatureMarks::compute_q_onset(pb.samples, fs, pb.rCol);
        const double jPointD = FeatureMarks::compute_j_point(pb.samples, fs, pb.rCol);
        if (std::isnan(qOnsetD) || std::isnan(jPointD)) return { -1, -1 };
        const int lo = static_cast<int>(std::lround(qOnsetD));
        const int hi = static_cast<int>(std::lround(jPointD));
        if (hi <= lo) return { -1, -1 };
        return { lo, hi };
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
            if (lo < 0) {
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
            if (lo < 0) {
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

    // feet: detected pulse onsets, ascending, in samples (find_foot_pulseox /
    // SegmentPPG output). peaks may be empty, in which case the systolic peak
    // is taken as the maximum of the foot-zeroed upstroke.
    //
    // ref: the subject's Global_Ref_person for this channel
    // (compute_pulse_global_ref for Option A, compute_pulse_global_ref_area for
    // Option B). Pass NaN to skip the PI scaling and keep foot-zeroed raw
    // units -- the alignment and zeroing still happen, so the ordering
    // guarantee above holds either way.
    //
    // tukeyK: fence width for the post-normalization amplitude/area exclusion,
    // in IQRs. 1.5 is the conventional Tukey fence and matches the pruning the
    // template path uses.
    inline PulseBeatSet segment_pulses_foot_anchored(
        const std::vector<double>& pulse,
        const std::vector<size_t>& feet,
        const std::vector<size_t>& peaks,
        double fs,
        double ref = std::numeric_limits<double>::quiet_NaN(),
        double tukeyK = 1.5)
    {
        PulseBeatSet out;
        out.refUsed = ref;
        const int64_t N = static_cast<int64_t>(pulse.size());
        if (N == 0 || feet.size() < 2 || !(fs > 0.0)) return out;

        // Expected beat count from the record's own median foot-to-foot, the
        // pulse analogue of the median-RR expectation on the ECG side.
        std::vector<double> pp;
        pp.reserve(feet.size() - 1);
        for (size_t i = 1; i < feet.size(); ++i)
            pp.push_back(static_cast<double>(feet[i] - feet[i - 1]));
        const double medPP = median_finite(pp);
        if (!(medPP > 0.0)) return out;
        out.nExpected = static_cast<int>(static_cast<double>(N) / medPP) + 1;

        // Uniform window from the median interval, so every beat lands on the
        // same column grid -- a per-beat window length would re-introduce the
        // misalignment this function exists to remove. A short beat is
        // zero-padded at the tail rather than stretched: resampling to a common
        // length would change the wave's duration, which is a feature here.
        const int W = static_cast<int>(std::lround(medPP));
        if (W < 4) return out;

        // ---- 1. SLICE AND ALIGN (foot at column 0) --------------------
        for (size_t i = 0; i + 1 < feet.size(); ++i) {
            const int64_t f0 = static_cast<int64_t>(feet[i]);
            if (f0 < 0 || f0 >= N) continue;

            PulseBeat pb;
            pb.samples.assign(static_cast<size_t>(W), std::nan(""));
            pb.footCol = 0;
            pb.ppLen = static_cast<int>(feet[i + 1] - feet[i]);
            for (int k = 0; k < W; ++k) {
                const int64_t src = f0 + k;
                if (src >= N) break;
                pb.samples[static_cast<size_t>(k)] = pulse[static_cast<size_t>(src)];
            }

            // ---- 2. NORMALIZE (foot-zero, then PI-scale) --------------
            const double footY = pb.samples[0];
            if (std::isnan(footY)) {
                pb.excluded = true;
                pb.exclusionReason = "no foot sample";
                out.beats.push_back(std::move(pb));
                continue;
            }
            pb.footBaseline = footY;
            for (double& y : pb.samples) if (!std::isnan(y)) y -= footY;
            if (!std::isnan(ref) && ref != 0.0) {
                // pulse_norm's own form, with the foot already removed above:
                // dividing the zeroed wave by the reference is what makes the
                // result comparable across bins and subjects.
                for (double& y : pb.samples) if (!std::isnan(y)) y /= ref;
            }

            // Peak located AFTER zeroing, on the normalized samples, so the
            // search is on the same data the exclusion test will see.
            if (i < peaks.size() && peaks[i] >= feet[i]) {
                const int64_t rel = static_cast<int64_t>(peaks[i]) - f0;
                if (rel >= 0 && rel < W) pb.peakCol = static_cast<int>(rel);
            }
            if (pb.peakCol < 0) {
                double best = -std::numeric_limits<double>::infinity();
                for (int k = 0; k < W; ++k) {
                    const double y = pb.samples[static_cast<size_t>(k)];
                    if (!std::isnan(y) && y > best) { best = y; pb.peakCol = k; }
                }
            }
            out.beats.push_back(std::move(pb));
        }
        if (out.beats.empty()) return out;

        // ---- 3. ONLY NOW REJECT --------------------------------------
        // Two statistics, both on the normalized, aligned samples: pulse
        // amplitude (peak above the foot) and pulse area. Amplitude alone
        // misses a beat with a normal peak and a collapsed or run-on
        // downstroke, which is what a movement artifact usually looks like
        // once the foot has been subtracted.
        std::vector<double> amp, area;
        amp.reserve(out.beats.size());
        area.reserve(out.beats.size());
        for (const PulseBeat& pb : out.beats) {
            if (pb.excluded) continue;
            const double a = (pb.peakCol >= 0)
                ? pb.samples[static_cast<size_t>(pb.peakCol)] : std::nan("");
            amp.push_back(a);
            area.push_back(segment_area(pb.samples, 0, W - 1, /*absolute=*/false));
        }
        auto fences = [tukeyK](std::vector<double> v) {
            std::vector<double> f;
            for (double x : v) if (!std::isnan(x)) f.push_back(x);
            if (f.size() < 4) return std::pair<double, double>{
                -std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity() };
            std::sort(f.begin(), f.end());
            const double q1 = f[f.size() / 4];
            const double q3 = f[(3 * f.size()) / 4];
            const double iqr = q3 - q1;
            return std::pair<double, double>{ q1 - tukeyK * iqr, q3 + tukeyK * iqr };
            };
        const auto [ampLo, ampHi] = fences(amp);
        const auto [arLo, arHi] = fences(area);

        size_t j = 0;
        for (PulseBeat& pb : out.beats) {
            if (pb.excluded) continue;
            const double a = amp[j];
            const double ar = area[j];
            ++j;
            if (std::isnan(a)) { pb.excluded = true; pb.exclusionReason = "no systolic peak"; continue; }
            if (a < ampLo || a > ampHi) { pb.excluded = true; pb.exclusionReason = "amplitude outlier"; continue; }
            if (std::isnan(ar) || ar < arLo || ar > arHi) { pb.excluded = true; pb.exclusionReason = "area outlier"; continue; }
        }

        out.nKept = 0;
        for (const PulseBeat& pb : out.beats) if (!pb.excluded) ++out.nKept;
        out.sufficient = (out.nExpected > 0)
            && (out.nKept >= static_cast<int>(0.5 * out.nExpected));
        return out;
    }

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

    inline PulseFeatureTimeSeries build_pulse_feature_time_series(const PulseBeatSet& set)
    {
        PulseFeatureTimeSeries out;
        out.sufficient = set.sufficient;
        out.length.reserve(set.beats.size());
        out.area.reserve(set.beats.size());
        out.amplitude.reserve(set.beats.size());
        for (const PulseBeat& pb : set.beats) {
            if (pb.excluded || pb.samples.empty()) {
                out.length.push_back(std::nan(""));
                out.area.push_back(std::nan(""));
                out.amplitude.push_back(std::nan(""));
                continue;
            }
            const int hi = static_cast<int>(pb.samples.size()) - 1;
            out.length.push_back(segment_length(pb.samples, 0, hi));
            out.area.push_back(segment_area(pb.samples, 0, hi, /*absolute=*/false));
            out.amplitude.push_back(pb.peakCol >= 0
                ? pb.samples[static_cast<size_t>(pb.peakCol)] : std::nan(""));
        }
        return out;
    }

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
        const std::string& subjectId, const Bins& bins, double ecgFs)
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
                segs[c] = normalize_features::segment_beats_proportional(*sig[c], *rp[c], ecgFs);

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