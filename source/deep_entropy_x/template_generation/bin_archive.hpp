#pragma once
//
// bin_archive.hpp
//
// Section 5.5.1 -- bin feature archival. Builds and writes one checkpoint
// row per bin (PQRST morphology, QRS-specific parameters, and the 22 bin
// quality parameters) to disk BEFORE any template deformation pass
// (Section 9.6) runs, so bin-level results survive even if deformation
// never runs, or fails partway through, on a later file.

#include "NormalizeFeatures.hpp"          // same folder (template_generation/)
#include "template_io.hpp"                // same folder
#include "template_marking_gui\feature_marks.hpp"   // FeatureMarks auto-detectors
#include "logging\sqi_ecg.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace bin_archive {

    inline constexpr int kNumEcgCh = 3;
    inline const double kNaN = std::numeric_limits<double>::quiet_NaN();

    // ---------------------------------------------------------------------
    // Per-channel morphology block
    // ---------------------------------------------------------------------
    struct ChannelArchive {
        // Amplitudes, in the template's own units (raw mV -- chN_raw.ecgTemplate).
        double p_amp = kNaN, q_amp = kNaN, r_amp = kNaN, s_amp = kNaN, t_amp = kNaN;

        // Intervals (ms), PER CHANNEL, from this channel's own auto-detected
        // onsets/offsets (see header note -- not the cross-lead global set).
        // See the JT note above.
        double pr_ms = kNaN, qrs_ms = kNaN, qt_ms = kNaN, jt_ms = kNaN;

        // Measured at this channel's own auto-detected J point.
        double st_level = kNaN;

        double t_wave_area = kNaN;   // segment_area, |ecg|, [t_begin, t_end]
        double qrs_area = kNaN;      // segment_area, |ecg|, [q_onset, j_point]

        // Per-sample STD (ddof=1) across the beats behind this channel's
        // template -- chN_raw.ecg_template_iqr, which despite the name holds
        // a STD (see NormalizeFeatures.hpp / create_ecg_templates.hpp step 7).
        std::vector<double> per_sample_std;

        // QRS-specific (5.5.1)
        double upstroke_slope_mv_per_s = kNaN;    // Q -> R
        double downstroke_slope_mv_per_s = kNaN;  // R -> S
        double q_to_r_ratio = kNaN;               // |Q| / |R|
        double s_to_r_ratio = kNaN;               // |S| / |R|
        double j_point_amp = kNaN;                // == st_level; named separately per spec wording
    };

    enum SqiField {
        SQI_TEMPLATE_CORR = 0, SQI_CHISQ0, SQI_CHISQ_ABS,
        SQI_CHISQ0_P, SQI_CHISQ0_QRS, SQI_CHISQ0_ST,
        SQI_BASELINE, SQI_NOISE, SQI_MOTION, SQI_COMPOSITE,
        SQI_FRAC_INCLUDE, kNumSqiFields
    };
    inline constexpr const char* kSqiFieldNames[kNumSqiFields] = {
        "template_corr", "chiSq0", "chiSqAbs",
        "chiSq0_P", "chiSq0_QRS", "chiSq0_ST",
        "baseline", "noise", "motion", "composite", "frac_include"
    };

    // ---------------------------------------------------------------------
    // Pulse channels (5.5.1), all four: PPG / ABP / ART / ART_PULM
    // ---------------------------------------------------------------------
    // Same shape as ChannelArchive, different physiology. The fiducials come
    // from FeatureMarks::detect_ppg_fiducials -- the SAME detector the viewer's
    // pulse glyphs and the markings CSV use, for the same reason
    // buildChannelArchive was changed to call detect_template_landmarks
    // instead of the six finders directly: an archive whose landmarks disagree
    // with what an operator sees for the same bin is worse than no archive.
    //
    // ALL AMPLITUDES ARE FOOT-REFERENCED. Every *_amp below is measured above
    // this template's own onset sample, never as an absolute level. A raw PPG
    // level is a function of LED gain and an arterial level of transducer zero,
    // so absolute amplitudes are not comparable across bins, let alone
    // subjects. This is the same rule that makes normalize_features::pulse_norm
    // take a foot argument, and it is why there is no pulse equivalent of the
    // ECG's st_level (an absolute offset from an isoelectric baseline that
    // pulse waves do not have).
    //
    // The arterial channels carry no construction-time fiducials in
    // BinTemplates -- only PPG has ppg_onset_col / ppg_peak_col -- so all four
    // are detected here on the averaged template. That is consistent rather
    // than a compromise: buildChannelArchive detects on the template too.
    struct PulseArchive {
        // Foot-referenced amplitudes, template units.
        double systolic_amp = kNaN;      // peak - foot; the per-bin PI numerator
        double dicrotic_amp = kNaN;      // notch - foot
        double peak2_amp = kNaN;         // diastolic/reflected peak - foot
        double perfusion_index = kNaN;   // 100 * (peak - foot) / |foot|, calculate_perfusion_index

        // Intervals (ms), from this channel's own detected fiducials.
        double foot_to_peak_ms = kNaN;   // systolic upstroke duration
        double foot_to_notch_ms = kNaN;  // ~ejection duration
        double pulse_width_ms = kNaN;    // foot -> end of wave
        double pw80_ms = kNaN;           // width at 80% of systolic amplitude

        // Shape. upstroke_slope is the pulse analogue of the QRS upstroke;
        // downstroke runs peak -> notch (or peak -> end when the notch is
        // absent), which is the interval a vasomotor change shows up in first.
        double upstroke_slope_per_s = kNaN;
        double downstroke_slope_per_s = kNaN;
        double pulse_area = kNaN;        // signed area above the foot, foot -> end
        double systolic_area = kNaN;     // foot -> notch
        double diastolic_area = kNaN;    // notch -> end
        double reflection_index = kNaN;  // amp(p2)/amp(p1), from the detector
        double aging_index = kNaN;       // APG (b-c-d-e)/a, from the detector
        int    dn_tier = 3;              // 1=IEM, 2=Windkessel, 3=absent
        double dn_confidence = 0.0;

        std::vector<double> per_sample_std;   // the channel's *Template_iqr
        uint64_t n_beats = 0;
    };

    // which: 0=PPG, 1=ABP, 2=ART, 3=ART_PULM -- normalize_features::pulseChanName
    // spells the column stems, so the two cannot drift.
    inline PulseArchive buildPulseArchive(const std::vector<double>& wave,
        const std::vector<double>& waveStd,
        uint64_t nBeats, double fs) {
        PulseArchive out;
        out.per_sample_std = waveStd;
        out.n_beats = nBeats;
        if (wave.empty() || !(fs > 0.0)) return out;

        const int W = static_cast<int>(wave.size());
        const FeatureMarks::PpgFiducials fid =
            FeatureMarks::detect_ppg_fiducials(wave, W, fs);

        const int foot = (fid.onset >= 0.0) ? static_cast<int>(std::lround(fid.onset)) : -1;
        if (foot < 0 || foot >= W) return out;
        const double footY = wave[static_cast<size_t>(foot)];

        auto ampAt = [&](double idx) -> double {
            if (idx < 0.0) return kNaN;
            const int k = static_cast<int>(std::lround(idx));
            if (k < 0 || k >= W) return kNaN;
            return wave[static_cast<size_t>(k)] - footY;   // FOOT-REFERENCED, always
            };
        auto msBetween = [&](double a, double b) -> double {
            if (a < 0.0 || b < 0.0) return kNaN;
            return (b - a) * 1000.0 / fs;
            };

        out.systolic_amp = ampAt(fid.peak);
        out.dicrotic_amp = fid.notch_found ? ampAt(fid.dicrotic) : kNaN;
        out.peak2_amp = ampAt(fid.peak2);
        if (fid.peak >= 0.0) {
            const int pk = static_cast<int>(std::lround(fid.peak));
            if (pk >= 0 && pk < W)
                out.perfusion_index =
                normalize_features::calculate_perfusion_index(wave[static_cast<size_t>(pk)], footY);
        }

        out.foot_to_peak_ms = msBetween(fid.onset, fid.peak);
        out.foot_to_notch_ms = fid.notch_found ? msBetween(fid.onset, fid.dicrotic) : kNaN;
        out.pulse_width_ms = msBetween(fid.onset, fid.end);
        out.pw80_ms = (fid.pw80 >= 0.0) ? fid.pw80 * 1000.0 / fs : kNaN;

        if (!std::isnan(out.systolic_amp) && out.foot_to_peak_ms > 0.0)
            out.upstroke_slope_per_s = out.systolic_amp / (out.foot_to_peak_ms / 1000.0);

        // Downstroke: peak -> notch when the notch resolved, else peak -> end.
        // Falling back to `end` keeps the column populated on the many bins
        // where no notch exists (dn_tier 3), and dn_tier records which of the
        // two definitions produced the number.
        {
            const double farIdx = fid.notch_found ? fid.dicrotic : fid.end;
            const double farAmp = ampAt(farIdx);
            const double dtMs = msBetween(fid.peak, farIdx);
            if (!std::isnan(out.systolic_amp) && !std::isnan(farAmp) && dtMs > 0.0)
                out.downstroke_slope_per_s = (farAmp - out.systolic_amp) / (dtMs / 1000.0);
        }

        // Areas, foot-zeroed and SIGNED -- not rectified. A pulse wave is
        // monophasic above its own foot, so |y| would fold diastolic
        // undershoot back in as if it were more pulse. (The ECG side rectifies
        // for the opposite reason: the QRS is biphasic and would cancel.)
        const int endIdx = (fid.end >= 0.0)
            ? std::min(W - 1, static_cast<int>(std::lround(fid.end))) : (W - 1);
        if (endIdx > foot) {
            std::vector<double> zeroed(static_cast<size_t>(endIdx - foot + 1));
            for (int k = foot; k <= endIdx; ++k)
                zeroed[static_cast<size_t>(k - foot)] = wave[static_cast<size_t>(k)] - footY;
            out.pulse_area = normalize_features::segment_area(
                zeroed, 0, static_cast<int>(zeroed.size()) - 1, /*absolute=*/false);
            if (fid.notch_found) {
                const int nk = static_cast<int>(std::lround(fid.dicrotic)) - foot;
                if (nk > 0 && nk < static_cast<int>(zeroed.size())) {
                    out.systolic_area = normalize_features::segment_area(zeroed, 0, nk, false);
                    out.diastolic_area = normalize_features::segment_area(
                        zeroed, nk, static_cast<int>(zeroed.size()) - 1, false);
                }
            }
        }

        out.reflection_index = fid.ri;
        out.aging_index = fid.agi;
        out.dn_tier = fid.dn_tier;
        out.dn_confidence = fid.dn_confidence;
        return out;
    }


    struct BinArchiveRow {
        int binIndex = -1;
        ChannelArchive ch[kNumEcgCh];
        double qrs_area_spatial = kNaN;   // 3-lead vector-magnitude integral over the QRS window (proxy; see header note)

        // All four pulse channels, indexed by normalize_features::pulseChanName
        // order (0=PPG, 1=ABP, 2=ART, 3=ART_PULM). A channel the recording does
        // not have leaves every field NaN and n_beats 0 -- present as columns,
        // empty as data, so the CSV's shape does not depend on which
        // transducers a given subject wore.
        //
        // NO SPATIAL PULSE ANALOGUE of qrs_area_spatial. That column is a
        // vector-magnitude integral over three orthogonal projections of one
        // cardiac dipole; the four pulse channels are four different arteries,
        // not three axes, so there is nothing to take a magnitude of. See the
        // note on build_pulse_feature_time_series.
        PulseArchive pulse[normalize_features::kNumPulseCh];

        double sqi_mean[kNumSqiFields];
        double sqi_std[kNumSqiFields];
        int    n_beats_scored = 0;        // pooled across ch1/ch2/ch3

        BinArchiveRow() {
            for (int k = 0; k < kNumSqiFields; ++k) { sqi_mean[k] = kNaN; sqi_std[k] = kNaN; }
        }
    };

    // ---------------------------------------------------------------------
    // Build: morphology -- all fiducials AUTO-DETECTED on the channel's own
    // averaged template waveform (bt.chN_raw.ecgTemplate), anchored at its
    // own r_col. No pre-placed markers (none exist at this pipeline stage).
    // ---------------------------------------------------------------------
    inline ChannelArchive buildChannelArchive(const template_io::ChannelMethodTemplate& chRaw,
        double fs) {
        ChannelArchive out;
        const std::vector<double>& ecg = chRaw.ecgTemplate;
        const int rPeak = chRaw.r_col;
        if (ecg.empty() || rPeak < 0 || fs <= 0.0) return out;

        // ONE DETECTOR, the same one the viewer and the bank columns use.
        //
        // This block used to call the six finders itself, and differed from the
        // viewer in TWO ways at once, so the archived landmarks did not match
        // what an operator saw for the same bin:
        //   1. it passed chRaw.r_col straight through, with no
        //      subsample_refine::symmetricExtremum refinement -- and the R
        //      anchor is the search origin for all six, so every one moved;
        //   2. it called compute_t_begin and compute_t_end WITHOUT their seed
        //      arguments, so each re-derived the J-point and T-onset internally
        //      by its own slightly different search instead of reusing the ones
        //      already computed here.
        // detect_template_landmarks does the refinement and the chained calls.
        const FeatureMarks::TemplateLandmarks lm =
            FeatureMarks::detect_template_landmarks(ecg, rPeak, fs);
        if (!lm.valid) return out;

        const double r_peak = lm.r_peak;
        const double qOnsetD = lm.q_begin;
        const double j_point = lm.s_end;
        const double tEndD = lm.t_end;
        const double pPeakD = lm.p_peak;
        const double q_peak = lm.q_peak;
        const double pBeginD = lm.p_begin;

        // Q-peak, S-peak and T-peak are not part of the shared six (nothing in
        // the marking path stores them), so they stay here -- but anchored on
        // the REFINED R and the shared T window, not on the raw r_col, or they
        // would reintroduce the same discrepancy one level down.
        const int rAnchor = static_cast<int>(r_peak);
        const double qIdxD = FeatureMarks::compute_q_peak(ecg, rAnchor, fs);
        const double sIdxD = FeatureMarks::compute_s_peak(ecg, rAnchor, fs);
        const double tPeakD = FeatureMarks::compute_t_peak(ecg, j_point, tEndD);

        // Amplitudes.
        out.r_amp = FeatureMarks::sample_at(ecg, r_peak);   // refined R, not r_col
        out.q_amp = FeatureMarks::sample_at(ecg, qIdxD);
        out.s_amp = FeatureMarks::sample_at(ecg, sIdxD);
        out.p_amp = FeatureMarks::sample_at(ecg, pPeakD);
        out.t_amp = FeatureMarks::sample_at(ecg, tPeakD);
        out.st_level = FeatureMarks::sample_at(ecg, j_point);
        out.j_point_amp = out.st_level;

        // Per-channel intervals (ms) from this channel's own landmarks.
        auto ms = [&](double a, double b) -> double {
            if (std::isnan(a) || std::isnan(b) || b <= a) return kNaN;
            return (b - a) / fs * 1000.0;
            };
        out.pr_ms = ms(pBeginD, qOnsetD);
        out.qrs_ms = ms(qOnsetD, j_point);
        out.qt_ms = ms(qOnsetD, tEndD);
        if (!std::isnan(out.qt_ms) && !std::isnan(out.qrs_ms))
            out.jt_ms = out.qt_ms - out.qrs_ms;

        // Areas (rounded fiducials for the integer window bounds).
        auto idx = [](double d) { return std::isnan(d) ? -1 : static_cast<int>(std::lround(d)); };
        const int qOnset = idx(qOnsetD), jPoint = idx(j_point);
        const int tEnd = idx(tEndD);
        out.t_wave_area = normalize_features::segment_area(ecg, j_point, tEnd, /*absolute=*/true);
        out.qrs_area = normalize_features::segment_area(ecg, qOnset, jPoint, /*absolute=*/true);

        // Slopes (mV/s) across the QRS limbs.
        // Measured against the REFINED R, consistent with r_amp and the
        // landmarks above. Mixing rPeakD into some intervals and the raw r_col
        // into others is how a table becomes internally inconsistent.
        if (!std::isnan(qIdxD) && qIdxD < r_peak) {
            const double dt = (r_peak - qIdxD) / fs;
            if (dt > 0.0) out.upstroke_slope_mv_per_s = (out.r_amp - out.q_amp) / dt;
        }
        if (!std::isnan(sIdxD) && sIdxD > r_peak) {
            const double dt = (sIdxD - r_peak) / fs;
            if (dt > 0.0) out.downstroke_slope_mv_per_s = (out.s_amp - out.r_amp) / dt;
        }
        if (!std::isnan(out.q_amp) && !std::isnan(out.r_amp) && out.r_amp != 0.0)
            out.q_to_r_ratio = std::abs(out.q_amp) / std::abs(out.r_amp);
        if (!std::isnan(out.s_amp) && !std::isnan(out.r_amp) && out.r_amp != 0.0)
            out.s_to_r_ratio = std::abs(out.s_amp) / std::abs(out.r_amp);

        // Per-sample STD (chN_raw.ecg_template_iqr holds a STD despite the name).
        out.per_sample_std = chRaw.ecg_template_iqr;
        return out;
    }

    // 3-lead vector-magnitude integral over the QRS window, built directly
    // from the three chN_raw template waveforms aligned on their own r_col.
    // A proxy for the SVD-basis spatial QRS area (see header note): the true
    // VCG loop needs the TemplateBin/vcg_signal_average machinery, which does
    // not exist at this checkpoint. Uses ch1's auto-detected Q-onset/J-point
    // as the window and each channel's r_col to co-register the three traces.
    inline double computeQrsAreaSpatial(const template_io::BinTemplates& bt, double fs) {
        const template_io::ChannelMethodTemplate* chs[3] = { &bt.ch1_raw, &bt.ch2_raw, &bt.ch3_raw };
        for (int c = 0; c < 3; ++c)
            if (chs[c]->ecgTemplate.empty() || chs[c]->r_col < 0) return kNaN;
        if (fs <= 0.0) return kNaN;

        // Window from ch1's auto-detected QRS onset/offset.
        const int r0 = chs[0]->r_col;
        const double qOnsetD = FeatureMarks::compute_q_onset(chs[0]->ecgTemplate, fs, r0);
        const double jPointD = FeatureMarks::compute_j_point(chs[0]->ecgTemplate, fs, r0);
        if (std::isnan(qOnsetD) || std::isnan(jPointD) || jPointD <= qOnsetD) return kNaN;

        // Integrate sqrt(sum ch_c(t)^2) over the window, each channel sampled
        // at the SAME offset-from-its-own-R so the three are co-registered.
        const int steps = static_cast<int>(std::lround(jPointD)) - static_cast<int>(std::lround(qOnsetD));
        if (steps < 1) return kNaN;
        const int off0 = static_cast<int>(std::lround(qOnsetD)) - r0;   // window start, relative to R

        double area = 0.0;
        bool any = false;
        auto magAt = [&](int rel) -> double {
            double sumsq = 0.0;
            for (int c = 0; c < 3; ++c) {
                const auto& w = chs[c]->ecgTemplate;
                const int k = chs[c]->r_col + rel;
                if (k < 0 || k >= static_cast<int>(w.size()) || std::isnan(w[k])) return kNaN;
                sumsq += w[k] * w[k];
            }
            return std::sqrt(sumsq);
            };
        for (int s = 0; s < steps; ++s) {
            const double a = magAt(off0 + s), b = magAt(off0 + s + 1);
            if (std::isnan(a) || std::isnan(b)) continue;
            area += 0.5 * (a + b);
            any = true;
        }
        return any ? area : kNaN;
    }

    // ---------------------------------------------------------------------
    // Build: 22 bin quality parameters, pooled across kept beats on every
    // ECG channel in this bin. Reuses buildSegments/computeEcgSQI
    // (sqi_ecg.hpp) directly -- the per-beat scoring itself is not
    // re-derived here, only reduced to bin level.
    // ---------------------------------------------------------------------
    inline void poolBinQuality(BinArchiveRow& row, const template_io::BinTemplates& bt,
        const template_io::BeatsFile& beats, double ecgFs) {
        struct ChannelSpec {
            const char* key;
            template_io::ChannelMethodTemplate template_io::BinTemplates::* raw;
            template_io::ChannelMethodTemplate template_io::BinTemplates::* absval;
        };
        // Same channel spec sqi_ecg.hpp's own writeEcgSQICsv driver uses.
        const ChannelSpec channels[] = {
            { "CH1", &template_io::BinTemplates::ch1_raw, &template_io::BinTemplates::ch1_absval },
            { "CH2", &template_io::BinTemplates::ch2_raw, &template_io::BinTemplates::ch2_absval },
            { "CH3", &template_io::BinTemplates::ch3_raw, &template_io::BinTemplates::ch3_absval },
        };
        constexpr int motionFlag = -1;   // same "unavailable at this pipeline stage" convention as sqi_ecg.hpp

        if (row.binIndex < 0 || bt.bad_segment) return;

        double accSum[kNumSqiFields] = { 0.0 };
        double accSumSq[kNumSqiFields] = { 0.0 };
        int n = 0;

        for (const ChannelSpec& ch : channels) {
            const auto& rawBlk = bt.*ch.raw;
            const auto& absBlk = bt.*ch.absval;
            if (rawBlk.ecgTemplate.empty() || rawBlk.r_col < 0) continue;

            const auto it = beats.per_channel_beats.find(ch.key);
            if (it == beats.per_channel_beats.end()
                || static_cast<size_t>(row.binIndex) >= it->second.size()) continue;
            const auto& binBeats = it->second[row.binIndex];
            if (binBeats.empty()) continue;

            const Segments seg = buildSegments(rawBlk.ecgTemplate, rawBlk.r_col, ecgFs);
            for (const auto& beat : binBeats) {
                const BeatSQI q = computeEcgSQI(beat, rawBlk.ecgTemplate, absBlk.ecgTemplate,
                    seg, motionFlag, ecgFs);
                const double vals[kNumSqiFields] = {
                    q.templateCorr, q.chiSq0, q.chiSqAbs, q.chiSq0_P, q.chiSq0_QRS, q.chiSq0_ST,
                    q.baseline, q.noise, static_cast<double>(q.motion), q.composite,
                    (q.handling == BeatSQI::INCLUDE) ? 1.0 : 0.0
                };
                for (int k = 0; k < kNumSqiFields; ++k) {
                    accSum[k] += vals[k];
                    accSumSq[k] += vals[k] * vals[k];
                }
                ++n;
            }
        }

        row.n_beats_scored = n;
        for (int k = 0; k < kNumSqiFields; ++k) {
            if (n < 1) { row.sqi_mean[k] = kNaN; row.sqi_std[k] = kNaN; continue; }
            const double m = accSum[k] / n;
            row.sqi_mean[k] = m;
            if (n < 2) { row.sqi_std[k] = kNaN; continue; }
            // Population variance -> ddof=1 via Bessel's correction, same
            // convention as raw_amplitude_iqr / segment_area's std siblings.
            const double var = std::max(0.0, accSumSq[k] / n - m * m)
                * static_cast<double>(n) / static_cast<double>(n - 1);
            row.sqi_std[k] = std::sqrt(var);
        }
    }

    // ---------------------------------------------------------------------
    // One bin's full archive row. `beats` is optional: pass it to populate
    // the 22 quality parameters; omit it to archive morphology/intervals/
    // areas only. All morphology is auto-detected on the bin's own R-pass
    // template waveforms (bt.chN_raw), so no anchor/marker input is needed.
    // ---------------------------------------------------------------------
    inline BinArchiveRow buildBinArchiveRow(const std::vector<template_io::BinTemplates>& bins, int i,
        double fs, const template_io::BeatsFile* beats = nullptr) {
        BinArchiveRow row;
        row.binIndex = i;
        if (i < 0 || static_cast<size_t>(i) >= bins.size()) return row;
        const template_io::BinTemplates& bt = bins[static_cast<size_t>(i)];
        if (bt.bad_segment) return row;

        const template_io::ChannelMethodTemplate* chs[3] = { &bt.ch1_raw, &bt.ch2_raw, &bt.ch3_raw };
        for (int c = 0; c < kNumEcgCh; ++c)
            row.ch[c] = buildChannelArchive(*chs[c], fs);

        row.qrs_area_spatial = computeQrsAreaSpatial(bt, fs);

        // ---- all four pulse channels ---------------------------------
        // Rate: the pulse channels are resampled onto the ECG grid upstream
        // (the templates are the same length as the ECG ones and share the
        // bin's sample rate), so `fs` is correct for all four. If a channel is
        // ever carried at its own native rate, this is the line that has to
        // learn about it -- the ms intervals inside buildPulseArchive all
        // divide by it.
        {
            const std::vector<double>* pw[normalize_features::kNumPulseCh] = {
                &bt.ppgTemplate, &bt.abpTemplate, &bt.artTemplate, &bt.artPulmTemplate };
            const std::vector<double>* ps[normalize_features::kNumPulseCh] = {
                &bt.ppg_template_iqr, &bt.abpTemplate_iqr,
                &bt.artTemplate_iqr, &bt.artPulmTemplate_iqr };
            // Only PPG carries a beat count in BinTemplates; the arterial
            // channels are foot-averaged over the same beats but their count
            // is not stored, so it is left 0 rather than guessed from ppg_n_beats.
            const uint64_t pn[normalize_features::kNumPulseCh] = { bt.ppg_n_beats, 0, 0, 0 };
            for (int q = 0; q < normalize_features::kNumPulseCh; ++q)
                row.pulse[q] = buildPulseArchive(*pw[q], *ps[q], pn[q], fs);
        }

        if (beats) poolBinQuality(row, bt, *beats, fs);

        return row;
    }

    inline std::vector<BinArchiveRow> buildBinArchive(const std::vector<template_io::BinTemplates>& bins,
        double fs, const template_io::BeatsFile* beats = nullptr) {
        std::vector<BinArchiveRow> rows;
        rows.reserve(bins.size());
        for (size_t i = 0; i < bins.size(); ++i)
            rows.push_back(buildBinArchiveRow(bins, static_cast<int>(i), fs, beats));
        return rows;
    }

    // ---------------------------------------------------------------------
    // CSV
    // ---------------------------------------------------------------------

    // Empty cell for NaN, so a missing value is never read as 0 downstream --
    // same convention as vcg_signal_average.hpp's num().
    inline std::string num(double v) {
        if (std::isnan(v)) return std::string();
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.6g", v);
        return std::string(buf);
    }

    // Per-sample std vector archived as one semicolon-joined cell, so the
    // row stays one line per bin -- a long-format table would need a
    // bin/channel/sample key on every row for no benefit here.
    inline std::string stdvec(const std::vector<double>& v) {
        std::string s;
        s.reserve(v.size() * 8);
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) s += ';';
            s += num(v[i]);
        }
        return s;
    }

    inline std::string csvHeader() {
        std::string h = "subject_id,bin_index,anchor,qrs_area_spatial";
        for (int c = 1; c <= 3; ++c) {
            const std::string p = "ch" + std::to_string(c) + "_";
            h += "," + p + "p_amp," + p + "q_amp," + p + "r_amp," + p + "s_amp," + p + "t_amp,"
                + p + "pr_ms," + p + "qrs_ms," + p + "qt_ms," + p + "jt_ms,"
                + p + "st_level," + p + "j_point_amp,"
                + p + "t_wave_area," + p + "qrs_area,"
                + p + "upstroke_slope," + p + "downstroke_slope,"
                + p + "q_to_r_ratio," + p + "s_to_r_ratio,"
                + p + "per_sample_std";
        }
        for (int q = 0; q < normalize_features::kNumPulseCh; ++q) {
            const std::string p = std::string(normalize_features::pulseChanName(q)) + "_";
            h += "," + p + "systolic_amp," + p + "dicrotic_amp," + p + "peak2_amp,"
                + p + "perfusion_index,"
                + p + "foot_to_peak_ms," + p + "foot_to_notch_ms,"
                + p + "pulse_width_ms," + p + "pw80_ms,"
                + p + "upstroke_slope," + p + "downstroke_slope,"
                + p + "pulse_area," + p + "systolic_area," + p + "diastolic_area,"
                + p + "reflection_index," + p + "aging_index,"
                + p + "dn_tier," + p + "dn_confidence,"
                + p + "n_beats," + p + "per_sample_std";
        }
        h += ",n_beats_scored";
        for (int k = 0; k < kNumSqiFields; ++k) h += std::string(",sqi_") + kSqiFieldNames[k] + "_mean";
        for (int k = 0; k < kNumSqiFields; ++k) h += std::string(",sqi_") + kSqiFieldNames[k] + "_std";
        return h;
    }

    /**
     * @brief Write one row per bin to <dir>/<subjectId>_bin_archive_<anchorLabel>.csv.
     * @return false if the file could not be opened.
     */
    inline bool writeBinArchiveCsv(const std::string& dir, const std::string& subjectId,
        const std::string& anchorLabel, const std::vector<BinArchiveRow>& rows) {
        const std::string path = dir + "/" + subjectId + "_bin_archive_" + anchorLabel + ".csv";
        std::ofstream f(path);
        if (!f) return false;
        f << csvHeader() << "\n";
        for (const BinArchiveRow& r : rows) {
            f << subjectId << ',' << r.binIndex << ',' << anchorLabel << ',' << num(r.qrs_area_spatial);
            for (int c = 0; c < kNumEcgCh; ++c) {
                const ChannelArchive& a = r.ch[c];
                f << ',' << num(a.p_amp) << ',' << num(a.q_amp) << ',' << num(a.r_amp)
                    << ',' << num(a.s_amp) << ',' << num(a.t_amp)
                    << ',' << num(a.pr_ms) << ',' << num(a.qrs_ms) << ',' << num(a.qt_ms) << ',' << num(a.jt_ms)
                    << ',' << num(a.st_level) << ',' << num(a.j_point_amp)
                    << ',' << num(a.t_wave_area) << ',' << num(a.qrs_area)
                    << ',' << num(a.upstroke_slope_mv_per_s) << ',' << num(a.downstroke_slope_mv_per_s)
                    << ',' << num(a.q_to_r_ratio) << ',' << num(a.s_to_r_ratio)
                    << ",\"" << stdvec(a.per_sample_std) << "\"";
            }
            for (int q = 0; q < normalize_features::kNumPulseCh; ++q) {
                const PulseArchive& pa = r.pulse[q];
                f << ',' << num(pa.systolic_amp) << ',' << num(pa.dicrotic_amp)
                    << ',' << num(pa.peak2_amp) << ',' << num(pa.perfusion_index)
                    << ',' << num(pa.foot_to_peak_ms) << ',' << num(pa.foot_to_notch_ms)
                    << ',' << num(pa.pulse_width_ms) << ',' << num(pa.pw80_ms)
                    << ',' << num(pa.upstroke_slope_per_s) << ',' << num(pa.downstroke_slope_per_s)
                    << ',' << num(pa.pulse_area) << ',' << num(pa.systolic_area)
                    << ',' << num(pa.diastolic_area)
                    << ',' << num(pa.reflection_index) << ',' << num(pa.aging_index)
                    << ',' << pa.dn_tier << ',' << num(pa.dn_confidence)
                    << ',' << pa.n_beats
                    << ",\"" << stdvec(pa.per_sample_std) << "\"";
            }
            f << ',' << r.n_beats_scored;
            for (int k = 0; k < kNumSqiFields; ++k) f << ',' << num(r.sqi_mean[k]);
            for (int k = 0; k < kNumSqiFields; ++k) f << ',' << num(r.sqi_std[k]);
            f << "\n";
        }
        return static_cast<bool>(f);
    }

    /**
     * @brief One call for the pipeline: build every bin's archive row and
     *        write the checkpoint CSV.
     *
     *        Call this BEFORE any template deformation pass (Section 9.6)
     *        runs, so bin-level results are on disk even if deformation
     *        never runs, or fails partway through.
     *
     * @param bins        job.tmpl.bins (the R-pass template build).
     * @param anchorLabel filename suffix only (e.g. "R") -- NOT a marker
     *                    source; all fiducials are auto-detected here.
     * @param beats       Optional. Pass to populate the 22 quality params;
     *                    omit for morphology / intervals / areas only.
     */
    inline bool writeBinFeatureArchive(const std::string& dir, const std::string& subjectId,
        const std::vector<template_io::BinTemplates>& bins, double fs,
        const std::string& anchorLabel,
        const template_io::BeatsFile* beats = nullptr) {
        const std::vector<BinArchiveRow> rows = buildBinArchive(bins, fs, beats);
        return writeBinArchiveCsv(dir, subjectId, anchorLabel, rows);
    }

}  // namespace bin_archive