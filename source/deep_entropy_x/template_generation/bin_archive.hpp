#pragma once
/*
* bin_archive.hpp
* Builds and writes one checkpoint
*  row per bin (PQRST morphology, QRS-specific parameters, and the 22 bin
*  quality parameters) to disk BEFORE any template deformation pass
*  (Section 9.6) runs, so bin-level results survive even if deformation
*  never runs, or fails partway through, on a later file.

* ---------------------------------------------------------------------------
* Two input models, matching how the rest of this codebase already splits
* them
* ---------------------------------------------------------------------------
   - TemplateBin (template_marking_gui/template_marking_bin_io.hpp) -- the
     viewer's marker-aware model. Morphology, intervals, ST level, and the
     area/spatial-volume features are all read from here, the same way
     global_intervals.hpp and vcg_signal_average.hpp already do.
   - template_io::TemplateFile / BeatsFile -- the lower-level per-beat data
     sqi_ecg.hpp's own driver (writeEcgSQICsv) consumes. Quality
     aggregation reuses sqi_ecg.hpp's buildSegments/computeEcgSQI directly
     (no re-derivation of per-beat scoring) and reduces it to bin level
     here -- mean+std across kept beats -- which sqi_ecg.hpp itself does
     not do.

 ASSUMPTION (flagged, not silently decided): `bins[i]` (viewer model) and
 `tmpl.bins[i]` / `beats.per_channel_beats[ch][i]` (template_io model) are
 assumed to be the SAME bin, by position -- the same positional
 correspondence sqi_ecg.hpp's own writeEcgSQICsv loop already relies on.
 If a caller's two inputs were ever built in different orders this breaks
 silently; nothing in the data itself can detect that.

 "22 bin quality parameters": the spec names this count without
 enumerating it. sqi_ecg.hpp's BeatSQI scores 10 numeric components per
 beat (templateCorr, chiSq0, chiSqAbs, chiSq0_P, chiSq0_QRS, chiSq0_ST,
 baseline, noise, motion, composite); adding the fraction of a bin's beats
 classified INCLUDE gives 11 measures, and mean+std of each -- pooled
* across every kept beat on every ECG channel in the bin, since the spec
* names ONE set of 22 per bin rather than one set per channel -- gives
* exactly 22. That pairing is this file's interpretation: swap
* kSqiFieldNames / poolBinQuality below if a different 22 was intended.

* "JT interval": global_intervals.hpp carries PR/QRS/QT but not JT
* directly. ASSUMPTION: JT here means QT minus QRS duration (the ST+T
* segment), the conventional definition, computed once QT and QRS are both
* available.

* "ST level": ASSUMPTION: measured AT this channel's own J point (s_end),
* not at a fixed offset past it (J+60ms / J+80ms is a common alternate
* convention in other protocols) -- flagged so a reviewer can shift it
* without re-deriving anything else here.
*/

#include "NormalizeFeatures.hpp"          // same folder (template_generation/)
#include "template_io.hpp"                // same folder
#include "template_marking_gui\template_marking_bin_io.hpp"
#include "template_marking_gui\global_intervals.hpp"
#include "template_marking_gui\vcg_signal_average.hpp"
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
        // Amplitudes, in whatever units the template carries (raw mV if the
        // caller passes the un-normalized ecgTemplate_raw, as buildBinArchive
        // below does).
        double p_amp = kNaN, q_amp = kNaN, r_amp = kNaN, s_amp = kNaN, t_amp = kNaN;

        // Intervals (ms), GLOBAL (union-of-leads, global_intervals.hpp) --
        // duplicated onto every channel's row for a flat CSV, not re-measured
        // per channel. See the JT ASSUMPTION note above.
        double pr_ms = kNaN, qrs_ms = kNaN, qt_ms = kNaN, jt_ms = kNaN;

        // See the ST-level ASSUMPTION note above.
        double st_level = kNaN;

        double t_wave_area = kNaN;   // segment_area, |ecg|, [t_begin, t_end]
        double qrs_area = kNaN;      // segment_area, |ecg|, [q_begin, s_end]

        // Per-sample STD (ddof=1) of the beats behind this channel's
        // template -- see NormalizeFeatures.hpp's note on the *_iqr naming;
        // archived here under its true statistical name.
        std::vector<double> per_sample_std;

        // QRS-specific (5.5.1)
        double upstroke_slope_mv_per_s = kNaN;    // Q -> R
        double downstroke_slope_mv_per_s = kNaN;  // R -> S
        double q_to_r_ratio = kNaN;               // |Q| / |R|
        double s_to_r_ratio = kNaN;               // |S| / |R|
        double j_point_amp = kNaN;                // == st_level; named separately per spec wording
    };

    // ---------------------------------------------------------------------
    // 22 bin quality parameters
    // ---------------------------------------------------------------------
    // Indices into sqi_mean/sqi_std, matching BeatSQI's fields (sqi_ecg.hpp,
    // global namespace -- that header has no namespace wrapper) plus the
    // derived "fraction INCLUDE".
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

    struct BinArchiveRow {
        int binIndex = -1;
        ChannelArchive ch[kNumEcgCh];
        double qrs_area_spatial = kNaN;   // segment_volume over the 3-lead vector magnitude, QRS window (Option C machinery)

        double sqi_mean[kNumSqiFields];
        double sqi_std[kNumSqiFields];
        int    n_beats_scored = 0;        // pooled across ch1/ch2/ch3

        BinArchiveRow() {
            for (int k = 0; k < kNumSqiFields; ++k) { sqi_mean[k] = kNaN; sqi_std[k] = kNaN; }
        }
    };

    // ---------------------------------------------------------------------
    // Build: morphology
    // ---------------------------------------------------------------------
    inline ChannelArchive buildChannelArchive(const TemplateBin& b, int ch, AnchorType anchor,
        double fs, const global_intervals::GlobalIntervals& g) {
        ChannelArchive out;
        const ChannelTemplateData* chs[3] = { &b.ch1, &b.ch2, &b.ch3 };
        const auto& ecg = chs[ch]->ecgTemplate_raw;
        if (ecg.empty()) return out;

        const TemplateBin::MarkerSet& umk = b.marks(anchor);
        const int pPeak = umk.p_peak_ch[ch];
        const int qBegin = umk.q_begin_ch[ch];
        const int rPeak = b.r_peak_ch[ch];      // auto-only, flat -- TemplateBin's own convention
        const int sEnd = umk.s_end_ch[ch];      // == J point (AnchorType::J_POINT comment)
        const int tBegin = umk.t_begin_ch[ch];
        const int tEnd = umk.t_end_ch[ch];

        out.p_amp = normalize_features::sample_y(ecg, pPeak);
        out.r_amp = normalize_features::sample_y(ecg, rPeak);

        // Q/S peaks inside the QRS: the SAME reactive finders
        // computeEcgFeatures uses for Option A's |R|+|S|, so this amplitude
        // and Option A's reference stay consistent with each other.
        const double qIdxD = FeatureMarks::compute_q_peak(ecg, rPeak, fs);
        const double sIdxD = FeatureMarks::compute_s_peak(ecg, rPeak, fs);
        out.q_amp = FeatureMarks::sample_at(ecg, qIdxD);
        out.s_amp = FeatureMarks::sample_at(ecg, sIdxD);

        // T peak: reactive, the same finder the CSV export / GUI glyph use.
        const double tPeakD = FeatureMarks::compute_t_peak(ecg, tBegin, tEnd);
        out.t_amp = FeatureMarks::sample_at(ecg, tPeakD);

        out.pr_ms = g.prInterval_ms;
        out.qrs_ms = g.qrsDuration_ms;
        out.qt_ms = g.qtInterval_ms;
        if (!std::isnan(out.qt_ms) && !std::isnan(out.qrs_ms))
            out.jt_ms = out.qt_ms - out.qrs_ms;

        out.st_level = normalize_features::sample_y(ecg, sEnd);
        out.j_point_amp = out.st_level;

        if (tBegin >= 0 && tEnd > tBegin)
            out.t_wave_area = normalize_features::segment_area(ecg, tBegin, tEnd, /*absolute=*/true);
        if (qBegin >= 0 && sEnd > qBegin)
            out.qrs_area = normalize_features::segment_area(ecg, qBegin, sEnd, /*absolute=*/true);

        if (!std::isnan(qIdxD) && rPeak >= 0 && rPeak > qIdxD && fs > 0.0) {
            const double dt = (rPeak - qIdxD) / fs;
            if (dt > 0.0) out.upstroke_slope_mv_per_s = (out.r_amp - out.q_amp) / dt;
        }
        if (!std::isnan(sIdxD) && rPeak >= 0 && sIdxD > rPeak && fs > 0.0) {
            const double dt = (sIdxD - rPeak) / fs;
            if (dt > 0.0) out.downstroke_slope_mv_per_s = (out.s_amp - out.r_amp) / dt;
        }
        if (!std::isnan(out.q_amp) && !std::isnan(out.r_amp) && out.r_amp != 0.0)
            out.q_to_r_ratio = std::abs(out.q_amp) / std::abs(out.r_amp);
        if (!std::isnan(out.s_amp) && !std::isnan(out.r_amp) && out.r_amp != 0.0)
            out.s_to_r_ratio = std::abs(out.s_amp) / std::abs(out.r_amp);

        // Already STD, ddof=1 -- see NormalizeFeatures.hpp's note on the
        // *_iqr naming.
        out.per_sample_std = chs[ch]->ecg_template_raw_iqr;
        return out;
    }

    // 3-lead spatial "volume" over the QRS window, reusing the SAME
    // cross-channel R-relative alignment vcg_signal_average.hpp's save-time
    // path already solves (loopFromTemplates), and NormalizeFeatures.hpp's
    // segment_volume for the actual trapezoidal integration.
    inline double computeQrsAreaSpatial(const TemplateBin& b,
        const global_intervals::GlobalIntervals& g, double fs, int marginSamples = 10) {
        if (!g.valid) return kNaN;
        int pre = static_cast<int>(std::ceil(-g.qrsOnset)) + marginSamples;
        int post = static_cast<int>(std::ceil(g.qrsOffset)) + marginSamples;
        if (pre <= 0)  pre = 40 + marginSamples;
        if (post <= 0) post = 60 + marginSamples;

        const vcg_avg::Loop loop = vcg_avg::loopFromTemplates(b, pre, post);
        if (loop.pts.empty()) return kNaN;
        const int qLo = loop.indexForOffset(g.qrsOnset);
        const int qHi = loop.indexForOffset(g.qrsOffset);
        if (qLo < 0 || qHi <= qLo) return kNaN;

        std::vector<double> x(loop.pts.size()), y(loop.pts.size()), z(loop.pts.size());
        for (size_t i = 0; i < loop.pts.size(); ++i) {
            x[i] = loop.pts[i].x; y[i] = loop.pts[i].y; z[i] = loop.pts[i].z;
        }
        return normalize_features::segment_volume(x, y, z, qLo, qHi);
    }

    // ---------------------------------------------------------------------
    // Build: 22 bin quality parameters, pooled across kept beats on every
    // ECG channel in this bin. Reuses buildSegments/computeEcgSQI
    // (sqi_ecg.hpp) directly -- the per-beat scoring itself is not
    // re-derived here, only reduced to bin level.
    // ---------------------------------------------------------------------
    inline void poolBinQuality(BinArchiveRow& row, const template_io::TemplateFile& tmpl,
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

        if (row.binIndex < 0 || static_cast<size_t>(row.binIndex) >= tmpl.bins.size()) return;
        const auto& bt = tmpl.bins[row.binIndex];
        if (bt.bad_segment) return;

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
    // One bin's full archive row. `tmpl`/`beats` are optional: pass both to
    // get the 22 quality parameters populated, or neither to archive
    // morphology/intervals/areas only (e.g. before SQI has been computed for
    // this file yet).
    // ---------------------------------------------------------------------
    inline BinArchiveRow buildBinArchiveRow(const std::vector<TemplateBin>& bins, int i,
        AnchorType anchor, double fs,
        const template_io::TemplateFile* tmpl = nullptr,
        const template_io::BeatsFile* beats = nullptr) {
        BinArchiveRow row;
        row.binIndex = i;
        if (i < 0 || static_cast<size_t>(i) >= bins.size()) return row;
        const TemplateBin& b = bins[static_cast<size_t>(i)];

        const global_intervals::GlobalIntervals g =
            global_intervals::computeGlobalIntervals(b, anchor, fs, global_intervals::MarkerSource::USER);

        for (int c = 0; c < kNumEcgCh; ++c)
            row.ch[c] = buildChannelArchive(b, c, anchor, fs, g);

        row.qrs_area_spatial = computeQrsAreaSpatial(b, g, fs);

        if (tmpl && beats) poolBinQuality(row, *tmpl, *beats, fs);

        return row;
    }

    inline std::vector<BinArchiveRow> buildBinArchive(const std::vector<TemplateBin>& bins,
        AnchorType anchor, double fs,
        const template_io::TemplateFile* tmpl = nullptr,
        const template_io::BeatsFile* beats = nullptr) {
        std::vector<BinArchiveRow> rows;
        rows.reserve(bins.size());
        for (size_t i = 0; i < bins.size(); ++i)
            rows.push_back(buildBinArchiveRow(bins, static_cast<int>(i), anchor, fs, tmpl, beats));
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
     *        never runs, or fails partway through -- the same
     *        "checkpoint before the next stage" placement sqi_ecg.hpp's
     *        writeEcgSQICsv already has relative to finalizeViewerJob's
     *        later steps.
     *
     * @param tmpl/beats  Optional. Pass both to populate the 22 quality
     *                    parameters; pass neither to archive morphology /
     *                    intervals / areas only.
     */
    inline bool writeBinFeatureArchive(const std::string& dir, const std::string& subjectId,
        const std::vector<TemplateBin>& bins, AnchorType anchor, double fs,
        const std::string& anchorLabel,
        const template_io::TemplateFile* tmpl = nullptr,
        const template_io::BeatsFile* beats = nullptr) {
        const std::vector<BinArchiveRow> rows = buildBinArchive(bins, anchor, fs, tmpl, beats);
        return writeBinArchiveCsv(dir, subjectId, anchorLabel, rows);
    }

}  // namespace bin_archive