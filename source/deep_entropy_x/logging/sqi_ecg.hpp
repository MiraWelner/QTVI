#pragma once
/**
 * @file   sqi_ecg.hpp
 * @brief  Per-beat ECG Signal Quality Index (SQI), scored against the
 *         bin/channel's own median template (and absolute-value template)
 *         that the template-generation pipeline already builds.
 *
 *         Wiring: writeEcgSQICsv() is called from finalizeViewerJob()
 *         (post_process.hpp) right after mergeTemplatesSlow() has produced
 *         the canonical job.tmpl / job.beats for a file. It writes one CSV
 *         per input file into cfg.quality_metric.
 *
 *         Segment boundaries (P/QRS/ST) are derived from FeatureMarks'
 *         existing auto-detectors -- the same ones that seed the viewer's
 *         movable markers -- anchored on the bin/channel's own r_col.
 *         The one boundary FeatureMarks doesn't expose directly (P onset)
 *         is estimated here; search "ASSUMPTION" below if that needs
 *         tightening.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "config_file_handling/config_entry.hpp"
#include "template_generation/template_io.hpp"
#include "template_marking_gui/feature_marks.hpp"   // FeatureMarks

 // P/QRS/ST sample ranges for one beat, in the beat's own sample coordinates
 // (the R-aligned template / kept-beat coordinate system).
struct Segments {
    int pLo = 0, pHi = 0;       // P wave
    int qrsLo = 0, qrsHi = 0;   // QRS complex
    int stLo = 0, stHi = 0;     // ST segment (J point -> T onset)
    int tHi = 0;                // isoelectric-window start = this beat's T-end
    int nextPLo = 0;            // isoelectric-window end = next beat's P-onset (TP segment)
};

// Builds a Segments from FeatureMarks' auto-detectors, anchored on r_col
// (the same R column the template was built around) at sample rate fs.
//
// The isoelectric window used for the noise metric is [tHi, nextPLo] = [this
// beat's T-end, the NEXT beat's P-onset] -- the TP segment, the true flat
// baseline between two consecutive beats. tHi is a direct FeatureMarks
// landmark; nextPLo has no direct detector (it belongs to a beat this
// function was never handed), so it's found the same way a normal p_begin
// is -- detect a peak, then FeatureMarks::compute_p_begin's anchor-fit -- just
// with the peak search restricted to [t_end, t_end+600ms]. If the array
// doesn't extend that far (no next beat in view), nextPLo falls back to tHi
// (a zero-width window, handled gracefully by the noise metric).
//
// ASSUMPTION: P onset (pLo) is estimated by reflecting the detected P-end
// around the detected P-peak (symmetric-P-wave assumption) -- FeatureMarks
// has no plain "P onset" detector (it's reactive / GUI-seeded elsewhere).
// Swap this for something more precise if you have a better source for it.
inline Segments buildSegments(const std::vector<double>& ecg, int r_col, double fs) {
    Segments s{};
    const int n = static_cast<int>(ecg.size());
    if (n == 0 || r_col < 0) return s;

    // Each landmark computed ONCE and reused: detect_p_end would otherwise
    // re-run the P seed, and the two T detectors would each re-run
    // compute_j_point (a full transitionAnchor fit).
    //
    // seed_p_peak, NOT the P landmark. The reported P peak is compute_p_peak,
    // bracketed by the P-onset and Q-onset bars -- but there are no bars here:
    // this scores a template with no operator marks, and all that is wanted is
    // the rough position that opens detect_p_end's search. Renamed from
    // detect_p_peak so a call site cannot mistake a seed for a measurement.
    const double pPeakD = FeatureMarks::seed_p_peak(ecg, r_col, fs);
    const int pPeak = (int)std::lround(pPeakD);
    const int pEnd = FeatureMarks::detect_p_end(ecg, r_col, fs, pPeakD);
    // compute_q_onset returns a sub-sample double and has an out-param this
    // caller does not need; -1 means no Q-onset and the fallback below applies.
    const double qBeginD = FeatureMarks::compute_q_onset(ecg, fs, r_col);
    const int qBegin = (qBeginD >= 0.0) ? (int)std::lround(qBeginD) : -1;
    const double jPointD = FeatureMarks::compute_j_point(ecg, fs, r_col);   // QRS end / J point
    const int jPoint = (int)std::lround(jPointD);
    const int tEnd = (int)std::lround(FeatureMarks::compute_t_end(ecg, fs, r_col, jPointD));

    auto clampIdx = [&](int v) { return std::max(0, std::min(n - 1, v)); };

    s.pHi = clampIdx(pEnd >= 0 ? pEnd : r_col);
    s.pLo = clampIdx((pPeak >= 0 && pEnd >= 0) ? (2 * pPeak - pEnd) : s.pHi);
    s.qrsLo = clampIdx(qBegin >= 0 ? qBegin : r_col);
    s.qrsHi = clampIdx(jPoint >= 0 ? jPoint : r_col);
    s.stLo = s.qrsHi;
    s.stHi = tEnd;
    // Isoelectric window for the noise metric: THIS beat's T-end -> the
    // NEXT beat's P-onset -- the TP segment, the true baseline between
    // consecutive beats.
    s.tHi = clampIdx(tEnd >= 0 ? tEnd : s.stHi);

    // Keep every range non-decreasing even if a detector fell back/failed.
    s.pHi = std::max(s.pHi, s.pLo);
    s.qrsHi = std::max(s.qrsHi, s.qrsLo);
    s.stHi = std::max(s.stHi, s.stLo);
    s.nextPLo = std::max(s.nextPLo, s.tHi);
    return s;
}

struct BeatSQI {
    double templateCorr = 0.0, chiSq0 = 0.0, chiSqAbs = 0.0;    // whole-beat
    double chiSq0_P = 0.0, chiSq0_QRS = 0.0, chiSq0_ST = 0.0;   // subsegmental
    double baseline = 0.0, noise = 0.0;
    int motion = 0;               // 1 clean, 0 motion, -1 unavailable
    double composite = 0.0;
    enum Handling { INCLUDE, SUBSTITUTE, EXCLUDE } handling = EXCLUDE;
};

// NaN-aware Pearson correlation (same convention as alignment.hpp's local
// pearson() lambda, factored out here so sqi_ecg.hpp has no dependency on it).
// pearsonSQI signature change: add lo/hi bounds (default = whole array,
// so any other caller is unaffected)
inline double pearsonSQI(const std::vector<double>& a, const std::vector<double>& b,
    int lo = 0, int hi = -1) {
    const int n = static_cast<int>(std::min(a.size(), b.size()));
    if (hi < 0 || hi > n) hi = n;
    lo = std::max(0, lo);
    if (hi - lo < 4) return 0.0;

    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
    int cnt = 0;
    for (int k = lo; k < hi; ++k) {
        const double av = a[k], bv = b[k];
        if (std::isnan(av) || std::isnan(bv)) continue;
        sa += av; sb += bv; saa += av * av; sbb += bv * bv; sab += av * bv;
        ++cnt;
    }
    if (cnt < 4) return 0.0;
    const double ma = sa / cnt, mb = sb / cnt;
    const double cov = sab / cnt - ma * mb;
    const double va = saa / cnt - ma * ma;
    const double vb = sbb / cnt - mb * mb;
    if (va <= 0.0 || vb <= 0.0) return 0.0;
    return cov / std::sqrt(va * vb);
}

inline double stddevSQI(const std::vector<double>& v, int lo, int hi) {
    lo = std::max(0, lo);
    hi = std::min(static_cast<int>(v.size()), hi);
    if (hi - lo < 2) return 0.0;
    double sum = 0.0; int n = 0;
    for (int i = lo; i < hi; ++i) if (!std::isnan(v[i])) { sum += v[i]; ++n; }
    if (n < 2) return 0.0;
    const double mean = sum / n;
    double sq = 0.0;
    for (int i = lo; i < hi; ++i) if (!std::isnan(v[i])) sq += (v[i] - mean) * (v[i] - mean);
    return std::sqrt(sq / (n - 1));
}

inline BeatSQI computeEcgSQI(const std::vector<double>& beat,
    const std::vector<double>& tmpl,      // median template, same length as beat
    const std::vector<double>& tmplAbs,   // absolute-value template
    const Segments& seg,                  // P/QRS/ST sample ranges
    int motionFlag,
    double fs) {
    BeatSQI q{};
    q.templateCorr = pearsonSQI(beat, tmpl, 0, seg.tHi + 0.050 * fs);
    auto chi = [&](const std::vector<double>& ref, int a, int b) {
        double s = 0.0;
        const int hi = std::min(b, static_cast<int>(std::min(beat.size(), ref.size())));
        for (int i = std::max(0, a); i < hi; ++i) {
            const double bi = beat[i], ri = ref[i];
            if (std::isnan(bi) || std::isnan(ri)) continue;
            const double d = bi - ri;
            s += d * d;
        }
        return s;
        };

    q.chiSq0 = chi(tmpl, 0, static_cast<int>(beat.size()));
    q.chiSqAbs = chi(tmplAbs, 0, static_cast<int>(beat.size()));
    q.chiSq0_P = chi(tmpl, seg.pLo, seg.pHi);
    q.chiSq0_QRS = chi(tmpl, seg.qrsLo, seg.qrsHi);
    q.chiSq0_ST = chi(tmpl, seg.stLo, seg.stHi);

    const int lastIdx = static_cast<int>(beat.size()) - 1;
    if (lastIdx >= 0) {
        const int pLo = std::clamp(seg.pLo, 0, lastIdx);
        const int pw = std::max(3, seg.pHi - seg.pLo);        // ~P-wave width
        const int wLo = std::max(0, pLo - pw);
        double preP;
        {
            std::vector<double> win;
            win.reserve(pLo - wLo);
            for (int i = wLo; i < pLo; ++i)
                if (!std::isnan(beat[i])) win.push_back(beat[i]);
            if (win.empty()) {
                preP = beat[pLo];                              // fallback: single sample
            }
            else {
                std::sort(win.begin(), win.end());
                const size_t m = win.size() / 2;
                preP = (win.size() % 2 == 0)
                    ? 0.5 * (win[m - 1] + win[m]) : win[m];
            }
        }

        const double postT = beat[std::clamp(seg.tHi, 0, lastIdx)];
        q.baseline = std::max(0.0, 1.0 - std::abs(preP - postT) / 0.5); // 0.5 mV
    }
    const int noiseHi = seg.tHi + static_cast<int>(std::lround(0.080 * fs));
    q.noise = std::max(0.0, 1.0 - stddevSQI(beat, seg.tHi, noiseHi) / 0.1);   // 0.1 mV budget

    q.motion = motionFlag;
    const double motionTerm = (motionFlag < 0) ? 1.0 : static_cast<double>(motionFlag);
    q.composite = q.templateCorr * q.baseline * q.noise * motionTerm;

    q.handling = q.composite >= 0.70 ? BeatSQI::INCLUDE
        : q.composite >= 0.50 ? BeatSQI::SUBSTITUTE : BeatSQI::EXCLUDE;
    return q;
}

// ---------------------------------------------------------------------
// File-level driver: scores every kept beat, on every ECG channel, in
// every bin, against that bin/channel's own raw + absval templates, and
// writes one row per beat to <cfg.quality_metric>/<stem>_quality.csv.
//
// Called from finalizeViewerJob() once job.tmpl/job.beats are final.
// motionFlag is left at -1 (unavailable): at this stage of the pipeline
// there's no finer per-beat motion signal available than bad_segment
// (already filtered out above), so composite falls back to the
// templateCorr*baseline*noise product with no motion penalty.
// ---------------------------------------------------------------------
inline void writeEcgSQICsv(const config_entry& cfg,
    const std::string& stem,
    const template_io::TemplateFile& tmpl,
    const template_io::BeatsFile& beats,
    double ecgFs) {
    const std::string outPath = cfg.quality_metric + "/" + stem + "_quality.csv";
    std::ofstream f(outPath);
    if (!f.is_open()) {
        std::cerr << "  WARNING: could not open " << outPath << " for SQI output\n";
        return;
    }

    f << "bin,channel,beat,template_corr,chiSq0,chiSqAbs,chiSq0_P,chiSq0_QRS,chiSq0_ST,"
        "baseline,noise,motion,composite,is_included\n";

    struct ChannelSpec {
        const char* key;
        template_io::ChannelMethodTemplate template_io::BinTemplates::* raw;
        template_io::ChannelMethodTemplate template_io::BinTemplates::* absval;
    };
    const ChannelSpec channels[] = {
        { "CH1", &template_io::BinTemplates::ch1_raw, &template_io::BinTemplates::ch1_absval },
        { "CH2", &template_io::BinTemplates::ch2_raw, &template_io::BinTemplates::ch2_absval },
        { "CH3", &template_io::BinTemplates::ch3_raw, &template_io::BinTemplates::ch3_absval },
    };
    static const char* const included_levels[] = { "INCLUDE", "SUBSTITUTE", "EXCLUDE" };

    constexpr int motionFlag = -1;   // see comment above

    for (size_t bin = 0; bin < tmpl.bins.size(); ++bin) {
        const auto& bt = tmpl.bins[bin];
        if (bt.bad_segment) continue;

        for (const ChannelSpec& ch : channels) {
            const auto& rawBlk = bt.*ch.raw;
            const auto& absBlk = bt.*ch.absval;
            if (rawBlk.ecgTemplate.empty() || rawBlk.r_col < 0) continue;

            const auto it = beats.per_channel_beats.find(ch.key);
            if (it == beats.per_channel_beats.end() || bin >= it->second.size()) continue;
            const auto& binBeats = it->second[bin];   // [beat][sample]
            if (binBeats.empty()) continue;

            const Segments seg = buildSegments(rawBlk.ecgTemplate, rawBlk.r_col, ecgFs);

            for (size_t bi = 0; bi < binBeats.size(); ++bi) {
                const BeatSQI q = computeEcgSQI(binBeats[bi], rawBlk.ecgTemplate,
                    absBlk.ecgTemplate, seg, motionFlag, ecgFs);
                f << bin << ',' << ch.key << ',' << bi << ','
                    << q.templateCorr << ',' << q.chiSq0 << ',' << q.chiSqAbs << ','
                    << q.chiSq0_P << ',' << q.chiSq0_QRS << ',' << q.chiSq0_ST << ','
                    << q.baseline << ',' << q.noise << ',' << q.motion << ','
                    << q.composite << ',' << included_levels[q.handling] << '\n';
            }
        }
    }
}