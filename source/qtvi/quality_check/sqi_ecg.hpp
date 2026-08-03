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
#include <string>
#include <vector>

#include "config_entry.hpp"
#include "template_generation/template_io.hpp"
#include "template_marking_gui/feature_marks.hpp"   // FeatureMarks

 // P/QRS/ST sample ranges for one beat, in the beat's own sample coordinates
 // (the R-aligned template / kept-beat coordinate system).
struct Segments {
    int pLo = 0, pHi = 0;       // P wave
    int qrsLo = 0, qrsHi = 0;   // QRS complex
    int stLo = 0, stHi = 0;     // ST segment (J point -> T onset)
    int tHi = 0;                // isoelectric-window start = QRS end / J point
    int nextPLo = 0;            // isoelectric-window end = T-wave onset
};

// Builds a Segments from FeatureMarks' auto-detectors, anchored on r_col
// (the same R column the template was built around) at sample rate fs.
//
// The isoelectric window used for the noise metric is [tHi, nextPLo] =
// [QRS end (J point), T-wave onset] -- the classic isoelectric segment
// right after the QRS, before the T wave starts. So tHi is the QRS-end
// landmark and nextPLo is the T-onset landmark, both taken directly from
// FeatureMarks (no estimation needed there).
//
// ASSUMPTION: P onset (pLo) is estimated by reflecting the detected P-end
// around the detected P-peak (symmetric-P-wave assumption) -- FeatureMarks
// has no plain "P onset" detector (it's reactive / GUI-seeded elsewhere).
// Swap this for something more precise if you have a better source for it.
inline Segments buildSegments(const std::vector<double>& ecg, int r_col, double fs) {
    Segments s{};
    const int n = static_cast<int>(ecg.size());
    if (n == 0 || r_col < 0) return s;

    const int pPeak = FeatureMarks::detect_p_peak(ecg, r_col);
    const int pEnd = FeatureMarks::detect_p_end(ecg, r_col);
    const int qBegin = FeatureMarks::detect_q_begin(ecg, r_col);
    const int sEnd = FeatureMarks::detect_s_end(ecg, r_col, fs);      // QRS end / J point
    const int tBegin = FeatureMarks::detect_t_begin(ecg, r_col, fs);    // T onset
    const int tEnd = FeatureMarks::detect_t_end(ecg, r_col, fs);

    std::cerr << "DEBUG n=" << n << " r_col=" << r_col
        << " pPeak=" << pPeak << " pEnd=" << pEnd
        << " qBegin=" << qBegin << " sEnd=" << sEnd
        << " tBegin=" << tBegin << " tEnd=" << tEnd << "\n";

    auto clampIdx = [&](int v) { return std::max(0, std::min(n - 1, v)); };

    s.pHi = clampIdx(pEnd >= 0 ? pEnd : r_col);
    s.pLo = clampIdx((pPeak >= 0 && pEnd >= 0) ? (2 * pPeak - pEnd) : s.pHi);
    s.qrsLo = clampIdx(qBegin >= 0 ? qBegin : r_col);
    s.qrsHi = clampIdx(sEnd >= 0 ? sEnd : r_col);
    s.stLo = s.qrsHi;
    s.stHi = clampIdx(tBegin >= 0 ? tBegin : s.qrsHi);
    // Isoelectric window for the noise metric: J point -> T onset.
    s.tHi = s.qrsHi;
    s.nextPLo = s.stHi;
    // tEnd isn't used for segment bounds (nothing here needs "after the T
    // wave"), but keep the variable named for clarity/future use.
    (void)tEnd;

    // Keep every range non-decreasing even if a detector fell back/failed.
    s.pHi = std::max(s.pHi, s.pLo);
    s.qrsHi = std::max(s.qrsHi, s.qrsLo);
    s.stHi = std::max(s.stHi, s.stLo);
    s.nextPLo = std::max(s.nextPLo, s.tHi);

    std::cerr << "  -> pLo=" << s.pLo << " pHi=" << s.pHi
        << " qrsLo=" << s.qrsLo << " qrsHi=" << s.qrsHi
        << " stLo=" << s.stLo << " stHi=" << s.stHi
        << " tHi=" << s.tHi << " nextPLo=" << s.nextPLo << "\n";
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
inline double pearsonSQI(const std::vector<double>& a, const std::vector<double>& b) {
    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0; int n = 0;
    const size_t len = std::min(a.size(), b.size());
    for (size_t k = 0; k < len; ++k) {
        if (std::isnan(a[k]) || std::isnan(b[k])) continue;
        sa += a[k]; sb += b[k]; saa += a[k] * a[k];
        sbb += b[k] * b[k]; sab += a[k] * b[k]; ++n;
    }
    if (n < 4) return 0.0;   // too few overlapping samples to trust
    const double ma = sa / n, mb = sb / n;
    const double cov = sab / n - ma * mb;
    const double va = saa / n - ma * ma, vb = sbb / n - mb * mb;
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
    int motionFlag) {
    BeatSQI q{};
    q.templateCorr = pearsonSQI(beat, tmpl);

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
        // Two isoelectric reference points: just before P, and just after
        // the QRS (J point) -- compares baseline level before/after the
        // beat's main deflection to catch drift.
        const double preP = beat[std::clamp(seg.pLo, 0, lastIdx)];
        const double postQRS = beat[std::clamp(seg.tHi, 0, lastIdx)];
        q.baseline = std::max(0.0, 1.0 - std::abs(preP - postQRS) / 0.5);    // 0.5 mV budget
    }
    q.noise = std::max(0.0, 1.0 - stddevSQI(beat, seg.tHi, seg.nextPLo) / 0.1);   // 0.1 mV budget

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
        "baseline,noise,motion,composite,handling\n";

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
    static const char* const handlingName[] = { "INCLUDE", "SUBSTITUTE", "EXCLUDE" };

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
                    absBlk.ecgTemplate, seg, motionFlag);
                f << bin << ',' << ch.key << ',' << bi << ','
                    << q.templateCorr << ',' << q.chiSq0 << ',' << q.chiSqAbs << ','
                    << q.chiSq0_P << ',' << q.chiSq0_QRS << ',' << q.chiSq0_ST << ','
                    << q.baseline << ',' << q.noise << ',' << q.motion << ','
                    << q.composite << ',' << handlingName[q.handling] << '\n';
            }
        }
    }
}