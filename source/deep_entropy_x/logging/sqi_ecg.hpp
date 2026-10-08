#pragma once
/**
 * @file   sqi_ecg.hpp
 * @brief  Per-beat ECG Signal Quality Index (SQI), scored against the
 *         beat's OWN TEMPLATE -- the bank slot whose members hold it -- as
 *         that template's median and absolute-value median, never against a
 *         bin-wide average.
 *
 *         Wiring: writeEcgSQICsv() is called from analysis_job::finalize()
 *         right after mergeTemplatesSlow() has produced the canonical
 *         job.tmpl / job.beats for a file. It writes one CSV per input file
 *         into cfg.training_log.
 *
 *         Segment boundaries (P/QRS/ST) are derived from FeatureMarks'
 *         existing auto-detectors -- the same ones that seed the viewer's
 *         movable markers -- run on each template's own median, anchored on
 *         the scored beat matrix's r_col.
 *         The one boundary FeatureMarks doesn't expose directly (P onset)
 *         is estimated here; search "ASSUMPTION" below if that needs
 *         tightening.
 *
 *         LEAD POLARITY COMES IN FROM THE CALLER. The Q- and S-side finders
 *         are defined by polarity (Q is the first NEGATIVE deflection before
 *         R), so they need the lead the right way up. That answer is the
 *         operator's per-channel "Lead Reversed" checkbox, carried here as a
 *         LeadPolarity and turned into a sign per channel -- NOT re-derived
 *         from the waveform, which is what the old qrs_positive_at did and
 *         which could disagree with the operator on a complex sitting near
 *         isoelectric.
 *
 *         MOTION AND TIMING COME IN FROM THE CALLER TOO. writeEcgSQICsv takes
 *         the job's beat times and accelerometer result. A beat's time gives
 *         (a) its motion flag -- 0 when any accelerometer sample over the
 *         beat has VM > 1.05 g (Task G) -- and (b) the next beat's P onset,
 *         which closes the TP-segment noise window. Either missing falls
 *         back: motion -1 (no penalty), noise over 80 ms after T-end.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "config_file_handling/config.hpp"
#include "template_generation/template_structs.hpp"
#include "fiducial_marker_finding/feature_marks.hpp"   // FeatureMarks, LeadPolarity
#include "stats_utils.hpp"   // pearson, PearsonResult
#include "prep_for_peakfinding/beat_times.hpp"   // BeatTimes
#include "accel/accel_pipeline.hpp"              // AccelResult (Task G)

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
// sgn is this channel's polarity sign, from LeadPolarity::sign(lead). It is a
// parameter rather than something derived here because this function is handed
// a bare trace with no channel index -- the caller is the only one that knows
// which lead this is.
//
// The isoelectric window used for the noise metric is [tHi, nextPLo] = [this
// beat's T-end, the NEXT beat's P-onset] -- the TP segment, the true flat
// baseline between two consecutive beats. tHi is a direct FeatureMarks
// landmark; nextPLo has no direct detector (it belongs to a beat this
// function was never handed), so it falls back to tHi (a zero-width window,
// handled gracefully by the noise metric).
//
// ASSUMPTION: P onset (pLo) is estimated by reflecting the detected P-end
// around the detected P-peak (symmetric-P-wave assumption). NOTE that
// FeatureMarks::find_p_begin is a real anchor-fit P-onset detector and would
// be a better source than this reflection -- swapping to it is a detector
// change and belongs on its own, but the "no plain P onset detector exists"
// justification this comment used to carry is no longer true.
inline Segments buildSegments(const std::vector<double>& ecg, int r_col, double fs, double sgn) {
    Segments s{};
    const int n = static_cast<int>(ecg.size());
    if (n == 0 || r_col < 0) return s;

    // Each landmark computed ONCE and reused: find_p_end would otherwise
    // re-run the P seed, and the two T detectors would each re-run
    // find_j_point (a full transitionAnchor fit).
    //
    // A P-PEAK SEED, NOT THE P LANDMARK. The reported P peak is find_p_peak
    // bracketed by the P-onset and Q-onset bars -- but there are no bars here:
    // this scores a template with no operator marks, and all that is wanted is
    // the rough position that opens find_p_end's search.
    //
    // find_p_peak and find_t_end take no sgn: both find their extremum by
    // distance from a local bracket baseline, in either direction, so an
    // inverted P or T on an upright lead is found on equal footing. Only the
    // Q/S-side finders need the lead sign.
    const double qOnsetD = FeatureMarks::find_q_onset(ecg, fs, r_col, sgn);
    const double pPeakD = FeatureMarks::find_p_peak(ecg, 0.0, qOnsetD, fs);
    const int pPeak = (int)std::lround(pPeakD);
    const int pEnd = FeatureMarks::find_p_end(ecg, r_col, fs, sgn, pPeakD);
    const int qBegin = (qOnsetD >= 0.0) ? (int)std::lround(qOnsetD) : -1;
    const double jPointD = FeatureMarks::find_j_point(ecg, fs, r_col, sgn);   // QRS end / J point
    const int jPoint = (int)std::lround(jPointD);
    const int tEnd = (int)std::lround(FeatureMarks::find_t_end(ecg, fs, r_col, jPointD));

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
    double chiSqAbs_P = 0.0, chiSqAbs_QRS = 0.0, chiSqAbs_ST = 0.0;   // subsegmental, vs tmplAbs
    double baseline = 0.0, noise = 0.0;
    int motion = 0;               // 1 clean, 0 motion, -1 unavailable
    double composite = 0.0;
    enum Handling { INCLUDE, SUBSTITUTE, EXCLUDE } handling = EXCLUDE;
};

// A WRAPPER over stats_utils::pearson, which is now the one implementation --
// this held a copy of the same five accumulators, and its own comment recorded
// that it was already a copy of alignment.hpp's local lambda.
//
// WHAT THIS KEEPS IS SQI'S POLICY, not the arithmetic: a 4-pair minimum, and
// 0.0 rather than NaN when the correlation is undefined. The zero is
// deliberate here and wrong elsewhere -- BeatSQI feeds a score that is summed
// and thresholded, so an unmeasurable correlation has to contribute a number,
// and "no evidence of a match" is the honest one to contribute. tbank's
// correlate() leaves it NaN for the opposite reason: there, a missing
// measurement must not read as a measured mismatch.
inline double pearsonSQI(const std::vector<double>& a, const std::vector<double>& b,
    int lo = 0, int hi = -1) {
    const PearsonResult pr = pearson(a, b, lo, hi);
    if (pr.n_overlap < 4 || !pr.defined()) return 0.0;
    return pr.r;
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
    const std::vector<double>& tmpl,      // the beat's own template's median, same frame
    const std::vector<double>& tmplAbs,   // that template's median of |beat|
    const Segments& seg,                  // P/QRS/ST sample ranges
    int motionFlag,
    double fs) {
    BeatSQI q{};
    q.templateCorr = std::max(0.0, pearsonSQI(beat, tmpl));   // whole beat, 0..1
    // Against tmplAbs (the median of |beat| over the template) the beat is compared
    // as |beat| too, so a beat identical to the template scores 0 on both.
    // Comparing the signed beat with the all-positive tmplAbs (the spec's code
    // as written) would score every negative deflection as a mismatch.
    auto chi = [&](const std::vector<double>& ref, int a, int b, bool absBeat = false) {
        double s = 0.0;
        const int hi = std::min(b, static_cast<int>(std::min(beat.size(), ref.size())));
        for (int i = std::max(0, a); i < hi; ++i) {
            const double bi = absBeat ? std::abs(beat[i]) : beat[i], ri = ref[i];
            if (std::isnan(bi) || std::isnan(ri)) continue;
            const double d = bi - ri;
            s += d * d;
        }
        return s;
        };

    q.chiSq0 = chi(tmpl, 0, static_cast<int>(beat.size()));
    q.chiSqAbs = chi(tmplAbs, 0, static_cast<int>(beat.size()), true);
    q.chiSq0_P = chi(tmpl, seg.pLo, seg.pHi);
    q.chiSq0_QRS = chi(tmpl, seg.qrsLo, seg.qrsHi);
    q.chiSq0_ST = chi(tmpl, seg.stLo, seg.stHi);
    q.chiSqAbs_P = chi(tmplAbs, seg.pLo, seg.pHi, true);
    q.chiSqAbs_QRS = chi(tmplAbs, seg.qrsLo, seg.qrsHi, true);
    q.chiSqAbs_ST = chi(tmplAbs, seg.stLo, seg.stHi, true);

    // Baseline: the sample at P onset against the sample at T end.
    const int lastIdx = static_cast<int>(beat.size()) - 1;
    if (lastIdx >= 0) {
        const double preP = beat[std::clamp(seg.pLo, 0, lastIdx)];
        const double postT = beat[std::clamp(seg.tHi, 0, lastIdx)];
        q.baseline = (std::isnan(preP) || std::isnan(postT)) ? 0.0
            : std::max(0.0, 1.0 - std::abs(preP - postT) / 0.5); // 0.5 mV
    }
    // Noise over the TP segment, [T end, next beat's P onset]. When the next
    // P onset is unknown (nextPLo not past tHi), 80 ms after T end instead.
    const int noiseHi = (seg.nextPLo > seg.tHi + 1) ? seg.nextPLo
        : seg.tHi + static_cast<int>(std::lround(0.080 * fs));
    q.noise = std::max(0.0, 1.0 - stddevSQI(beat, seg.tHi, noiseHi) / 0.1);   // 0.1 mV budget

    q.motion = motionFlag;
    const double motionTerm = (motionFlag < 0) ? 1.0 : static_cast<double>(motionFlag);
    q.composite = q.templateCorr * q.baseline * q.noise * motionTerm;

    q.handling = q.composite >= 0.70 ? BeatSQI::INCLUDE
        : q.composite >= 0.50 ? BeatSQI::SUBSTITUTE : BeatSQI::EXCLUDE;
    return q;
}

// ---------------------------------------------------------------------
// File-level driver: scores every beat against ITS OWN TEMPLATE -- the bank
// slot (ecg_bank[lead].templates[slot]) whose members hold the beat -- and
// writes one row per beat to <cfg.logs>/<stem>_quality.csv.
//
// PER TEMPLATE, NOT PER BIN. A bin can hold several morphologies (sinus, a
// PVC family, ...); scoring a PVC against the bin-wide median scored it
// against sinus, so every ectopic beat read as poor quality for being
// ectopic. Each template is now its own reference: its median and |x|
// median, its P/QRS/ST segments, and its R-R plausibility bound are all
// taken over that template's rows alone.
//
// THE REFERENCE IS REBUILT FROM THE SCORED MATRIX. per_channel_beats[CHc][bin]
// holds the beats aligned to the scoring anchor (build_bins, forScoring), and
// bank members are local rows of exactly that matrix. The reference is the
// column median of the template's AVERAGED rows (members_clean, else members)
// over that matrix -- so beat and template are co-framed by construction, and
// rawBlk.r_col (the matrix's R column, shared by every row) anchors the
// segments. The bank's own tmpl is not used: it was built before the anchor
// alignment and is in a different frame.
//
// WHICH ROWS ARE SCORED. Every member of a template, including members its
// fences left out of the average (in_template_average = 0), so the quality of
// an excluded beat is still on record. A row that no template claims
// (UNSCORABLE / SPAWN_LIMIT) has no reference and is not written.
//
// Called from analysis_job::finalize() once job.tmpl/job.beats are final.
// `pol` is built there from the job's per-channel inversion flags.
//
// times / accel are optional (nullptr = not available):
//   * motion: a beat spans [t - r_col/fs, t + (len - r_col)/fs] around its R
//     time t; the accelerometer answers 1 clean / 0 motion for that span, -1
//     with no accelerometer or no time.
//   * next P onset: this beat's pLo moved on by the R-R interval to the next
//     kept row, when that interval is plausible (under 1.5x the TEMPLATE's
//     median such interval, so a pruned beat in between does not stretch it);
//     else unknown. The next row may belong to another template -- it is the
//     next beat in time, which is what closes the TP segment.
// ---------------------------------------------------------------------
inline void writeEcgSQICsv(const config_entry& cfg,
    const std::string& stem,
    const template_structs::TemplateFile& tmpl,
    const template_structs::BeatsFile& beats,
    double ecgFs,
    const LeadPolarity& pol,
    const beat_times::BeatTimes* times = nullptr,
    const accel_pipeline::AccelResult* accel = nullptr) {
    const std::string outPath = cfg.logs + "/" + stem + "_quality.csv";
    std::ofstream f(outPath);
    if (!f.is_open()) {
        std::cerr << "  WARNING: could not open " << outPath << " for SQI output\n";
        return;
    }

    f << "bin,channel,template,slot,beat,in_template_average,"
        "template_corr,chiSq0,chiSqAbs,chiSq0_P,chiSq0_QRS,chiSq0_ST,"
        "chiSqAbs_P,chiSqAbs_QRS,chiSqAbs_ST,baseline,noise,motion,composite,is_included\n";

    // `lead` IS IN THE TABLE, not derived from the loop. The loop below is a
    // range-for over this array, so there is no counter to index pol with --
    // and adding one alongside would be a second thing to keep in step with
    // the key. One row, one channel, one polarity index (and one bank index).
    struct ChannelSpec {
        const char* key;
        int lead;   // index into LeadPolarity and ecg_bank; must match key
        template_structs::ChannelMethodTemplate template_structs::BinTemplates::* raw;
    };
    const ChannelSpec channels[] = {
        { "CH1", 0, &template_structs::BinTemplates::ch1_raw },
        { "CH2", 1, &template_structs::BinTemplates::ch2_raw },
        { "CH3", 2, &template_structs::BinTemplates::ch3_raw },
    };
    static const char* const included_levels[] = { "INCLUDE", "SUBSTITUTE", "EXCLUDE" };
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    for (size_t bin = 0; bin < tmpl.bins.size(); ++bin) {
        const auto& bt = tmpl.bins[bin];

        for (const ChannelSpec& ch : channels) {
            // The matrix's R column. Every row of the scored matrix shares it,
            // so it is every template's R column too.
            const int rCol = (bt.*ch.raw).r_col;
            if (rCol < 0) continue;

            const auto it = beats.per_channel_beats.find(ch.key);
            if (it == beats.per_channel_beats.end() || bin >= it->second.size()) continue;
            const auto& binBeats = it->second[bin];   // [row][sample]
            if (binBeats.empty()) continue;

            const tbank::TemplateBank& bank = bt.ecg_bank[ch.lead];
            if (bank.templates.empty()) continue;
            const std::vector<uint8_t> letters = tbank::letterRanks(bank);

            auto tAt = [&](size_t row) {
                return times ? times->at(ch.lead, bin, static_cast<uint32_t>(row)) : kNaN;
                };
            // R-R from a row to the next kept row; NaN where unknown.
            auto rrAfter = [&](size_t row) {
                return (row + 1 < binBeats.size()) ? tAt(row + 1) - tAt(row) : kNaN;
                };

            for (int slot = 0; slot < bank.size(); ++slot) {
                const tbank::template_of_all_signals& tp = bank.templates[slot];
                if (tp.members.empty()) continue;

                // THIS template's references, over its averaged rows.
                const std::vector<uint32_t>& avgRows = tbank::averagedRows(tp);
                const std::vector<double> ref = tbank::rowMedian(binBeats, avgRows, false);
                const std::vector<double> refAbs = tbank::rowMedian(binBeats, avgRows, true);
                if (ref.empty() || rCol >= static_cast<int>(ref.size())) continue;

                const Segments seg = buildSegments(ref, rCol, ecgFs, pol.sign(ch.lead));
                const std::string name = tbank::templateName(bank, slot, letters);

                std::vector<uint32_t> rows = tp.members;
                std::sort(rows.begin(), rows.end());
                rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
                std::vector<char> averaged(binBeats.size(), 0);
                for (const uint32_t r : avgRows) if (r < averaged.size()) averaged[r] = 1;

                // The template's median R-R to the next beat, over its members.
                double medRR = kNaN;
                {
                    std::vector<double> rr;
                    for (const uint32_t r : rows) {
                        const double d = rrAfter(r);
                        if (std::isfinite(d) && d > 0.0) rr.push_back(d);
                    }
                    if (!rr.empty()) {
                        std::nth_element(rr.begin(), rr.begin() + rr.size() / 2, rr.end());
                        medRR = rr[rr.size() / 2];
                    }
                }

                for (const uint32_t bi : rows) {
                    if (bi >= binBeats.size()) continue;   // stale index: drop, never clamp
                    const std::vector<double>& beat = binBeats[bi];
                    if (beat.empty()) continue;
                    const double t = tAt(bi);

                    Segments s = seg;
                    const double rr = rrAfter(bi);
                    if (std::isfinite(rr) && rr > 0.0 && std::isfinite(medRR) && rr < 1.5 * medRR) {
                        const int next = s.pLo + static_cast<int>(std::lround(rr * ecgFs));
                        s.nextPLo = std::min(next, static_cast<int>(beat.size()));
                    }

                    int motionFlag = -1;
                    if (accel && std::isfinite(t) && ecgFs > 0.0)
                        motionFlag = accel->sqiMotionFlag(t - rCol / ecgFs,
                            t + (static_cast<double>(beat.size()) - rCol) / ecgFs);

                    const BeatSQI q = computeEcgSQI(beat, ref, refAbs, s, motionFlag, ecgFs);
                    f << bin << ',' << ch.key << ',' << name << ',' << slot << ',' << bi << ','
                        << int(averaged[bi]) << ','
                        << q.templateCorr << ',' << q.chiSq0 << ',' << q.chiSqAbs << ','
                        << q.chiSq0_P << ',' << q.chiSq0_QRS << ',' << q.chiSq0_ST << ','
                        << q.chiSqAbs_P << ',' << q.chiSqAbs_QRS << ',' << q.chiSqAbs_ST << ','
                        << q.baseline << ',' << q.noise << ',' << q.motion << ','
                        << q.composite << ',' << included_levels[q.handling] << '\n';
                }
            }
        }
    }
}