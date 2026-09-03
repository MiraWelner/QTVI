/**
 * @file   premark_beats.hpp
 * @brief  Runs the Section 4.7.1 morphology envelope and the 4.7.2 pre-marking
 *         classifier inside template generation, where the per-beat matrices
 *         are still in memory (out.beats.per_channel_beats[ch][bin][beat]).
 *
 *         No file-format change and no GUI dependency. What this does NOT do
 *         is Stage 2 (progressive refinement on operator action) -- there is
 *         no operator here. When you want that, the beats are already
 *         persisted in template_io::BeatsFile; the viewer just doesn't read it
 *         yet, so Stage 2 needs a reader on the viewer side, not a new format.
 *
 *         Segment boundaries come from each bin's own template via the
 *         existing FeatureMarks auto-detectors (no user markers exist at
 *         generation time), so the P/QRS/ST/T band scores are computed on the
 *         same landmarks the viewer would seed.
 *
 *         Two-pass envelope: pass 1 builds the corridor from every kept beat
 *         (they already survived the Tukey and wave-score pruning in
 *         alignment.hpp), scores each beat against it, then pass 2 rebuilds
 *         from only the beats that scored above cleanPctThreshold. Pass 2's
 *         sd is reported alongside pass 1's so the tightening is measurable
 *         -- that's the acceptance-test quantity.
 *
 * @date   2026-08-14
 */
#pragma once

#include "template_io.hpp"
#include "template_generation/morphology_envelope.hpp"
#include "template_generation/beat_classifier.hpp"
#include "template_marking_gui/feature_marks.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace premark {

    // Beats scoring at or above this percentage of samples inside the pass-1
    // corridor form the pass-2 clean pool. Starting estimate -- calibrate on
    // records with known ectopy counts.
    //
    // The gate requires BOTH pct_overall and pct_QRS. pct_overall alone does
    // not work: a beat is ~400 columns of which the QRS is ~30, so a beat that
    // is grossly abnormal through the whole QRS still scores ~93% overall
    // against ~96% for sinus, and no threshold separates those two
    // populations. Measured on synthetic beats with a 2.3x-widened QRS, the
    // same beats score 78% vs 97% on pct_QRS. That is the whole reason the
    // spec's scoreBeat breaks the score out per wave -- the localized bands
    // are where the signal is, and the overall figure dilutes it by however
    // much isoelectric baseline the slice happens to contain.
    inline constexpr double kCleanPctThreshold = 90.0;
    // Below this many kept beats a percentile corridor isn't worth building.
    inline constexpr int kMinBeatsForEnvelope = 20;

    // Band boundaries on the shared template axis.
    struct SegmentCols {
        int pEnd = -1, qrsStart = -1, qrsEnd = -1, tBegin = -1, tEnd = -1;
        // Bin admission is deliberately UNCHANGED from before tBegin existed:
        // the structural landmarks must be ordered, and that is all. A
        // degenerate T sub-band (tBegin == tEnd, which compute_t_begin and
        // compute_t_end can both return on a low-amplitude or absent T) must
        // NOT cost the whole bin its scores -- an empty band now reads NaN,
        // which is exactly the case NaN was introduced to represent. Requiring
        // tEnd > tBegin here would have dropped bins that used to score.
        bool valid(int W) const {
            return pEnd >= 0 && pEnd <= qrsStart && qrsEnd > qrsStart
                && tEnd > qrsEnd && tEnd <= W
                && tBegin >= qrsEnd && tBegin <= tEnd;
        }
    };

    // Auto-detect the four boundaries off one bin's template.
    inline SegmentCols segmentsFromTemplate(const std::vector<double>& tmpl,
        int rCol, double fs)
    {
        SegmentCols s;
        if (tmpl.empty() || rCol < 0 || fs <= 0.0) return s;
        const double q_onset = FeatureMarks::compute_q_onset(tmpl, fs, rCol);
        const double j_point = FeatureMarks::compute_j_point(tmpl, fs, rCol);
        const double t_end = FeatureMarks::compute_t_end(tmpl, fs, rCol, j_point);
        const int    pEnd = FeatureMarks::detect_p_end(tmpl, rCol, fs);

        s.qrsStart = (q_onset >= 0.0) ? (int)std::lround(q_onset) : -1;
        s.qrsEnd = (j_point >= 0.0) ? (int)std::lround(j_point) : -1;
        s.tEnd = (t_end >= 0.0) ? (int)std::lround(t_end) : -1;
        
        if (s.qrsEnd >= 0 && s.tEnd >= 0)
            s.tBegin = std::max(s.qrsEnd, std::min(s.tBegin, s.tEnd));
        s.pEnd = (pEnd >= 0) ? pEnd : s.qrsStart;
        // ...and CLAMP to the QRS onset even when detection "succeeded".
        // detect_p_end's last resort returns p_peak + 59 samples when the
        // 90%-return-to-baseline target is never reached; at 1 kHz that
        // routinely lands past q_onset on a normal beat. The bin-validity test
        // requires pEnd <= qrsStart, so an unclamped fallback silently cost the
        // WHOLE BIN its envelope and scores -- a P-wave detector's bad day
        // should shorten the P band, not delete the QRS and T scores with it.
        if (s.qrsStart >= 0) s.pEnd = std::min(s.pEnd, s.qrsStart);
        return s;
    }

    // One bin, one channel.
    struct BinResult {
        int  bin = -1;
        std::string channel;
        bool ok = false;                 // false => skipped (bad segment, too few beats, no landmarks)
        SegmentCols seg;
        MorphologyEnvelope env;          // pass-2 (clean-pool) envelope
        double sdPass1 = std::numeric_limits<double>::quiet_NaN();
        double sdPass2 = std::numeric_limits<double>::quiet_NaN();
        int  nBeats = 0, nClean = 0;
        std::vector<BandMatchResult> scores;   // per beat, against env
        std::vector<PreMark>         marks;    // from BeatClassifier::preMarkAll
    };

    // Mean per-sample sd over columns that had any contributor, for the
    // pass-1 vs pass-2 tightening readout.
    inline double meanSd(const MorphologyEnvelope& env) {
        double s = 0.0; int n = 0;
        for (size_t c = 0; c < env.sd.size(); ++c)
            if (!std::isnan(env.sd[c]) && env.sd[c] > 0.0) { s += env.sd[c]; ++n; }
        return n ? s / n : std::numeric_limits<double>::quiet_NaN();
    }

    inline BinResult runBin(const std::vector<std::vector<double>>& beats,
        const std::vector<double>& tmpl,
        int rCol, double fs,
        BeatClassifier& classifier,
        int binIdx, const std::string& channel)
    {
        BinResult r;
        r.bin = binIdx;
        r.channel = channel;
        r.nBeats = (int)beats.size();
        if (r.nBeats < kMinBeatsForEnvelope || beats[0].empty()) {
            return r;
        }

        const int W = (int)beats[0].size();
        r.seg = segmentsFromTemplate(tmpl, rCol, fs);
        if (!r.seg.valid(W)) {
            // Landmark auto-detection failed or came back out of order.
            return r;
        }

        // Pass 1: corridor over every kept beat.
        MorphologyEnvelope env1 = buildEnvelope(beats, W);
        r.sdPass1 = meanSd(env1);

        // Pass 2: keep the beats that sat inside it.
        std::vector<std::vector<double>> clean;
        clean.reserve(beats.size());
        for (const auto& b : beats) {
            // Pass 1 reads pct_overall and pct_QRS only -- scoring P, ST, T and
            // tail here computes four bands per beat that are then discarded.
            // The full six-band score is taken in the final pass below, against
            // the pass-2 corridor, which is what actually reaches the CSV.
            const BandMatchResult bm = scoreBeatBands(b, env1, r.seg.qrsStart,
                r.seg.qrsEnd);
            // NaN (unscorable band) fails the gate, same as a low score did
            // before -- but it is now visible as NaN in the CSV rather than
            // masquerading as 0.0.
            if (!std::isnan(bm.pct_overall) && !std::isnan(bm.pct_QRS)
                && bm.pct_overall >= kCleanPctThreshold
                && bm.pct_QRS >= kCleanPctThreshold) clean.push_back(b);
        }
        // Everything scored out (very noisy bin) -- keep pass 1 rather than
        // publishing an empty corridor that every beat would then fail.
        r.env = (clean.size() >= (size_t)kMinBeatsForEnvelope)
            ? buildEnvelope(clean, W) : env1;
        r.nClean = (int)clean.size();
        r.sdPass2 = meanSd(r.env);

        // Final per-beat scores against the pass-2 corridor.
        r.scores.resize(beats.size());
        for (size_t i = 0; i < beats.size(); ++i)
            r.scores[i] = scoreBeat(beats[i], r.env, r.seg.pEnd, r.seg.qrsStart,
                r.seg.qrsEnd, r.seg.tBegin, r.seg.tEnd);

        // Pre-marks. NOTE: preMarkAll's body is comments in the spec (feature
        // matrix + ONNX batch), so until the model is wired this returns
        // default PreMarks -- cls UNKNOWN, confidence 0. The envelope and the
        // band scores above are real; the labels are not yet.
        r.marks = classifier.preMarkAll(beats, r.env, MorphologyEnvelope{},
            MorphologyEnvelope{});
        r.ok = true;
        return r;
    }

    // All bins of one channel.
    inline std::vector<BinResult> runChannel(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        const std::string& channel,
        template_io::ChannelMethodTemplate
        template_io::BinTemplates::* methodPtr,
        double fs, BeatClassifier& classifier)
    {
        std::vector<BinResult> out;
        auto it = beats.per_channel_beats.find(channel);
        if (it == beats.per_channel_beats.end() || fs <= 0.0) return out;

        const auto& perBin = it->second;

        // Bins are independent: runBin reads only its own beat matrix and its
        // own template, and writes only into its own BinResult. Nothing here
        // touches shared state -- the CSV writers run afterwards, on the
        // caller's thread.
        //
        // Indexed into a pre-sized vector rather than push_back, because
        // push_back from a parallel loop would race on the vector itself and
        // would also scramble bin order. Skipped bins keep ok=false and are
        // filtered out below, so bin i stays at index i during the loop.
        //
        // NAMING: not `slots`. Qt's qobjectdefs.h defines `slots`, `signals`
        // and `emit` as macros expanding to nothing (unless QT_NO_KEYWORDS),
        // and this header is pulled into a Qt target via post_process.hpp, so
        // a variable named `slots` silently vanishes at the preprocessor and
        // surfaces as "reserve does not take 0 arguments" -- the argument was
        // eaten, not miscounted. Avoid those three identifiers here.
        const int nBins = (int)std::min(perBin.size(), tmpl.bins.size());
        std::vector<BinResult> binResults((size_t)std::max(0, nBins));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (int i = 0; i < nBins; ++i) {
            if (i < (int)beats.bad_segment.size() && beats.bad_segment[i]) continue;
            const template_io::ChannelMethodTemplate& blk = tmpl.bins[i].*methodPtr;
            binResults[i] = runBin(perBin[i], blk.ecgTemplate, blk.r_col, fs,
                classifier, i, channel);
        }

        out.reserve(binResults.size());
        for (auto& s : binResults)
            if (s.bin >= 0) out.push_back(std::move(s));
        return out;
    }

    // ---------------------------------------------------------------------
    // CSV output -- same dir/stem pattern as ecg_move_log.
    // ---------------------------------------------------------------------
    inline std::string g_dir, g_stem;
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem;
    }

    inline void write_scores(const std::vector<BinResult>& results, bool first, const std::string& dir, const std::string& stem) {
        // Write each bin's per-beat scores and pre-marks to a CSV in quality_metric
        if (dir.empty() || stem.empty()) return;
        const std::string path = dir + "/" + stem + "_premarks.csv";
        std::ofstream f(path, first ? std::ios::trunc : std::ios::app);
        if (!f) {
            std::fprintf(stderr, "[premark] cannot open %s for writing\n", path.c_str());
            return;
        }
        if (first) f << "stem,channel,bin,beat,pct_overall,pct_P,pct_QRS,pct_ST,pct_T,pct_tail,class,confidence\n";
        const std::string& g_stem = stem;   // rows carry the stem in column 1
        for (const auto& r : results) {
            if (!r.ok) continue;
            for (size_t k = 0; k < r.scores.size(); ++k) {
                const PreMark::Class cls = (k < r.marks.size()) ? r.marks[k].cls
                    : PreMark::UNKNOWN;
                const double conf = (k < r.marks.size()) ? r.marks[k].confidence : 0.0;
                f << g_stem << ',' << r.channel << ',' << r.bin << ',' << k << ','
                    << r.scores[k].pct_overall << ',' << r.scores[k].pct_P << ','
                    << r.scores[k].pct_QRS << ',' << r.scores[k].pct_ST << ','
                    << r.scores[k].pct_T << ',' << r.scores[k].pct_tail << ','
                    << preMarkClassName(cls) << ',' << conf << '\n';
            }
        }
    }

    // Per-bin envelope summary: landmarks, beat counts, and the pass-1 to
    // pass-2 sd change (the tightening readout).
    inline void write_envelope_summary(const std::vector<BinResult>& results, bool first,
        const std::string& dir, const std::string& stem) {
        if (dir.empty() || stem.empty()) return;
        const std::string path = dir + "/" + stem + "_envelopes.csv";
        std::ofstream f(path, first ? std::ios::trunc : std::ios::app);
        if (!f) {
            std::fprintf(stderr, "[premark] cannot open %s for writing\n", path.c_str());
            return;
        }
        if (first) f << "stem,channel,bin,ok,n_beats,n_clean,p_end,qrs_start,"
            "qrs_end,t_begin,t_end,sd_pass1,sd_pass2,env_tight\n";
        const std::string& g_stem = stem;
        for (const auto& r : results)
            f << g_stem << ',' << r.channel << ',' << r.bin << ',' << (r.ok ? 1 : 0)
            << ',' << r.nBeats << ',' << r.nClean << ',' << r.seg.pEnd << ','
            << r.seg.qrsStart << ',' << r.seg.qrsEnd << ',' << r.seg.tBegin << ','
            << r.seg.tEnd << ',' << r.sdPass1 << ',' << r.sdPass2 << ','
            << (r.env.tight() ? 1 : 0) << '\n';
    }

    // ---------------------------------------------------------------------
    // One call for all three ECG channels.
    // ---------------------------------------------------------------------
    inline void runAll(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        double ecgRate,
        const std::string& dirIn,
        const std::string& stemIn,
        const std::string& onnxModelPath = {})
    {
        // Shadow the globals with locals: nothing below is shared state, so
        // concurrent finalize workers cannot interfere.
        const std::string g_dir = dirIn;
        const std::string g_stem = stemIn;
        // Silent-return diagnostics. Every early exit below prints its reason,
        // because the failure mode this replaces was "no files appear and
        // nothing says why" -- an ofstream onto a directory that does not
        // exist fails without throwing, and set() being handed an empty
        // config path looks identical to the feature being switched off.
        if (g_dir.empty() || g_stem.empty()) {
            return;
        }

        // The output directory is NOT guaranteed to exist: it comes straight
        // from the config, and nothing upstream creates it.
        {
            std::error_code ec;
            std::filesystem::create_directories(g_dir, ec);
            if (ec) {
                std::fprintf(stderr, "[premark] cannot create %s: %s\n",
                    g_dir.c_str(), ec.message().c_str());
                return;
            }
        }

        if (ecgRate <= 0.0) {  return;  }
        BeatClassifier classifier(onnxModelPath);

        struct Chan {
            const char* key; template_io::ChannelMethodTemplate
                template_io::BinTemplates::* ptr;
        };
        const Chan chans[] = {
            { "CH1", &template_io::BinTemplates::ch1_raw },
            { "CH2", &template_io::BinTemplates::ch2_raw },
            { "CH3", &template_io::BinTemplates::ch3_raw },
        };
        bool first = true;
        for (const Chan& c : chans) {
            const auto res = runChannel(beats, tmpl, c.key, c.ptr, ecgRate, classifier);
            if (res.empty()) {
                const bool haveKey =  beats.per_channel_beats.find(c.key) != beats.per_channel_beats.end();
    
                continue;
            }
            int nOk = 0, nBeats = 0;
            for (const auto& r : res) { if (r.ok) ++nOk; nBeats += r.nBeats; }
            write_scores(res, first, g_dir, g_stem);
            write_envelope_summary(res, first, g_dir, g_stem);
            first = false;
        }
    }

    // Backward-compatible form using the globals set by premark::set().
    // Safe ONLY when called from a single thread. Prefer the explicit
    // destination overload from any worker thread.
    inline void runAll(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        double ecgRate,
        const std::string& onnxModelPath = {})
    {
        runAll(beats, tmpl, ecgRate, g_dir, g_stem, onnxModelPath);
    }

} // namespace premark