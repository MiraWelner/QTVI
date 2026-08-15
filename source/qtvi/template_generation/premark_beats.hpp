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
        int pEnd = -1, qrsStart = -1, qrsEnd = -1, tEnd = -1;
        bool valid(int W) const {
            return pEnd >= 0 && qrsStart >= pEnd && qrsEnd > qrsStart
                && tEnd > qrsEnd && tEnd <= W;
        }
    };

    // Auto-detect the four boundaries off one bin's template.
    inline SegmentCols segmentsFromTemplate(const std::vector<double>& tmpl,
        int rCol, double fs)
    {
        SegmentCols s;
        if (tmpl.empty() || rCol < 0 || fs <= 0.0) return s;
        const double qOn = FeatureMarks::compute_q_onset(tmpl, fs, rCol);
        const double jPt = FeatureMarks::compute_j_point(tmpl, fs, rCol);
        const double tBeg = FeatureMarks::compute_t_begin(tmpl, fs, rCol, jPt);
        const double tEnd = FeatureMarks::compute_t_end(tmpl, fs, rCol, tBeg);
        const int    pEnd = FeatureMarks::detect_p_end(tmpl, rCol, fs);

        s.qrsStart = (qOn >= 0.0) ? (int)std::lround(qOn) : -1;
        s.qrsEnd = (jPt >= 0.0) ? (int)std::lround(jPt) : -1;
        s.tEnd = (tEnd >= 0.0) ? (int)std::lround(tEnd) : -1;
        // No P-end detected (absent or flat P) -- fall back to the QRS onset,
        // which makes the P band the whole pre-QRS span rather than dropping
        // the band entirely.
        s.pEnd = (pEnd >= 0) ? pEnd : s.qrsStart;
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
            if (env.sd[c] > 0.0) { s += env.sd[c]; ++n; }
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
            std::fprintf(stderr, "[premark]     bin %d %s: only %d beats "
                "(need %d)\n", binIdx, channel.c_str(), r.nBeats,
                kMinBeatsForEnvelope);
            return r;
        }

        const int W = (int)beats[0].size();
        r.seg = segmentsFromTemplate(tmpl, rCol, fs);
        if (!r.seg.valid(W)) {
            // Landmark auto-detection failed or came back out of order. Print
            // the columns so it is obvious which detector returned junk rather
            // than just losing the bin.
            std::fprintf(stderr, "[premark]     bin %d %s: bad landmarks "
                "pEnd=%d qrs=[%d,%d] tEnd=%d W=%d rCol=%d\n",
                binIdx, channel.c_str(), r.seg.pEnd, r.seg.qrsStart,
                r.seg.qrsEnd, r.seg.tEnd, W, rCol);
            return r;
        }

        // Pass 1: corridor over every kept beat.
        MorphologyEnvelope env1 = buildEnvelope(beats, W);
        r.sdPass1 = meanSd(env1);

        // Pass 2: keep the beats that sat inside it.
        std::vector<std::vector<double>> clean;
        clean.reserve(beats.size());
        for (const auto& b : beats) {
            const BandMatchResult bm = scoreBeat(b, env1, r.seg.pEnd, r.seg.qrsStart,
                r.seg.qrsEnd, r.seg.tEnd);
            if (bm.pct_overall >= kCleanPctThreshold
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
                r.seg.qrsEnd, r.seg.tEnd);

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
        out.reserve(perBin.size());
        for (size_t i = 0; i < perBin.size(); ++i) {
            if (i < beats.bad_segment.size() && beats.bad_segment[i]) continue;
            if (i >= tmpl.bins.size()) break;
            const template_io::ChannelMethodTemplate& blk = tmpl.bins[i].*methodPtr;
            out.push_back(runBin(perBin[i], blk.ecgTemplate, blk.r_col, fs,
                classifier, (int)i, channel));
        }
        return out;
    }

    // ---------------------------------------------------------------------
    // CSV output -- same dir/stem pattern as ecg_move_log.
    // ---------------------------------------------------------------------
    inline std::string g_dir, g_stem;
    inline void set(const std::string& dir, const std::string& stem) {
        g_dir = dir; g_stem = stem;
    }

    // Per-beat band scores + pre-mark. `first` truncates and writes the
    // header; later channels append.
    inline void write_scores(const std::vector<BinResult>& results, bool first) {
        if (g_dir.empty() || g_stem.empty()) return;
        const std::string path = g_dir + "/" + g_stem + "_premarks.csv";
        std::ofstream f(path, first ? std::ios::trunc : std::ios::app);
        if (!f) {
            std::fprintf(stderr, "[premark] cannot open %s for writing\n", path.c_str());
            return;
        }
        if (first) f << "stem,channel,bin,beat,pct_overall,pct_P,pct_QRS,pct_ST,"
            "pct_T,class,confidence\n";
        for (const auto& r : results) {
            if (!r.ok) continue;
            for (size_t k = 0; k < r.scores.size(); ++k) {
                const PreMark::Class cls = (k < r.marks.size()) ? r.marks[k].cls
                    : PreMark::UNKNOWN;
                const double conf = (k < r.marks.size()) ? r.marks[k].confidence : 0.0;
                f << g_stem << ',' << r.channel << ',' << r.bin << ',' << k << ','
                    << r.scores[k].pct_overall << ',' << r.scores[k].pct_P << ','
                    << r.scores[k].pct_QRS << ',' << r.scores[k].pct_ST << ','
                    << r.scores[k].pct_T << ',' << (int)cls << ',' << conf << '\n';
            }
        }
    }

    // Per-bin envelope summary: landmarks, beat counts, and the pass-1 to
    // pass-2 sd change (the tightening readout).
    inline void write_envelope_summary(const std::vector<BinResult>& results, bool first) {
        if (g_dir.empty() || g_stem.empty()) return;
        const std::string path = g_dir + "/" + g_stem + "_envelopes.csv";
        std::ofstream f(path, first ? std::ios::trunc : std::ios::app);
        if (!f) {
            std::fprintf(stderr, "[premark] cannot open %s for writing\n", path.c_str());
            return;
        }
        if (first) f << "stem,channel,bin,ok,n_beats,n_clean,p_end,qrs_start,"
            "qrs_end,t_end,sd_pass1,sd_pass2\n";
        for (const auto& r : results)
            f << g_stem << ',' << r.channel << ',' << r.bin << ',' << (r.ok ? 1 : 0)
            << ',' << r.nBeats << ',' << r.nClean << ',' << r.seg.pEnd << ','
            << r.seg.qrsStart << ',' << r.seg.qrsEnd << ',' << r.seg.tEnd << ','
            << r.sdPass1 << ',' << r.sdPass2 << '\n';
    }

    // ---------------------------------------------------------------------
    // One call for all three ECG channels.
    // ---------------------------------------------------------------------
    inline void runAll(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        double ecgRate,
        const std::string& onnxModelPath = {})
    {
        // Silent-return diagnostics. Every early exit below prints its reason,
        // because the failure mode this replaces was "no files appear and
        // nothing says why" -- an ofstream onto a directory that does not
        // exist fails without throwing, and set() being handed an empty
        // config path looks identical to the feature being switched off.
        if (g_dir.empty() || g_stem.empty()) {
            std::fprintf(stderr, "[premark] skipped: destination not set "
                "(dir='%s' stem='%s') -- check the config field passed to "
                "premark::set\n", g_dir.c_str(), g_stem.c_str());
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

        if (ecgRate <= 0.0) {
            std::fprintf(stderr, "[premark] skipped %s: ecg rate is %g\n",
                g_stem.c_str(), ecgRate);
            return;
        }
        std::fprintf(stderr, "[premark] %s -> %s | channels=%zu bins=%zu fs=%g\n",
            g_stem.c_str(), g_dir.c_str(), beats.per_channel_beats.size(),
            tmpl.bins.size(), ecgRate);

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
                const bool haveKey =
                    beats.per_channel_beats.find(c.key) != beats.per_channel_beats.end();
                std::fprintf(stderr, "[premark]   %s: no bins (key present=%d)\n",
                    c.key, haveKey ? 1 : 0);
                continue;
            }
            int nOk = 0, nBeats = 0;
            for (const auto& r : res) { if (r.ok) ++nOk; nBeats += r.nBeats; }
            std::fprintf(stderr, "[premark]   %s: %zu bins, %d scored, %d beats\n",
                c.key, res.size(), nOk, nBeats);
            write_scores(res, first);
            write_envelope_summary(res, first);
            first = false;
        }
    }

} // namespace premark