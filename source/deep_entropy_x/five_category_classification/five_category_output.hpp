#pragma once
/**
 * @file   five_category_output.hpp
 * @brief  Runs Task C over a record's beats and writes its results as CSV.
 *         Sections 4.5 to 4.7.
 *
 *         Task C's headers compute and return; none of them writes a file. This
 *         is the writer. It takes the beats and templates the pipeline already
 *         has in memory (post_process's job.beats / job.tmpl), runs the
 *         classification, the template bank, NSVT detection and the
 *         substitution over each bin, and emits six CSVs into
 *         cfg.five_category_output.
 *
 *         Called from post_process.hpp beside premark::runAll, which has the
 *         same inputs in scope.
 *
 *         WHAT IT WRITES, per record stem:
 *
 *           {stem}_beats.csv        one row per beat: rr, assigned template,
 *                                   composite SQI, template correlation, QRS
 *                                   width, and ONE flag -- NORMAL, PVC
 *                                   (premature) or VOTED_PVC (flagged by the
 *                                   5-of-8 vote on its neighbours). No
 *                                   category, label, handling or substitution
 *                                   marker: one flag, one place.
 *           {stem}_bins.csv         one row per bin: the five category
 *                                   percentages, PVC and PAC burden, the
 *                                   morphology-excluded share, the premature
 *                                   and voted-premature counts and shares, and
 *                                   the substitution tallies -- n_flagged, the
 *                                   substituted and skipped shares. Skipped is
 *                                   flagged-but-not-substitutable: no clean
 *                                   neighbour in reach, or an ill-formed row.
 *           {stem}_templates.csv    one row per template per bin: label, count,
 *                                   share
 *           {stem}_substituted.csv  every substituted beat, one row per
 *                                   sample: the beat, the two clean
 *                                   neighbours' mean that replaced it, and
 *                                   the blend, with both neighbour indices
 *
 *         A BIN IS A BIN. 4.6 keeps a template bank per bin and 4.5 reports
 *         category percentages per bin, so everything here is scoped to
 *         BeatsFile::per_channel_beats[ch][i], never the whole record. A bank
 *         spanning a record whose rate changes splits one morphology across
 *         templates for no physiological reason.
 *
 *         WHAT IT DOES NOT INVENT. The composite SQI comes from Task A
 *         (sqi_ecg.hpp) and the segment boundaries from Task B via
 *         buildSegments. Where a piece of evidence needs a measurement neither
 *         provides, the field is LEFT UNSET rather than estimated:
 *         BeatEvidence treats absent evidence as absent, and a fabricated
 *         number is worse than a missing one -- a stand-in SQI reading 0.34 on
 *         clean beats put a whole record in the noise categories during
 *         development.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "config_file_handling/config_entry.hpp"
#include "template_generation/template_io.hpp"
#include "five_category_classification/five_categories.hpp"
#include "five_category_classification/beat_substitute.hpp"
#include "logging/sqi_ecg.hpp"

namespace five_category_output {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    // A bin's beats are sliced heart-rate-proportionally: 0.3 RR before R and
    // 1.5 RR after, resampled to a common width. So a column is not a fixed
    // number of milliseconds -- it depends on the bin's own rate -- and the
    // boundary detectors need the slice's effective rate, not the record's.
    struct SliceGeometry {
        int    width = 0;
        int    rCol = -1;
        double msPerCol = kNaN;
        double sliceFs = kNaN;
    };

    inline SliceGeometry geometryOf(const template_io::ChannelMethodTemplate& t,
        double medianRrMs)
    {
        SliceGeometry g;
        g.width = static_cast<int>(t.ecgTemplate.size());
        g.rCol = t.r_col;
        if (g.width > 0 && medianRrMs > 0.0) {
            g.msPerCol = 1.8 * medianRrMs / g.width;
            g.sliceFs = 1000.0 / g.msPerCol;
        }
        return g;
    }

    // The one flag {stem}_beats.csv carries. Section 4.6 has two ways a beat
    // is called ectopic and they are not the same evidence, so the flag
    // distinguishes them rather than collapsing both to a boolean:
    //
    //   PVC        the prematurity filter fired on this beat directly --
    //              RR(t) < 0.80 * median of the trailing ten.
    //   VOTED_PVC  the beat is not itself premature; the 5-of-8 vote flagged
    //              it from its neighbours. This is the middle of a run, where
    //              the trailing median has itself gone short so an intra-run
    //              beat stops reading as premature on its own.
    //   NORMAL     neither fired.
    //
    // PVC wins when both apply: direct evidence over inferred.
    //
    // Beats the pruning passes rejected are not represented here at all --
    // they are discarded in alignment.hpp and never reach this file. A
    // rhythm-flagged beat is exempt from that pruning, which is the whole
    // reason the verdict is assigned before the Tukey passes run.
    enum class RhythmFlag { NORMAL = 0, PVC = 1, VOTED_PVC = 2 };

    inline const char* rhythmFlagName(RhythmFlag f) {
        switch (f) {
        case RhythmFlag::PVC:       return "PVC";
        case RhythmFlag::VOTED_PVC: return "VOTED_PVC";
        default:                    return "NORMAL";
        }
    }

    // ---------------------------------------------------------------------
    // WHERE THE BEATS ACTUALLY ARE
    // ---------------------------------------------------------------------
    // The beats live in per_channel_beats["CH1"|"CH2"|"CH3"|"PPG"], filled by
    // build_templates.hpp from TemplateInfo::kept_beats_by_channel. This used
    // to read BeatsFile::per_bin_beats -- a field nothing ever wrote -- which
    // made nBins collapse to zero and Task C write nothing on every record.
    // That field has been deleted.
    //
    // Take the first ECG channel that has beats, with THAT channel's templates
    // and rhythm codes: beats and the template they are scored against must
    // come from the same lead.
    using MethodPtr = template_io::ChannelMethodTemplate template_io::BinTemplates::*;

    struct BeatSource {
        const std::vector<std::vector<std::vector<double>>>* beats = nullptr;
        // Rhythm codes for the same channel, [bin][beat], or null when the
        // producer did not fill them.
        const std::vector<std::vector<uint8_t>>* rhythm = nullptr;
        MethodPtr   raw = &template_io::BinTemplates::ch1_raw;
        MethodPtr   absval = &template_io::BinTemplates::ch1_absval;
        const char* key = "";
    };

    inline BeatSource pickBeatSource(const template_io::BeatsFile& beats) {
        BeatSource s;
        struct Chan { const char* key; MethodPtr raw; MethodPtr absval; };
        static const Chan chans[] = {
            { "CH1", &template_io::BinTemplates::ch1_raw,
                     &template_io::BinTemplates::ch1_absval },
            { "CH2", &template_io::BinTemplates::ch2_raw,
                     &template_io::BinTemplates::ch2_absval },
            { "CH3", &template_io::BinTemplates::ch3_raw,
                     &template_io::BinTemplates::ch3_absval },
        };
        for (const Chan& c : chans) {
            const auto it = beats.per_channel_beats.find(c.key);
            if (it == beats.per_channel_beats.end() || it->second.empty()) continue;
            s.beats = &it->second;
            s.raw = c.raw;
            s.absval = c.absval;
            s.key = c.key;
            const auto rt = beats.per_channel_rhythm.find(c.key);
            if (rt != beats.per_channel_rhythm.end()) s.rhythm = &rt->second;
            return s;
        }
        return s;   // .beats == nullptr: caller reports and returns
    }

    inline double medianOf(std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return kNaN;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }

    // Per-beat measurements, plus ONE flag.
    //
    // The flag is the rhythm call and nothing else: no category, no template
    // label, no handling, no substitution marker. Section 4.6's two routes to
    // "premature" are the only thing distinguished here, and everything about
    // morphology stays where the morphology is.
    struct BeatRow {
        int    bin = 0, indexInBin = 0;
        double rrMs = kNaN;
        int    templateId = -1;
        double sqi = kNaN, templateCorr = kNaN, qrsWidthMs = kNaN;
        RhythmFlag flag = RhythmFlag::NORMAL;
    };

    struct BinResult {
        beatcls::BinCategoryReport report;
        TemplateBank bank;
        std::vector<double> templateWidthMs;
        // Substitution tallies, from beatsub::BinResult. nFlagged is the
        // rhythm-flagged count (PVC or VOTED_PVC); the other two partition it.
        // nSkipped is the part that matters and was previously discarded: a
        // beat flagged but NOT substitutable -- no clean neighbour within
        // reach inside a long run, or an ill-formed row. Without it,
        // "flagged minus substituted" is not derivable from any output.
        int nFlagged = 0, nSubstitutedBeats = 0, nSkipped = 0;
        // The premature/voted split, counted from the per-beat flags rather
        // than from the rhythm vector, so this file and _beats.csv cannot
        // disagree about how many beats were premature.
        int nPremature = 0, nVoted = 0;
        std::vector<BeatRow> rows;
        // One entry per substituted beat: the beat replaced, the two clean
        // neighbours blended to replace it, and the three waveforms.
        std::vector<beatsub::Substitution> substituted;
    };

    // QRS duration from the detected boundaries. Task B's, not an estimate.
    inline double qrsWidthFrom(const Segments& seg, double msPerCol) {
        const double cols = static_cast<double>(seg.qrsHi - seg.qrsLo);
        return (cols > 0.0 && std::isfinite(msPerCol)) ? cols * msPerCol : kNaN;
    }

    // RR intervals for a bin, used for the rate-proportional GEOMETRY only
    // (ms-per-column, the QRS-width scale). The beat matrix carries no timing,
    // so the fallback is the slice width, which spans 1.8 RR.
    //
    // The prematurity verdict does NOT come from here. It is assigned in
    // alignment.hpp before the pruning and travels with the beats -- see
    // BeatSource::rhythm. Deriving it from a reconstructed RR series is what
    // made every beat read NORMAL.
    inline std::vector<double> rrFromSliceWidths(
        const std::vector<std::vector<double>>& beats, double sliceRrMs)
    {
        return std::vector<double>(beats.size(), sliceRrMs);
    }

    // `rhythm` is the per-beat verdict assigned in alignment.hpp before the
    // pruning (0 NORMAL, 1 PVC, 2 VOTED_PVC), parallel to `beats`. Pass empty
    // when the producer did not fill it: every beat then reads NORMAL, which
    // is honest -- there is no way to recover the verdict here.
    inline BinResult runBin(const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rrMs,
        const std::vector<uint8_t>& rhythm,
        const template_io::ChannelMethodTemplate& tmplRaw,
        const template_io::ChannelMethodTemplate& tmplAbs,
        int binIndex)
    {
        BinResult out;
        if (beats.empty()) return out;

        // ---- TEMPORARY INSTRUMENTATION --------------------------------
        // This is the phase that runs between [fast-phases] and the
        // "five-category output:" line, i.e. where the freeze is.
        using fc_clk = std::chrono::steady_clock;
        auto fc_prev = fc_clk::now();
        const auto fc_t0 = fc_prev;
        auto fc_lap = [&fc_prev, binIndex](const char* what) {
            const auto now = fc_clk::now();
            std::fprintf(stderr, "    [5cat bin %d] %-22s %7lld ms\n", binIndex, what,
                (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                    now - fc_prev).count());
            std::fflush(stderr);
            fc_prev = now;
            };

        const double medRr = medianOf(rrMs);
        const SliceGeometry geo = geometryOf(tmplRaw, medRr);

        // ---- prematurity filter and 5-of-8 vote (4.6) -------------------
        // The verdict is READ, not recomputed. alignment.hpp assigned it after
        // the slice and before the Tukey passes; recomputing it here from the
        // surviving beats' reconstructed RR is exactly the mistake that made
        // the filter see a constant interval and never fire.
        const bool haveRhythm = (rhythm.size() == beats.size());
        const auto isPvc = [&](std::size_t i) {
            return haveRhythm && rhythm[i] == 1u;
            };
        const auto isVoted = [&](std::size_t i) {
            return haveRhythm && rhythm[i] == 2u;
            };

        // ---- template bank, seeded on this bin's own template (4.6) -----
        const int qrsLo = (geo.rCol >= 0) ? std::max(0, geo.rCol - 40) : -1;
        const int qrsHi = (geo.rCol >= 0) ? std::min(geo.width, geo.rCol + 40) : -1;
        out.bank = tmplRaw.ecgTemplate.empty()
            ? TemplateBank{}
        : beatcls::seedBank(tmplRaw.ecgTemplate, qrsLo, qrsHi);
        if (out.bank.size() == 0 && !beats.empty())
            out.bank = beatcls::seedBank(beats[0], qrsLo, qrsHi);

        fc_lap("seedBank");

        std::vector<int> templateId(beats.size(), -1);
        for (std::size_t i = 0; i < beats.size(); ++i)
            templateId[i] = assignToTemplate(beats[i], out.bank);
        fc_lap("assignToTemplate xN");

        // Boundaries from the TEMPLATE, once, shared by every beat in the bin.
        //
        // buildSegments runs four transitionAnchor fits (three candidate models
        // each, BIC-selected, 50 Newton iterations in the sigmoid), which costs
        // 13-15 ms. Per beat that was 92% of this whole phase -- ~166 s for
        // 11.5k beats. Every beat here sits on the template's axis: geo.rCol
        // comes from geometryOf(tmplRaw, ...) and the slice geometry puts each
        // beat's R at that column by construction. So the windows are the same
        // windows, and a sub-sample model fit per individual noisy beat is
        // template-grade work on beat-grade data.
        //
        // It is also the more defensible measurement: the per-beat SQI compares
        // beats against each other, and fixed windows are what make those
        // numbers comparable.
        std::vector<Segments> segOf(beats.size());
        {
            const Segments shared = tmplRaw.ecgTemplate.empty()
                ? buildSegments(beats[0], geo.rCol, geo.sliceFs)
                : buildSegments(tmplRaw.ecgTemplate, geo.rCol, geo.sliceFs);
            std::fill(segOf.begin(), segOf.end(), shared);
        }
        fc_lap("buildSegments xN");

        // Every assignment above marked its template's corridor dirty rather
        // than rebuilding it; rebuild them now, once, before the only consumer
        // (scoreAgainstTemplate, below) reads them.
        beatcls::refreshEnvs(out.bank);
        fc_lap("refreshEnvs");

        out.templateWidthMs.assign(static_cast<std::size_t>(out.bank.size()), kNaN);
        {
            std::vector<double> sum(static_cast<std::size_t>(out.bank.size()), 0.0);
            std::vector<int> n(static_cast<std::size_t>(out.bank.size()), 0);
            for (std::size_t i = 0; i < beats.size(); ++i) {
                if (templateId[i] < 0 || templateId[i] >= out.bank.size()) continue;
                const double w = qrsWidthFrom(segOf[i], geo.msPerCol);
                if (!std::isfinite(w)) continue;
                sum[static_cast<std::size_t>(templateId[i])] += w;
                ++n[static_cast<std::size_t>(templateId[i])];
            }
            double narrowest = std::numeric_limits<double>::infinity();
            for (int t = 0; t < out.bank.size(); ++t) {
                const double mean = n[static_cast<std::size_t>(t)]
                    ? sum[static_cast<std::size_t>(t)] / n[static_cast<std::size_t>(t)]
                    : kNaN;
                out.templateWidthMs[static_cast<std::size_t>(t)] = mean;
                if (n[static_cast<std::size_t>(t)] >= 2 && std::isfinite(mean))
                    narrowest = std::min(narrowest, mean);
            }
            if (!std::isfinite(narrowest)) narrowest = 90.0;   // nominal sinus QRS
            for (int t = 0; t < out.bank.size(); ++t) {
                const double mean = out.templateWidthMs[static_cast<std::size_t>(t)];
                const bool wide = std::isfinite(mean)
                    && mean > 110.0 && mean > 1.35 * narrowest;
                beatcls::confirmTemplateLabel(out.bank, t,
                    wide ? BeatClass::PVC_A : BeatClass::SINUS);
            }
        }

        fc_lap("template widths/labels");

        // ---- five-category classification (4.5) -------------------------
        std::vector<beatcls::BeatVerdict> verdicts(beats.size());
        std::vector<double> sqiOf(beats.size(), kNaN);
        std::vector<double> widthOf(beats.size(), kNaN);
        const int refT = beatcls::referenceTemplateIndex(out.bank);
        const double refWidth =
            (refT >= 0 && std::isfinite(out.templateWidthMs[(std::size_t)refT]))
            ? out.templateWidthMs[(std::size_t)refT] : 90.0;

        for (std::size_t i = 0; i < beats.size(); ++i) {
            const Segments& seg = segOf[i];
            beatcls::BeatEvidence e;

            // Task A: the composite SQI and the template correlation.
            if (!tmplRaw.ecgTemplate.empty()) {
                std::vector<double> absTmpl = tmplAbs.ecgTemplate;
                if (absTmpl.size() != tmplRaw.ecgTemplate.size()) {
                    absTmpl.resize(tmplRaw.ecgTemplate.size());
                    for (std::size_t c = 0; c < absTmpl.size(); ++c)
                        absTmpl[c] = std::fabs(tmplRaw.ecgTemplate[c]);
                }
                const BeatSQI q = computeEcgSQI(beats[i], tmplRaw.ecgTemplate,
                    absTmpl, seg, /*motionFlag=*/-1, geo.sliceFs);
                e.sqiComposite = q.composite;
                e.templateCorr = q.templateCorr;
                sqiOf[i] = q.composite;
            }
            else {
                e.templateCorr = beatcls::correlationToReference(beats[i], out.bank);
            }

            // Task B: QRS duration and whether a P segment was found.
            const double w = qrsWidthFrom(seg, geo.msPerCol);
            widthOf[i] = w;
            if (std::isfinite(w) && refWidth > 0.0) e.qrsWidthRatio = w / refWidth;
            e.pWavePresent = (seg.pHi > seg.pLo);

            // Band match against the beat's own template.
            const int tid = (templateId[i] >= 0) ? templateId[i] : refT;
            if (tid >= 0) {
                const BandMatchResult bm =
                    beatcls::scoreAgainstTemplate(beats[i], out.bank, tid);
                e.pctBandOverall = bm.pct_overall;
                e.pctBandQRS = beatcls::matchScoreOf(bm);
            }

            // Rhythm.
            // Rhythm. The vote WIDENS the prematurity flag rather than
            // filtering it (4.6): inside a run of ectopy the trailing median
            // has itself gone short, so an intra-run beat stops reading as
            // premature against its own neighbours and only the vote catches
            // it. Previously `confirmed` was computed, written to the CSV and
            // never acted on, which left exactly those mid-run beats
            // classified as sinus.
            e.premature = isPvc(i) || isVoted(i);
            {
                std::vector<double> w10;
                for (std::size_t k = (i >= 10 ? i - 10 : 0); k < i; ++k)
                    w10.push_back(rrMs[k]);
                const double med = medianOf(w10);
                if (std::isfinite(med) && med > 0.0) e.rrRatio = rrMs[i] / med;
            }
            // baselineDrift, noiseHfFrac, clipFrac, motion: not measured here.
            e.motion = -1;

            verdicts[i] = beatcls::classifyBeat(e);
        }
        fc_lap("SQI+bandmatch+classify xN");

        out.report = beatcls::summarizeBin(verdicts);

        // ---- substitution (4.6): temporal continuity ---------------------
        // Driven by the RHYTHM FLAG and nothing else: a beat the prematurity
        // filter or the 5-of-8 vote flagged is replaced by a blend of the
        // nearest clean beat either side, keeping 1/8 of itself. The SQI band
        // does not enter into it -- that is Section 4.3's handling decision,
        // and a mid-SQI beat with no flagged neighbour has no continuity
        // problem to repair.
        std::vector<char> rhythmFlag(beats.size(), 0);
        for (std::size_t i = 0; i < beats.size(); ++i)
            rhythmFlag[i] = (isPvc(i) || isVoted(i)) ? 1 : 0;

        const beatsub::BinResult sub = beatsub::applyToBin(beats, rhythmFlag);
        out.substituted = sub.subs;
        out.nFlagged = sub.nFlagged;
        out.nSubstitutedBeats = sub.nSubstituted;
        out.nSkipped = sub.nSkipped;
        // Which beats were substituted is reported by {stem}_substituted.csv
        // and nowhere else -- deliberately not carried back into the per-beat
        // rows.


        // ---- rows --------------------------------------------------------
        out.rows.reserve(beats.size());
        for (std::size_t i = 0; i < beats.size(); ++i) {
            BeatRow r;
            r.bin = binIndex;
            r.indexInBin = static_cast<int>(i);
            r.rrMs = rrMs[i];
            r.templateId = templateId[i];
            r.sqi = sqiOf[i];
            r.templateCorr = beatcls::correlationToReference(beats[i], out.bank);
            r.qrsWidthMs = widthOf[i];
            r.flag = isPvc(i) ? RhythmFlag::PVC
                : (isVoted(i) ? RhythmFlag::VOTED_PVC : RhythmFlag::NORMAL);
            if (r.flag == RhythmFlag::PVC)            ++out.nPremature;
            else if (r.flag == RhythmFlag::VOTED_PVC) ++out.nVoted;
            out.rows.push_back(std::move(r));
        }
        fc_lap("rows xN");
        std::fprintf(stderr, "    [5cat bin %d] %zu beats, TOTAL %lld ms\n",
            binIndex, beats.size(),
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                fc_clk::now() - fc_t0).count());
        std::fflush(stderr);
        return out;
    }

    // ---------------------------------------------------------------------
    // the entry point
    // ---------------------------------------------------------------------
    //
    // `rrPerBin` is optional: pass the real RR series per bin when the caller
    // has it. Without it the interval is taken as the bin's slice width over
    // 1.8, which is the only rate information the beat slices themselves
    // carry -- adequate for morphology, but the prematurity filter needs real
    // intervals to mean anything, so pass them when you have them.
    inline void runAll(const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        double ecgRate,
        const std::string& dir,
        const std::string& stem,
        const std::vector<std::vector<double>>& rrPerBin = {})
    {
        if (dir.empty() || stem.empty()) {
            std::fprintf(stderr, "five_category_output: empty dir or stem; "
                "nothing written\n");
            return;
        }
        // The output directory is not guaranteed to exist -- an ofstream onto a
        // missing directory fails silently, which looks identical to the
        // feature being switched off.
        {
            std::error_code ec;
            std::filesystem::create_directories(dir, ec);
            if (ec) {
                std::fprintf(stderr, "five_category_output: cannot create %s: %s\n",
                    dir.c_str(), ec.message().c_str());
                return;
            }
        }
        const BeatSource src = pickBeatSource(beats);
        if (!src.beats) {
            std::fprintf(stderr, "five_category_output: %s has no beats in "
                "per_channel_beats[CH1/CH2/CH3]; nothing "
                "written\n", stem.c_str());
            return;
        }
        const std::size_t nBins = std::min(src.beats->size(), tmpl.bins.size());
        if (nBins == 0) {
            std::fprintf(stderr, "five_category_output: no bins for %s "
                "(%zu beat bins from %s, %zu template bins)\n", stem.c_str(),
                src.beats->size(), src.key, tmpl.bins.size());
            return;
        }

        int nRhythmMissing = 0;
        std::vector<BinResult> results;
        results.reserve(nBins);
        // ---- TEMPORARY INSTRUMENTATION ------------------------------------
        const auto fc_all0 = std::chrono::steady_clock::now();
        for (std::size_t b = 0; b < nBins; ++b) {
            if (b < beats.bad_segment.size() && beats.bad_segment[b]) {
                results.push_back(BinResult{});
                continue;
            }
            const auto& bb = (*src.beats)[b];
            // Slice width over 1.8 gives the RR the slice was cut at.
            const double sliceRr = (!bb.empty() && !bb[0].empty() && ecgRate > 0.0)
                ? (bb[0].size() / 1.8) * (1000.0 / ecgRate) : 800.0;
            const std::vector<double> rr =
                (b < rrPerBin.size() && rrPerBin[b].size() == bb.size())
                ? rrPerBin[b] : rrFromSliceWidths(bb, sliceRr);
            // Beats and templates from the SAME lead -- see pickBeatSource.
            static const std::vector<uint8_t> kNoRhythm;
            const std::vector<uint8_t>& rh =
                (src.rhythm && b < src.rhythm->size()) ? (*src.rhythm)[b]
                : kNoRhythm;
            // ABSENT FLAGS ARE NOT "ALL NORMAL". runBin has to treat an empty
            // rhythm vector as no verdict, and no verdict renders as NORMAL --
            // which is indistinguishable in the CSV from a record that
            // genuinely has no ectopy. Say so, per bin, rather than emitting a
            // confident column of NORMAL built on nothing.
            if (rh.size() != bb.size()) {
                ++nRhythmMissing;
                std::fprintf(stderr,
                    "  *** five_category_output: bin %zu has %zu rhythm codes "
                    "for %zu beats -- the flag column will read NORMAL for "
                    "every beat in it. per_channel_rhythm[%s] is empty or "
                    "mismatched; the producers are alignment.hpp (assigns), "
                    "create_ecg_templates.hpp (captures), make_averaged_"
                    "templates.hpp and build_templates.hpp (carry).\n",
                    b, rh.size(), bb.size(), src.key);
            }
            results.push_back(runBin(bb, rr, rh, tmpl.bins[b].*src.raw,
                tmpl.bins[b].*src.absval, static_cast<int>(b)));
        }
        std::fprintf(stderr, "[5cat] all bins: %lld ms\n",
            (long long)std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - fc_all0).count());
        std::fflush(stderr);

        const std::string base = dir + "/" + stem;
        int nBeats = 0, nSub = 0;

        {
            std::ofstream f(base + "_beats.csv");
            // Measurements, plus the rhythm flag: NORMAL, PVC (premature) or
            // VOTED_PVC (flagged by the 5-of-8 vote on its neighbours).
            f << "bin,beat_in_bin,beat,rr_ms,template_id,sqi_composite,"
                "template_corr,qrs_width_ms,flag\n";
            int running = 0;
            for (const BinResult& r : results)
                for (const BeatRow& b : r.rows) {
                    f << b.bin << ',' << b.indexInBin << ',' << running++ << ','
                        << b.rrMs << ',' << b.templateId << ',' << b.sqi << ','
                        << b.templateCorr << ',' << b.qrsWidthMs << ','
                        << rhythmFlagName(b.flag) << '\n';
                    ++nBeats;
                }
        }
        {
            std::ofstream f(base + "_bins.csv");
            f << "bin,n_beats,pct_cat1,pct_cat2,pct_cat3,pct_cat4,pct_cat5,"
                "pvc_burden_pct,pac_burden_pct,"
                "morph_excluded_pct,n_premature,n_voted,pct_premature,"
                "pct_voted,n_flagged,substituted_pct,skipped_pct\n";
            for (std::size_t b = 0; b < results.size(); ++b) {
                const beatcls::BinCategoryReport& r = results[b].report;
                const BinResult& br = results[b];
                if (r.nBeats == 0) continue;
                const double d = 100.0 / r.nBeats;
                f << b << ',' << r.nBeats << ','
                    << r.pct(BeatCategory::BONAFIDE_PQRST) << ','
                    << r.pct(BeatCategory::ABNORMAL_RHYTHM) << ','
                    << r.pct(BeatCategory::ARTIFACTUAL_R) << ','
                    << r.pct(BeatCategory::MACHINE_NOISE) << ','
                    << r.pct(BeatCategory::HUMAN_NOISE) << ','
                    << r.pvcBurdenPct << ',' << r.pacBurdenPct << ','
                    << r.morphologyExcludedPct << ','
                    // Both shares are of ALL beats in the bin, not of the
                    // flagged ones, so they are comparable with the category
                    // percentages beside them and sum to n_flagged * d.
                    // n_premature: the prematurity filter fired directly.
                    // n_voted:     the 5-of-8 vote fired and the filter did
                    //              not -- the middle of a run.
                    << br.nPremature << ',' << br.nVoted << ','
                    << br.nPremature * d << ',' << br.nVoted * d << ','
                    << br.nFlagged << ','
                    << br.nSubstitutedBeats * d << ','
                    << br.nSkipped * d << '\n';
            }
        }
        {
            std::ofstream f(base + "_templates.csv");
            f << "bin,n_templates,template_id,label,count,share_pct,"
                "mean_qrs_width_ms\n";
            for (std::size_t b = 0; b < results.size(); ++b) {
                const TemplateBank& bank = results[b].bank;
                for (int t = 0; t < bank.size(); ++t)
                    f << b << ',' << bank.size() << ',' << t << ','
                    << beatClassName(bank.labels[(std::size_t)t]) << ','
                    << bank.counts[(std::size_t)t] << ',' << bank.share(t) << ','
                    << results[b].templateWidthMs[(std::size_t)t] << '\n';
            }
        }
        {
            std::ofstream f(base + "_substituted.csv");
            // Which beats were substituted, and which beats were blended to
            // substitute them. left_beat / right_beat are indices into the
            // same bin, so the source beats can be pulled back out of
            // _beats.csv. neighbour_mean is the avgOld handed to
            // substituteBeat; substituted is its output.
            f << "bin,beat_in_bin,left_beat,right_beat,note,column,"
                "original,neighbour_mean,substituted\n";
            for (std::size_t b = 0; b < results.size(); ++b)
                for (const beatsub::Substitution& s : results[b].substituted) {
                    const std::size_t n = std::min({ s.original.size(),
                        s.avgOld.size(), s.blend.size() });
                    for (std::size_t c = 0; c < n; ++c)
                        f << b << ',' << s.beat << ',' << s.leftBeat << ','
                        << s.rightBeat << ',' << s.note << ',' << c << ','
                        << s.original[c] << ',' << s.avgOld[c] << ','
                        << s.blend[c] << '\n';
                    ++nSub;
                }
        }

        std::fprintf(stderr, "  five-category output: %d beats from %s, %zu "
            "bins, %d substituted -> %s_*.csv\n",
            nBeats, src.key, nBins, nSub, base.c_str());
        if (nRhythmMissing > 0) {
            std::fprintf(stderr, "  *** %d of %zu bins had NO rhythm verdict. "
                "The flag column in %s_beats.csv is NOT evidence of a normal "
                "record.\n", nRhythmMissing, nBins, base.c_str());
        }
    }

    inline void runAll(const config_entry& cfg,
        const template_io::BeatsFile& beats,
        const template_io::TemplateFile& tmpl,
        const std::string& stem,
        const std::vector<std::vector<double>>& rrPerBin = {})
    {
        runAll(beats, tmpl, cfg.ecg_upsample_rate, cfg.five_category_output,
            stem, rrPerBin);
    }

} // namespace five_category_output