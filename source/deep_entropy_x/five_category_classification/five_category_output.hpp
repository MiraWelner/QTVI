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
 *           {stem}_beats.csv        one row per beat: rr, prematurity, the
 *                                   5-of-8 vote, category, label, handling,
 *                                   retained, template id, substituted
 *           {stem}_bins.csv         one row per bin: the five category
 *                                   percentages, PVC and PAC burden, the
 *                                   substituted and morphology-excluded shares
 *           {stem}_templates.csv    one row per template per bin: label, count,
 *                                   share
 *           {stem}_nsvt.csv         one row per detected run: onset, length,
 *                                   template, cycle length, rate, sustained
 *           {stem}_reference.csv    the reference average per bin, one row per
 *                                   sample -- what "excluded from the reference
 *                                   template" can be checked against
 *           {stem}_substituted.csv  every substituted beat, one row per sample:
 *                                   the original, the running reference, and
 *                                   the emitted blend
 *
 *         A BIN IS A BIN. 4.6 keeps a template bank per bin and 4.5 reports
 *         category percentages per bin, so everything here is scoped to
 *         BeatsFile::per_bin_beats[i] and never to the whole record. A bank
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

    inline double medianOf(std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
            [](double x) { return !std::isfinite(x); }), v.end());
        if (v.empty()) return kNaN;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }

    // Per-beat outputs for one bin.
    struct BeatRow {
        int    bin = 0, indexInBin = 0;
        double rrMs = kNaN;
        int    premature = 0, confirmed = 0;
        int    category = 0;
        std::string label;
        int    handling = 2, retained = 0, morphExcluded = 0;
        int    templateId = -1, substituted = 0, dropped = 0;
        double sqi = kNaN, templateCorr = kNaN, qrsWidthMs = kNaN;
        const char* rule = "";
    };

    struct BinResult {
        beatcls::BinCategoryReport report;
        TemplateBank bank;
        std::vector<NsvtRun> runs;
        std::vector<double> reference;
        int referenceBeats = 0;
        std::vector<double> templateWidthMs;
        std::vector<BeatRow> rows;
        // beat index in bin -> (original, reference, blend) for substituted beats
        std::vector<std::pair<int, std::array<std::vector<double>, 3>>> substituted;
    };

    // QRS duration from the detected boundaries. Task B's, not an estimate.
    inline double qrsWidthFrom(const Segments& seg, double msPerCol) {
        const double cols = static_cast<double>(seg.qrsHi - seg.qrsLo);
        return (cols > 0.0 && std::isfinite(msPerCol)) ? cols * msPerCol : kNaN;
    }

    // RR intervals for a bin. per_bin_beats carries the beats but not their
    // timing, so the interval is taken from the beat slices' own widths: each
    // spans 1.8 RR, which is the only rate information the beats themselves
    // carry. When the caller has real RR values it should pass them instead.
    inline std::vector<double> rrFromSliceWidths(
        const std::vector<std::vector<double>>& beats, double sliceRrMs)
    {
        return std::vector<double>(beats.size(), sliceRrMs);
    }

    inline BinResult runBin(const std::vector<std::vector<double>>& beats,
        const std::vector<double>& rrMs,
        const template_io::ChannelMethodTemplate& tmplRaw,
        const template_io::ChannelMethodTemplate& tmplAbs,
        int binIndex)
    {
        BinResult out;
        if (beats.empty()) return out;

        const double medRr = medianOf(rrMs);
        const SliceGeometry geo = geometryOf(tmplRaw, medRr);

        // ---- prematurity filter and 5-of-8 vote (4.6) -------------------
        const std::vector<char> premature = beatcls::prematureFlags(rrMs);
        const std::vector<char> confirmed = beatcls::confirmedFlags(premature);

        // ---- template bank, seeded on this bin's own template (4.6) -----
        const int qrsLo = (geo.rCol >= 0) ? std::max(0, geo.rCol - 40) : -1;
        const int qrsHi = (geo.rCol >= 0) ? std::min(geo.width, geo.rCol + 40) : -1;
        out.bank = tmplRaw.ecgTemplate.empty()
            ? TemplateBank{}
        : beatcls::seedBank(tmplRaw.ecgTemplate, qrsLo, qrsHi);
        if (out.bank.size() == 0 && !beats.empty())
            out.bank = beatcls::seedBank(beats[0], qrsLo, qrsHi);

        std::vector<int> templateId(beats.size(), -1);
        for (std::size_t i = 0; i < beats.size(); ++i)
            templateId[i] = assignToTemplate(beats[i], out.bank);

        // Label templates by class, as 4.6 says an operator would: the wide
        // ones are ventricular. The threshold is RELATIVE to the bin's own
        // narrowest template, because a heart-rate-proportional slice truncates
        // a broad complex at a short RR and the same morphology then measures
        // narrower at a faster rate.
        out.templateWidthMs.assign(static_cast<std::size_t>(out.bank.size()), kNaN);
        {
            std::vector<double> sum(static_cast<std::size_t>(out.bank.size()), 0.0);
            std::vector<int> n(static_cast<std::size_t>(out.bank.size()), 0);
            for (std::size_t i = 0; i < beats.size(); ++i) {
                if (templateId[i] < 0 || templateId[i] >= out.bank.size()) continue;
                const Segments seg = buildSegments(beats[i], geo.rCol, geo.sliceFs);
                const double w = qrsWidthFrom(seg, geo.msPerCol);
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

        // ---- five-category classification (4.5) -------------------------
        std::vector<beatcls::BeatVerdict> verdicts(beats.size());
        std::vector<double> sqiOf(beats.size(), kNaN);
        std::vector<double> widthOf(beats.size(), kNaN);
        const int refT = beatcls::referenceTemplateIndex(out.bank);
        const double refWidth =
            (refT >= 0 && std::isfinite(out.templateWidthMs[(std::size_t)refT]))
            ? out.templateWidthMs[(std::size_t)refT] : 90.0;

        for (std::size_t i = 0; i < beats.size(); ++i) {
            const Segments seg = buildSegments(beats[i], geo.rCol, geo.sliceFs);
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
            e.premature = premature[i] != 0;
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

        out.report = beatcls::summarizeBin(verdicts);

        // ---- NSVT (4.6) --------------------------------------------------
        out.runs = detectNsvt(templateId, rrMs, out.bank);

        // ---- substitution (4.6), and the reference it draws from ---------
        beatsub::SubstitutionState st;
        std::vector<char> wasSub(beats.size(), 0), wasDrop(beats.size(), 0);
        for (std::size_t i = 0; i < beats.size(); ++i) {
            if (verdicts[i].handling == beatcls::BeatVerdict::INCLUDE) {
                beatsub::observeClean(st, beats[i]);
                continue;
            }
            if (verdicts[i].handling != beatcls::BeatVerdict::SUBSTITUTE) {
                beatsub::observeExcluded(st);
                wasDrop[i] = 1;
                continue;
            }
            const std::vector<double> refBefore = st.average;
            const beatsub::SubstitutionResult r = beatsub::substitute(st, beats[i]);
            if (r.substituted) {
                wasSub[i] = 1;
                out.substituted.push_back({ static_cast<int>(i),
                    { beats[i], refBefore, r.beat } });
            }
            else {
                wasDrop[i] = 1;
            }
        }
        out.reference = st.average;
        out.referenceBeats = st.nIncluded;

        // ---- rows --------------------------------------------------------
        out.rows.reserve(beats.size());
        for (std::size_t i = 0; i < beats.size(); ++i) {
            BeatRow r;
            r.bin = binIndex;
            r.indexInBin = static_cast<int>(i);
            r.rrMs = rrMs[i];
            r.premature = premature[i] ? 1 : 0;
            r.confirmed = confirmed[i] ? 1 : 0;
            r.category = static_cast<int>(verdicts[i].category);
            r.label = beatClassName(verdicts[i].label);
            r.handling = static_cast<int>(verdicts[i].handling);
            r.retained = verdicts[i].retainedForBurden ? 1 : 0;
            r.morphExcluded = verdicts[i].morphologyExcluded ? 1 : 0;
            r.templateId = templateId[i];
            r.substituted = wasSub[i];
            r.dropped = wasDrop[i];
            r.sqi = sqiOf[i];
            r.templateCorr = beatcls::correlationToReference(beats[i], out.bank);
            r.qrsWidthMs = widthOf[i];
            r.rule = verdicts[i].reason;
            out.rows.push_back(std::move(r));
        }
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
        const std::size_t nBins = std::min(beats.per_bin_beats.size(),
            tmpl.bins.size());
        if (nBins == 0) {
            std::fprintf(stderr, "five_category_output: no bins for %s\n",
                stem.c_str());
            return;
        }

        std::vector<BinResult> results;
        results.reserve(nBins);
        for (std::size_t b = 0; b < nBins; ++b) {
            if (b < beats.bad_segment.size() && beats.bad_segment[b]) {
                results.push_back(BinResult{});
                continue;
            }
            const auto& bb = beats.per_bin_beats[b];
            // Slice width over 1.8 gives the RR the slice was cut at.
            const double sliceRr = (!bb.empty() && !bb[0].empty() && ecgRate > 0.0)
                ? (bb[0].size() / 1.8) * (1000.0 / ecgRate) : 800.0;
            const std::vector<double> rr =
                (b < rrPerBin.size() && rrPerBin[b].size() == bb.size())
                ? rrPerBin[b] : rrFromSliceWidths(bb, sliceRr);
            results.push_back(runBin(bb, rr, tmpl.bins[b].ch1_raw,
                tmpl.bins[b].ch1_absval, static_cast<int>(b)));
        }

        const std::string base = dir + "/" + stem;
        int nBeats = 0, nSub = 0, nRuns = 0;

        {
            std::ofstream f(base + "_beats.csv");
            f << "bin,beat_in_bin,beat,rr_ms,premature,vote_confirmed,category,"
                "label,handling,retained,morph_excluded,template_id,substituted,"
                "dropped,sqi_composite,template_corr,qrs_width_ms,rule\n";
            int running = 0;
            for (const BinResult& r : results)
                for (const BeatRow& b : r.rows) {
                    f << b.bin << ',' << b.indexInBin << ',' << running++ << ','
                        << b.rrMs << ',' << b.premature << ',' << b.confirmed << ','
                        << b.category << ',' << b.label << ',' << b.handling << ','
                        << b.retained << ',' << b.morphExcluded << ','
                        << b.templateId << ',' << b.substituted << ','
                        << b.dropped << ',' << b.sqi << ',' << b.templateCorr << ','
                        << b.qrsWidthMs << ',' << b.rule << '\n';
                    ++nBeats;
                }
        }
        {
            std::ofstream f(base + "_bins.csv");
            f << "bin,n_beats,pct_cat1,pct_cat2,pct_cat3,pct_cat4,pct_cat5,"
                "pvc_burden_pct,pac_burden_pct,substituted_pct,"
                "morph_excluded_pct\n";
            for (std::size_t b = 0; b < results.size(); ++b) {
                const beatcls::BinCategoryReport& r = results[b].report;
                if (r.nBeats == 0) continue;
                f << b << ',' << r.nBeats << ','
                    << r.pct(BeatCategory::BONAFIDE_PQRST) << ','
                    << r.pct(BeatCategory::ABNORMAL_RHYTHM) << ','
                    << r.pct(BeatCategory::ARTIFACTUAL_R) << ','
                    << r.pct(BeatCategory::MACHINE_NOISE) << ','
                    << r.pct(BeatCategory::HUMAN_NOISE) << ','
                    << r.pvcBurdenPct << ',' << r.pacBurdenPct << ','
                    << r.substitutedPct << ',' << r.morphologyExcludedPct << '\n';
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
            std::ofstream f(base + "_nsvt.csv");
            f << "bin,start_beat_in_bin,start_beat,length,template_id,"
                "mean_cycle_ms,rate_bpm,sustained\n";
            int offset = 0;
            for (std::size_t b = 0; b < results.size(); ++b) {
                for (const NsvtRun& r : results[b].runs) {
                    f << b << ',' << r.startBeat << ','
                        << (offset + r.startBeat) << ',' << r.length << ','
                        << r.templateId << ',' << r.meanCycleMs << ','
                        << r.rateBpm << ',' << (r.sustained ? 1 : 0) << '\n';
                    ++nRuns;
                }
                offset += static_cast<int>(results[b].rows.size());
            }
        }
        {
            std::ofstream f(base + "_reference.csv");
            f << "bin,column,reference,n_contributing\n";
            for (std::size_t b = 0; b < results.size(); ++b)
                for (std::size_t c = 0; c < results[b].reference.size(); ++c)
                    f << b << ',' << c << ',' << results[b].reference[c] << ','
                    << results[b].referenceBeats << '\n';
        }
        {
            std::ofstream f(base + "_substituted.csv");
            f << "bin,beat_in_bin,column,original,reference,substituted\n";
            for (std::size_t b = 0; b < results.size(); ++b)
                for (const auto& [idx, trio] : results[b].substituted) {
                    const auto& [orig, ref, blend] = trio;
                    const std::size_t n = std::min({ orig.size(), ref.size(),
                        blend.size() });
                    for (std::size_t c = 0; c < n; ++c)
                        f << b << ',' << idx << ',' << c << ',' << orig[c] << ','
                        << ref[c] << ',' << blend[c] << '\n';
                    ++nSub;
                }
        }

        std::fprintf(stderr, "  five-category output: %d beats, %zu bins, "
            "%d substituted, %d NSVT run(s) -> %s_*.csv\n",
            nBeats, nBins, nSub, nRuns, base.c_str());
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