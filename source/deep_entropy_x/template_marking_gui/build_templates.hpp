/**
* @file   build_templates.hpp
* @brief  Run the template-generation pipeline on peak-finding output
*         and produce a TemplateFile ready to write to disk.
*
*         Carries per-sample std for the ECG raw method (the one the
*         viewer displays) and for PPG, from the in-memory TemplateInfo
*         into the on-disk BinTemplates layout. The other three ECG
*         methods don't have std computed and write empty std vectors.
*/

#pragma once

#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <iostream>
#include <atomic>

#include "template_io.hpp"
#include "template_generation/make_averaged_templates.hpp"
#include "template_structs.hpp"
#include "template_generation/create_arterial_templates.hpp"
#include "template_marking_gui/alignment.hpp"   // align_beat_matrix, QAlignResult
#include "template_marking_gui/feature_marks.hpp"   // AnchorType, make_anchor_locator
#include "peak_finding/peakfinding_io.hpp"

namespace template_generation_detail {

    // Copy one template (+ optional std) into the on-disk method block.
    // The non-raw methods pass an empty tmpl_iqr and the on-disk std
    // field stays sz=0 (no payload). Both writer and reader handle that
    // uniformly, so there's only ever one code path.
    inline void copyMethod(template_io::ChannelMethodTemplate& dst,
        const std::vector<double>& tmpl,
        const std::vector<double>& tmpl_iqr,
        double alignment, int rCol)
    {
        dst.ecgTemplate = tmpl;
        dst.ecg_template_iqr = tmpl_iqr;
        dst.alignment_point = std::isnan(alignment) ? 0.0 : alignment;
        dst.r_col = rCol;
    }

    // Lift the per-channel banks out of TemplateInfo's string-keyed map into
    // BinTemplates' fixed array. The map keys mirror kept_beats_by_channel
    // ("CH1"/"CH2"/"CH3"); a channel absent from the map leaves an empty bank,
    // which reads downstream as "no bank for this channel" -- the same thing a
    // pre-v3 file produces, so nothing needs to distinguish the two cases.
    inline void copyBanks(template_io::BinTemplates& bt, const TemplateInfo& info)
    {
        static const char* kKeys[3] = { "CH1", "CH2", "CH3" };
        for (int c = 0; c < 3; ++c) {
            auto it = info.bank_by_channel.find(kKeys[c]);
            if (it != info.bank_by_channel.end())
                bt.ecg_bank[c] = it->second.bank;
        }
        // The PPG bank rides the SAME string-keyed map, under "PPG". Lifted
        // here rather than left in TemplateInfo for the same reason the ECG
        // banks are: a bank that never reaches BinTemplates dies with the
        // in-memory build, and nothing downstream can tell that apart from a
        // record with a single pulse morphology.
        auto pit = info.bank_by_channel.find("PPG");
        if (pit != info.bank_by_channel.end())
            bt.ppg_bank = pit->second.bank;
    }

    inline void packBin(template_io::BinTemplates& bt,
        const TemplateInfo& info, bool bad_segment)
    {
        bt.bad_segment = bad_segment;
        if (bad_segment) return;

        // Only the raw methods carry std; the other three pass an empty
        // vector via the default-constructed std::vector<double>{}.
        const std::vector<double> noStd;

        copyMethod(bt.ch1_raw, info.ch1.ecgTemplate_raw,
            info.ch1.ecgTemplate_raw_iqr,
            info.ch1.alignment_point_raw, info.ch1.r_col_raw);
        copyMethod(bt.ch1_squared, info.ch1.ecgTemplate_squared, noStd,
            info.ch1.alignment_point_squared, info.ch1.r_col_squared);
        copyMethod(bt.ch1_absval, info.ch1.ecgTemplate_absval, noStd,
            info.ch1.alignment_point_absval, info.ch1.r_col_absval);
        copyMethod(bt.ch1_unfiltered, info.ch1.ecgTemplate_unfiltered, noStd,
            info.ch1.alignment_point_unfiltered, info.ch1.r_col_unfiltered);

        copyMethod(bt.ch2_raw, info.ch2.ecgTemplate_raw,
            info.ch2.ecgTemplate_raw_iqr,
            info.ch2.alignment_point_raw, info.ch2.r_col_raw);
        copyMethod(bt.ch2_squared, info.ch2.ecgTemplate_squared, noStd,
            info.ch2.alignment_point_squared, info.ch2.r_col_squared);
        copyMethod(bt.ch2_absval, info.ch2.ecgTemplate_absval, noStd,
            info.ch2.alignment_point_absval, info.ch2.r_col_absval);
        copyMethod(bt.ch2_unfiltered, info.ch2.ecgTemplate_unfiltered, noStd,
            info.ch2.alignment_point_unfiltered, info.ch2.r_col_unfiltered);

        copyMethod(bt.ch3_raw, info.ch3.ecgTemplate_raw,
            info.ch3.ecgTemplate_raw_iqr,
            info.ch3.alignment_point_raw, info.ch3.r_col_raw);
        copyMethod(bt.ch3_squared, info.ch3.ecgTemplate_squared, noStd,
            info.ch3.alignment_point_squared, info.ch3.r_col_squared);
        copyMethod(bt.ch3_absval, info.ch3.ecgTemplate_absval, noStd,
            info.ch3.alignment_point_absval, info.ch3.r_col_absval);
        copyMethod(bt.ch3_unfiltered, info.ch3.ecgTemplate_unfiltered, noStd,
            info.ch3.alignment_point_unfiltered, info.ch3.r_col_unfiltered);

        bt.ppgTemplate = info.ppgTemplate;
        bt.ppg_template_iqr = info.ppg_template_iqr;

        // Per-channel + PPG slice counts (post drop-rules).
        bt.ch1_n_beats_raw = info.ch1.n_beats_raw;
        bt.ch2_n_beats_raw = info.ch2.n_beats_raw;
        bt.ch3_n_beats_raw = info.ch3.n_beats_raw;
        bt.ppg_n_beats = info.ppg_n_beats;
        bt.ppg_peak_col = info.ppg_peak_col;
        bt.ppg_onset_col = info.ppg_onset_col;

        // ---- Section 4.6 template banks ---------------------------------
        // TemplateInfo carries a whole ChannelOutput per channel (bank plus the
        // per-beat flags and assignment vectors); only the BANK reaches disk,
        // because the viewer needs templates and the per-beat vectors would
        // multiply the file size for data nothing in the GUI reads.
        //
        // Copied, not moved: packBin() and packBinFast() take info by const
        // reference, and the same TemplateInfo is also consumed by
        // bin_archive::writeBinFeatureArchive().
        copyBanks(bt, info);
    }

    // FAST pack: raw + unfiltered ECG blocks + PPG. Leaves the squared and
    // absval blocks default-empty for packBinSlow.
    inline void packBinFast(template_io::BinTemplates& bt,
        const TemplateInfo& info, bool bad_segment)
    {
        bt.bad_segment = bad_segment;
        if (bad_segment) return;

        copyMethod(bt.ch1_raw, info.ch1.ecgTemplate_raw,
            info.ch1.ecgTemplate_raw_iqr,
            info.ch1.alignment_point_raw, info.ch1.r_col_raw);
        copyMethod(bt.ch1_unfiltered, info.ch1.ecgTemplate_unfiltered, {},
            info.ch1.alignment_point_unfiltered, info.ch1.r_col_unfiltered);

        copyMethod(bt.ch2_raw, info.ch2.ecgTemplate_raw,
            info.ch2.ecgTemplate_raw_iqr,
            info.ch2.alignment_point_raw, info.ch2.r_col_raw);
        copyMethod(bt.ch2_unfiltered, info.ch2.ecgTemplate_unfiltered, {},
            info.ch2.alignment_point_unfiltered, info.ch2.r_col_unfiltered);

        copyMethod(bt.ch3_raw, info.ch3.ecgTemplate_raw,
            info.ch3.ecgTemplate_raw_iqr,
            info.ch3.alignment_point_raw, info.ch3.r_col_raw);
        copyMethod(bt.ch3_unfiltered, info.ch3.ecgTemplate_unfiltered, {},
            info.ch3.alignment_point_unfiltered, info.ch3.r_col_unfiltered);

        bt.ppgTemplate = info.ppgTemplate;
        bt.ppg_template_iqr = info.ppg_template_iqr;

        // Per-channel + PPG slice counts (post drop-rules).
        bt.ch1_n_beats_raw = info.ch1.n_beats_raw;
        bt.ch2_n_beats_raw = info.ch2.n_beats_raw;
        bt.ch3_n_beats_raw = info.ch3.n_beats_raw;
        bt.ppg_n_beats = info.ppg_n_beats;
        bt.ppg_peak_col = info.ppg_peak_col;
        bt.ppg_onset_col = info.ppg_onset_col;

        // ---- Section 4.6 template banks ---------------------------------
        // TemplateInfo carries a whole ChannelOutput per channel (bank plus the
        // per-beat flags and assignment vectors); only the BANK reaches disk,
        // because the viewer needs templates and the per-beat vectors would
        // multiply the file size for data nothing in the GUI reads.
        //
        // Copied, not moved: packBin() and packBinFast() take info by const
        // reference, and the same TemplateInfo is also consumed by
        // bin_archive::writeBinFeatureArchive().
        copyBanks(bt, info);
    }

    // SLOW pack: squared + absval blocks onto an already fast-packed bin.
    inline void packBinSlow(template_io::BinTemplates& bt,
        const TemplateInfo& info)
    {
        if (bt.bad_segment) return;
        const std::vector<double> noStd;
        copyMethod(bt.ch1_squared, info.ch1.ecgTemplate_squared, noStd,
            info.ch1.alignment_point_squared, info.ch1.r_col_squared);
        copyMethod(bt.ch1_absval, info.ch1.ecgTemplate_absval, noStd,
            info.ch1.alignment_point_absval, info.ch1.r_col_absval);
        copyMethod(bt.ch2_squared, info.ch2.ecgTemplate_squared, noStd,
            info.ch2.alignment_point_squared, info.ch2.r_col_squared);
        copyMethod(bt.ch2_absval, info.ch2.ecgTemplate_absval, noStd,
            info.ch2.alignment_point_absval, info.ch2.r_col_absval);
        copyMethod(bt.ch3_squared, info.ch3.ecgTemplate_squared, noStd,
            info.ch3.alignment_point_squared, info.ch3.r_col_squared);
        copyMethod(bt.ch3_absval, info.ch3.ecgTemplate_absval, noStd,
            info.ch3.alignment_point_absval, info.ch3.r_col_absval);
    }

}  // namespace template_generation_detail

// Carries the in-memory state from the fast build to the slow merge.
struct FastTemplateBuild {
    template_io::TemplateFile tmpl;   // raw/unfiltered/ppg blocks + their SAECG
    template_io::BeatsFile beats;     // ch1 raw kept beats (final)
    std::vector<TemplateInfo> info;   // retained so mergeTemplatesSlow can pack squared/absval
};

// FAST build: produces a TemplateFile the viewer can open immediately
// (raw/unfiltered/ppg + their SAECG averages) and the final beats file.
// The squared/absval per-bin blocks and their SAECG entries stay empty
// until mergeTemplatesSlow runs.
inline FastTemplateBuild
buildTemplatesAndBeatsFast(const std::vector<output_binfile_data>& peakResults,
    const SignalRates& rates)
{
    using namespace template_generation_detail;

    FastTemplateBuild out;
    out.info = GenerateTemplatesFast(peakResults, rates);

    out.tmpl.bins.resize(peakResults.size());
    for (size_t i = 0; i < peakResults.size(); ++i) {
        bool bad = peakResults[i].bad_segment;
        const TemplateInfo& info = (i < out.info.size()) ? out.info[i] : TemplateInfo{};
        packBinFast(out.tmpl.bins[i], info, bad);
        // (Per-channel n_beats fields are populated inside packBinFast from
        // the TemplateInfo's own counts -- no fallback needed here.)
    }

    // Arterial background-context templates (ABP / ART / ART_PULM). All
    // three are now R-anchored, same as PPG (CreatePulseTemplates --
    // borrows ch1.raw ECG R-peaks, same [R_i-pad, R_{i+1}+pad] slicer).
    // The old foot-anchored self-detection path (build_arterial_template_
    // foot_anchored / CreateArterialTemplates) is no longer called from
    // here for any pulse channel; it's kept in create_arterial_templates.hpp
    // only in case a future channel still needs self-detected anchoring.
    // Present-only; a channel with rate=0 in SignalRates yields an empty
    // result and is silently skipped when packed into the bins.
    {
        auto abp = CreatePulseTemplates(
            peakResults, &output_binfile_data::abpSignal, rates.ecg, rates.abp);
        auto art = CreatePulseTemplates(
            peakResults, &output_binfile_data::artSignal, rates.ecg, rates.art);
        auto artp = CreatePulseTemplates(
            peakResults, &output_binfile_data::artPulmSignal, rates.ecg, rates.artPulm);
        for (size_t i = 0; i < out.tmpl.bins.size(); ++i) {
            if (out.tmpl.bins[i].bad_segment) continue;
            if (i < abp.templates.size()) {
                out.tmpl.bins[i].abpTemplate = std::move(abp.templates[i]);
                out.tmpl.bins[i].abpTemplate_iqr = std::move(abp.iqrs[i]);
            }
            if (i < art.templates.size()) {
                out.tmpl.bins[i].artTemplate = std::move(art.templates[i]);
                out.tmpl.bins[i].artTemplate_iqr = std::move(art.iqrs[i]);
            }
            if (i < artp.templates.size()) {
                out.tmpl.bins[i].artPulmTemplate = std::move(artp.templates[i]);
                out.tmpl.bins[i].artPulmTemplate_iqr = std::move(artp.iqrs[i]);
            }
        }

        // Retain each arterial channel's snips into per_channel_beats.
        auto stashArt = [&](const char* name, PPGTemplatesResult& r) {
            auto& dst = out.beats.per_channel_beats[name];
            if (dst.size() < r.kept.size()) dst.resize(r.kept.size());
            for (size_t i = 0; i < r.kept.size(); ++i)
                if (i >= out.tmpl.bins.size() || !out.tmpl.bins[i].bad_segment)
                    dst[i] = std::move(r.kept[i]);
            };
        stashArt("ABP", abp);
        stashArt("ART", art);
        stashArt("ART_PULM", artp);
    }

    // Beats: assemble per-channel retained snips (CH1/CH2/CH3/PPG from the
    // TemplateInfo map; arterial channels already stashed above) so every
    const size_t nb = out.info.size();
    out.beats.bad_segment.resize(nb, false);
    auto ensureBins = [&](const std::string& ch)
        -> std::vector<std::vector<std::vector<double>>>&{
        auto& v = out.beats.per_channel_beats[ch];
        if (v.size() < nb) v.resize(nb);
        return v;
        };
    for (size_t i = 0; i < nb; ++i) {
        out.beats.bad_segment[i] = (i < peakResults.size())
            ? peakResults[i].bad_segment : false;
        if (out.beats.bad_segment[i]) continue;
        for (auto& kv : out.info[i].kept_beats_by_channel)
            ensureBins(kv.first)[i] = std::move(kv.second);
        for (auto& kv : out.info[i].kept_rhythm_by_channel) {
            auto& v = out.beats.per_channel_rhythm[kv.first];
            if (v.size() < nb) v.resize(nb);
            v[i] = std::move(kv.second);
        }
        for (auto& kv : out.info[i].ref_index_by_channel) {
            auto& v = out.beats.per_channel_ref_index[kv.first];
            if (v.size() < nb) v.resize(nb, -1);
            v[i] = kv.second;
        }
    }
    /*
    // DEBUG/TEST: corrupt ~50% of bin 0's beats with additive noise (remove when done).
    {
        const double sigma = 3.0;   // mV noise stddev
        static std::mt19937 rng(2025);
        std::normal_distribution<double> gauss(0.0, sigma);
        std::bernoulli_distribution coin(0.5);   // 50% of beats get hit
        const size_t targetBin = 0;
        for (auto& kv : out.beats.per_channel_beats) {
            auto& binsVec = kv.second;               // [bin][beat][sample]
            if (targetBin >= binsVec.size()) continue;
            for (auto& beat : binsVec[targetBin]) {
                if (!coin(rng)) continue;            // skip half, leave them clean
                for (double& s : beat)
                    if (!std::isnan(s)) s += gauss(rng);
            }
        }
    }
    */

    return out;
}

// Anchor name for the summary line below. Local to this file: post_process.hpp
// has its own anchorName, and this header must not depend on it.
inline const char* anchorName_bt(AnchorType a) {
    switch (a) {
    case AnchorType::R_PEAK:  return "R_PEAK";
    case AnchorType::P_ONSET: return "P_ONSET";
    case AnchorType::Q_ONSET: return "Q_ONSET";
    case AnchorType::J_POINT: return "J_POINT";
    }
    return "?";
}

inline void alignTemplatesFromCache(template_io::TemplateFile& tmpl, template_io::BeatsFile& beats, const SignalRates& rates, AnchorType anchor, bool forScoring = false)
{
    /* Re-align reuse PPG/arterial as-is; Re-align the R-pass beats.
    Because the alignment anchor is the median snippet, the median snippet doesn't move,
    so R stays on its column*/
    const double fs = rates.ecg;

    using MethodPtr =
        template_io::ChannelMethodTemplate template_io::BinTemplates::*;
    // chIdx addresses raw_anchors[tag][bin][chIdx]. absPtr is the matching
    // absval scalar block, used in scoring mode to store the co-framed
    // aligned absval template ( = column-median of |aligned raw beats| ).
    struct Chan { const char* key; MethodPtr ptr; MethodPtr absPtr; int chIdx; };
    const Chan channels[] = {
        { "CH1", &template_io::BinTemplates::ch1_raw, &template_io::BinTemplates::ch1_absval, 0 },
        { "CH2", &template_io::BinTemplates::ch2_raw, &template_io::BinTemplates::ch2_absval, 1 },
        { "CH3", &template_io::BinTemplates::ch3_raw, &template_io::BinTemplates::ch3_absval, 2 },
    };

    // This anchor's per-bin store. R_PEAK is the scalar base and is NOT put
    // here; only the re-aligned anchors accumulate. Sized to bins up front so
    // every (bin, channel) slot exists; unaligned slots stay default-empty
    // and readTemplateInfoBin falls back to the R base for them.
    const int anchorTag = static_cast<int>(anchor);
    std::vector<std::array<template_io::ChannelMethodTemplate, 3>>& store =
        tmpl.raw_anchors[anchorTag];
    if (store.size() != tmpl.bins.size())
        store.assign(tmpl.bins.size(), {});

    // The per-slot store for this anchor, sized the same way. Filled from the
    // SAME aligned beat matrix the whole-channel average comes from, so a slot
    // and its bin cannot disagree about what "this anchor" means.
    std::vector<std::array<std::vector<template_io::TemplateFile::BankSlotTemplate>, 3>>& slotStore =
        tmpl.bank_anchors[anchorTag];
    if (slotStore.size() != tmpl.bins.size())
        slotStore.assign(tmpl.bins.size(), {});

    // Per-anchor tallies, summarised in ONE line at the end instead of a
    // message per bin. Declared out here because the summary is printed after
    // the channel loop; atomic because the bin loop inside it is an omp
    // parallel-for.
    std::atomic<int> nRecovered{ 0 };   // target column taken from the template
    std::atomic<int> nNoAlign{ 0 };     // nothing moved: this bin equals R

    for (const auto& ch : channels) {
        auto it = beats.per_channel_beats.find(ch.key);
        if (it == beats.per_channel_beats.end()) continue;
        auto& perBin = it->second;
        const auto refIt = beats.per_channel_ref_index.find(ch.key);
        // Bins are independent (each iteration touches only store[i], bins[i],
        // perBin[i]), so align them in parallel -- this is the dominant cost
        // of an anchor step. int index for OpenMP.
        const int nBins = static_cast<int>(tmpl.bins.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(std::min(8, std::max(1, nBins)))
#endif
        for (int i = 0; i < nBins; ++i) {
            template_io::BinTemplates& bin = tmpl.bins[i];
            if (bin.bad_segment) continue;
            if ((size_t)i >= perBin.size() || perBin[i].empty()) continue;

            // The R base (blk) is the reference every anchor aligns FROM; it
            // is read-only here and never overwritten.
            const template_io::ChannelMethodTemplate& blk = bin.*(ch.ptr);
            if (blk.r_col < 0 || blk.ecgTemplate.empty()) continue;

            // Built here rather than after the reference selection, because
            // that selection now calls it: the fallback below needs to ask
            // whether the landmark is findable on this beat before deciding
            // which waveform to take the target column from.
            const int r_anchor = blk.r_col;
            AnchorLocator locate = make_anchor_locator(anchor, r_anchor, fs);

            int refIdx = (refIt != beats.per_channel_ref_index.end()
                && i < refIt->second.size()) ? refIt->second[i] : -1;
            const std::vector<double>* refP =
                (refIdx >= 0 && refIdx < (int)perBin[i].size())
                ? &perBin[i][refIdx] : &blk.ecgTemplate;

            // ---- LOCATE THE ANCHOR ON THE TEMPLATE WHEN THE BEAT FAILS ----
            //
            // The reference is ONE arbitrary beat, chosen for having the median
            // length. That is the right choice for framing -- it is why no beat
            // gets clipped -- but it is a poor signal to detect a P wave on: a
            // single beat carries the full noise, and compute_p_begin returned
            // -1 on roughly half the bins. align_beat_matrix wraps its whole
            // shift loop in `if (marker >= 0)`, so one failed detection on that
            // one beat silently turned the entire bin's alignment into a no-op
            // and produced an anchor template byte-identical to R.
            //
            // blk.ecgTemplate is the column-wise median of every beat in the
            // bin -- the same waveform, with the noise averaged down by sqrt(n).
            // If the landmark is findable anywhere it is findable there. So when
            // the beat-based reference fails, fall back to the template for the
            // TARGET COLUMN only; the beats themselves are still located
            // individually and shifted onto it.
            //
            // Only on failure, not always: the beat and the template can differ
            // in length, and the target column has to be inside the frame the
            // beats are padded to.
            if (!(locate(*refP) >= 0.0)
                && !blk.ecgTemplate.empty()
                && blk.ecgTemplate.size() == refP->size()
                && locate(blk.ecgTemplate) >= 0.0) {
                refP = &blk.ecgTemplate;
                ++nRecovered;   // recovered, not failed -- see the summary line
            }
            const std::vector<double>& ref_beat_of_median_length = *refP;


            alignment::aligned_beats q = alignment::align_beat_matrix(perBin[i], blk.r_col, fs, /*compute_iqr=*/true, ref_beat_of_median_length, locate);
            if (q.tmpl.empty()) continue;

            // DID ANY BEAT ACTUALLY MOVE? This is the only real failure
            // signal. align_beat_matrix wraps its whole shift loop in
            // `if (marker >= 0)`, and on failure still returns a valid-looking
            // template that is byte-identical to R -- so counting non-empty
            // templates cannot detect it.
            //
            // Not a median shift: the reference IS the median snippet, so the
            // median shift is zero by design, and a sub-sample shift never
            // appears in an integer median at all.
            {
                const double m0 = locate(ref_beat_of_median_length);
                bool anyMoved = false;
                for (const auto& bt : perBin[i]) {
                    const double mi = locate(bt);
                    if (mi >= 0.0 && std::abs(m0 - mi) >= 1e-3) { anyMoved = true; break; }
                }
                if (!anyMoved) ++nNoAlign;
            }

            // ROUNDING BOUNDARY. q.r_col is sub-sample (the anchor locators and
            // the R refinement in alignment.hpp both return doubles), but
            // ChannelMethodTemplate::r_col is int and is SERIALIZED as an int32
            // in the template .bin (see template_io.hpp's on-disk layout), so
            // widening it is a file-format change. Rounded here, once, visibly.
            const int alignedRcol = (q.r_col >= 0.0)
                ? static_cast<int>(std::lround(q.r_col)) : blk.r_col;

            // Store the aligned result in the per-anchor slot, leaving the R
            // base untouched so the NEXT anchor step still aligns from R.
            template_io::ChannelMethodTemplate& dst = store[i][ch.chIdx];
            dst.alignment_point = blk.alignment_point;
            dst.r_col = alignedRcol;

            if (forScoring) {
                // Project into the scalar ch*_raw (so writeEcgSQICsv reads the
                // anchor template) and write the co-framed aligned beats back
                // (so QC scores beats in the same frame). Caller owns copies.
                template_io::ChannelMethodTemplate& scalar = bin.*(ch.ptr);
                scalar.ecgTemplate = q.tmpl;          // copy: also stored below
                scalar.ecg_template_iqr = q.iqr;
                scalar.r_col = alignedRcol;

                // Co-framed aligned ABSVAL template: |shifted raw beat| IS the
                // shifted |raw| beat (abs is pointwise), so the column-wise
                // NaN-skipping median of |q.beats| equals what aligning the
                // absval beats on this anchor would produce -- exact, no
                // separate absval beat matrix needed. Write it into the absval
                // scalar so writeEcgSQICsv's chiSqAbs is scored co-framed.
                if (!q.beats.empty()) {
                    const size_t W = q.beats.front().size();
                    std::vector<double> absTmpl(W, std::numeric_limits<double>::quiet_NaN());
                    std::vector<double> col;
                    col.reserve(q.beats.size());
                    for (size_t c = 0; c < W; ++c) {
                        col.clear();
                        for (const auto& bt : q.beats) {
                            if (c < bt.size() && !std::isnan(bt[c])) col.push_back(std::abs(bt[c]));
                        }
                        if (col.empty()) continue;
                        std::sort(col.begin(), col.end());
                        const size_t m = col.size() / 2;
                        absTmpl[c] = (col.size() % 2 == 0)
                            ? 0.5 * (col[m - 1] + col[m]) : col[m];
                    }
                    template_io::ChannelMethodTemplate& absScalar = bin.*(ch.absPtr);
                    absScalar.ecgTemplate = std::move(absTmpl);
                    absScalar.r_col = alignedRcol;

                    perBin[i] = std::move(q.beats);
                }
            }

            // ---- PER-SLOT AVERAGES, FROM THESE SAME ALIGNED BEATS --------
            //
            // q.beats is the aligned matrix and was being discarded outside the
            // forScoring branch. The bank slots are row subsets of it, so their
            // aligned averages are column-wise reductions over those rows --
            // NOT a second alignment pass, which would give every slot its own
            // frame.
            //
            // BankTemplate::members is in this channel's LOCAL row space (see
            // projectToChannel), the same space q.beats rows are in, so no
            // index translation. An out-of-range member is skipped, not
            // clamped: a stale index must drop a beat rather than quietly
            // average the wrong one.
            //
            // Median plus a ddof=1 std, matching what *_iqr actually holds
            // despite its name. NaN-skipping per column, since aligned beats
            // are NaN-padded where a shift ran off the window.
            if (!q.beats.empty() && ch.chIdx >= 0 && ch.chIdx < 3) {
                const tbank::TemplateBank& bnk = bin.ecg_bank[ch.chIdx];
                auto& outSlots = slotStore[i][ch.chIdx];
                outSlots.assign(bnk.templates.size(),
                    template_io::TemplateFile::BankSlotTemplate{});

                const size_t W = q.beats.front().size();
                std::vector<double> col;
                for (size_t sl = 0; sl < bnk.templates.size(); ++sl) {
                    const auto& mem = bnk.templates[sl].members;
                    if (mem.empty()) continue;
                    auto& st = outSlots[sl];
                    st.n_members = static_cast<uint32_t>(mem.size());
                    st.tmpl.assign(W, std::numeric_limits<double>::quiet_NaN());
                    st.tmpl_iqr.assign(W, std::numeric_limits<double>::quiet_NaN());
                    for (size_t c2 = 0; c2 < W; ++c2) {
                        col.clear();
                        for (const uint32_t m : mem) {
                            if (m >= q.beats.size()) continue;
                            const std::vector<double>& bt = q.beats[m];
                            if (c2 < bt.size() && !std::isnan(bt[c2]))
                                col.push_back(bt[c2]);
                        }
                        if (col.empty()) continue;
                        std::sort(col.begin(), col.end());
                        const size_t mid = col.size() / 2;
                        st.tmpl[c2] = (col.size() % 2 == 0)
                            ? 0.5 * (col[mid - 1] + col[mid]) : col[mid];
                        if (col.size() >= 2) {
                            double mu = 0.0;
                            for (double x : col) mu += x;
                            mu /= static_cast<double>(col.size());
                            double ss = 0.0;
                            for (double x : col) ss += (x - mu) * (x - mu);
                            st.tmpl_iqr[c2] =
                                std::sqrt(ss / static_cast<double>(col.size() - 1));
                        }
                    }
                }
            }

            dst.ecgTemplate = std::move(q.tmpl);
            dst.ecg_template_iqr = std::move(q.iqr);
        }
    }

    // ---- ONE LINE PER ANCHOR --------------------------------------------
    //
    // recovered = the landmark was not findable on that bin's reference BEAT,
    // so the target column came from its median template instead and the
    // alignment ran normally. That is the fallback working, not a failure.
    //
    // no-align = nothing moved at all, so that bin's template is a copy of the
    // R one. This is the only number that means something went wrong.
    std::cerr << "  [anchors] " << anchorName_bt(anchor)
        << ": " << tmpl.bins.size() << " bins";
    if (nRecovered > 0) std::cerr << ", " << nRecovered << " recovered via template";
    if (nNoAlign > 0)   std::cerr << ", " << nNoAlign << " NOT ALIGNED (copy of R)";
    std::cerr << "\n";
}

// SLOW merge: fills the squared/absval per-bin blocks and their SAECG
// averages into a TemplateFile produced by buildTemplatesAndBeatsFast.
// peakResults must carry the squared/absval R-peaks + preprocessed signals
// (run augment_ecg_ppg_pairs_sqabs first, or load them from wave_markings).
inline void mergeTemplatesSlow(const std::vector<output_binfile_data>& peakResults,
    template_io::TemplateFile& tmpl,
    std::vector<TemplateInfo>& info,
    const SignalRates& rates)
{
    using namespace template_generation_detail;

    AugmentTemplatesSlow(peakResults, info, rates);

    for (size_t i = 0; i < tmpl.bins.size(); ++i) {
        const TemplateInfo& bi = (i < info.size()) ? info[i] : TemplateInfo{};
        packBinSlow(tmpl.bins[i], bi);
    }
}

inline std::pair<template_io::TemplateFile, template_io::BeatsFile>
buildTemplatesAndBeatsFromPeakResults(const std::vector<output_binfile_data>& peakResults,
    const SignalRates& rates)
{
    FastTemplateBuild fast = buildTemplatesAndBeatsFast(peakResults, rates);
    mergeTemplatesSlow(peakResults, fast.tmpl, fast.info, rates);
    return { std::move(fast.tmpl), std::move(fast.beats) };
}