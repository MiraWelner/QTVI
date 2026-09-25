/**
 * @file   make_averaged_templates.hpp
 * @brief  Orchestrate the full template generation pipeline.
 *         Port of GenerateTemplates.m
 *
 *         Under Patch B, PPG (and arterial) templates are built by the
 *         same [R_i - pad, R_{i+1} + pad] slicer as the ECG templates,
 *         driven by ch1.raw R-peaks. They come out R-anchored by
 *         construction, so the old find_foot -> AlignWaves -> NaN-strip
 *         PPG alignment pipeline is gone; PPG per-sample std rides through
 *         directly.
 *
 *         Rates for every channel arrive via a SignalRates struct
 *         (defined in TemplateTypes.hpp). A rate of 0 means the channel is
 *         absent from this dataset -- the slicer silently produces empty
 *         templates for those.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-08-20
 */
#pragma once

#include "template_structs.hpp"
#include "create_arterial_templates.hpp"
#include "create_ecg_templates.hpp"
#include "template_morphology_grouping/joint_bank.hpp"
#include "template_morphology_grouping/morphology_csv.hpp"
#include "template_morphology_grouping/nsvt_detect.hpp"
#include "noise_marking_gui/annotation_types.hpp"
 // noise_markings::Span / LoadResult / loadSpans / code_for_channel. The
 // operator class labels are an INPUT to template generation -- Section 4.6
 // partitions the bank by class BEFORE clustering -- so this pass reads the
 // noise-marking bin the GUI wrote.
#include "noise_marking_gui/user_annotation_handler.h"
#include <algorithm>   // std::nth_element for the per-bin median RR
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>      // numeric_limits<double>::quiet_NaN in the NSVT series
#include <memory>
#include <string>
#include <utility>


namespace morphology_writer {
    inline std::function<void()>& pending() {
        static std::function<void()> f;
        return f;
    }
    inline void runPending() {
        auto& f = pending();
        if (f) { f(); f = nullptr; }
    }
}

// NOT static. A template already has vague linkage, so `static` only gave
// every translation unit its own copy of an identical instantiation.
template <class PeakVec>
std::vector<uint8_t> sliceMarkCodes(
    const PeakVec& rPeaks, uint32_t n_slices,
    const std::vector<noise_markings::Span>& spans,
    uint64_t bin_index)
{
    std::vector<uint8_t> mark(n_slices, 0);
    if (rPeaks.empty() || spans.empty()) return mark;

    const uint8_t ppg = noise_markings::code_for_channel("PPG");

    for (const noise_markings::Span& s : spans) {
        if (s.channel_code == ppg) continue;
        // Scanned rather than testing codes 3 and 13 as literals: the whole
        // point of annotation_types is that renumbering an annotation needs one
        // edit, not two files.
        const annotation_types::AnnotationType* t = nullptr;
        for (const auto& row : annotation_types::noise_types)
            if (row.code == static_cast<int>(s.annotation_code)) { t = &row; break; }
        if (!t || t->paramEdit || t->invertEdit) continue;

        for (uint32_t k = 0; k < n_slices && k < rPeaks.size(); ++k) {
            const int64_t r = static_cast<int64_t>(rPeaks[k]);
            if (r < s.start_sample || r > s.end_sample) continue;
            if (mark[k] == 0) mark[k] = s.annotation_code;   // first in file order wins
        }
    }

    return mark;
}


inline vector<TemplateInfo> GenerateTemplatesFast(const vector<output_binfile_data>& wave_data,
    const SignalRates& rates,
    // <stem>_noise.bin from the noise-marking phase. DEFAULTED EMPTY so every
    // existing caller compiles and behaves exactly as before: no path means no
    // operator classes, every slice reads kUnlabeled, and the bank has one
    // partition -- which is the pre-partitioning behaviour.
    const std::string& noise_bin_path = {}) {
    size_t n = wave_data.size();
    // ---- phase timing: which part of the "fast" build is slow ----------
    auto _ms = [](std::chrono::steady_clock::time_point a,
        std::chrono::steady_clock::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count(); };

    // Check if any bin has PPG data.
    //
    // THE SIGNAL AND THE R-PEAKS, AND NOTHING ELSE. This also required
    // ppgMinAmps -- SegmentPPG's valley list -- which is a precondition Patch B
    // removed and this gate kept. The pulse templates are built by
    // build_pulse_template_pair_windowed, whose inputs are the signal, the
    // channel rate, ch1.raw and the ECG rate; neither it nor
    // extract_ppg_beats_and_align reads a valley. ppgMinAmps was what the old
    // find_foot -> AlignWaves aligner needed, and that aligner is gone.
    //
    // Requiring it meant any SegmentPPG throw -- whose catch(...) in
    // create_ecg_ppg_pairs clears the array -- skipped the entire pulse
    // channel, on a record whose PPG samples and R-peaks were sufficient to
    // build every template. The symptom is a channel that silently does not
    // exist: no pulse panels, no pulse columns, no message, and no response to
    // any pulse-QC setting because the QC never ran.
    //
    // ch1.raw is checked because the slicer needs at least two R-peaks to form
    // one [R_i - pad, R_i+1 + pad] window; it is the real precondition, and it
    // is the one CreatePulseTemplates itself tests per bin.
    bool has_ppg = false;
    for (size_t i = 0; i < n; ++i) {
        if (!wave_data[i].ppgSignal.empty() && wave_data[i].ch1.raw.size() >= 2) {
            has_ppg = true;
            break;
        }
    }


    // PPG templates (+ per-sample std, parallel shape). Under Patch B they
    // come out already R-anchored by construction (slice = [R_i-pad, R_i+1+pad]
    // at ppgRate, R_first at column pad*ppgRate), so the old find_foot ->
    // AlignWaves -> NaN-strip pipeline is unnecessary. We just hand the
    // templates through.
    vector<vector<double>> ppg_templates;
    vector<vector<double>> ppg_template_iqrs;
    vector<vector<vector<double>>> ppg_kept(n);
    vector<int> ppg_peak_cols(n, -1);
    vector<int> ppg_onset_cols(n, -1);
    // R-pair ordinals of the retained pulses, per bin. Hoisted out of the
    // ppg_res scope below because the joint bank needs them: without the
    // ordinal there is no way to say a pulse and a QRS are the same heartbeat.
    vector<vector<uint32_t>> ppg_kept_slices(n);
    vector<bool> template_good(n, false);

    const auto _ppg0 = std::chrono::steady_clock::now();
    if (has_ppg && rates.ppg > 0.0) {
        PPGTemplatesResult ppg_res = CreatePulseTemplates(wave_data, &output_binfile_data::ppgSignal, rates.ecg, rates.ppg);
        ppg_templates = std::move(ppg_res.templates);
        ppg_template_iqrs = std::move(ppg_res.iqrs);
        ppg_kept = std::move(ppg_res.kept);
        ppg_peak_cols = std::move(ppg_res.peakCol);
        ppg_onset_cols = std::move(ppg_res.footCol);
        ppg_kept_slices = std::move(ppg_res.keptSlices);

        for (size_t i = 0; i < n; ++i)
            for (double v : ppg_templates[i])
                if (!std::isnan(v)) { template_good[i] = true; break; }
    }

    // ECG templates -- FAST methods only (raw + unfiltered). The
    // squared/absval columns stay empty here; fill_channel copies those
    // empty vectors through harmlessly, and AugmentTemplatesSlow fills
    // them later.
    const auto _ppg1 = std::chrono::steady_clock::now();
    // NO SPAN PLUMBING HERE. An R-pair that is not a beat is one that straddles
    // a splice in the annealed bin, and process_channel_fast answers that from
    // bin.ecg_bin_indexs, which it already holds -- see
    // alignment::FragmentSeams. Nothing about annotations, coordinates or span
    // translation needs to reach the slicer.
    EcgTemplateResult ecg_res = CreateEcgTemplatesFast(wave_data, rates.ecg);
    const auto _ecg1 = std::chrono::steady_clock::now();
    std::fprintf(stderr, "[fast-phases] bins=%zu | PPG %8.1f  ECG(align+template) %8.1f ms\n", n, _ms(_ppg0, _ppg1), _ms(_ppg1, _ecg1));

    // Assemble TemplateInfo
    auto fill_channel = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_raw = src.ecgTemplates_raw[i];
        // Per-sample std for the raw method only -- the other three
        // methods are never displayed, so they don't have std computed.
        if (i < src.ecgTemplates_raw_iqr.size())
            dst.ecgTemplate_raw_iqr = src.ecgTemplates_raw_iqr[i];

        dst.ecgTemplate_squared = src.ecgTemplates_squared[i];
        dst.ecgTemplate_absval = src.ecgTemplates_absval[i];
        dst.ecgTemplate_unfiltered = src.ecgTemplates_unfiltered[i];

        dst.alignment_point_raw = std::isnan(src.ppg_alignment_point_raw[i]) ? 0.0 : src.ppg_alignment_point_raw[i];
        dst.alignment_point_squared = std::isnan(src.ppg_alignment_point_squared[i]) ? 0.0 : src.ppg_alignment_point_squared[i];
        dst.alignment_point_absval = std::isnan(src.ppg_alignment_point_absval[i]) ? 0.0 : src.ppg_alignment_point_absval[i];
        dst.alignment_point_unfiltered = std::isnan(src.ppg_alignment_point_unfiltered[i]) ? 0.0 : src.ppg_alignment_point_unfiltered[i];

        dst.r_col_raw = src.r_col_raw[i];
        dst.r_col_squared = src.r_col_squared[i];
        dst.r_col_absval = src.r_col_absval[i];
        dst.r_col_unfiltered = src.r_col_unfiltered[i];

        dst.n_beats_raw = (i < src.n_beats_raw.size()) ? src.n_beats_raw[i] : 0;
        };

    auto clear_channel = [](ChannelTemplates& dst) {
        dst.ecgTemplate_raw = {};
        dst.ecgTemplate_raw_iqr = {};
        dst.ecgTemplate_squared = {};
        dst.ecgTemplate_absval = {};
        dst.ecgTemplate_unfiltered = {};

        dst.alignment_point_raw = NaN;
        dst.alignment_point_squared = NaN;
        dst.alignment_point_absval = NaN;
        dst.alignment_point_unfiltered = NaN;

        dst.r_col_raw = -1;
        dst.r_col_squared = -1;
        dst.r_col_absval = -1;
        dst.r_col_unfiltered = -1;
        };

    vector<TemplateInfo> result(n);

    // ---- WHAT THE MORPHOLOGY ARCHIVE NEEDS, ACCUMULATED PER BIN ----------
    // Held here rather than rebuilt afterwards because the beat matrices are
    // MOVED into TemplateInfo inside the loop below: after the move, the only
    // place they exist is result[bin].kept_beats_by_channel, and the blocks
    // point there. Pointers, so nothing is copied -- a record's beats are the
    // largest thing in memory and duplicating them per channel to write a file
    // is what made the old beats writer look like a hang.
    static const char* kChanKeys[4] = { "CH1", "CH2", "CH3", "PPG" };
    std::array<morphology_csv::ChannelBlock, 4> mblocks;
    for (int c = 0; c < 4; ++c) {
        mblocks[c].channel = kChanKeys[c];
        mblocks[c].per_bin.assign(n, nullptr);
        mblocks[c].beats.assign(n, nullptr);
        mblocks[c].local_of_slice.assign(n, nullptr);
        mblocks[c].excluded_reason.assign(n, nullptr);
        mblocks[c].r_col.assign(n, -1);
    }
    // slice -> local row per channel per bin, kept alive for the writers. The
    // ChannelSet that produced it is a loop local, so the maps are copied out of
    // it here; they are one int32 per slice, which is nothing beside the beats.
    //
    // HEAP-OWNED, because the beats writer now runs after this function
    // returns. This is the ONE thing `blocks` points at that used to be a local
    // -- everything else points into `result`, whose addresses are stable -- so
    // the shared_ptr is the whole of what keeps the deferred closure valid. The
    // reference below means no other line in this function changes.
    auto local_of_slice_owned =
        std::make_shared<std::vector<std::array<std::vector<int32_t>, 4>>>(n);
    auto& local_of_slice = *local_of_slice_owned;

    // ---- OPERATOR CLASS LABELS, ONCE PER RECORD ------------------------
    //
    // Read here rather than per bin: the spans are record-wide and the file is
    // small, and re-opening it forty times would make a missing file forty log
    // lines instead of one.
    //
    // A MISSING FILE IS NOT AN ERROR. The record was never noise-marked, so
    // mark_code stays empty for every bin, jbank sees kUnlabeled throughout,
    // and the partition collapses to one -- exactly the behaviour before
    // partitioning existed. Said out loud, because "no operator classes" and
    // "the labels failed to load" produce the same partition and only one of
    // them is intended.
    noise_markings::LoadResult noiseSpans;
    if (!noise_bin_path.empty())
        noiseSpans = noise_markings::loadSpans(noise_bin_path);

    // One row per bin for <stem>_bins.csv: the 4.5 category census, the
    // partition's shape, and why beats left their group's average.
    std::vector<morphology_csv::BinRow> binRows;
    binRows.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        auto& info = result[i];

        const bool ppg_template_good = has_ppg && template_good[i];
        // ECG quality is independent of PPG availability. Previously this
        // function gated ECG fill on ppg_template_good, which cleared the
        // ECG templates for every bin of datasets without PPG (Bittium).
        // The top-of-file comment block documents the correct gate as
        // bad_segment, so we use that here.
        const bool ecg_good = (i < wave_data.size()) && !wave_data[i].bad_segment;

        // THE THREE ECG LEADS, BY INDEX. Every per-channel step below is the
        // same six statements three times over, and writing them out invited
        // the class of bug where ch3 keeps ch2's subscript.
        EcgChannelResult* const chRes[3] =
        { &ecg_res.ch1, &ecg_res.ch2, &ecg_res.ch3 };
        ChannelTemplates* const chDst[3] = { &info.ch1, &info.ch2, &info.ch3 };

        if (ecg_good) {
            const size_t nR = (i < wave_data.size()) ? wave_data[i].ch1.raw.size() : 0;
            if (nR >= 2) {
                jbank::BinBankInput ji;
                ji.n_slices = static_cast<uint32_t>(nR - 1);
                ji.bin_index = static_cast<uint64_t>(i);
                double medianRrSamples = 0.0;
                {
                    const auto& rp = wave_data[i].ch1.raw;
                    std::vector<double> rrs;
                    rrs.reserve(rp.size());
                    for (size_t sIdx = 0; sIdx + 1 < rp.size(); ++sIdx) {
                        const double d = (double)(rp[sIdx + 1] - rp[sIdx]);
                        if (d > 0.0) rrs.push_back(d);
                    }
                    if (!rrs.empty()) {
                        const size_t mid = rrs.size() / 2;
                        std::nth_element(rrs.begin(), rrs.begin() + mid, rrs.end());
                        medianRrSamples = rrs[mid];
                    }
                }
                if (medianRrSamples > 0.0) {
                    // The ECG window is already in ECG samples: the RR series
                    // above is measured on the ECG axis.
                    if (rates.ecg > 0.0 && rates.morph_halfwin_ecg_pct_rr > 0.0)
                        ji.ecg_corr_halfwin = static_cast<int>(
                            0.01 * rates.morph_halfwin_ecg_pct_rr
                            * medianRrSamples + 0.5);
                    // The pulse window is the same FRACTION OF THE SAME
                    // INTERVAL, re-expressed on the pulse axis -- via seconds,
                    // because the two channels are sampled at different rates
                    // and corr_halfwin is a sample count on its own channel.
                    if (rates.ecg > 0.0 && rates.ppg > 0.0
                        && rates.morph_halfwin_ppg_pct_rr > 0.0) {
                        const double rrSec = medianRrSamples / rates.ecg;
                        ji.ppg_corr_halfwin = static_cast<int>(
                            0.01 * rates.morph_halfwin_ppg_pct_rr
                            * rrSec * rates.ppg + 0.5);
                    }
                }

                for (int c = 0; c < 3; ++c) {
                    if (i >= chRes[c]->kept_beats_raw.size()) continue;
                    if (i >= ecg_res.kept_index[c].size()) continue;
                    ji.ecg_beats[c] = &chRes[c]->kept_beats_raw[i];
                    ji.ecg_forward[c] = &ecg_res.kept_index[c][i];
                    ji.ecg_r_col[c] = (i < chRes[c]->r_col_raw.size())
                        ? chRes[c]->r_col_raw[i] : -1;
                    if (i < chRes[c]->ecgTemplates_raw.size())
                        ji.ecg_phase1[c] = chRes[c]->ecgTemplates_raw[i];
                    if (i < chRes[c]->ecgTemplates_raw_iqr.size())
                        ji.ecg_phase1_spread[c] = chRes[c]->ecgTemplates_raw_iqr[i];
                }
                if (ppg_template_good && i < ppg_kept.size()
                    && i < ppg_kept_slices.size()) {
                    ji.ppg_beats = &ppg_kept[i];
                    ji.ppg_forward = &ppg_kept_slices[i];
                    ji.ppg_peak_col = ppg_peak_cols[i];
                    ji.ppg_phase1 = ppg_templates[i];
                    if (i < ppg_template_iqrs.size())
                        ji.ppg_phase1_spread = ppg_template_iqrs[i];
                }
                if (rates.ecg > 0.0) {
                    const auto& rp = wave_data[i].ch1.raw;
                    ji.rr_after_ms.assign(ji.n_slices, 0.0);
                    for (uint32_t sIdx = 0; sIdx + 1 < rp.size()
                        && sIdx < ji.n_slices; ++sIdx)
                        ji.rr_after_ms[sIdx] =
                        1000.0 * (double)(rp[sIdx + 1] - rp[sIdx]) / rates.ecg;
                }

                // ---- THE PARTITION KEY -------------------------------
                //
                // Empty when nothing was marked, which jbank reads as one
                // partition. See BeatGroup::partition.
                if (!noiseSpans.spans.empty())
                    ji.mark_code = sliceMarkCodes(wave_data[i].ch1.raw,
                        ji.n_slices, noiseSpans.spans, i);

                info.joint = jbank::buildBinBank(ji);
                info.joint_valid = true;

                // PROJECT IT INTO bank_by_channel, so the joint partition is
                // the ONLY partition. Nothing computes a per-channel bank
                // independently any more: these entries are channel views of
                // info.joint, template i of each being group i. Two independent
                // partitions of the same beats is the state that must not
                // exist, and this is what prevents it while the viewer and the
                // serializer still read the per-channel type.
                jbank::ChannelSet cs;
                for (int c = 0; c < 3; ++c)
                    if (ji.ecg_beats[c] && ji.ecg_forward[c])
                        jbank::setChannel(cs, c, *ji.ecg_beats[c],
                            *ji.ecg_forward[c], ji.n_slices, ji.ecg_r_col[c]);
                if (ji.ppg_beats && ji.ppg_forward)
                    jbank::setChannel(cs, jbank::kPpg, *ji.ppg_beats,
                        *ji.ppg_forward, ji.n_slices, ji.ppg_peak_col);

                for (int c = 0; c < 4; ++c) {
                    tbank::ChannelOutput co;
                    co.bank = jbank::projectToChannel(info.joint.bank, cs, c,
                        &info.joint.flags, &info.joint.rr_after_ms);

                    // BOTH IN SLICE SPACE, and the same length. flags and
                    // assignment used to be indexed by a channel's aligned row,
                    // which is why they could not be shared between channels;
                    // per slice they are one description of one set of
                    // heartbeats, and every block of the archive can key on it.
                    co.assignment = info.joint.group_of_slice;
                    co.flags = info.joint.flags;
                    co.pvc = info.joint.pvc;
                    co.counts.beats_detected = ji.n_slices;
                    co.counts.n_spawns = info.joint.counts.n_spawns;
                    co.counts.n_merges = info.joint.counts.n_merges;
                    co.counts.n_cap_raises = info.joint.counts.n_cap_raises;
                    co.counts.n_unscorable = info.joint.counts.n_unscorable;
                    info.bank_by_channel[kChanKeys[c]] = std::move(co);

                    // Copied out of the loop-local ChannelSet so the writers
                    // can still resolve slice -> row after this iteration ends.
                    local_of_slice[i][c] = cs[c].local_of_slice;
                    mblocks[c].local_of_slice[i] = &local_of_slice[i][c];
                    mblocks[c].excluded_reason[i] = &info.joint.excluded_reason;
                    mblocks[c].r_col[i] = cs[c].anchor_col;
                }

                // ---- the per-bin census row -----------------------------
                {
                    morphology_csv::BinRow row;
                    row.bin = static_cast<uint32_t>(i);
                    row.n_slices = ji.n_slices;
                    row.n_regular = info.joint.counts.n_regular;
                    row.n_ectopic = info.joint.counts.n_ectopic;
                    row.n_noise = info.joint.counts.n_noise;
                    row.n_groups = static_cast<uint32_t>(info.joint.bank.size());
                    row.n_spawns = info.joint.counts.n_spawns;
                    row.n_merges = info.joint.counts.n_merges;
                    row.n_cap_raises = info.joint.counts.n_cap_raises;
                    row.n_unscorable = info.joint.counts.n_unscorable;
                    row.ex_category = info.joint.clean.excluded_category;
                    row.ex_premature = info.joint.clean.excluded_premature;
                    row.ex_vote = info.joint.clean.excluded_vote;
                    row.ex_tukey = info.joint.clean.excluded_tukey;
                    row.n_kept = info.joint.clean.kept;
                    row.n_premature = info.joint.pvc.n_premature;
                    row.n_vote_only = info.joint.pvc.n_vote_only;
                    row.n_substituted = info.joint.subs.n_substituted;
                    row.n_sub_channel_blends = info.joint.subs.n_channel_blends;

                    for (uint32_t sIdx = 0; sIdx < ji.n_slices; ++sIdx) {
                        if (sIdx < info.joint.group_of_slice.size()
                            && info.joint.group_of_slice[sIdx] >= 0)
                            ++row.n_assigned;
                        if (sIdx < info.joint.excluded_reason.size()
                            && info.joint.excluded_reason[sIdx]
                            == static_cast<uint8_t>(jbank::ExcludeReason::NOT_A_MEMBER))
                            ++row.ex_not_member;
                        bool any = false;
                        for (int c = 0; c < 4 && !any; ++c)
                            if (sIdx < local_of_slice[i][c].size()
                                && local_of_slice[i][c][sIdx] >= 0) any = true;
                        if (any) ++row.n_became_beat;
                    }

                    // Operator-gated: before any marking this reports 0 and the
                    // writer prints none_confirmed, which is a different
                    // statement from monomorphic.
                    const jbank::PolymorphyVerdict pv =
                        jbank::polymorphyVerdict(info.joint.bank);
                    row.polymorphy_count = pv.count;
                    row.n_unconfirmed_groups =
                        static_cast<uint32_t>(pv.n_unconfirmed_groups);

                    row.seed_basis = (i < ecg_res.ch1.seed_basis_raw.size())
                        ? seed_pool::seedBasisName(
                            static_cast<seed_pool::SeedBasis>(
                                ecg_res.ch1.seed_basis_raw[i]))
                        : "";
                    binRows.push_back(row);
                }
            }

            // ---- PER-CHANNEL FILL, AND THE BEATS MOVED OUT ---------------
            //
            // bank_by_channel IS NOT WRITTEN HERE, and there is no longer
            // anything that could write it. It used to be assigned from a
            // per-channel bank built inside create_ecg_templates, AFTER the
            // joint projection above -- so the joint bank was built,
            // projected, and then silently overwritten by the very
            // per-channel partition it exists to abolish. Both that bank and
            // its last reader (the morphology writers, now on the projection)
            // are gone. The projection is the only writer.
            for (int c = 0; c < 3; ++c) {
                fill_channel(*chDst[c], *chRes[c], i);
                if (i >= chRes[c]->kept_beats_raw.size()) continue;
                info.kept_beats_by_channel[kChanKeys[c]] =
                    std::move(chRes[c]->kept_beats_raw[i]);
                if (i < chRes[c]->kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel[kChanKeys[c]] =
                    std::move(chRes[c]->kept_rhythm_raw[i]);
                info.ref_index_by_channel[kChanKeys[c]] =
                    (i < chRes[c]->ref_index_raw.size())
                    ? chRes[c]->ref_index_raw[i] : -1;
            }
        }
        else {
            for (int c = 0; c < 3; ++c) clear_channel(*chDst[c]);
        }

        if (ppg_template_good) {
            info.ppgTemplate = ppg_templates[i];
            info.ppg_template_iqr = ppg_template_iqrs[i];
            info.ppg_peak_col = ppg_peak_cols[i];
            info.ppg_onset_col = ppg_onset_cols[i];
            if (i < ppg_kept.size()) {
                info.ppg_n_beats = ppg_kept[i].size();

                info.kept_beats_by_channel["PPG"] = std::move(ppg_kept[i]);
            }
        }
        // else: info.ppgTemplate / ppg_template_iqr stay default-empty,
        // which the viewer already interprets as "no PPG for this bin".
    }

    //write templates.csv, beats.bin, templates.bin, and bins.csv.
    for (size_t i = 0; i < n; ++i) {
        for (int c = 0; c < 4; ++c) {
            const auto bit = result[i].bank_by_channel.find(kChanKeys[c]);
            if (bit != result[i].bank_by_channel.end())
                mblocks[c].per_bin[i] = &bit->second;
            const auto kit = result[i].kept_beats_by_channel.find(kChanKeys[c]);
            if (kit != result[i].kept_beats_by_channel.end())
                mblocks[c].beats[i] = &kit->second;
        }
    }
    const std::vector<morphology_csv::ChannelBlock> blocks(mblocks.begin(), mblocks.end());


    // =====================================================================
    // SECTION 4.6 NSVT: RECORD-LEVEL, ACROSS BIN BOUNDARIES
    // =====================================================================
    //
    // HERE and not inside the bin loop, because a run is not a per-bin object.
    // Three consecutive beats on one morphology can straddle a boundary, and
    // its beats carry a different GROUP INDEX on each side -- so scanning
    // per-bin indices both manufactures runs at every boundary (index 2 in bin
    // 7 and index 2 in bin 8 are unrelated morphologies) and misses the real
    // ones. Cross-bin global identity is the fix and it needs every bin's
    // groups in hand at once.
    std::vector<morphology_csv::NsvtRow> nsvtRows;
    uint32_t polyCandidates = 0;
    {
        // Per bin, the CH1 face of each group. Matching across bins is decided
        // on one channel; see nsvt::globalizeGroups.
        std::vector<std::vector<nsvt::GroupRef>> perBin(n);
        for (size_t i = 0; i < n; ++i) {
            if (!result[i].joint_valid) continue;
            const jbank::JointBank& jb = result[i].joint.bank;
            perBin[i].reserve(jb.groups.size());
            for (const jbank::BeatGroup& g : jb.groups) {
                nsvt::GroupRef gr;
                gr.tmpl = g.ch[jbank::kCh1].tmpl;
                gr.label_code = g.label_code;
                gr.subtype = g.subtype;
                gr.confirmed = g.confirmed();
                gr.n_members = static_cast<uint32_t>(g.memberCount());
                perBin[i].push_back(std::move(gr));
            }
        }
        const nsvt::JointGlobalMap gm = nsvt::globalizeGroups(perBin);

        // ---- the record-level per-beat series --------------------------
        // Beats in record order, which is bin order then slice order. RR comes
        // from the same per-slice series the prematurity filter used, so a run
        // rate and a prematurity verdict cannot disagree about an interval.
        nsvt::DetectInput di;
        for (size_t i = 0; i < n; ++i) {
            if (!result[i].joint_valid) continue;
            const auto& js = result[i].joint;
            const size_t ns = js.group_of_slice.size();
            for (size_t sIdx = 0; sIdx < ns; ++sIdx) {
                const int32_t grp = js.group_of_slice[sIdx];
                di.global_template.push_back(
                    (grp >= 0) ? gm.globalIdOf(i, grp) : -1);
                di.rr_after.push_back(
                    (sIdx < js.rr_after_ms.size()) ? js.rr_after_ms[sIdx]
                    : std::numeric_limits<double>::quiet_NaN());
                di.bin_of_beat.push_back(static_cast<uint32_t>(i));
            }
        }

        // NSVT IS OPERATOR-GATED: detectRuns only considers global templates
        // carrying a VENTRICULAR label, and labels come from marks. On an
        // unmarked record no run can be produced by any input, so an empty
        // result here is unfalsifiable rather than a clean verdict.
        const nsvt::GlobalMap flat = nsvt::asGlobalMap(gm);
        const std::vector<nsvt::NsvtRun> runs = nsvt::detectRuns(di, flat);
        polyCandidates = nsvt::countPolymorphicCandidates(di, flat);

        nsvtRows.reserve(runs.size());
        for (const nsvt::NsvtRun& r : runs) {
            morphology_csv::NsvtRow row;
            row.start_beat = r.start_beat;
            row.length = r.length;
            row.global_template = r.global_template;
            row.subtype = r.subtype;
            row.label_code = r.label_code;
            row.mean_cycle_ms = r.mean_cycle_ms;
            row.max_cycle_ms = r.max_cycle_ms;
            row.rate_bpm = r.rate_bpm;
            row.duration_ms = r.duration_ms;
            row.sustained = r.sustained ? 1u : 0u;
            row.crosses_bin = r.crosses_bin ? 1u : 0u;
            row.first_bin = r.first_bin;
            row.last_bin = r.last_bin;
            nsvtRows.push_back(row);
        }
        std::fprintf(stderr, "  [nsvt] morphologies=%zu beats=%zu runs=%zu polymorphic_candidates=%u\n",
            gm.morphologies.size(), di.global_template.size(),
            runs.size(), polyCandidates);
        std::fflush(stderr);
    }

    morphology_csv::writeBins(binRows);
    morphology_csv::writeNsvt(nsvtRows, polyCandidates);
    morphology_csv::writeTemplates(blocks);

    // DEFERRED, NOT WRITTEN. The ~3.8 s this used to cost lands on the
    // squared/absval pass instead -- see AugmentTemplatesSlow, which drains it.
    //
    // THE CLOSURE OUTLIVES THIS FUNCTION AND BORROWS FROM `result`. `blocks`
    // holds raw pointers into result[i].bank_by_channel and
    // kept_beats_by_channel -- std::map nodes, whose addresses survive the move
    // of `result` into the caller's variable, which is why this is sound today.
    // It is NOT sound if the caller drops (or copies and drops) that vector
    // before AugmentTemplatesSlow runs: the deferred write would then read
    // freed beats. local_of_slice_owned is a shared_ptr for exactly this
    // reason; the map interiors have no equivalent guard, only this contract.
    morphology_writer::pending() = [blocks, local_of_slice_owned] {
        morphology_csv::writeBeatsBin(blocks);
        };

    morphology_csv::writeTemplatesBin(blocks);

    return result;
}

// SLOW: fill the squared/absval ECG templates onto an existing
// vector<TemplateInfo> produced by GenerateTemplatesFast. Applies the same
// bad_segment gate as the fast pass, so squared/absval stay empty on bins
// the fast pass cleared.
inline void AugmentTemplatesSlow(const vector<output_binfile_data>& wave_data,
    vector<TemplateInfo>& templates,
    const SignalRates& rates)
{
    size_t n = wave_data.size();

    EcgTemplateResult ecg_res;
    init_channel_result(ecg_res.ch1, n);
    init_channel_result(ecg_res.ch2, n);
    init_channel_result(ecg_res.ch3, n);
    CreateEcgTemplatesSlow(wave_data, rates.ecg, ecg_res);

    auto fill_slow = [](ChannelTemplates& dst, const EcgChannelResult& src, size_t i) {
        dst.ecgTemplate_squared = src.ecgTemplates_squared[i];
        dst.ecgTemplate_absval = src.ecgTemplates_absval[i];
        dst.alignment_point_squared = std::isnan(src.ppg_alignment_point_squared[i]) ? 0.0 : src.ppg_alignment_point_squared[i];
        dst.alignment_point_absval = std::isnan(src.ppg_alignment_point_absval[i]) ? 0.0 : src.ppg_alignment_point_absval[i];
        dst.r_col_squared = src.r_col_squared[i];
        dst.r_col_absval = src.r_col_absval[i];
        };

    for (size_t i = 0; i < n && i < templates.size(); ++i) {
        const bool ecg_good = (i < wave_data.size()) && !wave_data[i].bad_segment;
        if (!ecg_good) continue;
        fill_slow(templates[i].ch1, ecg_res.ch1, i);
        fill_slow(templates[i].ch2, ecg_res.ch2, i);
        fill_slow(templates[i].ch3, ecg_res.ch3, i);
    }
    morphology_writer::runPending();//the beats.bin write takes ~3 so it is deferred with the abs and sqr processing
}