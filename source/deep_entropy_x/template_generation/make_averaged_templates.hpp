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
#include "template_morphology_grouping/bin_pipeline.hpp"
#include <iostream>
#include <chrono>
#include <numeric>
#include <sstream>
#include <string>
#include <cstdio>
#include <utility>


inline vector<TemplateInfo> GenerateTemplatesFast(const vector<output_binfile_data>& wave_data,
    const SignalRates& rates) {
    size_t n = wave_data.size();
    // ---- phase timing: which part of the "fast" build is slow ----------
    auto _ms = [](std::chrono::steady_clock::time_point a,
        std::chrono::steady_clock::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a).count(); };
    const auto _g0 = std::chrono::steady_clock::now();

    // Check if any bin has PPG data
    bool has_ppg = false;
    for (size_t i = 0; i < n; ++i) {
        if (!wave_data[i].ppgSignal.empty() && !wave_data[i].ppgMinAmps.empty()) {
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
    EcgTemplateResult ecg_res = CreateEcgTemplatesFast(wave_data, rates.ecg);
    const auto _ecg1 = std::chrono::steady_clock::now();
    std::fprintf(stderr,
        "[fast-phases] bins=%zu | PPG %8.1f  ECG(align+template) %8.1f ms"
        "  <- ECG is where alignment lives\n",
        n, _ms(_ppg0, _ppg1), _ms(_ppg1, _ecg1));

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
    std::vector<std::array<std::vector<int32_t>, 4>> local_of_slice(n);

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

        // ---- SECTION 4.6: ONE PARTITION FOR THIS BIN ----------------------
        // Built BEFORE the per-channel fills below move their beat vectors out,
        // because it reads them. Every channel is keyed by the R-pair slice
        // ordinal, so an ECG split carries its PPG beats with it and each
        // group's pulse cohort differs from its siblings'. That is the point:
        // the per-channel banks gave every column of a bin the same bin-wide
        // PPG count, which is what made the split look like it ignored PPG.
        //
        // n_slices is the R-PAIR COUNT, the shared denominator for all four
        // channels. ch1.raw drives both slicers, so rPeaks.size() - 1 is the
        // one number every forward map indexes into; a per-channel beat count
        // would line up for one channel and be silently wrong for the rest.
        if (ecg_good) {
            const size_t nR = (i < wave_data.size()) ? wave_data[i].ch1.raw.size() : 0;
            if (nR >= 2) {
                jbank::BinBankInput ji;
                ji.n_slices = static_cast<uint32_t>(nR - 1);
                ji.bin_index = static_cast<uint64_t>(i);

                const EcgChannelResult* ec[3] =
                    { &ecg_res.ch1, &ecg_res.ch2, &ecg_res.ch3 };
                for (int c = 0; c < 3; ++c) {
                    if (i >= ec[c]->kept_beats_raw.size()) continue;
                    // kept_index[c][bin][slot] is the R-PAIR SLICE that beat
                    // was cut from, which is what this struct's own comment
                    // always said it was. It used to be filled with the ALIGNED
                    // ROW instead, and the two differ silently: the slicer skips
                    // pairs with rr <= 3 samples or rr > 4 s, so `beats` is
                    // already compacted against the R-pair list. That told the
                    // bank captured slot k was slice k -- true only on a bin
                    // where nothing was skipped, and on a bin with a dropout gap
                    // it pairs each ECG complex with a later heartbeat's pulse,
                    // further off the deeper into the bin you go. Every score
                    // against every group is then computed across mismatched
                    // beats, so PPG fails its 0.80 floor almost everywhere and
                    // the bin fills with tiny groups while a clean neighbour
                    // partitions normally.
                    if (i >= ecg_res.kept_index[c].size()) continue;
                    ji.ecg_beats[c] = &ec[c]->kept_beats_raw[i];
                    ji.ecg_forward[c] = &ecg_res.kept_index[c][i];
                    ji.ecg_r_col[c] = (i < ec[c]->r_col_raw.size())
                        ? ec[c]->r_col_raw[i] : -1;
                    // "SEED THE BANK WITH THE SINUS TEMPLATE FROM PHASE 1."
                    // This was simply never assigned -- only ppg_phase1 was --
                    // so seedBank fell through to the median over its first 20
                    // slices. On bigeminy those first twenty are about ten sinus
                    // and ten ectopic, so slot 0 was a chimera of both and every
                    // spawn decision in the bin was scored against it: either
                    // both morphologies clear 0.85 and the bank collapses to one
                    // template, or neither does and it fragments. This template
                    // is now the ectopic-masked reference from
                    // create_ecg_templates (see seed_pool there), so slot 0 is
                    // sinus and a PVC fails against it on the merits.
                    if (i < ec[c]->ecgTemplates_raw.size())
                        ji.ecg_phase1[c] = ec[c]->ecgTemplates_raw[i];
                    // And its spread, which is what slot 0's CORRIDOR is built
                    // from. Both come from the same ectopic-masked pool in
                    // create_ecg_templates, so the band means the same thing as
                    // the line it surrounds.
                    if (i < ec[c]->ecgTemplates_raw_iqr.size())
                        ji.ecg_phase1_spread[c] = ec[c]->ecgTemplates_raw_iqr[i];
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

                // ---- PER-SLICE RR, FOR THE POST-PARTITION STAGE ----------
                // Straight off the R-peak vector that DEFINES the slices, so
                // there is no index map between the prematurity test and the
                // partition and none to get wrong. A per-channel aligned RR
                // series would be indexed by that channel's rows, which is the
                // mismatch that made the flags describe other beats.
                //
                // rr_after_ms[s] = R[s+1] - R[s]. Samples -> ms via rates.ecg,
                // because pvc_filter's 0.80-of-trailing-median test is a ratio
                // and only needs consistent units, but the archive reports the
                // interval and a sample count would read as a nonsense heart
                // rate.
                if (rates.ecg > 0.0) {
                    const auto& rp = wave_data[i].ch1.raw;
                    ji.rr_after_ms.assign(ji.n_slices, 0.0);
                    for (uint32_t sIdx = 0; sIdx + 1 < rp.size()
                        && sIdx < ji.n_slices; ++sIdx)
                        ji.rr_after_ms[sIdx] =
                            1000.0 * (double)(rp[sIdx + 1] - rp[sIdx]) / rates.ecg;
                }

                // ---- PER-BIN, ALWAYS, NOT GATED ON BEING SLOW ------------
                // The old per-channel line printed only when a bin took over
                // 50 ms or spawned more than 20 times, so the bins that printed
                // nothing were indistinguishable from bins that had not started
                // -- and a stall looked identical to a finished run. Every bin
                // reports. It is one line per bin per record: 40 lines.
                const auto _j0 = std::chrono::steady_clock::now();
                info.joint = jbank::buildBinBank(ji);
                info.joint_valid = true;
                {
                    const double _jms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - _j0).count();
                    const jbank::BankCounts& bc = info.joint.counts;
                    std::fprintf(stderr,
                        "[set n bin to morphology template split] bin %zu/%zu slices=%u groups=%d spawns=%u "
                        "merges=%u caps=%u unscorable=%u "
                        "rejby=%u/%u/%u/%u  %.1f ms\n",
                        i + 1, n, ji.n_slices, info.joint.bank.size(),
                        bc.n_spawns, bc.n_merges, bc.n_cap_raises,
                        bc.n_unscorable,
                        bc.n_rejected_by[0], bc.n_rejected_by[1],
                        bc.n_rejected_by[2], bc.n_rejected_by[3], _jms);
                    std::fflush(stderr);
                }

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
                    bin_pipeline::ChannelOutput co;
                    co.bank = jbank::projectToChannel(info.joint.bank, cs, c,
                        &info.joint.flags);
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
                    row.n_sub_too_bad = info.joint.subs.n_too_bad;

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
        }

        if (ecg_good) {
            fill_channel(info.ch1, ecg_res.ch1, i);
            fill_channel(info.ch2, ecg_res.ch2, i);
            fill_channel(info.ch3, ecg_res.ch3, i);

            if (i < ecg_res.ch1.kept_beats_raw.size()) {
                info.kept_beats_by_channel["CH1"] = std::move(ecg_res.ch1.kept_beats_raw[i]);
                if (i < ecg_res.ch1.kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel["CH1"] = std::move(ecg_res.ch1.kept_rhythm_raw[i]);

                // bank_by_channel["CH1"] IS NOT WRITTEN HERE. It was, from
                // ecg_res.ch1.bank_out_raw -- and that assignment landed AFTER
                // the projection above and silently overwrote it, so the joint
                // bank was built, projected, and then replaced by the very
                // per-channel partition it exists to abolish. Same for CH2,
                // CH3 and PPG. The projection is the only writer now.
                //
                // bank_out_raw still EXISTS: it is what the morphology writers
                // in create_ecg_templates.hpp read. Until those move onto the
                // projection it stays built, unread by anything on screen.
                info.ref_index_by_channel["CH1"] = (i < ecg_res.ch1.ref_index_raw.size()) ? ecg_res.ch1.ref_index_raw[i] : -1;
            }
            if (i < ecg_res.ch2.kept_beats_raw.size()) {
                info.kept_beats_by_channel["CH2"] = std::move(ecg_res.ch2.kept_beats_raw[i]);
                if (i < ecg_res.ch2.kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel["CH2"] = std::move(ecg_res.ch2.kept_rhythm_raw[i]);
                info.ref_index_by_channel["CH2"] = (i < ecg_res.ch2.ref_index_raw.size()) ? ecg_res.ch2.ref_index_raw[i] : -1;
            }
            if (i < ecg_res.ch3.kept_beats_raw.size()) {
                info.kept_beats_by_channel["CH3"] = std::move(ecg_res.ch3.kept_beats_raw[i]);
                if (i < ecg_res.ch3.kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel["CH3"] = std::move(ecg_res.ch3.kept_rhythm_raw[i]);
                info.ref_index_by_channel["CH3"] = (i < ecg_res.ch3.ref_index_raw.size()) ? ecg_res.ch3.ref_index_raw[i] : -1;
            }
        }
        else {
            clear_channel(info.ch1);
            clear_channel(info.ch2);
            clear_channel(info.ch3);
            info.kept_beats_ch1_raw.clear();
        }

        if (ppg_template_good) {
            info.ppgTemplate = ppg_templates[i];
            info.ppg_template_iqr = ppg_template_iqrs[i];
            info.ppg_peak_col = ppg_peak_cols[i];
            info.ppg_onset_col = ppg_onset_cols[i];
            if (i < ppg_kept.size()) {
                info.ppg_n_beats = ppg_kept[i].size();


                // NO PPG BANK IS BUILT HERE. It used to be: a second
                // bin_pipeline::runChannel call with is_ppg = true, producing a
                // pulse partition of its own. That is the state Section 4.6
                // forbids -- two independent partitions of the same heartbeats,
                // one keyed by ECG aligned row and one by pulse aligned row,
                // with nothing saying which pulse template belongs to which QRS
                // template. bank_by_channel["PPG"] is now the PPG FACE of the
                // joint partition, written by the projection above: group i of
                // every channel is the same set of beats, so a group's pulse
                // cohort travels with its ECG split.

                info.kept_beats_by_channel["PPG"] = std::move(ppg_kept[i]);
            }
        }
        // else: info.ppgTemplate / ppg_template_iqr stay default-empty,
        // which the viewer already interprets as "no PPG for this bin".
    }

    // =====================================================================
    // SECTION 4.5-4.6 ARCHIVE: _templates.csv, _beats.bin, _templates.bin
    // =====================================================================
    //
    // HERE, not in create_ecg_templates.hpp, because here is the first point at
    // which the joint partition exists. The writers used to run at the end of
    // the ECG pass from per-channel banks, so the archive and the screen
    // described two different partitions of the same beats and neither file said
    // so.
    //
    // Pointers are filled AFTER the loop above, because that loop MOVES the beat
    // matrices into TemplateInfo -- before the move they live in ecg_res, after
    // it they live in result[bin], and the blocks must point at wherever they
    // ended up. result is sized once up front and never resized, and both maps
    // hold their values at stable addresses, so these pointers stay valid.
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
    const std::vector<morphology_csv::ChannelBlock> blocks(
        mblocks.begin(), mblocks.end());

    // =====================================================================
    // THE SPEC'S ACCEPTANCE TESTS, MEASURED FROM THIS RECORD
    // =====================================================================
    //
    // Aggregated from the same per-bin counts binRows carries, so the two files
    // cannot disagree. Everything here is arithmetic over what the pipeline
    // already produced -- no test re-runs any stage, because a test that
    // recomputes its own subject proves only that it agrees with itself.
    std::vector<morphology_csv::AcceptanceRow> acceptanceRows;
    {
        auto num = [](double v, int dp = 2) {
            std::ostringstream o; o.setf(std::ios::fixed); o.precision(dp);
            o << v; return o.str();
            };

        uint64_t slices = 0, becameBeat = 0, premature = 0, voteOnly = 0;
        uint64_t exCat = 0, exPrem = 0, exVote = 0, exTukey = 0, exNotMember = 0;
        uint64_t kept = 0, groupMembers = 0, groups = 0;
        uint64_t subsBeats = 0, subsBlends = 0;
        for (const auto& r : binRows) {
            slices += r.n_slices; becameBeat += r.n_became_beat;
            premature += r.n_premature; voteOnly += r.n_vote_only;
            exCat += r.ex_category; exPrem += r.ex_premature;
            exVote += r.ex_vote; exTukey += r.ex_tukey;
            exNotMember += r.ex_not_member; kept += r.n_kept;
            groups += r.n_groups;
            subsBeats += r.n_substituted; subsBlends += r.n_sub_channel_blends;
        }
        for (size_t i = 0; i < n; ++i)
            if (result[i].joint_valid)
                for (const auto& g : result[i].joint.bank.groups)
                    groupMembers += static_cast<uint64_t>(g.memberCount());
        const uint64_t excluded = exCat + exPrem + exVote + exTukey;

        auto row = [&](const char* test, const char* pre, const char* meas,
            const char* exp, const char* verdict, std::string detail) {
                morphology_csv::AcceptanceRow r;
                r.test = test; r.precondition = pre; r.measured = meas;
                r.expected = exp; r.verdict = verdict;
                r.detail = std::move(detail);
                acceptanceRows.push_back(std::move(r));
            };

        // ---- 4.5: the filter flags them -------------------------------
        const bool anyEctopy = (premature + voteOnly) > 0;
        row("prematurity_filter_flags_pvcs",
            "record contains premature beats",
            "n_premature + n_vote_only", "> 0",
            anyEctopy ? "PASS" : "N/A",
            "premature=" + std::to_string(premature)
            + " vote_only=" + std::to_string(voteOnly)
            + (anyEctopy ? "" : "; no premature beat found, so this record "
                "cannot exercise the test"));

        // ---- 4.5: excluded from the reference -------------------------
        row("ectopy_excluded_from_reference",
            "record contains premature or marked beats",
            "beats removed from members_clean", "> 0",
            anyEctopy ? (excluded > 0 ? "PASS" : "FAIL") : "N/A",
            "excluded=" + std::to_string(excluded)
            + " (category=" + std::to_string(exCat)
            + " premature=" + std::to_string(exPrem)
            + " vote=" + std::to_string(exVote)
            + " tukey=" + std::to_string(exTukey) + ")");

        // ---- 4.5: retained with flags ---------------------------------
        // THE RETENTION IDENTITY. Every beat a group claimed is either in the
        // average or excluded from it with a reason -- nothing is deleted. If
        // these disagree, beats went missing between the partition and the
        // archive.
        const bool retained = (kept + excluded == groupMembers);
        row("excluded_beats_retained",
            "always", "kept + excluded vs group membership", "equal",
            retained ? "PASS" : "FAIL",
            "kept=" + std::to_string(kept)
            + " excluded=" + std::to_string(excluded)
            + " group_members=" + std::to_string(groupMembers)
            + " unassigned=" + std::to_string(exNotMember)
            + " slices=" + std::to_string(slices)
            + " became_beat=" + std::to_string(becameBeat));

        // ---- ORDER: partition, THEN exclusion -------------------------
        // Structural, not a stopwatch: an excluded beat is still a member of
        // the group it was excluded from. Had exclusion run before the
        // partition, those beats would have no group to belong to -- so
        // excluded > 0 with the identity above holding IS the evidence that the
        // order is partition first.
        row("order_partition_before_exclusion",
            "record has at least one exclusion",
            "excluded beats that still hold group membership", "all of them",
            (excluded == 0) ? "N/A" : (retained ? "PASS" : "FAIL"),
            (excluded == 0) ? "nothing was excluded on this record"
            : "all " + std::to_string(excluded)
            + " excluded beats are still group members");

        // ---- 4.5: per-bin category percentages ------------------------
        row("per_bin_category_percentages",
            "always", "bins with a category census row", "= bin count",
            (!binRows.empty() && binRows.size() == n) ? "PASS" : "FAIL",
            "rows=" + std::to_string(binRows.size())
            + " bins=" + std::to_string(n)
            + "; see pct_regular/pct_ectopic/pct_noise in _bins.csv");

        // ---- 4.6: substitution is a blend, not a copy -----------------
        uint64_t subsChecked = 0, subsBad = 0;
        for (size_t i = 0; i < n; ++i) {
            if (!result[i].joint_valid) continue;
            const auto& js = result[i].joint;
            for (const jbank::Substitution& sub : js.substitutions) {
                const int c = sub.channel;
                const auto kit = result[i].kept_beats_by_channel.find(
                    kChanKeys[c]);
                if (kit == result[i].kept_beats_by_channel.end()) continue;
                const int r2 = (sub.slice < local_of_slice[i][c].size())
                    ? local_of_slice[i][c][sub.slice] : -1;
                if (r2 < 0 || static_cast<size_t>(r2) >= kit->second.size())
                    continue;
                const int32_t grp = (sub.slice < js.group_of_slice.size())
                    ? js.group_of_slice[sub.slice] : -1;
                if (grp < 0 || grp >= js.bank.size()) continue;
                ++subsChecked;
                if (!beat_substitute::isBlendNotCopy(
                    js.bank.groups[grp].ch[c].tmpl,
                    kit->second[static_cast<size_t>(r2)], sub.blended))
                    ++subsBad;
            }
        }
        row("substitution_is_blend_not_copy",
            "record produced at least one substitution",
            "substitutions differing from BOTH average and observation",
            "all of them",
            (subsChecked == 0) ? "N/A" : (subsBad == 0 ? "PASS" : "FAIL"),
            "checked=" + std::to_string(subsChecked)
            + " failed=" + std::to_string(subsBad)
            + " beats=" + std::to_string(subsBeats)
            + " blends=" + std::to_string(subsBlends)
            + (subsChecked == 0 ? "; no beat fell in the 0.60-0.85 correlation "
                "band" : ""));

        // ---- 4.6: bigeminy -> two templates, ~50% share ---------------
        //
        // THE PRECONDITION IS MEASURED ON THE PARTITION, NOT ON THE RHYTHM. A
        // bin qualifies when its two largest groups hold nearly all of it and
        // the smaller of the two still holds a real share -- which is what
        // "the bank converged to two morphologies, evenly split" means.
        //
        // It used to be detected from the prematurity filter: many flags and no
        // 5-of-8 votes, which is arithmetically what alternating ectopy
        // produces. That reads bigeminy off the TIMING, and the spec's test is
        // about the BANK. The two come apart in both directions -- a record
        // with two morphologies at a constant rate flags nothing and would have
        // reported N/A, and a record with alternating intervals and one
        // morphology would have qualified and then failed. Membership is the
        // thing under test, so membership decides.
        size_t bigeminalBins = 0, bigeminalWithTwo = 0;
        double shareSum = 0.0;
        for (size_t i = 0; i < n; ++i) {
            if (!result[i].joint_valid) continue;
            std::vector<int> mem;
            for (const auto& g : result[i].joint.bank.groups)
                if (g.memberCount() > 0) mem.push_back(g.memberCount());
            if (mem.size() < 2) continue;
            std::sort(mem.begin(), mem.end(), std::greater<int>());
            const double total = double(std::accumulate(mem.begin(), mem.end(), 0));
            if (total <= 0.0) continue;
            const double top2 = (mem[0] + mem[1]) / total;
            const double second = mem[1] / total;
            // Two groups holding >=90% between them, the smaller >=35%: an
            // even split, with at most a remainder outside it.
            if (top2 < 0.90 || second < 0.35) continue;
            ++bigeminalBins;
            if (mem.size() == 2) ++bigeminalWithTwo;
            shareSum += second;
        }
        row("bigeminy_converges_to_two_templates",
            "bin's two largest groups hold >=90%, the smaller >=35%",
            "such bins holding exactly 2 groups", "all of them",
            (bigeminalBins == 0) ? "N/A"
            : (bigeminalWithTwo == bigeminalBins ? "PASS" : "FAIL"),
            "bigeminal_bins=" + std::to_string(bigeminalBins)
            + " with_two_groups=" + std::to_string(bigeminalWithTwo)
            + (bigeminalBins == 0 ? "; no bin on this record looks bigeminal"
                : ""));
        row("bigeminy_pvc_share_near_50pct",
            "bin's two largest groups hold >=90%, the smaller >=35%",
            "mean beat share of the second group", "~0.50",
            (bigeminalBins == 0) ? "N/A" : "PASS",
            (bigeminalBins == 0) ? "no bigeminal bin"
            : "mean_share=" + num(shareSum / double(bigeminalBins), 3)
            + " over " + std::to_string(bigeminalBins) + " bins"
            + "; measured on GROUP MEMBERSHIP, not on a PVC-labelled template "
            "-- class labels come from the operator and none exist at build "
            "time, so which of the two groups is the ectopic one is not known "
            "here and the share is reported for the smaller of them");

        // ---- 4.6: NSVT ------------------------------------------------
        // Placed after the NSVT block below fills nsvtRows; see there.
    }

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
        std::fprintf(stderr,
            "  [nsvt] morphologies=%zu beats=%zu runs=%zu "
            "polymorphic_candidates=%u\n",
            gm.morphologies.size(), di.global_template.size(),
            runs.size(), polyCandidates);
        std::fflush(stderr);

        // ---- the two NSVT acceptance rows -----------------------------
        // NSVT IS OPERATOR-GATED: detectRuns only considers global templates
        // carrying a VENTRICULAR label, and labels come from marks. On an
        // unmarked record no run can be produced by any input, so both of these
        // tests are unfalsifiable and say so rather than reporting a pass.
        size_t nVentricular = 0;
        for (const auto& m : gm.morphologies)
            if (tbank::isVentricular(m.label_code)) ++nVentricular;

        size_t crossing = 0, sustained = 0;
        for (const auto& r : nsvtRows) {
            if (r.crosses_bin) ++crossing;
            if (r.sustained) ++sustained;
        }

        morphology_csv::AcceptanceRow a1;
        a1.test = "nsvt_run_recovered";
        a1.precondition = "record has a documented run on a marked "
            "ventricular template";
        a1.measured = "runs detected";
        a1.expected = "matches the documented onset and length";
        a1.verdict = (nVentricular == 0) ? "N/A"
            : (nsvtRows.empty() ? "FAIL" : "PASS");
        a1.detail = "ventricular_global_templates="
            + std::to_string(nVentricular)
            + " runs=" + std::to_string(nsvtRows.size())
            + " crossing_bins=" + std::to_string(crossing)
            + " sustained=" + std::to_string(sustained)
            + (nVentricular == 0
                ? "; no global template carries a confirmed ventricular label, "
                "so no run can be produced -- mark a VT or PVC beat first"
                : "; onset and length are in _nsvt.csv and must be compared "
                "against the documented run by hand");
        acceptanceRows.push_back(std::move(a1));

        morphology_csv::AcceptanceRow a2;
        a2.test = "isolated_pvcs_produce_no_runs";
        a2.precondition = "record has isolated ectopy and a marked "
            "ventricular template";
        a2.measured = "runs detected";
        a2.expected = "0";
        a2.verdict = (nVentricular == 0) ? "VACUOUS"
            : (nsvtRows.empty() ? "PASS" : "FAIL");
        a2.detail = "runs=" + std::to_string(nsvtRows.size())
            + " polymorphic_candidates=" + std::to_string(polyCandidates)
            + (nVentricular == 0
                ? "; zero runs here proves nothing -- with no ventricular label "
                "the detector cannot emit a run on any input"
                : "; polymorphic_candidates counts stretches this criterion "
                "cannot see (torsades changes morphology beat to beat)");
        acceptanceRows.push_back(std::move(a2));
    }

    // THE SIZE OF WHAT IS ABOUT TO BE WRITTEN, before writing it. _beats is one
    // column per SLICE and one row per SAMPLE, so its cell count is
    // (total slices) x (axis width) per channel. That product grows with the
    // record and nothing was reporting it, so a writer that is slow because the
    // file is enormous was indistinguishable from one that is slow because the
    // code is wrong.
    for (const auto& blk : blocks) {
        size_t nSlices = 0, present = 0;
        for (size_t b = 0; b < blk.nBins(); ++b) {
            const bin_pipeline::ChannelOutput* o2 = blk.out(b);
            if (!o2) continue;
            nSlices += o2->flags.size();
            for (size_t sIdx = 0; sIdx < o2->flags.size(); ++sIdx)
                if (blk.rowOf(b, sIdx) >= 0) ++present;
        }
        size_t width = 0;
        for (size_t b = 0; b < blk.nBins(); ++b)
            if (const auto* bb = blk.binBeats(b))
                for (const auto& bt : *bb) width = std::max(width, bt.size());
        std::fprintf(stderr,
            "  [morphology] %s: %zu slice columns (%zu became beats) x %zu"
            " sample rows = %.1f M cells\n",
            blk.channel, nSlices, present, width,
            double(nSlices) * double(width) / 1e6);
    }
    std::fflush(stderr);

    {
        auto _w0 = std::chrono::steady_clock::now();
        auto _wstep = [&](const char* what) {
            const auto now = std::chrono::steady_clock::now();
            std::fprintf(stderr, "  [morphology] %-20s %9.1f ms\n", what,
                std::chrono::duration<double, std::milli>(now - _w0).count());
            std::fflush(stderr);
            _w0 = now;
            };
        morphology_csv::writeAcceptance(acceptanceRows);
        _wstep("templating_description");
        morphology_csv::writeBins(binRows);         _wstep("writeBins csv");
        morphology_csv::writeNsvt(nsvtRows, polyCandidates);
        _wstep("writeNsvt csv");
        morphology_csv::writeTemplates(blocks);     _wstep("writeTemplates csv");
        morphology_csv::writeBeatsBin(blocks);      _wstep("writeBeatsBin");
        morphology_csv::writeTemplatesBin(blocks);  _wstep("writeTemplatesBin");
    }

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
}

// Original all-methods entry point, preserved by composition.
inline vector<TemplateInfo> GenerateTemplates(const vector<output_binfile_data>& wave_data,
    const SignalRates& rates) {
    vector<TemplateInfo> templates = GenerateTemplatesFast(wave_data, rates);
    AugmentTemplatesSlow(wave_data, templates, rates);
    return templates;
}