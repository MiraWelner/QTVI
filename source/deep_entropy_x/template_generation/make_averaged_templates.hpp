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
#include "template_morphology_grouping/bin_pipeline.hpp"
#include <iostream>
#include <chrono>
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
                    if (i >= ecg_res.kept_index[c].size()) continue;
                    ji.ecg_beats[c] = &ec[c]->kept_beats_raw[i];
                    ji.ecg_forward[c] = &ecg_res.kept_index[c][i];
                    ji.ecg_r_col[c] = (i < ec[c]->r_col_raw.size())
                        ? ec[c]->r_col_raw[i] : -1;
                }
                if (ppg_template_good && i < ppg_kept.size()
                    && i < ppg_kept_slices.size()) {
                    ji.ppg_beats = &ppg_kept[i];
                    ji.ppg_forward = &ppg_kept_slices[i];
                    ji.ppg_peak_col = ppg_peak_cols[i];
                    ji.ppg_phase1 = ppg_templates[i];
                }

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

                static const char* kKeys[4] = { "CH1", "CH2", "CH3", "PPG" };
                for (int c = 0; c < 4; ++c) {
                    bin_pipeline::ChannelOutput co;
                    co.bank = jbank::projectToChannel(info.joint.bank, cs, c);
                    co.assignment = info.joint.group_of_slice;
                    co.counts.n_spawns = info.joint.counts.n_spawns;
                    co.counts.n_merges = info.joint.counts.n_merges;
                    co.counts.n_cap_raises = info.joint.counts.n_cap_raises;
                    co.counts.n_unscorable = info.joint.counts.n_unscorable;
                    info.bank_by_channel[kKeys[c]] = std::move(co);
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
                // Section 4.6 bank. Without this the bank dies with
                // EcgChannelResult: create_ecg_templates builds it, and nothing
                // downstream ever sees it -- which renders in the viewer as
                // exactly one column per bin, indistinguishable from a record
                // with no ectopy.
                if (i < ecg_res.ch1.bank_out_raw.size())
                    info.bank_by_channel["CH1"] =
                        std::move(ecg_res.ch1.bank_out_raw[i]);
                info.ref_index_by_channel["CH1"] = (i < ecg_res.ch1.ref_index_raw.size()) ? ecg_res.ch1.ref_index_raw[i] : -1;
            }
            if (i < ecg_res.ch2.kept_beats_raw.size()) {
                info.kept_beats_by_channel["CH2"] = std::move(ecg_res.ch2.kept_beats_raw[i]);
                if (i < ecg_res.ch2.kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel["CH2"] = std::move(ecg_res.ch2.kept_rhythm_raw[i]);
                // Section 4.6 bank. Without this the bank dies with
                // EcgChannelResult: create_ecg_templates builds it, and nothing
                // downstream ever sees it -- which renders in the viewer as
                // exactly one column per bin, indistinguishable from a record
                // with no ectopy.
                if (i < ecg_res.ch2.bank_out_raw.size())
                    info.bank_by_channel["CH2"] =
                        std::move(ecg_res.ch2.bank_out_raw[i]);
                info.ref_index_by_channel["CH2"] = (i < ecg_res.ch2.ref_index_raw.size()) ? ecg_res.ch2.ref_index_raw[i] : -1;
            }
            if (i < ecg_res.ch3.kept_beats_raw.size()) {
                info.kept_beats_by_channel["CH3"] = std::move(ecg_res.ch3.kept_beats_raw[i]);
                if (i < ecg_res.ch3.kept_rhythm_raw.size())
                    info.kept_rhythm_by_channel["CH3"] = std::move(ecg_res.ch3.kept_rhythm_raw[i]);
                // Section 4.6 bank. Without this the bank dies with
                // EcgChannelResult: create_ecg_templates builds it, and nothing
                // downstream ever sees it -- which renders in the viewer as
                // exactly one column per bin, indistinguishable from a record
                // with no ectopy.
                if (i < ecg_res.ch3.bank_out_raw.size())
                    info.bank_by_channel["CH3"] =
                        std::move(ecg_res.ch3.bank_out_raw[i]);
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

                // ---- Section 4.6 bank, ON PPG -------------------------
                // The spec names two morphology thresholds, 0.85 for ECG and
                // 0.80 for PPG, so it asks for a bank on BOTH channel types.
                // Only the ECG banks were ever built: tbank::kMatchFloorPpg,
                // TemplateBank::matchFloorFor(is_ppg) and BeatFlags::
                // template_id_ppg all existed, and every caller passed
                // is_ppg = false, so the 0.80 floor was unreachable and a PPG
                // that changed morphology averaged into one template per bin --
                // the exact variance inflation 4.6 exists to remove, on the one
                // channel nobody split.
                //
                // NOTHING IN bin_pipeline::runChannel IS ECG-SPECIFIC. It takes
                // is_ppg and threads it to the floor, so seeding, best-match
                // assignment, spawn-at-floor, the cap, the merge, the
                // confirmed-member cap raise and the census are the SAME code
                // for both. That is why this is a call site rather than a
                // second implementation -- a parallel PPG bank would be a
                // second place for the seven rules to drift.
                bin_pipeline::ChannelInput pin;
                pin.beats = &ppg_kept[i];
                pin.width = ppg_kept[i].empty()
                    ? 0 : static_cast<int>(ppg_kept[i].front().size());
                // Column of the systolic peak on the shared axis. The pulse
                // slicer anchors every beat the same way, so one column serves
                // the whole bin -- the same contract r_col carries for ECG.
                pin.r_col = ppg_peak_cols[i];
                pin.is_ppg = true;                 // <- selects kMatchFloorPpg
                pin.bin_index = static_cast<uint64_t>(i);
                // NOT 0..2. Those index BeatFlags::template_id_ecg, and a PPG
                // assignment written into an ECG slot would read downstream as
                // an ECG template id. runChannel writes template_id_ppg from
                // is_ppg instead, and skips the ECG array when channel is out
                // of 0..2, so 3 is the value that keeps the two apart.
                pin.channel = 3;
                pin.max_templates_per_bin = 0;     // 0 => kDefaultMaxTemplatesPerBin
                // Seed slot 0 with the Phase 1 pulse template for this bin,
                // which is what ppg_templates[i] is. Same caveat as the ECG
                // seed and for the same reason: it is a column-wise median with
                // no rhythm test, so it is the sinus pulse by name rather than
                // by content -- hence not claiming otherwise here.
                pin.phase1_template = ppg_templates[i];
                pin.phase1_is_verified_sinus = false;
                // rr_after LEFT EMPTY, deliberately. The PVC filter is a timing
                // test on consecutive RR, and the PPG kept set is filtered
                // independently of the ECG kept set, so a beat at PPG index k
                // is not in general the beat at ECG index k. Supplying ECG RR
                // here would attach one beat's timing to another's morphology
                // and every premature verdict in the bin would be plausible and
                // wrong. Absent RR makes n_premature_members 0, so a PPG
                // template's presumed category rests on reproducibility alone,
                // which is honest about what is actually known.

                bin_pipeline::ChannelOutput pout = bin_pipeline::runChannel(pin);
                info.bank_by_channel["PPG"] = std::move(pout);

                info.kept_beats_by_channel["PPG"] = std::move(ppg_kept[i]);
            }
        }
        // else: info.ppgTemplate / ppg_template_iqr stay default-empty,
        // which the viewer already interprets as "no PPG for this bin".
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