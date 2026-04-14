#pragma once

#include "ppg_features.hpp"
#include <fstream>
#include <stdexcept>
#include <cstring>

namespace ppg {

    // ═══════════════════════════════════════════════════════════════════════════
    // Annealed Segments .bin
    // ─────────────────────────────────────────────────────────────────────────
    // Written by 3_anneal_segments.cpp:
    //   Header:  [uint64 nSegments] [double ppgSR] [double ecgSR] [double epochSec]
    //   Per seg: ppg_bin_indexs  → [uint64 nPairs] then nPairs × [uint64, uint64]
    //            ecg_bin_indexs  → same layout
    //            ppg             → [uint64 len] [len × double]
    //            ecg1            → same
    //            ecg2            → same
    //            ecg3            → same
    //            sleep_stages    → same
    // ═══════════════════════════════════════════════════════════════════════════

    inline std::vector<AnnealedSegment> read_annealed_bin(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open())
            throw std::runtime_error("Cannot open annealed bin: " + path);

        f.seekg(0, std::ios::end);
        uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        auto r64 = [&]() -> uint64_t {
            uint64_t v = 0; f.read(reinterpret_cast<char*>(&v), 8);
            if (!f.good()) throw std::runtime_error("Unexpected EOF in annealed at pos " + std::to_string(f.tellg()));
            return v;
            };
        auto rd = [&]() -> double {
            double v; f.read(reinterpret_cast<char*>(&v), 8); return v;
            };
        auto readVec = [&]() -> std::vector<double> {
            uint64_t n = r64();
            std::vector<double> v(n);
            if (n > 0) f.read(reinterpret_cast<char*>(v.data()), n * 8);
            return v;
            };
        auto readPairs = [&]() -> std::vector<std::vector<int>> {
            uint64_t n = r64();
            std::vector<std::vector<int>> pairs(n, std::vector<int>(2));
            for (uint64_t i = 0; i < n; i++) {
                pairs[i][0] = static_cast<int>(r64());
                pairs[i][1] = static_cast<int>(r64());
            }
            return pairs;
            };

        uint64_t nSegments = r64();
        double ppgSR = 0, ecgSR = 0, epochSec = 0;
        if (nSegments > 0) {
            ppgSR = rd();
            ecgSR = rd();
            epochSec = rd();
        }

        std::vector<AnnealedSegment> segs(nSegments);
        for (uint64_t i = 0; i < nSegments; i++) {
            segs[i].ppgSampleRate = ppgSR;
            segs[i].ecgSampleRate = ecgSR;
            segs[i].ppg_bin_indexs = readPairs();

            // Skip ecg_bin_indexs (not used in feature generation)
            {
                uint64_t n = r64();
                f.seekg(n * 16, std::ios::cur);
            }

            segs[i].po = readVec();  // PPG signal

            // Skip ecg channels (not used in feature generation)
            for (int ch = 0; ch < 3; ch++) {
                uint64_t n = r64();
                f.seekg(n * 8, std::ios::cur);
            }

            segs[i].sleep_stages = readVec();
        }
        return segs;
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Wave Data .bin
    // ─────────────────────────────────────────────────────────────────────────
    // Layout:
    //   Header:  [uint64 nBins]
    //   Per bin:
    //     -- 9 ECG R-peak index arrays (ch1 raw/sq/abs, ch2 raw/sq/abs, ch3 raw/sq/abs)
    //        each: [uint64 len] then len × [int64]
    //     -- 2 PPG index arrays (maxAmps, minAmps)
    //        each: [uint64 len] then len × [double]
    //     -- 4 raw signal arrays (ppg, ecg1, ecg2, ecg3)
    //        each: [uint64 len] then len × [double]
    //     -- 6 preprocessed signal arrays (ch1 sq/abs, ch2 sq/abs, ch3 sq/abs)
    //        each: [uint64 len] then len × [double]
    //     -- 9 noise flags (1 byte each, uint8)
    //     -- pairs array: [uint64 nPairs] then nPairs × [int64 ppgIdx, int64 ecgIdx]
    //     -- ppg_bin_indexs: [uint64 n] then n × [uint64, uint64]
    //     -- ecg_bin_indexs: [uint64 n] then n × [uint64, uint64]
    // ═══════════════════════════════════════════════════════════════════════════

    inline std::vector<WaveData> read_wave_data_bin(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open())
            throw std::runtime_error("Cannot open wave_data bin: " + path);

        // File size for bounds checking
        f.seekg(0, std::ios::end);
        uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        auto r64 = [&]() -> uint64_t {
            uint64_t v = 0;
            f.read(reinterpret_cast<char*>(&v), 8);
            if (!f.good()) throw std::runtime_error("Unexpected EOF in wave_data at pos " + std::to_string(f.tellg()));
            return v;
            };

        // Read an array of uint64 indices, converting from 1-based to 0-based
        auto readIdxArr = [&](const char* label) -> std::vector<int> {
            uint64_t n = r64();
            if (n > fileSize / 8) throw std::runtime_error(
                std::string("Implausible array size ") + std::to_string(n) + " for " + label
                + " at pos " + std::to_string(f.tellg()));
            std::vector<int> v(n);
            for (uint64_t i = 0; i < n; i++) {
                uint64_t val;
                f.read(reinterpret_cast<char*>(&val), 8);
                v[i] = static_cast<int>(val > 0 ? val - 1 : 0);  // 1-based → 0-based
            }
            return v;
            };

        // Read an array of doubles
        auto readDblArr = [&](const char* label) -> std::vector<double> {
            uint64_t n = r64();
            if (n > fileSize / 8) throw std::runtime_error(
                std::string("Implausible array size ") + std::to_string(n) + " for " + label
                + " at pos " + std::to_string(f.tellg()));
            std::vector<double> v(n);
            if (n > 0) f.read(reinterpret_cast<char*>(v.data()), n * 8);
            return v;
            };

        // Skip a double array without storing
        auto skipDblArr = [&](const char* label) {
            uint64_t n = r64();
            if (n > fileSize / 8) throw std::runtime_error(
                std::string("Implausible skip size ") + std::to_string(n) + " for " + label
                + " at pos " + std::to_string(f.tellg()));
            if (n > 0) f.seekg(n * 8, std::ios::cur);
            };

        // Read pairs of uint64
        auto readPairArr = [&]() -> std::vector<std::vector<int>> {
            uint64_t n = r64();
            std::vector<std::vector<int>> pairs(n, std::vector<int>(2));
            for (uint64_t i = 0; i < n; i++) {
                pairs[i][0] = static_cast<int>(r64());
                pairs[i][1] = static_cast<int>(r64());
            }
            return pairs;
            };

        uint64_t nBins = r64();
        if (nBins > 100000) throw std::runtime_error(
            "Implausible nBins=" + std::to_string(nBins) + " — wrong file format?");
        std::vector<WaveData> segs(nBins);

        for (uint64_t b = 0; b < nBins; b++) {
            auto& seg = segs[b];
            if (b < 2 || b == nBins - 1 || b % 100 == 0)
            {
            } // progress (removed verbose output)

       // ── 9 ECG R-peak index arrays ──
            seg.ch1_raw_idx = readIdxArr("ch1_raw_idx");
            seg.ch1_sq_idx = readIdxArr("ch1_sq_idx");
            seg.ch1_abs_idx = readIdxArr("ch1_abs_idx");
            seg.ch2_raw_idx = readIdxArr("ch2_raw_idx");
            seg.ch2_sq_idx = readIdxArr("ch2_sq_idx");
            seg.ch2_abs_idx = readIdxArr("ch2_abs_idx");
            seg.ch3_raw_idx = readIdxArr("ch3_raw_idx");
            seg.ch3_sq_idx = readIdxArr("ch3_sq_idx");
            seg.ch3_abs_idx = readIdxArr("ch3_abs_idx");

            // ── 2 PPG index arrays (not used in feature gen — skip) ──
            { uint64_t n = r64(); f.seekg(n * 8, std::ios::cur); }  // ppgMaxAmps
            { uint64_t n = r64(); f.seekg(n * 8, std::ios::cur); }  // ppgMinAmps

            // ── 4 raw signal arrays (not needed for feature gen — skip) ──
            skipDblArr("ppg_raw");
            skipDblArr("ecg1_raw");
            skipDblArr("ecg2_raw");
            skipDblArr("ecg3_raw");

            // ── 6 preprocessed signal arrays (skip) ──
            skipDblArr("ch1_sq");
            skipDblArr("ch1_abs");
            skipDblArr("ch2_sq");
            skipDblArr("ch2_abs");
            skipDblArr("ch3_sq");
            skipDblArr("ch3_abs");

            // ── 9 noise flags (1 byte each) ──
            f.read(reinterpret_cast<char*>(seg.noise_flags), 9);

            // ── pairs array (int64, 1-based → 0-based, -1 stays as -1) ──
            {
                uint64_t nPairs = r64();
                seg.pairs.data.resize(nPairs, std::vector<int>(2));
                for (uint64_t p = 0; p < nPairs; p++) {
                    int64_t a, b_val;
                    f.read(reinterpret_cast<char*>(&a), 8);
                    f.read(reinterpret_cast<char*>(&b_val), 8);
                    seg.pairs.data[p][0] = (a == -1) ? -1 : static_cast<int>(a - 1);
                    seg.pairs.data[p][1] = (b_val == -1) ? -1 : static_cast<int>(b_val - 1);
                }
            }

            // ── ppg_bin_indexs ──
            seg.ppg_bin_indexs = readPairArr();

            // ── ecg_bin_indexs ──
            seg.ecg_bin_indexs = readPairArr();

            // Derive bad_segment: bad if no pairs were detected
            seg.bad_segment = seg.pairs.data.empty() ? 1 : 0;
        }
        return segs;
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Template Info .bin  (written by your template generation step)
    // ─────────────────────────────────────────────────────────────────────────
    // Expected layout:
    //   Header:  [uint64 nSegments]
    //   Per seg: [uint64 ppgTemplateLen] [ppgTemplateLen × double]
    //            [uint64 ecgTemplateLen] [ecgTemplateLen × double]
    //            [double alignment_point]
    //
    // This file carries the computed templates but NOT the manual markings.
    // ═══════════════════════════════════════════════════════════════════════════

    inline std::vector<TemplateInfo> read_template_info_bin(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open())
            throw std::runtime_error("Cannot open template_info bin: " + path);

        f.seekg(0, std::ios::end);
        uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        auto r64 = [&]() -> uint64_t {
            uint64_t v = 0; f.read(reinterpret_cast<char*>(&v), 8);
            if (!f.good()) throw std::runtime_error("EOF in template at pos " + std::to_string(f.tellg()));
            return v;
            };
        auto rd = [&]() -> double {
            double v = 0; f.read(reinterpret_cast<char*>(&v), 8); return v;
            };
        auto readVecDbl = [&]() -> std::vector<double> {
            uint64_t n = r64();
            if (n > fileSize / 8) throw std::runtime_error(
                "vec too large " + std::to_string(n) + " at pos " + std::to_string(f.tellg()));
            std::vector<double> v(n);
            if (n > 0) f.read(reinterpret_cast<char*>(v.data()), n * 8);
            return v;
            };
        auto skipVecDbl = [&]() {
            uint64_t n = r64();
            if (n > fileSize / 8) throw std::runtime_error(
                "skip too large " + std::to_string(n) + " at pos " + std::to_string(f.tellg()));
            if (n > 0) f.seekg(n * 8, std::ios::cur);
            };
        auto skipPairs = [&]() {
            uint64_t n = r64();
            if (n > fileSize / 16) throw std::runtime_error(
                "pairs too large " + std::to_string(n) + " at pos " + std::to_string(f.tellg()));
            if (n > 0) f.seekg(n * 16, std::ios::cur);
            };

        uint64_t nT = r64();
        if (nT > 100000) throw std::runtime_error("nTemplates too large");

        std::vector<TemplateInfo> infos(nT);

        for (uint64_t i = 0; i < nT; i++) {
            uint64_t idx = r64();
            (void)idx;

            skipPairs();
            skipPairs();

            uint8_t bad = 0;
            f.read(reinterpret_cast<char*>(&bad), 1);
            infos[i].TemplateBad = bad ? 1 : 0;

            for (int ch = 0; ch < 3; ch++) {
                skipVecDbl();
                skipVecDbl();
                skipVecDbl();
                rd(); rd(); rd();
                rd(); rd(); rd();
            }

            infos[i].ppgTemplate = readVecDbl();

            infos[i].bad_r_templates = bad ? 1 : 0;
            infos[i].bad_ppg_templates = infos[i].ppgTemplate.empty() ? 1 : 0;
            infos[i].Dicrotic = NaN;
            infos[i].Onset = NaN;
            infos[i].Peak = NaN;
            infos[i].End = NaN;
        }
        return infos;
    }
    // ═══════════════════════════════════════════════════════════════════════════
    // Template Markings .bin  (written by writeTemplateMarkingsBin)
    // ─────────────────────────────────────────────────────────────────────────
    // Layout:
    //   [uint64 nBins]
    //   Per bin: [uint64 index]
    //            [uint8  bad_r_ch0] [uint8 bad_r_ch1] [uint8 bad_r_ch2]
    //            [uint8  ppg_issue]
    //            [int32  dicrotic] [int32 onset] [int32 peak] [int32 end_idx]
    // ═══════════════════════════════════════════════════════════════════════════

    inline void apply_template_markings_bin(const std::string& path,
        std::vector<TemplateInfo>& infos) {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) return;  // No markings file is OK (unreviewed)

        auto r64 = [&]() -> uint64_t {
            uint64_t v; f.read(reinterpret_cast<char*>(&v), 8); return v;
            };

        uint64_t nBins = r64();
        for (uint64_t i = 0; i < nBins; i++) {
            uint64_t index = r64();

            uint8_t bad_r0, bad_r1, bad_r2, ppg_issue;
            f.read(reinterpret_cast<char*>(&bad_r0), 1);
            f.read(reinterpret_cast<char*>(&bad_r1), 1);
            f.read(reinterpret_cast<char*>(&bad_r2), 1);
            f.read(reinterpret_cast<char*>(&ppg_issue), 1);

            int32_t dicrotic, onset, peak, end_idx;
            f.read(reinterpret_cast<char*>(&dicrotic), 4);
            f.read(reinterpret_cast<char*>(&onset), 4);
            f.read(reinterpret_cast<char*>(&peak), 4);
            f.read(reinterpret_cast<char*>(&end_idx), 4);

            if (index < infos.size()) {
                auto& t = infos[index];

                // All 3 ECG channels bad → bad R templates
                t.bad_r_templates = (bad_r0 && bad_r1 && bad_r2) ? 1 : 0;
                t.bad_ppg_templates = ppg_issue ? 1 : 0;
                t.TemplateBad = (t.bad_r_templates || t.bad_ppg_templates) ? 1 : 0;

                t.Dicrotic = (dicrotic >= 0) ? static_cast<double>(dicrotic) : NaN;
                t.Onset = (onset >= 0) ? static_cast<double>(onset) : NaN;
                t.Peak = (peak >= 0) ? static_cast<double>(peak) : NaN;
                t.End = (end_idx >= 0) ? static_cast<double>(end_idx) : NaN;
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Feature Output .bin  (replaces save to .mat)
    // ─────────────────────────────────────────────────────────────────────────
    // Layout:
    //   [uint64 nBeats]
    //   [uint64 nFields]
    //   Per field: [uint64 nameLen] [chars...] [nBeats × double]
    //
    // Also writes ppg_wout_noise as a separate trailing block:
    //   [uint64 ppgLen] [ppgLen × double]
    // ═══════════════════════════════════════════════════════════════════════════

    inline void write_feature_output_bin(const std::string& path,
        const BeatsFlattened& flat,
        int nBeats) {
        std::ofstream f(path, std::ios::binary);
        if (!f.is_open())
            throw std::runtime_error("Cannot write output: " + path);

        auto w64 = [&](uint64_t v) { f.write(reinterpret_cast<const char*>(&v), 8); };
        auto writeVec = [&](const std::vector<double>& v) {
            uint64_t n = v.size();
            w64(n);
            if (n > 0)
                f.write(reinterpret_cast<const char*>(v.data()), n * 8);
            };
        auto writeField = [&](const std::string& name, const std::vector<double>& data) {
            uint64_t nameLen = name.size();
            w64(nameLen);
            f.write(name.data(), nameLen);
            writeVec(data);
            };

        w64(static_cast<uint64_t>(nBeats));

        // Collect all named fields into a list for writing
        std::vector<std::pair<std::string, const std::vector<double>*>> all_fields;

#define ADD_FIELD(member) all_fields.push_back({#member, &flat.member})

        // Index fields
        ADD_FIELD(idx_begin); ADD_FIELD(idx_end); ADD_FIELD(idx_foot);
        ADD_FIELD(idx_pos_slope); ADD_FIELD(idx_systolic);
        ADD_FIELD(idx_neg_slope_b4); ADD_FIELD(idx_neg_slope_after);
        ADD_FIELD(idx_diastolic); ADD_FIELD(idx_dnotch);

        // tP/tR x-values
        ADD_FIELD(tP_20_x); ADD_FIELD(tP_50_x); ADD_FIELD(tP_80_x);
        ADD_FIELD(tP_20_x_inv); ADD_FIELD(tP_50_x_inv); ADD_FIELD(tP_80_x_inv);
        ADD_FIELD(tR_20_x); ADD_FIELD(tR_50_x); ADD_FIELD(tR_80_x);
        ADD_FIELD(tR_20_x_inv); ADD_FIELD(tR_50_x_inv); ADD_FIELD(tR_80_x_inv);

        // Scalars
        ADD_FIELD(sleep_stages); ADD_FIELD(area_baselined); ADD_FIELD(area);
        ADD_FIELD(amp_delta_systolic); ADD_FIELD(abs_amp_foot); ADD_FIELD(abs_amp_peak);
        ADD_FIELD(msec_beat_length);

        // tP50-relative timings
        ADD_FIELD(msec_tP_50_2_first_valley); ADD_FIELD(msec_tP_50_2_foot);
        ADD_FIELD(msec_tP_50_2_tP_20); ADD_FIELD(msec_tP_50_2_tP_80);
        ADD_FIELD(msec_tP_50_2_tP_20_inv); ADD_FIELD(msec_tP_50_2_tP_80_inv);
        ADD_FIELD(msec_tP_50_2_pos_slope); ADD_FIELD(msec_tP_50_2_systolic_peak);
        ADD_FIELD(msec_tP_50_2_negslopes_pre_dnotch); ADD_FIELD(msec_tP_50_2_dicrotic_notch);
        ADD_FIELD(msec_tP_50_2_diastolic_peak); ADD_FIELD(msec_tP_50_2_negslopes_post_dnotch);
        ADD_FIELD(msec_tP_50_2_second_valley);
        ADD_FIELD(msec_tP_50_2_tR_20); ADD_FIELD(msec_tP_50_2_tR_50); ADD_FIELD(msec_tP_50_2_tR_80);
        ADD_FIELD(msec_tP_50_2_tR_20_inv); ADD_FIELD(msec_tP_50_2_tR_50_inv); ADD_FIELD(msec_tP_50_2_tR_80_inv);
        ADD_FIELD(msec_total_duration_20); ADD_FIELD(msec_total_duration_50); ADD_FIELD(msec_total_duration_80);
        ADD_FIELD(msec_total_duration_tR_20); ADD_FIELD(msec_total_duration_tR_50); ADD_FIELD(msec_total_duration_tR_80);

        // R-peak timings
        ADD_FIELD(msec_R_2_first_valley); ADD_FIELD(msec_R_2_foot);
        ADD_FIELD(msec_R_2_tP_20); ADD_FIELD(msec_R_2_tP_50); ADD_FIELD(msec_R_2_tP_80);
        ADD_FIELD(msec_R_2_tP_20_inv); ADD_FIELD(msec_R_2_tP_50_inv); ADD_FIELD(msec_R_2_tP_80_inv);
        ADD_FIELD(msec_R_2_pos_slope); ADD_FIELD(msec_R_2_systolic_peak);
        ADD_FIELD(msec_R_2_negslopes_pre_dnotch); ADD_FIELD(msec_R_2_dicrotic_notch);
        ADD_FIELD(msec_R_2_tR_20); ADD_FIELD(msec_R_2_tR_50); ADD_FIELD(msec_R_2_tR_80);
        ADD_FIELD(msec_R_2_tR_20_inv); ADD_FIELD(msec_R_2_tR_50_inv); ADD_FIELD(msec_R_2_tR_80_inv);
        ADD_FIELD(msec_R_2_diastolic_peak); ADD_FIELD(msec_R_2_negslopes_post_dnotch);
        ADD_FIELD(msec_R_2_second_valley);

        // Raw amplitudes
        ADD_FIELD(amp_raw_vallies); ADD_FIELD(amp_raw_feets);
        ADD_FIELD(amp_raw_tP_20); ADD_FIELD(amp_raw_tP_50); ADD_FIELD(amp_raw_tP_80);
        ADD_FIELD(amp_raw_tP_20_inv); ADD_FIELD(amp_raw_tP_50_inv); ADD_FIELD(amp_raw_tP_80_inv);
        ADD_FIELD(amp_raw_tR_20); ADD_FIELD(amp_raw_tR_50); ADD_FIELD(amp_raw_tR_80);
        ADD_FIELD(amp_raw_tR_20_inv); ADD_FIELD(amp_raw_tR_50_inv); ADD_FIELD(amp_raw_tR_80_inv);
        ADD_FIELD(amp_raw_pos_slopes); ADD_FIELD(amp_raw_systolic_peaks);
        ADD_FIELD(amp_raw_neg_slopes_pre_dnotch); ADD_FIELD(amp_raw_dicrotic_notches);
        ADD_FIELD(amp_raw_diastolic_peaks); ADD_FIELD(amp_raw_neg_slopes_after_dnotch);

        // Baselined amplitudes
        ADD_FIELD(amp_baselined_feets);
        ADD_FIELD(amp_baselined_tP_20); ADD_FIELD(amp_baselined_tP_50); ADD_FIELD(amp_baselined_tP_80);
        ADD_FIELD(amp_baselined_tP_20_inv); ADD_FIELD(amp_baselined_tP_50_inv); ADD_FIELD(amp_baselined_tP_80_inv);
        ADD_FIELD(amp_baselined_tR_20); ADD_FIELD(amp_baselined_tR_50); ADD_FIELD(amp_baselined_tR_80);
        ADD_FIELD(amp_baselined_tR_20_inv); ADD_FIELD(amp_baselined_tR_50_inv); ADD_FIELD(amp_baselined_tR_80_inv);
        ADD_FIELD(amp_baselined_pos_slopes); ADD_FIELD(amp_baselined_systolic_peaks);
        ADD_FIELD(amp_baselined_neg_slopes_pre_dnotch); ADD_FIELD(amp_baselined_dicrotic_notches);
        ADD_FIELD(amp_baselined_diastolic_peaks); ADD_FIELD(amp_baselined_neg_slopes_after_dnotch);
        ADD_FIELD(proportional_pulse_amp);

        // SQI
        ADD_FIELD(sqi_mean_corr_dtw); ADD_FIELD(sqi_corrcoff_direct);
        ADD_FIELD(sqi_corrcoff_interp); ADD_FIELD(sqi_dtw); ADD_FIELD(sqi_frechet);

        // Inter-beat intervals
        ADD_FIELD(sec_valley_2_valley); ADD_FIELD(sec_foot_2_foot);
        ADD_FIELD(sec_tP_20_2_tP_20); ADD_FIELD(sec_tP_50_2_tP_50); ADD_FIELD(sec_tP_80_2_tP_80);
        ADD_FIELD(sec_tP_20_inv_2_tP_20_inv); ADD_FIELD(sec_tP_50_inv_2_tP_50_inv); ADD_FIELD(sec_tP_80_inv_2_tP_80_inv);
        ADD_FIELD(sec_pos_slope_2_pos_slope); ADD_FIELD(sec_systolic_2_systolic);
        ADD_FIELD(sec_neg_slope_b4_2_neg_slope_b4); ADD_FIELD(sec_neg_slope_after_2_neg_slope_after);
        ADD_FIELD(sec_diastolic_2_diastolic); ADD_FIELD(sec_dnotch_2_dnotch);
        ADD_FIELD(sec_tR_20_2_tR_20); ADD_FIELD(sec_tR_50_2_tR_50); ADD_FIELD(sec_tR_80_2_tR_80);
        ADD_FIELD(sec_tR_20_inv_2_tR_20_inv); ADD_FIELD(sec_tR_50_inv_2_tR_50_inv); ADD_FIELD(sec_tR_80_inv_2_tR_80_inv);

        // Special fields
        ADD_FIELD(adjusted_sleep_state);
        ADD_FIELD(corrected_time_sec);
        ADD_FIELD(sec_to_first_onset_of_sleep);
        ADD_FIELD(sec_from_last_onset_of_sleep);
        ADD_FIELD(ppg_flat_time_msec);

#undef ADD_FIELD

        w64(static_cast<uint64_t>(all_fields.size() + 1)); // +1 for edge_beat_mask

        for (auto& [name, data] : all_fields) {
            writeField(name, *data);
        }

        // edge_beat_mask as doubles
        {
            std::vector<double> edge(flat.edge_beat_mask.begin(), flat.edge_beat_mask.end());
            writeField("edge_beat_mask", edge);
        }

        // ppg_wout_noise as trailing block
        writeVec(flat.ppg_wout_noise);
    }

    // ═══════════════════════════════════════════════════════════════════════════
    // Setup: find matching files across directories
    // ═══════════════════════════════════════════════════════════════════════════

    struct BinPaths {
        std::string uuid;
        std::string anneal_path;      // *_annealed.bin  or *.bin
        std::string wave_path;        // *_wave_data.bin
        std::string template_path;    // *_template_info.bin  (empty if not using templates)
        std::string marking_path;     // *_template_markings.bin
    };

    inline std::vector<BinPaths> find_analysis_files(
        const std::string& annealed_dir,
        const std::string& wave_dir,
        const std::string& template_dir,
        const std::string& marking_dir,
        bool use_templates)
    {
        namespace fs = std::filesystem;

        // Collect .bin stems from a directory, optionally stripping a suffix.
        // If suffix is empty, the full stem is the ID.
        // If suffix is non-empty, only files ending with that suffix match,
        // and the suffix is stripped to produce the ID.
        auto collect_ids = [](const std::string& dir,
            const std::string& suffix = "")
            -> std::vector<std::string>
            {
                std::vector<std::string> ids;
                if (!fs::exists(dir)) return ids;
                for (auto& entry : fs::recursive_directory_iterator(dir)) {
                    if (!entry.is_regular_file()) continue;
                    if (entry.path().extension() != ".bin" &&
                        entry.path().extension() != ".BIN") continue;
                    std::string stem = entry.path().stem().string();
                    if (suffix.empty()) {
                        ids.push_back(stem);
                    }
                    else {
                        auto pos = stem.rfind(suffix);
                        if (pos != std::string::npos &&
                            pos + suffix.size() == stem.size()) {
                            ids.push_back(stem.substr(0, pos));
                        }
                    }
                }
                std::sort(ids.begin(), ids.end());
                ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                return ids;
            };

        // ── Annealed: <id>.bin (full stem = id) ──
        auto anneal_ids = collect_ids(annealed_dir);
        std::cout << "  Annealed dir: " << annealed_dir
            << "  ->  " << anneal_ids.size() << " .bin files\n";

        // ── Wave: try known suffixes, fall back to bare <id>.bin ──
        auto wave_ids = collect_ids(wave_dir, "_wave_data");
        std::string wave_suffix = "_wave_data";
        if (wave_ids.empty()) {
            wave_ids = collect_ids(wave_dir, "_wave_markings");
            wave_suffix = "_wave_markings";
        }
        if (wave_ids.empty()) {
            wave_ids = collect_ids(wave_dir);
            wave_suffix = "";
        }
        std::cout << "  Wave dir:     " << wave_dir
            << "  ->  " << wave_ids.size() << " .bin files"
            << (wave_suffix.empty() ? " (bare stems)"
                : (" (*" + wave_suffix + ")").c_str()) << "\n";

        // Intersect annealed ∩ wave
        std::vector<std::string> common;
        std::set_intersection(anneal_ids.begin(), anneal_ids.end(),
            wave_ids.begin(), wave_ids.end(),
            std::back_inserter(common));
        std::cout << "  Annealed ∩ Wave: " << common.size() << " matched IDs\n";

        if (use_templates) {
            auto tmpl_ids = collect_ids(template_dir, "_template_info");
            std::cout << "  Template dir: " << template_dir
                << "  ->  " << tmpl_ids.size() << " *_template_info.bin files\n";
            std::vector<std::string> common2;
            std::set_intersection(common.begin(), common.end(),
                tmpl_ids.begin(), tmpl_ids.end(),
                std::back_inserter(common2));
            common = common2;
            std::cout << "  After template intersect: " << common.size() << " matched IDs\n";
        }

        // Print a few sample IDs for debugging
        if (!common.empty()) {
            int show = std::min(3, (int)common.size());
            std::cout << "  Sample IDs: ";
            for (int i = 0; i < show; i++) std::cout << common[i] << "  ";
            std::cout << (common.size() > 3 ? "..." : "") << "\n";
        }
        else if (!anneal_ids.empty() && !wave_ids.empty()) {
            std::cout << "  WARNING: No ID overlap. Sample annealed ID: " << anneal_ids[0]
                << "  Sample wave ID: " << wave_ids[0] << "\n";
        }

        std::vector<BinPaths> result;
        for (auto& id : common) {
            BinPaths bp;
            bp.uuid = id;
            bp.anneal_path = (fs::path(annealed_dir) / (id + ".bin")).string();
            bp.wave_path = (fs::path(wave_dir) / (id + wave_suffix + ".bin")).string();
            if (use_templates) {
                bp.template_path = (fs::path(template_dir) / (id + "_template_info.bin")).string();
                bp.marking_path = (fs::path(marking_dir) / (id + "_template_markings.bin")).string();
            }
            result.push_back(bp);
        }
        return result;
    }

} // namespace ppg