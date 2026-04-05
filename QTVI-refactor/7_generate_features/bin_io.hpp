#pragma once

#include "data_types.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <iostream>

/**
 * @file bin_io.hpp
 * @brief Binary file readers and writers for every stage of the pipeline.
 *
 * Formats are derived directly from:
 *   - 3_anneal_segments.cpp  → write_output_bin       (annealed segments)
 *   - main.cpp               → write_output_binfile   (wave markings)
 *   - 5_generate_templates.cpp → write_template_info_binfile (templates)
 *
 * All multi-byte values are little-endian (x86 native).
 * Index arrays written by main.cpp are **1-based** for MATLAB compatibility;
 * readers here convert back to 0-based.
 */

namespace ppg {
    namespace io {

        // ═══════════════════════════════════════════════════════════════════════════
        //  Low-level helpers
        // ═══════════════════════════════════════════════════════════════════════════

        namespace detail {

            inline void read_bytes(std::ifstream& f, void* dst, size_t n) {
                if (!f.read(reinterpret_cast<char*>(dst), static_cast<std::streamsize>(n)))
                    throw std::runtime_error("bin_io: unexpected end of file");
            }

            inline uint64_t read_u64(std::ifstream& f) {
                uint64_t v = 0;
                read_bytes(f, &v, 8);
                return v;
            }

            inline double read_f64(std::ifstream& f) {
                double v = 0.0;
                read_bytes(f, &v, 8);
                return v;
            }

            inline uint8_t read_u8(std::ifstream& f) {
                uint8_t v = 0;
                read_bytes(f, &v, 1);
                return v;
            }

            /**
             * @brief Read a length-prefixed double array.
             * @param f    Input stream.
             * @param dest Destination vector (resized).
             */
            inline void read_vec_f64(std::ifstream& f, std::vector<double>& dest) {
                uint64_t sz = read_u64(f);
                if (sz > 50'000'000) throw std::runtime_error("bin_io: vector size too large");
                dest.resize(sz);
                if (sz > 0) read_bytes(f, dest.data(), sz * 8);
            }

            /**
             * @brief Read a length-prefixed index array (stored as uint64, 1-based) → 0-based int vector.
             * @param f    Input stream.
             * @param dest Destination vector.
             */
            inline void read_idx_array(std::ifstream& f, std::vector<int>& dest) {
                uint64_t sz = read_u64(f);
                if (sz > 50'000'000) throw std::runtime_error("bin_io: index array too large");
                dest.resize(sz);
                if (sz > 0) {
                    std::vector<uint64_t> tmp(sz);
                    read_bytes(f, tmp.data(), sz * 8);
                    for (uint64_t i = 0; i < sz; ++i)
                        dest[i] = static_cast<int>(tmp[i] > 0 ? tmp[i] - 1 : 0);
                }
            }

            /**
             * @brief Read a length-prefixed pair<uint64,uint64> array.
             * @param f    Input stream.
             * @param dest Destination vector.
             */
            inline void read_pair_vec(std::ifstream& f, std::vector<std::pair<int, int>>& dest) {
                uint64_t sz = read_u64(f);
                if (sz > 50'000'000) throw std::runtime_error("bin_io: pair vector too large");
                dest.resize(sz);
                if (sz > 0) {
                    std::vector<uint64_t> buf(sz * 2);
                    read_bytes(f, buf.data(), sz * 16);
                    for (uint64_t i = 0; i < sz; ++i) {
                        dest[i].first = static_cast<int>(buf[i * 2]);
                        dest[i].second = static_cast<int>(buf[i * 2 + 1]);
                    }
                }
            }

            // ── Write helpers ────────────────────────────────────────────────────────

            inline void write_bytes(std::ofstream& f, const void* src, size_t n) {
                f.write(reinterpret_cast<const char*>(src), static_cast<std::streamsize>(n));
            }

            inline void write_u64(std::ofstream& f, uint64_t v) { write_bytes(f, &v, 8); }
            inline void write_f64(std::ofstream& f, double v) { write_bytes(f, &v, 8); }
            inline void write_u8(std::ofstream& f, uint8_t v) { write_bytes(f, &v, 1); }

            inline void write_vec_f64(std::ofstream& f, const std::vector<double>& v) {
                write_u64(f, v.size());
                if (!v.empty()) write_bytes(f, v.data(), v.size() * 8);
            }

            /**
             * @brief Write 0-based int indices as 1-based uint64 array (MATLAB compat).
             */
            inline void write_idx_array(std::ofstream& f, const std::vector<int>& v) {
                write_u64(f, v.size());
                if (!v.empty()) {
                    std::vector<uint64_t> tmp(v.size());
                    for (size_t i = 0; i < v.size(); ++i)
                        tmp[i] = static_cast<uint64_t>(v[i]) + 1;
                    write_bytes(f, tmp.data(), v.size() * 8);
                }
            }

            inline void write_pair_vec_u64(
                std::ofstream& f,
                const std::vector<std::pair<int, int>>& v)
            {
                write_u64(f, v.size());
                for (const auto& p : v) {
                    uint64_t a = static_cast<uint64_t>(p.first);
                    uint64_t b = static_cast<uint64_t>(p.second);
                    write_bytes(f, &a, 8);
                    write_bytes(f, &b, 8);
                }
            }

        } // namespace detail


        // ═══════════════════════════════════════════════════════════════════════════
        //  1. Annealed Segments (.bin)
        //
        //  Written by 3_anneal_segments.cpp → write_output_bin:
        //    [u64  nSegments]
        //    [f64  ppgSR]  [f64  ecgSR]  [f64  epochSec]
        //    per segment:
        //      ppg_bin_indexs (pair vec)
        //      ecg_bin_indexs (pair vec)
        //      ppg   (f64 vec)
        //      ecg1  (f64 vec)
        //      ecg2  (f64 vec)
        //      ecg3  (f64 vec)
        //      sleep (f64 vec)
        // ═══════════════════════════════════════════════════════════════════════════

        /**
         * @brief Read annealed segments from a .bin file.
         * @param path  File path.
         * @return Vector of AnnealedSegment structs.
         */
        inline std::vector<AnnealedSegment> read_annealed_bin(const std::string& path) {
            std::ifstream f(path, std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open annealed bin: " + path);

            uint64_t n_seg = detail::read_u64(f);

            double ppg_sr = 256.0, ecg_sr = 256.0, epoch_sec = 30.0;
            if (n_seg > 0) {
                ppg_sr = detail::read_f64(f);
                ecg_sr = detail::read_f64(f);
                epoch_sec = detail::read_f64(f);
            }

            std::vector<AnnealedSegment> segs(n_seg);
            for (uint64_t i = 0; i < n_seg; ++i) {
                auto& s = segs[i];
                s.ppg_sample_rate = ppg_sr;
                s.ecg_sample_rate = ecg_sr;

                detail::read_pair_vec(f, s.ppg_bin_indices);

                // ecg_bin_indices (read but stored loosely — not needed downstream)
                std::vector<std::pair<int, int>> ecg_bins;
                detail::read_pair_vec(f, ecg_bins);

                detail::read_vec_f64(f, s.po);

                // ecg1, ecg2, ecg3 — skip (not used by feature extraction)
                std::vector<double> tmp;
                detail::read_vec_f64(f, tmp);
                detail::read_vec_f64(f, tmp);
                detail::read_vec_f64(f, tmp);

                detail::read_vec_f64(f, s.sleep_stages);
            }
            return segs;
        }


        // ═══════════════════════════════════════════════════════════════════════════
        //  2. Wave Markings (.bin)
        //
        //  Written by main.cpp → write_output_binfile:
        //    [u64 numBins]
        //    per bin:
        //      9 index arrays  (ch1 raw/sq/abs, ch2 raw/sq/abs, ch3 raw/sq/abs)
        //      2 PPG idx arrays (maxAmps, minAmps)
        //      4 raw signals    (ppg, ecg1, ecg2, ecg3)
        //      6 preproc sigs   (ch1 sq/abs, ch2 sq/abs, ch3 sq/abs)
        //      9 noise flags    (1 byte each)
        //      pairs:           [u64 count] [int64 × count × 2]
        //      ppg_bin_indexs   (pair vec)
        //      ecg_bin_indexs   (pair vec)
        // ═══════════════════════════════════════════════════════════════════════════

        /**
         * @brief Read wave markings from a .bin file.
         * @param path  File path.
         * @return Vector of WaveData structs.
         */
        inline std::vector<WaveData> read_wave_data_bin(const std::string& path) {
            std::ifstream f(path, std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open wave bin: " + path);

            uint64_t n_bins = detail::read_u64(f);
            std::vector<WaveData> result(n_bins);

            for (uint64_t b = 0; b < n_bins; ++b) {
                auto& wd = result[b];

                // 9 R-peak index arrays — skip (not used directly by feature extraction)
                std::vector<int> skip_idx;
                for (int k = 0; k < 9; ++k) detail::read_idx_array(f, skip_idx);

                // 2 PPG index arrays — skip
                for (int k = 0; k < 2; ++k) detail::read_idx_array(f, skip_idx);

                // 4 raw signals — skip
                std::vector<double> skip_sig;
                for (int k = 0; k < 4; ++k) detail::read_vec_f64(f, skip_sig);

                // 6 preprocessed signals — skip
                for (int k = 0; k < 6; ++k) detail::read_vec_f64(f, skip_sig);

                // 9 noise flags
                uint8_t flags[9];
                detail::read_bytes(f, flags, 9);

                // Pairs: [u64 count] [int64 × count × 2]
                // 1-based, -1 = unmatched
                uint64_t num_pairs = detail::read_u64(f);
                wd.pairs.resize(num_pairs);
                if (num_pairs > 0) {
                    std::vector<int64_t> pair_buf(num_pairs * 2);
                    detail::read_bytes(f, pair_buf.data(), num_pairs * 16);
                    for (uint64_t i = 0; i < num_pairs; ++i) {
                        // Convert 1-based → 0-based; -1 stays as -1
                        int ppg_idx = (pair_buf[i * 2] == -1)
                            ? -1 : static_cast<int>(pair_buf[i * 2] - 1);
                        int ecg_idx = (pair_buf[i * 2 + 1] == -1)
                            ? -1 : static_cast<int>(pair_buf[i * 2 + 1] - 1);
                        wd.pairs[i] = { ppg_idx, ecg_idx };
                    }
                }

                // ppg/ecg bin index pairs — skip
                std::vector<std::pair<int, int>> skip_pairs;
                detail::read_pair_vec(f, skip_pairs);
                detail::read_pair_vec(f, skip_pairs);

                // Determine bad_segment heuristically
                // (matches logic in 5_generate_templates.cpp)
                wd.bad_segment = (num_pairs == 0);
            }
            return result;
        }


        // ═══════════════════════════════════════════════════════════════════════════
        //  3. Template Info (.bin)
        //
        //  Written by 5_generate_templates.cpp → write_template_info_binfile:
        //    [u64 nInfos]
        //    per info:
        //      [u64 index]
        //      ppg_bin_indexs  (pair vec)
        //      ecg_bin_indexs  (pair vec)
        //      [u8  bad_segment]
        //      3 × channel template:
        //        ecgTemplate_raw     (f64 vec)
        //        ecgTemplate_squared (f64 vec)
        //        ecgTemplate_absval  (f64 vec)
        //        alignment_point_raw     (f64)
        //        alignment_point_squared (f64)
        //        alignment_point_absval  (f64)
        //        avg_r_expand_raw        (f64)
        //        avg_r_expand_squared    (f64)
        //        avg_r_expand_absval     (f64)
        //      ppgTemplate     (f64 vec)
        // ═══════════════════════════════════════════════════════════════════════════

        /**
         * @brief Read template info from a .bin file produced by 5_generate_templates.
         *
         * Note: The step-5 binary does NOT contain Dicrotic/Onset/Peak/End annotations.
         * Those come from the manual marking step (step 6). If a marking .bin exists
         * for the same subject, call read_template_markings_bin() to overlay those
         * values onto the TemplateInfo vector.
         *
         * @param path  File path.
         * @return Vector of TemplateInfo structs.
         */
        inline std::vector<TemplateInfo> read_template_info_bin(const std::string& path) {
            std::ifstream f(path, std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open template bin: " + path);

            uint64_t n = detail::read_u64(f);
            std::vector<TemplateInfo> result(n);

            for (uint64_t i = 0; i < n; ++i) {
                auto& ti = result[i];

                /* index */ detail::read_u64(f);

                // bin index pairs — skip
                std::vector<std::pair<int, int>> skip_pairs;
                detail::read_pair_vec(f, skip_pairs);
                detail::read_pair_vec(f, skip_pairs);

                // bad_segment flag
                ti.template_bad = (detail::read_u8(f) != 0);

                // 3 channel templates — skip ECG template data,
                // but read through to keep stream position correct
                for (int ch = 0; ch < 3; ++ch) {
                    std::vector<double> skip_vec;
                    detail::read_vec_f64(f, skip_vec);  // ecgTemplate_raw
                    detail::read_vec_f64(f, skip_vec);  // ecgTemplate_squared
                    detail::read_vec_f64(f, skip_vec);  // ecgTemplate_absval

                    double ap_raw = detail::read_f64(f);
                    /* ap_sq  */ detail::read_f64(f);
                    /* ap_abs */ detail::read_f64(f);
                    /* avg_r_raw */ detail::read_f64(f);
                    /* avg_r_sq  */ detail::read_f64(f);
                    /* avg_r_abs */ detail::read_f64(f);

                    // Use channel 1's alignment point as the main one
                    if (ch == 0) ti.alignment_point = ap_raw;
                }

                // PPG template
                detail::read_vec_f64(f, ti.ppg_template);

                // Derive flags from template content
                ti.bad_ppg_templates = ti.ppg_template.empty();
                ti.bad_r_templates = ti.template_bad;

                // Dicrotic/Onset/Peak/End default to NaN — must be overlaid from
                // marking files if available (see read_template_markings_bin).
            }
            return result;
        }

        // ═══════════════════════════════════════════════════════════════════════════
        //  3b. Template Markings (.bin) — from step 6 GUI
        //
        //  Written by writeTemplateMarkingsBin() in TemplateBinIO.hpp:
        //    [u64  numBins]
        //    per bin:
        //      [u64  index]
        //      [u8   bad_r_ch1]
        //      [u8   bad_r_ch2]
        //      [u8   bad_r_ch3]
        //      [u8   ppg_issue]     0=ok, 1=bad, 2=no ppg
        //      [i32  dicrotic]      -1 = NaN
        //      [i32  onset]         -1 = NaN
        //      [i32  peak]          -1 = NaN
        //      [i32  end_idx]       -1 = NaN
        // ═══════════════════════════════════════════════════════════════════════════

        /**
         * @brief Read manually-reviewed template markings and overlay onto TemplateInfo.
         *
         * Format matches TemplateBinIO.hpp → writeTemplateMarkingsBin().
         * If the file does not exist, this is a no-op and the pipeline uses
         * automatic dicrotic notch estimation (NaN pathway).
         *
         * @param path       Path to _template_markings.bin.
         * @param templates  TemplateInfo vector to update in-place.
         */
        inline void read_template_markings_bin(
            const std::string& path,
            std::vector<TemplateInfo>& templates)
        {
            std::ifstream f(path, std::ios::binary);
            if (!f) return;

            uint64_t count = detail::read_u64(f);

            for (uint64_t i = 0; i < count && i < templates.size(); ++i) {
                /* index */ detail::read_u64(f);

                uint8_t bad_r_ch1 = detail::read_u8(f);
                uint8_t bad_r_ch2 = detail::read_u8(f);
                uint8_t bad_r_ch3 = detail::read_u8(f);
                uint8_t ppg_issue = detail::read_u8(f);

                int32_t dicrotic_i32 = 0, onset_i32 = 0, peak_i32 = 0, end_i32 = 0;
                detail::read_bytes(f, &dicrotic_i32, 4);
                detail::read_bytes(f, &onset_i32, 4);
                detail::read_bytes(f, &peak_i32, 4);
                detail::read_bytes(f, &end_i32, 4);

                // Any R channel flagged bad → bad_r_templates
                templates[i].bad_r_templates = (bad_r_ch1 != 0 || bad_r_ch2 != 0 || bad_r_ch3 != 0);

                // ppg_issue: 0=ok, 1=bad, 2=no ppg
                templates[i].bad_ppg_templates = (ppg_issue != 0);
                templates[i].template_bad = (ppg_issue != 0)
                    || (bad_r_ch1 != 0 && bad_r_ch2 != 0 && bad_r_ch3 != 0);

                // -1 sentinel → NaN
                templates[i].dicrotic = (dicrotic_i32 == -1) ? kNaN : static_cast<double>(dicrotic_i32);
                templates[i].onset = (onset_i32 == -1) ? kNaN : static_cast<double>(onset_i32);
                templates[i].peak = (peak_i32 == -1) ? kNaN : static_cast<double>(peak_i32);
                templates[i].end = (end_i32 == -1) ? kNaN : static_cast<double>(end_i32);
            }
        }


        // ═══════════════════════════════════════════════════════════════════════════
        //  4. Feature Output (.bin)
        //
        //  Custom format for the pipeline output (replaces the .mat output).
        //  Layout:
        //    [u64  numBeats]
        //    [f64  ppgSampleRate]
        //    per beat:
        //      All BeatFeatures fields as consecutive f64 values.
        //      [u64 sqi_count] [f64 × sqi_count] for SQI vector.
        //    Trailing arrays (one f64 per beat):
        //      sec_valley_to_valley
        //      sec_foot_to_foot
        //      adjusted_sleep_state
        //      corrected_time_sec
        //    [u64 ppg_length] [f64 × ppg_length] for ppg_wout_noise
        // ═══════════════════════════════════════════════════════════════════════════

        /**
         * @brief Write flattened feature output to a .bin file.
         * @param path  Output file path.
         * @param flat  Flattened beat features.
         */
        inline void write_feature_output_bin(const std::string& path, const FlattenedBeats& flat) {
            std::ofstream f(path, std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open output bin: " + path);

            char buf[1 << 16];
            f.rdbuf()->pubsetbuf(buf, sizeof(buf));

            const uint64_t n = flat.beats.size();
            detail::write_u64(f, n);
            detail::write_f64(f, flat.ppg_sample_rate);

            // Write each beat's scalar features as a fixed-width record
            for (const auto& b : flat.beats) {
                auto w = [&](double v) { detail::write_f64(f, v); };
                auto wi = [&](int v) { detail::write_f64(f, static_cast<double>(v)); };

                // Indices
                wi(b.idx_begin); wi(b.idx_end); wi(b.idx_foot);
                wi(b.idx_pos_slope); wi(b.idx_systolic);
                wi(b.idx_neg_slope_b4); wi(b.idx_dnotch);
                wi(b.idx_diastolic); wi(b.idx_neg_slope_after);

                // tP positions
                w(b.tP_20_x); w(b.tP_50_x); w(b.tP_80_x);
                w(b.tP_20_x_inv); w(b.tP_50_x_inv); w(b.tP_80_x_inv);

                // tR positions
                w(b.tR_20_x); w(b.tR_50_x); w(b.tR_80_x);
                w(b.tR_20_x_inv); w(b.tR_50_x_inv); w(b.tR_80_x_inv);

                // tP amplitudes
                w(b.tP_20_y); w(b.tP_50_y); w(b.tP_80_y);
                w(b.tP_20_y_inv); w(b.tP_50_y_inv); w(b.tP_80_y_inv);

                // tR amplitudes
                w(b.tR_20_y); w(b.tR_50_y); w(b.tR_80_y);
                w(b.tR_20_y_inv); w(b.tR_50_y_inv); w(b.tR_80_y_inv);

                // tP50 timings
                w(b.msec_tP50_to_valley1); w(b.msec_tP50_to_foot);
                w(b.msec_tP50_to_tP20); w(b.msec_tP50_to_tP80);
                w(b.msec_tP50_to_tP20_inv); w(b.msec_tP50_to_tP80_inv);
                w(b.msec_tP50_to_pos_slope); w(b.msec_tP50_to_systolic);
                w(b.msec_tP50_to_neg_pre); w(b.msec_tP50_to_dnotch);
                w(b.msec_tP50_to_diastolic); w(b.msec_tP50_to_neg_post);
                w(b.msec_tP50_to_tR20); w(b.msec_tP50_to_tR50); w(b.msec_tP50_to_tR80);
                w(b.msec_tP50_to_tR20_inv); w(b.msec_tP50_to_tR50_inv); w(b.msec_tP50_to_tR80_inv);

                // Durations
                w(b.msec_beat_length);
                w(b.msec_dur_tP20); w(b.msec_dur_tP50); w(b.msec_dur_tP80);
                w(b.msec_dur_tR20); w(b.msec_dur_tR50); w(b.msec_dur_tR80);

                // Raw amplitudes
                w(b.amp_raw_valley); w(b.amp_raw_foot);
                w(b.amp_raw_tP20); w(b.amp_raw_tP50); w(b.amp_raw_tP80);
                w(b.amp_raw_tP20_inv); w(b.amp_raw_tP50_inv); w(b.amp_raw_tP80_inv);
                w(b.amp_raw_pos_slope); w(b.amp_raw_systolic);
                w(b.amp_raw_neg_pre); w(b.amp_raw_dnotch);
                w(b.amp_raw_diastolic); w(b.amp_raw_neg_post);
                w(b.amp_raw_tR20); w(b.amp_raw_tR50); w(b.amp_raw_tR80);
                w(b.amp_raw_tR20_inv); w(b.amp_raw_tR50_inv); w(b.amp_raw_tR80_inv);

                // Baselined amplitudes
                w(b.amp_bl_foot);
                w(b.amp_bl_tP20); w(b.amp_bl_tP50); w(b.amp_bl_tP80);
                w(b.amp_bl_tP20_inv); w(b.amp_bl_tP50_inv); w(b.amp_bl_tP80_inv);
                w(b.amp_bl_pos_slope); w(b.amp_bl_systolic);
                w(b.amp_bl_neg_pre); w(b.amp_bl_dnotch);
                w(b.amp_bl_diastolic); w(b.amp_bl_neg_post);
                w(b.amp_bl_tR20); w(b.amp_bl_tR50); w(b.amp_bl_tR80);
                w(b.amp_bl_tR20_inv); w(b.amp_bl_tR50_inv); w(b.amp_bl_tR80_inv);

                // Area & misc
                w(b.area); w(b.area_baselined);
                w(b.amp_delta_systolic); w(b.proportional_pulse_amp);
                w(b.abs_amp_foot); w(b.abs_amp_peak);

                // R-peak timings
                w(b.msec_R_to_valley1); w(b.msec_R_to_foot);
                w(b.msec_R_to_tP20); w(b.msec_R_to_tP50); w(b.msec_R_to_tP80);
                w(b.msec_R_to_tP20_inv); w(b.msec_R_to_tP50_inv); w(b.msec_R_to_tP80_inv);
                w(b.msec_R_to_pos_slope); w(b.msec_R_to_systolic);
                w(b.msec_R_to_neg_pre); w(b.msec_R_to_dnotch);
                w(b.msec_R_to_diastolic); w(b.msec_R_to_neg_post);
                w(b.msec_R_to_valley2);
                w(b.msec_R_to_tR20); w(b.msec_R_to_tR50); w(b.msec_R_to_tR80);
                w(b.msec_R_to_tR20_inv); w(b.msec_R_to_tR50_inv); w(b.msec_R_to_tR80_inv);

                // Sleep stage
                w(b.sleep_stage);

                // SQI vector
                detail::write_u64(f, b.sqi.size());
                for (double s : b.sqi) detail::write_f64(f, s);
            }

            // Trailing per-beat arrays
            auto write_arr = [&](const std::vector<double>& v) {
                detail::write_vec_f64(f, v);
                };
            write_arr(flat.sec_valley_to_valley);
            write_arr(flat.sec_foot_to_foot);
            write_arr(flat.adjusted_sleep_state);
            write_arr(flat.corrected_time_sec);

            // Edge beat mask as doubles
            {
                std::vector<double> mask(flat.edge_beat_mask.begin(), flat.edge_beat_mask.end());
                write_arr(mask);
            }

            // Concatenated clean PPG
            write_arr(flat.ppg_wout_noise);
        }

        /**
         * @brief Read flattened feature output from a .bin file.
         * @param path  Input file path.
         * @return FlattenedBeats struct.
         */
        inline FlattenedBeats read_feature_output_bin(const std::string& path) {
            std::ifstream f(path, std::ios::binary);
            if (!f) throw std::runtime_error("Cannot open feature bin: " + path);

            FlattenedBeats flat;
            uint64_t n = detail::read_u64(f);
            flat.ppg_sample_rate = detail::read_f64(f);
            flat.beats.resize(n);

            for (uint64_t i = 0; i < n; ++i) {
                auto& b = flat.beats[i];
                auto r = [&]() { return detail::read_f64(f); };
                auto ri = [&]() { return static_cast<int>(detail::read_f64(f)); };

                b.idx_begin = ri(); b.idx_end = ri(); b.idx_foot = ri();
                b.idx_pos_slope = ri(); b.idx_systolic = ri();
                b.idx_neg_slope_b4 = ri(); b.idx_dnotch = ri();
                b.idx_diastolic = ri(); b.idx_neg_slope_after = ri();

                b.tP_20_x = r(); b.tP_50_x = r(); b.tP_80_x = r();
                b.tP_20_x_inv = r(); b.tP_50_x_inv = r(); b.tP_80_x_inv = r();
                b.tR_20_x = r(); b.tR_50_x = r(); b.tR_80_x = r();
                b.tR_20_x_inv = r(); b.tR_50_x_inv = r(); b.tR_80_x_inv = r();

                b.tP_20_y = r(); b.tP_50_y = r(); b.tP_80_y = r();
                b.tP_20_y_inv = r(); b.tP_50_y_inv = r(); b.tP_80_y_inv = r();
                b.tR_20_y = r(); b.tR_50_y = r(); b.tR_80_y = r();
                b.tR_20_y_inv = r(); b.tR_50_y_inv = r(); b.tR_80_y_inv = r();

                b.msec_tP50_to_valley1 = r(); b.msec_tP50_to_foot = r();
                b.msec_tP50_to_tP20 = r(); b.msec_tP50_to_tP80 = r();
                b.msec_tP50_to_tP20_inv = r(); b.msec_tP50_to_tP80_inv = r();
                b.msec_tP50_to_pos_slope = r(); b.msec_tP50_to_systolic = r();
                b.msec_tP50_to_neg_pre = r(); b.msec_tP50_to_dnotch = r();
                b.msec_tP50_to_diastolic = r(); b.msec_tP50_to_neg_post = r();
                b.msec_tP50_to_tR20 = r(); b.msec_tP50_to_tR50 = r(); b.msec_tP50_to_tR80 = r();
                b.msec_tP50_to_tR20_inv = r(); b.msec_tP50_to_tR50_inv = r(); b.msec_tP50_to_tR80_inv = r();

                b.msec_beat_length = r();
                b.msec_dur_tP20 = r(); b.msec_dur_tP50 = r(); b.msec_dur_tP80 = r();
                b.msec_dur_tR20 = r(); b.msec_dur_tR50 = r(); b.msec_dur_tR80 = r();

                b.amp_raw_valley = r(); b.amp_raw_foot = r();
                b.amp_raw_tP20 = r(); b.amp_raw_tP50 = r(); b.amp_raw_tP80 = r();
                b.amp_raw_tP20_inv = r(); b.amp_raw_tP50_inv = r(); b.amp_raw_tP80_inv = r();
                b.amp_raw_pos_slope = r(); b.amp_raw_systolic = r();
                b.amp_raw_neg_pre = r(); b.amp_raw_dnotch = r();
                b.amp_raw_diastolic = r(); b.amp_raw_neg_post = r();
                b.amp_raw_tR20 = r(); b.amp_raw_tR50 = r(); b.amp_raw_tR80 = r();
                b.amp_raw_tR20_inv = r(); b.amp_raw_tR50_inv = r(); b.amp_raw_tR80_inv = r();

                b.amp_bl_foot = r();
                b.amp_bl_tP20 = r(); b.amp_bl_tP50 = r(); b.amp_bl_tP80 = r();
                b.amp_bl_tP20_inv = r(); b.amp_bl_tP50_inv = r(); b.amp_bl_tP80_inv = r();
                b.amp_bl_pos_slope = r(); b.amp_bl_systolic = r();
                b.amp_bl_neg_pre = r(); b.amp_bl_dnotch = r();
                b.amp_bl_diastolic = r(); b.amp_bl_neg_post = r();
                b.amp_bl_tR20 = r(); b.amp_bl_tR50 = r(); b.amp_bl_tR80 = r();
                b.amp_bl_tR20_inv = r(); b.amp_bl_tR50_inv = r(); b.amp_bl_tR80_inv = r();

                b.area = r(); b.area_baselined = r();
                b.amp_delta_systolic = r(); b.proportional_pulse_amp = r();
                b.abs_amp_foot = r(); b.abs_amp_peak = r();

                b.msec_R_to_valley1 = r(); b.msec_R_to_foot = r();
                b.msec_R_to_tP20 = r(); b.msec_R_to_tP50 = r(); b.msec_R_to_tP80 = r();
                b.msec_R_to_tP20_inv = r(); b.msec_R_to_tP50_inv = r(); b.msec_R_to_tP80_inv = r();
                b.msec_R_to_pos_slope = r(); b.msec_R_to_systolic = r();
                b.msec_R_to_neg_pre = r(); b.msec_R_to_dnotch = r();
                b.msec_R_to_diastolic = r(); b.msec_R_to_neg_post = r();
                b.msec_R_to_valley2 = r();
                b.msec_R_to_tR20 = r(); b.msec_R_to_tR50 = r(); b.msec_R_to_tR80 = r();
                b.msec_R_to_tR20_inv = r(); b.msec_R_to_tR50_inv = r(); b.msec_R_to_tR80_inv = r();

                b.sleep_stage = r();

                uint64_t sqi_n = detail::read_u64(f);
                b.sqi.resize(sqi_n);
                for (uint64_t j = 0; j < sqi_n; ++j) b.sqi[j] = detail::read_f64(f);
            }

            detail::read_vec_f64(f, flat.sec_valley_to_valley);
            detail::read_vec_f64(f, flat.sec_foot_to_foot);
            detail::read_vec_f64(f, flat.adjusted_sleep_state);
            detail::read_vec_f64(f, flat.corrected_time_sec);

            {
                std::vector<double> mask;
                detail::read_vec_f64(f, mask);
                flat.edge_beat_mask.resize(mask.size());
                for (size_t i = 0; i < mask.size(); ++i)
                    flat.edge_beat_mask[i] = static_cast<int>(mask[i]);
            }

            detail::read_vec_f64(f, flat.ppg_wout_noise);
            return flat;
        }

    }
} // namespace ppg::io