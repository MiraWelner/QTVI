#pragma once

#include "common.hpp"
#include "data_types.hpp"
#include "beat_features.hpp"
#include "ppg_sqi.hpp"
#include "flatten_beats.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <future>

/**
 * @file generate_features.hpp
 * @brief Top-level feature generation pipeline.
 *
 * Port of GenerateFeatures.m and GenerateFeaturesRunner.m.
 * Uses a single shared thread pool for all parallelism.
 * Work is chunked to avoid flooding the pool with tiny tasks.
 */

namespace ppg {

    struct BinResult {
        std::vector<BeatFeatures> beats;
        std::vector<BeatSqi>      sqi;
        BinMarks                  marks;
    };

    /**
     * @brief Process a single segment: compute SQI and extract beat features.
     */
    inline BinResult process_bin(
        const AnnealedSegment& segment,
        const WaveData& wd,
        const TemplateInfo& tmpl,
        double sqi_threshold,
        int window_length,
        ThreadPool& pool)
    {
        BinResult result;
        result.marks.could_not_identify_ppg = wd.bad_segment;
        result.marks.ppg_template_excluded = tmpl.bad_ppg_templates;
        result.marks.ecg_template_excluded = tmpl.bad_r_templates;
        result.marks.poor_template_excluded = tmpl.template_bad;

        const auto& ppg = segment.po;
        const auto& pairs = wd.pairs;
        const double ppg_fs = segment.ppg_sample_rate;
        const double ecg_fs = segment.ecg_sample_rate;

        if (pairs.size() < 2) return result;

        // Beat onset indices
        std::vector<int> ann_times;
        ann_times.reserve(pairs.size());
        for (auto& p : pairs) ann_times.push_back(p.first);

        // ── SQI (parallel via pool) ──
        int win_samples = window_length * static_cast<int>(ppg_fs);
        result.sqi = ppg_sqi(ppg, ann_times, tmpl.ppg_template, win_samples, ppg_fs, &pool);

        // Dicrotic notch ratio from template markings
        double dn_ratio = kNaN;
        double size_peak_to_end = tmpl.end - tmpl.peak;
        if (!std::isnan(size_peak_to_end) && size_peak_to_end >= 1.0)
            dn_ratio = std::abs(tmpl.dicrotic - tmpl.peak) / size_peak_to_end;

        // ── Identify valid beats ──
        const int num_beats = static_cast<int>(ann_times.size()) - 1;

        std::vector<int> valid_indices;
        valid_indices.reserve(num_beats);
        for (int i = 0; i < num_beats; ++i) {
            double sqi_val = (i < static_cast<int>(result.sqi.size()))
                ? result.sqi[i].mean_corr : 0.0;
            if (sqi_val < sqi_threshold)
                valid_indices.push_back(i);
        }

        const int nValid = static_cast<int>(valid_indices.size());
        if (nValid == 0) return result;

        // ── Feature extraction (chunked parallel) ──
        result.beats.resize(nValid);

        constexpr int CHUNK = 32;
        int nChunks = (nValid + CHUNK - 1) / CHUNK;

        std::vector<std::future<void>> futures;
        futures.reserve(nChunks);

        for (int ch = 0; ch < nChunks; ++ch) {
            int kStart = ch * CHUNK;
            int kEnd = std::min(kStart + CHUNK, nValid);

            futures.push_back(pool.enqueue(
                [&, kStart, kEnd]() {
                    for (int k = kStart; k < kEnd; ++k) {
                        int i = valid_indices[k];
                        int r_idx = (i < static_cast<int>(pairs.size())) ? pairs[i].second : -1;
                        int bb = ann_times[i];
                        int be = ann_times[i + 1];

                        auto feat = extract_beat_features(
                            ppg, bb, be, segment.sleep_stages,
                            r_idx, ppg_fs, ecg_fs, dn_ratio);

                        if (i < static_cast<int>(result.sqi.size())) {
                            auto& s = result.sqi[i];
                            feat.sqi = { s.mean_corr, s.corr_direct, s.corr_interp,
                                         s.corr_dtw, s.frechet };
                        }
                        result.beats[k] = std::move(feat);
                    }
                }));
        }

        for (auto& f : futures) f.get();

        return result;
    }

    /**
     * @brief Run the full feature-extraction pipeline across all segments.
     */
    inline FlattenedBeats generate_features(
        const std::vector<AnnealedSegment>& segments,
        const std::vector<WaveData>& wave_data,
        const std::vector<TemplateInfo>& template_info = {},
        int window_length = 30,
        double sqi_threshold = kInf,
        unsigned num_threads = 0)
    {
        const int num_bins = static_cast<int>(segments.size());
        bool has_templates = !template_info.empty();

        if (static_cast<int>(wave_data.size()) != num_bins)
            throw std::runtime_error("generate_features: segment/wave_data size mismatch");
        if (has_templates && static_cast<int>(template_info.size()) != num_bins)
            throw std::runtime_error("generate_features: segment/template_info size mismatch");

        // Single shared thread pool for ALL work
        ThreadPool pool(num_threads);

        std::vector<std::vector<BeatFeatures>> bins_beats(num_bins);

        for (int t = 0; t < num_bins; ++t) {
            std::cout << "  Section " << (t + 1) << "/" << num_bins << std::flush;

            TemplateInfo ti;
            if (has_templates) ti = template_info[t];

            auto result = process_bin(
                segments[t], wave_data[t], ti,
                sqi_threshold, window_length, pool);

            bins_beats[t] = std::move(result.beats);
            std::cout << " — " << bins_beats[t].size() << " beats\n";
        }

        std::cout << "Flattening...\n";
        return flatten_beat_idx(bins_beats, segments);
    }

} // namespace ppg