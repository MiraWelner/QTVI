/**
 * @file   AnnealSegments.hpp
 * @brief  The core annealing algorithm (header-only).
 *
 * ECG noise is only excluded when ALL 3 ECG channels have overlapping
 * noise for a given time region. PPG noise is excluded independently.
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-03-22
 */
#pragma once

#include "common.hpp"
#include "getBinBreaksAndCount.hpp"
#include "RoundToClosestBin.hpp"
#include "splitOverlappingBins.hpp"
#include "markForMovement.hpp"
#include "mergeIntervals.hpp"
#include "MergeSegments.hpp"

namespace detail {

    inline void mergeOverlapping(std::vector<std::pair<double, double>>& intervals) {
        if (intervals.empty()) return;
        std::sort(intervals.begin(), intervals.end());
        std::vector<std::pair<double, double>> merged;
        merged.push_back(intervals[0]);
        for (size_t i = 1; i < intervals.size(); ++i) {
            if (intervals[i].first <= merged.back().second)
                merged.back().second = std::max(merged.back().second, intervals[i].second);
            else
                merged.push_back(intervals[i]);
        }
        intervals = std::move(merged);
    }

    inline std::vector<std::pair<double, double>> intersectIntervals(
        const std::vector<std::pair<double, double>>& a,
        const std::vector<std::pair<double, double>>& b)
    {
        std::vector<std::pair<double, double>> result;
        size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            double lo = std::max(a[i].first, b[j].first);
            double hi = std::min(a[i].second, b[j].second);
            if (lo < hi)
                result.push_back({ lo, hi });
            if (a[i].second < b[j].second)
                ++i;
            else
                ++j;
        }
        return result;
    }

    inline std::vector<std::pair<double, double>> computeExclusions(const NoiseMarkings& markings) {
        auto ecg1 = markings.ecg1;
        auto ecg2 = markings.ecg2;
        auto ecg3 = markings.ecg3;
        auto ppg = markings.ppg;

        mergeOverlapping(ecg1);
        mergeOverlapping(ecg2);
        mergeOverlapping(ecg3);
        mergeOverlapping(ppg);

        std::vector<std::pair<double, double>> ecgExclusion;
        if (!ecg1.empty() && !ecg2.empty() && !ecg3.empty()) {
            auto ecg12 = intersectIntervals(ecg1, ecg2);
            ecgExclusion = intersectIntervals(ecg12, ecg3);
        }

        std::vector<std::pair<double, double>> combined;
        combined.insert(combined.end(), ecgExclusion.begin(), ecgExclusion.end());
        combined.insert(combined.end(), ppg.begin(), ppg.end());
        mergeOverlapping(combined);

        return combined;
    }

} // namespace detail


inline std::vector<FinalSegment> AnnealSegments(
    const RawData& data,
    const NoiseMarkings& noiseMarkings,
    double targetLenMins)
{
    const bool hasPpg = data.ppg.size() > 1;
    const bool hasEcg = data.ecg1.size() > 1;

    if (!hasPpg && !hasEcg) return {};

    const bool ecgOnly = !hasPpg && hasEcg;
    const std::vector<double>& primarySignal = ecgOnly ? data.ecg1 : data.ppg;
    const double primarySR = ecgOnly ? data.ecgSR : data.ppgSR;

    const uint64_t bin_size_samples = static_cast<uint64_t>(primarySR * 60.0 * targetLenMins);
    const double min_bin_size_mins = targetLenMins / 2.0;
    const double min_exclusion_bin_size_seconds = 5.0;
    const uint64_t total_len = primarySignal.size();

    // 1. Bin breaks and count
    auto bbc = getBinBreaksAndCount(total_len, bin_size_samples, primarySR, min_bin_size_mins);
    if (bbc.bin_count <= 0) return {};
    const int bin_count = bbc.bin_count;
    const auto& bin_breaks = bbc.bin_breaks;

    // 2. Compute combined exclusions
    auto combinedExclusions = detail::computeExclusions(noiseMarkings);

    // 3. Filter exclusions by minimum length
    std::vector<std::pair<double, double>> exclusions_seconds;
    for (const auto& seg : combinedExclusions) {
        if ((seg.second - seg.first) >= min_exclusion_bin_size_seconds)
            exclusions_seconds.push_back(seg);
    }

    // 4. Convert exclusion times to primary sample indices (1-based)
    std::vector<uint64_t> exclusions_indexs_flat;
    for (const auto& seg : exclusions_seconds) {
        exclusions_indexs_flat.push_back(closest_idx(seg.first, primarySR));
        exclusions_indexs_flat.push_back(closest_idx(seg.second, primarySR));
    }

    // 5. Determine which bin each exclusion boundary falls in
    std::vector<int> exclusions_bin = RoundToClosestBin(bin_breaks, exclusions_indexs_flat);

    std::vector<Exclusion> exclusions;
    for (size_t i = 0; i < exclusions_seconds.size(); ++i) {
        exclusions.push_back({
            exclusions_indexs_flat[i * 2],
            exclusions_indexs_flat[i * 2 + 1],
            exclusions_bin[i * 2],
            exclusions_bin[i * 2 + 1]
            });
    }

    // 6. Split exclusions that span multiple bins
    splitOverlappingBins(exclusions, bin_breaks);

    // 7. Build good_sections
    std::set<int> update_bins_set;
    for (const auto& ex : exclusions) update_bins_set.insert(ex.bin_start);

    std::vector<Section> good_sections;

    for (int b = 1; b <= bin_count; ++b) {
        if (update_bins_set.find(b) == update_bins_set.end()) {
            uint64_t bin_begin = bin_breaks[b - 1] - bin_size_samples;
            uint64_t bin_end = bin_breaks[b - 1];
            good_sections.push_back({ bin_begin, bin_end, 0, 0 });
        }
    }

    auto marked = markForMovement(exclusions, bin_breaks, bin_count,
        primarySR, bin_size_samples, min_bin_size_mins);
    good_sections.insert(good_sections.end(), marked.begin(), marked.end());

    std::sort(good_sections.begin(), good_sections.end(),
        [](const Section& a, const Section& b) { return a.begin < b.begin; });

    // 8. Merge adjacent moving sections
    mergeIntervals(good_sections, primarySR, min_bin_size_mins);

    // 9. Assign sections to final bins
    struct BinIdx {
        std::vector<std::pair<uint64_t, uint64_t>> po;
    };
    std::vector<BinIdx> final_bin_idx;
    final_bin_idx.push_back(BinIdx());

    int current_bin = 0;
    std::vector<std::pair<uint64_t, uint64_t>> temp_bin;

    for (const auto& sec : good_sections) {
        if (sec.flag) {
            if (sec.dir == 1) {
                if (current_bin > 0)
                    final_bin_idx[current_bin - 1].po.push_back({ sec.begin, sec.end });
                else
                    final_bin_idx[current_bin].po.push_back({ sec.begin, sec.end });
            }
            else {
                temp_bin.push_back({ sec.begin, sec.end });
            }
        }
        else {
            for (const auto& tb : temp_bin)
                final_bin_idx[current_bin].po.push_back(tb);
            final_bin_idx[current_bin].po.push_back({ sec.begin, sec.end });
            temp_bin.clear();
            current_bin++;
            final_bin_idx.push_back(BinIdx());
        }
    }

    while (!final_bin_idx.empty() && final_bin_idx.back().po.empty())
        final_bin_idx.pop_back();

    // 10. Merge overlapping segments within each bin
    for (auto& fb : final_bin_idx) {
        if (fb.po.size() > 1)
            fb.po = MergeSegments(fb.po);
    }

    // 11. Correct overlaps between adjacent bins
    for (size_t i = 0; i + 1 < final_bin_idx.size(); ++i) {
        if (!final_bin_idx[i].po.empty() && !final_bin_idx[i + 1].po.empty()) {
            if (final_bin_idx[i].po.back().second == final_bin_idx[i + 1].po.front().first)
                final_bin_idx[i].po.back().second--;
        }
    }

    // 12. Compute secondary index pairs
    const double secondarySR = ecgOnly ? data.ppgSR : data.ecgSR;
    const bool hasSecondary = ecgOnly ? hasPpg : hasEcg;

    struct BinIdxFull {
        std::vector<std::pair<uint64_t, uint64_t>> primary;
        std::vector<std::pair<uint64_t, uint64_t>> secondary;
    };
    std::vector<BinIdxFull> final_bins(final_bin_idx.size());

    for (size_t i = 0; i < final_bin_idx.size(); ++i) {
        final_bins[i].primary = final_bin_idx[i].po;
        if (hasSecondary && secondarySR > 0) {
            for (const auto& seg : final_bin_idx[i].po) {
                double time_start = (double)(seg.first - 1) / primarySR;
                double time_end = (double)(seg.second - 1) / primarySR;
                uint64_t sec_start = closest_idx(time_start, secondarySR);
                uint64_t sec_end = closest_idx(time_end, secondarySR);
                final_bins[i].secondary.push_back({ sec_start, sec_end });
            }
        }
    }

    // 13. Collect data from indices
    std::vector<double> sleep_stage_times(data.sleepStages.size());
    for (size_t s = 0; s < data.sleepStages.size(); ++s)
        sleep_stage_times[s] = (s + 1) * data.scoringEpochSec;

    auto extractSamples = [](const std::vector<double>& signal,
        const std::vector<std::pair<uint64_t, uint64_t>>& idxPairs,
        std::vector<double>& out) {
            for (const auto& seg : idxPairs) {
                for (uint64_t k = seg.first; k <= seg.second && k <= signal.size(); ++k)
                    out.push_back(signal[k - 1]);
            }
        };

    std::vector<FinalSegment> results(final_bins.size());
    for (size_t i = 0; i < final_bins.size(); ++i) {
        results[i].ppgSampleRate = data.ppgSR;
        results[i].ecgSampleRate = data.ecgSR;
        results[i].scoring_epoch_size_sec = data.scoringEpochSec;

        if (ecgOnly) {
            results[i].ecg_bin_indexs = final_bins[i].primary;
            results[i].ppg_bin_indexs = final_bins[i].secondary;
        }
        else {
            results[i].ppg_bin_indexs = final_bins[i].primary;
            results[i].ecg_bin_indexs = final_bins[i].secondary;
        }

        const auto& ecgIdxPairs = ecgOnly ? final_bins[i].primary : final_bins[i].secondary;
        const auto& ppgIdxPairs = ecgOnly ? final_bins[i].secondary : final_bins[i].primary;

        if (!data.ppg.empty() && !ppgIdxPairs.empty())
            extractSamples(data.ppg, ppgIdxPairs, results[i].ppg);

        if (!data.ecg1.empty() && !ecgIdxPairs.empty())
            extractSamples(data.ecg1, ecgIdxPairs, results[i].ecg1);
        if (!data.ecg2.empty() && !ecgIdxPairs.empty())
            extractSamples(data.ecg2, ecgIdxPairs, results[i].ecg2);
        if (!data.ecg3.empty() && !ecgIdxPairs.empty())
            extractSamples(data.ecg3, ecgIdxPairs, results[i].ecg3);

        // Sleep stages
        for (size_t w = 0; w < final_bins[i].primary.size(); ++w) {
            uint64_t pri_start = final_bins[i].primary[w].first;
            uint64_t pri_end = final_bins[i].primary[w].second;
            double time_first = (double)(pri_start - 1) / primarySR;
            double time_last = (double)(pri_end - 1) / primarySR;
            for (size_t s = 0; s < data.sleepStages.size(); ++s) {
                if (sleep_stage_times[s] >= time_first && sleep_stage_times[s] <= time_last)
                    results[i].sleep_stages.push_back(data.sleepStages[s]);
            }
        }
    }

    return results;
}