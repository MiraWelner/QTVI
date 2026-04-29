/**
 * @file   AnnealSegments.hpp
 * @brief  Header-only annealing algorithm.
 *
 *         Takes a full-night recording and noise markings, splits into
 *         fixed-length bins, excises noisy regions, and redistributes
 *         leftover fragments to neighbouring bins when possible.
 *
 *         ECG noise is only excluded when ALL 3 ECG channels have
 *         overlapping noise for a given time region.
 *         PPG noise is excluded independently.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */
#pragma once

#include "common.hpp"

namespace anneal {

    // ============================================================================
    // Interval helpers
    // ============================================================================

    /// Merge overlapping or touching intervals in-place (sorts first).
    inline void mergeOverlapping(std::vector<std::pair<double, double>>& iv) {
        if (iv.empty()) return;
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<double, double>> m;
        m.push_back(iv[0]);
        for (size_t i = 1; i < iv.size(); ++i) {
            if (iv[i].first <= m.back().second)
                m.back().second = std::max(m.back().second, iv[i].second);
            else
                m.push_back(iv[i]);
        }
        iv = std::move(m);
    }

    /// Return the intersection of two sorted, non-overlapping interval lists.
    inline std::vector<std::pair<double, double>> intersectIntervals(
        const std::vector<std::pair<double, double>>& a,
        const std::vector<std::pair<double, double>>& b)
    {
        std::vector<std::pair<double, double>> out;
        size_t i = 0, j = 0;
        while (i < a.size() && j < b.size()) {
            double lo = std::max(a[i].first, b[j].first);
            double hi = std::min(a[i].second, b[j].second);
            if (lo < hi) out.push_back({ lo, hi });
            (a[i].second < b[j].second) ? ++i : ++j;
        }
        return out;
    }

    /// Merge overlapping or touching uint64 index pairs.
    inline std::vector<std::pair<uint64_t, uint64_t>> mergeSegments(
        std::vector<std::pair<uint64_t, uint64_t>> segs)
    {
        if (segs.empty()) return {};
        std::sort(segs.begin(), segs.end());
        std::vector<std::pair<uint64_t, uint64_t>> m;
        m.push_back(segs[0]);
        for (size_t i = 1; i < segs.size(); ++i) {
            if (segs[i].first <= m.back().second)
                m.back().second = std::max(m.back().second, segs[i].second);
            else
                m.push_back(segs[i]);
        }
        return m;
    }

    // ============================================================================
    // Noise exclusion logic
    // ============================================================================

    /// Combine ECG (3-channel intersection) and PPG noise into one exclusion list.
    inline std::vector<std::pair<double, double>> computeExclusions(
        const NoiseMarkings& nm)
    {
        auto e1 = nm.ecg1, e2 = nm.ecg2, e3 = nm.ecg3, pp = nm.ppg;
        mergeOverlapping(e1);
        mergeOverlapping(e2);
        mergeOverlapping(e3);
        mergeOverlapping(pp);

        std::vector<std::pair<double, double>> ecg;
        if (!e1.empty() && !e2.empty() && !e3.empty())
            ecg = intersectIntervals(intersectIntervals(e1, e2), e3);

        std::vector<std::pair<double, double>> combined;
        combined.insert(combined.end(), ecg.begin(), ecg.end());
        combined.insert(combined.end(), pp.begin(), pp.end());
        mergeOverlapping(combined);
        return combined;
    }

    // ============================================================================
    // Bin geometry
    // ============================================================================

    /// Compute 1-based bin break positions and total bin count.
    inline BinBreaksResult getBinBreaksAndCount(
        uint64_t total_len, uint64_t bin_size, double sr, double min_mins)
    {
        BinBreaksResult r;
        double remainder_mins = (double)(total_len % bin_size) / sr / 60.0;
        r.bin_count = (remainder_mins < min_mins)
            ? (int)std::floor((double)total_len / bin_size)
            : (int)std::ceil((double)total_len / bin_size);

        // Stride is (bin_size + 1), not bin_size, to match MATLAB's
        // annealer. MATLAB places bin n+1 starting one sample after
        // bin n's last sample, so consecutive bins do not share their
        // boundary sample. With stride bin_size, C++ would have bins
        // overlap by 1 sample on the boundary; the old code masked this
        // with the per-bin --second fixup at step 11. With that fixup
        // disabled (to keep the inclusive [first,second] range at
        // bin_size+1 samples like MATLAB), the only way to also match
        // MATLAB's bin start times is to advance breaks by bin_size+1
        // each iteration.
        for (uint64_t b = bin_size + 1; b <= total_len; b += bin_size + 1)
            r.bin_breaks.push_back(b);

        if ((int)r.bin_breaks.size() < r.bin_count)
            r.bin_breaks.push_back(total_len);
        else if (!r.bin_breaks.empty())
            r.bin_breaks.back() = total_len;

        return r;
    }

    /// For each index, return the 1-based bin number it falls in.
    inline std::vector<int> roundToClosestBin(
        const std::vector<uint64_t>& breaks,
        const std::vector<uint64_t>& indices)
    {
        std::vector<int> out(indices.size());
        for (size_t i = 0; i < indices.size(); ++i) {
            out[i] = (int)breaks.size(); // default: last bin
            for (size_t b = 0; b < breaks.size(); ++b) {
                if (indices[i] <= breaks[b]) { out[i] = (int)b + 1; break; }
            }
        }
        return out;
    }

    // ============================================================================
    // Exclusion splitting
    // ============================================================================

    /// Split any exclusion spanning multiple bins so each piece sits in one bin.
    /// MATLAB quirk: only iterates over the original size; appended entries
    /// are not revisited.
    inline void splitOverlappingBins(
        std::vector<Exclusion>& ex,
        const std::vector<uint64_t>& breaks)
    {
        size_t orig = ex.size();
        for (size_t i = 0; i < orig; ++i) {
            if (ex[i].bin_start == ex[i].bin_end) continue;

            uint64_t saved_end = ex[i].idx_end;
            int saved_bin_end = ex[i].bin_end;

            for (int b = ex[i].bin_start; b <= saved_bin_end; ++b) {
                if (b == ex[i].bin_start) {
                    ex[i].idx_end = breaks[b - 1];
                    ex[i].bin_end = b;
                }
                else if (b == saved_bin_end) {
                    ex.push_back({ breaks[b - 2], saved_end, b, b });
                }
                else {
                    ex.push_back({ breaks[b - 2], breaks[b - 1], b, b });
                }
            }
        }
        std::sort(ex.begin(), ex.end(),
            [](const Exclusion& a, const Exclusion& b) {
                return a.idx_start < b.idx_start;
            });
    }

    // ============================================================================
    // Good-interval extraction
    // ============================================================================

    /// Given interval [a, b] and a set of exclusion pairs within it,
    /// return the complementary "good" intervals.
    inline std::vector<std::pair<uint64_t, uint64_t>> getGoodIntervals(
        uint64_t a, uint64_t b,
        const std::vector<std::pair<uint64_t, uint64_t>>& breaks)
    {
        std::vector<uint64_t> flat;
        flat.push_back(a);
        for (const auto& br : breaks) {
            flat.push_back(br.first);
            flat.push_back(br.second);
        }
        flat.push_back(b);

        std::vector<std::pair<uint64_t, uint64_t>> good;
        for (size_t j = 0; j + 1 < flat.size(); j += 2) {
            if (flat[j + 1] > flat[j])
                good.push_back({ flat[j], flat[j + 1] });
        }
        return good;
    }

    // ============================================================================
    // Movement marking
    // ============================================================================

    /// For each bin containing exclusions, determine which good fragments
    /// should be moved left, right, or kept in place.
    inline std::vector<Section> markForMovement(
        const std::vector<Exclusion>& exclusions,
        const std::vector<uint64_t>& breaks,
        int bin_count, double sr, uint64_t bin_size, double min_mins)
    {
        std::set<int> affected;
        for (const auto& ex : exclusions) affected.insert(ex.bin_start);

        std::vector<Section> marked;

        for (int cur : affected) {
            uint64_t bin_begin = breaks[cur - 1] - bin_size;
            uint64_t bin_end = breaks[cur - 1];
            double   bin_half = (double)(bin_end - (bin_size / 2));

            // Collect exclusions belonging to this bin
            std::vector<std::pair<uint64_t, uint64_t>> bin_excl;
            for (const auto& ex : exclusions)
                if (ex.bin_start == cur)
                    bin_excl.push_back({ ex.idx_start, ex.idx_end });

            auto good = getGoodIntervals(bin_begin, bin_end, bin_excl);

            for (const auto& g : good) {
                double m = std::max((double)g.first - bin_half,
                    (double)g.second - bin_half);
                double dur_mins = ((double)(g.second - g.first) / sr) / 60.0;
                bool too_small = dur_mins < min_mins;

                int dir;
                if (cur == 1)         dir = (m <= 0) ? 2 : 0;
                else if (cur == bin_count)  dir = (m >= 0) ? 1 : 0;
                else                        dir = (m <= 0) ? 1 : 0;

                int flag;
                if (!too_small) { dir = 0; flag = 0; }
                else { flag = 1; }

                // Filter out unmovable edge fragments
                if (cur == 1 && m <= 0 && flag == 1) continue;
                if (cur == bin_count && m >= 0 && flag == 1) continue;

                marked.push_back({ g.first, g.second, dir, flag });
            }
        }
        return marked;
    }

    // ============================================================================
    // Section merging
    // ============================================================================

    /// Merge adjacent flagged sections that share a boundary.
    /// MATLAB quirk: index always advances after a merge, skipping the next pair.
    inline void mergeAdjacentSections(
        std::vector<Section>& secs, double sr, double min_mins)
    {
        size_t i = 0;
        while (i + 1 < secs.size()) {
            if (secs[i].end == secs[i + 1].begin &&
                secs[i].flag != 0 && secs[i + 1].flag != 0)
            {
                double t1 = ((double)(secs[i].end - secs[i].begin) / sr) / 60.0;
                double t2 = ((double)(secs[i + 1].end - secs[i + 1].begin) / sr) / 60.0;

                if (t1 + t2 >= min_mins) {
                    secs[i].dir = 0;
                    secs[i].flag = 0;
                }
                else {
                    int idx = (t1 >= t2) ? 0 : 1;
                    secs[i].dir = secs[i + idx].dir;
                    secs[i].flag = 1;
                }
                secs[i].end = secs[i + 1].end;
                secs.erase(secs.begin() + (int64_t)i + 1);
            }
            ++i; // always advance (MATLAB behaviour)
        }
    }

} // namespace anneal


// ============================================================================
// Main entry point
// ============================================================================

inline std::vector<FinalSegment> AnnealSegments(
    const RawData& data,
    const NoiseMarkings& noiseMarkings,
    double targetLenMins)
{
    using namespace anneal;

    const bool hasPpg = data.ppg.size() > 1;
    const bool hasEcg = data.ecg1.size() > 1;
    if (!hasPpg && !hasEcg) return {};

    const bool ecgOnly = !hasPpg && hasEcg;
    const auto& primarySignal = ecgOnly ? data.ecg1 : data.ppg;
    const double primarySR = ecgOnly ? data.ecgSR : data.ppgSR;

    const uint64_t bin_size = static_cast<uint64_t>(primarySR * 60.0 * targetLenMins);
    const double   min_mins = targetLenMins / 2.0;
    const double   min_excl_s = 5.0;
    const uint64_t total_len = primarySignal.size();

    // 1. Bin geometry
    auto bbc = getBinBreaksAndCount(total_len, bin_size, primarySR, min_mins);
    if (bbc.bin_count <= 0) return {};
    const int bin_count = bbc.bin_count;
    const auto& breaks = bbc.bin_breaks;

    // 2. Compute combined noise exclusions
    auto combined = computeExclusions(noiseMarkings);

    // 3. Filter short exclusions
    std::vector<std::pair<double, double>> excl_sec;
    for (const auto& seg : combined)
        if ((seg.second - seg.first) >= min_excl_s)
            excl_sec.push_back(seg);

    // 4. Convert exclusion times to sample indices
    std::vector<uint64_t> excl_flat;
    for (const auto& seg : excl_sec) {
        excl_flat.push_back(closest_idx(seg.first, primarySR));
        excl_flat.push_back(closest_idx(seg.second, primarySR));
    }

    // 5. Map exclusion boundaries to bins
    auto excl_bins = roundToClosestBin(breaks, excl_flat);

    std::vector<Exclusion> exclusions;
    for (size_t i = 0; i < excl_sec.size(); ++i) {
        exclusions.push_back({
            excl_flat[i * 2], excl_flat[i * 2 + 1],
            excl_bins[i * 2], excl_bins[i * 2 + 1]
            });
    }

    // 6. Split multi-bin exclusions
    splitOverlappingBins(exclusions, breaks);

    // 7. Build good sections: unaffected bins + fragments from affected bins
    std::set<int> affected;
    for (const auto& ex : exclusions) affected.insert(ex.bin_start);

    std::vector<Section> good;
    for (int b = 1; b <= bin_count; ++b) {
        if (affected.count(b) == 0)
            good.push_back({ breaks[b - 1] - bin_size, breaks[b - 1], 0, 0 });
    }

    auto marked = markForMovement(exclusions, breaks, bin_count,
        primarySR, bin_size, min_mins);
    good.insert(good.end(), marked.begin(), marked.end());

    std::sort(good.begin(), good.end(),
        [](const Section& a, const Section& b) { return a.begin < b.begin; });

    // 8. Merge adjacent flagged sections
    mergeAdjacentSections(good, primarySR, min_mins);

    // 9. Assign sections to final bins
    struct BinIdx { std::vector<std::pair<uint64_t, uint64_t>> po; };
    std::vector<BinIdx> final_idx(1);

    int cur = 0;
    std::vector<std::pair<uint64_t, uint64_t>> temp;

    for (const auto& sec : good) {
        if (sec.flag) {
            if (sec.dir == 1) {
                auto& target = (cur > 0) ? final_idx[cur - 1] : final_idx[cur];
                target.po.push_back({ sec.begin, sec.end });
            }
            else {
                temp.push_back({ sec.begin, sec.end });
            }
        }
        else {
            for (const auto& t : temp)
                final_idx[cur].po.push_back(t);
            final_idx[cur].po.push_back({ sec.begin, sec.end });
            temp.clear();
            cur++;
            final_idx.push_back(BinIdx());
        }
    }

    while (!final_idx.empty() && final_idx.back().po.empty())
        final_idx.pop_back();

    // 10. Merge overlapping segments within each bin
    for (auto& fb : final_idx)
        if (fb.po.size() > 1)
            fb.po = mergeSegments(fb.po);

    // 11. Fix shared boundaries between adjacent bins
    //
    // DISABLED to match the MATLAB annealer, which does NOT apply this
    // fixup. With it enabled, every bin except the last got its end
    // sample decremented, producing 60000-sample bins at 1000 Hz. The
    // MATLAB annealer leaves the inclusive [first, second] range alone,
    // so its bins are 60001 samples. Disabling the fixup here makes the
    // two pipelines produce the same bin lengths and identical R-peak
    // detection times.
    // for (size_t i = 0; i + 1 < final_idx.size(); ++i) {
    //     if (!final_idx[i].po.empty() && !final_idx[i + 1].po.empty()) {
    //         if (final_idx[i].po.back().second == final_idx[i + 1].po.front().first)
    //             final_idx[i].po.back().second--;
    //     }
    // }

    // 12. Compute secondary (cross-modality) index pairs
    const double secondarySR = ecgOnly ? data.ppgSR : data.ecgSR;
    const bool hasSecondary = ecgOnly ? hasPpg : hasEcg;

    struct BinIdxFull {
        std::vector<std::pair<uint64_t, uint64_t>> primary;
        std::vector<std::pair<uint64_t, uint64_t>> secondary;
    };
    std::vector<BinIdxFull> final_bins(final_idx.size());

    for (size_t i = 0; i < final_idx.size(); ++i) {
        final_bins[i].primary = final_idx[i].po;
        if (hasSecondary && secondarySR > 0) {
            for (const auto& seg : final_idx[i].po) {
                double t0 = (double)(seg.first - 1) / primarySR;
                double t1 = (double)(seg.second - 1) / primarySR;
                final_bins[i].secondary.push_back({
                    closest_idx(t0, secondarySR),
                    closest_idx(t1, secondarySR)
                    });
            }
        }
    }

    // 13. Extract signal samples and build output segments
    std::vector<double> sleep_times(data.sleepStages.size());
    for (size_t s = 0; s < data.sleepStages.size(); ++s)
        sleep_times[s] = (s + 1) * data.scoringEpochSec;

    auto extract = [](const std::vector<double>& sig,
        const std::vector<std::pair<uint64_t, uint64_t>>& idx,
        std::vector<double>& out) {
            for (const auto& seg : idx)
                for (uint64_t k = seg.first; k <= seg.second && k <= sig.size(); ++k)
                    out.push_back(sig[k - 1]);
        };

    std::vector<FinalSegment> results(final_bins.size());
    for (size_t i = 0; i < final_bins.size(); ++i) {
        auto& r = results[i];
        r.ppgSampleRate = data.ppgSR;
        r.ecgSampleRate = data.ecgSR;
        r.scoring_epoch_size_sec = data.scoringEpochSec;

        if (ecgOnly) {
            r.ecg_bin_indexs = final_bins[i].primary;
            r.ppg_bin_indexs = final_bins[i].secondary;
        }
        else {
            r.ppg_bin_indexs = final_bins[i].primary;
            r.ecg_bin_indexs = final_bins[i].secondary;
        }

        const auto& ecgIdx = ecgOnly ? final_bins[i].primary : final_bins[i].secondary;
        const auto& ppgIdx = ecgOnly ? final_bins[i].secondary : final_bins[i].primary;

        if (!data.ppg.empty() && !ppgIdx.empty()) extract(data.ppg, ppgIdx, r.ppg);
        if (!data.ecg1.empty() && !ecgIdx.empty()) extract(data.ecg1, ecgIdx, r.ecg1);
        if (!data.ecg2.empty() && !ecgIdx.empty()) extract(data.ecg2, ecgIdx, r.ecg2);
        if (!data.ecg3.empty() && !ecgIdx.empty()) extract(data.ecg3, ecgIdx, r.ecg3);

        // Sleep stages falling within this bin's time range
        for (const auto& seg : final_bins[i].primary) {
            double t0 = (double)(seg.first - 1) / primarySR;
            double t1 = (double)(seg.second - 1) / primarySR;
            for (size_t s = 0; s < data.sleepStages.size(); ++s)
                if (sleep_times[s] >= t0 && sleep_times[s] <= t1)
                    r.sleep_stages.push_back(data.sleepStages[s]);
        }
    }

    return results;
}