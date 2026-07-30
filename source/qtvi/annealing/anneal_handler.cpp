/**
 * @file   anneal_handler.cpp
 * @brief  Step 3 of the QTVi pipeline. Reads a v2 data .bin (from
 *         file_to_bin) plus an optional noise-markings .bin (from the
 *         marking GUI), runs the annealing algorithm, and writes an
 *         annealed .bin with all 40 channels preserved.
 *
 *         File-private structure:
 *           - Data types (RawData, FinalSegment, NoiseMarkings, ...)
 *           - Algorithm helpers (in anneal:: namespace)
 *           - AnnealSegments() entry point for the algorithm
 *           - File I/O (read_noise_bin, read_data_bin, write_output_bin)
 *           - annealOneFile() public dispatch
 *
 *         Algorithm summary: split the recording into fixed-length bins,
 *         excise noisy regions, redistribute leftover fragments to
 *         neighbouring bins when possible.
 *
 *         ECG noise is excluded only when EVERY ECG channel that the
 *         user marked agrees on the noisy region. Channels with no
 *         marks at all do not constrain (a missing channel does not
 *         mean "always clean"). Single-lead datasets like MESA where
 *         the user only marks ECG1 still get their noise excluded.
 *         PPG noise is excluded independently.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-22
 */

#include "anneal_handler.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

    // ============================================================================
    // Algorithm data structures
    // ============================================================================

    struct RawData {
        std::vector<double> ppg, ecg1, ecg2, ecg3, sleepStages;
        double ppgSR = 0, ecgSR = 0, scoringEpochSec = 0;
    };

    struct FinalSegment {
        std::vector<double> ppg;
        std::vector<double> ecg1, ecg2, ecg3;
        std::vector<double> sleep_stages;
        std::vector<std::pair<uint64_t, uint64_t>> ppg_bin_indexs;
        std::vector<std::pair<uint64_t, uint64_t>> ecg_bin_indexs;
        double ppgSampleRate = 0;
        double ecgSampleRate = 0;
        double scoring_epoch_size_sec = 0;
    };

    /// Per-channel noise markings as time intervals (seconds).
    /// ECG noise is only excluded when all 3 channels overlap.
    /// PPG noise is excluded independently.
    struct NoiseMarkings {
        std::vector<std::pair<double, double>> ecg1, ecg2, ecg3, ppg;
        std::vector<std::pair<double, double>> accel, abp, art, art_pulm;
    };

    struct Exclusion {
        uint64_t idx_start, idx_end;
        int bin_start, bin_end;
    };

    /// A contiguous good section of signal.
    ///   dir:  0 = stationary, 1 = move left, 2 = move right
    ///   flag: 1 = too small, needs merging with a neighbour
    struct Section {
        uint64_t begin, end;
        int dir = 0;
        int flag = 0;
    };

    struct BinBreaksResult {
        std::vector<uint64_t> bin_breaks;
        int bin_count = 0;
    };

    // ============================================================================
    // MATLAB-compatible rounding
    // ============================================================================

    /// Banker's rounding (round half to even) to match MATLAB's round().
    double matlab_round(double x) {
        double r = std::round(x);
        double frac = x - std::floor(x);
        if (std::abs(frac - 0.5) < 1e-12) {
            double f = std::floor(x);
            r = (std::fmod(std::abs(f), 2.0) < 0.5) ? f : f + 1.0;
        }
        return r;
    }

    /// Convert a time in seconds to a 1-based sample index.
    uint64_t closest_idx(double time_sec, double sr) {
        double rounded = matlab_round(time_sec * sr);
        if (rounded < 0.0) rounded = 0.0;
        return static_cast<uint64_t>(rounded) + 1;
    }

    // ============================================================================
    // Algorithm helpers
    // ============================================================================

    namespace anneal {

        // ------------------------------------------------------------------------
        // Interval helpers
        // ------------------------------------------------------------------------

        /// Merge overlapping or touching intervals in-place (sorts first).
        void mergeOverlapping(std::vector<std::pair<double, double>>& iv) {
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
        std::vector<std::pair<double, double>> intersectIntervals(
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
        std::vector<std::pair<uint64_t, uint64_t>> mergeSegments(
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


        std::vector<std::pair<double, double>> computeExclusions(const NoiseMarkings& nm)
        {

            // Intersect over the channels the user actually marked. An empty list =
            // "not present / not reviewed" and does NOT constrain -- so datasets that
            // don't have (or don't mark) ABP/ART/ART_PULM, i.e. MESA and BITTIUM, are
            // completely unaffected. A region is excised only where every MARKED
            // channel agrees it is noisy.
            auto e1 = nm.ecg1, e2 = nm.ecg2, e3 = nm.ecg3, pp = nm.ppg;
            auto ab = nm.abp, ar = nm.art, ap = nm.art_pulm;
            for (auto* v : { &e1, &e2, &e3, &pp, &ab, &ar, &ap })
                mergeOverlapping(*v);
            std::vector<std::vector<std::pair<double, double>>*> present;
            for (auto* v : { &e1, &e2, &e3, &pp, &ab, &ar, &ap })
                if (!v->empty()) present.push_back(v);

            if (present.empty()) return {};

            std::vector<std::pair<double, double>> common = *present[0];
            for (size_t i = 1; i < present.size(); ++i)
                common = intersectIntervals(common, *present[i]);
            mergeOverlapping(common);
            return common;
        }

        // ------------------------------------------------------------------------
        // Bin geometry
        // ------------------------------------------------------------------------

        // Compute 1-based bin break positions and total bin count.
        BinBreaksResult getBinBreaksAndCount(uint64_t total_len, uint64_t bin_size, double sr, double min_mins)
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
        std::vector<int> roundToClosestBin(
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

        // ------------------------------------------------------------------------
        // Exclusion splitting
        // ------------------------------------------------------------------------

        /// Split any exclusion spanning multiple bins so each piece sits in one bin.
        /// MATLAB quirk: only iterates over the original size; appended entries
        /// are not revisited.
        void splitOverlappingBins(
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

        // ------------------------------------------------------------------------
        // Good-interval extraction
        // ------------------------------------------------------------------------

        /// Given interval [a, b] and a set of exclusion pairs within it,
        /// return the complementary "good" intervals.
        std::vector<std::pair<uint64_t, uint64_t>> getGoodIntervals(
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

        // ------------------------------------------------------------------------
        // Movement marking
        // ------------------------------------------------------------------------

        /// For each bin containing exclusions, determine which good fragments
        /// should be moved left, right, or kept in place.
        std::vector<Section> markForMovement(
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

        // ------------------------------------------------------------------------
        // Section merging
        // ------------------------------------------------------------------------

        /// Merge adjacent flagged sections that share a boundary.
        /// MATLAB quirk: index always advances after a merge, skipping the next pair.
        void mergeAdjacentSections(
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
    // Algorithm entry point
    // ============================================================================

    std::vector<FinalSegment> AnnealSegments(
        const RawData& data,
        const NoiseMarkings& noiseMarkings,
        double targetLenMins)
    {
        using namespace anneal;

        const bool hasPpg = data.ppg.size() > 1 && data.ppgSR > 0.0;
        const bool hasEcg = data.ecg1.size() > 1 && data.ecgSR > 0.0;
        if (!hasPpg && !hasEcg) {
            std::cerr << "  AnnealSegments empty: no usable channel (ppg n="
                << data.ppg.size() << " sr=" << data.ppgSR
                << ", ecg n=" << data.ecg1.size() << " sr=" << data.ecgSR << ")\n";
            return {};
        }

        const bool ecgOnly = !hasPpg && hasEcg;
        const auto& primarySignal = ecgOnly ? data.ecg1 : data.ppg;
        const double primarySR = ecgOnly ? data.ecgSR : data.ppgSR;

        const uint64_t bin_size = static_cast<uint64_t>(primarySR * 60.0 * targetLenMins);
        const double   min_mins = targetLenMins / 2.0;
        const double   min_excl_s = 5.0;
        const uint64_t total_len = primarySignal.size();

        // 1. Bin geometry
        auto bbc = getBinBreaksAndCount(total_len, bin_size, primarySR, min_mins);
        if (bbc.bin_count <= 0) {
            std::cerr << "  AnnealSegments empty: zero bins (total_len=" << total_len
                << ", bin_size=" << bin_size << ", primarySR=" << primarySR
                << ", min_mins=" << min_mins << ")\n";
            return {};
        }
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
                // Reserve once for all segments so multi-segment bins don't
                // repeatedly reallocate/copy as `out` grows.
                size_t total = 0;
                for (const auto& seg : idx) {
                    if (seg.first > sig.size()) continue;
                    const uint64_t last = std::min<uint64_t>(seg.second, sig.size());
                    if (last >= seg.first) total += static_cast<size_t>(last - seg.first + 1);
                }
                out.reserve(out.size() + total);

                // Bulk range-copy per segment instead of one push_back per
                // sample -- same clamping behaviour as the original
                // (k <= sig.size()) bound, just done as one contiguous copy.
                for (const auto& seg : idx) {
                    if (seg.first > sig.size()) continue;
                    const uint64_t last = std::min<uint64_t>(seg.second, sig.size());
                    if (last < seg.first) continue;
                    out.insert(out.end(),
                        sig.begin() + (seg.first - 1),
                        sig.begin() + last);
                }
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

    // ============================================================================
    // File I/O
    // ============================================================================

    constexpr int NUM_CHANNELS = 36;

    struct Extras {
        // NB: parenthesis-init via assignment to invoke the count constructor.
        // Brace-init {NUM_CHANNELS} on a vector<vector<double>> hits the
        // initializer_list ctor, producing a 1-element outer vector instead
        // of a NUM_CHANNELS-element one -- which silently corrupts memory the
        // moment we index past slot 0.
        std::vector<std::vector<double>> upsampled =
            std::vector<std::vector<double>>(NUM_CHANNELS);
        std::vector<std::vector<double>> rawFlat =
            std::vector<std::vector<double>>(NUM_CHANNELS);
        std::vector<float> nativeRates = std::vector<float>(NUM_CHANNELS, 0.0f);
        uint32_t signal_rate = 0, boolean_rate = 0, pacemaker_rate = 0, sleep_rate = 0;
    };

    // Noise reader
    //
    // Format: [uint64 count] then count x [6 doubles per row]
    //   Row: startSample, endSample, startSec, endSec, labelId, typeId
    //   labelId: 0=unknown, 1=PPG, 2=ECG1, 3=ECG2, 4=ECG3, 5=ABP
    // (typeId is unused here -- we only care about which channel was marked.)
    NoiseMarkings read_noise_bin(const std::filesystem::path& path) {
        NoiseMarkings m;
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) return m;

        uint64_t count = 0;
        f.read(reinterpret_cast<char*>(&count), 8);

        for (uint64_t i = 0; i < count; ++i) {
            double row[6];
            f.read(reinterpret_cast<char*>(row), 48);
            std::pair<double, double> iv = { row[2], row[3] };
            switch (static_cast<int>(row[4])) {
            case 1: m.ppg.push_back(iv);      break;
            case 2: m.ecg1.push_back(iv);     break;
            case 3: m.ecg2.push_back(iv);     break;
            case 4: m.ecg3.push_back(iv);     break;
            case 5: m.abp.push_back(iv);      break;
            case 7: m.art.push_back(iv);      break;
            case 8: m.art_pulm.push_back(iv); break;
                // case 6 (ACCEL) intentionally excluded from annealing
            default: break;
            }
        }
        return m;
    }

    // Header layout MUST match write_header_and_close() in file_to_bin.cpp:
//   [uint32 header_version]
//   [uint32 n_channels]
//   [uint32 sleep_state_len]
//   [NUM_CH x uint32 sizes_up]
//   [NUM_CH x uint32 sizes_raw]
//   [NUM_CH x float32 native_rates]
//   [NUM_CH x float32 up_rates]
//   [uint32 sleep_size]
// NUM_CH = 36 (file_to_bin NUM_CHANNELS). Channel slots per file_to_bin
// ChannelIdx: 0=timestamp, 1=ECG1, 2=ECG2, 3=ECG3, 4=PPG, ...
    void read_data_bin(const std::filesystem::path& path,
        RawData& data, Extras& extras)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f.is_open()) throw std::runtime_error("cannot open: " + path.string());

        char read_buf[1 << 16];
        f.rdbuf()->pubsetbuf(read_buf, sizeof(read_buf));

        f.seekg(0, std::ios::end);
        const uint64_t fileSize = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        constexpr int    NCH = 36;                       // file_to_bin NUM_CHANNELS
        constexpr size_t NHF = 4 + 4 * NCH;              // version+n_channels+sleep_state_len+sleep_size + 4*36 = 148 fields
        constexpr size_t HDR = NHF * 4;                  // = 592 bytes
        if (fileSize < HDR)
            throw std::runtime_error("bin too small (need " + std::to_string(HDR)
                + "-byte header): " + path.string());

        uint32_t header_version = 0, n_channels = 0;
        uint32_t sleep_state_len = 0, sleep_size = 0;
        std::vector<uint32_t> sizes_up(NCH), sizes_raw(NCH);
        std::vector<float>    native_rates(NCH), up_rates(NCH);

        f.read(reinterpret_cast<char*>(&header_version), 4);
        f.read(reinterpret_cast<char*>(&n_channels), 4);
        if (header_version != 1 || n_channels != static_cast<uint32_t>(NCH))
            throw std::runtime_error("unsupported .bin header (version=" +
                std::to_string(header_version) + ", n_channels=" + std::to_string(n_channels) +
                "); regenerate with file_to_bin: " + path.string());
        f.read(reinterpret_cast<char*>(&sleep_state_len), 4);
        f.read(reinterpret_cast<char*>(sizes_up.data()), NCH * 4);
        f.read(reinterpret_cast<char*>(sizes_raw.data()), NCH * 4);
        f.read(reinterpret_cast<char*>(native_rates.data()), NCH * 4);
        f.read(reinterpret_cast<char*>(up_rates.data()), NCH * 4);
        f.read(reinterpret_cast<char*>(&sleep_size), 4);

        // Per-channel upsample rates: the anneal consumes the UPSAMPLED
        // blocks, so use up_rates (not native_rates). ECG1 = slot 1, PPG = slot 4.
        constexpr int SLOT_ECG1 = 1, SLOT_ECG2 = 2, SLOT_ECG3 = 3, SLOT_PPG = 4;
        data.ecgSR = static_cast<double>(up_rates[SLOT_ECG1]);
        data.ppgSR = static_cast<double>(up_rates[SLOT_PPG]);
        data.scoringEpochSec = static_cast<double>(sleep_state_len);

        if (data.ecgSR <= 0.0 && data.ppgSR <= 0.0)
            throw std::runtime_error("input .bin has zero ECG and PPG upsample rate: "
                + path.string());

        // Body: NCH x { upsampled block, raw (t,v) block }, then sleep stages.
        auto safeReadDoubles = [&](std::vector<double>& dest, uint64_t count) {
            const uint64_t pos = static_cast<uint64_t>(f.tellg());
            const uint64_t avail = (fileSize > pos) ? (fileSize - pos) / 8 : 0;
            const uint64_t actual = std::min<uint64_t>(count, avail);
            dest.resize(actual);
            if (actual > 0)
                f.read(reinterpret_cast<char*>(dest.data()),
                    static_cast<std::streamsize>(actual * 8));
            };

        // Keep Extras sized to NCH so the writer half stays consistent.
        extras.upsampled.assign(NCH, {});
        extras.rawFlat.assign(NCH, {});
        extras.nativeRates.assign(NCH, 0.0f);
        for (int i = 0; i < NCH; ++i) {
            safeReadDoubles(extras.upsampled[i], sizes_up[i]);
            safeReadDoubles(extras.rawFlat[i], static_cast<uint64_t>(sizes_raw[i]) * 2);
            extras.nativeRates[i] = native_rates[i];
        }
        safeReadDoubles(data.sleepStages, sleep_size);

        data.ecg1 = extras.upsampled[SLOT_ECG1];
        data.ecg2 = extras.upsampled[SLOT_ECG2];
        data.ecg3 = extras.upsampled[SLOT_ECG3];
        data.ppg = extras.upsampled[SLOT_PPG];
    }

    // Output writer
    //
    // Layout:
    //   Header: [uint64 nSegments][double ppgSR][double ecgSR][double epochSec]
    //           [uint32 nChannels=36][36 x float32 nativeRates]
    //           [uint8 ecg1_inverted][uint8 ecg2_inverted][uint8 ecg3_inverted]
    //   Per segment:
    //     ppg_bin_indexs, ecg_bin_indexs, ppg, ecg1, ecg2, ecg3, sleep,
    //     then 36 x {upsampled_slice, raw_slice}.
    //
    // The 36 trailing per-channel slices preserve every input channel sliced
    // to the segment's time window:
    //   - Upsampled slice: indices proportional to ecg_bin_indexs. ECG-index
    //     i in an N-sample ECG block maps to channel-X-index ceil(i*M/N) in
    //     an M-sample block of channel X. Avoids per-channel rate bookkeeping.
    //   - Raw slice: filter (t,v) pairs whose t falls inside the ECG time
    //     window. Timestamps stay in absolute seconds-from-recording-start.
    void write_output_bin(const std::filesystem::path& path,
        const std::vector<FinalSegment>& segs,
        const Extras& extras,
        bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted)
    {
        std::ofstream out(path, std::ios::binary);
        if (!out.is_open())
            throw std::runtime_error("cannot create: " + path.string());

        // Buffered I/O for throughput -- default iostream buffering issues
        // many small writes for a file this size; a 64KB buffer cuts that
        // substantially. Matches the read-side buffering already used in
        // run_find_r_peaks.hpp's read_input_binfile.
        char write_buf[1 << 16];
        out.rdbuf()->pubsetbuf(write_buf, sizeof(write_buf));

        uint64_t n = segs.size();
        out.write(reinterpret_cast<char*>(&n), 8);

        if (n > 0) {
            out.write(reinterpret_cast<const char*>(&segs[0].ppgSampleRate), 8);
            out.write(reinterpret_cast<const char*>(&segs[0].ecgSampleRate), 8);
            out.write(reinterpret_cast<const char*>(&segs[0].scoring_epoch_size_sec), 8);
        }
        const uint32_t nch = NUM_CHANNELS;
        out.write(reinterpret_cast<const char*>(&nch), 4);
        out.write(reinterpret_cast<const char*>(extras.nativeRates.data()),
            NUM_CHANNELS * sizeof(float));

        // Per-channel "Inverted Lead?" checkbox state (noise_marking_gui.ui:
        // ecg_1_reverse / ecg_2_reverse / ecg_3_reverse), one value for the
        // whole file. Must land here -- immediately after nativeRates,
        // before any per-segment data -- to match read_input_binfile() in
        // run_find_r_peaks.hpp.
        const uint8_t inv1 = ecg1_inverted ? 1 : 0;
        const uint8_t inv2 = ecg2_inverted ? 1 : 0;
        const uint8_t inv3 = ecg3_inverted ? 1 : 0;
        out.write(reinterpret_cast<const char*>(&inv1), 1);
        out.write(reinterpret_cast<const char*>(&inv2), 1);
        out.write(reinterpret_cast<const char*>(&inv3), 1);

        auto writePairs = [&](const std::vector<std::pair<uint64_t, uint64_t>>& v) {
            uint64_t sz = v.size();
            out.write(reinterpret_cast<char*>(&sz), 8);
            for (const auto& p : v) {
                out.write(reinterpret_cast<const char*>(&p.first), 8);
                out.write(reinterpret_cast<const char*>(&p.second), 8);
            }
            };
        auto writeVec = [&](const std::vector<double>& v) {
            uint64_t sz = v.size();
            out.write(reinterpret_cast<char*>(&sz), 8);
            if (!v.empty())
                out.write(reinterpret_cast<const char*>(v.data()), sz * 8);
            };

        const size_t ecgN = extras.upsampled[1].size();   // slot 1 = ECG1
        const double sr = (extras.signal_rate > 0)
            ? static_cast<double>(extras.signal_rate) : 1000.0;

        // Per-channel forward cursors into rwSrc, preserved across segments.
        // Raw pairs are time-sorted and segments are time-ordered, so we never
        // revisit pairs we've already passed. Without this, MESA files (~500
        // segments, dense raw blocks) re-scanned each raw vector from the
        // start for every segment -- quadratic in the number of segments.
        std::vector<size_t> rawCursor(NUM_CHANNELS, 0);

        for (const auto& s : segs) {
            writePairs(s.ppg_bin_indexs);
            writePairs(s.ecg_bin_indexs);
            writeVec(s.ppg);
            writeVec(s.ecg1);
            writeVec(s.ecg2);
            writeVec(s.ecg3);
            writeVec(s.sleep_stages);

            for (int ch = 0; ch < NUM_CHANNELS; ++ch) {
                const auto& upSrc = extras.upsampled[ch];
                const auto& rwSrc = extras.rawFlat[ch];

                // ---- Upsampled slice via proportional indexing ----
                // Sentinel passthrough: file_to_bin writes a single -1.0 for
                // missing channels. Pass it through so the schema is uniform.
                std::vector<double> upOut;
                if (upSrc.size() == 1 && upSrc[0] == -1.0) {
                    upOut = upSrc;
                }
                else if (!upSrc.empty() && ecgN > 0) {
                    const double scale = static_cast<double>(upSrc.size()) /
                        static_cast<double>(ecgN);
                    size_t totalSamples = 0;
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        uint64_t a = static_cast<uint64_t>((p.first - 1) * scale) + 1;
                        uint64_t b = static_cast<uint64_t>((p.second - 1) * scale) + 1;
                        if (a > upSrc.size()) continue;
                        if (b > upSrc.size()) b = upSrc.size();
                        totalSamples += static_cast<size_t>(b - a + 1);
                    }
                    upOut.reserve(totalSamples);
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        uint64_t a = static_cast<uint64_t>((p.first - 1) * scale) + 1;
                        uint64_t b = static_cast<uint64_t>((p.second - 1) * scale) + 1;
                        if (a > upSrc.size()) continue;
                        if (b > upSrc.size()) b = upSrc.size();
                        upOut.insert(upOut.end(),
                            upSrc.begin() + (a - 1),
                            upSrc.begin() + b);
                    }
                }
                writeVec(upOut);

                // ---- Raw (t, v) slice via time-window filter ----
                // Sentinel passthrough: missing raw channel is a single (-1, -1).
                std::vector<double> rwOut;
                const bool rawIsSentinel = (rwSrc.size() == 2 &&
                    rwSrc[0] == -1.0 && rwSrc[1] == -1.0);
                if (rawIsSentinel) {
                    rwOut = { -1.0, -1.0 };
                }
                else if (!rwSrc.empty()) {
                    size_t k = rawCursor[ch];
                    const size_t N = rwSrc.size();
                    for (const auto& p : s.ecg_bin_indexs) {
                        if (p.second < p.first) continue;
                        const double t0 = static_cast<double>(p.first - 1) / sr;
                        const double t1 = static_cast<double>(p.second - 1) / sr;
                        while (k + 1 < N && rwSrc[k] < t0) k += 2;
                        while (k + 1 < N && rwSrc[k] <= t1) {
                            rwOut.push_back(rwSrc[k]);
                            rwOut.push_back(rwSrc[k + 1]);
                            k += 2;
                        }
                    }
                    rawCursor[ch] = k;
                }
                const uint64_t nPairs = rwOut.size() / 2;
                out.write(reinterpret_cast<const char*>(&nPairs), 8);
                if (!rwOut.empty())
                    out.write(reinterpret_cast<const char*>(rwOut.data()),
                        rwOut.size() * 8);
            }
        }
    }

}   // anonymous namespace

// ============================================================================
// Public entry point
// ============================================================================

bool annealOneFile(const std::filesystem::path& binPath,
    const std::filesystem::path& noisePath,
    const std::filesystem::path& outPath,
    double binLengthMin,
    bool ecg1_inverted, bool ecg2_inverted, bool ecg3_inverted)
{
    try {
        RawData raw;
        Extras  extras;
        read_data_bin(binPath, raw, extras);

        NoiseMarkings noise;
        if (std::filesystem::exists(noisePath))
            noise = read_noise_bin(noisePath);
        else
            std::cerr << "  no noise file at " << noisePath
            << " -- annealing with no exclusions\n";

        auto results = AnnealSegments(raw, noise, binLengthMin);
        write_output_bin(outPath, results, extras, ecg1_inverted, ecg2_inverted, ecg3_inverted);

        std::cerr << "  -> " << results.size() << " bins -> "
            << outPath.filename() << "\n";
        return true;
    }
    catch (const std::exception& e) {
        std::cerr << "  ERROR: " << e.what() << std::endl;
        return false;
    }
    catch (...) {
        std::cerr << "  ERROR: unknown exception type\n" << std::flush;
        return false;
    }
}