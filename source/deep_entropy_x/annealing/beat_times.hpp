#pragma once
/**
 * @file   beat_times.hpp
 * @brief  Recording time, in seconds from the start of the record, of every
 *         kept ECG beat -- indexed the way bank members are, [lead][bin][row].
 *
 * ---------------------------------------------------------------------------
 * HOW A MEMBER ROW BECOMES A TIME
 * ---------------------------------------------------------------------------
 *
 * ecg_bank[c].templates[slot].members are LOCAL ROWS of lead c's kept-beat
 * matrix for that bin (per_channel_beats[CHc][bin]). A row is not a heartbeat
 * ordinal -- each lead prunes independently -- so the chain is:
 *
 *   row  --kept_index[c][bin]-->            slice s (the R-pair ordinal)
 *   s    --peakResults[bin].ch1.raw[s]-->   R sample in the bin's SPLICED signal
 *   p    --peakResults[bin].ecg_bin_indexs--> sample in the ORIGINAL record
 *   idx  --/ ECG1 upsample rate-->          seconds
 *
 * Slice s is the pair (ch1.raw[s], ch1.raw[s+1]) and the beat is cut around
 * ch1.raw[s] on EVERY lead (ch1.raw drives all three slicers -- bin_ecg.hpp),
 * so a beat's time is its CH1 R peak whichever lead's row it is.
 *
 * THE SPLICE. The annealer builds a bin by concatenating the record's good
 * stretches: ecg_bin_indexs holds 1-based inclusive (first, last) ranges into
 * the upsampled ECG1 array, copied in order, each clamped to the array's end
 * and skipped if it starts past it (anneal_handler.cpp, extract). Splice below
 * inverts exactly that, so a beat after a removed noise gap gets its real
 * time rather than one short by the gap.
 *
 * Nothing here knows about sleep. record_sleep.hpp turns a time into a stage.
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "peak_finding/peakfinding_structs.hpp"

namespace beat_times {

    // [lead][bin][kept row] -> seconds; NaN where the row has no usable time.
    // Empty inner vectors where a lead has no beats in that bin.
    struct BeatTimes {
        std::array<std::vector<std::vector<double>>, 3> time_s;

        // NaN for anything out of range, so callers need no bounds checks.
        double at(int lead, size_t bin, uint32_t row) const {
            constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
            if (lead < 0 || lead > 2 || bin >= time_s[lead].size()) return kNaN;
            const std::vector<double>& v = time_s[lead][bin];
            return (row < v.size()) ? v[row] : kNaN;
        }
    };

    // One bin's splice, inverted: spliced position -> original 0-based sample.
    struct Splice {
        struct Run { uint64_t concat0, orig0, len; };
        std::vector<Run> runs;

        Splice(const std::vector<std::pair<uint64_t, uint64_t>>& idx, uint64_t sigLen) {
            uint64_t at = 0;
            for (const auto& seg : idx) {
                // extract(): skip a range starting past the array, clamp its
                // end to the array, skip an inverted one. first is 1-based.
                if (seg.first == 0 || seg.first > sigLen) continue;
                const uint64_t last = std::min<uint64_t>(seg.second, sigLen);
                if (last < seg.first) continue;
                const uint64_t len = last - seg.first + 1;
                runs.push_back({ at, seg.first - 1, len });
                at += len;
            }
        }
        // False when p lies past the spliced signal.
        bool original(uint64_t p, uint64_t& out) const {
            for (const Run& r : runs)
                if (p >= r.concat0 && p < r.concat0 + r.len) {
                    out = r.orig0 + (p - r.concat0);
                    return true;
                }
            return false;
        }
    };

    // sliceOfRow[bin][lead] is that lead's kept_index for the bin (row -> slice).
    // ecgRateHz / ecgLength are the UPSAMPLED ECG1 rate and sample count the
    // annealer cut the bins from; pass max() for the length when unknown.
    inline BeatTimes build(const std::vector<output_binfile_data>& peakResults,
        const std::vector<std::array<std::vector<size_t>, 3>>& sliceOfRow,
        double ecgRateHz,
        uint64_t ecgLength = std::numeric_limits<uint64_t>::max())
    {
        BeatTimes bt;
        const size_t nb = std::min(peakResults.size(), sliceOfRow.size());
        for (int c = 0; c < 3; ++c) bt.time_s[c].assign(nb, {});
        if (!(ecgRateHz > 0.0)) return bt;
        constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

        for (size_t b = 0; b < nb; ++b) {
            const output_binfile_data& pr = peakResults[b];
            const Splice sp(pr.ecg_bin_indexs, ecgLength);
            const std::vector<std::size_t>& rp = pr.ch1.raw;
            for (int c = 0; c < 3; ++c) {
                const std::vector<size_t>& rows = sliceOfRow[b][c];
                std::vector<double>& ts = bt.time_s[c][b];
                ts.assign(rows.size(), kNaN);
                for (size_t k = 0; k < rows.size(); ++k) {
                    const size_t s = rows[k];
                    if (s >= rp.size()) continue;
                    uint64_t orig = 0;
                    if (!sp.original(static_cast<uint64_t>(rp[s]), orig)) continue;
                    ts[k] = static_cast<double>(orig) / ecgRateHz;
                }
            }
        }
        return bt;
    }

}   // namespace beat_times