#pragma once
/**
 * @file   ppg_realign.hpp
 * @brief  Re-stack one pulse template about an operator-corrected foot.
 *
 *         WHAT THIS IS FOR. The pulse template a column draws is a column-wise
 *         median over the beats the build assigned to that morphology, stacked
 *         on whatever fiducial alignment.hpp was told to use. When the foot
 *         detector is wrong it is usually wrong the SAME WAY on every beat in
 *         the bin -- it lands on the dicrotic notch, or on the reflected wave's
 *         upstroke -- and the averaged template then smears about the wrong
 *         landmark. One operator correction, dragged on the average, is enough
 *         to say where the foot really is; this re-detects each beat's own foot
 *         near that column and re-medians the stack about it.
 *
 *         WHAT IT DELIBERATELY DOES NOT DO. It does not re-run the pulse QC
 *         filter, and it does not touch membership: the survivor set is exactly
 *         the one the build chose. Re-filtering here would change which beats
 *         are in the group, which changes the morphology split, which is the
 *         pipeline's job and not a mouse-release's. So this is an honest
 *         re-average of a fixed cohort and nothing more -- and the template it
 *         produces is anchored one way while its cohort was selected under
 *         another, which the caller should say out loud rather than hide.
 *
 *         ONE NUMBER CANNOT BE AN ALIGNMENT. The dragged bar is a single column
 *         on the AVERAGE; a re-stack needs a foot per beat. Applying the
 *         operator's column as a rigid offset to every beat would translate the
 *         template and change nothing about its shape, which is not what the
 *         gesture means. So the column is used as a SEARCH HINT: each beat's
 *         own trough within +-search_halfwin of it. That is why the correction
 *         can be large (notch to foot, 200 ms) while the window stays small --
 *         the window only has to cover beat-to-beat scatter of the true foot,
 *         not the size of the correction.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "template_morphology_grouping/morphology_csv.hpp"
#include "fiducial_marker_finding/feature_marks.hpp"
#include "template_generation/normalize_template_amplitude.hpp"

namespace ppg_realign {

    // One bin's pulse beats, in the channel's LOCAL ROW SPACE -- the space
    // tbank::BankTemplate::members is in (see the join note on loadBin).
    struct BinBeats {
        int width = 0;
        std::vector<std::vector<double>> rows;
        bool empty() const { return rows.empty() || width <= 0; }
    };

    // ------------------------------------------------------------------
    // Reading ONE bin's pulse block out of <stem>_beats.bin
    // ------------------------------------------------------------------
    //
    // NOT morphology_csv::readBeatsBin, and this is the one place in the
    // codebase that reads the format without it. readBeatsBin materialises
    // EVERY block: nCols x width doubles per channel, four channels, and on a
    // full-night record that is hundreds of megabytes to answer a question
    // about one bin. This walks the same layout instead, seeking over the
    // blocks and columns it does not want, and allocates only the bin asked
    // for.
    //
    // The LAYOUT still has exactly one definition: morphology_csv::BeatRecord
    // and morphology_csv::bin_version are used directly, so a change to either
    // reaches this reader as a compile-time or version-check failure rather
    // than as silent misparsing. The traversal is duplicated; the format is
    // not.
    //
    // ---- THE JOIN, WHICH IS THE PART THAT IS EASY TO GET WRONG ----------
    //
    // _beats.bin has ONE COLUMN PER SLICE (an R-pair), present whether or not
    // this channel produced a beat for it; became_beat says which. But
    // BankTemplate::members holds LOCAL ROWS -- jbank::projectToChannel
    // translates each group member from slice to row via
    // ChannelBeats::local_of_slice before the archive ever sees it. The two
    // spaces are NOT the same, despite what morphology_csv::BinBlock::members'
    // own comment says.
    //
    // They are reconcilable because row order is slice order: the pulse slicer
    // appends beats in increasing slice index and every prune downstream
    // preserves that order, so local row k is the k-th became_beat column of
    // that bin. Counting them in file order reproduces local_of_slice exactly,
    // which is what the loop below does -- rows.size() is the next row index.
    inline BinBeats loadBin(const std::string& beatsPath, uint32_t bin,
        const char* channel = "PPG")
    {
        BinBeats out;
        std::ifstream f(beatsPath, std::ios::binary);
        if (!f) return out;

        auto rd = [&f](void* dst, std::streamsize n) -> bool {
            return static_cast<bool>(f.read(static_cast<char*>(dst), n));
            };

        uint32_t ver = 0, nBlocks = 0;
        if (!rd(&ver, 4) || !rd(&nBlocks, 4)) return out;
        if (ver > morphology_csv::bin_version) return out;   // newer than this build

        const std::streamsize recSize =
            static_cast<std::streamsize>(sizeof(morphology_csv::BeatRecord));

        for (uint32_t bi = 0; bi < nBlocks; ++bi) {
            uint32_t len = 0;
            if (!rd(&len, 4) || len > 64) return out;
            std::string name(len, '\0');
            if (len && !rd(name.data(), len)) return out;
            uint32_t width = 0;
            uint64_t nCols = 0;
            if (!rd(&width, 4) || !rd(&nCols, 8)) return out;

            const std::streamsize sampleBytes =
                static_cast<std::streamsize>(width) * 8;

            if (name != channel) {
                // Whole block skipped in one seek. A block is a fixed stride:
                // every column is one record plus `width` doubles, with no
                // per-column length prefix, which is what makes this possible.
                f.seekg((recSize + sampleBytes) * static_cast<std::streamoff>(nCols),
                    std::ios::cur);
                if (!f) return out;
                continue;
            }

            out.width = static_cast<int>(width);
            for (uint64_t k = 0; k < nCols; ++k) {
                morphology_csv::BeatRecord rec;
                if (!rd(&rec, recSize)) { out.rows.clear(); return out; }
                const bool want = (rec.bin == bin) && (rec.became_beat != 0);
                if (!want || width == 0) {
                    if (sampleBytes) f.seekg(sampleBytes, std::ios::cur);
                    if (!f) { out.rows.clear(); return out; }
                    continue;
                }
                std::vector<double> row(width);
                if (!rd(row.data(), sampleBytes)) { out.rows.clear(); return out; }
                out.rows.push_back(std::move(row));
            }
            return out;   // the block we came for; nothing after it matters
        }
        return out;
    }

    // ------------------------------------------------------------------
    // The re-stack
    // ------------------------------------------------------------------

    struct Result {
        bool ok = false;
        std::string why;                 // set when !ok, for the status bar
        std::vector<double> tmpl;        // re-medianed waveform
        std::vector<double> iqr;         // local-ratio IQR about anchor_col
        int anchor_col = -1;             // where every beat's foot now sits
        int n_members = 0;               // cohort offered
        int n_used = 0;                  // cohort with a foot in the window
        int max_shift = 0;               // largest |shift| applied, samples
        double median_shift = 0.0;       // median signed shift, samples
    };

    // `rows` is the member list IN LOCAL ROW SPACE -- pass members_clean when
    // it is non-empty, since that is the set the existing waveform was averaged
    // over and the point here is to re-average THE SAME SET.
    //
    // search_halfwin is in samples. The caller sets it from the pulse rate; see
    // kDefaultSearchSeconds.
    inline Result realign(const BinBeats& beats,
        const std::vector<uint32_t>& rows,
        double operatorFootCol,
        int search_halfwin)
    {
        Result out;
        out.n_members = static_cast<int>(rows.size());
        if (beats.empty()) { out.why = "no beat matrix for this bin"; return out; }
        const int w = beats.width;
        const int target = static_cast<int>(std::lround(operatorFootCol));
        if (target < 0 || target >= w) { out.why = "foot bar outside the beat window"; return out; }
        if (search_halfwin < 1) search_halfwin = 1;

        // ANCHORED ON THE OPERATOR'S OWN COLUMN, not on the median of the feet
        // found below. The bar stays where they put it, the re-stacked template
        // lines up under it, and the marks already stored against this slot do
        // not need translating -- whereas anchoring to the median would slide
        // the whole pulse out from under every bar on the panel.
        out.anchor_col = target;

        std::vector<std::vector<double>> shifted;
        shifted.reserve(rows.size());
        std::vector<int> shifts;
        shifts.reserve(rows.size());
        const double kNaN = std::numeric_limits<double>::quiet_NaN();

        for (const uint32_t r : rows) {
            if (r >= beats.rows.size()) continue;       // stale member list
            const std::vector<double>& src = beats.rows[r];
            // Each beat's OWN trough, bounded to the hint window. trough_in
            // clamps its own bounds and returns -1 when the window is all NaN,
            // which is a beat whose samples do not reach here -- skipped, not
            // forced to a column.
            const int foot = FeatureMarks::trough_in(src,
                target - search_halfwin, target + search_halfwin);
            if (foot < 0) continue;

            const int shift = target - foot;
            std::vector<double> dst(static_cast<size_t>(w), kNaN);
            for (int k = 0; k < (int)src.size(); ++k) {
                const int d = k + shift;
                if (d >= 0 && d < w) dst[static_cast<size_t>(d)] = src[k];
            }
            shifted.push_back(std::move(dst));
            shifts.push_back(shift);
            if (std::abs(shift) > out.max_shift) out.max_shift = std::abs(shift);
        }

        out.n_used = static_cast<int>(shifted.size());
        // TWO IS NOT A TEMPLATE. Refusing here leaves the built waveform in
        // place, which is the right outcome: a median over one or two beats
        // would replace a real template with noise and look like a successful
        // re-alignment.
        if (out.n_used < 3) {
            out.why = "fewer than 3 beats have a foot near that column";
            return out;
        }

        {
            std::vector<int> s = shifts;
            std::sort(s.begin(), s.end());
            out.median_shift = s[s.size() / 2];
        }

        // Column-wise NaN-skipping median, the same statistic (and the same
        // nth_element form) the build uses in build_pulse_template_pair_
        // windowed. A different average here would mean the re-stacked
        // template and the built one were not comparable.
        out.tmpl.assign(static_cast<size_t>(w), kNaN);
        std::vector<double> col;
        col.reserve(shifted.size());
        for (int c = 0; c < w; ++c) {
            col.clear();
            for (const auto& b : shifted)
                if (!std::isnan(b[static_cast<size_t>(c)]))
                    col.push_back(b[static_cast<size_t>(c)]);
            if (col.empty()) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            out.tmpl[static_cast<size_t>(c)] = (nc % 2)
                ? hi
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi);
        }

        // The spread comes from the same primitive the build uses, about the
        // new anchor: per-beat local perfusion ratios, then the cross-beat IQR
        // of those. Recomputed rather than carried over, because the stored
        // tmpl_iqr describes the OLD stacking and would draw a band that no
        // longer belongs to the trace inside it.
        out.iqr = normalize_features::local_ratio_iqr(shifted, out.anchor_col);

        out.ok = true;
        return out;
    }

    // 100 ms either side of the operator's column, the default search window.
    // Wide enough for the beat-to-beat scatter of a real foot; narrow enough
    // that it cannot reach the dicrotic notch, which sits well past 150 ms
    // after systole at any plausible rate.
    inline constexpr double kDefaultSearchSeconds = 0.10;

    inline int searchHalfWinFor(double pulseRateHz) {
        if (!(pulseRateHz > 0.0)) return 5;
        return std::max(3, static_cast<int>(
            std::lround(kDefaultSearchSeconds * pulseRateHz)));
    }

}   // namespace ppg_realign#pragma once
