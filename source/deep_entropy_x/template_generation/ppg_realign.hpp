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

    // ---- HOW MANY ROWS A COLUMN NEEDS BEFORE IT GETS A VALUE ------------
    //
    // The beat matrix is padded to the bin's LONGEST slice, so its far tail is
    // populated by only the few longest beats -- the same far-tail problem
    // build_bins.hpp documents for the ECG, where a dropped R detection is the
    // only row with samples out there and the column median becomes its later
    // complexes. A median over three rows out of a thousand is not an average
    // and the spread of those three is not a corridor, but both come out
    // finite and get drawn, which is how a re-levelled template ended up with
    // a band the height of the panel.
    //
    // HALF THE USED ROWS. The bank slots the viewer draws are already trimmed
    // to roughly where their own membership runs out (two slots of one bin came
    // back 747 and 819 columns wide on a 1466-column matrix), so this
    // reproduces that trim from the same evidence rather than guessing a width.
    inline int minColumnRows(int n_used) {
        return std::max(3, n_used / 2);
    }


    // ---- THE SPREAD THE BANK SLOT FIELD ACTUALLY HOLDS ------------------
    //
    // RAW-AMPLITUDE q3 - q1 PER COLUMN, in the same units as tmpl. Copied
    // from template_assign.hpp, which is what writes tbank::BankTemplate::
    // tmpl_iqr at build time -- and the units are the whole point:
    //
    //   pulseTraceForSlot draws the band as
    //     scale_pulse_spread_by_ref(slot.tmpl_iqr, sample_y(tmpl, foot), ref)
    //
    // so the DISPLAY divides the stored spread by the foot amplitude. A
    // perfusion-ratio spread (local_ratio_iqr) has already divided by that
    // foot, so storing one here divides by it TWICE. On a pulse whose stored
    // baseline is near zero -- and it is: -0.065 on a real record -- that is a
    // factor of ~1500, and the band came out at 6.4 where the build's was
    // 0.0049 and filled the entire panel.
    //
    // local_ratio_iqr is the right statistic for the BIN template
    // (build_pulse_template_pair_windowed uses it, and nothing divides that
    // one again). It is the wrong statistic for a bank slot. The two fields
    // share a name and not a unit.
    inline std::vector<double> rawIqrColumns(
        const std::vector<std::vector<double>>& rows, int n_used)
    {
        std::size_t w = 0;
        for (const auto& r : rows) w = std::max(w, r.size());
        std::vector<double> out(w, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> col;
        col.reserve(rows.size());
        for (std::size_t c = 0; c < w; ++c) {
            col.clear();
            for (const auto& r : rows)
                if (c < r.size() && !std::isnan(r[c])) col.push_back(r[c]);
            if ((int)col.size() < minColumnRows(n_used)) continue;
            const std::size_t n = col.size();
            const std::size_t iq1 = n / 4;
            const std::size_t iq3 = (3 * n) / 4;
            std::nth_element(col.begin(), col.begin() + iq1, col.end());
            const double q1 = col[iq1];
            std::nth_element(col.begin() + iq1, col.begin() + iq3, col.end());
            const double q3 = col[iq3];
            out[c] = q3 - q1;
        }
        return out;
    }

    struct Result {
        bool ok = false;
        std::string why;                 // set when !ok, for the status bar
        std::vector<double> tmpl;        // re-medianed waveform
        std::vector<double> iqr;         // local-ratio IQR about foot_col
        int anchor_col = -1;             // where every beat's ANCHOR now sits
        // WHERE THE FOOT ENDED UP, which is NOT anchor_col once pct > 0. The
        // two were one field while the only alignment was the foot, and the
        // difference matters to the two things that read it: the perfusion
        // ratio underneath the IQR divides by the foot's amplitude (see
        // local_ratio_iqr), and the operator's foot bar is drawn at the foot.
        // Anchoring 10% up the upstroke and then dividing by the value there
        // would scale the band by a level that is not the baseline.
        int foot_col = -1;
        double pct = 0.0;                // % up the upstroke this stacked on
        int n_members = 0;               // cohort offered
        int n_used = 0;                  // cohort with an anchor in the window
        int max_shift = 0;               // largest |shift| applied, samples
        double median_shift = 0.0;       // median signed shift, samples
    };

    // ------------------------------------------------------------------
    // WHERE A PULSE IS `pct` PERCENT UP ITS OWN UPSTROKE
    // ------------------------------------------------------------------
    //
    // IN AMPLITUDE, NOT IN TIME. The column returned is the first one on the
    // rising edge at which the pulse reaches
    //
    //     foot_y + (pct/100) * (peak_y - foot_y)
    //
    // so 0 IS the foot -- returned exactly, with no search -- and 100 is the
    // systolic peak. A time fraction would be a different physiological
    // instant on every beat, because upstroke DURATION varies with rate and
    // contractility while the fraction of the rise does not; stacking on a
    // time fraction would therefore smear the very feature it was meant to
    // sharpen.
    //
    // WHY ANY OF THIS IS WORTH OFFERING OVER THE FOOT. The foot is the worst
    // -conditioned landmark on a pulse: it is a turning point, the slope
    // through it is zero by definition, and a millivolt of noise moves it by
    // tens of milliseconds. Partway up the upstroke the slope is at its
    // steepest, so the same noise moves the crossing by almost nothing.
    //
    // THE BUILD ALREADY AGREES. extract_ppg_beats_and_align anchors every beat
    // to the median up50 column -- its 50% foot-to-peak crossing -- and
    // discards any beat whose up50 cannot be found. So this is not a new idea
    // being introduced at the viewer: it is the pipeline's own alignment,
    // exposed with the fraction as a control instead of frozen at one half.
    //
    // ---- IT IS THE BUILD'S OWN MEASUREMENT, NOT A SECOND ONE ------------
    //
    // FeatureMarks::first_crossing(v, foot, peak, frac) is the primitive
    // extract_ppg_beats_and_align calls with frac = 0.50 to find the up50
    // column it anchors every beat to -- so the build's pulse alignment IS
    // this function at 50, and "Auto" in the viewer's radio group is that same
    // alignment as the pipeline left it.
    //
    // Which is exactly why this delegates rather than computing a crossing of
    // its own. A private copy would be a SECOND definition of "% up the
    // upstroke", and the one thing worse than the viewer disagreeing with the
    // build about where 50% is would be its doing so by a sub-sample, where
    // nothing on screen could show it. first_crossing interpolates without
    // clamping the straddle fraction, deliberately (see its comment: clamping
    // moved a handful of up50 columns and shifted the shared width by one) --
    // reproducing that by hand is not the kind of detail that survives being
    // reproduced.
    //
    // THE PEAK IS THE BRACKET, and it comes from the upstroke rather than
    // argmax for the reason the slicer gives: on a pulse whose reflected wave
    // exceeds systole, argmax is the SECOND peak, the crossing gets searched
    // on the notch-to-P2 rise instead of the systolic one, and there is often
    // no clean single crossing there at all.
    //
    // Sub-sample. The caller rounds if it needs an integer; baking a
    // half-sample bias in here would put it into every beat's shift in the
    // same direction. -1 when the pulse has no resolvable upstroke.
    inline double upstrokePctCol(const std::vector<double>& v, int foot,
        double pct)
    {
        const int n = static_cast<int>(v.size());
        if (foot < 0 || foot >= n) return -1.0;
        if (!std::isfinite(v[static_cast<size_t>(foot)])) return -1.0;
        // 0 IS THE FOOT, AS ITS OWN CASE AND NOT AS A LIMIT. first_crossing
        // needs v[i-1] to sit strictly BELOW the target to call it a crossing,
        // and at frac = 0 the target is the foot's own amplitude -- so it
        // reports no crossing and returns -1, which here would read as "this
        // pulse has no upstroke" rather than "you asked for the foot".
        if (!(pct > 0.0)) return static_cast<double>(foot);
        if (pct > 100.0) pct = 100.0;

        const int peak = FeatureMarks::detect_ppg_upstroke_peak(v, foot, n);
        if (peak <= foot || peak >= n) return -1.0;

        return FeatureMarks::first_crossing(v, foot, peak, pct / 100.0);
    }

    // `rows` is the member list IN LOCAL ROW SPACE -- pass members_clean when
    // it is non-empty, since that is the set the existing waveform was averaged
    // over and the point here is to re-average THE SAME SET.
    //
    // search_halfwin is in samples. The caller sets it from the pulse rate; see
    // kDefaultSearchSeconds.
    //
    // ---- THE TWO COLUMNS, AND WHY THEY ARE BOTH PARAMETERS --------------
    //
    // operatorFootCol is WHERE THE FOOT IS: the column the operator dragged
    // the foot bar to, and the hint every beat's own trough is searched
    // around.
    //
    // anchorCol is WHERE THE STACK LINES UP. At pct = 0 the two are the same
    // column and this is the original foot re-stack exactly. At pct > 0 the
    // caller measures the percent point ON THE SLOT'S CURRENT AVERAGE (via
    // upstrokePctCol, from that same foot) and passes it here, so the foot
    // still lands on the operator's column while the beats are lined up on
    // their steepest, best-conditioned instant. Deriving anchorCol in here
    // instead would mean measuring it on a beat rather than on the average,
    // and one beat is exactly the noisy measurement the percent alignment
    // exists to avoid.
    //
    // ONE NUMBER STILL CANNOT BE AN ALIGNMENT (see the header note): anchorCol
    // is only the column the stack is anchored AT. Every beat contributes its
    // own crossing, so the shape sharpens; a rigid offset would only translate
    // the template.
    inline Result realignAt(const BinBeats& beats,
        const std::vector<uint32_t>& rows,
        double operatorFootCol,
        int search_halfwin,
        double pct,
        double anchorCol)
    {
        Result out;
        out.n_members = static_cast<int>(rows.size());
        out.pct = (pct > 0.0) ? std::min(pct, 100.0) : 0.0;
        if (beats.empty()) { out.why = "no beat matrix for this bin"; return out; }
        const int w = beats.width;
        const int foot_target = static_cast<int>(std::lround(operatorFootCol));
        if (foot_target < 0 || foot_target >= w) { out.why = "foot bar outside the beat window"; return out; }
        if (search_halfwin < 1) search_halfwin = 1;

        const int target = (out.pct > 0.0 && anchorCol >= 0.0)
            ? static_cast<int>(std::lround(anchorCol)) : foot_target;
        if (target < 0 || target >= w) {
            out.why = "alignment point outside the beat window";
            return out;
        }

        // ANCHORED ON THE OPERATOR'S OWN COLUMN, not on the median of the
        // anchors found below. The bar stays where they put it, the re-stacked
        // template lines up under it, and the marks already stored against this
        // slot do not need translating -- whereas anchoring to the median would
        // slide the whole pulse out from under every bar on the panel.
        out.anchor_col = target;
        out.foot_col = foot_target;

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
                foot_target - search_halfwin, foot_target + search_halfwin);
            if (foot < 0) continue;

            // AND THEN ITS OWN CROSSING, FROM ITS OWN FOOT. A beat whose
            // upstroke cannot be resolved is skipped rather than aligned on its
            // foot as a fallback: a stack in which most rows are held at the
            // 10% level and a few at their feet is aligned on neither, and the
            // n_used the caller reports would no longer mean what it says.
            int anchor = foot;
            if (out.pct > 0.0) {
                const double c = upstrokePctCol(src, foot, out.pct);
                if (!(c >= 0.0)) continue;
                anchor = static_cast<int>(std::lround(c));
            }

            const int shift = target - anchor;
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
            out.why = (out.pct > 0.0)
                ? "fewer than 3 beats have a resolvable upstroke near that column"
                : "fewer than 3 beats have a foot near that column";
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
            // LEFT NaN, NOT ZERO, when too few rows reach here: NaN is how
            // every reader downstream already spells "no average at this
            // column", and it is what the tail trim keys on.
            if ((int)col.size() < minColumnRows(out.n_used)) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            out.tmpl[static_cast<size_t>(c)] = (nc % 2)
                ? hi
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi);
        }

        // The spread comes from the same primitive the build uses: per-beat
        // local perfusion ratios, then the cross-beat IQR of those.
        // Recomputed rather than carried over, because the stored tmpl_iqr
        // describes the OLD stacking and would draw a band that no longer
        // belongs to the trace inside it.
        //
        // RAW-AMPLITUDE UNITS, because that is what this field is drawn as.
        // See rawIqrColumns -- a perfusion-ratio spread here is divided by the
        // foot twice and fills the panel.
        out.iqr = rawIqrColumns(shifted, out.n_used);

        out.ok = true;
        return out;
    }

    // The foot alignment, which is realignAt with no percent and no separate
    // anchor. Kept as its own name because it is what a foot-bar drag means
    // and because every existing caller passes exactly these four arguments.
    inline Result realign(const BinBeats& beats,
        const std::vector<uint32_t>& rows,
        double operatorFootCol,
        int search_halfwin)
    {
        return realignAt(beats, rows, operatorFootCol, search_halfwin,
            /*pct=*/0.0, /*anchorCol=*/operatorFootCol);
    }

    // ======================================================================
    // THE VERTICAL SIBLING: RE-LEVEL A COHORT ON THE OPERATOR'S FOOT
    // ======================================================================
    //
    // A DIFFERENT AXIS FROM realignAt, AND THAT IS THE WHOLE DISTINCTION. The
    // "Align PPG Horizontal" group re-TIMES the stack -- it decides which
    // instant every beat is shifted onto. The foot bar does not: the foot is
    // the pulse's vertical reference, in the build (pass 2 of
    // extract_ppg_beats_and_align DC-matches each beat's baseline to the
    // reference beat's) and in the viewer (normalize_ppg_or_similar computes
    // 100*(sample - footY)/footY, so the foot column is subtracted AND divided
    // by, and the band is scaled about the same column).
    //
    // So "re-do the foot alignment" means: level every member row so they
    // agree at the column the operator put the bar on, then re-median. No
    // sample moves sideways.
    //
    // ONE SHARED COLUMN, NOT A SEARCH. The bar says "the foot is HERE", so
    // that column is read on every row directly -- no per-row trough hunt,
    // which would re-derive an answer the operator has just overridden. It is
    // a slightly stronger claim than the detector makes (the detector returns
    // a foot per beat), and it is a fair one here because the rows arrive
    // already up50 time-aligned, so their feet very nearly coincide in column
    // anyway.
    //
    // EACH ROW GETS ITS OWN OFFSET, which is what makes this more than a
    // cosmetic shift of the average: the rows' relative vertical positions
    // change, so the column median changes at every other column too. The
    // waveform's SHAPE moves, not just its height.
    //
    // ---- AND IT LEVELS TO A VALUE, NOT TO ZERO --------------------------
    //
    // Subtracting each row's own level outright would leave tmpl[foot] == 0 --
    // and that is precisely the divisor normalize_ppg_or_similar uses. A
    // near-zero divisor is the failure maybeNotchTrace's comment describes: a
    // tiny DC change becomes a huge relative one and the normalized trace
    // collapses to the bottom of the axis, looking like the channel vanished.
    //
    // So the rows are brought to a COMMON BASELINE VALUE -- the median of
    // their own levels, a real amplitude from this cohort -- exactly as the
    // build matches each beat to its reference beat rather than to zero.
    struct RelevelResult {
        bool ok = false;
        std::string why;
        std::vector<double> tmpl;     // re-medianed waveform
        std::vector<double> iqr;      // local-ratio IQR about foot_col
        int    foot_col = -1;         // the operator's column
        double baseline = 0.0;        // the common level every row was brought to
        int    n_members = 0;         // cohort offered
        int    n_used = 0;            // cohort with a readable level there
        double median_level = 0.0;    // median signed offset applied
        double max_level = 0.0;       // largest |offset| applied
    };

    // mean_halfwin is in samples: the level for a row is the MEAN of the
    // finite samples in [foot - hw, foot + hw], not the single sample at the
    // foot. One noisy sample would otherwise set that row's entire offset, and
    // the build takes a windowed mean for the same reason ("avoids
    // order-statistic bias").
    inline RelevelResult relevelAt(const BinBeats& beats,
        const std::vector<uint32_t>& rows,
        double footCol,
        int mean_halfwin)
    {
        RelevelResult out;
        out.n_members = static_cast<int>(rows.size());
        if (beats.empty()) { out.why = "no beat matrix for this bin"; return out; }
        const int w = beats.width;
        const int target = static_cast<int>(std::lround(footCol));
        if (target < 0 || target >= w) {
            out.why = "foot bar outside the beat window";
            return out;
        }
        if (mean_halfwin < 0) mean_halfwin = 0;
        out.foot_col = target;

        // ---- pass one: each row's own level at that column ---------------
        std::vector<const std::vector<double>*> src;
        std::vector<double> level;
        src.reserve(rows.size());
        level.reserve(rows.size());
        for (const uint32_t r : rows) {
            if (r >= beats.rows.size()) continue;       // stale member list
            const std::vector<double>& row = beats.rows[r];
            const int lo = std::max(0, target - mean_halfwin);
            const int hi = std::min<int>(static_cast<int>(row.size()) - 1,
                target + mean_halfwin);
            double sum = 0.0; int n = 0;
            for (int k = lo; k <= hi; ++k)
                if (!std::isnan(row[static_cast<size_t>(k)])) {
                    sum += row[static_cast<size_t>(k)]; ++n;
                }
            // ALL NaN THERE means this beat's samples do not reach the
            // operator's column -- a clipped row, skipped rather than levelled
            // on a guess.
            if (n == 0) continue;
            src.push_back(&row);
            level.push_back(sum / static_cast<double>(n));
        }

        out.n_used = static_cast<int>(src.size());
        // TWO IS NOT A TEMPLATE -- the same refusal as realignAt, and for the
        // same reason: leaving the built waveform in place beats replacing it
        // with a median over two rows and reporting success.
        if (out.n_used < 3) {
            out.why = "fewer than 3 beats have a readable baseline at that column";
            return out;
        }

        // The common baseline: the median of the levels themselves. Robust to
        // a row whose baseline is an outlier, and guaranteed to be a value
        // this cohort actually exhibits.
        {
            std::vector<double> s = level;
            std::sort(s.begin(), s.end());
            out.baseline = s[s.size() / 2];
        }

        // ---- pass two: bring every row to it ----------------------------
        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> levelled;
        levelled.reserve(src.size());
        std::vector<double> offs;
        offs.reserve(src.size());
        for (std::size_t i = 0; i < src.size(); ++i) {
            const double off = out.baseline - level[i];
            std::vector<double> dst(static_cast<size_t>(w), kNaN);
            const int n = std::min<int>(w, static_cast<int>(src[i]->size()));
            for (int k = 0; k < n; ++k) {
                const double v = (*src[i])[static_cast<size_t>(k)];
                if (!std::isnan(v)) dst[static_cast<size_t>(k)] = v + off;
            }
            levelled.push_back(std::move(dst));
            offs.push_back(off);
            if (std::abs(off) > std::abs(out.max_level)) out.max_level = off;
        }
        {
            std::vector<double> s = offs;
            std::sort(s.begin(), s.end());
            out.median_level = s[s.size() / 2];
        }

        // Column-wise NaN-skipping median, the same statistic in the same
        // nth_element form realignAt and the build both use. A different
        // average here would make the re-levelled template incomparable with
        // the one it replaces.
        out.tmpl.assign(static_cast<size_t>(w), kNaN);
        std::vector<double> col;
        col.reserve(levelled.size());
        for (int c = 0; c < w; ++c) {
            col.clear();
            for (const auto& b : levelled)
                if (!std::isnan(b[static_cast<size_t>(c)]))
                    col.push_back(b[static_cast<size_t>(c)]);
            if ((int)col.size() < minColumnRows(out.n_used)) continue;
            const size_t nc = col.size();
            const size_t mid = nc / 2;
            std::nth_element(col.begin(), col.begin() + mid, col.end());
            const double hi = col[mid];
            out.tmpl[static_cast<size_t>(c)] = (nc % 2)
                ? hi
                : 0.5 * (*std::max_element(col.begin(), col.begin() + mid) + hi);
        }

        // ---- ON THE LEVELLED ROWS, IN RAW UNITS ---------------------------
        //
        // THE LEVELLED ONES, because the corridor is the EVIDENCE the levelling
        // took: if the beats really do share a baseline at the operator's
        // column then their spread must be small there and grow away from it,
        // and that is exactly what the operator needs to see. The spread at the
        // bar comes out at 0 because every row was brought to the same value
        // there -- that is the signal, not a defect.
        //
        // RAW AMPLITUDE, not a perfusion ratio: the display divides this field
        // by the foot amplitude itself. See rawIqrColumns.
        out.iqr = rawIqrColumns(levelled, out.n_used);

        out.ok = true;
        return out;
    }

    // 10 ms either side of the foot, the baseline-mean window. Short enough
    // that it cannot climb the upstroke at any plausible rate, long enough to
    // average a few samples of noise out of the level.
    inline constexpr double kDefaultLevelSeconds = 0.010;

    inline int levelHalfWinFor(double pulseRateHz) {
        if (!(pulseRateHz > 0.0)) return 2;
        return std::max(1, static_cast<int>(
            std::lround(kDefaultLevelSeconds * pulseRateHz)));
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

}   // namespace ppg_realign