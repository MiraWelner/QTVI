#pragma once
/**
 * @file   ecg_realign.hpp
 * @brief  Re-stack one ECG bank slot about an operator-corrected landmark.
 *
 *         THE GESTURE. Drag the P-onset bar to where the P onset really is,
 *         release. Each member beat's own P onset is re-found, the stack is
 *         re-medianed about the operator's column, and the P_ONSET per-slot
 *         average is replaced in place.
 *
 *         WHY IT IS NEEDED AT ALL. The P-aligned average a panel draws was
 *         built by alignTemplatesFromCache, whose target column came from
 *         make_anchor_locator run on ONE reference beat -- or, when that
 *         failed, on the bin's whole-bin median. Its own comment records that
 *         compute_p_begin returned -1 on roughly half the bins, and the
 *         [anchors] summary line counts how many bins ended up NOT ALIGNED
 *         (byte-identical to R). So the P alignment on screen is frequently
 *         anchored on a column nobody would have chosen, and until now the
 *         operator could move the BAR but not the STACK: the bar moved, the
 *         waveform under it did not, and the two disagreed about where P was.
 *
 *         THE SAME SHAPE AS ppg_realign, AND THE SAME LIMITS. The cohort is
 *         not re-selected and the correlation filter is not re-run: the beats
 *         averaged are exactly the ones the build assigned to this slot. A
 *         re-filter would change group membership, which changes the
 *         morphology split, which is the pipeline's decision and not a mouse
 *         release's. Say so out loud rather than hide it.
 *
 *         ONE NUMBER CANNOT BE AN ALIGNMENT. The bar is a single column on the
 *         AVERAGE. Applying it as a rigid offset to every beat would translate
 *         the template and change nothing about its shape, which is not what
 *         the gesture means. So the column is a TARGET, and each beat
 *         contributes its own P onset, found by the build's own locator.
 *
 *         THE ROWS ARE R-FRAMED. _beats.bin holds the ECG slices as
 *         extract_beats_and_align produced them, R at a fixed column -- NOT
 *         the per-anchor matrices, which alignTemplatesFromCache builds on the
 *         fly and discards. So a bar stored in the P_ONSET frame has to be
 *         converted before it can index a row: see TemplateBin::frameShift,
 *         and the caller does that conversion, not this header.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */

#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "fiducial_marker_finding/feature_marks.hpp"
#include "template_generation/ppg_realign.hpp"

namespace ecg_realign {

    // ONE READER FOR THE FILE, shared with the pulse path. ppg_realign::loadBin
    // already takes a channel name and walks the same layout -- "CH1", "CH2",
    // "CH3" are blocks in the same file as "PPG", and its local-row-space join
    // note applies verbatim to them. A second traversal here would be a second
    // place the format is assumed.
    using ppg_realign::BinBeats;

    inline const char* channelKey(int lead) {
        switch (lead) {
        case 0: return "CH1";
        case 1: return "CH2";
        case 2: return "CH3";
        }
        return nullptr;
    }

    struct Result {
        bool ok = false;
        std::string why;                 // set when !ok, for the status bar
        std::vector<double> tmpl;        // re-medianed waveform
        std::vector<double> iqr;         // per-column SD (ddof=1), as the build
        int anchor_col = -1;             // where every beat's P onset now sits
        int n_members = 0;               // cohort offered
        int n_used = 0;                  // cohort whose landmark was findable
        int max_shift = 0;               // largest |shift| applied, samples
        double median_shift = 0.0;       // median signed shift, samples
    };

    // `rows` is the member list IN THIS CHANNEL'S LOCAL ROW SPACE -- pass
    // members_clean when non-empty, since that is the set the waveform on
    // screen was averaged over and the point is to re-average THE SAME SET.
    //
    // targetCol is the operator's column ALREADY CONVERTED TO THE ROW FRAME
    // (R-framed). rCol and fs are what make_anchor_locator needs; anchor is
    // normally AnchorType::P_ONSET but nothing here is specific to P, so the
    // other three bars can use this the day they want to.
    //
    // search_halfwin bounds how far from the operator's column a beat's own
    // landmark may be found. Unlike the pulse foot there is no trough to lock
    // onto: the locator returns a column or -1, so the window is a sanity
    // fence rather than a search hint -- a beat whose P onset comes back 300 ms
    // from where the operator says P is has not been corrected, it has been
    // mis-detected, and averaging it in would smear the result.
    inline Result realignAt(const BinBeats& beats,
        const std::vector<uint32_t>& rows,
        double targetCol,
        int rCol,
        double fs,
        AnchorType anchor,
        int search_halfwin)
    {
        Result out;
        out.n_members = static_cast<int>(rows.size());
        if (beats.empty()) { out.why = "no beat matrix for this bin/lead"; return out; }
        const int w = beats.width;
        const int target = static_cast<int>(std::lround(targetCol));
        if (target < 0 || target >= w) {
            out.why = "bar outside the beat window";
            return out;
        }
        if (!(fs > 0.0)) { out.why = "no sample rate"; return out; }
        if (search_halfwin < 1) search_halfwin = 1;
        out.anchor_col = target;

        // THE BUILD'S OWN LOCATOR, not a private P finder. alignTemplatesFromCache
        // locates every beat with this; a second detector here would mean the
        // re-stacked average and the one it replaces were anchored by two
        // different definitions of the same landmark.
        AnchorLocator locate = make_anchor_locator(anchor, rCol, fs);

        const double kNaN = std::numeric_limits<double>::quiet_NaN();
        std::vector<std::vector<double>> shifted;
        std::vector<int> shifts;
        shifted.reserve(rows.size());
        shifts.reserve(rows.size());

        // ---- PASS 1: LOCATE, ACROSS CORES -------------------------------
        //
        // THIS IS THE COST OF THE WHOLE GESTURE. locate() is the real landmark
        // detector -- compute_p_begin and friends -- run once per member beat,
        // and a bin carries on the order of a thousand of them. Serially that
        // is seconds of frozen UI on every mouse-up, which is what made moving
        // the bar feel slow even though nothing runs during the drag itself.
        //
        // SAFE TO PARALLELISE BECAUSE IT IS PURE: each iteration reads one row
        // and writes one slot of `own`, and AnchorLocator carries no state
        // across calls (make_anchor_locator closes over r_col and fs, both
        // by value). -1 marks a row to skip, decided here and acted on in the
        // serial pass below -- so the surviving ROW ORDER does not depend on
        // thread scheduling, which a push_back from inside the loop would.
        //
        // alignTemplatesFromCache parallelises the same work one level up, over
        // bins; it cannot help here because this is one bin.
        const int nRows = static_cast<int>(rows.size());
        std::vector<int> own(static_cast<size_t>(nRows), -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(std::min(8, std::max(1, nRows)))
#endif
        for (int i = 0; i < nRows; ++i) {
            const uint32_t r = rows[static_cast<size_t>(i)];
            if (r >= beats.rows.size()) continue;       // stale member list
            const double c = locate(beats.rows[r]);
            // NOT FOUND, OR FOUND SOMEWHERE ELSE ENTIRELY: left at -1 and
            // skipped, not forced onto the target. A P wave that genuinely is
            // not there (the case BankTemplate's own comment is about) must
            // drop out of the average rather than contribute a guessed column.
            if (!(c >= 0.0)) continue;
            const int col_i = static_cast<int>(std::lround(c));
            if (std::abs(col_i - target) > search_halfwin) continue;
            own[static_cast<size_t>(i)] = col_i;
        }

        // ---- PASS 2: SHIFT, SERIALLY ------------------------------------
        // Pure memory movement, and it appends to shared vectors -- cheap
        // enough that parallelising it would cost more in synchronisation than
        // it saved.
        for (int i = 0; i < nRows; ++i) {
            if (own[static_cast<size_t>(i)] < 0) continue;
            const std::vector<double>& src = beats.rows[rows[static_cast<size_t>(i)]];
            const int shift = target - own[static_cast<size_t>(i)];
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
        // TWO IS NOT A TEMPLATE -- refusing leaves the built average in place,
        // which beats replacing a real waveform with noise and reporting
        // success. Same floor as the pulse path.
        if (out.n_used < 3) {
            out.why = "fewer than 3 beats have that landmark near the bar";
            return out;
        }
        {
            std::vector<int> s = shifts;
            std::sort(s.begin(), s.end());
            out.median_shift = s[s.size() / 2];
        }

        // MEDIAN PLUS A ddof=1 SD, NaN-skipping per column -- exactly the
        // statistic alignTemplatesFromCache computes for a per-slot average,
        // despite the field being named _iqr. A different average here would
        // make the re-stacked slot incomparable with its siblings.
        out.tmpl.assign(static_cast<size_t>(w), kNaN);
        out.iqr.assign(static_cast<size_t>(w), kNaN);
        std::vector<double> col;
        col.reserve(shifted.size());
        for (int c = 0; c < w; ++c) {
            col.clear();
            for (const auto& b : shifted)
                if (!std::isnan(b[static_cast<size_t>(c)]))
                    col.push_back(b[static_cast<size_t>(c)]);
            // THE SAME FAR-TAIL GATE THE PULSE PATH USES, and for the reason
            // build_bins.hpp spells out for this very channel: the matrix is
            // padded to the bin's longest slice, so its tail is populated by
            // only the few longest beats -- a dropped R detection being the
            // classic case, whose "RR" spans several cardiac cycles, and which
            // is then the ONLY row with samples out there. A median over three
            // rows out of a thousand is not an average and the SD of those
            // three is not a corridor, but both come out finite and get drawn.
            if ((int)col.size() < ppg_realign::minColumnRows(out.n_used)) continue;
            std::sort(col.begin(), col.end());
            const size_t mid = col.size() / 2;
            out.tmpl[static_cast<size_t>(c)] = (col.size() % 2 == 0)
                ? 0.5 * (col[mid - 1] + col[mid]) : col[mid];
            if (col.size() >= 2) {
                double mu = 0.0;
                for (double x : col) mu += x;
                mu /= static_cast<double>(col.size());
                double ss = 0.0;
                for (double x : col) ss += (x - mu) * (x - mu);
                out.iqr[static_cast<size_t>(c)] =
                    std::sqrt(ss / static_cast<double>(col.size() - 1));
            }
        }

        out.ok = true;
        return out;
    }

    // 80 ms either side of the operator's column. Wide enough for the
    // beat-to-beat scatter of a real P onset at any plausible rate, narrow
    // enough that it cannot reach into the QRS on one side or the preceding
    // T wave on the other -- the two places a failed P detection lands.
    inline constexpr double kDefaultFenceSeconds = 0.080;

    inline int fenceHalfWinFor(double ecgRateHz) {
        if (!(ecgRateHz > 0.0)) return 20;
        return std::max(4, static_cast<int>(
            std::lround(kDefaultFenceSeconds * ecgRateHz)));
    }

}   // namespace ecg_realign