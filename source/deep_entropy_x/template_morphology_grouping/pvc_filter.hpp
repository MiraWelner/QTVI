#pragma once
/**
 * @file   pvc_filter.hpp
 * @brief  Section 4.6 prematurity filter and 5-of-8 voting, run across ALL
 *         detected beats as a standalone step.
 *
 *         THE FUNCTIONS ARE UNMODIFIED. The spec's isPremature() indexes rr[t]
 *         as beat t's OWN interval -- the interval PRECEDING beat t -- and takes
 *         its median over rr[t-10 .. t-1]. Your data is not stored that way:
 *         alignment.hpp's rr_lens[i] holds the interval FOLLOWING beat i
 *         (R[i+1] - R[i]), because that is the span the slice covers.
 *
 *         Rather than rewrite the spec function to chase the offset, the INPUT
 *         is converted once by toPrecedingIntervals(). Two reasons that is the
 *         right way round: the function stays byte-identical to the document so
 *         it can be diffed against it, and the convention shift lives in exactly
 *         one place instead of being smeared through index arithmetic at every
 *         use.
 *
 *         With the adapter in place this reproduces alignment.hpp's internal
 *         filter exactly: its median over rr_lens[t-11 .. t-2] compared against
 *         rr_lens[t-1] is, in preceding-interval terms, the median over
 *         rr_pre[t-10 .. t-1] against rr_pre[t] -- which is what the spec's
 *         function computes.
 *
 *         WHY THIS EXISTS SEPARATELY FROM alignment.hpp. alignment.hpp already
 *         implements both tests (lines ~247-273), correctly and in the right
 *         place: after slicing, before the Tukey passes, because Tukey rejects
 *         on RR length at 1.5*IQR and a premature beat is short by definition,
 *         so flagging after pruning would find nothing left to flag. But those
 *         flags live on `out.premature` / `out.voted`, which are compacted by
 *         every apply_mask alongside `beats` -- so a beat that Tukey rejects
 *         loses its flag with it.
 *
 *         The pipeline order requires the filter to run across all beats, and
 *         the archive requires three flags for EVERY beat, including the ones
 *         alignment discarded. A beat rejected by Tukey and premature is one of
 *         the more interesting rows in the file: it is the case where the two
 *         gates agree that something is wrong, and it is invisible if the flag
 *         was compacted away with the beat.
 *
 *         So this runs on the detected R-peak list before any pruning, keyed on
 *         the ORIGINAL beat index, and its output is never compacted. It is not
 *         a reimplementation competing with alignment's copy -- alignment needs
 *         its own internal flags to drive the line-175 pruning exemption, and
 *         those must agree with these. runFilter() is written to produce
 *         identical verdicts on the beats the two share; see kVoteWindow.
 *
 *         TWO PATHS, DISJOINT COMPETENCE, COUNTED SEPARATELY. The raw
 *         prematurity test catches isolated ectopy. The vote only rescues beats
 *         inside consecutive runs. Alternating patterns fall between them and
 *         neither works:
 *
 *           - Perfect bigeminy: exactly 4 of any 8 consecutive beats are
 *             ectopic, the count never reaches 5, THE VOTE NEVER FIRES.
 *           - Trigeminy: ~2.7 of 8. Never fires.
 *           - VT run / salvo: the neighborhood is densely flagged, the vote
 *             fires, and it is the only thing that catches the mid-run beats
 *             whose trailing-ten median has already collapsed.
 *
 *         Which means the COUNT PATTERN IDENTIFIES THE RHYTHM. High raw with
 *         zero vote is bigeminy. High vote is a run. That is why the two are
 *         never summed into one number.
 *
 *         AND THE RAW TEST DECAYS INSIDE A RUN. RR(t) < 0.80 * median of the
 *         trailing ten stops firing once the trailing ten are themselves short
 *         -- by beat 11 of a VT run the ratio approaches 1.0 and the test goes
 *         quiet on exactly the beats that matter most. The vote is the patch
 *         for that, and it is a timing-domain patch for a timing-domain blind
 *         spot: it votes over prematurity flags and never looks at a waveform.
 *         Morphology has no equivalent decay, which is why the bank and this
 *         filter are independent gates rather than one gate twice.
 */

#include <algorithm>
#include <cstdint>
#include <vector>

#include "template_bank.hpp"   // tbank::PvcFilter

namespace pvc_filter {

    // ---------------------------------------------------------------------
    // Section 4.6, VERBATIM from the specification. Do not edit these two to
    // fix a convention mismatch -- adapt the input instead, via
    // toPrecedingIntervals() below.
    // ---------------------------------------------------------------------

    inline bool isPremature(const std::vector<double>& rr, int t) { // rr in ms
        if (t < 10) return false;
        std::vector<double> w(rr.begin() + t - 10, rr.begin() + t);
        std::sort(w.begin(), w.end());
        double med = w[w.size() / 2];
        return rr[t] < 0.80 * med;   // RR(t) < 0.80 * median
    }

    // 5-of-8 voting: flag beat t if >=5 of the surrounding 8 beats are flagged
    inline bool voteFlag(const std::vector<char>& flag, int t) {
        int lo = std::max(0, t - 4), hi = std::min((int)flag.size(), t + 4), c = 0;
        for (int i = lo; i < hi; ++i) c += flag[i];
        return c >= 5;
    }

    // ---------------------------------------------------------------------
    // Convention adapter
    // ---------------------------------------------------------------------

    // rr_after[i] = R[i+1] - R[i], alignment.hpp's rr_lens convention.
    // Returns rr_pre[i] = the interval PRECEDING beat i, which is what the
    // spec's isPremature() reads at index t.
    //
    // A beat is premature when the interval BEFORE it is short. Feeding rr_lens
    // straight in tests the interval AFTER each beat, which flags the beat
    // PRECEDING every PVC instead of the PVC -- and since the beat before a PVC
    // looks normal while the beat after is a compensatory pause, that detector
    // fires at roughly the right rate in roughly the right places. It looks
    // exactly like a detector that nearly works.
    //
    // rr_pre[0] has no predecessor and is NaN, not 0.0, so nothing can mistake
    // "no preceding beat" for "a zero-length interval". isPremature() returns
    // false for t < 10 anyway, so it is never read there.
    inline std::vector<double> toPrecedingIntervals(
        const std::vector<double>& rr_after)
    {
        const size_t n = rr_after.size();
        std::vector<double> rr_pre(n, std::numeric_limits<double>::quiet_NaN());
        for (size_t i = 1; i < n; ++i) rr_pre[i] = rr_after[i - 1];
        return rr_pre;
    }

    // Section 4.6 constants, for reporting only -- the values live inside the
    // verbatim functions above. Named here so counts and comments can refer to
    // them without a second definition drifting from the first.
    inline constexpr double kPrematureRatio = 0.80;
    inline constexpr int    kTrailingWindow = 10;
    inline constexpr int    kVoteWindow     = 8;
    inline constexpr int    kVoteRequired   = 5;
    inline constexpr size_t kMinBeatsForFilter = 12;

    struct FilterResult {
        // Parallel to the detected beat list, never compacted.
        std::vector<tbank::PvcFilter> verdict;

        // The ratio the test compared, per beat: the continuous form of a binary
        // decision. A run of beats sitting at 0.82 is a different situation from
        // a run at 0.40, and only the ratio shows it. NaN where no trailing
        // median was available.
        std::vector<double> ratio;

        uint32_t n_premature  = 0;
        uint32_t n_vote_only  = 0;   // voted but not itself premature
        uint32_t n_no_median  = 0;   // t < 10, no trailing window
        uint32_t n_vote_blind = 0;   // vote window clamped below 5 beats
    };

    // rr_after: alignment.hpp's rr_lens convention, in milliseconds.
    inline FilterResult runFilter(const std::vector<double>& rr_after) {
        const size_t n = rr_after.size();
        FilterResult out;
        out.verdict.assign(n, tbank::PvcFilter::NONE);
        out.ratio.assign(n, std::numeric_limits<double>::quiet_NaN());
        if (n == 0) return out;
        if (n < kMinBeatsForFilter) {
            // alignment.hpp requires nb >= 12 before running either test.
            // Matched rather than relaxed: a filter that fires on short records
            // here and not there would make the two disagree exactly where a bin
            // is too short to trust anyway.
            out.n_no_median = static_cast<uint32_t>(n);
            return out;
        }

        const std::vector<double> rr = toPrecedingIntervals(rr_after);

        // ---- pass A: the spec's isPremature(), per beat -------------------
        std::vector<char> flag(n, 0);
        for (size_t t = 0; t < n; ++t) {
            if (t < static_cast<size_t>(kTrailingWindow)) { ++out.n_no_median; continue; }
            if (std::isnan(rr[t])) continue;

            // Ratio recorded over the same window the spec's function uses, so
            // the two cannot disagree about which beats were in the median.
            std::vector<double> w(rr.begin() + t - kTrailingWindow, rr.begin() + t);
            w.erase(std::remove_if(w.begin(), w.end(),
                [](double v) { return std::isnan(v); }), w.end());
            if (!w.empty()) {
                std::sort(w.begin(), w.end());
                const double med = w[w.size() / 2];
                if (med > 0.0) out.ratio[t] = rr[t] / med;
            }

            if (isPremature(rr, static_cast<int>(t))) flag[t] = 1;
        }

        // ---- pass B: the spec's voteFlag(), over pass A's flags -----------
        // A SECOND PASS over the completed array, never interleaved: a vote that
        // could see votes would propagate a flag along the whole record from a
        // single dense neighbourhood.
        for (size_t t = 0; t < n; ++t) {
            const long lo = std::max<long>(0, static_cast<long>(t) - 4);
            const long hi = std::min<long>(static_cast<long>(n),
                static_cast<long>(t) + 4);
            // Near the record edges the spec's own clamping leaves fewer than
            // kVoteRequired beats in the window, so the vote cannot fire there
            // however ectopic the neighbourhood is. Preserved rather than
            // "fixed" -- alignment.hpp clamps identically, and a silently
            // different edge rule between the two would be worse than a known
            // blind spot -- but counted, so the blind region shows up in the
            // archive instead of being discovered later.
            if (hi - lo < kVoteRequired) ++out.n_vote_blind;

            if (flag[t]) {
                out.verdict[t] = tbank::PvcFilter::PREMATURE;
                ++out.n_premature;
            }
            else if (voteFlag(flag, static_cast<int>(t))) {
                // Premature wins over voted: direct evidence over inferred.
                // Same precedence as create_ecg_templates.hpp lines 170-172.
                out.verdict[t] = tbank::PvcFilter::VOTE;
                ++out.n_vote_only;
            }
        }

        return out;
    }

    // Diagnostic for the bigeminy blind spot, read off the verdicts rather than
    // the rhythm. A substantial premature count with zero votes is alternating
    // ectopy -- the pattern where the vote is arithmetically incapable of
    // firing. A substantial vote count means runs.
    inline bool looksAlternating(const FilterResult& r) {
        return r.n_premature >= 8 && r.n_vote_only == 0;
    }

}  // namespace pvc_filter
