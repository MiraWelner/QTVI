#pragma once
/**
 * @file   beat_substitute.hpp
 * @brief  Section 4.6 beat substitution: EWMA with alpha = 1/8.
 *
 *         SPEC FIDELITY. substituteBeat is transcribed verbatim and is the
 *         only arithmetic here. Everything else answers the two questions the
 *         spec leaves to the caller -- WHICH beat gets substituted, and what
 *         `avgOld` is.
 *
 *         WHAT SUBSTITUTION IS FOR: TEMPORAL CONTINUITY. A beat flagged by the
 *         prematurity filter or by the 5-of-8 vote cannot go into the average,
 *         but deleting it leaves a hole in the sequence. So it is REPLACED by
 *         a blend of the nearest clean beat either side of it, keeping 1/8 of
 *         its own samples. The row stays, the shape becomes its neighbours',
 *         and the trace does not step.
 *
 *           avgOld  = the mean of the nearest clean beat on the left and the
 *                     nearest clean beat on the right
 *           current = the flagged beat itself
 *           out     = 0.875 * avgOld + 0.125 * current
 *
 *         WHAT IT IS NOT FOR, and what was removed to get here:
 *
 *           - NOT a running reference. The previous SubstitutionState fed
 *             every clean beat through substituteBeat to maintain an EWMA
 *             average. That made alpha = 1/8 a forgetting factor on the
 *             reference, which is a different job with a different meaning,
 *             and it is gone. A reference average, if one is wanted, is a mean
 *             over the clean beats -- not this function.
 *           - NOT quality repair. The previous driver substituted beats whose
 *             composite SQI fell in the 0.50-0.70 band. SQI handling is
 *             Section 4.3's business; adjacency to a rhythm flag is 4.6's, and
 *             a mid-SQI sinus beat next to nothing at all has no continuity
 *             problem to fix.
 *
 *         NEIGHBOURS ARE NEVER SUBSTITUTED, AND NEVER CHAINED. Only flagged
 *         beats are replaced. The beats used as source material are read from
 *         the ORIGINAL matrix and must themselves be unflagged, so a run of
 *         ectopy cannot bootstrap itself: in bigeminy, every other beat is
 *         flagged and each one reaches past its flagged neighbours to real
 *         sinus beats on both sides. If substitutions were allowed to feed
 *         later substitutions, one flagged beat's blend would propagate along
 *         the whole run and the output would converge on a single shape.
 *
 *         maxSearch BOUNDS THAT REACH. Inside a long VT run the nearest clean
 *         beat may be many beats away, and a blend against a beat ten seconds
 *         distant is not continuity, it is fabrication. Past maxSearch the
 *         beat is left unsubstituted and reported as such -- an honest gap
 *         beats a manufactured beat.
 *
 * @date   2026-08-25
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

 // Section 4.6, EWMA with alpha = 1/8 -- VERBATIM.
inline std::vector<double> substituteBeat(const std::vector<double>& avgOld,
    const std::vector<double>& current, double alpha = 0.125)
{
    std::vector<double> out(avgOld.size());
    for (size_t j = 0; j < out.size(); ++j) out[j] = (1 - alpha) * avgOld[j] + alpha * current[j];
    return out;                                  // preserves temporal continuity
}

namespace beatsub {

    inline constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    inline constexpr double kAlpha = 0.125;      // 1/8, as specified

    struct Params {
        double alpha = kAlpha;
        // How far to look for a clean beat on each side. 8 covers a couplet,
        // a triplet and short bigeminy; a long VT run deliberately exceeds it.
        int    maxSearch = 8;
        // Accept a one-sided blend when only one side has a clean beat within
        // reach (start of a bin, end of a bin, a run against the edge).
        bool   allowOneSided = true;
    };

    // One substituted beat, and the beats that made it. This is what
    // {stem}_substituted.csv reports.
    struct Substitution {
        int  beat = -1;   // index in the bin: the beat replaced
        int  leftBeat = -1;   // index of the clean beat on the left, or -1
        int  rightBeat = -1;   // index of the clean beat on the right, or -1
        std::vector<double> avgOld;     // the neighbour mean fed to substituteBeat
        std::vector<double> original;   // the flagged beat as it was
        std::vector<double> blend;      // what replaces it
        bool substituted = false;
        const char* note = "";
    };

    struct BinResult {
        std::vector<Substitution> subs;   // one per substituted beat
        int nFlagged = 0, nSubstituted = 0, nSkipped = 0;
    };

    // A row usable as neighbour material or as `current`.
    inline bool wellFormed(const std::vector<double>& beat, std::size_t width) {
        if (beat.empty() || beat.size() != width) return false;
        for (double v : beat) if (!std::isfinite(v)) return false;
        return true;
    }

    // Element-wise mean of the two neighbours. NaN-aware per column: where one
    // side is missing a sample the other carries it alone, because the aligned
    // matrix is NaN-padded at the edges and discarding the column outright
    // would put a hole in the middle of the blend.
    inline std::vector<double> neighbourMean(const std::vector<double>* left,
        const std::vector<double>* right,
        std::size_t width)
    {
        std::vector<double> out(width, kNaN);
        for (std::size_t c = 0; c < width; ++c) {
            const double a = (left && c < left->size()) ? (*left)[c] : kNaN;
            const double b = (right && c < right->size()) ? (*right)[c] : kNaN;
            const bool fa = std::isfinite(a), fb = std::isfinite(b);
            if (fa && fb) out[c] = 0.5 * (a + b);
            else if (fa)  out[c] = a;
            else if (fb)  out[c] = b;
        }
        return out;
    }

    // `flagged[i]` is the rhythm flag: premature OR voted premature (4.6).
    // Substitution is driven by that and by nothing else. Neighbour material
    // is any unflagged beat.
    inline BinResult applyToBin(const std::vector<std::vector<double>>& beats,
        const std::vector<char>& flagged,
        const Params& p = {})
    {
        BinResult out;
        const std::size_t n = beats.size();
        if (n == 0 || flagged.size() != n) return out;

        // Width comes from the first well-formed beat: the aligned matrix
        // shares one axis, so a row of a different length is a bad row, not a
        // different geometry.
        std::size_t width = 0;
        for (const std::vector<double>& b : beats)
            if (!b.empty()) { width = b.size(); break; }
        if (width == 0) return out;

        const auto usableSource = [&](std::size_t k) { return flagged[k] == 0; };

        for (std::size_t i = 0; i < n; ++i) {
            if (!flagged[i]) continue;
            ++out.nFlagged;

            Substitution s;
            s.beat = static_cast<int>(i);

            if (!wellFormed(beats[i], width)) {
                ++out.nSkipped;
                continue;
            }

            // Nearest UNFLAGGED, well-formed beat each side, from the original
            // matrix. Never a substituted one -- see the header note on runs.
            const std::vector<double>* left = nullptr;
            const std::vector<double>* right = nullptr;
            for (int k = 1; k <= p.maxSearch; ++k) {
                const long li = static_cast<long>(i) - k;
                if (!left && li >= 0
                    && usableSource(static_cast<std::size_t>(li))
                    && wellFormed(beats[static_cast<std::size_t>(li)], width)) {
                    left = &beats[static_cast<std::size_t>(li)];
                    s.leftBeat = static_cast<int>(li);
                }
                const std::size_t ri = i + static_cast<std::size_t>(k);
                if (!right && ri < n && usableSource(ri)
                    && wellFormed(beats[ri], width)) {
                    right = &beats[ri];
                    s.rightBeat = static_cast<int>(ri);
                }
                if (left && right) break;
            }

            if (!left && !right) { ++out.nSkipped; continue; }
            if ((!left || !right) && !p.allowOneSided) {
                ++out.nSkipped;
                continue;
            }

            s.avgOld = neighbourMean(left, right, width);
            s.original = beats[i];
            s.blend = substituteBeat(s.avgOld, s.original, p.alpha);
            s.substituted = true;
            s.note = (left && right) ? "both_sides" : "one_sided";

            ++out.nSubstituted;
            out.subs.push_back(std::move(s));
        }
        return out;
    }

} // namespace beatsub