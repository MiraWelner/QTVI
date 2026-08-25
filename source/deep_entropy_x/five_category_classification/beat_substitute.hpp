/**
 * @file   beat_substitute.hpp
 * @brief  Beat substitution by exponentially weighted moving average,
 *         alpha = 1/8. Spec Section 4.6.
 *
 *         SPEC FIDELITY. substituteBeat is transcribed verbatim. Its
 *         preconditions are therefore the spec's preconditions and are stated
 *         below rather than defended against in code:
 *
 *           - avgOld and current must be the SAME LENGTH. The output is sized
 *             from avgOld and indexes current at the same j, so a shorter
 *             `current` reads past its end. Callers must pass equal-length
 *             beats; the bin beat matrices out of Phase 1 are uniform width,
 *             so this holds for the pipeline path.
 *           - avgOld must be non-empty. An empty running average returns an
 *             empty beat.
 *           - Both must be free of non-finite samples. EWMA has infinite
 *             memory, so one NaN entering a column makes that column NaN for
 *             the remainder of the record.
 *
 *         applyToBin below enforces those preconditions at the boundary --
 *         seeding on the first clean beat and skipping ill-formed rows --
 *         without altering substituteBeat itself.
 *
 * @date   2026-08-24
 */
#pragma once

#include "five_categories.hpp"   // beatcls::BeatVerdict

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

    // Per-bin state: the running average that substituted beats are drawn
    // from, and the counts feeding the Task I max_correction_pct gate.
    struct SubstitutionState {
        std::vector<double> average;
        double alpha = kAlpha;
        int nIncluded = 0, nSubstituted = 0, nExcluded = 0;
        bool ready() const { return !average.empty() && nIncluded > 0; }
    };

    // Is this row a legal argument to substituteBeat against `average`?
    inline bool wellFormed(const std::vector<double>& beat,
        const std::vector<double>& average)
    {
        if (beat.empty() || beat.size() != average.size()) return false;
        for (double v : beat) if (!std::isfinite(v)) return false;
        return true;
    }

    // The running average is updated from INCLUDE-quality beats only. A
    // substituted beat is a read of the average, never a write to it --
    // otherwise the reference becomes an average of the beats it replaced.
    inline bool observeClean(SubstitutionState& st, const std::vector<double>& beat) {
        if (beat.empty()) return false;
        if (st.average.empty()) {                    // seed
            bool ok = true;
            for (double v : beat) if (!std::isfinite(v)) ok = false;
            if (!ok) return false;
            st.average = beat;
            ++st.nIncluded;
            return true;
        }
        if (!wellFormed(beat, st.average)) return false;
        st.average = substituteBeat(st.average, beat, st.alpha);
        ++st.nIncluded;
        return true;
    }

    struct SubstitutionResult {
        std::vector<double> beat;
        bool substituted = false;
        bool dropped = false;
        const char* note = "";
    };

    inline SubstitutionResult substitute(SubstitutionState& st,
        const std::vector<double>& beat)
    {
        SubstitutionResult r;
        if (!st.ready()) {
            r.dropped = true; r.note = "no_reference_yet"; ++st.nExcluded; return r;
        }
        if (!wellFormed(beat, st.average)) {
            r.dropped = true; r.note = "ill_formed_beat"; ++st.nExcluded; return r;
        }
        r.beat = substituteBeat(st.average, beat, st.alpha);
        r.substituted = true;
        r.note = "ewma_blend";
        ++st.nSubstituted;
        return r;
    }

    inline void observeExcluded(SubstitutionState& st) { ++st.nExcluded; }

    struct SubstitutionSummary {
        int nBeats = 0, nIncluded = 0, nSubstituted = 0, nExcluded = 0;
        double substitutedPct = kNaN;
        bool overCorrectionFlag = false;    // Task I max_correction_pct
    };

    inline SubstitutionSummary summarize(const SubstitutionState& st,
        double maxPct = 1.0)
    {
        SubstitutionSummary s;
        s.nIncluded = st.nIncluded;
        s.nSubstituted = st.nSubstituted;
        s.nExcluded = st.nExcluded;
        s.nBeats = st.nIncluded + st.nSubstituted + st.nExcluded;
        if (s.nBeats > 0) {
            s.substitutedPct = 100.0 * st.nSubstituted / s.nBeats;
            s.overCorrectionFlag = (s.substitutedPct > maxPct);
        }
        return s;
    }

    // One pass over a bin. Dropped beats come back as empty rows so the row
    // index still equals the beat index.
    struct BinSubstitutionResult {
        std::vector<std::vector<double>> beats;
        std::vector<char> wasSubstituted, wasDropped;
        SubstitutionSummary summary;
    };

    inline BinSubstitutionResult applyToBin(const std::vector<std::vector<double>>& beats,
        const std::vector<beatcls::BeatVerdict>& verdicts,
        double alpha = kAlpha,
        double maxCorrectionPct = 1.0)
    {
        BinSubstitutionResult out;
        const std::size_t n = std::min(beats.size(), verdicts.size());
        out.beats.resize(beats.size());
        out.wasSubstituted.assign(beats.size(), 0);
        out.wasDropped.assign(beats.size(), 0);

        SubstitutionState st;
        st.alpha = alpha;

        for (std::size_t i = 0; i < n; ++i) {
            switch (verdicts[i].handling) {
            case beatcls::BeatVerdict::INCLUDE:
                if (observeClean(st, beats[i])) out.beats[i] = beats[i];
                else                            out.wasDropped[i] = 1;
                break;
            case beatcls::BeatVerdict::SUBSTITUTE: {
                const SubstitutionResult r = substitute(st, beats[i]);
                if (r.substituted) { out.beats[i] = r.beat; out.wasSubstituted[i] = 1; }
                else { out.wasDropped[i] = 1; }
                break;
            }
            default:
                observeExcluded(st);
                out.wasDropped[i] = 1;
                break;
            }
        }
        out.summary = summarize(st, maxCorrectionPct);
        return out;
    }

} // namespace beatsub