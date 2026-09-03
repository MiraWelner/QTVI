#pragma once
/**
 * @file   beat_substitute.hpp
 * @brief  Section 4.6 beat substitution. EWMA with alpha = 1/8.
 *
 *         WHAT IT IS FOR. A rejected beat leaves a hole in every beat-to-beat
 *         series, and a hole is worse than a smoothed estimate for anything
 *         that measures variability: an RR series with gaps reports variance
 *         that depends on which beats were dropped. So a borderline beat is
 *         replaced by a blend of the running average and itself, preserving
 *         temporal continuity.
 *
 *         WHERE IT RUNS. AFTER bank assignment, never before. The spec's
 *         substituteBeat() blends against "avgOld", which has one meaning when
 *         a bin holds one template and six meanings when it holds a bank -- it
 *         must be the assigned template's OWN average, or a borderline PVC gets
 *         blended toward sinus and the substitution manufactures a morphology
 *         that never occurred.
 *
 *         WHAT MUST BE FLAGGED. A substituted beat is not an observation. Every
 *         one is recorded in BeatFlags::substituted so downstream variance,
 *         corridor, and IQR calculations can exclude them -- a corridor fitted
 *         partly to synthetic values understates its own spread, and nothing
 *         downstream could detect that from the numbers alone.
 *
 *         NOTE ON THE MEAN/MEDIAN MISMATCH. This is a mean recursion feeding a
 *         median pipeline: templates are column-wise NaN-skipping medians
 *         (create_ecg_templates.hpp), while an EWMA is a running mean. That is
 *         the spec's design and is kept, but it is why substituted beats are
 *         excluded from template recomputation by default -- injecting blended
 *         values into a median's input set makes the median partly a function
 *         of its own history.
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "template_bank.hpp"

namespace beat_substitute {

    // Section 4.6: alpha = 1/8.
    inline constexpr double kAlpha = 0.125;

    // A beat is "borderline" when it clears neither the reject floor nor the
    // accept ceiling: good enough to be worth keeping in the series, not good
    // enough to feed a reference. Below the reject floor a beat is not
    // borderline, it is bad, and substituting it would smuggle an artifact into
    // the series wearing the average's shape.
    inline constexpr double kBorderlineLo = 0.60;
    // The top of the band IS the 4.6 morphology threshold, so it tracks the
    // configured ECG floor rather than restating 0.85. Above the floor a beat
    // is a good match and needs no substitution; the band runs from the reject
    // floor up to it.
    inline double borderlineHi() { return tbank::matchFloorEcg(); }

    // ---------------------------------------------------------------------
    // Section 4.6, VERBATIM from the specification.
    // ---------------------------------------------------------------------

    inline std::vector<double> substituteBeat(const std::vector<double>& avgOld,
        const std::vector<double>& current, double alpha = 0.125) {
        std::vector<double> out(avgOld.size());
        for (size_t j = 0; j < out.size(); ++j)
            out[j] = (1 - alpha) * avgOld[j] + alpha * current[j];
        return out;   // preserves temporal continuity
    }

    // NaN-safe variant, used internally by run().
    //
    // The verbatim function above is correct arithmetic on clean vectors, but
    // the shared beat axis is NaN-PADDED at both ends -- alignment lays every
    // beat on one axis and fills the overhang with NaN. Multiply NaN by
    // (1-alpha) and it stays NaN, so the verbatim blend propagates whichever
    // input had MORE padding into the result, silently shortening the
    // substituted beat. It also reads current[j] past its end when the two
    // vectors differ in length, which they can when a beat sits near a bin edge.
    //
    // So: where one side is NaN the other passes through unchanged, and where
    // both are, the output is NaN. Everything else is the same expression.
    inline std::vector<double> substituteBeatNaNSafe(
        const std::vector<double>& avgOld,
        const std::vector<double>& current,
        double alpha = kAlpha)
    {
        const size_t n = std::max(avgOld.size(), current.size());
        std::vector<double> out(n, std::numeric_limits<double>::quiet_NaN());
        for (size_t j = 0; j < n; ++j) {
            const double a = (j < avgOld.size()) ? avgOld[j]
                : std::numeric_limits<double>::quiet_NaN();
            const double c = (j < current.size()) ? current[j]
                : std::numeric_limits<double>::quiet_NaN();
            const bool na = std::isnan(a), nc = std::isnan(c);
            if (na && nc) continue;
            if (na) { out[j] = c; continue; }
            if (nc) { out[j] = a; continue; }
            out[j] = (1.0 - alpha) * a + alpha * c;
        }
        return out;
    }

    struct SubstitutionResult {
        // Replacement waveforms, keyed by beat index. Only borderline beats
        // appear; everything else is left exactly as observed.
        std::vector<std::pair<uint32_t, std::vector<double>>> replacements;
        uint32_t n_substituted = 0;
        uint32_t n_too_bad = 0;     // below kBorderlineLo: rejected, not blended
        uint32_t n_no_average = 0;  // assigned template had no usable average
    };

    // Walks the beats in order, substituting the borderline ones against the
    // running average of the template each is assigned to.
    //
    // The running average is per template and is updated with the SUBSTITUTED
    // value, which is what makes this a recursion rather than a fixed offset --
    // consecutive borderline beats drift toward the average rather than each
    // being pulled the same distance from it.
    //
    // `score[i]` is the beat's correlation against its assigned template, i.e.
    // the same quantity the bank used to assign it.
    inline SubstitutionResult run(const std::vector<std::vector<double>>& beats,
        const std::vector<int32_t>& assignment,
        const std::vector<double>& score,
        const tbank::TemplateBank& bank,
        std::vector<tbank::BeatFlags>* flags = nullptr,
        double alpha = kAlpha,
        double lo = kBorderlineLo,
        double hi = borderlineHi())
    {
        SubstitutionResult out;

        // Seed each template's running average from its median. Starting from
        // the template rather than from the first member means the first
        // borderline beat is blended against the whole population, not against
        // whichever beat happened to arrive first.
        std::vector<std::vector<double>> avg(bank.size());
        for (int t = 0; t < bank.size(); ++t) avg[t] = bank.templates[t].tmpl;

        const size_t n = std::min(beats.size(), assignment.size());
        for (uint32_t i = 0; i < n; ++i) {
            const int32_t t = assignment[i];
            if (t < 0 || t >= static_cast<int32_t>(avg.size())) continue;
            if (i >= score.size() || std::isnan(score[i])) continue;

            if (score[i] >= hi) continue;          // good enough as observed
            if (score[i] < lo) { ++out.n_too_bad; continue; }
            if (avg[t].empty()) { ++out.n_no_average; continue; }

            std::vector<double> blended =
                substituteBeatNaNSafe(avg[t], beats[i], alpha);
            avg[t] = blended;                      // recursion, not a fixed offset
            out.replacements.emplace_back(i, std::move(blended));
            ++out.n_substituted;

            if (flags && i < flags->size()) (*flags)[i].substituted = true;
        }
        return out;
    }

    // A substituted beat is a blend, never a copy: it must differ from both the
    // average it was blended against and the beat it replaced. With alpha = 1/8
    // it sits one eighth of the way from the average toward the observation.
    // Exposed so the acceptance test can assert it rather than assume it.
    inline bool isBlendNotCopy(const std::vector<double>& avgOld,
        const std::vector<double>& current,
        const std::vector<double>& result,
        double tol = 1e-12)
    {
        bool differsFromAvg = false, differsFromCurrent = false;
        const size_t n = result.size();
        for (size_t j = 0; j < n; ++j) {
            if (std::isnan(result[j])) continue;
            if (j < avgOld.size() && !std::isnan(avgOld[j])
                && std::fabs(result[j] - avgOld[j]) > tol) differsFromAvg = true;
            if (j < current.size() && !std::isnan(current[j])
                && std::fabs(result[j] - current[j]) > tol) differsFromCurrent = true;
        }
        return differsFromAvg && differsFromCurrent;
    }

}  // namespace beat_substitute
