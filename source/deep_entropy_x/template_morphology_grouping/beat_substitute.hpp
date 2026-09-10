#pragma once
/**
 * @file   beat_substitute.hpp
 * @brief  when a beat is premature (ectopic) or voted ectopic due to being surrounded by ectopic beats, it cannot
 *         be left as a hole, thus it is substituted with a blend of the previous and subsequent beats
 */

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "template_bank.hpp"

namespace beat_substitute {
    inline constexpr double ewma_alpha = 0.125; //1/8

    struct SubstitutionResult {
        std::vector<std::pair<uint32_t, std::vector<double>>> replacements;
        uint32_t n_no_average = 0;
    };

    inline std::vector<double> substituteBeat(const std::vector<double>& avgOld, const std::vector<double>& current) {
        //substitute premature beat (current) with the exponentially weighted moving average of the average old beat and the current beat
        std::vector<double> out(avgOld.size());
        for (size_t j = 0; j < out.size(); ++j)
            out[j] = (1 - ewma_alpha) * avgOld[j] + ewma_alpha * current[j];
        return out;   // preserves temporal continuity
    }


    inline std::vector<double> substituteBeatNaNSafe( const std::vector<double>& avgOld, const std::vector<double>& current)
    {
        //wraps above but NAN safe
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        const size_t n = std::max(avgOld.size(), current.size());
        std::vector<double> a(n, NaN), c(n, NaN);
        for (size_t j = 0; j < n; ++j) {
            const double av = (j < avgOld.size()) ? avgOld[j] : NaN;
            const double cv = (j < current.size()) ? current[j] : NaN;
            const bool na = std::isnan(av), nc = std::isnan(cv);
            if (na && nc) continue;          // both NaN -> blend yields NaN
            a[j] = na ? cv : av;
            c[j] = nc ? av : cv;
        }
        return substituteBeat(a, c);
    }

    inline SubstitutionResult run(const std::vector<std::vector<double>>& beats, const std::vector<int32_t>& assignment,  const std::vector<char>& premature, const tbank::TemplateBank& bank, std::vector<tbank::BeatFlags>* flags = nullptr)
    {
        /* Walks the beats in order, substituting the premature ones against the
        running average of the template each is assigned to.*/
        SubstitutionResult out;

        // Seed each template's running average from its median. Starting from
        // the template rather than from the first member means the first
        // substituted beat is blended against the whole population, not against
        // whichever beat happened to arrive first.
        std::vector<std::vector<double>> avg(bank.size());
        for (int t = 0; t < bank.size(); ++t) avg[t] = bank.templates[t].tmpl;

        const size_t n = std::min(beats.size(), assignment.size());
        for (uint32_t i = 0; i < n; ++i) {
            const int32_t t = assignment[i];
            if (t < 0 || t >= static_cast<int32_t>(avg.size())) continue;
            if (i >= premature.size() || !premature[i]) continue;

            // The one guard the spec implies: substituteBeat reads avgOld[j],
            // so there has to be an avgOld. A template with no median has no
            // population to blend toward.
            if (avg[t].empty()) { ++out.n_no_average; continue; }

            std::vector<double> blended = substituteBeatNaNSafe(avg[t], beats[i]);
            avg[t] = blended;                      // recursion, not a fixed offset
            out.replacements.emplace_back(i, std::move(blended));

            if (flags && i < flags->size()) (*flags)[i].substituted = true;
        }
        return out;
    }
}  // namespace beat_substitute
