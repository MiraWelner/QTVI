#pragma once
/**
 * @file   state_stratified.hpp
 * @brief  Task J (state-stratified features) -- STUB.
 *
 *         Task H's record summary carries a StateStratifiedFeatures member, so
 *         the type and the entry point exist now; the stratification itself
 *         does not. stratify() returns available == false and a note, and the
 *         summary page prints that note in the panel Task J will fill.
 *
 *         WHEN TASK J LANDS, replace the body of stratify() and add fields to
 *         StateStratifiedFeatures. Keep `available` and `note`: the summary
 *         writer and the page renderer read only those two to decide whether
 *         to draw the panel, so neither has to change when the fields do.
 *
 *         The inputs are what a per-state split will need -- each template's
 *         member times and the record's staging -- passed now so the call site
 *         in record_summary.hpp does not have to change either.
 */

#include <string>
#include <vector>

#include "template_generation/template_structs.hpp"
#include "prep_for_peakfinding/beat_times.hpp"
#include "prep_for_peakfinding/record_sleep.hpp"

namespace state_stratified {

    struct StateStratifiedFeatures {
        bool available = false;
        std::string note = "Task J (state-stratified features) not implemented";
        // Task J: per-state feature tables go here.
    };

    // STUB. See the file comment.
    inline StateStratifiedFeatures stratify(
        const template_structs::TemplateFile& /*tmpl*/,
        const beat_times::BeatTimes* /*times*/,
        const record_sleep::RecordSleep* /*sleep*/)
    {
        return StateStratifiedFeatures{};
    }

}  // namespace state_stratified#pragma once
