#pragma once
/**
 * @file   bin_pipeline.hpp
 * @brief  Per-bin driver enforcing the Section 4.5-4.6 order of operations.
 *
 *              1. Detect all beats.
 *              2. Classify every beat into the three categories. Remove
 *                 nothing.
 *              3. Run the PVC filter across all beats. Flag category 2.
 *              4. Build templates from category 1 only. Apply Tukey INSIDE
 *                 that group.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "template_bank.hpp"
#include "template_assign.hpp"
#include "pvc_filter.hpp"
#include "seed_pool.hpp"
#include "template_marking_gui/alignment.hpp"

namespace bin_pipeline {

    // ---------------------------------------------------------------------
    // Step 2: classification from operator marks
    // ---------------------------------------------------------------------
    // categoryForCode / postEligibleCode / categoryFromMarks MOVED to
    // tbank::categoriesFromMarks (template_bank.hpp), where the annotation code
    // constants live. They had no caller left once runChannel went, and the copy
    // here spelled the codes as bare literals (case 4: case 5: case 9:) with
    // nothing tying them to annotation_types -- so a renumbered annotation would
    // have reclassified beats here and not there.
    // ---------------------------------------------------------------------

    // ---------------------------------------------------------------------
    // Step 4a: Tukey -- CONSUMED, NOT RUN
    // ---------------------------------------------------------------------
    //
    // tukeyInPool() USED TO LIVE HERE AND IS GONE. It computed quartiles and
    // fences over the seed pool and rejected on RR length, which meant every
    // beat passed through a Tukey fence twice: once in alignment.hpp, which runs
    // four passes (RR length, amplitude, R location, wave score) over the sliced
    // beats, and again here over the survivors of those four. Two fences, two
    // populations, two answers, and no rule for which one a downstream reader
    // was looking at.
    //
    // It existed for a reason that was real: apply_mask() discards rejected
    // beats, so nothing downstream could ask what happened to a beat that did
    // not survive -- only survivors existed. Re-running the fence was the only
    // way to obtain a per-beat verdict. That gap is now closed at the source:
    // alignment records its verdict against each beat's ORIGINAL slice index in
    // a vector it never compacts, and reports each pass's quartiles and fences
    // alongside. There is one Tukey in this codebase, it is in alignment.hpp,
    // and this file reads what it decided.
    //
    // The two enums are mirrored rather than shared, because alignment.hpp is
    // upstream of template_bank.hpp and must not depend on it. Mirroring is the
    // failure annotation_types.hpp exists to prevent, so it is checked rather
    // than trusted -- exactly as annotation_code_check.hpp argued, but at
    // compile time.
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::NOT_ELIGIBLE)
        == static_cast<uint8_t>(tbank::TukeyOutcome::NOT_ELIGIBLE), "TukeyOutcome drift");
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::KEPT)
        == static_cast<uint8_t>(tbank::TukeyOutcome::KEPT), "TukeyOutcome drift");
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::REJ_RR_LENGTH)
        == static_cast<uint8_t>(tbank::TukeyOutcome::REJ_RR_LENGTH), "TukeyOutcome drift");
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::REJ_AMPLITUDE)
        == static_cast<uint8_t>(tbank::TukeyOutcome::REJ_AMPLITUDE), "TukeyOutcome drift");
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::REJ_R_LOCATION)
        == static_cast<uint8_t>(tbank::TukeyOutcome::REJ_R_LOCATION), "TukeyOutcome drift");
    static_assert(static_cast<uint8_t>(alignment::TukeyOutcome::REJ_WAVE_SCORE)
        == static_cast<uint8_t>(tbank::TukeyOutcome::REJ_WAVE_SCORE), "TukeyOutcome drift");

    // THE CONVERTERS BELOW HAVE NO CALLER and stay anyway. The static_asserts
    // are the point: they are the only thing checking that the two mirrored
    // enumerations still agree, and jbank writes tbank::TukeyOutcome values that
    // alignment's fences decided. Deleting the functions for tidiness would take
    // the check with them.
    inline tbank::TukeyOutcome fromAlignment(alignment::TukeyOutcome o) {
        return static_cast<tbank::TukeyOutcome>(static_cast<uint8_t>(o));
    }

    inline tbank::TukeyPassCounts fromAlignment(const alignment::TukeyStats& s) {
        tbank::TukeyPassCounts c;
        c.beats_in = s.beats_in;
        c.beats_out = s.beats_out;
        c.rejected = s.rejected;
        c.q1 = s.q1;
        c.q3 = s.q3;
        c.fence_lo = s.fence_lo;
        c.fence_hi = s.fence_hi;
        return c;
    }

    struct ChannelOutput {
        tbank::TemplateBank      bank;
        std::vector<int32_t>     assignment;   // per beat
        seed_pool::SeedSelection seed;
        std::vector<tbank::CapRaiseEvent> cap_raises;
        tbank::BinCounts         counts;
        std::vector<tbank::BeatFlags> flags;
        pvc_filter::FilterResult pvc;
    };

}  // namespace bin_pipeline
