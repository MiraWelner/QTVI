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
 *
 *         WHY STEP 4 IS THE WHOLE POINT. Tukey assumes one dominant
 *         population. Run over sinus + ectopic + noise it computes an IQR over
 *         a MIXTURE, which fattens it, which makes the fences too wide to catch
 *         anything. Restricted to category 1 the IQR describes sinus beats
 *         only, so the fences mean "unusual for a sinus beat".
 *
 *         alignment.hpp currently runs its three Tukey passes over all sliced
 *         beats and exempts flagged beats from rejection at line 175. Exemption
 *         is NOT the same thing as restriction: an exempted beat is still in
 *         the sample that computes the IQR, so it widens the fences it is then
 *         excused from. That is why the exemption becomes unnecessary once the
 *         passes move inside the pool -- ectopic beats were never in the sample
 *         to begin with.
 *
 *         AND WHY MOVING IT IS NOT A FIX, ONLY A RELOCATION. With ectopy in the
 *         minority both quartiles sit in the sinus cluster and every ectopic
 *         beat falls outside the fence. As ectopy approaches half the bin the
 *         IQR widens across both clusters and Tukey stops rejecting anything,
 *         including real artifact. Restricting the pool does not remove that
 *         failure; it makes Tukey's correctness CONDITIONAL ON CLASSIFICATION
 *         RECALL, one layer deeper and quieter, with no backstop behind it --
 *         Tukey WAS the backstop. The two mechanisms also degrade together
 *         rather than independently: at high burden the trailing-ten median is
 *         dragged down by ectopic RRs, so the prematurity test fires less often
 *         at exactly the burden where the fences have widened past usefulness.
 *
 *         Hence every pass reports its Q1, Q3 and fences alongside its
 *         rejection count. Counts give the outcome; the IQR gives the cause,
 *         and it is continuous where counts are discrete. The signature to
 *         watch is LOW REJECTION RATE WITH WIDE IQR, which is
 *         indistinguishable from "this bin is genuinely clean" if you only
 *         have the counts. This defect class produces output that looks fine.
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
    // ---------------------------------------------------------------------
    // runChannel() AND ChannelInput ARE GONE.
    //
    // runChannel(ChannelInput) -> ChannelOutput built a bank over ONE channel's
    // aligned beats: seed pool, slot-0 seeding, both assignment passes, the
    // merge, the cap and the census. Section 4.6 has exactly one partition,
    // shared by the three ECG leads and the pulse, and it is
    // jbank::buildBinBank -- driven per bin from make_averaged_templates.hpp.
    //
    // Its last two consumers are gone with it: the morphology archive now reads
    // the joint projection, and envelope_database.hpp was a stale copy of
    // create_ecg_templates.hpp that nothing included.
    //
    // IT WAS THE LAST PLACE THE ORDER RAN BACKWARDS. Slot 0 came from
    // seed_pool::selectSeedPool filtered by the Tukey verdict, so prematurity
    // and Tukey chose the seed the whole partition was scored against. The
    // required order is partition, then remove premature, then Tukey on what is
    // left, and there is now no code that does it the other way.
    //
    // ChannelOutput STAYS. It is the type the viewer, the serializer and the
    // morphology writers read, and jbank::projectToChannel fills it -- a channel
    // VIEW of the one partition rather than a partition of its own.
    // ---------------------------------------------------------------------

    struct ChannelOutput {
        tbank::TemplateBank      bank;
        std::vector<int32_t>     assignment;   // per beat
        seed_pool::SeedSelection seed;
        std::vector<tbank::CapRaiseEvent> cap_raises;
        tbank::BinCounts         counts;
        std::vector<tbank::BeatFlags> flags;
        pvc_filter::FilterResult pvc;
    };


    // ---------------------------------------------------------------------
    // WHAT USED TO BE HERE: runBin() AND ITS TYPES. Deleted as unreachable.
    //
    // runBin(BinInput) -> BinResult was a THREE-CHANNEL entry point that did the
    // same job as runChannel three times over, with its own copies of the seed
    // pool, the Tukey read, the slot-0 seeding, both bank passes and the census,
    // plus ChannelResult, BinResult, BinInput and RecordAccumulator. Nothing
    // outside this header ever called any of them -- the whole external surface
    // of this file is runChannel, ChannelInput and ChannelOutput -- so it was a
    // second implementation of the partition that no record was ever built by,
    // and it could not be kept correct by anything except reading it.
    //
    // It is also where tbank::polymorphicVerdict was computed, as a max over
    // three per-channel banks. That verdict has to be rebuilt on the joint
    // groups, where there is one partition to count, rather than restored here.
    // ---------------------------------------------------------------------

}  // namespace bin_pipeline
