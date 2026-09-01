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

    // Annotation codes, read from annotation_types.hpp rather than restated.
    // The spec text says "PVC is annotation type 5"; the table says PVC is 4
    // and PAC is 5. The table wins -- it is the declared single source of truth
    // for marking codes, and taking the spec's number would label every PVC
    // template as PAC.
    //
    // Current placement, per operator decision:
    //   ECTOPIC : 4 PVC, 5 PAC, 9 VT
    //   NOISE   : 2 Minor Noise
    //   REGULAR : everything else, including 6 Cond. Delay, 7 AF, 8 SVT
    //
    // AF and SVT sit in REGULAR because those beats are morphologically normal
    // -- they would correlate ~0.97 against the sinus template anyway, so
    // calling them ectopic would mean fighting the metric. The consequence to
    // watch: they enter the Tukey IQR sample, and AF's defining feature is
    // irregular RR, so an AF span widens the RR-length fences for the whole
    // bin. The template stays clean while the INTERVAL features degrade, which
    // is a failure invisible in the template. If it shows up, the fix is a
    // timing-based exclusion from interval features only, not a category
    // change.
    inline tbank::Category categoryForCode(uint8_t code) {
        switch (code) {
        case 4: case 5: case 9: return tbank::Category::ECTOPIC;
        case 2:                 return tbank::Category::NOISE;
        default:                return tbank::Category::REGULAR;
        }
    }

    // Whether the beat FOLLOWING a marked beat inherits its code. This mirrors
    // annotation_types.hpp's postEligible flag, which is true for exactly
    // 4, 5, 7, 8, 9 -- the compensatory beat after an ectopic beat is not a
    // normal beat either and must not feed the sinus reference.
    inline bool postEligibleCode(uint8_t code) {
        return code == 4 || code == 5 || code == 7 || code == 8 || code == 9;
    }

    // mark_code[i] is the annotation code covering beat i, 0 where unmarked.
    // The absence of a mark IS the regular verdict, so this never leaves a beat
    // unclassified.
    inline std::vector<tbank::Category> categoryFromMarks(
        const std::vector<uint8_t>& mark_code)
    {
        const size_t n = mark_code.size();
        std::vector<tbank::Category> out(n, tbank::Category::REGULAR);
        for (size_t i = 0; i < n; ++i) {
            out[i] = categoryForCode(mark_code[i]);
            // Post-beat inheritance, applied forward so an unmarked beat
            // following an ectopic one is excluded from the sinus reference.
            // Does not overwrite a beat that carries its own mark.
            if (i > 0 && mark_code[i] == 0 && postEligibleCode(mark_code[i - 1])
                && categoryForCode(mark_code[i - 1]) == tbank::Category::ECTOPIC)
                out[i] = tbank::Category::ECTOPIC;
        }
        return out;
    }

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
    // Single-channel entry point
    //
    // build_ecg_template_for_method() in create_ecg_templates.hpp handles ONE
    // channel and returns one SingleMethodResult, so runBin()'s three-channel
    // shape does not fit it. Assignment is genuinely per channel and belongs
    // there; the polymorphic verdict is a max ACROSS channels and cannot be
    // computed until all three exist, so it is deferred to whoever assembles
    // them (make_averaged_templates / build_templates).
    //
    // NOTE ON MARKS. This runs in the template-generation pipeline, before any
    // operator review, so `mark_code` is normally empty and every beat reads
    // REGULAR. That is the design working rather than a gap: morphology does
    // the sorting, and operator marks only ever supply LABELS. The bank
    // separates the PVC morphology before anyone has called it a PVC. Marks
    // arrive later and propagate through tbank::propagateLabel().
    // ---------------------------------------------------------------------

    struct ChannelInput {
        const std::vector<std::vector<double>>* beats = nullptr;  // aligned, shared axis
        int    width = 0;
        int    r_col = -1;
        std::vector<double>  rr_after;      // ms, R[i+1]-R[i], parallel to beats
        std::vector<uint8_t> baseline_ok;   // 0 where baseline_source == NONE
        std::vector<uint8_t> mark_code;     // usually empty at build time
        bool     is_ppg = false;
        uint64_t bin_index = 0;
        int      channel = 0;

        // ---- TWO INDEX SPACES, AND THE MAP BETWEEN THEM -------------------
        //
        // `beats` is alignment's KEPT set, so everything in this struct that is
        // parallel to it -- rr_after, baseline_ok, assignment, flags -- lives in
        // ALIGNED space. tukey_outcome cannot: it has to describe beats that
        // were pruned, and a pruned beat has no aligned index. It is therefore
        // indexed by ORIGINAL SLICE index, and original_index is the map.
        //
        // Getting this wrong is silent. An aligned index used directly against
        // tukey_outcome reads some other beat's verdict, and the further into
        // the bin you go the further off it is -- every flag plausible, all of
        // them wrong. Hence the map rather than an assumption that the two
        // spaces coincide (they only do when nothing was pruned).
        //
        // Note what this makes visible: aligned.beats DOES contain rejected
        // beats, because apply_mask exempts rhythm-flagged ones. A beat that is
        // both premature and a Tukey RR outlier is retained, and this is the
        // only path by which the archive can say so -- pvc_filter.hpp calls that
        // the most interesting row in the file.
        std::vector<size_t> original_index;             // parallel to *beats
        std::vector<alignment::TukeyOutcome> tukey_outcome;   // by ORIGINAL index
        alignment::TukeyStats tukey_rr, tukey_amplitude,
            tukey_r_location, tukey_wave_score;

        // THE PHASE 1 SINUS TEMPLATE, which seeds slot 0. Section 4.6: "Seed
        // the bank with the sinus template from Phase 1."
        //
        // Supplied rather than re-derived. It used to be rebuilt here as the
        // median of the seed pool, which is a different waveform: features
        // measured in Phase 1 and features measured against slot 0 then used two
        // references and nothing in the archive said which. Seeding closes that.
        //
        // READ THIS BEFORE TRUSTING IT. Phase 1's template is a column-wise
        // median over `usable`, and `usable` is filtered on one condition --
        // baseline_source != NONE. There is NO rhythm test, so ectopic beats are
        // in it (see seed_pool.hpp, which was written about exactly this). It is
        // therefore the sinus template by name and not necessarily by content,
        // and seeding from it propagates that contamination into every spawn
        // decision in the bin, because every decision is scored against slot 0.
        // The fix belongs upstream, in create_ecg_templates.hpp: apply the seed
        // pool there and Phase 1's template becomes what this clause assumes it
        // already is. Until then, seed_fallback_used below is the thing to
        // watch, and phase1_is_verified_sinus records whether the caller can
        // vouch for it.
        std::vector<double> phase1_template;
        bool phase1_is_verified_sinus = false;

        // max_templates_per_bin, from cfg. 0 = use the built-in default.
        int32_t max_templates_per_bin = 0;
    };

    struct ChannelOutput {
        tbank::TemplateBank      bank;
        std::vector<int32_t>     assignment;   // per beat
        seed_pool::SeedSelection seed;
        std::vector<tbank::CapRaiseEvent> cap_raises;
        tbank::BinCounts         counts;
        std::vector<tbank::BeatFlags> flags;
        pvc_filter::FilterResult pvc;
    };

    inline ChannelOutput runChannel(const ChannelInput& in)
    {
        ChannelOutput out;
        if (!in.beats || in.beats->empty() || in.width <= 0) return out;
        const size_t n = in.beats->size();

        out.flags.assign(n, tbank::BeatFlags{});
        out.counts.beats_detected = static_cast<uint32_t>(n);

        // ---- step 2: classify -------------------------------------------
        std::vector<tbank::Category> cat;
        if (in.mark_code.size() == n) {
            cat = categoryFromMarks(in.mark_code);
        }
        else {
            // No marks yet. Every beat is REGULAR, which is correct: the
            // absence of a mark IS the regular verdict, and at build time no
            // mark has been made.
            cat.assign(n, tbank::Category::REGULAR);
        }
        for (size_t i = 0; i < n; ++i) {
            out.flags[i].category = cat[i];
            switch (cat[i]) {
            case tbank::Category::REGULAR: ++out.counts.n_regular; break;
            case tbank::Category::ECTOPIC: ++out.counts.n_ectopic; break;
            case tbank::Category::NOISE:   ++out.counts.n_noise;   break;
            }
        }

        // ---- step 3: PVC filter -- MEASURED NOW, APPLIED AFTER PARTITIONING
        //
        // ORDER: align -> partition and merge -> remove and flag premature ->
        // Tukey on the clean beats. The filter is computed here because it is
        // pure arithmetic on rr_after and cheap, but its verdicts must NOT
        // shape the partition: a premature beat is a beat, and excluding it
        // before the bank runs is what deleted the ectopic morphologies before
        // they could be separated. out.pvc holds the verdicts; they are written
        // onto the flags and acted on in the post-partition stage below.
        out.pvc = pvc_filter::runFilter(in.rr_after);
        for (size_t i = 0; i < n && i < out.pvc.verdict.size(); ++i)
            out.flags[i].pvc = out.pvc.verdict[i];
        out.counts.n_premature = out.pvc.n_premature;
        out.counts.n_vote_only = out.pvc.n_vote_only;

        // SEED POOL WITHOUT RHYTHM OR TUKEY. Both used to steer it, and both
        // are now decided after the partition exists. Seeding on category and
        // baseline alone is the only self-consistent choice: anything else
        // would reintroduce the inversion this reordering removes.
        std::vector<uint8_t> rhythm(n, 0u);

        // ---- step 4: seed from category 1, Tukey inside that group ------
        std::vector<uint32_t> candidates;
        candidates.reserve(n);
        for (uint32_t i = 0; i < n; ++i) {
            if (i < in.baseline_ok.size() && in.baseline_ok[i] == 0u) continue;
            candidates.push_back(i);
        }

        out.seed = seed_pool::selectSeedPool(candidates, rhythm, cat);
        out.assignment.assign(n, -1);

        // ---- step 4a: read alignment's Tukey verdict --------------------
        //
        // NOT RECOMPUTED. in.tukey_outcome is what alignment's four passes
        // decided about each beat, indexed by its original slice index and
        // never compacted, and in.tukey_* carry those passes' quartiles and
        // fences. Copied through, not re-derived, so the flags in the archive
        // and the beats in the templates were decided by the same fence.
        //
        // A beat with no verdict (an input that predates the recording of one)
        // reads NOT_ELIGIBLE, which is honest: no pass tested it here. It is
        // NOT silently promoted to KEPT -- that would assert a decision nobody
        // made, and is the direction that hides a missing plumbing step.
        out.counts.tukey_rr = fromAlignment(in.tukey_rr);
        out.counts.tukey_amplitude = fromAlignment(in.tukey_amplitude);
        out.counts.tukey_r_location = fromAlignment(in.tukey_r_location);
        out.counts.tukey_wave_score = fromAlignment(in.tukey_wave_score);

        std::vector<uint32_t> seeded;
        seeded.reserve(out.seed.members.size());
        for (size_t j = 0; j < out.seed.members.size(); ++j) {
            const uint32_t bi = out.seed.members[j];
            // Aligned index -> original slice index -> verdict. Never index
            // tukey_outcome with bi directly; see ChannelInput.
            const size_t orig = (bi < in.original_index.size())
                ? in.original_index[bi] : static_cast<size_t>(bi);
            const alignment::TukeyOutcome verdict =
                (orig < in.tukey_outcome.size()) ? in.tukey_outcome[orig]
                : alignment::TukeyOutcome::NOT_ELIGIBLE;
            out.flags[bi].tukey = fromAlignment(verdict);
            if (verdict != alignment::TukeyOutcome::KEPT) {
                if (verdict == alignment::TukeyOutcome::REJ_RR_LENGTH
                    && out.flags[bi].category == tbank::Category::REGULAR)
                    ++out.counts.n_regular_rejected_on_rr;
                continue;
            }
            seeded.push_back(bi);
        }

        // ---- slot 0: THE PHASE 1 SINUS TEMPLATE -------------------------
        // Seeded, not rebuilt. Section 4.6's first bullet. Slot 0 starts with
        // the supplied waveform and NO members, and every beat in the bin --
        // including the seed pool's own -- is then scored against it by the
        // second bullet ("for each beat ... assign it to the best match").
        // Nothing is pre-assigned, so no beat reaches a template without having
        // matched it.
        //
        // The seed pool and its Tukey pass are still computed above. They are
        // Section 4.5's requirement, they populate the counters, and they are
        // what SHOULD be feeding Phase 1's template upstream. They just no
        // longer hand slot 0 its membership.
        //
        // FALLBACK. With no Phase 1 template supplied there is nothing to seed
        // from and a bank cannot be built without a slot 0, so the old seed-pool
        // median is used and counted. A record where this fires is a record
        // whose slot 0 is not the Phase 1 template, which is exactly the
        // discrepancy the clause exists to remove -- so it must be visible in
        // the archive rather than inferred.
        if (in.max_templates_per_bin > 0)
            out.bank.setCap(in.max_templates_per_bin);

        tbank::BankTemplate slot0;
        slot0.spawn_seq = out.bank.next_spawn_seq++;
        slot0.r_col = in.r_col;

        if (!in.phase1_template.empty()) {
            slot0.tmpl = in.phase1_template;
            slot0.tmpl.resize(in.width,
                std::numeric_limits<double>::quiet_NaN());
            // A seeded slot 0 has no members yet, so it has no corridor of its
            // own. Derive one from the seed pool's spread -- the same beats
            // Phase 1 should have built the template from -- so the first beat
            // scored has something to be scored against. Once beats accumulate,
            // recomputeAll() replaces it with slot 0's measured corridor.
            tbank::BankTemplate spread;
            spread.members = seeded;
            tbank::recomputeTemplate(spread, *in.beats, in.width);
            slot0.band_lo.assign(in.width,
                std::numeric_limits<double>::quiet_NaN());
            slot0.band_hi.assign(in.width,
                std::numeric_limits<double>::quiet_NaN());
            for (int c = 0; c < in.width; ++c) {
                if (std::isnan(slot0.tmpl[c])) continue;
                double half = std::numeric_limits<double>::quiet_NaN();
                if (c < (int)spread.band_lo.size()
                    && !std::isnan(spread.band_lo[c])
                    && !std::isnan(spread.band_hi[c]))
                    half = 0.5 * (spread.band_hi[c] - spread.band_lo[c]);
                if (std::isnan(half) || half <= 0.0) continue;
                slot0.band_lo[c] = slot0.tmpl[c] - half;
                slot0.band_hi[c] = slot0.tmpl[c] + half;
            }
            out.counts.seed_from_phase1 = 1u;
            out.counts.seed_verified_sinus = in.phase1_is_verified_sinus ? 1u : 0u;
        }
        else {
            slot0.members = seeded;
            tbank::recomputeTemplate(slot0, *in.beats, in.width);
            out.counts.seed_fallback_used = 1u;
            for (uint32_t bi : seeded) out.assignment[bi] = 0;
        }

        out.bank.templates.push_back(std::move(slot0));
        out.bank.assigned_beats =
            static_cast<uint32_t>(out.bank.templates[0].members.size());

        // ---- bank pass 1, then pass 2 -----------------------------------
        // Pass 1 defers template recomputation (see assignBeat); everything is
        // rebuilt once here, before pass 2 scores against it.
        for (uint32_t bi = 0; bi < n; ++bi) {
            if (out.assignment[bi] >= 0) continue;
            const tbank::AssignOutcome ao = tbank::assignBeat(
                out.bank, bi, *in.beats, in.width, in.is_ppg,
                in.bin_index, in.channel, &out.cap_raises, &out.counts);
            out.assignment[bi] = ao.template_id;
        }
        tbank::recomputeAll(out.bank, *in.beats, in.width);
        tbank::refinePass(out.bank, *in.beats, out.assignment,
            in.width, in.is_ppg, &out.counts);

        // ==================================================================
        // POST-PARTITION: REMOVE AND FLAG PREMATURE, THEN TUKEY PER GROUP
        //
        // This is the stage the pipeline was missing. The order is align ->
        // partition and merge -> remove and flag premature beats -> Tukey on
        // the clean beats, and both of the last two steps happen HERE, after
        // the bank exists, because both are judgements about a beat relative to
        // ITS OWN morphology.
        //
        // Tukey used to run inside align_beat_matrix, before any partition. Its
        // amplitude and wave-score passes reject beats that do not resemble the
        // population, and before partitioning the population is every
        // morphology at once -- so an ectopic beat is an outlier by
        // construction and was deleted before the bank could separate it.
        // alignment.hpp now measures and records without pruning
        // (kTukeyPrunesInAlignment == false) and this stage does the excluding.
        //
        // EXCLUDED FROM THE TEMPLATE, RETAINED WITH FLAGS. A flagged beat keeps
        // its membership record in out.assignment and its verdict on
        // out.flags -- it is still that group's beat and the archive must be
        // able to say so -- but it is dropped from the member list the waveform
        // is averaged from. Deleting it instead would make the count
        // unreconstructable; leaving it in the average is the contamination the
        // section exists to prevent.
        {
            const double NaN = std::numeric_limits<double>::quiet_NaN();

            for (int t = 0; t < out.bank.size(); ++t) {
                tbank::BankTemplate& tp = out.bank.templates[t];

                // ---- 1: premature removal, per group --------------------
                std::vector<uint32_t> clean;
                clean.reserve(tp.members.size());
                for (const uint32_t m : tp.members) {
                    const bool premature = (m < out.flags.size())
                        && (out.flags[m].pvc != tbank::PvcFilter::NONE);
                    if (!premature) clean.push_back(m);
                }

                // A GROUP THAT IS ENTIRELY PREMATURE IS NOT CLEANED, IT IS
                // KEPT AS IT IS. That is the ectopic morphology itself, and
                // emptying it would delete exactly what 4.6 exists to
                // preserve -- the template would lose its waveform and the
                // operator would have nothing to confirm.
                if (clean.empty()) clean = tp.members;

                // ---- 2: Tukey INSIDE the group, THREE PASSES ------------
                //
                // R-LOCATION, AMPLITUDE, WAVE-SCORE. All three measured on the
                // SAME post-premature set and rejected as a union, NOT
                // sequentially. alignment.hpp ran them in series, each pass
                // computing its quartiles over the previous pass's survivors,
                // which is far more aggressive: three sequential 1.5*IQR passes
                // on a 12-member group can cascade to almost nothing, and small
                // groups are exactly what this pipeline produces. Simultaneous
                // fences are also order-independent, so nobody has to know
                // which pass ran first to interpret the result.
                //
                // RR-LENGTH IS NOT ONE OF THEM. It is a rhythm test, not a
                // morphology test, and premature removal above has already
                // taken the short intervals. What would be left for it to
                // reject is mostly long ones -- post-ectopic pauses and dropped
                // detections -- and inside a group of ectopics it would compute
                // ectopic RR fences and reject on those, which is a different
                // statement again. It stays a bin-level rhythm flag.
                //
                // Wave-score is the pass that gains most from the reorder: it
                // now scores a beat against ITS OWN group's template instead of
                // against a template built from every morphology at once, which
                // is what the test was always supposed to mean.
                std::vector<double> ampV(clean.size(), NaN);
                std::vector<double> rLoc(clean.size(), NaN);
                std::vector<double> wave(clean.size(), NaN);

                for (size_t k = 0; k < clean.size(); ++k) {
                    const uint32_t m = clean[k];
                    if (m >= in.beats->size()) continue;
                    const std::vector<double>& b = (*in.beats)[m];

                    double lo = std::numeric_limits<double>::infinity();
                    double hi = -std::numeric_limits<double>::infinity();
                    int argmax = -1; double peak = -1.0;
                    for (size_t j = 0; j < b.size(); ++j) {
                        const double v = b[j];
                        if (std::isnan(v)) continue;
                        lo = std::min(lo, v);
                        hi = std::max(hi, v);
                        const double a = std::abs(v);
                        if (a > peak) { peak = a; argmax = static_cast<int>(j); }
                    }
                    if (hi >= lo) ampV[k] = hi - lo;

                    // Offset of this beat's own dominant deflection from the
                    // frame's R column. A beat whose R did not land where the
                    // frame says it did corrupts every shape measurement made
                    // on it, which is why this is measured rather than assumed.
                    if (argmax >= 0 && tp.r_col >= 0)
                        rLoc[k] = static_cast<double>(argmax - tp.r_col);

                    // Shape against this group's own template. NaN when the
                    // template has no waveform yet; correlate() returns NaN r
                    // for a flat vector, which keep_within_tukey treats as
                    // unjudgeable rather than as an outlier.
                    if (!tp.tmpl.empty())
                        wave[k] = tbank::correlate(b, tp.tmpl).r;
                }

                std::vector<uint32_t> kept = clean;
                if (clean.size() >= 8) {
                    alignment::TukeyStats sA, sR, sW;
                    const std::vector<bool> keepA =
                        alignment::keep_within_tukey(ampV, 1.5, &sA);
                    const std::vector<bool> keepR =
                        alignment::keep_within_tukey(rLoc, 1.5, &sR);
                    const std::vector<bool> keepW =
                        alignment::keep_within_tukey(wave, 1.5, &sW);

                    std::vector<uint32_t> k2;
                    k2.reserve(clean.size());
                    for (size_t k = 0; k < clean.size(); ++k) {
                        // A metric that could not be measured does NOT reject.
                        // NaN means unjudgeable, and rejecting on it would prune
                        // beats for being unmeasurable rather than for being
                        // outliers -- the same conflation bandMatch avoids by
                        // excluding incomparable columns from its denominator.
                        const bool okA = std::isnan(ampV[k])
                            || (k < keepA.size() && keepA[k]);
                        const bool okR = std::isnan(rLoc[k])
                            || (k < keepR.size() && keepR[k]);
                        const bool okW = std::isnan(wave[k])
                            || (k < keepW.size() && keepW[k]);
                        if (okA && okR && okW) { k2.push_back(clean[k]); continue; }

                        const uint32_t m = clean[k];
                        if (m < out.flags.size()) {
                            // FIRST failing pass wins the label, so the reason
                            // is a single value. Ordered R-location, amplitude,
                            // wave-score: an R that landed wrong explains a bad
                            // amplitude and a bad shape, so it is the more
                            // informative attribution when several fire.
                            // fromAlignment(), not a direct assignment:
                            // BeatFlags::tukey is tbank::TukeyOutcome while the
                            // reason names here are alignment::TukeyOutcome. The
                            // two enumerations share values today, which is
                            // exactly why the conversion is explicit -- a
                            // static_cast would keep compiling if either list
                            // ever gained a member.
                            out.flags[m].tukey = fromAlignment(
                                !okR ? alignment::TukeyOutcome::REJ_R_LOCATION
                                : (!okA ? alignment::TukeyOutcome::REJ_AMPLITUDE
                                    : alignment::TukeyOutcome::REJ_WAVE_SCORE));
                        }
                    }
                    if (!k2.empty()) kept = std::move(k2);   // never empty a group

                    // out.counts.tukey_* are tbank::TukeyPassCounts, not
                    // alignment::TukeyStats. Accumulating COUNTS across groups
                    // is meaningful; accumulating FENCES is not -- q1/q3 and the
                    // fences describe one group's distribution, and summing or
                    // overwriting them across groups yields a quartile that
                    // describes nothing. Only the three counts are added; the
                    // fence fields stay NaN at bin level.
                    auto acc = [&](tbank::TukeyPassCounts& dst,
                        const alignment::TukeyStats& src) {
                            dst.beats_in += src.beats_in;
                            dst.beats_out += src.beats_out;
                            dst.rejected += src.rejected;
                        };
                    acc(out.counts.tukey_amplitude, sA);
                    acc(out.counts.tukey_r_location, sR);
                    acc(out.counts.tukey_wave_score, sW);
                }

                // ---- 3: record the clean subset, KEEP the full membership
                // members is NOT overwritten. An excluded beat stays this
                // template's beat so the archive can write it with its reason;
                // members_clean is what the waveform and the screen use.
                tp.members_clean = std::move(kept);
            }

            // ---- rebuild the waveforms from the CLEAN members ------------
            // recomputeAll averages over `members`, so the clean lists are
            // swapped in for the duration and swapped back after. Swapped
            // rather than passed as a parameter because recomputeTemplate is
            // shared with the pre-partition paths and changing its signature
            // would silently alter what those average over.
            //
            // Corridors depend on membership, so this happens once after every
            // group is cleaned rather than per group -- slot 0's corridor is
            // what the young groups inherit and must be final first.
            std::vector<std::vector<uint32_t>> full(out.bank.size());
            for (int t = 0; t < out.bank.size(); ++t) {
                full[t] = out.bank.templates[t].members;
                if (!out.bank.templates[t].members_clean.empty())
                    out.bank.templates[t].members = out.bank.templates[t].members_clean;
            }
            tbank::recomputeAll(out.bank, *in.beats, in.width);

            // ---- SEED-QUALITY CORRECTION: reassign against clean templates --
            //
            // The seed pool can no longer be steered by rhythm or Tukey
            // verdicts, because under this order neither exists when it runs --
            // both are decided after the partition. Slot 0's initial corridor
            // is therefore derived from a less filtered population than before,
            // and the first pass assigned beats against it.
            //
            // Fixed by reassigning now, against templates rebuilt from the
            // CLEAN members. This is the first moment a group's waveform
            // reflects only the beats that belong in it, so it is the first
            // moment assignment can be trusted -- doing it earlier would just
            // repeat the first pass against the same unfiltered reference.
            //
            // Runs with the clean lists still swapped in, deliberately: the
            // corridors a beat is scored against must be the clean ones. New
            // members land in `members`, which at this instant IS the clean
            // list, and the union with the full list is taken below so no beat
            // is lost from the archive.
            tbank::refinePass(out.bank, *in.beats, out.assignment,
                in.width, in.is_ppg, &out.counts);
            tbank::recomputeAll(out.bank, *in.beats, in.width);

            for (int t = 0; t < out.bank.size(); ++t) {
                tbank::BankTemplate& tp = out.bank.templates[t];
                // What refinePass produced is the clean membership.
                tp.members_clean = tp.members;
                // The archive keeps everything: the pre-clean membership plus
                // anything reassignment moved in. Union, not replacement -- a
                // beat excluded from this group's average is still a beat this
                // group is responsible for reporting.
                std::vector<uint32_t>& f = full[t];
                f.insert(f.end(), tp.members_clean.begin(), tp.members_clean.end());
                std::sort(f.begin(), f.end());
                f.erase(std::unique(f.begin(), f.end()), f.end());
                tp.members = std::move(f);
            }
        }


        // ---- per-template census ----------------------------------------
        // Aggregate the per-beat verdicts onto the templates their beats landed
        // in. This is what makes a template's PRESUMED category computable
        // before any operator mark exists, and it is also the "premature" row of
        // templates.csv -- which is deliberately NOT the same thing as the class
        // row: a template named PVC_A and a template whose members are mostly
        // premature will usually agree, but 4.6 never reassigns a category, it
        // only excludes, so the two must be reported side by side rather than
        // reconciled.
        for (auto& t : out.bank.templates) {
            t.n_premature_members = 0;
            t.n_voted_members = 0;
            t.n_noise_members = 0;
            for (uint32_t m : t.members) {
                if (m >= out.flags.size()) continue;
                switch (out.flags[m].pvc) {
                case tbank::PvcFilter::PREMATURE: ++t.n_premature_members; break;
                case tbank::PvcFilter::VOTE:      ++t.n_voted_members;     break;
                default: break;
                }
                if (out.flags[m].category == tbank::Category::NOISE)
                    ++t.n_noise_members;
            }
        }

        for (uint32_t bi = 0; bi < n; ++bi) {
            if (in.channel >= 0 && in.channel < 3)
                out.flags[bi].template_id_ecg[in.channel] = out.assignment[bi];
            if (in.is_ppg) out.flags[bi].template_id_ppg = out.assignment[bi];
        }

        // Every spawned template inherits slot 0's R column: the bank works on
        // the shared axis, where alignment placed every beat's detected R at
        // r_aligned_col regardless of morphology.
        for (auto& t : out.bank.templates) if (t.r_col < 0) t.r_col = in.r_col;

        return out;
    }

    // ---------------------------------------------------------------------
    // Per-bin, per-channel result
    // ---------------------------------------------------------------------

    struct ChannelResult {
        tbank::TemplateBank        bank;
        std::vector<int32_t>       assignment;   // per beat, -1 / kUnscorable
        seed_pool::SeedSelection   seed;
        std::vector<tbank::CapRaiseEvent> cap_raises;
    };

    struct BinResult {
        std::vector<tbank::BeatFlags> flags;      // per beat, all three axes
        pvc_filter::FilterResult      pvc;
        tbank::BinCounts              counts;
        std::array<ChannelResult, 3>  ecg;
        tbank::PolymorphicVerdict     verdict;
    };

    struct BinInput {
        // Per beat, all indexed identically and none of them compacted.
        std::vector<uint8_t> mark_code;    // 0 = unmarked
        std::vector<double>  rr_after;     // R[i+1] - R[i], ms
        std::vector<uint8_t> baseline_ok;  // 0 where baseline_source == NONE

        // Aligned beat waveforms per ECG channel, on the shared axis.
        std::array<std::vector<std::vector<double>>, 3> beats;
        std::array<int, 3> width = { 0, 0, 0 };
        std::array<bool, 3> channel_present = { false, false, false };

        // Phase 1 sinus template per ECG channel, on the same shared axis as
        // beats[c]. Empty = not supplied; the pipeline falls back to the seed
        // pool median and counts it. See ChannelInput::phase1_template for the
        // caveat that matters.
        std::array<std::vector<double>, 3> phase1_template;
        std::array<bool, 3> phase1_is_verified_sinus = { false, false, false };

        // max_templates_per_bin, from cfg. 0 = use the built-in default.
        int32_t max_templates_per_bin = 0;

        // ---- THE SINGLE TUKEY'S VERDICT, from alignment.hpp ---------------
        // Per beat, indexed by ORIGINAL slice index, uncompacted -- copy
        // ecg_beat_set::tukey_outcome straight in. Empty means no verdict was
        // supplied, and every beat then reads NOT_ELIGIBLE rather than being
        // promoted to KEPT: asserting a decision nobody made is the direction
        // that hides a missing plumbing step.
        std::vector<alignment::TukeyOutcome> tukey_outcome;

        // Map from aligned index to original slice index, per channel, parallel
        // to beats[c]. See ChannelInput for why tukey_outcome needs it.
        std::array<std::vector<size_t>, 3> original_index;

        // The four passes' quartiles, fences and counts, copied from the same
        // ecg_beat_set. Nothing here recomputes them.
        alignment::TukeyStats tukey_rr, tukey_amplitude,
            tukey_r_location, tukey_wave_score;

        // R column on the shared axis, per channel. Every beat's detected R
        // sits here by construction (alignment put it there), so it is a
        // property of the bin rather than of any beat. ChannelInput already
        // carried this; BinInput did not, which is one reason the morphology
        // writers had no way to emit their r_col descriptor row.
        std::array<int, 3> r_col = { -1, -1, -1 };

        uint64_t bin_index = 0;
    };

    // ---------------------------------------------------------------------
    // The driver
    // ---------------------------------------------------------------------

    inline BinResult runBin(const BinInput& in)
    {
        BinResult res;
        const size_t n = in.mark_code.size();
        res.flags.assign(n, tbank::BeatFlags{});
        res.counts.beats_detected = static_cast<uint32_t>(n);

        // ---- step 2: classify. Remove nothing. --------------------------
        const std::vector<tbank::Category> cat = categoryFromMarks(in.mark_code);
        for (size_t i = 0; i < n; ++i) {
            res.flags[i].category = cat[i];
            switch (cat[i]) {
            case tbank::Category::REGULAR: ++res.counts.n_regular; break;
            case tbank::Category::ECTOPIC: ++res.counts.n_ectopic; break;
            case tbank::Category::NOISE:   ++res.counts.n_noise;   break;
            }
        }

        // ---- step 3: PVC filter across ALL beats ------------------------
        // Not restricted to any category. Running it over everything is what
        // makes the timing-vs-shape agreement matrix fall out of the archive
        // for free: an ECTOPIC beat with no prematurity flag is a direct
        // measurement of the timing test's blind spot.
        res.pvc = pvc_filter::runFilter(in.rr_after);
        for (size_t i = 0; i < n && i < res.pvc.verdict.size(); ++i)
            res.flags[i].pvc = res.pvc.verdict[i];
        res.counts.n_premature = res.pvc.n_premature;
        res.counts.n_vote_only = res.pvc.n_vote_only;

        // Rhythm in create_ecg_templates' kept_rhythm encoding, for the seed
        // pool: 0 NORMAL, 1 PVC (premature), 2 VOTED_PVC.
        std::vector<uint8_t> rhythm(n, 0u);
        for (size_t i = 0; i < n; ++i)
            rhythm[i] = (res.flags[i].pvc == tbank::PvcFilter::PREMATURE) ? 1u
            : (res.flags[i].pvc == tbank::PvcFilter::VOTE) ? 2u : 0u;

        // ---- step 4: build from category 1, Tukey inside that group -----
        std::vector<uint32_t> candidates;
        candidates.reserve(n);
        for (uint32_t i = 0; i < n; ++i) {
            if (i < in.baseline_ok.size() && in.baseline_ok[i] == 0u) continue;
            candidates.push_back(i);
        }

        for (int c = 0; c < 3; ++c) {
            if (!in.channel_present[c]) continue;
            ChannelResult& cr = res.ecg[c];
            cr.assignment.assign(n, -1);

            // The seed pool: category 1 with the rhythm flags applied. This is
            // the one place a rhythm flag should hard-gate anything, because
            // every spawn decision in the bank is scored against slot 0 and a
            // seed sitting between sinus and PVC lets PVCs clear 0.85 against
            // it -- after which no ectopic template ever opens and the whole
            // section is silently inert.
            cr.seed = seed_pool::selectSeedPool(candidates, rhythm, cat);

            // ---- step 4a: read alignment's Tukey verdict ----------------
            //
            // NOT RECOMPUTED. All four passes ran once, in alignment.hpp, over
            // the sliced beats; in.tukey_outcome is what they decided, indexed
            // by original slice index and never compacted. Running a fence again
            // here -- which is what tukeyInPool did -- put every beat through
            // two fences over two different populations, so the flags in the
            // archive and the beats in the templates were decided by different
            // quartiles.
            //
            // The counts are copied from the passes that made the decisions, so
            // the low-rejection-with-wide-IQR signature this file's header is
            // about now describes a pass that actually pruned something.
            res.counts.tukey_rr = fromAlignment(in.tukey_rr);
            res.counts.tukey_amplitude = fromAlignment(in.tukey_amplitude);
            res.counts.tukey_r_location = fromAlignment(in.tukey_r_location);
            res.counts.tukey_wave_score = fromAlignment(in.tukey_wave_score);

            std::vector<uint32_t> seeded;
            seeded.reserve(cr.seed.members.size());
            for (size_t j = 0; j < cr.seed.members.size(); ++j) {
                const uint32_t bi = cr.seed.members[j];
                const size_t orig = (bi < in.original_index[c].size())
                    ? in.original_index[c][bi] : static_cast<size_t>(bi);
                const alignment::TukeyOutcome verdict =
                    (orig < in.tukey_outcome.size()) ? in.tukey_outcome[orig]
                    : alignment::TukeyOutcome::NOT_ELIGIBLE;
                res.flags[bi].tukey = fromAlignment(verdict);
                if (verdict != alignment::TukeyOutcome::KEPT) {
                    // A REGULAR beat rejected for being SHORT is the cheapest
                    // estimate of classifier recall available: very likely
                    // ectopy that classification missed. This count rising
                    // while total rejections fall is the high-burden failure
                    // announcing itself.
                    if (verdict == alignment::TukeyOutcome::REJ_RR_LENGTH
                        && res.flags[bi].category == tbank::Category::REGULAR)
                        ++res.counts.n_regular_rejected_on_rr;
                    continue;
                }
                seeded.push_back(bi);
            }

            // ---- slot 0: THE PHASE 1 SINUS TEMPLATE ---------------------
            // See the single-channel path above for why this is seeded rather
            // than rebuilt, and for the caveat about what Phase 1's template
            // actually contains.
            if (in.max_templates_per_bin > 0)
                cr.bank.setCap(in.max_templates_per_bin);

            tbank::BankTemplate slot0;
            slot0.spawn_seq = cr.bank.next_spawn_seq++;

            if (!in.phase1_template[c].empty()) {
                slot0.tmpl = in.phase1_template[c];
                slot0.tmpl.resize(in.width[c],
                    std::numeric_limits<double>::quiet_NaN());
                tbank::BankTemplate spread;
                spread.members = seeded;
                tbank::recomputeTemplate(spread, in.beats[c], in.width[c]);
                slot0.band_lo.assign(in.width[c],
                    std::numeric_limits<double>::quiet_NaN());
                slot0.band_hi.assign(in.width[c],
                    std::numeric_limits<double>::quiet_NaN());
                for (int k = 0; k < in.width[c]; ++k) {
                    if (std::isnan(slot0.tmpl[k])) continue;
                    double half = std::numeric_limits<double>::quiet_NaN();
                    if (k < (int)spread.band_lo.size()
                        && !std::isnan(spread.band_lo[k])
                        && !std::isnan(spread.band_hi[k]))
                        half = 0.5 * (spread.band_hi[k] - spread.band_lo[k]);
                    if (std::isnan(half) || half <= 0.0) continue;
                    slot0.band_lo[k] = slot0.tmpl[k] - half;
                    slot0.band_hi[k] = slot0.tmpl[k] + half;
                }
                res.counts.seed_from_phase1 = 1u;
                res.counts.seed_verified_sinus =
                    in.phase1_is_verified_sinus[c] ? 1u : 0u;
            }
            else {
                slot0.members = seeded;
                tbank::recomputeTemplate(slot0, in.beats[c], in.width[c]);
                res.counts.seed_fallback_used = 1u;
                for (uint32_t bi : seeded) cr.assignment[bi] = 0;
            }

            cr.bank.templates.push_back(std::move(slot0));
            cr.bank.assigned_beats =
                static_cast<uint32_t>(cr.bank.templates[0].members.size());

            // ---- bank pass 1: causal assignment over every OTHER beat ---
            // Every beat, not just the ectopic ones. Noise beats are routed
            // here too and separate themselves by failing to match anything:
            // ectopy is reproducible morphology and accumulates members, noise
            // is non-reproducible by definition and sits alone. Member count is
            // therefore the discriminator, and no classifier is needed to
            // produce it.
            for (uint32_t bi = 0; bi < n; ++bi) {
                if (cr.assignment[bi] >= 0) continue;              // in slot 0
                if (bi >= in.beats[c].size()) continue;
                const tbank::AssignOutcome ao = tbank::assignBeat(
                    cr.bank, bi, in.beats[c], in.width[c],
                    /*is_ppg=*/false, in.bin_index, c,
                    &cr.cap_raises, &res.counts);
                cr.assignment[bi] = ao.template_id;
                res.flags[bi].template_id_ecg[c] = ao.template_id;
            }

            // ---- bank pass 2: refinement --------------------------------
            tbank::recomputeAll(cr.bank, in.beats[c], in.width[c]);
            tbank::refinePass(cr.bank, in.beats[c], cr.assignment,
                in.width[c], /*is_ppg=*/false, &res.counts);
            for (uint32_t bi = 0; bi < n; ++bi)
                res.flags[bi].template_id_ecg[c] = cr.assignment[bi];
        }

        // ---- bin verdict ------------------------------------------------
        // Max over ECG channels of CONFIRMED PVC templates. Because only
        // confirmed templates count, this reads 1 until the operator confirms a
        // beat in a second template: polymorphy is an operator-gated finding,
        // never an algorithmic assertion. Do not report the count as an
        // algorithmic output.
        res.verdict = tbank::polymorphicVerdict(
            res.ecg[0].bank, res.ecg[1].bank, res.ecg[2].bank);

        return res;
    }


    // ---------------------------------------------------------------------
    // Record accumulator
    // ---------------------------------------------------------------------
    //
    // WHY THIS EXISTS. morphology_csv's ChannelBlock wants, per channel, a
    // vector of ChannelOutput over the whole record plus the beat matrices and
    // R columns. runChannel() produces a ChannelOutput; runBin() does not -- it
    // produces ChannelResult per channel (bank, assignment, seed, cap_raises)
    // and keeps flags, pvc and counts at the BIN level, because the
    // classification and the PVC filter run once over all beats and not three
    // times. So the two sides never fit together and nothing could call the
    // writers. That is the whole reason no _beats.csv or _templates.csv has ever
    // been produced: the writers are complete and were unreachable.
    //
    // WHY THE BIN-LEVEL FIELDS ARE COPIED INTO ALL THREE CHANNELS rather than
    // split. Category, prematurity and the Tukey verdict are per BEAT, and a
    // beat is one event recorded on three leads -- it is premature or it is not,
    // regardless of which lead you look at. Copying is therefore not
    // duplication of unrelated data, it is the same fact stated in each
    // channel's view, and the writer needs it in that view because it emits one
    // column block per channel. What genuinely differs per channel -- bank,
    // assignment, cap raises, seed -- comes from ChannelResult.
    //
    // Beat matrices are COPIED, not referenced. The caller's BinInput is
    // typically a loop variable that is refilled for the next bin, so holding
    // pointers into it would leave the writers reading the last bin's samples
    // for every column. That is a real cost -- a full record of aligned beats
    // in memory at once -- and it is the reason to call writeAll() per record
    // and then clear(), not to accumulate across an entire study.
    struct RecordAccumulator {
        std::array<std::vector<ChannelOutput>, 3> per_bin;
        std::array<std::vector<std::vector<std::vector<double>>>, 3> beats;
        std::array<std::vector<int>, 3> r_col;

        void addBin(const BinInput& in, const BinResult& res) {
            for (int c = 0; c < 3; ++c) {
                if (!in.channel_present[c]) continue;

                ChannelOutput out;
                out.bank = res.ecg[c].bank;
                out.assignment = res.ecg[c].assignment;
                out.seed = res.ecg[c].seed;
                out.cap_raises = res.ecg[c].cap_raises;
                out.flags = res.flags;     // per beat, shared across leads
                out.pvc = res.pvc;
                out.counts = res.counts;

                per_bin[c].push_back(std::move(out));
                beats[c].push_back(in.beats[c]);
                r_col[c].push_back(in.r_col[c]);
            }
        }

        // A channel with no bins is omitted from the output rather than emitted
        // empty, so a two-lead record does not produce a CH3 block of nothing.
        bool has(int c) const { return c >= 0 && c < 3 && !per_bin[c].empty(); }

        void clear() {
            for (int c = 0; c < 3; ++c) {
                per_bin[c].clear(); beats[c].clear(); r_col[c].clear();
            }
        }
    };

}  // namespace bin_pipeline
