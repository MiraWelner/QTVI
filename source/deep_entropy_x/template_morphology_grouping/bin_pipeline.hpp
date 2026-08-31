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
    // Step 4a: Tukey, restricted to a pool
    // ---------------------------------------------------------------------

    struct TukeyPassResult {
        tbank::TukeyPassCounts counts;
        std::vector<uint8_t>   rejected;   // parallel to `pool`
    };

    // Quartiles and fences computed over `pool` ONLY. Beats outside the pool
    // are neither tested nor counted -- they were never eligible, which is a
    // different state from "tested and kept" and is recorded as such in
    // BeatFlags::tukey.
    inline TukeyPassResult tukeyInPool(const std::vector<double>& values,
        const std::vector<uint32_t>& pool,
        double k = 1.5)
    {
        TukeyPassResult out;
        out.rejected.assign(pool.size(), 0u);
        out.counts.beats_in = static_cast<uint32_t>(pool.size());

        std::vector<double> v;
        v.reserve(pool.size());
        for (uint32_t i : pool)
            if (i < values.size() && !std::isnan(values[i])) v.push_back(values[i]);

        if (v.size() < 4) {          // too few order statistics to fence on
            out.counts.beats_out = out.counts.beats_in;
            return out;
        }
        std::sort(v.begin(), v.end());
        const size_t n = v.size();
        out.counts.q1 = v[n / 4];
        out.counts.q3 = v[std::min(n - 1, (3 * n) / 4)];
        const double iqr = out.counts.q3 - out.counts.q1;
        out.counts.fence_lo = out.counts.q1 - k * iqr;
        out.counts.fence_hi = out.counts.q3 + k * iqr;

        for (size_t j = 0; j < pool.size(); ++j) {
            const uint32_t i = pool[j];
            if (i >= values.size() || std::isnan(values[i])) continue;
            const double x = values[i];
            if (x < out.counts.fence_lo || x > out.counts.fence_hi) {
                out.rejected[j] = 1u;
                ++out.counts.rejected;
            }
        }
        out.counts.beats_out = out.counts.beats_in - out.counts.rejected;
        return out;
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

        // ---- step 3: PVC filter across ALL beats ------------------------
        out.pvc = pvc_filter::runFilter(in.rr_after);
        for (size_t i = 0; i < n && i < out.pvc.verdict.size(); ++i)
            out.flags[i].pvc = out.pvc.verdict[i];
        out.counts.n_premature = out.pvc.n_premature;
        out.counts.n_vote_only = out.pvc.n_vote_only;

        std::vector<uint8_t> rhythm(n, 0u);
        for (size_t i = 0; i < n; ++i)
            rhythm[i] = (out.flags[i].pvc == tbank::PvcFilter::PREMATURE) ? 1u
            : (out.flags[i].pvc == tbank::PvcFilter::VOTE) ? 2u : 0u;

        // ---- step 4: seed from category 1, Tukey inside that group ------
        std::vector<uint32_t> candidates;
        candidates.reserve(n);
        for (uint32_t i = 0; i < n; ++i) {
            if (i < in.baseline_ok.size() && in.baseline_ok[i] == 0u) continue;
            candidates.push_back(i);
        }

        out.seed = seed_pool::selectSeedPool(candidates, rhythm, cat);
        out.assignment.assign(n, -1);

        const TukeyPassResult t_rr =
            tukeyInPool(in.rr_after, out.seed.members, 1.5);
        out.counts.tukey_rr = t_rr.counts;

        std::vector<uint32_t> seeded;
        seeded.reserve(out.seed.members.size());
        for (size_t j = 0; j < out.seed.members.size(); ++j) {
            const uint32_t bi = out.seed.members[j];
            if (t_rr.rejected[j]) {
                out.flags[bi].tukey = tbank::TukeyOutcome::REJ_RR_LENGTH;
                if (out.flags[bi].category == tbank::Category::REGULAR)
                    ++out.counts.n_regular_rejected_on_rr;
                continue;
            }
            out.flags[bi].tukey = tbank::TukeyOutcome::KEPT;
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

            // Tukey INSIDE the seed pool. Only the RR-length pass is shown
            // here; the amplitude and R-location passes take the same shape
            // with their own value vectors, and each keeps its own fence state.
            const TukeyPassResult t_rr =
                tukeyInPool(in.rr_after, cr.seed.members, 1.5);
            res.counts.tukey_rr = t_rr.counts;

            std::vector<uint32_t> seeded;
            seeded.reserve(cr.seed.members.size());
            for (size_t j = 0; j < cr.seed.members.size(); ++j) {
                const uint32_t bi = cr.seed.members[j];
                if (t_rr.rejected[j]) {
                    res.flags[bi].tukey = tbank::TukeyOutcome::REJ_RR_LENGTH;
                    // A REGULAR beat rejected for being SHORT is the cheapest
                    // estimate of classifier recall available: very likely
                    // ectopy that classification missed. This count rising
                    // while total rejections fall is the high-burden failure
                    // announcing itself.
                    if (res.flags[bi].category == tbank::Category::REGULAR)
                        ++res.counts.n_regular_rejected_on_rr;
                    continue;
                }
                res.flags[bi].tukey = tbank::TukeyOutcome::KEPT;
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

}  // namespace bin_pipeline
