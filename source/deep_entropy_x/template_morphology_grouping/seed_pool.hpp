#pragma once
/**
 * @file   seed_pool.hpp
 * @brief  The ectopic mask that was specified but never built: selects the
 *         beats allowed to form bank slot 0, so that "seed the bank with the
 *         sinus template" is true rather than aspirational.
 *
 *         WHAT IS ACTUALLY IN THE TREE TODAY. alignment.hpp assigns the
 *         prematurity and 5-of-8 vote flags before pruning and exempts flagged
 *         beats from every apply_mask (line 175), so ectopic beats survive to
 *         the template stage -- and the comment at line 171 says excluding
 *         them is "the ectopic mask's job (create_ecg_templates.hpp)". That
 *         mask does not exist. create_ecg_templates.hpp:155 builds the median
 *         over `usable`, and `usable` is filtered on exactly one condition:
 *         baseline_source == NONE. No rhythm test. The flags then travel from
 *         kept_rhythm through make_averaged_templates into
 *         beats.per_channel_rhythm and are read by nothing.
 *
 *         So the current net effect of the flags is to INCREASE ectopic
 *         contamination: they rescue beats from RR-length pruning and then
 *         nothing excludes them from the median.
 *
 *         WHY SLOT 0 IS THE ONE PLACE A FLAG SHOULD HARD-GATE ANYTHING. Every
 *         spawn decision in the bank is scored against slot 0. A seed whose
 *         shape sits between sinus and PVC lets PVCs correlate above 0.85
 *         against it, no ectopic template ever opens, and Section 4.6 becomes
 *         silently inert. Elsewhere the flags should only inform, not route:
 *         routing by flag would build a timing-based classifier with a
 *         shape-based storage layer and inherit every prematurity blind spot
 *         (PACs look like sinus, late ventricular beats are not premature,
 *         mid-run beats sit where the collapsed trailing median puts them).
 *
 *         ASYMMETRIC COSTS ARE WHY THIS GATE IS STRICT. A misclassified beat
 *         is one wrong number in the burden statistics. A contaminated seed
 *         damages the template, the corridor built from the same pool, and
 *         every feature for the whole bin -- and it compounds, because a wider
 *         corridor admits the next ectopic beat more easily. False rejection
 *         costs one beat's contribution to a median over hundreds. Reject on
 *         any doubt.
 *
 *         SCOPE. This selects the pool. Relocating Tukey to run INSIDE the
 *         selected pool is a separate change in alignment.hpp, where the three
 *         Tukey passes currently run over all sliced beats -- see the note on
 *         tukey_relocation_pending below.
 */

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "template_bank.hpp"   // tbank::Category

namespace seed_pool {

    // Below this many sinus beats the median rests on too few contributors to
    // be a usable reference, and the fallback ladder engages rather than
    // producing a seed nobody should trust. Deliberately lower than
    // MorphologyEnvelope::kTightMinBeats (40): a median is far more tolerant
    // of small n than a 2.5/97.5 corridor is.
    inline constexpr int kMinSeedBeats = 12;

    // How the pool was chosen. Recorded per bin per channel, because a bin
    // whose seed fell back is a bin whose features are measured against a
    // reference that is not purely sinus, and that has to be visible in the
    // archive rather than inferred later.
    enum class SeedBasis : uint8_t {
        EMPTY = 0,   // no beats at all
        SINUS_ONLY = 1,   // the intended path
        SINUS_PLUS_VOTED = 2,   // too few strictly-unflagged beats; voted
        //   beats readmitted (inferred ectopy, weaker
        //   evidence than a direct prematurity flag)
        ALL_USABLE = 3    // ectopy is the majority or nearly so; the
        //   seed is contaminated and says so
    };

    inline const char* seedBasisName(SeedBasis b) {
        switch (b) {
        case SeedBasis::SINUS_ONLY:       return "sinus_only";
        case SeedBasis::SINUS_PLUS_VOTED: return "sinus_plus_voted";
        case SeedBasis::ALL_USABLE:       return "all_usable";
        default:                          return "empty";
        }
    }

    struct SeedSelection {
        // Indices into whatever beat vector the caller passed. Parallel to
        // nothing else -- callers index their own arrays with these.
        std::vector<uint32_t> members;

        SeedBasis basis = SeedBasis::EMPTY;

        // Population accounting for BinCounts. n_excluded_* are the reasons a
        // beat was kept out, and they are not mutually exclusive with each
        // other in principle, so each is counted where it was decided.
        uint32_t n_candidates = 0;
        uint32_t n_excluded_premature = 0;
        uint32_t n_excluded_voted = 0;
        uint32_t n_excluded_category = 0;
        uint32_t n_selected = 0;

        // Ectopic fraction of the candidate population, computed BEFORE the
        // fallback ladder. This is the number that predicts whether Tukey
        // still works in this bin: as it approaches 0.5 the RR-length IQR
        // widens across both clusters and the fences stop rejecting anything,
        // including real artifact. Pair it with the fence IQR from BinCounts.
        double ectopic_fraction = 0.0;

        // True while the three Tukey passes still run upstream over all
        // sliced beats rather than inside this pool. Left as a field rather
        // than a comment so the condition is queryable from the archive: the
        // fences that trimmed these beats were computed over a mixed
        // population, so they were wider than they should have been.
        bool tukey_relocation_pending = true;

        bool contaminated() const { return basis == SeedBasis::ALL_USABLE; }
    };

    // rhythm: 0 NORMAL, 1 PVC (premature), 2 VOTED_PVC -- the kept_rhythm
    //         encoding from create_ecg_templates.hpp. Pass an empty vector
    //         when unavailable; the selection then rests on `category` alone.
    // category: tbank::Category per beat, derived from operator marks. Always
    //         populated in normal use, because the absence of a mark is the
    //         REGULAR verdict -- there is no unknown state. An empty vector is
    //         accepted only for the case where marks have not been loaded at
    //         all, and then every beat reads REGULAR and the selection rests
    //         on the rhythm flags alone.
    //
    // `candidates` are the indices already cleared by upstream filters -- in
    // create_ecg_templates.hpp's terms, usableIdx after the baseline_source
    // test.
    inline SeedSelection selectSeedPool(
        const std::vector<uint32_t>& candidates,
        const std::vector<uint8_t>& rhythm,
        const std::vector<tbank::Category>& category)
    {
        SeedSelection out;
        out.n_candidates = static_cast<uint32_t>(candidates.size());
        if (candidates.empty()) return out;

        const bool haveRhythm = !rhythm.empty();
        const bool haveCat = !category.empty();

        auto rhythmOf = [&](uint32_t i) -> uint8_t {
            return (haveRhythm && i < rhythm.size()) ? rhythm[i] : 0u;
            };
        auto categoryOf = [&](uint32_t i) -> tbank::Category {
            return (haveCat && i < category.size()) ? category[i]
                : tbank::Category::REGULAR;
            };

        std::vector<uint32_t> strict;      // category REGULAR, no rhythm flag
        std::vector<uint32_t> plus_voted;  // strict + voted-only beats
        uint32_t n_ectopic_evidence = 0;

        for (uint32_t i : candidates) {
            const tbank::Category cat = categoryOf(i);
            if (cat == tbank::Category::ECTOPIC || cat == tbank::Category::NOISE) {
                ++out.n_excluded_category;
                ++n_ectopic_evidence;
                continue;
            }
            const uint8_t rh = rhythmOf(i);
            if (rh == 1u) {                 // direct prematurity evidence
                ++out.n_excluded_premature;
                ++n_ectopic_evidence;
                continue;
            }
            if (rh == 2u) {                 // inferred by the 5-of-8 vote
                ++out.n_excluded_voted;
                ++n_ectopic_evidence;
                plus_voted.push_back(i);
                continue;
            }
            strict.push_back(i);
            plus_voted.push_back(i);
        }

        out.ectopic_fraction = static_cast<double>(n_ectopic_evidence)
            / static_cast<double>(candidates.size());

        // Fallback ladder. Never return an empty pool when candidates exist:
        // an empty seed means no template at all for the bin, which is a worse
        // outcome than a contaminated one that is labeled as contaminated.
        //
        // Voted beats are readmitted before unflagged-only is abandoned
        // because the vote is inferred ectopy -- weaker evidence than a direct
        // prematurity flag. Note the vote's own blind spot: in perfect
        // bigeminy exactly 4 of any 8 beats are ectopic, the count never
        // reaches 5, and the vote never fires, so in the bins where this
        // fallback matters most the voted set is often empty anyway.
        if (static_cast<int>(strict.size()) >= kMinSeedBeats) {
            out.members = std::move(strict);
            out.basis = SeedBasis::SINUS_ONLY;
        }
        else if (static_cast<int>(plus_voted.size()) >= kMinSeedBeats) {
            out.members = std::move(plus_voted);
            out.basis = SeedBasis::SINUS_PLUS_VOTED;
        }
        else if (!strict.empty()) {
            // Fewer than kMinSeedBeats either way, but some unflagged beats
            // exist: prefer the clean minority over the contaminated majority.
            out.members = std::move(strict);
            out.basis = SeedBasis::SINUS_ONLY;
        }
        else {
            out.members = candidates;
            out.basis = SeedBasis::ALL_USABLE;
        }

        out.n_selected = static_cast<uint32_t>(out.members.size());
        return out;
    }

}  // namespace seed_pool
