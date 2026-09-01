#pragma once
/**
 * @file   joint_bank.hpp
 * @brief  Section 4.6 morphology segregation as ONE partition shared by all
 *         four channels (CH1, CH2, CH3, PPG), rather than four independent
 *         banks.
 *
 *         WHY THIS REPLACES THE PER-CHANNEL BANK. A beat has one morphology.
 *         Splitting each channel separately produces four unrelated groupings
 *         of the same beats, and then nothing in the record says which of CH1's
 *         templates corresponds to which of CH2's -- so a PVC is separated
 *         three times over with no link between the three, and its pulse is
 *         separated a fourth time for unrelated perfusion reasons. Slot B of
 *         CH1 and slot B of PPG were different sets of beats that happened to
 *         share an index. This file makes them the same set by construction:
 *         a template IS a group of beats, and every channel contributes its own
 *         average over that one group.
 *
 *         ------------------------------------------------------------------
 *         THE SHARED KEY IS THE SLICE INDEX, NOT A CHANNEL'S BEAT INDEX
 *         ------------------------------------------------------------------
 *         Each channel prunes independently -- alignment's Tukey pass on ECG,
 *         the fit-error threshold on PPG -- so beat 12 of CH1 and beat 12 of
 *         PPG are different beats. A group therefore holds SLICE indices (the
 *         R-pair ordinal every slicer is driven by), and each channel supplies
 *         `local_of_slice`, mapping a slice index to that channel's own row or
 *         -1 when the channel dropped it.
 *
 *         A MEMBER NEED NOT BE PRESENT ON EVERY CHANNEL, and this is the common
 *         case, not an edge case: a beat whose ECG survived pruning but whose
 *         PPG did not is a full member of its group with no PPG sample. Its
 *         ECG average includes it and its PPG average does not. Assuming
 *         otherwise -- dropping such beats, or treating absence as a mismatch
 *         -- would either discard most of the record or split it on pruning
 *         luck rather than morphology.
 *
 *         ------------------------------------------------------------------
 *         SPAWN RULE: EITHER CHANNEL FAILS
 *         ------------------------------------------------------------------
 *         A beat joins a group only if it clears the floor on EVERY channel
 *         where both it and the group are scorable -- 0.85 for the three ECG
 *         leads, 0.80 for PPG, the two numbers the spec names. One failing
 *         channel is enough to reject the group, and a beat that no group
 *         accepts opens a new one.
 *
 *         STILL ONE THRESHOLD PER CHANNEL. The spec forbids a second, looser
 *         assignment threshold, and there isn't one: the same floor decides
 *         "this beat belongs to that group" and "no new group is needed". What
 *         is new is the conjunction across channels, which adds no number.
 *
 *         UNSCORABLE IS NOT FAILURE. A channel with no corridor yet, or too
 *         little overlap, abstains rather than voting no -- kUnscorable and
 *         "matched nothing" have to stay distinct, or a beat that merely could
 *         not be compared spawns a template. A group is a candidate only if at
 *         least one channel actually scored it; a beat no channel can score is
 *         returned as unscorable and assigned nowhere.
 *
 *         Consequence worth stating plainly: the conjunction splits more
 *         readily than any single channel would, and PPG is the noisiest of the
 *         four. Expect the cap to be reached more often than it was with ECG
 *         alone, which means the merge path and the confirmed-member cap raise
 *         matter more here than they did before, not less.
 *
 *         ------------------------------------------------------------------
 *         WHAT IS REUSED RATHER THAN REWRITTEN
 *         ------------------------------------------------------------------
 *         tbank::bandMatch and tbank::recomputeTemplate are per-channel and
 *         already correct, so each channel's face of a group is carried as a
 *         tbank::BankTemplate purely to hold {tmpl, tmpl_iqr, band_lo, band_hi}
 *         and be fed to those two functions.
 *
 *         THE PER-CHANNEL BankTemplate's IDENTITY FIELDS ARE NOT AUTHORITATIVE.
 *         label_code, confirmed_by_operator, subtype and spawn_seq exist on it
 *         because it is a whole BankTemplate, but a group's class is a property
 *         of the GROUP -- there is one label for one morphology, not four.
 *         BeatGroup carries them, and nothing should read them off ch[i].
 *         Reading the copy is how the four channels would come to disagree
 *         about what class a beat is.
 */

#include "template_bank.hpp"
#include "template_assign.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <vector>

namespace jbank {

    // Channel indices. PPG is the fourth, deliberately after the ECG leads so
    // that a loop over 0..2 is still "the ECG leads" and reads the same as it
    // does everywhere else in the pipeline.
    inline constexpr int kCh1 = 0, kCh2 = 1, kCh3 = 2, kPpg = 3;
    inline constexpr int kNumChannels = 4;
    inline constexpr int kNumEcgCh = 3;      // 0..2; kPpg is the fourth

    // The two floors Section 4.6 names. Indexed by channel so the assignment
    // loop never has to ask "is this the pulse channel".
    inline constexpr double kFloor[kNumChannels] = {
        tbank::kMatchFloorEcg, tbank::kMatchFloorEcg, tbank::kMatchFloorEcg,
        tbank::kMatchFloorPpg
    };

    inline const char* channelName(int c) {
        switch (c) {
        case kCh1: return "CH1";
        case kCh2: return "CH2";
        case kCh3: return "CH3";
        case kPpg: return "PPG";
        }
        return "?";
    }

    // ---------------------------------------------------------------------
    // One channel's beats, plus the map from the shared key into them
    // ---------------------------------------------------------------------
    struct ChannelBeats {
        const std::vector<std::vector<double>>* beats = nullptr;  // aligned, shared axis
        int width = 0;
        int anchor_col = -1;   // r_col for ECG, systolic peak column for PPG

        // local_of_slice[slice] = row in *beats, or -1 if this channel dropped
        // that slice. Sized to the slice count, NOT to beats->size().
        //
        // THE WHOLE DESIGN RESTS ON THIS VECTOR. Without it there is no way to
        // say that a PPG row and an ECG row are the same heartbeat, and the
        // partition cannot be joint. create_arterial_templates.hpp does not
        // currently produce one for PPG -- it discards the R-pair ordinal when
        // it filters on fit error -- so that has to be retained before this
        // file can be driven for real.
        std::vector<int32_t> local_of_slice;

        bool present() const { return beats && !beats->empty() && width > 0; }

        int localFor(uint32_t slice) const {
            if (slice >= local_of_slice.size()) return -1;
            return local_of_slice[slice];
        }
        const std::vector<double>* beatFor(uint32_t slice) const {
            const int k = localFor(slice);
            if (k < 0 || !beats || static_cast<size_t>(k) >= beats->size())
                return nullptr;
            return &(*beats)[static_cast<size_t>(k)];
        }
    };

    using ChannelSet = std::array<ChannelBeats, kNumChannels>;

    // ---------------------------------------------------------------------
    // Building local_of_slice from what the pipeline already produces
    // ---------------------------------------------------------------------
    //
    // BOTH SIDES ALREADY CARRY THE MAP, JUST INVERTED. alignment's BeatSet has
    // original_index[alignedRow] = slice, and the pulse path now has
    // PPGTemplatesResult::keptSlices[bin][keptRow] = slice. Either is the
    // forward direction; ChannelBeats wants the reverse, slice -> local row,
    // because assignment walks slices and asks each channel what it has.
    //
    // n_slices is the R-pair count for the bin (rPeaks.size() - 1), which is
    // the SAME denominator for every channel -- that is what makes the slice
    // index a shared key rather than four unrelated local ones. Passing a
    // per-channel beat count here instead would silently produce a map that
    // works for one channel and misaligns the rest.
    //
    // A slice no channel retained simply has -1 everywhere, which assignment
    // skips. A slice appearing twice cannot happen (a row maps to one slice),
    // but if a stale map ever produced it the LAST row wins, which is why the
    // loop does not break early -- a silent duplicate should not depend on
    // iteration order.
    template <typename IndexVec>
    inline std::vector<int32_t> localOfSlice(const IndexVec& forward,
        uint32_t n_slices)
    {
        std::vector<int32_t> out(n_slices, -1);
        for (size_t row = 0; row < forward.size(); ++row) {
            const uint64_t slice = static_cast<uint64_t>(forward[row]);
            if (slice < n_slices) out[static_cast<size_t>(slice)] =
                static_cast<int32_t>(row);
        }
        return out;
    }

    // Convenience: fill one channel from a forward map. anchor_col is r_col for
    // the ECG leads and the systolic peak column for PPG.
    template <typename IndexVec>
    inline void setChannel(ChannelSet& set, int channel,
        const std::vector<std::vector<double>>& beats,
        const IndexVec& forward, uint32_t n_slices, int anchor_col)
    {
        if (channel < 0 || channel >= kNumChannels) return;
        ChannelBeats& cb = set[channel];
        cb.beats = &beats;
        cb.width = beats.empty() ? 0 : static_cast<int>(beats.front().size());
        cb.anchor_col = anchor_col;
        cb.local_of_slice = localOfSlice(forward, n_slices);
    }

    // ---------------------------------------------------------------------
    // A group: one set of beats, four faces
    // ---------------------------------------------------------------------
    struct BeatGroup {
        // Shared membership, in SLICE indices. This is the partition.
        std::vector<uint32_t> members;

        // Per-channel waveform + corridor over `members`. Only the waveform
        // fields are meaningful; see the header note on identity fields.
        std::array<tbank::BankTemplate, kNumChannels> ch;

        // ---- group identity: ONE class for ONE morphology ----------------
        uint8_t  label_code = tbank::kUnlabeled;
        bool     confirmed_by_operator = false;   // never inferred from label_code
        int32_t  subtype = -1;
        uint32_t spawn_seq = 0;

        std::map<int32_t, tbank::BankMarkerSet> markers_by_anchor;

        tbank::BankMarkerSet& marks(int32_t a) { return markers_by_anchor[a]; }

        bool confirmed() const { return confirmed_by_operator; }
        bool isSeed()    const { return spawn_seq == 0; }
        int  memberCount() const { return static_cast<int>(members.size()); }
        bool earnsColumn() const {
            return memberCount() >= tbank::kMinMembersForColumn;
        }
        bool isJunk() const { return memberCount() <= tbank::kMaxJunkMembers; }

        // Members this channel actually has. Lower than memberCount() whenever
        // the channel pruned some of them, which is normal.
        int memberCountOn(int c) const {
            return (c >= 0 && c < kNumChannels)
                ? static_cast<int>(ch[c].members.size()) : 0;
        }
    };

    // ---------------------------------------------------------------------
    // The bank
    // ---------------------------------------------------------------------
    struct JointBank {
        std::vector<BeatGroup> groups;

        int32_t configured_cap = tbank::kDefaultMaxTemplatesPerBin;
        int32_t effective_cap = tbank::kDefaultMaxTemplatesPerBin;
        uint32_t next_spawn_seq = 0;
        uint32_t assigned_beats = 0;

        void setCap(int32_t cap) {
            if (cap < 1) cap = 1;          // slot 0 always exists
            configured_cap = cap;
            effective_cap = cap;
        }

        int  size()  const { return static_cast<int>(groups.size()); }
        bool atCap() const { return size() >= effective_cap; }

        int nConfirmed() const {
            int n = 0;
            for (const auto& g : groups) if (g.confirmed()) ++n;
            return n;
        }

        // Which group holds this slice, or -1.
        int findBySlice(uint32_t slice) const {
            for (int i = 0; i < size(); ++i)
                for (uint32_t m : groups[i].members)
                    if (m == slice) return i;
            return -1;
        }

        // Next subtype index for a class, in order of first appearance. Same
        // rule as the per-channel bank: the operator never supplies one.
        int32_t nextSubtypeFor(uint8_t label_code) const {
            int32_t hi = 0;
            for (const auto& g : groups)
                if (g.label_code == label_code && g.subtype > hi) hi = g.subtype;
            return hi + 1;
        }
    };

    // ---------------------------------------------------------------------
    // Recomputation
    // ---------------------------------------------------------------------

    // Rebuild one channel's face of a group: project the group's slice members
    // into that channel's local rows, then hand the result to the existing
    // per-channel median/corridor builder.
    inline void recomputeGroupChannel(BeatGroup& g, int c,
        const ChannelSet& chans,
        const std::vector<double>* floor_corridor = nullptr)
    {
        const ChannelBeats& cb = chans[c];
        tbank::BankTemplate& t = g.ch[c];
        t.members.clear();
        if (!cb.present()) {
            t.tmpl.clear(); t.tmpl_iqr.clear();
            t.band_lo.clear(); t.band_hi.clear();
            return;
        }
        t.members.reserve(g.members.size());
        for (uint32_t slice : g.members) {
            const int k = cb.localFor(slice);
            if (k >= 0) t.members.push_back(static_cast<uint32_t>(k));
        }
        t.r_col = cb.anchor_col;
        tbank::recomputeTemplate(t, *cb.beats, cb.width, floor_corridor);
    }

    inline void recomputeGroup(BeatGroup& g, const ChannelSet& chans,
        const std::array<std::vector<double>, kNumChannels>* seed_corridors = nullptr)
    {
        for (int c = 0; c < kNumChannels; ++c) {
            const std::vector<double>* fc =
                (seed_corridors && !(*seed_corridors)[c].empty())
                ? &(*seed_corridors)[c] : nullptr;
            recomputeGroupChannel(g, c, chans, fc);
        }
    }

    inline void recomputeAll(JointBank& bank, const ChannelSet& chans)
    {
        // Slot 0 first and with no inherited corridor: it IS the corridor every
        // young group inherits, so it has to exist before they are rebuilt.
        if (!bank.groups.empty())
            recomputeGroup(bank.groups[0], chans, nullptr);

        std::array<std::vector<double>, kNumChannels> seeds;
        if (!bank.groups.empty()) {
            for (int c = 0; c < kNumChannels; ++c) {
                const tbank::BankTemplate& s0 = bank.groups[0].ch[c];
                const size_t w = std::min(s0.band_lo.size(), s0.band_hi.size());
                seeds[c].assign(w, std::numeric_limits<double>::quiet_NaN());
                for (size_t k = 0; k < w; ++k)
                    if (!std::isnan(s0.band_lo[k]) && !std::isnan(s0.band_hi[k]))
                        seeds[c][k] = 0.5 * (s0.band_hi[k] - s0.band_lo[k]);
            }
        }
        for (int i = 1; i < bank.size(); ++i)
            recomputeGroup(bank.groups[i], chans, &seeds);
    }

    // ---------------------------------------------------------------------
    // Scoring one beat against one group
    // ---------------------------------------------------------------------

    struct JointScore {
        double mean_score = std::numeric_limits<double>::quiet_NaN();
        int    n_scored = 0;      // channels that could be compared
        int    failing_channel = -1;   // first channel below its floor, or -1
        double per_channel[kNumChannels] = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN() };

        // Accepted only if something scored it and nothing failed.
        bool accepts() const { return n_scored > 0 && failing_channel < 0; }
    };

    // EVERY channel is scored even after one fails, rather than breaking out at
    // the first failure. The cost is small (four band matches on a few hundred
    // columns) and the diagnostic value is the point: "rejected on PPG while
    // passing all three ECG leads at 0.97" is the row that tells you the pulse
    // channel is driving your spawn rate, and a short-circuit would throw that
    // away. failing_channel keeps the FIRST failure so the summary stays a
    // single number.
    inline JointScore scoreAgainst(const BeatGroup& g, uint32_t slice,
        const ChannelSet& chans)
    {
        JointScore js;
        double sum = 0.0;
        for (int c = 0; c < kNumChannels; ++c) {
            const std::vector<double>* beat = chans[c].beatFor(slice);
            if (!beat) continue;                 // channel dropped this beat
            const tbank::BandResult br = tbank::bandMatch(*beat, g.ch[c]);
            if (!br.scorable()) continue;        // abstains; not a failure
            js.per_channel[c] = br.score;
            sum += br.score;
            ++js.n_scored;
            if (br.score < kFloor[c] && js.failing_channel < 0)
                js.failing_channel = c;
        }
        if (js.n_scored > 0) js.mean_score = sum / js.n_scored;
        return js;
    }

    // ---------------------------------------------------------------------
    // Merge selection, with the confirmed-member rule
    // ---------------------------------------------------------------------

    struct MergeCandidate {
        int    a = -1, b = -1;
        double closeness = -1.0;
        bool   both_confirmed = false;   // true => raise the cap instead
        bool   garbage_pair = false;
        bool   valid() const { return a >= 0 && b >= 0; }
    };

    // Closeness between two groups: the mean over channels of the symmetrized
    // band match between their medians. Symmetrized by max because band-match
    // is directional -- a wide group swallows a narrow one asymmetrically --
    // and averaged over channels because the partition is joint, so "the two
    // closest" has to mean closest overall rather than closest on one lead.
    inline double groupCloseness(const BeatGroup& x, const BeatGroup& y)
    {
        double sum = 0.0; int n = 0;
        for (int c = 0; c < kNumChannels; ++c) {
            if (x.ch[c].tmpl.empty() || y.ch[c].tmpl.empty()) continue;
            const double ab = tbank::bandMatch(x.ch[c].tmpl, y.ch[c]).score;
            const double ba = tbank::bandMatch(y.ch[c].tmpl, x.ch[c]).score;
            double s;
            if (std::isnan(ab) && std::isnan(ba)) continue;
            else if (std::isnan(ab)) s = ba;
            else if (std::isnan(ba)) s = ab;
            else s = std::max(ab, ba);
            sum += s; ++n;
        }
        return (n > 0) ? sum / n : std::numeric_limits<double>::quiet_NaN();
    }

    // A pair containing ANY confirmed group is blocked, not merely a pair where
    // both are confirmed. The spec's words are "both", and this is stricter on
    // purpose: merging a confirmed group into an unconfirmed one absorbs
    // unlabeled beats into a labeled class, which is the label-by-similarity
    // inference the section forbids elsewhere. The cost is more cap raises than
    // the spec's literal rule would produce, which is the safe direction.
    inline MergeCandidate findMergePair(const JointBank& bank)
    {
        MergeCandidate best_any, best_blocked;
        const int n = bank.size();
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double s = groupCloseness(bank.groups[i], bank.groups[j]);
                if (std::isnan(s)) continue;
                if (bank.groups[i].confirmed() || bank.groups[j].confirmed()) {
                    if (s > best_blocked.closeness) {
                        best_blocked.a = i; best_blocked.b = j;
                        best_blocked.closeness = s;
                        best_blocked.both_confirmed = true;
                    }
                    continue;
                }
                if (s > best_any.closeness) {
                    best_any.a = i; best_any.b = j;
                    best_any.closeness = s;
                    best_any.garbage_pair = bank.groups[i].isJunk()
                        && bank.groups[j].isJunk();
                }
            }
        }
        if (best_any.valid()) return best_any;
        return best_blocked;
    }

    // Fold b into a and drop b. NO LABEL IS COPIED: findMergePair blocks any
    // pair touching a confirmed group, so this is unreachable for one, and the
    // guard stays anyway because a caller arriving by another route must not
    // acquire a class it was never confirmed as.
    inline void mergeGroups(JointBank& bank, int a, int b,
        const ChannelSet& chans)
    {
        if (a < 0 || b < 0 || a >= bank.size() || b >= bank.size() || a == b)
            return;
        if (a > b) std::swap(a, b);

        BeatGroup& dst = bank.groups[a];
        BeatGroup& src = bank.groups[b];
        dst.members.insert(dst.members.end(),
            src.members.begin(), src.members.end());
        std::sort(dst.members.begin(), dst.members.end());
        // The LOWER spawn_seq survives, so seed identity (spawn_seq 0) survives
        // merge history without a stored flag -- same contract as the
        // per-channel bank's isSeed().
        dst.spawn_seq = std::min(dst.spawn_seq, src.spawn_seq);

        bank.groups.erase(bank.groups.begin() + b);
        recomputeAll(bank, chans);
    }

    // ---------------------------------------------------------------------
    // Assignment
    // ---------------------------------------------------------------------

    struct AssignOutcome {
        int    group_id = tbank::kNoMatch;
        double mean_score = std::numeric_limits<double>::quiet_NaN();
        int    failing_channel = -1;   // on the best REJECTED group, for diagnostics
        bool   spawned = false;
        bool   merged = false;
        bool   cap_raised = false;
    };

    struct BankCounts {
        uint32_t slices_seen = 0;
        uint32_t n_spawns = 0;
        uint32_t n_merges = 0;
        uint32_t n_cap_raises = 0;
        uint32_t n_unscorable = 0;
        // Which channel rejected the winning group, per channel. This is the
        // number that says whether the conjunction is splitting on cardiac
        // morphology or on pulse noise, and it is the first thing to look at if
        // the cap is being hit on ordinary records.
        uint32_t n_rejected_by[kNumChannels] = { 0, 0, 0, 0 };
    };

    inline AssignOutcome assignSlice(JointBank& bank, uint32_t slice,
        const ChannelSet& chans,
        uint64_t bin_index = 0,
        std::vector<tbank::CapRaiseEvent>* events = nullptr,
        BankCounts* counts = nullptr)
    {
        AssignOutcome out;
        if (counts) ++counts->slices_seen;

        // ---- best accepting group ---------------------------------------
        int best = -1; double best_mean = -1.0;
        int anyScored = 0;
        JointScore bestRejected;
        for (int i = 0; i < bank.size(); ++i) {
            const JointScore js = scoreAgainst(bank.groups[i], slice, chans);
            if (js.n_scored > 0) ++anyScored;
            if (js.accepts()) {
                if (js.mean_score > best_mean) { best_mean = js.mean_score; best = i; }
            }
            else if (js.n_scored > 0
                && (std::isnan(bestRejected.mean_score)
                    || js.mean_score > bestRejected.mean_score)) {
                bestRejected = js;
            }
        }

        if (best >= 0) {
            bank.groups[best].members.push_back(slice);
            ++bank.assigned_beats;
            out.group_id = best;
            out.mean_score = best_mean;
            return out;
        }

        // A beat NO group could score is not a new morphology, it is an
        // unmeasurable beat. Spawning on it is how a bin fills with
        // single-member templates built from noise.
        if (anyScored == 0 && bank.size() > 0) {
            if (counts) ++counts->n_unscorable;
            out.group_id = tbank::kUnscorable;
            return out;
        }

        out.failing_channel = bestRejected.failing_channel;
        if (counts && out.failing_channel >= 0
            && out.failing_channel < kNumChannels)
            ++counts->n_rejected_by[out.failing_channel];

        // ---- spawn ------------------------------------------------------
        if (bank.atCap()) {
            const MergeCandidate mc = findMergePair(bank);
            if (mc.valid() && !mc.both_confirmed) {
                mergeGroups(bank, mc.a, mc.b, chans);
                out.merged = true;
                if (counts) ++counts->n_merges;
            }
            else {
                const int32_t old_cap = bank.effective_cap;
                ++bank.effective_cap;
                out.cap_raised = true;
                if (counts) ++counts->n_cap_raises;

                // LOGGED EVEN WHEN THERE IS NO CANDIDATE PAIR. mc can come back
                // invalid -- groupCloseness is NaN for every pair when no two
                // groups share a scorable channel, which a single-channel
                // dataset can reach early in a bin. The cap still has to rise,
                // and an unexplained cap rise is worse than a partially
                // populated event: the census would show effective_cap above
                // configured_cap with nothing saying why. The pair fields stay
                // -1/NaN, which is what "no pair was in the running" looks like.
                if (events) {
                    tbank::CapRaiseEvent e;
                    e.bin_index = bin_index;
                    e.channel = -1;   // joint: not attributable to one channel
                    e.old_cap = old_cap;
                    e.new_cap = bank.effective_cap;
                    if (mc.valid()) {
                        e.template_a = mc.a;
                        e.template_b = mc.b;
                        e.closeness = mc.closeness;
                        e.label_a = bank.groups[mc.a].label_code;
                        e.label_b = bank.groups[mc.b].label_code;
                    }
                    events->push_back(e);
                }
            }
        }

        BeatGroup g;
        g.spawn_seq = bank.next_spawn_seq++;
        g.members.push_back(slice);
        bank.groups.push_back(std::move(g));
        ++bank.assigned_beats;
        out.group_id = bank.size() - 1;
        out.spawned = true;
        if (counts) ++counts->n_spawns;

        // Rebuilt through recomputeAll so the new group inherits slot 0's
        // corridor on every channel. A one-member group has no spread of its
        // own -- lo == hi at every column -- so without the inherited corridor
        // the next beat of the same morphology scores ~0 and spawns again.
        recomputeAll(bank, chans);
        return out;
    }

    // ---------------------------------------------------------------------
    // Seed + drive
    // ---------------------------------------------------------------------

    // Slot 0 from Phase 1. One group, seeded per channel from that channel's
    // Phase 1 template: the ECG leads from chN's sinus template, PPG from the
    // bin's pulse template. `seed_slices` are the slices Phase 1 built those
    // templates from, and they become slot 0's initial membership.
    inline void seedBank(JointBank& bank, const ChannelSet& chans,
        const std::array<std::vector<double>, kNumChannels>& phase1,
        const std::vector<uint32_t>& seed_slices,
        int32_t max_templates_per_bin = 0)
    {
        bank.groups.clear();
        bank.next_spawn_seq = 0;
        bank.assigned_beats = 0;
        if (max_templates_per_bin > 0) bank.setCap(max_templates_per_bin);

        BeatGroup g;
        g.spawn_seq = bank.next_spawn_seq++;   // 0: the seed
        g.members = seed_slices;
        bank.groups.push_back(std::move(g));

        recomputeAll(bank, chans);

        // Phase 1's waveform overrides the recomputed median where it exists,
        // so features measured in Phase 1 and features measured against slot 0
        // use ONE reference. The corridor stays as recomputed: Phase 1 supplies
        // a waveform, not a spread.
        for (int c = 0; c < kNumChannels; ++c)
            if (!phase1[c].empty())
                bank.groups[0].ch[c].tmpl = phase1[c];

        bank.assigned_beats =
            static_cast<uint32_t>(bank.groups[0].members.size());
    }

    // Every slice in order. `n_slices` is the R-pair count for the bin, which
    // is the length local_of_slice was sized to on every channel.
    inline void runBank(JointBank& bank, const ChannelSet& chans,
        uint32_t n_slices,
        std::vector<int32_t>& out_group_of_slice,
        uint64_t bin_index = 0,
        std::vector<tbank::CapRaiseEvent>* events = nullptr,
        BankCounts* counts = nullptr)
    {
        out_group_of_slice.assign(n_slices, tbank::kNoMatch);
        for (const uint32_t m : bank.groups.empty()
            ? std::vector<uint32_t>{} : bank.groups[0].members)
            if (m < n_slices) out_group_of_slice[m] = 0;

        for (uint32_t s = 0; s < n_slices; ++s) {
            if (out_group_of_slice[s] >= 0) continue;   // seeded
            // A slice no channel retained is not a beat as far as this bank is
            // concerned; skipped rather than counted as unscorable, which would
            // inflate that diagnostic with beats nobody tried to measure.
            bool anyPresent = false;
            for (int c = 0; c < kNumChannels && !anyPresent; ++c)
                if (chans[c].beatFor(s)) anyPresent = true;
            if (!anyPresent) continue;

            const AssignOutcome ao = assignSlice(bank, s, chans,
                bin_index, events, counts);
            // Recorded, but NOT trusted -- see the rebuild below.
            out_group_of_slice[s] = ao.group_id;
        }
        recomputeAll(bank, chans);

        // ---- REBUILD THE ASSIGNMENT FROM MEMBERSHIP -----------------------
        // The ids returned by assignSlice go STALE. mergeGroups erases an
        // element, so every group after it shifts down by one, and any id
        // recorded before that merge now points at the wrong group. On a bin
        // that merges even once, the earlier half of the assignment vector is
        // silently off -- and it is off by a variable amount, so it looks like
        // a plausible partition rather than an error.
        //
        // Membership is the single source of truth (it is what merging
        // actually moves), so the vector is rederived from it after the run
        // rather than patched at each merge. Anything not claimed by a group
        // keeps whatever assignSlice returned, which is where kUnscorable and
        // kNoMatch survive.
        for (int gi = 0; gi < bank.size(); ++gi)
            for (const uint32_t m : bank.groups[gi].members)
                if (m < n_slices) out_group_of_slice[m] = gi;
    }

    // ---------------------------------------------------------------------
    // Label propagation
    // ---------------------------------------------------------------------

    // ONE call covers every channel, which is the point of a joint partition:
    // "the class label propagates to the template that beat is assigned to, and
    // from there to every other beat assigned to the same template" is now a
    // single statement rather than four that have to be kept consistent.
    //
    // Takes a SLICE index, not a channel's beat index. An operator clicking a
    // beat on the CH2 panel is confirming a slice; the caller maps CH2's local
    // row back through local_of_slice before calling.
    struct PropagationResult {
        int32_t group = -1;
        int32_t subtype = -1;
        int     beats_relabeled = 0;
    };

    inline PropagationResult propagateLabel(JointBank& bank, uint32_t slice,
        uint8_t label_code)
    {
        PropagationResult out;
        if (label_code == tbank::kUnlabeled) return out;
        const int gi = bank.findBySlice(slice);
        if (gi < 0) return out;

        BeatGroup& g = bank.groups[gi];
        // Subtype issued once, then immutable: a re-confirmation, or a
        // confirmation of a different beat in an already-labeled group, must
        // not mint a new index.
        if (g.subtype < 0 || g.label_code != label_code) {
            if (g.label_code != label_code)
                g.subtype = bank.nextSubtypeFor(label_code);
            g.label_code = label_code;
        }
        g.confirmed_by_operator = true;

        out.group = gi;
        out.subtype = g.subtype;
        out.beats_relabeled = g.memberCount();
        return out;
    }

    // ---------------------------------------------------------------------
    // Presumed category, and the display letter
    // ---------------------------------------------------------------------

    // Presumed, NOT confirmed. A confirmed label overrides the presumption
    // outright -- the operator's verdict is not a hypothesis to be re-derived.
    // Unlabeled is PQRST, not a third "unknown" state.
    inline tbank::Category presumedCategory(const BeatGroup& g,
        uint32_t n_premature_members = 0, uint32_t n_voted_members = 0)
    {
        if (g.confirmed_by_operator) {
            if (g.label_code == tbank::kCodeMinorNoise) return tbank::Category::NOISE;
            if (g.label_code == tbank::kCodePvc || g.label_code == tbank::kCodePac
                || g.label_code == tbank::kCodeVt) return tbank::Category::ECTOPIC;
            return tbank::Category::REGULAR;
        }
        if (!g.earnsColumn()) return tbank::Category::NOISE;
        const uint32_t n = static_cast<uint32_t>(g.members.size());
        if (n == 0) return tbank::Category::NOISE;
        if ((n_premature_members + n_voted_members) * 2 > n)
            return tbank::Category::ECTOPIC;
        return tbank::Category::REGULAR;
    }

    // Contiguous letters over the groups that exist, ordered by spawn_seq --
    // the joint counterpart of tbank::letterRanks, and for the same reason: a
    // merge erases an element and shifts everything after it, so slot position
    // is not stable between runs while spawn order is. One letter per group now
    // names the morphology on ALL FOUR channels, which is what the joint
    // partition buys.
    inline std::vector<uint8_t> letterRanks(const JointBank& bank) {
        const int n = bank.size();
        std::vector<int> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i)
            if (bank.groups[i].memberCount() > 0) order.push_back(i);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return bank.groups[a].spawn_seq < bank.groups[b].spawn_seq;
            });

        std::vector<uint8_t> letter(n, 0);
        int rank = 0;
        for (int i : order) {
            const BeatGroup& g = bank.groups[i];
            const int idx = (g.label_code != tbank::kUnlabeled && g.subtype > 0)
                ? g.subtype - 1 : rank;
            letter[i] = static_cast<uint8_t>(idx % 26);
            ++rank;
        }
        return letter;
    }

    // ---------------------------------------------------------------------
    // THE DRIVER: one joint bank per bin
    // ---------------------------------------------------------------------
    //
    // Takes exactly what the pipeline already has, per bin:
    //   ecg_beats[c]   res.chN.kept_beats_raw[bin]      (aligned, shared axis)
    //   ecg_forward[c] keptIdx[c][bin]                  (alignedRow -> slice)
    //   ecg_r_col[c]   res.chN.r_col_raw[bin]
    //   ppg_beats      ppg_res.kept[bin]
    //   ppg_forward    ppg_res.keptSlices[bin]
    //   ppg_peak_col   ppg_res.peakCol[bin]
    //   n_slices       rPeaks.size() - 1 for this bin
    //
    // n_slices IS THE SHARED DENOMINATOR and must be the R-pair count, not any
    // channel's beat count. Every forward map indexes into the same R-peak
    // vector; substituting a per-channel count here produces a map that lines
    // up for one channel and is silently wrong for the rest, which is the whole
    // failure mode the slice key exists to prevent.
    //
    // A channel with no beats contributes nothing and abstains from every
    // score -- that is how a 1-lead, no-PPG record (SHHS) degrades to the
    // single-channel behaviour without a special case.
    struct BinBankInput {
        const std::vector<std::vector<double>>* ecg_beats[kNumEcgCh] = { nullptr, nullptr, nullptr };
        const std::vector<size_t>* ecg_forward[kNumEcgCh] = { nullptr, nullptr, nullptr };
        int ecg_r_col[kNumEcgCh] = { -1, -1, -1 };
        std::vector<double> ecg_phase1[kNumEcgCh];

        const std::vector<std::vector<double>>* ppg_beats = nullptr;
        const std::vector<uint32_t>* ppg_forward = nullptr;
        int ppg_peak_col = -1;
        std::vector<double> ppg_phase1;

        uint32_t n_slices = 0;
        int32_t  max_templates_per_bin = 0;   // 0 => kDefaultMaxTemplatesPerBin
        uint64_t bin_index = 0;
    };

    struct BinBankOutput {
        JointBank bank;
        std::vector<int32_t> group_of_slice;   // sized n_slices; -1/kNoMatch/kUnscorable
        BankCounts counts;
        std::vector<tbank::CapRaiseEvent> cap_raises;
    };

    inline BinBankOutput buildBinBank(const BinBankInput& in)
    {
        BinBankOutput out;

        ChannelSet chans;
        for (int c = 0; c < kNumEcgCh; ++c) {
            if (!in.ecg_beats[c] || !in.ecg_forward[c]) continue;
            setChannel(chans, c, *in.ecg_beats[c], *in.ecg_forward[c],
                in.n_slices, in.ecg_r_col[c]);
        }
        if (in.ppg_beats && in.ppg_forward)
            setChannel(chans, kPpg, *in.ppg_beats, *in.ppg_forward,
                in.n_slices, in.ppg_peak_col);

        std::array<std::vector<double>, kNumChannels> phase1;
        for (int c = 0; c < kNumEcgCh; ++c) phase1[c] = in.ecg_phase1[c];
        phase1[kPpg] = in.ppg_phase1;

        // SEED MEMBERSHIP: the slices the seed can actually be scored against.
        // A slice present on at least one channel is eligible; a slice no
        // channel retained is not a beat as far as this bank is concerned.
        //
        // Seeded from the FIRST eligible slices rather than a sample across the
        // bin, matching the per-channel bank's behaviour, so slot 0's initial
        // corridor comes from a contiguous stretch. runBank then walks every
        // slice including these, and recomputeAll replaces the seed corridor
        // with slot 0's measured one once members accumulate.
        std::vector<uint32_t> seedSlices;
        const uint32_t kSeedTarget = 20;
        for (uint32_t s = 0; s < in.n_slices && seedSlices.size() < kSeedTarget; ++s) {
            bool any = false;
            for (int c = 0; c < kNumChannels && !any; ++c)
                if (chans[c].beatFor(s)) any = true;
            if (any) seedSlices.push_back(s);
        }

        seedBank(out.bank, chans, phase1, seedSlices, in.max_templates_per_bin);
        runBank(out.bank, chans, in.n_slices, out.group_of_slice,
            in.bin_index, &out.cap_raises, &out.counts);
        return out;
    }

    // ---------------------------------------------------------------------
    // PROJECTION: the joint partition, expressed in the existing structures
    // ---------------------------------------------------------------------
    //
    // THERE IS ONLY ONE MERGE, AND IT HAPPENS IN THIS FILE. The viewer, the
    // serializer and the morphology writers all read tbank::TemplateBank per
    // channel. Rather than migrate all of them at once -- and rather than leave
    // two independent partitions coexisting, which is the thing that must not
    // happen -- the joint groups are PROJECTED into those per-channel banks.
    //
    // So `ecg_bank[c]` keeps its type and every consumer keeps working, but it
    // is no longer a partition in its own right: it is channel c's view of the
    // one partition. Template i of every channel is group i. That is what makes
    // the columns line up across channels, which independent per-channel banks
    // never did.
    //
    // MEMBERS ARE IN THE CHANNEL'S OWN LOCAL INDEX SPACE, because that is what
    // every existing consumer expects of BankTemplate::members. The group's
    // slice membership is the authority; this is a translation of it, and a
    // slice the channel dropped simply does not appear.
    inline tbank::TemplateBank projectToChannel(const JointBank& bank,
        const ChannelSet& chans, int channel)
    {
        tbank::TemplateBank out;
        out.configured_cap = bank.configured_cap;
        out.effective_cap = bank.effective_cap;
        out.next_spawn_seq = bank.next_spawn_seq;
        if (channel < 0 || channel >= kNumChannels) return out;

        const ChannelBeats& cb = chans[channel];
        out.templates.reserve(bank.groups.size());
        for (const BeatGroup& g : bank.groups) {
            tbank::BankTemplate t = g.ch[channel];   // waveform + corridor
            t.members.clear();
            for (const uint32_t slice : g.members) {
                const int k = cb.localFor(slice);
                if (k >= 0) t.members.push_back(static_cast<uint32_t>(k));
            }
            // GROUP identity, not the per-channel copy's. The per-channel
            // BankTemplate carries label/confirmed/subtype fields because it is
            // a whole BankTemplate; only the group's are authoritative, and
            // this is the one place they are written onto it.
            t.label_code = g.label_code;
            t.confirmed_by_operator = g.confirmed_by_operator;
            t.subtype = g.subtype;
            t.spawn_seq = g.spawn_seq;
            t.markers_by_anchor = g.markers_by_anchor;
            t.n_ppg_members = static_cast<int32_t>(g.memberCountOn(kPpg));
            out.templates.push_back(std::move(t));
        }
        out.assigned_beats = bank.assigned_beats;
        return out;
    }

    // All four at once, rebuilding the ChannelSet the bank was run against.
    // The caller must pass the SAME ChannelSet, or the local index spaces will
    // not match the ones the members were resolved from.
    inline std::array<tbank::TemplateBank, kNumChannels> projectAll(
        const JointBank& bank, const ChannelSet& chans)
    {
        std::array<tbank::TemplateBank, kNumChannels> out;
        for (int c = 0; c < kNumChannels; ++c)
            out[c] = projectToChannel(bank, chans, c);
        return out;
    }

}  // namespace jbank
