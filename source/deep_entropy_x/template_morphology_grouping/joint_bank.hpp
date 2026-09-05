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
// keep_within_tukey and TukeyStats. THE ONLY TUKEY FENCE IN THE CODEBASE lives
// in alignment.hpp and this file consumes it rather than reimplementing the
// quartile arithmetic -- two fences would be two things to keep in agreement,
// and they would disagree first on the small groups that matter most.
#include "template_marking_gui/alignment.hpp"
// runFilter: the verbatim isPremature + 5-of-8 vote. Driven here from the
// per-slice RR series, so the verdict is indexed the same way the partition is.
#include "pvc_filter.hpp"
// substituteBeatNaNSafe + the 4.6 alpha and borderline band.
#include "beat_substitute.hpp"

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
    // Per-channel match floor, read at the point of use.
    //
    // WAS A constexpr ARRAY initialised from tbank's constants, which froze
    // both values at compile time -- so a config-set floor would have been
    // honoured by everything except the conjunction that actually decides the
    // partition. A function cannot go stale that way.
    inline double floorFor(int channel) {
        return (channel == kPpg) ? tbank::matchFloorPpg()
                                 : tbank::matchFloorEcg();
    }

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

        // ---- THE CLEAN SUBSET, ALSO IN SLICE INDICES --------------------
        //
        // members minus the premature and minus the Tukey-rejected, filled by
        // cleanGroups() after the partition is final. It is what every
        // channel's waveform is averaged from, what is drawn and marked, and
        // what "kept" means in the outputs. `members` keeps everything, because
        // an excluded beat is still this group's beat and the archive has to say
        // so and say why.
        //
        // ONE LIST FOR ALL FOUR CHANNELS, which is the point of a joint
        // partition: a beat is excluded from the group, not from CH2. Per
        // channel it would be four different subsets of one morphology, and the
        // three leads and the pulse would each average a different set of
        // heartbeats.
        //
        // Empty means "not yet cleaned" -- consumers fall back to `members`
        // rather than treating the group as empty.
        std::vector<uint32_t> members_clean;

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
        int  cleanCount() const {
            return static_cast<int>(members_clean.empty()
                ? members.size() : members_clean.size());
        }
        int  excludedCount() const {
            return static_cast<int>(members_clean.empty()
                ? 0 : members.size() - members_clean.size());
        }
        // The clean count, not the full one. A group whose two members were
        // both excluded has nothing behind its waveform and must not present as
        // a markable morphology.
        bool earnsColumn() const {
            return cleanCount() >= tbank::kMinMembersForColumn;
        }
        // Which list the waveform is built from.
        const std::vector<uint32_t>& averagedMembers() const {
            return members_clean.empty() ? members : members_clean;
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
        // averagedMembers(), NOT members. Before cleanGroups() runs the two are
        // the same vector, so this costs nothing then; afterwards it is the
        // whole mechanism by which an excluded beat leaves the waveform and the
        // corridor while staying in the group. No swapping of member lists
        // around the recompute, and no second recompute path that could average
        // a different set than the one on screen.
        const std::vector<uint32_t>& src = g.averagedMembers();
        t.members.reserve(src.size());
        for (uint32_t slice : src) {
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

    // Slot 0's corridor half-widths per channel, which every young group
    // inherits. Factored out because rebuilding ONE group needs them just as
    // much as rebuilding all of them does, and recomputing slot 0 to get them
    // is the entire cost this avoids.
    inline std::array<std::vector<double>, kNumChannels> inheritedCorridors(
        const JointBank& bank)
    {
        std::array<std::vector<double>, kNumChannels> seeds;
        if (bank.groups.empty()) return seeds;
        for (int c = 0; c < kNumChannels; ++c) {
            const tbank::BankTemplate& s0 = bank.groups[0].ch[c];
            const size_t w = std::min(s0.band_lo.size(), s0.band_hi.size());
            seeds[c].assign(w, std::numeric_limits<double>::quiet_NaN());
            for (size_t k = 0; k < w; ++k)
                if (!std::isnan(s0.band_lo[k]) && !std::isnan(s0.band_hi[k]))
                    seeds[c][k] = 0.5 * (s0.band_hi[k] - s0.band_lo[k]);
        }
        return seeds;
    }

    // ---- REBUILD ONE GROUP, NOT THE WHOLE BANK ---------------------------
    //
    // THIS IS THE DIFFERENCE BETWEEN A BIN AND AN AFTERNOON. A spawn changes
    // exactly one group -- the new one, with one member -- and leaves slot 0's
    // membership untouched, so its corridor is unchanged and nothing else needs
    // rebuilding. recomputeAll on every spawn instead rebuilt slot 0 over its
    // entire membership: a column-wise median over up to a thousand beats, on
    // four channels, across an 1800-column axis, per spawn. A bin that spawns a
    // couple of hundred times paid that a couple of hundred times, and the
    // profile looks like a hang rather than like slow code.
    //
    // Callers must NOT use this when slot 0 itself changed -- see mergeGroups.
    inline void recomputeOne(JointBank& bank, const ChannelSet& chans, int idx)
    {
        if (idx < 0 || idx >= bank.size()) return;
        if (idx == 0) { recomputeGroup(bank.groups[0], chans, nullptr); return; }
        const std::array<std::vector<double>, kNumChannels> seeds =
            inheritedCorridors(bank);
        recomputeGroup(bank.groups[idx], chans, &seeds);
    }

    inline void recomputeAll(JointBank& bank, const ChannelSet& chans)
    {
        // Slot 0 first and with no inherited corridor: it IS the corridor every
        // young group inherits, so it has to exist before they are rebuilt.
        if (!bank.groups.empty())
            recomputeGroup(bank.groups[0], chans, nullptr);

        const std::array<std::vector<double>, kNumChannels> seeds =
            inheritedCorridors(bank);
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
            if (br.score < floorFor(c) && js.failing_channel < 0)
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

    // Closeness between two groups: the mean over channels of the match between
    // their medians, averaged over channels because the partition is joint, so
    // "the two closest" has to mean closest overall rather than closest on one
    // lead.
    //
    // THE SYMMETRISATION IS NOW REDUNDANT AND IS KEPT ANYWAY. bandMatch scored
    // the fraction of one waveform inside the other's corridor, which is
    // directional -- a wide group swallowed a narrow one asymmetrically -- so
    // the two directions were taken and the larger used. Correlation is
    // symmetric, so ab and ba are equal to within floating-point noise and the
    // max is a no-op. Left in place rather than halved: it costs one extra
    // correlate per candidate pair, only when the bank is at its cap, and
    // removing it would silently change the merge choice if the metric ever
    // becomes directional again.
    inline double groupCloseness(const BeatGroup& x, const BeatGroup& y)
    {
        double sum = 0.0; int n = 0;
        for (int c = 0; c < kNumChannels; ++c) {
            if (x.ch[c].tmpl.empty() || y.ch[c].tmpl.empty()) continue;
            // ONE DIRECTION. The old score was the fraction of one waveform
            // inside the other's corridor, which is directional -- a wide group
            // swallowed a narrow one asymmetrically -- so both were taken and
            // the larger used. The score is a correlation now, and correlation
            // is symmetric, so the second call returned the same number and
            // doubled the cost of every candidate pair. findMergePair runs on
            // EVERY spawn once the bank is at cap, which on a thrashing bin is
            // hundreds of times.
            const double s = tbank::bandMatch(x.ch[c].tmpl, y.ch[c]).score;
            if (std::isnan(s)) continue;
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

        // ---- ONLY THE SURVIVOR IS REBUILT, EVEN WHEN IT IS SLOT 0 --------
        //
        // This used to call recomputeAll when a == 0, because slot 0's corridor
        // is what every young group inherits, so moving it moved theirs. That
        // reasoning died with the metric: bandMatch scores a CORRELATION
        // against tmpl now, and band_lo/band_hi decide nothing. The inherited
        // corridor is a display and QC quantity, and rebuilding six groups on
        // four channels to keep it exact during a merge storm is not worth what
        // it costs.
        //
        // WHAT IT COST. A record with 592 merges per bin, most of them landing
        // on slot 0: each one re-medianed slot 0's ~600 members over a
        // 2344-column axis on four channels, then every other group. 29 seconds
        // per bin, 48 bins.
        //
        // The corridors are exact again at the end of the pass: runBank calls
        // recomputeAll once when the assignment loop finishes, and cleanGroups
        // calls it again after excluding. Nothing reads a corridor in between.
        recomputeOne(bank, chans, a);
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

        // ---- SECTION 4.5 PER-BIN CATEGORY CENSUS -------------------------
        // Counts, never percentages: a percentage cannot be re-aggregated
        // across bins and cannot say how many beats it rested on. The writer
        // divides.
        //
        // These were filled by bin_pipeline::runChannel, per channel, over a
        // per-channel partition. A category is a property of the beat, so
        // counting it per channel gave three answers to a one-answer question.
        // Counted here once, over slices.
        uint32_t n_regular = 0;
        uint32_t n_ectopic = 0;
        uint32_t n_noise = 0;
    };

    // Distinct CONFIRMED templates carrying one class label. One means
    // monomorphic, two or more polymorphic.
    //
    // Counted over GROUPS, with no max across channels. tbank::polymorphicVerdict
    // takes three banks and reports the largest per-channel count, because there
    // used to be three partitions that could disagree about how many
    // morphologies a bin held. There is one partition now, so the count is just
    // a count, and "which channel drove it" is not a question that exists.
    //
    // ONLY CONFIRMED GROUPS COUNT, so the algorithm can propose a second
    // morphology but cannot declare polymorphy on its own authority -- it stays
    // an operator-gated finding.
    struct PolymorphyVerdict {
        int count = 0;
        int n_unconfirmed_groups = 0;   // proposed, not yet ruled on
        bool polymorphic() const { return count >= 2; }
        bool monomorphic() const { return count == 1; }
    };

    inline PolymorphyVerdict polymorphyVerdict(const JointBank& bank,
        uint8_t code = tbank::kCodePvc)
    {
        PolymorphyVerdict v;
        for (const BeatGroup& g : bank.groups) {
            const tbank::BankTemplate& t = g.ch[kCh1];
            if (!t.confirmed()) {
                if (g.memberCount() > 0) ++v.n_unconfirmed_groups;
                continue;
            }
            if (t.label_code == code) ++v.count;
        }
        return v;
    }

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
        if (bank.next_spawn_seq > 200) { out.group_id = tbank::kUnscorable; return out; }
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

        // Rebuilt so the new group inherits slot 0's corridor on every
        // channel. A one-member group has no spread of its own -- lo == hi at
        // every column -- so without the inherited corridor the next beat of
        // the same morphology scores ~0 and spawns again.
        //
        // ONLY the new group: a spawn does not touch slot 0's membership, so
        // slot 0's corridor is the same corridor it was a moment ago and every
        // other group's inheritance from it is unchanged.
        recomputeOne(bank, chans, bank.size() - 1);
        return out;
    }

    // ---------------------------------------------------------------------
    // Seed + drive
    // ---------------------------------------------------------------------

    // Slot 0 from Phase 1. One group, seeded per channel from that channel's
    // Phase 1 template: the ECG leads from chN's sinus template, PPG from the
    // bin's pulse template. `seed_slices` are the slices Phase 1 built those
    // templates from, and they become slot 0's initial membership.
    // How wide slot 0's corridor is when Phase 1 supplies the spread, in
    // standard deviations either side of the Phase 1 waveform.
    //
    // 1.96 because recomputeTemplate builds its own corridors from the 2.5 and
    // 97.5 percentiles, and +/- 1.96 sigma is the same interval for a normal
    // column -- so slot 0's band means the same thing as every other group's
    // band rather than being a differently-defined object that happens to sit
    // in the same field.
    //
    // THIS NUMBER SETS HOW READILY EVERY BIN SPAWNS. bandMatch scores the
    // fraction of columns inside the corridor against 0.85, so a wider band
    // admits more morphologies into slot 0 and a narrower one splits sooner.
    // It is the one free parameter in the seeding and it is named here so it can
    // be changed in one place.
    inline constexpr double kPhase1BandSigma = 1.96;

    inline void seedBank(JointBank& bank, const ChannelSet& chans,
        const std::array<std::vector<double>, kNumChannels>& phase1,
        const std::vector<uint32_t>& seed_slices,
        int32_t max_templates_per_bin = 0,
        // Per-sample spread of the beats Phase 1's waveform was built from, per
        // channel. Empty entries leave that channel's corridor as recomputed.
        const std::array<std::vector<double>, kNumChannels>* phase1_spread = nullptr)
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
        // use ONE reference.
        for (int c = 0; c < kNumChannels; ++c)
            if (!phase1[c].empty())
                bank.groups[0].ch[c].tmpl = phase1[c];

        // ---- AND THE CORRIDOR, WHICH IS THE PART THAT ACTUALLY DECIDES -----
        //
        // It used to stay as recomputed, on the reasoning that Phase 1 supplies
        // a waveform and not a spread. But bandMatch scores against the
        // CORRIDOR, not the waveform: a beat is assigned by what fraction of its
        // columns fall between band_lo and band_hi. The recomputed corridor
        // comes from seed_slices, which is the first twenty slices present on
        // any channel -- and on bigeminy that is about ten sinus beats and ten
        // ectopic ones, so the band spans the gap BETWEEN the two morphologies
        // and admits both. A sinus waveform on slot 0 does not help when the
        // band drawn around it is wide enough to swallow a PVC. The bank
        // converged to one template on a record that has two, which is the
        // failure Section 4.6 exists to prevent.
        //
        // So Phase 1 supplies both, from the same masked pool: the median is the
        // waveform, the per-sample spread is the band. No rhythm flag touches
        // the partition -- the mask was applied upstream, when the reference was
        // built -- so the required order still holds: partition first, then
        // remove premature, then Tukey on what is left.
        if (phase1_spread)
            for (int c = 0; c < kNumChannels; ++c) {
                const std::vector<double>& sd = (*phase1_spread)[c];
                std::vector<double>& tm = bank.groups[0].ch[c].tmpl;
                if (sd.empty() || tm.empty()) continue;
                tbank::BankTemplate& t = bank.groups[0].ch[c];
                t.band_lo.assign(tm.size(),
                    std::numeric_limits<double>::quiet_NaN());
                t.band_hi = t.band_lo;
                for (size_t j = 0; j < tm.size(); ++j) {
                    if (j >= sd.size() || std::isnan(tm[j])) continue;
                    // A zero-spread column is not a zero-width corridor: it
                    // means the reference beats agreed exactly there, usually
                    // because only one contributed. Left NaN, which bandMatch
                    // treats as an unscorable column and excludes from its
                    // denominator, rather than as a band nothing can fall
                    // inside.
                    if (!(sd[j] > 0.0)) continue;
                    t.band_lo[j] = tm[j] - kPhase1BandSigma * sd[j];
                    t.band_hi[j] = tm[j] + kPhase1BandSigma * sd[j];
                }
                t.corridor_inherited = false;
            }

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
    // POST-PARTITION: REMOVE AND FLAG PREMATURE, THEN TUKEY ON WHAT IS LEFT
    // ---------------------------------------------------------------------
    //
    // ORDER. align -> partition and merge -> remove and flag premature -> Tukey
    // on the clean beats. This function is the last two steps, and it runs
    // AFTER runBank() because both are judgements about a beat relative to ITS
    // OWN morphology, which does not exist until the partition does.
    //
    // Tukey used to run inside align_beat_matrix, before any partition. Its
    // amplitude and wave-score passes reject beats that do not resemble the
    // population, and before partitioning the population is every morphology at
    // once -- so an ectopic beat is an outlier by construction and was deleted
    // before the bank could separate it. alignment.hpp now measures and records
    // without pruning; this is where the excluding happens.
    //
    // EXCLUDED FROM THE TEMPLATE, RETAINED WITH FLAGS. An excluded beat keeps
    // its membership in `members` and gains a reason in `excluded_reason`. It is
    // dropped only from `members_clean`, which is what the waveform averages
    // over. Deleting it would make the count unreconstructable; leaving it in
    // the average is the contamination Section 4.6 exists to prevent.

    // Why a beat is not in its group's average. One value per slice, so the
    // archive can state the reason rather than only the fact.
    enum class ExcludeReason : uint8_t {
        KEPT = 0,
        NOT_A_MEMBER = 1,   // no group claimed the slice at all
        PREMATURE = 2,   // the prematurity test fired on this beat
        VOTE = 3,   // the 5-of-8 vote fired, the beat itself did not
        TUKEY_R_LOCATION = 4,
        TUKEY_AMPLITUDE = 5,
        TUKEY_WAVE_SCORE = 6,
        // The operator marked this beat ECTOPIC or NOISE. Section 4.5: only
        // category 1 feeds averaging. Distinct from PREMATURE because the
        // evidence is different in kind -- a mark is a human verdict on the
        // morphology, prematurity is an arithmetic verdict on an interval, and
        // they disagree in both directions (a PAC is marked and not premature;
        // a sinus beat after a pause is premature and unmarked).
        CATEGORY = 7
    };

    inline const char* excludeReasonName(uint8_t r) {
        switch (static_cast<ExcludeReason>(r)) {
        case ExcludeReason::KEPT:             return "kept";
        case ExcludeReason::NOT_A_MEMBER:     return "not_a_member";
        case ExcludeReason::PREMATURE:        return "premature";
        case ExcludeReason::VOTE:             return "vote";
        case ExcludeReason::TUKEY_R_LOCATION: return "tukey_r_location";
        case ExcludeReason::TUKEY_AMPLITUDE:  return "tukey_amplitude";
        case ExcludeReason::TUKEY_WAVE_SCORE: return "tukey_wave_score";
        case ExcludeReason::CATEGORY:         return "category";
        }
        return "?";
    }

    // ExcludeReason -> the per-beat Tukey verdict the archive already carries.
    // BeatFlags::tukey is tbank::TukeyOutcome, a separate enumeration from the
    // reason codes here, and the mapping is written out rather than cast: the
    // two lists share no values and a static_cast would compile into nonsense.
    inline tbank::TukeyOutcome tukeyOutcomeFor(uint8_t reason) {
        switch (static_cast<ExcludeReason>(reason)) {
        case ExcludeReason::TUKEY_R_LOCATION: return tbank::TukeyOutcome::REJ_R_LOCATION;
        case ExcludeReason::TUKEY_AMPLITUDE:  return tbank::TukeyOutcome::REJ_AMPLITUDE;
        case ExcludeReason::TUKEY_WAVE_SCORE: return tbank::TukeyOutcome::REJ_WAVE_SCORE;
        default: break;
        }
        return tbank::TukeyOutcome::KEPT;
    }

    struct CleanCounts {
        uint32_t members_total = 0;
        uint32_t excluded_premature = 0;
        uint32_t excluded_vote = 0;
        uint32_t excluded_tukey = 0;
        uint32_t excluded_category = 0;
        uint32_t kept = 0;
        // Groups whose every member was premature, so nothing was removed. That
        // is the ectopic morphology itself and it must keep its waveform;
        // counted because a bin where this fires on the SEED is a bin whose
        // reference is ectopic.
        uint32_t groups_all_premature = 0;
        // Groups where every eligible member failed a Tukey fence. Also left
        // intact, for the same reason, and also worth seeing.
        uint32_t groups_all_rejected = 0;
    };

    // Metrics for one member on one channel. Amplitude and R-location are
    // measured on the beat; wave-score is the beat against ITS OWN group's
    // template, which is the pass that gains most from running after the
    // partition instead of before it.
    namespace detail {

        struct MemberMetrics {
            std::vector<double> amp;    // max - min over non-NaN samples
            std::vector<double> rloc;   // argmax|v| - anchor column
            std::vector<double> wave;   // correlation against the group template
        };

        inline MemberMetrics measureOnChannel(const BeatGroup& g, int c,
            const ChannelSet& chans,
            const std::vector<uint32_t>& cand)
        {
            const double NaN = std::numeric_limits<double>::quiet_NaN();
            MemberMetrics m;
            m.amp.assign(cand.size(), NaN);
            m.rloc.assign(cand.size(), NaN);
            m.wave.assign(cand.size(), NaN);

            const ChannelBeats& cb = chans[c];
            const tbank::BankTemplate& tp = g.ch[c];
            if (!cb.present()) return m;

            for (size_t k = 0; k < cand.size(); ++k) {
                const std::vector<double>* bp = cb.beatFor(cand[k]);
                if (!bp) continue;                 // channel dropped this slice
                const std::vector<double>& b = *bp;

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
                if (hi >= lo) m.amp[k] = hi - lo;

                // Offset of this beat's own dominant deflection from the
                // frame's anchor column. A beat whose R (or systolic peak) did
                // not land where the frame says it did corrupts every shape
                // measurement made on it, which is why this is measured rather
                // than assumed.
                if (argmax >= 0 && cb.anchor_col >= 0)
                    m.rloc[k] = static_cast<double>(argmax - cb.anchor_col);

                if (!tp.tmpl.empty())
                    m.wave[k] = tbank::correlate(b, tp.tmpl).r;
            }
            return m;
        }

    }  // namespace detail

    // `flags` is indexed by SLICE and must already carry the prematurity
    // verdict (BeatFlags::pvc). This function reads that and writes
    // BeatFlags::tukey, so after it returns the flag vector is the complete
    // per-beat record: category from the operator, prematurity from timing,
    // Tukey from morphology inside the beat's own group.
    inline void cleanGroups(JointBank& bank, const ChannelSet& chans,
        std::vector<tbank::BeatFlags>& flags,
        std::vector<uint8_t>& excluded_reason,
        CleanCounts* counts = nullptr)
    {
        const uint32_t n_slices = static_cast<uint32_t>(flags.size());
        excluded_reason.assign(n_slices,
            static_cast<uint8_t>(ExcludeReason::NOT_A_MEMBER));

        for (BeatGroup& g : bank.groups) {
            if (counts) counts->members_total +=
                static_cast<uint32_t>(g.members.size());

            // ---- 1: premature removal, per group -------------------------
            // Counted per group and folded into `counts` only if this group is
            // actually cleaned. The restore path below has to be able to undo
            // its own accounting, and it cannot do that against a running total
            // that other groups have already contributed to -- an earlier
            // version zeroed the record-wide counters whenever any one group
            // turned out to be entirely premature, so a single all-ectopic
            // group erased every exclusion count in the bin.
            uint32_t ex_cat = 0, ex_prem = 0, ex_vote = 0;
            std::vector<uint32_t> clean;
            clean.reserve(g.members.size());
            for (const uint32_t m : g.members) {
                const tbank::Category cat = (m < flags.size())
                    ? flags[m].category : tbank::Category::REGULAR;
                // CATEGORY BEFORE PREMATURITY, so the reason reported is the
                // operator's when both apply. A marked PVC is usually also
                // premature, and "the operator called this ectopic" is the
                // stronger and more auditable statement of the two.
                if (cat != tbank::Category::REGULAR) {
                    excluded_reason[m] =
                        static_cast<uint8_t>(ExcludeReason::CATEGORY);
                    ++ex_cat;
                    continue;
                }
                const tbank::PvcFilter v = (m < flags.size())
                    ? flags[m].pvc : tbank::PvcFilter::NONE;
                if (v == tbank::PvcFilter::NONE) { clean.push_back(m); continue; }
                excluded_reason[m] = static_cast<uint8_t>(
                    v == tbank::PvcFilter::PREMATURE ? ExcludeReason::PREMATURE
                    : ExcludeReason::VOTE);
                if (v == tbank::PvcFilter::PREMATURE) ++ex_prem; else ++ex_vote;
            }

            // A GROUP THAT IS ENTIRELY PREMATURE IS NOT CLEANED. That is the
            // ectopic morphology itself, and emptying it would delete exactly
            // what 4.6 exists to preserve -- the template would lose its
            // waveform and the operator would have nothing to confirm.
            if (clean.empty()) {
                clean = g.members;
                for (const uint32_t m : clean)
                    excluded_reason[m] = static_cast<uint8_t>(ExcludeReason::KEPT);
                ex_cat = ex_prem = ex_vote = 0;   // this group only
                if (counts) ++counts->groups_all_premature;
            }
            if (counts) {
                counts->excluded_category += ex_cat;
                counts->excluded_premature += ex_prem;
                counts->excluded_vote += ex_vote;
            }

            // ---- 2: Tukey INSIDE the group, THREE PASSES, ALL CHANNELS ---
            //
            // R-LOCATION, AMPLITUDE, WAVE-SCORE, measured on the SAME
            // post-premature set and rejected as a UNION, not sequentially.
            // alignment.hpp ran them in series, each pass computing its
            // quartiles over the previous pass's survivors, which is far more
            // aggressive: three sequential 1.5*IQR passes on a 12-member group
            // can cascade to almost nothing, and small groups are exactly what
            // this pipeline produces. Simultaneous fences are also
            // order-independent, so nobody has to know which pass ran first to
            // interpret the result.
            //
            // RR-LENGTH IS NOT ONE OF THEM. It is a rhythm test, not a
            // morphology test, and premature removal above has already taken
            // the short intervals. What would be left for it to reject is
            // mostly long ones -- post-ectopic pauses and dropped detections --
            // and inside a group of ectopics it would compute ectopic RR fences
            // and reject on those, which is a different statement again.
            //
            // FENCES ARE PER CHANNEL, THE VERDICT IS NOT. Amplitude on CH1 and
            // amplitude on PPG are different quantities in different units, so
            // each channel gets its own quartiles. A failure on ANY channel
            // excludes the beat from the group, because the group is one set of
            // beats -- excluding per channel would give the three leads and the
            // pulse four different averages of one morphology.
            std::vector<uint32_t> kept = clean;
            if (clean.size() >= 8) {
                std::vector<uint8_t> reject(clean.size(), 0);   // ExcludeReason

                for (int c = 0; c < kNumChannels; ++c) {
                    if (!chans[c].present()) continue;
                    const detail::MemberMetrics mm =
                        detail::measureOnChannel(g, c, chans, clean);

                    alignment::TukeyStats sA, sR, sW;
                    const std::vector<bool> keepA =
                        alignment::keep_within_tukey(mm.amp, 1.5, &sA);
                    const std::vector<bool> keepR =
                        alignment::keep_within_tukey(mm.rloc, 1.5, &sR);
                    const std::vector<bool> keepW =
                        alignment::keep_within_tukey(mm.wave, 1.5, &sW);

                    for (size_t k = 0; k < clean.size(); ++k) {
                        // A metric that could not be measured does NOT reject.
                        // NaN means unjudgeable, and rejecting on it would prune
                        // beats for being unmeasurable rather than for being
                        // outliers -- the same conflation bandMatch avoids by
                        // excluding incomparable columns from its denominator.
                        // This is also what makes a channel that dropped the
                        // beat abstain instead of voting against it.
                        const bool okA = std::isnan(mm.amp[k])
                            || (k < keepA.size() && keepA[k]);
                        const bool okR = std::isnan(mm.rloc[k])
                            || (k < keepR.size() && keepR[k]);
                        const bool okW = std::isnan(mm.wave[k])
                            || (k < keepW.size() && keepW[k]);
                        if (okA && okR && okW) continue;
                        if (reject[k]) continue;   // first channel to fail owns it

                        // FIRST failing pass wins the reason, ordered
                        // R-location, amplitude, wave-score: an R that landed
                        // wrong explains a bad amplitude and a bad shape, so it
                        // is the more informative attribution when several fire.
                        reject[k] = static_cast<uint8_t>(!okR
                            ? ExcludeReason::TUKEY_R_LOCATION
                            : (!okA ? ExcludeReason::TUKEY_AMPLITUDE
                                : ExcludeReason::TUKEY_WAVE_SCORE));
                    }
                }

                std::vector<uint32_t> k2;
                k2.reserve(clean.size());
                for (size_t k = 0; k < clean.size(); ++k)
                    if (!reject[k]) k2.push_back(clean[k]);

                if (!k2.empty()) {
                    for (size_t k = 0; k < clean.size(); ++k)
                        if (reject[k]) {
                            excluded_reason[clean[k]] = reject[k];
                            if (clean[k] < flags.size())
                                flags[clean[k]].tukey = tukeyOutcomeFor(reject[k]);
                            if (counts) ++counts->excluded_tukey;
                        }
                    kept = std::move(k2);
                }
                else if (counts) {
                    // Every eligible member failed. The group keeps them all
                    // rather than losing its waveform, and nothing is flagged,
                    // because "excluded" would then describe beats that are
                    // still in the average.
                    ++counts->groups_all_rejected;
                }
            }

            for (const uint32_t m : kept)
                excluded_reason[m] = static_cast<uint8_t>(ExcludeReason::KEPT);
            if (counts) counts->kept += static_cast<uint32_t>(kept.size());

            g.members_clean = std::move(kept);
        }

        // Corridors depend on membership, so this happens once after every
        // group is cleaned rather than per group: slot 0's corridor is what the
        // young groups inherit and it has to be final first.
        recomputeAll(bank, chans);
    }

    // ---------------------------------------------------------------------
    // SECTION 4.6 BEAT SUBSTITUTION
    // ---------------------------------------------------------------------
    //
    // A borderline beat is replaced by an EWMA blend of its group's running
    // average and itself, alpha = 1/8, so that beat-to-beat series have no
    // holes. A hole is worse than a smoothed estimate for anything measuring
    // variability: an RR or feature series with gaps reports a variance that
    // depends on which beats were dropped.
    //
    // WHICH SCORE DECIDES "BORDERLINE", and why it is not the assignment score.
    // beat_substitute's band is 0.60 to 0.85 and 0.85 is the assignment floor,
    // so if the score were the one the bank assigned on, THE BAND WOULD BE
    // EMPTY BY CONSTRUCTION -- every assigned beat cleared 0.85 and everything
    // below it spawned its own group instead. The two 0.85s are different
    // quantities: bandMatch's score is the FRACTION OF COLUMNS inside the
    // corridor, while the spec's "morphology correlation threshold r < 0.85" is
    // a CORRELATION. So this uses tbank::correlate against the beat's own
    // group template, which is what the spec's sentence says, and which has a
    // real spread among beats the corridor accepted.
    //
    // RUNS AFTER cleanGroups, so a beat already excluded for prematurity, for a
    // category or by a Tukey fence is not also substituted -- substituting a
    // beat that is not in the average is work with no consumer, and flagging it
    // twice would double-count it in any burden statistic.
    //
    // NOTHING IS OVERWRITTEN. The blend is stored beside the beat, never in
    // place of it: a substituted beat is not an observation. members_clean is
    // untouched, so the group's median stays a median of observations -- feeding
    // blends back into it would make the median partly a function of its own
    // history, which is the mean/median mismatch beat_substitute.hpp warns
    // about.
    struct Substitution {
        uint32_t slice = 0;
        uint8_t  channel = 0;
        double   score = 0.0;          // correlation that put it in the band
        std::vector<double> blended;
    };

    struct SubstitutionCounts {
        uint32_t n_substituted = 0;    // beats, not beat-channels
        uint32_t n_channel_blends = 0;
        uint32_t n_too_bad = 0;        // below the band: rejected, not blended
    };

    inline void substituteBorderline(const JointBank& bank,
        const ChannelSet& chans,
        const std::vector<uint8_t>& excluded_reason,
        std::vector<tbank::BeatFlags>& flags,
        std::vector<Substitution>& out,
        SubstitutionCounts* counts = nullptr)
    {
        for (const BeatGroup& g : bank.groups) {
            // Running average per channel, seeded from the group's template.
            // Seeding from the median rather than from the first member means
            // the first borderline beat is blended against the whole
            // population, not against whichever beat happened to arrive first.
            std::array<std::vector<double>, kNumChannels> avg;
            for (int c = 0; c < kNumChannels; ++c) avg[c] = g.ch[c].tmpl;

            // IN SLICE ORDER, because this is a recursion in time: the average
            // is updated with the substituted value, so consecutive borderline
            // beats drift toward it rather than each being pulled the same
            // distance from it. members is sorted, so iterating it is time
            // order.
            for (const uint32_t slice : g.members) {
                if (slice < excluded_reason.size()
                    && excluded_reason[slice]
                    != static_cast<uint8_t>(ExcludeReason::KEPT)) continue;

                bool blended_any = false;
                for (int c = 0; c < kNumChannels; ++c) {
                    const std::vector<double>* beat = chans[c].beatFor(slice);
                    if (!beat || avg[c].empty()) continue;
                    const tbank::CorrResult cr =
                        tbank::correlate(*beat, g.ch[c].tmpl);
                    if (!cr.scorable()) continue;

                    if (cr.r >= beat_substitute::borderlineHi()) continue;
                    if (cr.r < beat_substitute::kBorderlineLo) {
                        // Below the floor a beat is not borderline, it is bad,
                        // and blending it would smuggle an artifact into the
                        // series wearing the average's shape.
                        if (counts) ++counts->n_too_bad;
                        continue;
                    }

                    Substitution sub;
                    sub.slice = slice;
                    sub.channel = static_cast<uint8_t>(c);
                    sub.score = cr.r;
                    sub.blended = beat_substitute::substituteBeatNaNSafe(
                        avg[c], *beat, beat_substitute::kAlpha);
                    avg[c] = sub.blended;      // recursion, not a fixed offset
                    out.push_back(std::move(sub));
                    blended_any = true;
                    if (counts) ++counts->n_channel_blends;
                }

                // FLAGGED ON THE BEAT, not per channel. BeatFlags is per slice
                // and the flag's job is to let downstream variance, corridor
                // and IQR calculations exclude synthetic values -- a beat
                // blended on one lead is not a clean observation on that lead,
                // and one flag per beat is what the archive can carry.
                if (blended_any) {
                    if (slice < flags.size()) flags[slice].substituted = true;
                    if (counts) ++counts->n_substituted;
                }
            }
        }
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
        // Per-sample spread of the beats ecg_phase1 was built from
        // (ecgTemplates_raw_iqr, which is a std with ddof=1 despite the name).
        // Supplies slot 0's corridor; see seedBank.
        std::vector<double> ecg_phase1_spread[kNumEcgCh];

        const std::vector<std::vector<double>>* ppg_beats = nullptr;
        const std::vector<uint32_t>* ppg_forward = nullptr;
        int ppg_peak_col = -1;
        std::vector<double> ppg_phase1;
        std::vector<double> ppg_phase1_spread;

        uint32_t n_slices = 0;
        int32_t  max_templates_per_bin = 0;   // 0 => kDefaultMaxTemplatesPerBin
        uint64_t bin_index = 0;

        // ---- PER-SLICE INPUTS FOR THE POST-PARTITION STAGE ---------------
        //
        // rr_after_ms[s] = R[s+1] - R[s] in milliseconds, one entry per slice.
        // Taken straight off the R-peak vector that defines the slices, so it
        // needs no index map and cannot desynchronise from them -- which is
        // exactly what went wrong when the prematurity test was fed a
        // per-channel aligned RR series and asked about a joint partition.
        //
        // Empty is legal and means no prematurity verdict is available: every
        // beat then reads NONE and nothing is excluded for rhythm. Honest, and
        // it degrades a bin rather than the record.
        std::vector<double> rr_after_ms;

        // Operator marks per slice, 0 = unmarked. Normally empty at build time:
        // marking happens after templates exist, so morphology does the sorting
        // and marks only ever supply LABELS later.
        std::vector<uint8_t> mark_code;
    };

    struct BinBankOutput {
        JointBank bank;
        std::vector<int32_t> group_of_slice;   // sized n_slices; -1/kNoMatch/kUnscorable
        BankCounts counts;
        std::vector<tbank::CapRaiseEvent> cap_raises;

        // ---- PER-SLICE, AND SHARED BY ALL FOUR CHANNELS ------------------
        //
        // One verdict per heartbeat, not one per channel-row. category comes
        // from the operator's mark, pvc from the timing test, tukey from the
        // fences inside the beat's own group. Sized n_slices, so it is indexed
        // identically to group_of_slice and to excluded_reason -- the three of
        // them together are the per-beat record the archive writes.
        std::vector<tbank::BeatFlags> flags;

        // jbank::ExcludeReason per slice. KEPT means the beat is in its group's
        // average; NOT_A_MEMBER means no group claimed the slice at all, which
        // is a different statement from being excluded by a test.
        std::vector<uint8_t> excluded_reason;

        CleanCounts clean;
        pvc_filter::FilterResult pvc;

        // 4.6 substitutions: the blends, and what they cost. Stored, not
        // applied -- see substituteBorderline.
        std::vector<Substitution> substitutions;
        SubstitutionCounts subs;

        // The per-slice RR series this bin was scored with, carried through so
        // the record-level NSVT pass uses the SAME intervals the prematurity
        // filter did. Recomputing them downstream from the R-peaks would be a
        // second derivation of one quantity, and the two would differ first on
        // the bins with dropout gaps -- exactly the bins where a run matters.
        std::vector<double> rr_after_ms;
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
        std::array<std::vector<double>, kNumChannels> spread;
        for (int c = 0; c < kNumEcgCh; ++c) spread[c] = in.ecg_phase1_spread[c];
        spread[kPpg] = in.ppg_phase1_spread;

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

        seedBank(out.bank, chans, phase1, seedSlices, in.max_templates_per_bin,
            &spread);
        runBank(out.bank, chans, in.n_slices, out.group_of_slice,
            in.bin_index, &out.cap_raises, &out.counts);

        // ---- PER-SLICE FLAGS --------------------------------------------
        //
        // Built AFTER the partition and consumed only by the stage below, which
        // is the required order: the bank is rhythm-blind and Tukey-blind, and
        // both verdicts act only once each beat has a morphology to be judged
        // against. Computing them here rather than before runBank() is what
        // makes that structural instead of a convention.
        out.flags.assign(in.n_slices, tbank::BeatFlags{});
        // categoriesFromMarks, not categoryForLabelCode per slice: the beat
        // AFTER an ectopic one inherits ECTOPIC when it carries no mark of its
        // own, and that rule needs the neighbouring mark, so it cannot be
        // applied one beat at a time.
        const std::vector<tbank::Category> cats =
            in.mark_code.empty() ? std::vector<tbank::Category>{}
                                 : tbank::categoriesFromMarks(in.mark_code);
        for (uint32_t s = 0; s < in.n_slices; ++s) {
            out.flags[s].category = (s < cats.size())
                ? cats[s] : tbank::Category::REGULAR;
            const int32_t g = (s < out.group_of_slice.size())
                ? out.group_of_slice[s] : tbank::kNoMatch;
            for (int c = 0; c < kNumEcgCh; ++c)
                out.flags[s].template_id_ecg[c] = g;
            out.flags[s].template_id_ppg = g;
        }
        if (!in.rr_after_ms.empty()) {
            out.pvc = pvc_filter::runFilter(in.rr_after_ms);
            for (uint32_t s = 0; s < in.n_slices
                && s < out.pvc.verdict.size(); ++s)
                out.flags[s].pvc = out.pvc.verdict[s];
        }

        // The category census, over the flags just built. Before cleanGroups so
        // it describes the population the partition saw, not the survivors.
        for (uint32_t s = 0; s < in.n_slices; ++s) {
            switch (out.flags[s].category) {
            case tbank::Category::ECTOPIC: ++out.counts.n_ectopic; break;
            case tbank::Category::NOISE:   ++out.counts.n_noise;   break;
            default:                       ++out.counts.n_regular; break;
            }
        }

        cleanGroups(out.bank, chans, out.flags, out.excluded_reason, &out.clean);

        substituteBorderline(out.bank, chans, out.excluded_reason, out.flags,
            out.substitutions, &out.subs);
        out.rr_after_ms = in.rr_after_ms;
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
        const ChannelSet& chans, int channel,
        // Per-SLICE flags, when available. Used only to fill the per-template
        // census (n_premature_members and friends), which templates.csv reports
        // beside the class row and presumedCategory() reads. Left null the
        // census stays zero, and every template presents as non-ectopic -- so
        // this is not optional in practice, only in signature.
        const std::vector<tbank::BeatFlags>* flags = nullptr)
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
            // BOTH LISTS ARE TRANSLATED. `members` is everything the partition
            // assigned, which is what the archive writes; `members_clean` is
            // what the waveform in t.tmpl was actually averaged over. Carrying
            // only one of them is how the two descriptions of a template come
            // apart: t.memberCount() would then either overstate what is behind
            // the trace or lose the excluded beats from the record entirely.
            t.members.clear();
            for (const uint32_t slice : g.members) {
                const int k = cb.localFor(slice);
                if (k >= 0) t.members.push_back(static_cast<uint32_t>(k));
            }
            t.members_clean.clear();
            if (!g.members_clean.empty())
                for (const uint32_t slice : g.members_clean) {
                    const int k = cb.localFor(slice);
                    if (k >= 0) t.members_clean.push_back(static_cast<uint32_t>(k));
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

            // ---- census, over the FULL membership --------------------------
            // Deliberately the full list: "this template holds mostly premature
            // beats" is a statement about what landed in the morphology, and
            // counting only the clean members would make it read zero by
            // construction, since premature beats are exactly what cleaning
            // removed. 4.6 never reassigns a category, it only excludes, so the
            // premature row and the class row are reported side by side.
            t.n_premature_members = 0;
            t.n_voted_members = 0;
            t.n_noise_members = 0;
            if (flags) {
                for (const uint32_t slice : g.members) {
                    if (slice >= flags->size()) continue;
                    const tbank::BeatFlags& bf = (*flags)[slice];
                    switch (bf.pvc) {
                    case tbank::PvcFilter::PREMATURE: ++t.n_premature_members; break;
                    case tbank::PvcFilter::VOTE:      ++t.n_voted_members;     break;
                    default: break;
                    }
                    if (bf.category == tbank::Category::NOISE) ++t.n_noise_members;
                    if (bf.tukey != tbank::TukeyOutcome::NOT_ELIGIBLE
                        && bf.tukey != tbank::TukeyOutcome::KEPT)
                        ++t.n_tukey_members;
                }
            }
            out.templates.push_back(std::move(t));
        }
        out.assigned_beats = bank.assigned_beats;
        return out;
    }

    // All four at once, rebuilding the ChannelSet the bank was run against.
    // The caller must pass the SAME ChannelSet, or the local index spaces will
    // not match the ones the members were resolved from.
    inline std::array<tbank::TemplateBank, kNumChannels> projectAll(
        const JointBank& bank, const ChannelSet& chans,
        const std::vector<tbank::BeatFlags>* flags = nullptr)
    {
        std::array<tbank::TemplateBank, kNumChannels> out;
        for (int c = 0; c < kNumChannels; ++c)
            out[c] = projectToChannel(bank, chans, c, flags);
        return out;
    }

}  // namespace jbank
