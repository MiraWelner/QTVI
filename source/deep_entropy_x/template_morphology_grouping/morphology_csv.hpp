#pragma once
/**
 * @file   morphology_csv.hpp
 * @brief  The descriptions of the outputs, per tempalte and per bin. beats.csv was too large so its just beats.bin
 *
 *           <stem>_bins.csv       one row per BIN          (text)
 *           <stem>_nsvt.csv       one row per RUN          (text)
 *           <stem>_templates.csv  one column per TEMPLATE  (text)
 *           <stem>_beats.bin      one column per SLICE     (binary)
 *           <stem>_templates.bin  one column per TEMPLATE  (binary)
 *
 *         _templates.bin is the only one of these that is READ BACK: it is
 *         what carries a morphology split forward, so a config or algorithm
 *         change does not repartition a record that has already been
 *         partitioned. See bank_reload.hpp. Everything else is write-only.
 */

#include <algorithm>   // std::max, for the per-bin run length
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "bin_pipeline.hpp"

namespace morphology_csv {

    inline std::string g_dir;
    inline std::string g_stem;

    inline void set(const std::string& dir, const std::string& stem,
        const std::string& subdir = "") {
        g_stem = stem;
        if (dir.empty()) { g_dir.clear(); return; }
        std::filesystem::path p(dir);
        if (!subdir.empty()) p /= subdir;
        std::error_code ec;
        std::filesystem::create_directories(p, ec);
        // ofstream on a missing directory fails silently and the writer would
        // return having produced nothing, so the destination is created here
        // rather than assumed.
        g_dir = ec ? dir : p.string();
    }

    namespace detail {

        // EVERY enumerator named explicitly, no catch-all. `default: return
        // "pqrst"` would map any future enumerator to the word for NORMAL,
        // silently, so a new category would enter the archive as sinus without
        // anyone deciding that. A missing case is now a compiler warning instead
        // of a plausible-looking column.
        inline const char* categoryWord(tbank::Category c) {
            switch (c) {
            case tbank::Category::REGULAR:     return "pqrst";
            case tbank::Category::ECTOPIC:     return "ectopic";
            case tbank::Category::NOISE:       return "noise";
            }
            return "pqrst";   // unreachable; every enumerator is named above
        }

        inline std::vector<uint8_t> letterRanks(const tbank::TemplateBank& bank) {
            return tbank::letterRanks(bank);
        }

        inline std::string templateName(const tbank::BankTemplate& t,
            int letter_rank = -1) {
            std::string cls = "PQRST";
            int letterIdx = -1;

            if (t.label_code != tbank::kUnlabeled) {
                switch (t.label_code) {
                case tbank::kCodePvc:        cls = "PVC";   break;
                case tbank::kCodePac:        cls = "PAC";   break;
                case tbank::kCodeVt:         cls = "VT";    break;
                case tbank::kCodeMinorNoise: cls = "NOISE"; break;
                default: cls = "CODE" + std::to_string(t.label_code); break;
                }
                if (t.subtype > 0) letterIdx = t.subtype - 1;
            }
            // Ranked letter when the caller supplied one; spawn_seq only as a
            // last resort. spawn_seq counts every template ever created in the
            // bin, so on a bank that merged 460 times it names surviving
            // templates A, C, E, ... with the gaps being the merged ones.
            if (letterIdx < 0)
                letterIdx = (letter_rank >= 0) ? letter_rank
                : static_cast<int>(t.spawn_seq);
            const char letter = static_cast<char>('A' + (letterIdx % 26));
            return cls + "_" + std::string(1, letter);
        }

        inline const char* prematureWord(tbank::PvcFilter p) {
            switch (p) {
            case tbank::PvcFilter::PREMATURE: return "premature";
            case tbank::PvcFilter::VOTE:      return "vote";
            default:                          return "no";
            }
        }
        inline std::string prematurePctAgg(const tbank::BankTemplate& t) {
            const uint32_t n = (uint32_t)t.members.size();
            if (n == 0) return "na";
            char b[16]; std::snprintf(b, sizeof b, "%.1f",
                100.0 * double(t.n_premature_members) / double(n)); return b;
        }
        inline std::string votedPctAgg(const tbank::BankTemplate& t) {
            const uint32_t n = (uint32_t)t.members.size();
            if (n == 0) return "na";
            char b[16]; std::snprintf(b, sizeof b, "%.1f",
                100.0 * double(t.n_voted_members) / double(n)); return b;
        }
        inline std::string tukeyPctAgg(const tbank::BankTemplate& t) {
            const int32_t n = (int32_t)t.members.size();
            const int32_t flagged = t.n_premature_members + t.n_voted_members;
            const int32_t elig = (n > flagged) ? (n - flagged) : 0;
            if (elig <= 0) return "na";
            char b[16]; std::snprintf(b, sizeof b, "%.1f",
                100.0 * double(t.n_tukey_members) / double(elig)); return b;
        }
        inline std::string blendPctAgg(const tbank::BankTemplate& t) {
            const uint32_t n = (uint32_t)t.members.size();
            if (n == 0) return "na";
            char b[16]; std::snprintf(b, sizeof b, "%.1f",
                100.0 * double(t.n_blended_members) / double(n)); return b;
        }

        // Instantaneous rate implied by this template's member R-R intervals.
        inline std::string bpmPerTemplate(const tbank::BankTemplate& t) {
            if (t.mean_rr_ms <= 0.0) return "na";
            char b[16]; std::snprintf(b, sizeof b, "%.1f", 60000.0 / t.mean_rr_ms);
            return b;
        }


        // Row 5. Tukey runs only on beats NOT flagged premature: a premature
        // beat is already excluded from the reference set, so there is nothing
        // for Tukey to decide about it. `not_eligible` is therefore a real and
        // common state, distinct from "tested and kept" -- and the specific
        // rejection reason is kept because a beat classified pqrst and rejected
        // on RR LENGTH is very likely ectopy the classifier missed, which is the
        // cheapest estimate of classifier recall available.
        inline const char* tukeyWord(tbank::TukeyOutcome t) {
            switch (t) {
            case tbank::TukeyOutcome::KEPT:           return "kept";
            case tbank::TukeyOutcome::REJ_RR_LENGTH:  return "removed_rr";
            case tbank::TukeyOutcome::REJ_AMPLITUDE:  return "removed_amplitude";
            case tbank::TukeyOutcome::REJ_R_LOCATION: return "removed_r_location";
            case tbank::TukeyOutcome::REJ_WAVE_SCORE: return "removed_wave_score";
            default:                                  return "not_eligible";
            }
        }

    }  // namespace detail

    struct ChannelBlock {
        const char* channel = "CH1";

        // ---- EVERYTHING HERE IS PER BIN, IN BIN ORDER, AND BY POINTER ------
        //
        // The projected bank, the beat matrix and the slice map all live inside
        // the caller's TemplateInfo, one object per bin. Pointers rather than
        // copies because the alternative duplicates every beat waveform in the
        // record once per channel block.
        //
        // A null entry means this channel has nothing for that bin -- a two-lead
        // record, or a bin the pulse channel found no usable beat in. Its slices
        // still get columns, with no samples: a slice that produced no beat on
        // this channel is a fact the archive has to be able to state.
        std::vector<const bin_pipeline::ChannelOutput*> per_bin;

        // This channel's captured beats for the bin, [row][sample], on the
        // shared axis. Indexed by ROW, which is not a slice; local_of_slice
        // below is the only route from a column to a waveform.
        std::vector<const std::vector<std::vector<double>>*> beats;

        // slice -> row in beats[bin], or -1 when this channel produced no beat
        // for that slice. This is jbank::ChannelBeats::local_of_slice handed
        // straight through -- the same map the partition resolved its members
        // with, so a column's descriptors and its samples cannot end up
        // describing different heartbeats.
        std::vector<const std::vector<int32_t>*> local_of_slice;

        // jbank::ExcludeReason per slice, shared by every channel because the
        // exclusion is a property of the beat, not of the lead.
        std::vector<const std::vector<uint8_t>*> excluded_reason;

        // Anchor column on the shared axis, per bin: r_col for the ECG leads,
        // the systolic peak column for PPG.
        std::vector<int> r_col;

        bool empty() const { return per_bin.empty(); }
        size_t nBins() const { return per_bin.size(); }

        const bin_pipeline::ChannelOutput* out(size_t b) const {
            return (b < per_bin.size()) ? per_bin[b] : nullptr;
        }
        const std::vector<std::vector<double>>* binBeats(size_t b) const {
            return (b < beats.size()) ? beats[b] : nullptr;
        }
        const std::vector<uint8_t>* reasons(size_t b) const {
            return (b < excluded_reason.size()) ? excluded_reason[b] : nullptr;
        }
        int rCol(size_t b) const { return (b < r_col.size()) ? r_col[b] : -1; }

        // slice -> row, or -1. Identity when no map was supplied, which is what
        // a caller whose rows and slices coincide has.
        int rowOf(size_t b, size_t slice) const {
            const std::vector<int32_t>* m =
                (b < local_of_slice.size()) ? local_of_slice[b] : nullptr;
            if (!m) return static_cast<int>(slice);
            return (slice < m->size()) ? (*m)[slice] : -1;
        }
    };

    // One column per beat. Channels are emitted as successive row blocks, each
    // preceded by a `channel` row naming it, so a reader can split on that row
    // rather than parsing three files.
    // ---- SlotMap IS GONE ---------------------------------------------------
    //
    // It inverted kept_idx to turn an ALIGNED BEAT INDEX into a row of the
    // captured beat matrix, because a column of _beats was one aligned beat.
    // A column is now one SLICE, and the slice -> row map already exists as
    // jbank::ChannelBeats::local_of_slice -- built by the partition itself,
    // from the same forward maps, before any of this runs. ChannelBlock
    // carries it directly (rowOf), so there is nothing left to invert and no
    // second inversion that could disagree with the first.

    // =====================================================================
    // <stem>_bins.csv -- ONE ROW PER BIN
    // =====================================================================
    //
    // The 4.5 acceptance test asks for "per-bin category percentages", and
    // nothing produced them: BinCounts held n_regular / n_ectopic / n_noise and
    // no writer ever read them. Everything else in this file is per template or
    // per beat, so a per-bin fact had nowhere to go.
    //
    // COUNTS ARE PASSED IN, PERCENTAGES ARE COMPUTED HERE, from those same
    // counts on the same row. A caller that computed both could hand over a
    // percentage that disagrees with its own numerator, and nothing downstream
    // could tell.
    //
    // ONE ROW PER BIN, not per bin per channel: every quantity on it is a
    // property of the beats -- which category they are, which group they landed
    // in, why they left the average -- and the whole point of the joint
    // partition is that those have one answer per bin rather than four.
    struct BinRow {
        uint32_t bin = 0;
        uint32_t n_slices = 0;      // R-pairs in the bin
        uint32_t n_became_beat = 0; // slices with a beat on at least one channel

        // 4.5 categories, over slices.
        uint32_t n_regular = 0;
        uint32_t n_ectopic = 0;
        uint32_t n_noise = 0;

        // Partition.
        uint32_t n_groups = 0;
        uint32_t n_assigned = 0;
        uint32_t n_unscorable = 0;
        uint32_t n_spawns = 0;
        uint32_t n_merges = 0;
        uint32_t n_cap_raises = 0;

        // Why beats are not in their group's average, by jbank::ExcludeReason.
        uint32_t ex_not_member = 0;
        uint32_t ex_category = 0;
        uint32_t ex_premature = 0;
        uint32_t ex_vote = 0;
        uint32_t ex_tukey = 0;
        uint32_t n_kept = 0;

        // Prematurity, from the filter rather than from the marks. The PAIR is
        // the diagnostic: high premature with zero vote is alternating ectopy
        // (in perfect bigeminy exactly 4 of any 8 beats are ectopic, the count
        // never reaches 5, and the vote is arithmetically unable to fire); a
        // high vote count means runs.
        uint32_t n_premature = 0;
        uint32_t n_vote_only = 0;

        // 4.6 substitutions, over the beats PREMATURE or VOTE rejected.
        // n_substituted counts BEATS; n_sub_blends counts beat-channels, and is
        // the larger of the two whenever a beat was substituted on more than
        // one lead.
        //
        // n_sub_too_bad is gone with the correlation band it counted. There is
        // no 0.60 floor any more: 4.6 substitutes every flagged beat, because a
        // beat too bad to blend still leaves a hole in the RR series.
        uint32_t n_substituted = 0;
        uint32_t n_sub_channel_blends = 0;

        // -1 = not computed. 0 = no confirmed PVC template in this bin, which
        // is the normal state before marking and is NOT "monomorphic".
        int32_t  polymorphy_count = -1;
        uint32_t n_unconfirmed_groups = 0;

        // seed_pool::seedBasisName for the Phase 1 reference on CH1. A bin whose
        // basis is not sinus_only has a reference that is not purely sinus.
        const char* seed_basis = "";
    };

    inline bool writeBins(const std::vector<BinRow>& rows) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_bins.csv", std::ios::trunc);
        if (!f) return false;

        f << "bin,n_slices,n_became_beat,"
            "n_regular,n_ectopic,n_noise,"
            "pct_regular,pct_ectopic,pct_noise,"
            "n_groups,n_assigned,n_unscorable,n_spawns,n_merges,n_cap_raises,"
            "ex_not_member,ex_category,ex_premature,ex_vote,ex_tukey,n_kept,"
            "pct_excluded,"
            "n_premature,n_vote_only,looks_alternating,"
            "n_substituted,n_sub_blends,"
            "polymorphy_count,polymorphy,n_unconfirmed_groups,seed_basis\n";

        f.setf(std::ios::fixed);
        for (const BinRow& r : rows) {
            const double den = (r.n_slices > 0) ? double(r.n_slices) : 1.0;
            const uint32_t members = r.ex_category + r.ex_premature
                + r.ex_vote + r.ex_tukey + r.n_kept;
            const double mden = (members > 0) ? double(members) : 1.0;

            // Percentages of the SLICE count, not of the assigned count, so the
            // three add to 100 and a bin where nothing was assigned still
            // reports its categories instead of dividing by zero.
            f << r.bin << ',' << r.n_slices << ',' << r.n_became_beat << ','
                << r.n_regular << ',' << r.n_ectopic << ',' << r.n_noise << ',';
            f.precision(2);
            f << 100.0 * r.n_regular / den << ','
                << 100.0 * r.n_ectopic / den << ','
                << 100.0 * r.n_noise / den << ',';
            f << r.n_groups << ',' << r.n_assigned << ',' << r.n_unscorable
                << ',' << r.n_spawns << ',' << r.n_merges << ','
                << r.n_cap_raises << ','
                << r.ex_not_member << ',' << r.ex_category << ','
                << r.ex_premature << ',' << r.ex_vote << ',' << r.ex_tukey
                << ',' << r.n_kept << ',';
            // Of the beats that HAVE a group: what fraction left the average.
            f.precision(2);
            f << 100.0 * double(members - r.n_kept) / mden << ',';
            f << r.n_premature << ',' << r.n_vote_only << ','
                << ((r.n_premature >= 8 && r.n_vote_only == 0) ? 1 : 0) << ',';
            f << r.n_substituted << ',' << r.n_sub_channel_blends << ',';
            f << r.polymorphy_count << ',';
            if (r.polymorphy_count < 0)      f << "unknown";
            else if (r.polymorphy_count >= 2) f << "polymorphic";
            else if (r.polymorphy_count == 1) f << "monomorphic";
            else                              f << "none_confirmed";
            f << ',' << r.n_unconfirmed_groups << ',' << r.seed_basis << '\n';
        }
        return static_cast<bool>(f);
    }

    // =====================================================================
    // <stem>_nsvt.csv -- ONE ROW PER RUN
    // =====================================================================
    //
    // Runs are RECORD-LEVEL, not per bin, which is the whole reason cross-bin
    // global template identity exists: a run that straddles a bin boundary is
    // invisible if you scan per-bin group ids, because its beats carry two
    // different local numbers. crosses_bin is recorded because such a run is
    // assembled from two bins' groups that a correlation judged identical, and
    // that judgement is worth being able to audit.
    //
    // SUSTAINED RUNS ARE IN THIS FILE TOO, with the flag set. The spec says a
    // run of 30 s or more "is sustained VT and is escalated rather than logged
    // as NSVT", but a run silently absent from every file is the worst outcome
    // available -- an escalation path can filter on the column.
    struct NsvtRow {
        uint32_t start_beat = 0;
        uint32_t length = 0;        // beats
        int32_t  global_template = -1;
        int32_t  subtype = -1;
        uint8_t  label_code = 0;
        double   mean_cycle_ms = 0.0;
        double   max_cycle_ms = 0.0;
        double   rate_bpm = 0.0;
        double   duration_ms = 0.0;
        uint8_t  sustained = 0;
        uint8_t  crosses_bin = 0;
        uint32_t first_bin = 0;
        uint32_t last_bin = 0;
    };

    // Written even when empty -- a header-only file says the detector ran and
    // found nothing, which is a different statement from no file at all. The
    // 4.6 acceptance test "a record with isolated unifocal PVCs produces no
    // runs" is only evidence if the run can be distinguished from the detector
    // never having executed.
    inline bool writeNsvt(const std::vector<NsvtRow>& runs,
        uint32_t polymorphic_candidates = 0)
    {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_nsvt.csv", std::ios::trunc);
        if (!f) return false;

        // The blind spot, as a number rather than a comment. The criterion is
        // three or more consecutive beats on THE SAME template, so polymorphic
        // VT and torsades -- which change morphology beat to beat -- scatter
        // across templates and never form a run. A non-zero count here on a
        // record with no detected runs is the signal that the same-template
        // criterion is costing something real.
        f << "# polymorphic_candidates," << polymorphic_candidates << '\n';
        f << "start_beat,length,global_template,subtype,label_code,"
            "mean_cycle_ms,max_cycle_ms,rate_bpm,duration_ms,"
            "sustained,crosses_bin,first_bin,last_bin\n";
        f.setf(std::ios::fixed);
        f.precision(1);
        for (const NsvtRow& r : runs)
            f << r.start_beat << ',' << r.length << ',' << r.global_template
            << ',' << r.subtype << ',' << int(r.label_code) << ','
            << r.mean_cycle_ms << ',' << r.max_cycle_ms << ','
            << r.rate_bpm << ',' << r.duration_ms << ','
            << int(r.sustained) << ',' << int(r.crosses_bin) << ','
            << r.first_bin << ',' << r.last_bin << '\n';
        return static_cast<bool>(f);
    }

    namespace detail {

        // Row 4, aggregated. `vote` and `mixed` are kept distinct from
        // `premature` because the two filter paths have disjoint competence: the
        // raw test catches isolated ectopy, the 5-of-8 vote only fires inside
        // runs, and in perfect bigeminy exactly 4 of any 8 beats are ectopic so
        // the vote never reaches 5 at all. A template reading `vote` was
        // therefore found by a mechanism that cannot see alternating rhythms.
        inline std::string prematureWordAgg(const tbank::BankTemplate& t) {
            const uint32_t n = static_cast<uint32_t>(t.members.size());
            if (n == 0) return "na";
            const uint32_t prem = t.n_premature_members;
            const uint32_t vote = t.n_voted_members;
            if (prem * 2 > n) return "premature";
            if ((prem + vote) * 2 > n) return "vote";
            if (prem + vote > 0) return "mixed";
            return "no";
        }

        // Row 5, aggregated. `not_eligible` dominates on ectopic templates by
        // construction: Tukey runs only on beats NOT flagged premature, since a
        // premature beat is already excluded from the reference set and there is
        // nothing left for Tukey to decide about it.
        // FROM THE TEMPLATE'S OWN COUNTS, NOT BY INDEXING THE FLAGS.
        //
        // This used to walk t.members and read out.flags[m]. Under the joint
        // partition those are different index spaces -- members is in the
        // channel's local row space, flags is in slice space -- so every lookup
        // returned some other beat's Tukey verdict, plausibly and wrongly. The
        // counts are computed in jbank::projectToChannel, where the group's
        // slice membership and the per-slice flags are both in hand.
        //
        // ELIGIBLE means "not premature": Tukey only ever ran on the beats
        // premature removal left behind, so a premature member has no Tukey
        // verdict to report and must not be counted as one that passed.
        inline std::string tukeyWordAgg(const tbank::BankTemplate& t,
            const bin_pipeline::ChannelOutput&)
        {
            const int32_t n = static_cast<int32_t>(t.members.size());
            const int32_t flagged = t.n_premature_members + t.n_voted_members;
            const int32_t eligible = (n > flagged) ? (n - flagged) : 0;
            if (eligible <= 0) return "not_eligible";
            if (t.n_tukey_members <= 0) return "kept";
            if (t.n_tukey_members >= eligible) return "removed";
            return "partial";
        }

        // -----------------------------------------------------------------
        // Does this template belong in the TEMPLATES file?
        // -----------------------------------------------------------------
        //
        // Two exclusions, both on the templates file ONLY. The beats file still
        // carries every one of these beats with its full descriptors, so nothing
        // is lost -- the question here is narrower: is this row a MORPHOLOGY
        // worth a column, or is it a record of beats the pipeline already
        // decided against?
        //
        //   entirely Tukey-removed  every eligible member failed a Tukey pass.
        //                           The pipeline excluded all of them, so the
        //                           median is an average of rejected beats and
        //                           has no standing as a template.
        //
        //   entirely premature      every member is premature. Not majority --
        //                           ENTIRELY. A majority-premature template is
        //                           already reported as ectopic and denied
        //                           landmarks, and it is a real morphology worth
        //                           a column; one whose every member is
        //                           premature is a run, not a shape the operator
        //                           needs to see averaged here.
        //
        // Slot 0 is exempt. It is the Phase 1 sinus seed and the reference every
        // other column is implicitly compared against, so dropping it would
        // leave a bin's templates file describing a bin with no sinus template
        // -- which reads as "this bin had no normal beats" rather than "the seed
        // had a bad night".
        inline bool belongsInTemplatesFile(const tbank::BankTemplate& tp,
            const bin_pipeline::ChannelOutput& out, const char** why)
        {
            if (tp.isSeed()) return true;
            if (tukeyWordAgg(tp, out) == "removed") {
                if (why) *why = "all members Tukey-removed";
                return false;
            }
            const uint32_t n = static_cast<uint32_t>(tp.members.size());
            if (n > 0 && tp.n_premature_members == n) {
                if (why) *why = "all members premature";
                return false;
            }
            return true;
        }

    }  // namespace detail

    inline bool writeTemplates(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_templates.csv", std::ios::trunc);
        if (!f) return false;

        // ---- ONE TABLE, A COLUMN PER (bin, template, channel) -------------
        //
        // This wrote one BLOCK PER CHANNEL, stacked vertically: every CH1
        // column, then a fresh `channel` row, then CH2, then CH3, then PPG at
        // the bottom. Two things were wrong with that.
        //
        // A GROUP'S FOUR FACES WERE IN FOUR SEPARATE TABLES. A group is ONE set
        // of heartbeats seen on three leads and the pulse, and reading one meant
        // finding its column in each table and joining them by (bin, template)
        // by hand -- and NOT by position, because the blocks had different
        // column counts. A channel with no usable beats in a bin contributed no
        // column at all, so the nth column of the CH1 block and the nth of the
        // PPG block were generally different groups.
        //
        // AND AN ABSENT CHANNEL VANISHED ENTIRELY. A block whose every template
        // had an empty waveform produced no rows -- not even a header -- so
        // "this record has no usable pulse" and "the pulse never got written"
        // were indistinguishable in the file. On a record whose pulse filter
        // rejected all but a few hundred beats that cost several rounds of
        // looking in the wrong place.
        //
        // Now: columns are ordered by (bin, template, channel), so a group is a
        // contiguous run -- CH1, CH2, CH3, PPG for bin 0 template A, then bin 0
        // template B, and so on. EVERY channel gets a column in every run, even
        // when it has no waveform: the metadata rows say what happened and the
        // waveform rows are simply empty. A gap in a run is now a stated fact
        // rather than a missing column.
        struct Col {
            std::string category, bin, name, premature, voted, tukey, blended,
                bpm, members, excluded, share, too_few;
            const std::vector<double>* wave = nullptr;   // null => no waveform
            const char* channel = "";
            const char* channel_type = "";              // "ecg" or "ppg"
            int  r_col = -1;                            // -1 => none for this channel
            size_t binIdx = 0;
            int  tmplIdx = 0;
            int  chOrder = 0;                           // block order, the tiebreak
        };

        std::vector<Col> cols;
        size_t dropped = 0;   // reported, so a thin table is explicable

        // How many templates any channel holds for a bin. The banks are allowed
        // to disagree -- a morphology separable on CH1 may not be on CH2 -- so
        // the run length for a bin is the widest of them, and a channel that
        // runs short gets a placeholder column rather than the run going ragged.
        auto slotsInBin = [&](size_t b) -> int {
            int n = 0;
            for (const ChannelBlock& blk : blocks) {
                if (blk.empty()) continue;
                const bin_pipeline::ChannelOutput* op = blk.out(b);
                if (!op) continue;
                n = std::max(n, static_cast<int>(op->bank.size()));
            }
            return n;
            };

        size_t nBins = 0;
        for (const ChannelBlock& blk : blocks)
            if (!blk.empty()) nBins = std::max(nBins, blk.nBins());

        for (size_t b = 0; b < nBins; ++b) {
            const int nSlots = slotsInBin(b);
            for (int t = 0; t < nSlots; ++t) {
                for (size_t bi = 0; bi < blocks.size(); ++bi) {
                    const ChannelBlock& blk = blocks[bi];
                    // Which minimum applies. Compared on the block's channel
                    // name rather than its index, because the caller decides
                    // the block order and a positional assumption here would
                    // silently apply the ECG minimum to the pulse channel on
                    // any reordering.
                    const bool isPpgBlock =
                        (blk.channel && std::string(blk.channel) == "PPG");

                    Col c;
                    c.channel = blk.channel ? blk.channel : "";
                    c.channel_type = isPpgBlock ? "ppg" : "ecg";
                    c.chOrder = static_cast<int>(bi);
                    c.binIdx = b;
                    c.tmplIdx = t;
                    c.bin = std::to_string(b);

                    const bin_pipeline::ChannelOutput* op =
                        blk.empty() ? nullptr : blk.out(b);

                    // ---- THE PLACEHOLDER CASES -----------------------------
                    // A column is emitted for every (bin, template, channel),
                    // so the reasons a face is missing have to be stated rather
                    // than implied by absence. Each `name` below is that
                    // statement, and every other cell stays blank.
                    if (!op) {
                        c.name = "(no data this bin)";
                        cols.push_back(std::move(c));
                        continue;
                    }
                    if (t >= op->bank.size()) {
                        c.name = "(no slot on this channel)";
                        cols.push_back(std::move(c));
                        continue;
                    }
                    const tbank::BankTemplate& tp = op->bank.templates[t];
                    if (tp.tmpl.empty()) {
                        // The case that used to make the whole PPG block
                        // disappear: a slot that exists with no waveform, which
                        // is what a bin whose pulse filter produced nothing
                        // looks like.
                        c.name = "(no waveform)";
                        c.members = std::to_string(tp.memberCount());
                        cols.push_back(std::move(c));
                        continue;
                    }
                    const char* why = nullptr;
                    if (!detail::belongsInTemplatesFile(tp, *op, &why)) {
                        // Excluded from THIS file only; the beats file keeps
                        // them. Stated in place rather than counted and
                        // omitted, now that there is a column to state it in.
                        ++dropped;
                        c.name = why ? why : "(excluded)";
                        c.members = std::to_string(tp.memberCount());
                        cols.push_back(std::move(c));
                        continue;
                    }

                    const std::vector<uint8_t> letters = detail::letterRanks(op->bank);
                    c.r_col = blk.r_col.empty() ? -1 : blk.rCol(b);
                    c.category = detail::categoryWord(tp.presumedCategory());
                    c.name = detail::templateName(tp,
                        (static_cast<size_t>(t) < letters.size()) ? letters[t] : 0);
                    c.premature = detail::prematurePctAgg(tp);
                    c.voted = detail::votedPctAgg(tp);
                    c.tukey = detail::tukeyPctAgg(tp);
                    c.blended = detail::blendPctAgg(tp);
                    c.members = std::to_string(tp.memberCount());
                    c.excluded = std::to_string(tp.excludedCount());
                    c.bpm = detail::bpmPerTemplate(tp);
                    // FLAGGED, NOT OMITTED, and now flagged where the sibling
                    // channels can be seen beside it: a CH1 column reading "no"
                    // next to its PPG column reading "yes" describes a panel
                    // the operator never saw as valid, and that disagreement
                    // used to require joining two tables to notice.
                    c.too_few = tp.tooFewBeats(isPpgBlock) ? "yes" : "no";
                    c.share = std::to_string(op->bank.beatShare(t));
                    c.wave = &tp.tmpl;
                    cols.push_back(std::move(c));
                }
            }
        }

        if (dropped)
            std::fprintf(stderr,
                "  [morphology] %zu template(s) marked excluded in _templates "
                "(all members Tukey-removed or all premature); their beats "
                "remain in _beats\n", dropped);

        if (cols.empty()) return true;

        auto row = [&](const char* label, std::string Col::* field) {
            f << label;
            for (auto& c : cols) f << ',' << (c.*field);
            f << '\n';
            };

        // r_col first: it is the axis every ECG column's samples are on, and
        // blank for the pulse, which is foot-anchored and has no R column.
        {
            bool anyR = false;
            for (auto& c : cols) if (c.r_col >= 0) { anyR = true; break; }
            if (anyR) {
                f << "r_col";
                for (auto& c : cols) { f << ','; if (c.r_col >= 0) f << c.r_col; }
                f << '\n';
            }
        }

        row("category", &Col::category);
        row("bin", &Col::bin);
        row("template", &Col::name);
        row("premature_pct", &Col::premature);
        row("voted_pvc_pct", &Col::voted);
        row("bpm", &Col::bpm);
        row("tukey_pct", &Col::tukey);
        row("blended_pct", &Col::blended);
        // NO `confirmed` ROW, AND NO `marking` ROW.
        //
        // confirmed could never be right here: this file is written by
        // GenerateTemplatesFast, inside prepareViewerJob, BEFORE the marking
        // window exists -- so every template read "presumed" in every record no
        // matter how much marking followed, which looked like a bug in the
        // confirmation logic rather than a promise this file cannot keep.
        //
        // Nor can it be fixed by rewriting afterwards: the ChannelBlocks this
        // function takes hold raw pointers into local_of_slice,
        // excluded_reason and TemplateInfo::bank_by_channel, all owned by
        // GenerateTemplatesFast's own frame, so none of them exist once it
        // returns.
        //
        // The VIEWER writes <id>_template_confirmations.csv instead, from the
        // bins it owns and after the operator is done. Join on
        // (bin, channel, template).
        row("too_few_beats", &Col::too_few);
        row("n_members", &Col::members);
        row("n_excluded", &Col::excluded);
        row("beat_share", &Col::share);

        // ---- THE CHANNEL ROWS, LAST OF THE METADATA ---------------------
        // After beat_share and before the waveform, not at the head of the
        // table: the waveform below is thousands of rows tall, because an ECG
        // axis is framed on the bin's LONGEST RR. A channel name at the top
        // scrolls away long before the samples run out.
        f << "channel";
        for (auto& c : cols) f << ',' << c.channel;
        f << '\n';
        // ecg OR ppg. Not redundant with `channel`: that names the LEAD,
        // this names the KIND. Everything that has to treat the pulse
        // differently -- which minimum applies, the perfusion-index
        // normalization, which landmarks exist -- was deciding it by
        // string-matching "PPG" against the channel name, independently in
        // several places.
        f << "channel_type";
        for (auto& c : cols) f << ',' << c.channel_type;
        f << '\n';

        // Waveform: the template median, one row per sample index. Empty cell
        // for NaN -- the shared axis is NaN-padded at both ends, and 0 is a
        // real amplitude -- and empty for a whole column that has no waveform.
        //
        // THE WIDTH IS THE MAX ACROSS CHANNELS. An ECG axis is framed on the
        // bin's longest RR and a pulse axis on its own foot-to-foot interval,
        // so they differ by thousands of samples; the shorter columns run out
        // and go empty, which is what a reader keyed on the sample index
        // expects.
        size_t width = 0;
        for (auto& c : cols)
            if (c.wave) width = std::max(width, c.wave->size());
        for (size_t sIdx = 0; sIdx < width; ++sIdx) {
            f << sIdx;
            for (auto& c : cols) {
                f << ',';
                if (c.wave && sIdx < c.wave->size()
                    && !std::isnan((*c.wave)[sIdx]))
                    f << (*c.wave)[sIdx];
            }
            f << '\n';
        }
        return true;
    }

    // =====================================================================
    // Binary forms
    // =====================================================================
    //
    // Layout, both files:
    //
    //   [uint32 version][uint32 nChannelBlocks]
    //   per block:
    //     [uint32 nameLen][char[nameLen] channel]
    //     [uint32 width]              samples per column, 0 = no waveform
    //     [uint64 nColumns]
    //     per column: a fixed descriptor record, then `width` doubles
    //                 _templates.bin only, then:
    //                   [uint32 n][n x uint32]  members
    //                   [uint32 n][n x uint32]  members_clean
    //                   [TemplateTrailer]       64 bytes, raw
    //                   [uint32 n][n x double]  tmpl_iqr
    //
    // NaN is written as NaN, not as a sentinel. The CSV writes an empty cell
    // because a literal would force consumers to special-case a string in a
    // numeric column; binary has no such problem, and IEEE NaN survives a
    // read/write pair exactly.

    // NO MAGIC. There was a kBeatsMagic / kTemplatesMagic pair, on the argument
    // that a wrong-type file should fail immediately rather than be read as a
    // plausible column count. The failure it was guarding against turned out to
    // be a filename collision in another writer, fixed at the source -- and the
    // two readers here take an explicit path and are named for their file, so
    // feeding one the other's path is a typo, not a class of failure.
    //
    // ONE FORMAT VERSION, AND IT STAYS 1. There was a v1..v3 history here
    // describing migrations that never happened -- nothing has been released,
    // so there are no old files to be compatible with. A field gets added and
    // every archive is regenerated. The number is kept only so a future release
    // has somewhere to start counting.
    inline constexpr uint32_t bin_version = 1;

    // One record per beat column. Fixed size, so a reader can stride over
    // descriptors without parsing them.
    // ONE RECORD PER SLICE. A slice is an R-pair -- one heartbeat interval --
    // and it exists whether or not any channel produced a usable waveform for
    // it. That is the whole reason the column key changed: the old per-beat
    // column could not represent a dropped beat at all, so a slice the pulse
    // channel lost simply vanished from the PPG block and column k of CH1 and
    // column k of PPG were different heartbeats. Now column k is slice k in
    // every block, and became_beat says whether this channel has samples for it.
    struct BeatRecord {
        uint32_t bin = 0;
        uint8_t  category = 0;   // tbank::Category
        uint8_t  premature = 0;   // tbank::PvcFilter
        uint8_t  tukey = 0;   // tbank::TukeyOutcome
        uint8_t  confirmed = 0;   // 0 presumed, 1 confirmed, 2 n/a
        uint8_t  label_code = 0;   // annotation code, 0 = unlabeled
        uint8_t  letter = 0;   // separation index: 0 = A
        // 1 when this channel produced a beat for this slice, 0 when it did not
        // -- dropped by the slicer, by the baseline filter, or by a failed
        // fiducial. A 0 column carries a full row of NaN and is NOT an error.
        uint8_t  became_beat = 0;
        // jbank::ExcludeReason. 0 kept, 1 not_a_member, 2 premature, 3 vote,
        // 4/5/6 Tukey r_location / amplitude / wave_score. This is the field
        // that says a beat was excluded from its template's average WITHOUT
        // being removed from the record.
        uint8_t  excluded = 0;
        int32_t  template_id = -1;  // bank slot, -2 unscorable, -1 unassigned
        int32_t  r_col = -1;
    };

    // One record per template column.
    struct TemplateRecord {
        uint32_t bin = 0;
        uint8_t  category = 0;   // presumed tbank::Category
        uint8_t  premature = 0;   // 0 no, 1 mixed, 2 vote, 3 premature
        uint8_t  tukey = 0;   // 0 not_eligible, 1 kept, 2 partial, 3 removed
        uint8_t  confirmed = 0;
        uint8_t  label_code = 0;
        uint8_t  letter = 0;
        uint8_t  landmark_marked = 0;   // 1 = category 1, gets landmarks
        int32_t  template_id = -1;
        int32_t  r_col = -1;
        uint32_t n_members = 0;
        // Members excluded from the average: premature plus Tukey-rejected.
        // n_members stays the FULL membership, so n_members - n_excluded is what
        // the waveform in this column was actually built from. Reported as two
        // numbers rather than one because the difference is the interesting part.
        uint32_t n_excluded = 0;
        // 1 when cleanCount() is below the configured minimum for this
        // channel (tbank::minBeatsEcg / minBeatsPpg). The template is still
        // written in full -- flagged, not omitted -- because a suppressed row
        // would leave its beats unaccounted for and give a reader no way to
        // distinguish a thin template from one that never existed.
        uint8_t  too_few_beats = 0;
        double   beat_share = 0.0;
    };

    // THE RECORDS ARE WRITTEN RAW, so their padding is part of the on-disk
    // format. Both are pinned here: a compiler that packs them differently
    // fails the build instead of silently emitting a file no reader can parse.
    // read_morphology_bin.py hardcodes these same layouts, and the padding
    // offsets are the reason it must -- see the struct format strings there.
    //
    //   BeatRecord      uint32 bin, 8x uint8, int32 template_id,
    //                   int32 r_col                                      = 20
    //   TemplateRecord  uint32 bin, 7x uint8, 1 pad, int32 template_id,
    //                   int32 r_col, uint32 n_members, uint32 n_excluded,
    //                   4 pad, double beat_share                         = 40
    //
    // THE TWO GREW DIFFERENTLY, AND THIS IS THE DANGEROUS PART. BeatRecord's
    // two new uint8 fields fit in padding it already had, so its SIZE IS
    // UNCHANGED at 20: a v1 reader on a v2 beats file strides correctly through
    // every record and silently reports padding as became_beat. Nothing about
    // the file's shape reveals the change, so the version field is the only
    // thing that can, and a reader that does not check it is wrong without
    // failing. TemplateRecord had no spare room, so it grew 32 -> 40 and a v1
    // reader walks off alignment on the second record -- the loud failure, and
    // the preferable one.
    static_assert(sizeof(BeatRecord) == 20,
        "BeatRecord layout changed; update read_morphology_bin.py");
    static_assert(sizeof(TemplateRecord) == 40,
        "TemplateRecord layout changed; update read_morphology_bin.py");

    namespace detail {

        // NOT NAMED w / r. They were, and writeBins has a `const BinRow& r` in
        // its row loop -- so inside that loop `r` is the row and detail::r is
        // shadowed. Harmless today because the two are never used in one scope,
        // and a trap the first time anyone adds a read there.
        inline void writeRaw(std::ofstream& f, const void* p, size_t n) {
            f.write(reinterpret_cast<const char*>(p), static_cast<std::streamsize>(n));
        }
        inline bool readRaw(std::ifstream& f, void* p, size_t n) {
            return static_cast<bool>(
                f.read(reinterpret_cast<char*>(p), static_cast<std::streamsize>(n)));
        }

        // ---- LENGTH-PREFIXED VECTORS -----------------------------------
        //
        // [uint32 n][n x element]. Bounded on read, so a corrupt or misaligned
        // length fails here rather than reaching the allocator with garbage.
        // The cap is generous against any real bin: a 15-minute bin holds low
        // thousands of beats.
        inline constexpr uint32_t kMaxVecLen = 4u * 1024u * 1024u;

        inline void writeVecU32(std::ofstream& f, const std::vector<uint32_t>& v) {
            const uint32_t n = static_cast<uint32_t>(v.size());
            writeRaw(f, &n, 4);
            if (n) writeRaw(f, v.data(), static_cast<size_t>(n) * 4);
        }
        inline bool readVecU32(std::ifstream& f, std::vector<uint32_t>& v) {
            uint32_t n = 0;
            if (!readRaw(f, &n, 4) || n > kMaxVecLen) return false;
            v.resize(n);
            return n ? readRaw(f, v.data(), static_cast<size_t>(n) * 4) : true;
        }
        inline void writeVecD(std::ofstream& f, const std::vector<double>& v) {
            const uint32_t n = static_cast<uint32_t>(v.size());
            writeRaw(f, &n, 4);
            if (n) writeRaw(f, v.data(), static_cast<size_t>(n) * 8);
        }
        inline bool readVecD(std::ifstream& f, std::vector<double>& v) {
            uint32_t n = 0;
            if (!readRaw(f, &n, 4) || n > kMaxVecLen) return false;
            v.resize(n);
            return n ? readRaw(f, v.data(), static_cast<size_t>(n) * 8) : true;
        }

        // ---- EVERYTHING TemplateRecord HAS NO ROOM FOR -----------------
        //
        // Written raw and appended after the waveform, because
        // TemplateRecord's size is pinned by static_assert and
        // read_morphology_bin.py strides on it -- growing that struct would
        // break the python reader on every record.
        //
        // This is what lets the archive RESTORE a partition rather than only
        // describe one. n_members says how many beats a template holds; the
        // census says what kind they were, which is what presumedCategory()
        // reads; and subtype / spawn_seq are what letters resolve against.
        //
        // THE LAST FOUR ARE PER BANK, NOT PER TEMPLATE. Records are flat, one
        // per template with no per-bin header, so there is nowhere to hang
        // per-bank fields. 16 bytes on each record, and a reader takes them
        // from any record of the bank -- they agree by construction.
        struct TemplateTrailer {
            int32_t  subtype = -1;
            uint32_t spawn_seq = 0;
            int32_t  n_ppg_members = -1;
            int32_t  n_tukey_members = 0;
            uint32_t n_blended_members = 0;
            uint32_t n_premature_members = 0;
            uint32_t n_voted_members = 0;
            uint32_t n_noise_members = 0;
            uint8_t  marked_invalid_template = 0;
            uint8_t  operator_state = 0;
            uint8_t  confirmed_by_operator = 0;
            uint8_t  pad0 = 0;
            int32_t  configured_cap = 0;
            int32_t  effective_cap = 0;
            uint32_t next_spawn_seq = 0;
            uint32_t assigned_beats = 0;
            uint32_t pad1 = 0;            // keeps mean_rr_ms 8-byte aligned
            double   mean_rr_ms = 0.0;
        };
        // Raw-written, so the padding is part of the format. Pinned for the same
        // reason the two record structs are.
        static_assert(sizeof(TemplateTrailer) == 64,
            "TemplateTrailer layout changed; update read_morphology_bin.py");

        inline void writeBlockHeader(std::ofstream& f, const char* channel,
            uint32_t width, uint64_t nCols) {
            const uint32_t len = static_cast<uint32_t>(std::string(channel).size());
            writeRaw(f, &len, 4);
            writeRaw(f, channel, len);
            writeRaw(f, &width, 4);
            writeRaw(f, &nCols, 8);
        }

        // Letter index: the subtype for a confirmed template, spawn order for an
        // unlabeled one. They differ because subtypes are issued per class in
        // order of first confirmation while spawn order is order of first
        // appearance -- and spawn_seq is used rather than bank position because a
        // merge erases an element and shifts everything after it, which would
        // otherwise renumber templates between runs.

        inline uint8_t prematureCode(const std::string& word) {
            if (word == "premature") return 3;
            if (word == "vote")      return 2;
            if (word == "mixed")     return 1;
            return 0;
        }
        inline uint8_t tukeyAggCode(const std::string& word) {
            if (word == "removed") return 3;
            if (word == "partial") return 2;
            if (word == "kept")    return 1;
            return 0;
        }

    }  // namespace detail

    // ---- THE BEATS FILE IS THE BIG ONE, AND IT IS SWITCHABLE --------------
    //
    // ON BY DEFAULT. It is the only per-beat output there is: which group each
    // heartbeat landed in, whether it became a beat on each channel, and why it
    // left its group's average. Nothing else records any of that -- the other
    // two files are per template and per bin -- so with this off a partition
    // cannot be audited after the fact.
    //
    // It is also, by a wide margin, the largest thing the pipeline writes:
    // (slices x axis width) doubles per channel, four channels. Most of that is
    // NaN padding, because `width` is the widest beat in the whole block, so one
    // long RR anywhere in a channel pads every column in it. Per-bin width is
    // the fix and has not been done.
    //
    // So: a switch rather than a commented-out call. Set it false at the top of
    // a debug run and the writer returns immediately, leaving an empty file
    // rather than a stale one from a previous run -- a half-written file with
    // yesterday's partition in it is worse than no file.
    inline bool g_write_beats_bin = true;
    inline void setWriteBeatsBin(bool on) { g_write_beats_bin = on; }

    inline bool writeBeatsBin(const std::vector<ChannelBlock>& blocks) {
        if (!g_write_beats_bin) {
            // Truncate to nothing, so the stem's beats file is unmistakably
            // empty rather than left over from an earlier run.
            if (!g_dir.empty() && !g_stem.empty())
                std::ofstream(g_dir + "/" + g_stem + "_beats.bin",
                    std::ios::binary | std::ios::trunc);
            return true;
        }
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_beats.bin",
            std::ios::binary | std::ios::trunc);
        if (!f) return false;

        uint32_t nBlocks = 0;
        for (const auto& b : blocks) if (!b.empty()) ++nBlocks;

        detail::writeRaw(f, &bin_version, 4);
        detail::writeRaw(f, &nBlocks, 4);

        for (const auto& blk : blocks) {
            if (blk.empty()) continue;

            // ---- width, and why it is per BLOCK and not per beat -----------
            // Every column in a block is `width` doubles, so the block's width
            // is the widest beat in it and the rest are NaN-padded out to it.
            // NO BEAT IS EVER TRUNCATED to make a file smaller -- storage keeps
            // every sample, and the plot width is what gets clipped to the
            // median beat length instead.
            //
            // nColumns is the SLICE COUNT summed over bins, taken from the flag
            // vector because that is the one array in ChannelOutput sized to the
            // slice count by construction. A slice with no beat on this channel
            // still gets its column.
            uint32_t width = 0;
            uint64_t nCols = 0;
            for (size_t b = 0; b < blk.nBins(); ++b) {
                const bin_pipeline::ChannelOutput* out = blk.out(b);
                if (!out) continue;
                nCols += out->flags.size();
                if (const auto* bb = blk.binBeats(b))
                    for (const auto& bt : *bb)
                        width = std::max(width, static_cast<uint32_t>(bt.size()));
            }
            detail::writeBlockHeader(f, blk.channel, width, nCols);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            // ONE BUFFER, REUSED. Every sample used to be a separate 8-byte
            // ofstream::write; a record is (slices x width) of them, tens of
            // millions, and the syscall-per-double dominated the writer's time.
            // One write per column instead.
            std::vector<double> row;
            row.reserve(width);

            for (size_t b = 0; b < blk.nBins(); ++b) {
                const bin_pipeline::ChannelOutput* outp = blk.out(b);
                if (!outp) continue;
                const bin_pipeline::ChannelOutput& out = *outp;

                // Letters for this bin's bank, contiguous over the surviving
                // templates. Computed once per bin, not per column.
                const std::vector<uint8_t> letters = detail::letterRanks(out.bank);
                const std::vector<uint8_t>* why = blk.reasons(b);
                const std::vector<std::vector<double>>* binBeats = blk.binBeats(b);

                for (size_t slice = 0; slice < out.flags.size(); ++slice) {
                    BeatRecord rec;
                    rec.bin = static_cast<uint32_t>(b);
                    rec.category = static_cast<uint8_t>(out.flags[slice].category);
                    rec.premature = static_cast<uint8_t>(out.flags[slice].pvc);
                    rec.tukey = static_cast<uint8_t>(out.flags[slice].tukey);
                    rec.r_col = blk.rCol(b);
                    rec.excluded = (why && slice < why->size()) ? (*why)[slice] : 0u;

                    const int32_t t = (slice < out.assignment.size())
                        ? out.assignment[slice] : -1;
                    rec.template_id = t;
                    if (t >= 0 && t < out.bank.size()) {
                        const tbank::BankTemplate& tp = out.bank.templates[t];
                        rec.label_code = tp.label_code;
                        rec.letter = letters[t];
                        rec.confirmed = tp.confirmed() ? 1 : 0;
                    }
                    else {
                        rec.confirmed = 2;   // n/a: unassigned or unscorable
                    }

                    // ---- did this slice become a beat on THIS channel? -----
                    const int r = blk.rowOf(b, slice);
                    const std::vector<double>* bt =
                        (binBeats && r >= 0
                            && static_cast<size_t>(r) < binBeats->size())
                        ? &(*binBeats)[static_cast<size_t>(r)] : nullptr;
                    rec.became_beat = bt ? 1u : 0u;

                    detail::writeRaw(f, &rec, sizeof(rec));

                    row.assign(width, nan);
                    if (bt) {
                        const size_t n = std::min<size_t>(width, bt->size());
                        for (size_t k = 0; k < n; ++k) row[k] = (*bt)[k];
                    }
                    if (width) detail::writeRaw(f, row.data(),
                        static_cast<size_t>(width) * sizeof(double));
                }
            }
        }
        return static_cast<bool>(f);
    }

    inline bool writeTemplatesBin(const std::vector<ChannelBlock>& blocks) {
        if (g_dir.empty() || g_stem.empty()) return false;
        std::ofstream f(g_dir + "/" + g_stem + "_templates.bin",
            std::ios::binary | std::ios::trunc);
        if (!f) return false;

        uint32_t nBlocks = 0;
        for (const auto& b : blocks) if (!b.empty()) ++nBlocks;

        detail::writeRaw(f, &bin_version, 4);
        detail::writeRaw(f, &nBlocks, 4);

        for (const auto& blk : blocks) {
            if (blk.empty()) continue;

            const bool isPpgBin =
                (blk.channel && std::string(blk.channel) == "PPG");

            // COUNT AND WRITE MUST APPLY THE SAME FILTER. nCols goes in the
            // block header, so a template excluded from the payload but counted
            // here leaves the reader expecting one record more than exists --
            // it then reads the next block's header as a record and produces
            // plausible garbage rather than an error. Identical predicate,
            // identical order, both loops.
            uint32_t width = 0;
            uint64_t nCols = 0;
            for (const bin_pipeline::ChannelOutput* op : blk.per_bin) {
                if (!op) continue;
                const bin_pipeline::ChannelOutput& out = *op;
                for (int t = 0; t < out.bank.size(); ++t) {
                    const tbank::BankTemplate& tp = out.bank.templates[t];
                    if (tp.tmpl.empty()) continue;
                    if (!detail::belongsInTemplatesFile(tp, out, nullptr)) continue;
                    ++nCols;
                    width = std::max(width,
                        static_cast<uint32_t>(tp.tmpl.size()));
                }
            }
            detail::writeBlockHeader(f, blk.channel, width, nCols);

            const double nan = std::numeric_limits<double>::quiet_NaN();
            for (size_t b = 0; b < blk.nBins(); ++b) {
                const bin_pipeline::ChannelOutput* op = blk.out(b);
                if (!op) continue;
                const bin_pipeline::ChannelOutput& out = *op;
                // Letters for this bin's bank, contiguous over the surviving
                // templates. Computed once per bin, not per template.
                const std::vector<uint8_t> letters = detail::letterRanks(out.bank);
                for (int t = 0; t < out.bank.size(); ++t) {
                    const tbank::BankTemplate& tp = out.bank.templates[t];
                    if (tp.tmpl.empty()) continue;
                    if (!detail::belongsInTemplatesFile(tp, out, nullptr)) continue;

                    TemplateRecord rec;
                    rec.bin = static_cast<uint32_t>(b);
                    rec.category = static_cast<uint8_t>(tp.presumedCategory());
                    rec.premature = detail::prematureCode(detail::prematureWordAgg(tp));
                    rec.tukey = detail::tukeyAggCode(detail::tukeyWordAgg(tp, out));
                    rec.confirmed = tp.confirmed() ? 1 : 0;
                    rec.label_code = tp.label_code;
                    rec.letter = letters[t];
                    rec.landmark_marked = tp.wantsLandmarkMarking() ? 1 : 0;
                    rec.template_id = t;
                    rec.r_col = blk.rCol(b);
                    rec.n_members = static_cast<uint32_t>(tp.memberCount());
                    rec.n_excluded = static_cast<uint32_t>(tp.excludedCount());
                    rec.too_few_beats = tp.tooFewBeats(isPpgBin) ? 1u : 0u;
                    rec.beat_share = out.bank.beatShare(t);
                    detail::writeRaw(f, &rec, sizeof(rec));

                    for (uint32_t s = 0; s < width; ++s) {
                        const double v = (s < tp.tmpl.size()) ? tp.tmpl[s] : nan;
                        detail::writeRaw(f, &v, 8);
                    }

                    // ---- WHICH beats, not just how many ----------------
                    //
                    // BOTH LISTS. members is everything the partition assigned;
                    // members_clean is what tmpl was actually averaged over.
                    // Carrying only the first would make a reader that rebuilds
                    // the average include the premature and Tukey-rejected
                    // beats that cleaning removed; only the second would lose
                    // those beats from the record entirely. n_members and
                    // n_excluded above are the counts of these two, so the four
                    // agree by construction.
                    detail::writeVecU32(f, tp.members);
                    detail::writeVecU32(f, tp.members_clean);

                    detail::TemplateTrailer tr;
                    tr.subtype = tp.subtype;
                    tr.spawn_seq = tp.spawn_seq;
                    tr.n_ppg_members = tp.n_ppg_members;
                    tr.n_tukey_members = tp.n_tukey_members;
                    tr.n_blended_members = tp.n_blended_members;
                    tr.n_premature_members = tp.n_premature_members;
                    tr.n_voted_members = tp.n_voted_members;
                    tr.n_noise_members = tp.n_noise_members;
                    tr.mean_rr_ms = tp.mean_rr_ms;
                    tr.marked_invalid_template =
                        tp.marked_invalid_template ? 1u : 0u;
                    tr.operator_state = tp.operator_state;
                    tr.confirmed_by_operator =
                        tp.confirmed_by_operator ? 1u : 0u;
                    // Per bank, repeated per record. See TemplateTrailer.
                    tr.configured_cap = out.bank.configured_cap;
                    tr.effective_cap = out.bank.effective_cap;
                    tr.next_spawn_seq = out.bank.next_spawn_seq;
                    tr.assigned_beats = out.bank.assigned_beats;
                    detail::writeRaw(f, &tr, sizeof(tr));

                    // LENGTH-PREFIXED, not padded out to `width`: the spread is
                    // the same length as tmpl by construction, and padding it
                    // would put NaNs into a reloaded tmpl_iqr that were never
                    // measured.
                    detail::writeVecD(f, tp.tmpl_iqr);
                }
            }
        }
        return static_cast<bool>(f);
    }

    // ---------------------------------------------------------------------
    // Readers. Present so the format has exactly one definition rather than
    // one here and another in whatever consumes it -- a binary layout with a
    // writer and no reader is a layout nobody can verify.
    // ---------------------------------------------------------------------

    template <class Rec>
    struct BinBlock {
        std::string channel;
        uint32_t width = 0;
        std::vector<Rec> records;
        std::vector<double> samples;   // records.size() * width, column-major

        // ---- TEMPLATES ONLY, parallel to `records` ---------------------
        //
        // Empty on a beats file, because readBin only fills them when the
        // caller says the format carries them. So a consumer must check the
        // size rather than assume it matches `records` -- which is exactly what
        // bank_reload::readSplit does to decide whether a reload can be exact.
        //
        // members is the partition itself, in the channel's LOCAL ROW SPACE as
        // projectToChannel left it -- the same space _beats.bin's columns are
        // in, so the two files join on it without translation.
        std::vector<std::vector<uint32_t>> members;
        std::vector<std::vector<uint32_t>> members_clean;
        std::vector<detail::TemplateTrailer> trailers;
        std::vector<std::vector<double>> tmpl_iqr;

        const double* column(size_t k) const {
            return (k < records.size() && width) ? &samples[k * width] : nullptr;
        }
    };

    // hasExtras: this format appends the per-template member lists, trailer and
    // spread after each column's waveform. True only for the templates file. A
    // flag rather than a test on the record type, because the two are told
    // apart by which entry point was called and nothing else.
    template <class Rec>
    inline bool readBin(const std::string& path,
        std::vector<BinBlock<Rec>>& out, bool hasExtras = false)
    {
        std::ifstream f(path, std::ios::binary);
        if (!f) return false;
        uint32_t ver = 0, nBlocks = 0;
        if (!detail::readRaw(f, &ver, 4) || !detail::readRaw(f, &nBlocks, 4))
            return false;
        if (ver > bin_version) return false;   // newer than this build understands

        out.clear();
        for (uint32_t bi = 0; bi < nBlocks; ++bi) {
            BinBlock<Rec> blk;
            uint32_t len = 0;
            if (!detail::readRaw(f, &len, 4) || len > 64) return false;
            blk.channel.resize(len);
            if (len && !detail::readRaw(f, blk.channel.data(), len)) return false;
            uint64_t nCols = 0;
            if (!detail::readRaw(f, &blk.width, 4) || !detail::readRaw(f, &nCols, 8)) return false;

            blk.records.resize(nCols);
            blk.samples.assign(static_cast<size_t>(nCols) * blk.width, 0.0);
            if (hasExtras) {
                blk.members.resize(nCols);
                blk.members_clean.resize(nCols);
                blk.trailers.resize(nCols);
                blk.tmpl_iqr.resize(nCols);
            }
            for (uint64_t k = 0; k < nCols; ++k) {
                if (!detail::readRaw(f, &blk.records[k], sizeof(Rec))) return false;
                if (blk.width && !detail::readRaw(f,
                    &blk.samples[static_cast<size_t>(k) * blk.width],
                    static_cast<size_t>(blk.width) * 8)) return false;
                if (hasExtras) {
                    if (!detail::readVecU32(f, blk.members[k])) return false;
                    if (!detail::readVecU32(f, blk.members_clean[k])) return false;
                    if (!detail::readRaw(f, &blk.trailers[k],
                        sizeof(detail::TemplateTrailer))) return false;
                    if (!detail::readVecD(f, blk.tmpl_iqr[k])) return false;
                }
            }
            out.push_back(std::move(blk));
        }
        return true;
    }

    inline bool readBeatsBin(const std::string& path,
        std::vector<BinBlock<BeatRecord>>& out) {
        return readBin<BeatRecord>(path, out);
    }
    // MEMBER LISTS, TRAILER AND SPREAD INCLUDED. Returned in BinBlock::members
    // / ::members_clean / ::trailers / ::tmpl_iqr, parallel to ::records.
    inline bool readTemplatesBin(const std::string& path,
        std::vector<BinBlock<TemplateRecord>>& out) {
        return readBin<TemplateRecord>(path, out, /*hasExtras=*/true);
    }

}  // namespace morphology_csv