#pragma once
/**
 * @file   template_bank.hpp
 * @brief  Multi-template morphology segregation (Spec Section 4.6).
 *
 *        This file handles the bank of templates which each new beat is compared to.
 *        If it resembles an existing template, it is assigned to that template.
 *        If it does not resemble any existing template, a new template is spawned.
 *
 */

#include <algorithm>
#include <array>
#include "noise_marking_gui/annotation_types.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace tbank {
    namespace correlation_floors {
        inline double g_ecg = 0.0;   //these are set later by the config
        inline double g_ppg = 0.0;
    }

    inline double matchFloorEcg() { return correlation_floors::g_ecg; }
    // ---- BankTemplate::split_source encoding ------------------------------
    // Channel index PLUS ONE, so that zero stays free for "no split / not
    // recorded". Writing the raw channel index would make CH1 indistinguishable
    // from an unset field, and every seed template would read as "split by
    // CH1".
    inline constexpr uint8_t kSplitUnknown = 0;
    inline constexpr uint8_t kSplitCh1 = 1;
    inline constexpr uint8_t kSplitCh2 = 2;
    inline constexpr uint8_t kSplitCh3 = 3;
    inline constexpr uint8_t kSplitPpg = 4;


    inline const char* splitSourceLabel(uint8_t s) {
        switch (s) {
        case kSplitCh1:
        case kSplitCh2:
        case kSplitCh3: return "ECG";
        case kSplitPpg: return "PPG";
        default:        return nullptr;
        }
    }

    inline double matchFloorPpg() { return correlation_floors::g_ppg; }

    namespace minimum_beats {
        inline int g_ecg = 0.0; //these are also set later by the config
        inline int g_ppg = 0.0;
    }

    inline int minBeatsEcg() { return minimum_beats::g_ecg; }
    inline int minBeatsPpg() { return minimum_beats::g_ppg; }

    // Negative is refused and changes nothing.
    inline bool setMinBeats(int ecg, int ppg) {
        if (ecg < 0 || ppg < 0) return false;
        minimum_beats::g_ecg = ecg;
        minimum_beats::g_ppg = ppg;
        return true;
    }

    inline bool setMatchFloors(double ecg, double ppg) {
        if (!(ecg > 0.0 && ecg <= 1.0)) return false;
        if (!(ppg > 0.0 && ppg <= 1.0)) return false;
        correlation_floors::g_ecg = ecg;
        correlation_floors::g_ppg = ppg;
        return true;
    }

    inline constexpr int max_templates_per_bin = 6;//num templates in bank before merging

    // Minimum overlapping non-NaN columns for a correlation to mean anything.
    // alignment.hpp's pearson() returns 0.0 below its own floor, and 0.0 is
    // below every match floor, so an unscorable beat would spawn a template.
    // Beats below this get kUnscorable instead of an assignment.
    inline constexpr int kMinOverlapColumns = 8;

    inline constexpr int min_members_per_column = 2; //how many beats must be in a template to have a column at all

    // Merge-eligibility ceiling for the garbage tier, deliberately ABOVE
    // kMinMembersForColumn. Tying the two together is a trap: merging two
    // 1-member noise templates yields a 2-member template, which earns a column
    // and is therefore no longer garbage -- so junk gets promoted out of the
    // tier that exists to collect it, accumulates as 2-member zombies holding
    // slots, and total churn rises rather than falls. A ceiling of 3 lets a
    // merged pair remain collectable once more.
    inline constexpr int kMaxJunkMembers = 3;

    // Members below which a template's own 2.5/97.5 corridor is not an estimate
    // of anything and slot 0's spread is inherited instead (design note 1).
    // Four is the fewest that gives the 2.5 and 97.5 percentiles two distinct
    // values to interpolate between at all; below it they collapse onto the
    // min and max, which is a range, not a percentile.
    inline constexpr int kMinMembersForCorridor = 4;

    // Corridor half-width when even slot 0 has no usable spread at a column --
    // an all-NaN column, or a bin so sparse that slot 0 itself is degenerate.
    // Expressed as a fraction of the template's own peak-to-peak amplitude so
    // it carries no unit and survives normalization. Deliberately generous: a
    // corridor that is too wide admits a beat that should have spawned, which
    // pass 2 can still move; one that is too narrow spawns a template per beat,
    // which nothing downstream can undo.
    inline constexpr double kFallbackCorridorFrac = 0.15;

    inline double corridorInflation(int n_members) {
        if (n_members < 1) return 1.0;
        return std::sqrt(1.0 + 1.0 / static_cast<double>(n_members));
    }

    inline constexpr int kUnscorable = -2;   // returned by assignment
    inline constexpr int kNoMatch = -1;   // spawn required

    inline constexpr uint8_t kUnlabeled = 0;   // label_code sentinel

    // Class codes, READ FROM annotation_types.hpp rather than mirrored.
    //
    // These used to be literals duplicating that table, because the table
    // header pulled in QString and QColor and everything in this namespace --
    // bin_pipeline, seed_pool, nsvt_detect -- runs on the template-generation
    // side, which is Qt-free and should stay that way. The duplication was a
    // known drift risk and annotation_code_check.hpp existed to catch it at
    // startup. It was never called from anywhere, so for as long as the mirror
    // existed the drift was entirely unguarded.
    //
    // annotation_table.hpp is the table with the Qt helpers split off, so there
    // is nothing left to mirror: each constant is the table entry, resolved at
    // compile time by label. A renamed or removed label is a build failure on
    // the static_assert below rather than a 0 that silently becomes "unlabeled",
    // and there is no second copy for anyone to update out of step.
    //
    // Note that the Section 4.6 addendum asserts "PVC is annotation type 5",
    // while the table has 4) PVC at code 4 and 5) PAC at code 5. The table wins,
    // and now it wins structurally: taking the spec's number would mean writing
    // a literal here again.
    inline constexpr uint8_t kCodeMinorNoise = annotation_types::code_for_label("2) Minor Noise");
    inline constexpr uint8_t kCodePvc = annotation_types::code_for_label("4) PVC");
    inline constexpr uint8_t kCodePac = annotation_types::code_for_label("5) PAC");
    inline constexpr uint8_t kCodeCondDelay = annotation_types::code_for_label("6) Cond. Delay");
    inline constexpr uint8_t kCodeAf = annotation_types::code_for_label("7) AF");
    inline constexpr uint8_t kCodeSvt = annotation_types::code_for_label("8) SVT");
    inline constexpr uint8_t kCodeVt = annotation_types::code_for_label("9) VT");

    // The whole of what verifyAnnotationCodes() did, at compile time. 0 is the
    // "no such label" return and is also kUnlabeled, so a missing entry would
    // otherwise turn every template of that class into an unlabeled one.
    static_assert(kCodeMinorNoise != 0, "annotation_types: no \"2) Minor Noise\" entry");
    static_assert(kCodePvc != 0, "annotation_types: no \"4) PVC\" entry");
    static_assert(kCodePac != 0, "annotation_types: no \"5) PAC\" entry");
    static_assert(kCodeCondDelay != 0, "annotation_types: no \"6) Cond. Delay\" entry");
    static_assert(kCodeAf != 0, "annotation_types: no \"7) AF\" entry");
    static_assert(kCodeSvt != 0, "annotation_types: no \"8) SVT\" entry");
    static_assert(kCodeVt != 0, "annotation_types: no \"9) VT\" entry");
    static_assert(kCodePvc != kCodePac, "annotation_types: PVC and PAC share a code");

    // Ventricular origin, for NSVT run detection. PAC, AF and SVT are
    // supraventricular and must not qualify.
    inline bool isVentricular(uint8_t label_code) {
        return label_code == kCodePvc || label_code == kCodeVt;
    }

    //the only categories taken into account are regular - eventually there will be 5
    enum class Category : uint8_t {
        REGULAR = 1,   // unmarked; also 6) Cond. Delay, 7) AF, 8) SVT for now
        ECTOPIC = 2,   // marks 4) PVC, 5) PAC, 9) VT, and the postEligible
        //   beat following each of them
        NOISE = 3    // mark 2) Minor Noise
    };

    // An operator mark code -> the 4.5 category it implies. ONE definition, and
    // it is here because this is where the codes are named: bin_pipeline had its
    // own copy written as bare literals (case 4: case 5: case 9:), which is the
    // same table with nothing tying it to annotation_types, so a renumbered
    // annotation would silently reclassify beats there and not here.
    //
    // Absence of a mark is the REGULAR verdict, not an unknown one -- there is
    // no third state, which is why the default arm returns REGULAR rather than
    // failing.
    inline Category categoryForLabelCode(uint8_t code) {
        // NOISE covers "2) Minor Noise" AND any type flagged
        // suppressesDetection in annotation_types -- today that is "1) R Peak
        // Noise". Both mean the beat's morphology is corrupted and must not be
        // averaged into a template; keying the second condition off the table
        // flag (rather than hardcoding code 1) keeps the rule in one place.
        // Previously only Minor Noise was recognized here, so R-Peak-Noise
        // beats fell through to REGULAR and were averaged into the plotted
        // template despite being marked.
        if (code == kCodeMinorNoise
            || annotation_types::code_suppresses_detection(code))
            return Category::NOISE;
        if (code == kCodePvc || code == kCodePac || code == kCodeVt)
            return Category::ECTOPIC;
        return Category::REGULAR;
    }

    // Per-beat categories from a mark vector, WITH THE POST-ECTOPIC RULE.
    //
    // The beat FOLLOWING an ectopic one inherits ECTOPIC when it carries no mark
    // of its own. That beat is the compensatory pause -- its RR is long, its
    // baseline is still recovering -- and admitting it to the sinus reference
    // contaminates the template with post-ectopic morphology. It is a separate
    // rule from categoryForLabelCode and cannot be derived from a single mark,
    // which is why this takes the whole vector.
    //
    // A MARKED BEAT IS NEVER OVERWRITTEN: the operator's verdict on beat i wins
    // over what beat i-1 implies about it.
    //
    // Moved here from bin_pipeline::categoryFromMarks, which died with
    // runChannel. Its postEligibleCode() also listed AF and SVT, but they only
    // mattered in conjunction with an ECTOPIC verdict on the same mark and both
    // map to REGULAR -- so the set that can actually trigger inheritance is
    // PVC / PAC / VT, and that is what this tests. Same behaviour, one fewer
    // table.
    inline std::vector<Category> categoriesFromMarks(
        const std::vector<uint8_t>& mark_code)
    {
        const size_t n = mark_code.size();
        std::vector<Category> out(n, Category::REGULAR);
        for (size_t i = 0; i < n; ++i) {
            out[i] = categoryForLabelCode(mark_code[i]);
            if (i > 0 && mark_code[i] == 0
                && categoryForLabelCode(mark_code[i - 1]) == Category::ECTOPIC)
                out[i] = Category::ECTOPIC;
        }
        return out;
    }



    // The two paths are recorded separately because their competence is
    // disjoint. The raw flag catches isolated ectopy. The 5-of-8 vote only
    // rescues beats inside CONSECUTIVE runs -- in perfect bigeminy exactly 4
    // of any 8 beats are ectopic, the count never reaches 5, and the vote
    // never fires. High raw / zero vote is bigeminy; high vote is a run.
    enum class PvcFilter : uint8_t {
        NONE = 0,
        PREMATURE = 1,   // RR(t) < 0.80 * median of trailing ten
        VOTE = 2    // 5-of-8 over the prematurity flags
    };

    // Tukey runs INSIDE the classified category-1 population, so most beats
    // were never eligible for it. Recording which pass rejected a beat costs
    // nothing and pays for itself: a REGULAR beat rejected on RR length is a
    // strong hint of ectopy the classifier missed, and that count rising while
    // total rejections fall is the high-burden failure announcing itself.
    enum class TukeyOutcome : uint8_t {
        NOT_ELIGIBLE = 0,   // not category 1
        KEPT = 1,
        REJ_RR_LENGTH = 2,
        REJ_AMPLITUDE = 3,
        REJ_R_LOCATION = 4,
        REJ_WAVE_SCORE = 5
    };

    struct BeatFlags {
        Category     category = Category::REGULAR;
        PvcFilter    pvc = PvcFilter::NONE;
        TukeyOutcome tukey = TukeyOutcome::NOT_ELIGIBLE;

        // Bank assignment per ECG channel, and PPG. -1 = unassigned,
        // kUnscorable = too little axis overlap to score.
        std::array<int32_t, 3> template_id_ecg = { -1, -1, -1 };
        int32_t                template_id_ppg = -1;

        // True when this beat's waveform is an EWMA substitution rather than
        // an observation (4.6 beat substitution, alpha = 1/8, run AFTER
        // assignment against the assigned template's own average). Variance
        // measures must not treat these as observations.
        bool substituted = false;
    };

    // ---------------------------------------------------------------------
    // Per-bin counts. This defect class produces output that looks fine and
    // shows up only in the counts, so they are first-class state, not
    // debug logging.
    //
    // Tukey assumes one dominant population. With ectopy in the minority both
    // quartiles sit in the sinus cluster and every ectopic beat falls outside
    // the fence. As ectopy approaches half the bin the IQR widens across both
    // clusters and Tukey stops rejecting anything, including real artifact.
    // Moving Tukey inside the classified pool does not remove that failure --
    // it makes Tukey's correctness CONDITIONAL ON THE CLASSIFIER'S RECALL,
    // one layer deeper and quieter, with no backstop behind it. The two
    // mechanisms also degrade together rather than independently: at high
    // burden the trailing-ten median is dragged down by ectopic RRs, so
    // RR(t) < 0.80*median fires less often at exactly the burden where the
    // fences have widened past usefulness.
    //
    // Hence fence_iqr alongside the rejection counts. Counts give the
    // outcome; the IQR gives the cause, and it is continuous where the counts
    // are discrete. The signature to watch for is LOW REJECTION RATE WITH
    // WIDE IQR -- indistinguishable from "the bin is genuinely clean" if you
    // only have the counts.
    // ---------------------------------------------------------------------

    struct TukeyPassCounts {
        uint32_t beats_in = 0;
        uint32_t beats_out = 0;
        uint32_t rejected = 0;
        double   q1 = std::numeric_limits<double>::quiet_NaN();
        double   q3 = std::numeric_limits<double>::quiet_NaN();
        double   fence_lo = std::numeric_limits<double>::quiet_NaN();
        double   fence_hi = std::numeric_limits<double>::quiet_NaN();
        double   iqr() const { return q3 - q1; }
    };

    struct BinCounts {
        uint32_t beats_detected = 0;

        // Classification, before anything is removed.
        uint32_t n_regular = 0;
        uint32_t n_ectopic = 0;
        uint32_t n_noise = 0;

        // PVC filter over ALL beats, paths kept separate.
        uint32_t n_premature = 0;
        uint32_t n_vote_only = 0;   // voted but not itself premature

        // Tukey inside category 1, one entry per pass.
        TukeyPassCounts tukey_rr, tukey_amplitude, tukey_r_location,
            tukey_wave_score;

        // Category-1 beats Tukey rejected on RR length: the cheapest estimate
        // of classifier recall you will get. A sinus-labeled beat rejected for
        // being short is very likely ectopy that was missed.
        uint32_t n_regular_rejected_on_rr = 0;

        // Bank churn. Merges are mostly garbage collection of lone noise
        // templates, so this doubles as a noise metric: heavy merge activity
        // means a high noise fraction in the bin.
        uint32_t n_spawns = 0;
        uint32_t n_merges = 0;
        // Merges split by tier. n_merges_garbage collected two templates that
        // both sat below kMinMembersForColumn; n_merges_real absorbed at least
        // one template with real membership, which is the case worth watching --
        // a bin with many real merges lost morphologies to cap pressure and its
        // template count understates what was there.
        uint32_t n_merges_garbage = 0;
        uint32_t n_merges_real = 0;
        uint32_t n_cap_raises = 0;
        uint32_t n_unscorable = 0;
        uint32_t n_reassigned_pass2 = 0;

        // --- how slot 0 was obtained ------------------------------------
        // Exactly one of the first two is 1. seed_fallback_used == 1 means no
        // Phase 1 template reached this bin and slot 0 is a locally rebuilt
        // seed-pool median, so this bin does NOT satisfy Section 4.6's first
        // bullet and its features are measured against a different reference
        // from Phase 1's. seed_verified_sinus == 0 with seed_from_phase1 == 1
        // means the template was supplied but the caller could not vouch that
        // it excludes ectopy -- which, given create_ecg_templates.hpp applies no
        // rhythm test, is the expected state today.
        uint32_t seed_from_phase1 = 0;
        uint32_t seed_fallback_used = 0;
        uint32_t seed_verified_sinus = 0;

        // NOT SERIALIZED, deliberately. writeCounts() lays BinCounts out
        // field-by-field with writeSeed() immediately after it inside the bank
        // block, so appending here shifts everything downstream and needs a
        // kVersionBank -> v2 bump plus a version branch in readCounts(). These
        // three are build-time diagnostics and surface in the per-bin stderr
        // log. Persist them with the next format bump, not before one.
    };

    // ---------------------------------------------------------------------
    // One template in a bank
    // ---------------------------------------------------------------------

    struct BankMarkerSet {
        // THE FOUR DRAGGABLE BARS. Operator-owned: -1 means untouched, and an
        // untouched bar falls back to the auto value below.
        double p_begin = -1, q_onset = -1, s_end = -1, t_end = -1;

        // THE DETECTOR'S ANSWER, STORED. seed_bank_template already computed
        // every one of these and then threw all but the bars away, so each
        // reader re-ran the detection to recover them -- the widget, the focus
        // panel and the export each independently, and a landmark detected
        // several times is a landmark that can be several values.
        //
        // Same treatment BankPulseMarkerSet already gives the pulse glyphs
        // (onset_auto, peak_auto, ...), which is why the pulse side never had
        // the bar-vs-glyph disagreements the ECG side did.
        //
        // NOT SERIALIZED, like the pulse autos and for the same reason:
        // template_bank_serialize.hpp stores nothing derived, and
        // markers_by_anchor is deliberately absent from the file. These are
        // recomputed at load, so adding them changes no format.
        // FIXED LANDMARKS ONLY. No p_peak_auto and no t_peak_auto: those two
        // are reactive glyphs, bracketed by the bars and recomputed by
        // FeatureMarks::reactive_ecg at every read. p_peak was deliberately
        // removed from this struct once before, because a stored copy drifted
        // from the X on screen as soon as a bracket bar moved; putting it back
        // would undo that.
        double p_begin_auto = -1;
        double q_onset_auto = -1, q_peak_auto = -1;
        double r_peak_auto = -1;
        double s_peak_auto = -1, s_end_auto = -1;
        double t_end_auto = -1;
        bool   q_onset_found_auto = false;

        // WHETHER THE DETECTOR HAS RUN, which is not the same question as
        // "does this set hold a bar". hasDetectedMarks used to answer the
        // second: an alignment owning no bar (R and J, once the J-point and
        // T-end bars moved onto Q) seeded an all -1 set, reported "never
        // seeded" forever, and re-detected on EVERY display -- R being the
        // alignment the operator starts on. seed_bank_template sets this
        // whether or not it had a bar to store.
        bool   seeded = false;

        bool isUnset() const {
            return p_begin < 0 && q_onset < 0 && s_end < 0 && t_end < 0;
        }
    };
    struct BankPulseMarkerSet {
        // The three bars: the pulse markers markerAtX will hand out.
        double onset = -1, dicrotic = -1, end = -1;

        // Frozen detector output, straight from FeatureMarks::PpgFiducials.
        // NOT derived from the bars and NOT persisted -- recomputed each load,
        // and read by BinPlotWidget::overridePulseGlyphs to paint the glyphs.
        double onset_auto = -1.0, peak_auto = -1.0, dicrotic_auto = -1.0;
        double peak2_auto = -1.0, end_auto = -1.0;
        bool notch_found = false;

        bool isUnset() const {
            return onset < 0 && dicrotic < 0 && end < 0;
        }
    };

    struct BankTemplate {
        // Column-wise NaN-skipping median over members, on the bin's shared
        // axis. Recomputed whenever membership changes (design note 3).
        std::vector<double> tmpl;
        std::vector<double> tmpl_iqr;      // per-sample spread
        int                 r_col = -1;
        // +-samples around r_col that the morphology-split correlation is
        // restricted to (0.5 s at the channel rate). <= 0 means unrestricted.
        // Set by recomputeGroupChannel from the channel's ChannelBeats; not
        // serialized (a reloaded bank scores full-width, which is the old
        // behaviour).
        int                 corr_halfwin = -1;
        uint8_t label_code = kUnlabeled;
        bool confirmed_by_operator = false;
        int32_t subtype = -1; //PVC_2, etc

        // Spawn order within the bin, assigned in pass 1 and NOT recomputed
        // in pass 2. This is what "order of first appearance" resolves
        // against, so that reassignment cannot renumber subtypes between
        // runs.
        uint32_t spawn_seq = 0;

        // PPG members of THIS template, when the template is a projection of a
        // joint (four-channel) group. -1 = unknown, which is what a
        // per-channel-built bank leaves it as.
        //
        // Exists because the pulse cohort is a property of the GROUP, not of the
        // bin. Every column of a bin used to report the same bin-wide
        // ppg_n_beats, which made a split look like it had ignored PPG
        // entirely. A group's PPG members are the members it actually has on
        // that channel, and they differ between siblings.
        int32_t n_ppg_members = -1;

        int32_t n_tukey_members = 0;// how many members were rejected due to tukey? this is the FINAL rejection step after morphology split and premature/voting split
        uint32_t n_blended_members = 0;//how many members were blended with the ones beside them due to being PVC or voted PVC
        double mean_rr_ms = 0.0;// mean R-R over member slices, ms; 60000/this = bpm

        // ---- WHICH CHANNEL SPLIT THIS TEMPLATE OFF ----------------------
        //
        // The partition is joint over three ECG leads and the pulse, so a new
        // template is spawned when the best existing group is rejected -- and
        // ONE channel is the one that rejected it. That channel is the answer
        // to "is this bin splitting on cardiac morphology or on pulse noise",
        // which was previously only available as a per-bin tally
        // (jbank::BankCounts::n_rejected_by) and so could not be shown beside
        // the template it explains.
        //
        // 0 IS "NOT A SPLIT", which covers both the seed (spawn_seq 0, nothing
        // was rejected to create it) and an archive written before the field
        // existed. Readers must not render 0 as a channel; use
        // tbank::splitSourceLabel, which returns nullptr for it.
        uint8_t split_source = kSplitUnknown;

        bool marked_invalid_template = false;// operator/pipeline flagged this template as invalid

        // ---- THE OPERATOR'S QUALITY VERDICT ON THIS PANEL -----------------
        //
        // 0 = good, 1 = bad R detection, 2 = bad pulse - determined by right click
        uint8_t operator_state = 0;

        std::vector<uint32_t> members;//the templates, including the ones excluded by the morphology split
        std::vector<uint32_t> members_clean;//the templates, excluded the ones excluded by the morphology split

        int cleanCount() const {
            return static_cast<int>(members_clean.empty()
                ? members.size() : members_clean.size());
        }
        int excludedCount() const {
            return static_cast<int>(members_clean.empty()
                ? 0 : members.size() - members_clean.size());
        }

        //the minimum beats to be a template is stored in the config, i have it as 8 right now but it can change
        //(I know this comment is terrible and begging to be stale but i don't think i'm changing it from 8)
        bool tooFewBeats(bool is_ppg) const {
            const int lim = is_ppg ? minBeatsPpg() : minBeatsEcg();
            return lim > 0 && cleanCount() < lim;
        }

        // Landmarks per alignment anchor. A ventricular template has no P
        // wave, so p_begin/p_peak stay -1 legitimately and every P-dependent
        // feature must come out NaN rather than 0. Downstream extraction has
        // to treat -1 here as a valid state, not an error, or it will report
        // a PR interval measured from a P wave that does not exist.
        std::map<int32_t, BankMarkerSet> markers_by_anchor;   // key: AnchorType
        BankPulseMarkerSet pulse_marks;   // meaningful only on ppg_bank slots
        bool hasDetectedPulseMarks() const { return !pulse_marks.isUnset(); }

        // NOTE: INSERTS. markers_by_anchor[a] default-constructs an all -1 set
        // when the key is absent, so a read through this overload is a write.
        // Use hasDetectedMarks() to ask whether a set is populated.
        BankMarkerSet& marks(int32_t a) { return markers_by_anchor[a]; }

        // Non-inserting, and asks about CONTENT rather than key presence.
        // ASKS WHETHER THE DETECTOR RAN, not whether a bar survived it. The
        // !isUnset() form re-ran the detection forever on any alignment that
        // owns no bar; see BankMarkerSet::seeded.
        bool hasDetectedMarks(int32_t a) const {
            auto it = markers_by_anchor.find(a);
            return it != markers_by_anchor.end() && it->second.seeded;
        }
        const BankMarkerSet& marks(int32_t a) const {
            static const BankMarkerSet kEmpty;
            auto it = markers_by_anchor.find(a);
            return (it == markers_by_anchor.end()) ? kEmpty : it->second;
        }

        bool confirmed() const { return confirmed_by_operator; }

        // Per-member census, filled after assignment settles. These are what
        // make a template's PRESUMED category computable before any operator
        // mark exists: prematurity is a timing verdict available at build time,
        // and non-reproducibility shows up as a member count that never grew.
        uint32_t n_premature_members = 0;
        uint32_t n_voted_members = 0;
        uint32_t n_noise_members = 0;

        // Presumed, NOT confirmed. Decides which templates go in front of an
        // operator for landmark marking: only category 1 templates get
        // landmarks, because only category 1 beats feed feature extraction -- a
        // P-onset on a PVC template has nothing downstream that consumes it,
        // and a PVC's QT is not comparable to a sinus QT.
        //
        // UNLABELED IS PQRST, NOT A THIRD STATE. Section 4.6: "A template with
        // no confirmed member stays unlabeled. Unlabeled means not yet
        // confirmed, NOT unknown class." That second sentence forbids exactly
        // what an UNCONFIRMED category would be -- an unknown bucket that gets
        // withheld from display. A template is displayed as PQRST until an
        // operator marks it otherwise; being unconfirmed is a statement about
        // the operator's progress, not about the morphology.
        //
        // What the clause DOES forbid is inferring a label from morphological
        // similarity to a labeled template, and that prohibition is enforced
        // where the inference actually happened: mergeTemplates() no longer
        // copies a confirmed label onto the template it absorbs, and
        // findMergePair() blocks any pair containing a confirmed template. That
        // is the operation the spec names. This function is not it.
        //
        // A confirmed label overrides the presumption outright; the operator's
        // verdict is not a hypothesis to be re-derived.
        Category presumedCategory() const {
            if (confirmed_by_operator) {
                if (label_code == kCodeMinorNoise
                    || annotation_types::code_suppresses_detection(label_code))
                    return Category::NOISE;
                if (label_code == kCodePvc || label_code == kCodePac
                    || label_code == kCodeVt) return Category::ECTOPIC;
                return Category::REGULAR;
            }
            return Category::REGULAR; //unless marked PVC, PAC, VT or noise, presume regular
        }

        // The Phase 1 sinus seed. spawn_seq 0 is issued once, to slot 0, by the
        // pipeline before pass 1 begins, and mergeTemplates() keeps the LOWER
        // spawn_seq of a merged pair, so this survives merge history without a
        // stored flag -- which also keeps the serialized layout byte-identical.
        bool isSeed() const { return spawn_seq == 0; }

        bool wantsLandmarkMarking() const {
            return presumedCategory() == Category::REGULAR;
        }
        // Per-column 2.5/97.5 corridor over the members, and the flag saying
        // whether it is the template's OWN spread or slot 0's inherited one.
        // Derived from members exactly as tmpl is, so neither is serialized --
        // template_bank_serialize.hpp stores nothing derived, on purpose.
        std::vector<double> band_lo;
        std::vector<double> band_hi;
        bool corridor_inherited = false;

        int  memberCount() const { return static_cast<int>(members.size()); }
        bool earnsColumn() const { return memberCount() >= min_members_per_column; }

        // Eligible for the garbage-collection merge tier. Not the negation of
        // earnsColumn() -- see kMaxJunkMembers.
        bool isJunk() const { return memberCount() <= kMaxJunkMembers; }
    };

    // ---------------------------------------------------------------------
    // The bank: one per channel per bin
    // ---------------------------------------------------------------------

    struct TemplateBank {
        std::vector<BankTemplate> templates;

        // Per-bin, because the confirmed-member rule raises it. The addendum's
        // assignToTemplate() takes maxTemplates by value and returns a single
        // int, so it can neither persist a raise nor report one for logging;
        // the cap lives here instead.
        // max_templates_per_bin. A CONFIG VALUE with a default, which is what
        // the spec says it is: "Cap the bank at max_templates_per_bin (default
        // 6)". It was previously reachable only by editing
        // kDefaultMaxTemplatesPerBin and recompiling, so the "default" was the
        // only value the program had. Set from cfg at pipeline entry
        // (BinInput::max_templates_per_bin); the constant is now only the
        // fallback when nothing supplies one.
        //
        // configured_cap is the value the operator asked for and never changes.
        // effective_cap starts equal to it and only ever RISES, by the
        // confirmed-member rule. Keeping both means a raise is visible as a
        // difference rather than having to be reconstructed from the log.
        int32_t configured_cap = max_templates_per_bin;
        int32_t effective_cap = max_templates_per_bin;

        void setCap(int32_t cap) {
            if (cap < 1) cap = 1;   // slot 0 always exists
            configured_cap = cap;
            effective_cap = cap;
        }

        uint32_t next_spawn_seq = 0;

        // Total beats routed into this bank, for beat-share reporting.
        uint32_t assigned_beats = 0;

        double matchFloorFor(bool is_ppg) const {
            return is_ppg ? matchFloorPpg() : matchFloorEcg();
        }

        int  size() const { return static_cast<int>(templates.size()); }
        bool atCap() const { return size() >= effective_cap; }

        int nConfirmed() const {
            int c = 0;
            for (const auto& t : templates) if (t.confirmed()) ++c;
            return c;
        }

        // Distinct templates confirmed as `code`, i.e. the highest subtype
        // index reached for that class. One PVC template is monomorphic; two
        // or more is polymorphic.
        int countLabeled(uint8_t code) const {
            int c = 0;
            for (const auto& t : templates) if (t.label_code == code) ++c;
            return c;
        }

        // Next subtype index for a class, by order of first appearance. Reads
        // the highest index already issued rather than counting current
        // members, so a merge that removes a labeled template cannot cause
        // the next one to reuse a retired index.
        int32_t nextSubtypeFor(uint8_t code) const {
            int32_t hi = 0;
            for (const auto& t : templates)
                if (t.label_code == code) hi = std::max(hi, t.subtype);
            return hi + 1;
        }

        int findByBeat(uint32_t beat_idx) const {
            for (int i = 0; i < size(); ++i)
                for (uint32_t m : templates[i].members)
                    if (m == beat_idx) return i;
            return -1;
        }

        double beatShare(int template_idx) const {
            if (template_idx < 0 || template_idx >= size() || assigned_beats == 0)
                return std::numeric_limits<double>::quiet_NaN();
            return static_cast<double>(templates[template_idx].memberCount())
                / static_cast<double>(assigned_beats);
        }
    };

    // ---------------------------------------------------------------------
    // Cap-raise events. "Log the event" is a requirement, not bookkeeping: a
    // bin that raised its cap has three or more confirmed morphologies
    // competing for six slots, which makes it the most clinically interesting
    // bin in the record and the one most likely to expose a bug.
    // ---------------------------------------------------------------------

    struct CapRaiseEvent {
        uint64_t bin_index = 0;
        int      channel = 0;
        int32_t  old_cap = 0;
        int32_t  new_cap = 0;
        int32_t  template_a = -1;   // the two closest, both confirmed
        int32_t  template_b = -1;
        double   closeness = 0.0;   // r between them
        uint8_t  label_a = kUnlabeled;
        uint8_t  label_b = kUnlabeled;
    };

    // (PolymorphicVerdict / polymorphicVerdict REMOVED. It took the max over
    //  three channels and reported which one drove it, because three
    //  per-channel partitions could disagree about how many morphologies a bin
    //  held. There is one partition now, so the count is just a count and
    //  "which channel drove it" is not a question that exists;
    //  jbank::polymorphyVerdict counts groups and is what the pipeline calls.)

    // ---------------------------------------------------------------------
    // Label propagation
    //
    // Operator confirms ONE beat -> the class label attaches to the template
    // that beat is assigned to -> and from there to every other beat in that
    // template. One click labels a whole morphology cluster.
    //
    // KEYED ON THE SLOT. Under the joint partition template i is group i on
    // every channel, so the slot is the shared handle and there is one
    // function below rather than a beat-keyed family beside it.
    // ---------------------------------------------------------------------

    struct PropagationResult {
        std::array<int32_t, 3> labeled_template = { -1, -1, -1 };
        std::array<int32_t, 3> subtype = { -1, -1, -1 };
        int beats_relabeled = 0;

        int32_t ppg_labeled_template = -1;
        int32_t ppg_subtype = -1;
        int     ppg_beats_relabeled = 0;
    };

    inline PropagationResult propagateLabelBySlot(std::array<TemplateBank, 3>& ecg_banks, TemplateBank& ppg_bank, int slot, uint8_t label_code)
    {
        PropagationResult out;
        if (label_code == kUnlabeled || slot < 0) return out;

        auto label = [&](TemplateBank& bk, int32_t* which, int32_t* sub,
            int* relabeled) {
                if (slot >= bk.size()) return;
                BankTemplate& t = bk.templates[slot];
                // Subtype is issued once and then immutable: a re-confirmation,
                // or a confirmation of another beat in an already-labeled
                // template, must not mint a second index for one class.
                if (t.subtype < 0 || t.label_code != label_code) {
                    if (t.label_code != label_code)
                        t.subtype = bk.nextSubtypeFor(label_code);
                    t.label_code = label_code;
                }
                t.confirmed_by_operator = true;
                if (which) *which = slot;
                if (sub) *sub = t.subtype;
                if (relabeled) *relabeled += t.memberCount();
            };

        for (int c = 0; c < 3; ++c)
            label(ecg_banks[c], &out.labeled_template[c], &out.subtype[c],
                &out.beats_relabeled);
        label(ppg_bank, &out.ppg_labeled_template, &out.ppg_subtype,
            &out.ppg_beats_relabeled);
        return out;
    }

    inline std::vector<uint8_t> letterRanks(const TemplateBank& bank) {
        const int n = bank.size();
        std::vector<int> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i)
            if (!bank.templates[i].tmpl.empty()) order.push_back(i);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return bank.templates[a].spawn_seq < bank.templates[b].spawn_seq;
            });

        std::vector<uint8_t> letter(n, 0);
        int rank = 0;
        for (int i : order) {
            const BankTemplate& t = bank.templates[i];
            const int idx = (t.label_code != kUnlabeled && t.subtype > 0)
                ? t.subtype - 1 : rank;
            letter[i] = static_cast<uint8_t>(idx % 26);
            ++rank;
        }
        return letter;
    }

}  // namespace tbank