#pragma once
/**
 * @file   bank_structs.hpp
 * @brief  In order to split a temporal bin into a morpohologically split template, it is neccecary for there to be a 'bank'
 *         of templates which each new split is compared to. This file controls the structure for this bank
 */

#include <algorithm>
#include <array>
#include "noise_marking_gui/annotation_types.hpp"
#include "stats_utils.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace tbank {
    // ---- MORPHOLOGY SPLIT THRESHOLDS -------------------------------------
    //
    // The band-match score a beat must reach against a group's template to
    // join it. Separate per side because the two channels are not comparable:
    // an ECG lead is far more repeatable beat to beat than a pulse-ox trace,
    // so one number that suits both suits neither.
    //
    // COMPILED DEFAULTS, NOT ZERO. These were 0.0 with the config expected to
    // fill them, and 0.0 is not an inert value -- it is a threshold every beat
    // clears, so a missing config cell put every beat in one group and no bin
    // was ever split. The symptom was a partition that silently did nothing.
    // 0.85 and 0.80 are the Section 4.6 values; a config cell overrides them.
    namespace correlation_floors {
        inline double g_ecg = 0.85;
        inline double g_ppg = 0.80;
    }

    inline double morphThresholdEcg() { return correlation_floors::g_ecg; }
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

    inline double morphThresholdPpg() { return correlation_floors::g_ppg; }

    namespace minimum_beats {
        inline int g_ecg = 0; //these are also set later by the config
        inline int g_ppg = 0;
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

    // ---- SET INDEPENDENTLY, ONE SIDE AT A TIME ---------------------------
    //
    // WAS ONE setMatchFloors(ecg, ppg) TAKING BOTH OR NEITHER: it returned
    // false and set NOTHING if either value was outside (0, 1], so a config
    // with a good ECG threshold and a blank PPG cell lost both. That coupling
    // made the two thresholds one setting with two numbers in it.
    //
    // Each returns false and changes nothing when its own value is outside
    // (0, 1]; the other side is unaffected either way.
    // ---- THE BEAT-AFTER-ECTOPIC SWITCH (config exclude_beat_after_ectopic) --
    namespace post_ectopic_rule { inline bool g_on = true; }
    inline void setExcludePostEctopic(bool on) { post_ectopic_rule::g_on = on; }
    inline bool excludePostEctopic() { return post_ectopic_rule::g_on; }

    inline bool setMorphThresholdEcg(double r) {
        if (!(r > 0.0 && r <= 1.0)) return false;
        correlation_floors::g_ecg = r;
        return true;
    }

    inline bool setMorphThresholdPpg(double r) {
        if (!(r > 0.0 && r <= 1.0)) return false;
        correlation_floors::g_ppg = r;
        return true;
    }

    inline constexpr int max_templates_per_bin = 6;//num templates in bank before merging
    inline constexpr int kMinOverlapColumns = 8;
    inline constexpr int kMaxJunkMembers = 3;
    inline constexpr int kMinMembersForCorridor = 4;
    inline constexpr double kFallbackCorridorFrac = 0.15;

    inline double corridorInflation(int n_members) {
        if (n_members < 1) return 1.0;
        return std::sqrt(1.0 + 1.0 / static_cast<double>(n_members));
    }

    inline constexpr int kUnscorable = -2;   // returned by assignment
    inline constexpr int kNoMatch = -1;   // spawn required

    inline constexpr uint8_t kUnlabeled = 0;   // label_code sentinel

    // ---- WHERE A TEMPLATE'S LABEL CAME FROM ------------------------------
    //
    // A label_code alone cannot say this, and two DIFFERENT operator inputs
    // can produce one:
    //
    //   kLabelFromPartition -- the group was built only from beats the
    //       operator marked in the NOISE-MARKING stage. The mark reaches here
    //       as the partition key, so every member carries that class by
    //       construction: cross-class assignment is structurally impossible.
    //   kLabelFromOperator  -- the operator confirmed the morphology in the
    //       TEMPLATE VIEWER, via propagateLabelBySlot.
    //
    // WHY THE DISTINCTION IS KEPT RATHER THAN COLLAPSED. Both are the
    // operator's judgement, so both are legitimate grounds for a label -- but
    // they are not the same claim. A noise-stage PVC mark says "this beat is
    // ectopic"; a viewer confirmation says "this is a distinct morphology in
    // this bin". polymorphyVerdict wants the second and only the second, and
    // it reads confirmed_by_operator, which stays false for a partition label.
    // nsvt_detect wants either, because it asks whether the template is
    // ventricular, not whether the bin is polymorphic.
    inline constexpr uint8_t kLabelNone = 0;
    inline constexpr uint8_t kLabelFromPartition = 1;
    inline constexpr uint8_t kLabelFromOperator = 2;

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
        ECTOPIC = 2,   // marks 4) PVC, 5) PAC, 9) VT. NOT the beat after one:
        //   that beat is sinus and stays REGULAR; see BeatFlags::post_ectopic
        NOISE = 3    // mark 2) Minor Noise
    };

    // An operator mark code -> the 4.5 category it implies. ONE definition, and
    // it is here because this is where the codes are named: a downstream copy had its
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
        if (code == kCodeMinorNoise || annotation_types::code_suppresses_detection(code))
            return Category::NOISE;
        if (code == kCodePvc || code == kCodePac || code == kCodeVt)
            return Category::ECTOPIC;
        return Category::REGULAR;
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

        // THE BEAT AFTER AN ECTOPIC ONE (baseqt_gissipheno.m's beat-after-PVC
        // rule): the previous beat was marked PVC / PAC / VT or flagged
        // PREMATURE, and this one is otherwise good -- unmarked, not premature,
        // not voted. A sinus beat, so it stays in the sinus partition and its RR
        // stays in the series; but its repolarisation is disturbed by the
        // preceding compensatory pause, so it is left out of its template's
        // average (ExcludeReason::POST_ECTOPIC) and so out of QT. Set only when
        // excludePostEctopic() is on.
        bool post_ectopic = false;
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

        // Bank churn. Merges are mostly garbage collection of lone noise
        // templates, so this doubles as a noise metric: heavy merge activity
        // means a high noise fraction in the bin.
        uint32_t n_spawns = 0;
        uint32_t n_merges = 0;
        // Merges split by tier. n_merges_garbage collected two templates that
        // both sat in the junk tier; n_merges_real absorbed at least
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

        // The height the Foot variant's rows were levelled at, 0 at the foot to
        // 100 at the apex. PERSISTED, because Auto resolves it per column from
        // a threshold and re-resolving it on reload rebuilds _F at a different
        // height. -1 means not resolved yet: decide it, then keep it.
        double foot_pct = -1.0;

        bool isUnset() const {
            return onset < 0 && dicrotic < 0 && end < 0;
        }
    };

    // ---- THE TWO PULSE ALIGNMENTS ----------------------------------------
    //
    // A pulse average exists in two versions of itself, the way an ECG average
    // exists once per AnchorType. Both are the column-wise median of the SAME
    // member beats over the SAME time axis; they differ only in what each beat
    // was levelled on before the median was taken:
    //
    //   Foot  each beat's own trough (or pct up its own upstroke), located
    //         from the operator's foot column. Rebuilt when that column or the
    //         percentage changes -- it is the operator's alignment.
    //   Peak  each beat's own systolic apex. No bar is among its inputs, so
    //         there is exactly one answer per slot, for the life of the record.
    //
    // WHICH LANDMARK IS MEASURED ON WHICH. The foot is measured on Foot: it
    // defines that alignment, so that is where it is sharp. The dicrotic notch
    // and the end of cycle are measured on Peak: both sit on the decay, which
    // a foot-levelled stack smears precisely because the beats were brought
    // together at the other end of the pulse.
    //
    // VERTICAL ONLY, which is what makes the pair usable together.
    // ppg_realign::relevelAtOwnCrossing adds a per-row offset and moves no row
    // sideways, so both variants share one width and one time axis and a column
    // index means the same instant in each. A foot read off Foot and a notch
    // read off Peak are therefore directly comparable in time; only their
    // AMPLITUDES belong to different waveforms, so a y value must always be
    // read from the variant its column was measured on.
    enum class PulseAnchor : int32_t { Foot = 0, Peak = 1 };

    inline const char* pulseAnchorSuffix(PulseAnchor v) {
        return (v == PulseAnchor::Peak) ? "_P" : "_F";
    }
    inline constexpr std::array<PulseAnchor, 2> pulse_anchor_array = {
        PulseAnchor::Foot, PulseAnchor::Peak
    };

    struct PulseVariant {
        std::vector<double> tmpl;
        std::vector<double> tmpl_iqr;
        BankPulseMarkerSet  marks;
        // FALSE UNTIL THE MEDIAN HAS ACTUALLY BEEN TAKEN. An empty tmpl is not
        // the same state: a slot whose re-level was refused (too few rows
        // reaching the columns, no beat matrix) has to be distinguishable from
        // one nobody has asked for yet, or every paint retries the refusal.
        bool built = false;
        bool ok() const { return built && !tmpl.empty(); }

        // ---- WHAT THE BUILD DID TO EACH ROW, for <id>_pulse_moves.csv ----
        //
        // Parallel arrays over the cohort rows that were offered (row_ids are
        // local rows of per_channel_beats["PPG"][bin]); NaN where the re-level
        // skipped the row. See ppg_realign::RelevelResult for the meaning of
        // each. Set by every build, cleared by a refused one; NOT serialized.
        double pct = -1.0;                    // % up the upstroke it levelled at
        std::vector<uint32_t> row_ids;
        std::vector<double>   row_anchor_col; // samples, in the beat frame
        std::vector<double>   row_v_shift;    // amount subtracted from the row
    };

    struct template_of_all_signals {
        // Column-wise NaN-skipping median over members, on the bin's shared
        // axis. Recomputed whenever membership changes (design note 3).
        std::vector<double> tmpl;
        std::vector<double> tmpl_std;      // per-sample spread
        int                 r_col = -1;
        // +-samples around r_col that the morphology-split correlation is
        // restricted to (0.5 s at the channel rate). <= 0 means unrestricted.
        // Set by recomputeGroupChannel from the channel's ChannelBeats; not
        // serialized (a reloaded bank scores full-width, which is the old
        // behaviour).
        int                 corr_halfwin = -1;
        uint8_t label_code = kUnlabeled;
        // kLabelNone / kLabelFromPartition / kLabelFromOperator. See the note
        // beside those constants: label_code alone does not say which operator
        // input produced it, and the two are different claims.
        uint8_t label_source = kLabelNone;
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
        // A BITFIELD, independent bits: the right-click cycle has a both-bad
        // step and the markings CSV reports bad_ecg and bad_ppg as two columns.
        //
        // The ECG bit is per lead, so it lives on ecg_bank[lead].templates[slot];
        // the pulse bit is shared across a panel's leads and lives on
        // ppg_bank.templates[slot]. Both banks serialize this field.
        //
        // PER TEMPLATE, NO SPECIAL CASE FOR ANY SLOT INDEX. TemplateBin's
        // bin-shaped bad_r_ch / bad_ppg are the pipeline's; they seed these
        // bits once at load and are not read for display or export after.
        static constexpr uint8_t kOperatorGood = 0u;
        static constexpr uint8_t kOperatorBadEcg = 1u << 0;
        static constexpr uint8_t kOperatorBadPulse = 1u << 1;
        uint8_t operator_state = kOperatorGood;

        bool badEcgMarked()   const { return (operator_state & kOperatorBadEcg) != 0; }
        bool badPulseMarked() const { return (operator_state & kOperatorBadPulse) != 0; }
        void setBadEcg(bool bad) {
            if (bad) operator_state |= kOperatorBadEcg;
            else     operator_state &= static_cast<uint8_t>(~kOperatorBadEcg);
        }
        void setBadPulse(bool bad) {
            if (bad) operator_state |= kOperatorBadPulse;
            else     operator_state &= static_cast<uint8_t>(~kOperatorBadPulse);
        }

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
        // THE COMPOSED SET THE PANEL DRAWS, and the one every existing reader
        // already knows about: the foot from the Foot variant, the dicrotic
        // notch and the end from Peak, the derived crossings from those. It is
        // written by composePulseMarks and is not itself a source of truth --
        // pulse_by_variant is.
        BankPulseMarkerSet pulse_marks;   // meaningful only on ppg_bank slots

        // ---- THE TWO VARIANTS, KEYED LIKE markers_by_anchor ---------------
        //
        // NOT SERIALIZED IN THE BANK FILE. Both are derivable from the beat
        // matrix and the foot column, so an archive carries neither; they are
        // rebuilt on demand. What DOES have to persist is each variant's
        // operator marks, which travel in the markings file.
        std::map<int32_t, PulseVariant> pulse_by_variant;

        PulseVariant& pulseVariant(PulseAnchor v) {
            return pulse_by_variant[static_cast<int32_t>(v)];
        }
        // NON-INSERTING, unlike the overload above -- a read through
        // operator[] is a write, the same trap marks() documents.
        const PulseVariant& pulseVariant(PulseAnchor v) const {
            static const PulseVariant kEmpty;
            const auto it = pulse_by_variant.find(static_cast<int32_t>(v));
            return (it == pulse_by_variant.end()) ? kEmpty : it->second;
        }
        bool hasPulseVariant(PulseAnchor v) const {
            return pulseVariant(v).ok();
        }

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

        // Eligible for the garbage-collection merge tier. Deliberately not the
        // complement of "has enough beats to keep" -- see kMaxJunkMembers. The
        // ceiling sits ABOVE the keep-threshold on purpose: merging two
        // 1-member noise templates yields a 2-member one, and if that were
        // instantly no longer junk, junk would be promoted out of the tier
        // that exists to collect it and accumulate as zombies holding slots.
        bool isJunk() const { return memberCount() <= kMaxJunkMembers; }
    };

    // ---------------------------------------------------------------------
    // The bank: one per channel per bin
    // ---------------------------------------------------------------------

    struct TemplateBank {
        std::vector<template_of_all_signals> templates;


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


        int  size() const { return static_cast<int>(templates.size()); }
        bool atCap() const { return size() >= effective_cap; }

        int nConfirmed() const {
            int c = 0;
            for (const auto& t : templates) if (t.confirmed()) ++c;
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

    // STAMPS kLabelFromOperator. This is the viewer-confirmation path, so it
    // is the stronger of the two provenances and overwrites a partition label
    // on the same group -- the operator looked at the waveform and agreed.
    inline PropagationResult propagateLabelBySlot(std::array<TemplateBank, 3>& ecg_banks, TemplateBank& ppg_bank, int slot, uint8_t label_code)
    {
        PropagationResult out;
        if (label_code == kUnlabeled || slot < 0) return out;

        auto label = [&](TemplateBank& bk, int32_t* which, int32_t* sub,
            int* relabeled) {
                if (slot >= bk.size()) return;
                template_of_all_signals& t = bk.templates[slot];
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
            const template_of_all_signals& t = bank.templates[i];
            const int idx = (t.label_code != kUnlabeled && t.subtype > 0)
                ? t.subtype - 1 : rank;
            letter[i] = static_cast<uint8_t>(idx % 26);
            ++rank;
        }
        return letter;
    }

    // =====================================================================
    // SCORING PRIMITIVES
    // =====================================================================
    //
    // correlate(), bandMatch() and recomputeTemplate(): measure similarity,
    // measure band containment, rebuild a centroid and its corridor. These
    // are the three operations every assignment decision is made of.
    //
    // They were template_assign.hpp, a file named for a job it no longer did:
    // assignment, spawning, merge/cap-raise and pass-2 refinement over a
    // PER-CHANNEL bank all moved into jbank::runBank when the partition became
    // joint across CH1/CH2/CH3/PPG, because a per-channel assigner is a second
    // partition of the same beats. What was left was three functions in
    // namespace tbank operating on tbank state, so they are here with the
    // state. correlate() is also called by nsvt_detect.hpp to match
    // morphologies ACROSS bins, which is why it is a free function rather than
    // a member of the bank.
    // ---------------------------------------------------------------------
    // Correlation
    // ---------------------------------------------------------------------

    struct CorrResult {
        double r = std::numeric_limits<double>::quiet_NaN();
        int    n_overlap = 0;
        bool   scorable() const { return n_overlap >= kMinOverlapColumns && !std::isnan(r); }
    };

    // A WRAPPER NOW. The five accumulators moved to stats_utils::pearson, the
    // one implementation; what stays here is this bank's two conventions,
    // which are the part that is actually tbank's business:
    //
    //   INCLUSIVE hi, converted to pearson's half-open [lo, hi+1). Callers
    //   pass column indices, and nsvt_detect matches morphologies across bins
    //   with them, so changing the convention would be a silent one-column
    //   change at every call site.
    //
    //   kMinOverlapColumns, the point below which two templates are not
    //   comparable at all. A short overlap leaves r at NaN and reports
    //   n_overlap, so scorable() can tell "too little evidence" from "measured
    //   and different" -- the distinction the partition's diagnostics rest on.
    inline CorrResult correlate(const std::vector<double>& a,
        const std::vector<double>& b, int lo = 0, int hi = std::numeric_limits<int>::max())
    {
        CorrResult out;
        const int w = static_cast<int>(std::min(a.size(), b.size()));
        // hi is inclusive here; guard the +1 against the sentinel default.
        const int hiOpen = (hi >= w) ? w : (hi < 0 ? 0 : hi + 1);
        const PearsonResult pr = pearson(a, b, lo, hiOpen);

        out.n_overlap = pr.n_overlap;
        if (pr.n_overlap < kMinOverlapColumns) return out;   // r stays NaN
        out.r = pr.r;                                        // NaN if undefined
        return out;
    }

    // ---------------------------------------------------------------------
    // The band-match score (Section 4.6's metric)
    // ---------------------------------------------------------------------

    struct BandResult {
        // THE CORRELATION, in [-1, 1]. This is what every floor is compared
        // against; see bandMatch below for why it stopped being a fraction.
        double score = std::numeric_limits<double>::quiet_NaN();
        int    n_overlap = 0;
        double r = std::numeric_limits<double>::quiet_NaN();      // == score
        // Fraction of comparable columns inside the corridor: what `score`
        // used to be. Kept as a DIAGNOSTIC and filled only for the winning
        // template, by assignBeat. Nothing routes on it.
        double frac_in_band = std::numeric_limits<double>::quiet_NaN();
        bool   scorable() const {
            return n_overlap >= kMinOverlapColumns && !std::isnan(score);
        }
    };
    inline BandResult bandMatch(const std::vector<double>& beat,
        const template_of_all_signals& t)
    {
        BandResult out;
        if (t.tmpl.empty()) return out;

        CorrResult cr;
        if (t.corr_halfwin > 0 && t.r_col >= 0) {
            cr = correlate(beat, t.tmpl,
                t.r_col - t.corr_halfwin, t.r_col + t.corr_halfwin);
        }
        else {
            cr = correlate(beat, t.tmpl);
        }
        out.n_overlap = cr.n_overlap;
        out.score = cr.r;
        out.r = cr.r;
        // frac_in_band stays NaN: it costs another pass over the overlap for
        // every beat against every template, and only the winner's value is
        // ever read. assignBeat fills it there.
        return out;
    }

    inline void recomputeTemplate(template_of_all_signals& t,
        const std::vector<std::vector<double>>& beats,
        int width,
        const std::vector<double>* floor_corridor = nullptr)
    {
        const double NaN = std::numeric_limits<double>::quiet_NaN();
        t.tmpl.assign(width, NaN);
        t.tmpl_std.assign(width, NaN);
        t.band_lo.assign(width, NaN);
        t.band_hi.assign(width, NaN);
        t.corridor_inherited = false;
        if (t.members.empty() || width <= 0) return;

        // Own corridor only when there are enough members for percentiles to
        // mean something; otherwise widths come from the floor below.
        const bool own_corridor = t.memberCount() >= kMinMembersForCorridor;
        t.corridor_inherited = !own_corridor;

        std::vector<double> col;
        col.reserve(t.members.size());
        for (int c = 0; c < width; ++c) {
            col.clear();
            for (uint32_t m : t.members) {
                if (m >= beats.size()) continue;
                const auto& b = beats[m];
                if (c < static_cast<int>(b.size()) && !std::isnan(b[c]))
                    col.push_back(b[c]);
            }
            if (col.empty()) continue;
            const size_t n = col.size();

            // nth_element, not sort. One order statistic is needed per column
            // (the median), and this runs once per column per template per
            // recompute -- on a 991-member slot 0 across a 200 column axis that
            // is the single hottest loop in the pass. Partial selection is O(n)
            // against sort's O(n log n) and measured roughly 3x faster here at
            // these sizes.
            //
            // Q1 AND Q3 ARE NO LONGER SELECTED. They existed only to form the
            // spread as q3 - q1; the spread is now the per-sample SD, which
            // needs no order statistic at all. Two of the three nth_element
            // calls are therefore gone.
            const size_t imid = n / 2;
            std::nth_element(col.begin(), col.begin() + imid, col.end());
            const double hi_mid = col[imid];

            if (n % 2) {
                t.tmpl[c] = hi_mid;
            }
            else {
                // Even n: the median averages the two central values, and the
                // lower one is the max of everything below imid -- already
                // partitioned there by the nth_element calls above, so no
                // further selection is needed.
                const double lo_mid =
                    *std::max_element(col.begin(), col.begin() + imid);
                t.tmpl[c] = 0.5 * (lo_mid + hi_mid);
            }
            // PER-SAMPLE SD (ddof = 1), NOT AN IQR, despite the field name.
            // Every spread carried alongside a template in this pipeline is an
            // SD -- create_ecg_templates and normalize_template_amplitude
            // already computed one, and this column used to be a q3 - q1 while
            // feeding the same field, so the field's contents depended on which
            // producer filled it. For a normal column IQR is about 1.35 SD, so
            // the two were roughly 35% apart in the same units and any
            // threshold compared against the field was right for only one of
            // them (see kAutoFootIqrMax).
            //
            // n < 2 leaves 0.0: one observation has no spread, and a NaN here
            // would propagate into every band drawn from it.
            if (n >= 2) {
                double mean = 0.0;
                for (double v : col) mean += v;
                mean /= static_cast<double>(n);
                double sumsq = 0.0;
                for (double v : col) sumsq += (v - mean) * (v - mean);
                t.tmpl_std[c] =
                    std::sqrt(sumsq / static_cast<double>(n - 1));
            }

            // --- the 2.5/97.5 corridor -------------------------------------
            if (own_corridor) {

                if (n < static_cast<size_t>(kMinMembersForCorridor)) {
                    double half = std::numeric_limits<double>::quiet_NaN();
                    if (floor_corridor && c < static_cast<int>(floor_corridor->size()))
                        half = (*floor_corridor)[c];
                    if (std::isnan(half) || half <= 0.0) {
                        double plo = std::numeric_limits<double>::infinity();
                        double phi = -std::numeric_limits<double>::infinity();
                        for (int k = 0; k < width; ++k) {
                            if (std::isnan(t.tmpl[k])) continue;
                            plo = std::min(plo, t.tmpl[k]);
                            phi = std::max(phi, t.tmpl[k]);
                        }
                        const double ptp = (phi > plo) ? (phi - plo) : 0.0;
                        half = kFallbackCorridorFrac * ptp;
                    }
                    if (half > 0.0) {
                        half *= corridorInflation(static_cast<int>(n));
                        t.band_lo[c] = t.tmpl[c] - half;
                        t.band_hi[c] = t.tmpl[c] + half;
                    }
                    continue;
                }

                // Percentiles by position with linear interpolation, over the
                // column's values. nth_element, NOT sort -- for the reason the
                // median block above states and which the first version of this
                // block ignored. A std::sort here runs per column, per
                // template, per recompute, inside the hottest loop in the pass.
                // Four extra order statistics, ordered ascending so each
                // partitions only what the previous left.
                const double x_lo = 0.025 * (static_cast<double>(n) - 1.0);
                const double x_hi = 0.975 * (static_cast<double>(n) - 1.0);
                const size_t i_lo = static_cast<size_t>(std::floor(x_lo));
                const size_t i_lo1 = std::min(n - 1, i_lo + 1);
                const size_t i_hi = static_cast<size_t>(std::floor(x_hi));
                const size_t i_hi1 = std::min(n - 1, i_hi + 1);
                const double f_lo = x_lo - static_cast<double>(i_lo);
                const double f_hi = x_hi - static_cast<double>(i_hi);
                std::nth_element(col.begin(), col.begin() + i_lo, col.end());
                const double p2a = col[i_lo];
                std::nth_element(col.begin() + i_lo, col.begin() + i_lo1, col.end());
                const double p2b = col[i_lo1];
                std::nth_element(col.begin() + i_lo1, col.begin() + i_hi, col.end());
                const double p97a = col[i_hi];
                std::nth_element(col.begin() + i_hi, col.begin() + i_hi1, col.end());
                const double p97b = col[i_hi1];
                const double p2 = p2a * (1.0 - f_lo) + p2b * f_lo;
                const double p97 = p97a * (1.0 - f_hi) + p97b * f_hi;
                const double mid = 0.5 * (p2 + p97);
                const double half = 0.5 * (p97 - p2)
                    * corridorInflation(t.memberCount());
                t.band_lo[c] = mid - half;
                t.band_hi[c] = mid + half;
            }
        }

        // --- inherited corridor ---------------------------------------------
        // Centred on THIS template's own median, widened by slot 0's spread.
        // Centring on the inheritor and not the donor is the whole point: the
        // shape being scored against is this template's, only the tolerance is
        // borrowed. Centring on slot 0 would make every young template a
        // restatement of sinus and no beat would ever fail to match it.
        if (!own_corridor) {
            for (int c = 0; c < width; ++c) {
                if (std::isnan(t.tmpl[c])) continue;
                double half = std::numeric_limits<double>::quiet_NaN();
                if (floor_corridor && c < static_cast<int>(floor_corridor->size()))
                    half = (*floor_corridor)[c];
                if (std::isnan(half) || half <= 0.0) {
                    // Neither this template nor slot 0 has a spread here. Fall
                    // back to a fraction of the template's own amplitude, which
                    // is unitless and survives normalization.
                    double lo = std::numeric_limits<double>::infinity();
                    double hi = -std::numeric_limits<double>::infinity();
                    for (int k = 0; k < width; ++k) {
                        if (std::isnan(t.tmpl[k])) continue;
                        lo = std::min(lo, t.tmpl[k]);
                        hi = std::max(hi, t.tmpl[k]);
                    }
                    const double ptp = (hi > lo) ? (hi - lo) : 0.0;
                    half = kFallbackCorridorFrac * ptp;
                    if (half <= 0.0) continue;   // flat template: leave NaN
                }
                half *= corridorInflation(t.memberCount());
                t.band_lo[c] = t.tmpl[c] - half;
                t.band_hi[c] = t.tmpl[c] + half;
            }
        }
    }

}  // namespace tbank